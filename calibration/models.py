"""
Interface to the MongoDB persistent backend.

"""

__all__ = ["Model", "SDSSModel"]

import os
import logging

import numpy as np
import pymongo

from db import Database
from data import SDSSRun
from patch import Patch
from conversions import *

_db = Database(name=os.environ.get("MONGO_DB", "sdss"))

# Ensure indices.
_db.stars.ensure_index([("coords", pymongo.GEO2D)])

_db.runs.ensure_index("decMin")
_db.runs.ensure_index("decMax")

_db.photometry.ensure_index([("position", pymongo.GEO2D)])

class Model(object):
    """
    The base class for the database access model. The initializer takes a
    set of keyword arguments defining the object. If the keywords _fully_
    specify the object (i.e. all of the default fields and an `_id` field
    are provided) then the database is not queried. If the document is
    fully specified except for an `_id`, the database is queried to see
    if such an object already exists and if it doesn't a new document is
    created in memory but not committed until `save`. Finally, if the
    keyword arguments do not fully specify the model then a query is run to
    find the document satisfying the keyword arguments.

    """
    # The name of the collection that is associated with the model is given
    # by the attribute `cname` that must be implemented by subclasses.
    cname  = None

    # A list of the default fields to always load for every instance.
    fields = []

    # For geospatially indexed collections, `coords` gives the name of the
    # field that is geospatially indexed.
    coords = None

    def __init__(self, **kwargs):
        if all([f in kwargs for f in self.fields]):
            # The provided document is _fully specified_.
            self.doc = kwargs
            if "_id" not in kwargs:
                # Query the database to find the `_id` of this object if it
                # exists. If not, just leave `_id` blank and let it be
                # generated on `save`.
                d = self.collection.find_one(kwargs, {"_id": 1})
                if d is not None:
                    self.doc["_id"] = d["_id"]
        else:
            # Query the database using the provided document as the search.
            fields = dict([(k, 1) for k in self.fields])
            self.doc = self.collection.find_one(kwargs, fields)
            assert self.doc is not None

    @property
    def collection(self):
        return _db[self.cname]

    def save(self):
        """
        Commit the current state of the object to the database and update
        the `_id` if necessary.

        """
        self.doc["_id"] = self.collection.save(self.doc)

    def __getitem__(self, k):
        return self.get(k)

    def __getattr__(self, k):
        try:
            return self.get(k)
        except KeyError:
            raise AttributeError("'%s' object has not attribute '%s'"
                    % (type(self), k))

    def get(self, k):
        """
        Get the value of a given attribute for the object and query the
        database if needed.

        ## Arguments

        * `k` (str): The attribute to get.

        """
        try:
            return self.doc[k]
        except KeyError:
            d = self.collection.find_one({"_id": self._id}, {"_id":0, k: 1})
            self.doc[k] = d[k]
            return d[k]

    def get_multiple(self, f):
        """
        Get the values for a set of attributes of the object and query the
        database if needed.

        ## Arguments

        * `f` (list): The attributes to get.

        ## Returns

        * `values` (dict): A dictionary of the requested attributes.

        """
        result = {}
        get_f  = {"_id": 0}
        for k in f:
            try:
                result[k] = self.doc[k]
            except KeyError:
                get_f[k]  = 1
        if len(get_f) > 1:
            d = self.collection.find_one({"_id": self._id}, get_f)
            for k in get_f:
                self.doc[k] = d[k]
                result[k] = d[k]
        return result

    @classmethod
    def find(cls, q, limit=None):
        """
        Run a query on the model collection and return the resulting objects.

        ## Arguments

        * `q` (dict): The query to run.

        ## Keyword Arguments

        * `limit` (int): How many documents should the results be limited to?

        ## Returns

        * `results` (list): A list of `Model` objects returned by the query
          or `None` if nothing was found.

        """
        f = dict([(k, 1) for k in cls.fields])
        c = _db[cls.cname].find(q, f)
        if c is None:
            return None
        if limit is not None:
            c = c.limit(limit)
        return [cls(**d) for d in c]

    @classmethod
    def sphere(cls, center, radius=None, q={}, **kwargs):
        """
        For a collection with a geospatial index, search in spherical
        coordinates (as specified by the `coords` attribute).

        ## Arguments

        * `center` (tuple): The `(ra, dec)` coordinates for the center of the
          search.

        ## Keyword Arguments

        * `radius` (float): The radius (in arcminutes) to search within.
        * `q` (dict): Any extra query elements.

        ## Returns

        * `results` (list): A list of `Model` objects returned by the query
          or `None` if nothing was found.

        """
        assert cls.coords is not None

        if radius is not None:
            radius = np.radians(radius/60.0)
            while center[0] > 180.:
                center[0] -= 360.0
            q[cls.coords] = {"$within": {"$centerSphere": [center,radius]}}
        else:
            q[cls.coords] ={"$nearSphere": center}
        return cls.find(q, **kwargs)

class Star(Model):
    cname  = "stars"
    fields = ["ra", "dec"] + "u g r i z".split()
    coords = "coords"

class Run(Model):
    cname  = "runs"
    fields = ["run", "camcol", "band", "raMin", "raMax", "decMin", "decMax"]

    @classmethod
    def point(cls, pt, q={}, **kwargs):
        """
        Find all of the runs that overlap with a given point.

        WARNING: This doesn't actually match the R.A. right now because the
            `preprocess` script is _broken_ in how it deals with `raMin` and
            `raMax`.

        ## Arguments

        * `pt` (tuple): The coordinates `(ra, dec)` to search for.

        ## Keyword Arguments

        * `q` (dict): Any extra query elements.

        ## Returns

        * `results` (list): A list of `Model` objects returned by the query
          or `None` if nothing was found.

        """
        ra, dec = pt
        q["decMin"] = {"$lt": dec}
        q["decMax"] = {"$gt": dec}
        return cls.find(q, **kwargs)

    @property
    def data(self):
        """
        Lazily access the interface to the HDF5 data file.

        ## Returns

        * `data` (SDSSRun): The data access interface object.

        """
        try:
            return self._data
        except AttributeError:
            self._data = SDSSRun(self.run, self.camcol, self.band)
            return self._data

    def get_stars(self):
        """
        Get all the stars within the bounds of this run.

        WARNING: This only matches the `dec` value.

        ## Returns

        * `stars` (list): The list of `Star` objects within the run.

        """
        q = {"dec": {"$gt": self.decMin, "$lt": self.decMax}}
        stars = Star.find(q)
        return stars

    def do_photometry(self):
        """
        Do the photometry for all of the stars in a given run.

        """
        s0 = self.get_stars()
        sids = [s._id for s in s0]
        done = Measurement.find({"star": {"$in": sids}, "run": self._id})
        dids = [m.star for m in done]
        stars = []
        for s in s0:
            if s._id not in dids:
                stars.append(s)
        print "Photometering %d stars in run (%d, %d, %s)"\
                % (len(stars), self["run"], self["camcol"], self["band"])
        for star in stars:
            try:
                m = Measurement.measure(self, star)
            except Exception as e:
                logging.warn("Run (%d, %d, %s) failed with: "
                        %(self["run"], self["camcol"], self["band"]) + str(e))
                break
            m.save()

class Measurement(Model):
    cname  = "photometry"
    fields = ["star", "position", "run", "tai", "flux", "bg", "dx", "dy"]
    coords = "position"

    @classmethod
    def measure(cls, run, star, clobber=False):
        doc = {"star": star._id, "run": run._id}

        ra, dec = star.ra, star.dec
        while ra < 0:
            ra += 360.

        doc[cls.coords] = [star.ra, star.dec]
        doc["tai"] = run.data.get_tai(ra, dec)

        try:
            val, var = run.data.photometry(ra, dec)
        except IndexError:
            doc["out_of_bounds"] = True
            doc["flux"] = {"value": 0, "ivar": 0}
            doc["bg"]   = {"value": 0, "ivar": 0}
            doc["dx"]   = {"value": 0, "ivar": 0}
            doc["dy"]   = {"value": 0, "ivar": 0}
            return cls(**doc)

        bg, flux, fx, fy = val
        bg_var, flux_var, fx_var, fy_var = var

        doc["dx"] = {"value": fx/flux,
                      "ivar": 1./(fx_var/flux**2 + flux_var*(fx/flux**2)**2)}
        doc["dy"] = {"value": fy/flux,
                      "ivar": 1./(fy_var/flux**2 + flux_var*(fy/flux**2)**2)}
        doc["flux"] = {"value": flux, "ivar": 1./flux_var}
        doc["bg"]   = {"value": bg, "ivar": 1./bg_var}

        return cls(**doc)

class CalibPatch(Model):
    cname  = "patches"
    fields = []
    coords = "position"

    @classmethod
    def calibrate(cls, band, ra, dec, radius=None, limit=None):
        # FIXME: Add band specification in the query here.
        ms = Measurement.sphere([ra, dec], radius=radius, limit=limit,
                q={"out_of_bounds": {"$exists": False}})

        stars = set([])
        runs  = set([])

        # Get the unique stars and runs.
        for m in ms:
            stars.add(m.star)
            runs.add(m.run)
        stars, runs = list(stars), list(runs)

        # Build the data arrays.
        flux = np.zeros((len(runs), len(stars)))
        ivar = np.zeros_like(flux)

        for m in ms:
            i = runs.index(m.run)
            j = stars.index(m.star)
            flux[i, j] = m.flux["value"]
            ivar[i, j] = m.flux["ivar"]

        # Get the priors.
        star_prior = Star.find({"_id": {"$in": stars}})
        fp = np.zeros(len(stars))
        ivp = 10**(-0.4 * (0.5 - 22.5))*np.ones_like(fp)

        for s in star_prior:
            i = stars.index(s._id)
            mag = s[band]
            f0  = mag2nmgy(mag)
            sig = 0.5 * (mag2nmgy(mag+0.01)-mag2nmgy(mag-0.01))

            fp[i] = f0
            ivp[i] = 1./(0.05*f0)**2

        print fp

        patch = Patch(flux, ivar)
        p = patch.optimize(fp, ivp)

        print patch.f0
        print patch.fs
        print np.sqrt(patch.w2)
        print np.sqrt(patch.s2)

        import matplotlib.pyplot as pl

        for si, s in enumerate(star_prior):
            pl.figure()
            pl.plot(np.arange(len(flux[:,si])), flux[:,si]/patch.f0, ".k")
            pl.gca().axhline(fp[si])
            pl.gca().axhline(patch.fs[si])
            pl.ylim(min(pl.gca().get_ylim()[0], 0), 2*fp[si])
            pl.title("%.3e"%np.sqrt(patch.w2[si]))
            pl.savefig("lc/%d.png"%si)

# def _do_photo(doc):
#     run = Run(**doc)
#     run.do_photometry()

if __name__ == "__main__":
    import sys
    from multiprocessing import Pool

    # if "--photo" in sys.argv:
    #     runs = [r.doc for r in Run.find({"band": "g"})]
    #     pool = Pool(10)
    #     pool.map(_do_photo, runs)

    p = CalibPatch.calibrate("g", -5, 0.5, radius=1)

