"""
Interface to the MongoDB persistent backend.

"""

__all__ = ["Model", "SDSSModel"]

import os

import numpy as np
import pymongo

from db import Database

_db = Database(name=os.environ.get("MONGO_DB", "sdss"))

# Ensure indices.
_db.stars.ensure_index([("coords", pymongo.GEO2D)])

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
        print self.doc.get("_id")
        self.doc["_id"] = self.collection.save(self.doc)
        print self.doc.get("_id")

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

