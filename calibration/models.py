"""
The interface to all yer models.

"""

__all__ = ["Model"]

import os
import time
import logging
import cPickle as pickle

import numpy as np
import psycopg2

from data import SDSSRun
from patch import Patch
from conversions import mag2nmgy, nmgy2mag


# Set up the database connection
_connection = None


def get_db_connection():
    """
    Get the existing database connection or create a new one if needed.

    """
    global _connection
    if _connection is None:
        _connection = psycopg2.connect("dbname='{0}'".format(
                    os.environ.get("SDSS_DB_NAME", "sdss")))
    return _connection


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
    # The name of the table that is associated with the model.
    table_name = None

    # An ordered list of the column names in the table.
    columns = []

    def __init__(self, **kwargs):
        if all([f in kwargs for f in self.columns]):
            # The provided document is _fully specified_.
            self.doc = kwargs
        else:
            logging.warn("{0} is not fully specified by {1}"
                    .format(type(self), kwargs.keys()))

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
        Get the value of a given attribute.

        :param k:
            The name of the attribute to get.

        """
        return self.doc[k]

    @classmethod
    def find(cls, q="", args=[], order=None, limit=None):
        """
        Find a list of objects in the table that match a particular set of
        query parameters.

        This is just a shortcut to run ``SELECT * FROM {{ table }} WHERE
        {{ q }}``.

        :param q: (optional)
            The SQL query to run.

        :param args: (optional)
            The positional arguments for the query.

        :param order: (optional)
            Which column should the results be sorted on?

        :param limit: (optional)
            How many objects should the results be limited to?

        :returns results:
            A list of ``Model`` objects (or subclasses thereof) returned by
            the query. This will return an empty list if nothing was found.

        """
        connection = get_db_connection()
        cursor = connection.cursor()
        cmd = "SELECT {0} FROM {1}".format(",".join(cls.columns),
                                             cls.table_name)
        if q != "":
            cmd += " WHERE " + q
        if order is not None:
            cmd += " ORDER BY {0}".format(order)
        if limit is not None:
            cmd += " LIMIT {0}".format(limit)
        cursor.execute(cmd, args)
        return [cls(**dict(zip(cls.columns, d))) for d in cursor.fetchall()]

    @classmethod
    def find_one(cls, **kwargs):
        """
        Takes the same arguments as :func:`find` but only returns a single
        result. This will return ``None`` if nothing matches the query.

        """
        kwargs["limit"] = 1
        try:
            return cls.find(**kwargs)[0]
        except IndexError:
            return None


class Star(Model):
    table_name = "stars"
    columns = ["id", "has_lightcurve", "ra", "dec"] + "u g r i z".split()

    # def get_lightcurve(self, band=None):
    #     measurements = Measurement.find({"star": self._id,
    #         "out_of_bounds": {"$exists": False}, "band": "g"})

    #     N = len(measurements)
    #     tai = np.empty(N)
    #     flux = np.empty(N)
    #     ferr = np.empty(N)
    #     mask = np.ones(N, dtype=bool)

    #     for i, m in enumerate(measurements):
    #         try:
    #             cal = m.calibrated
    #         except AttributeError:
    #             mask[i] = False
    #         else:
    #             flist = [c["flux"] for c in cal]
    #             flux[i] = np.mean(flist)
    #             ferr[i] = np.sqrt(np.mean([c["ferr"] for c in cal]) ** 2
    #                     + np.var(flist))
    #             tai[i] = m.tai

    #     return tai[mask], flux[mask], ferr[mask]


class Run(Model):
    table_name = "runs"
    columns = ["id", "run", "camcol", "field_min", "field_max", "band",
              "ramin", "ramax", "decmin", "decmax"]

    @classmethod
    def point(cls, ra, dec, q="", **kwargs):
        """
        Find all of the runs that overlap with a given point.

        :param ra:
            The R.A. coordinate to search for.

        :param dec:
            The Dec. coordinate to search for.

        ## Keyword Arguments

        * `q` (dict): Any extra query elements.

        ## Returns

        * `results` (list): A list of `Model` objects returned by the query
          or `None` if nothing was found.

        """
        if q != "":
            q += " AND "
        q += "{0} BETWEEN ramin AND ramax AND {1} BETWEEN decmin AND decmax" \
                .format(ra, dec)
        return cls.find(q=q, **kwargs)

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

        ## Returns

        * `stars` (list): The list of `Star` objects within the run.

        """
        q = """ra BETWEEN {ramin} AND {ramax}
            AND dec BETWEEN {decmin} AND {decmax}""".format(**self.doc)
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
        logging.info(
            "Photometering {0} stars in run ({1.run}, {1.camcol}, {1.band})"\
                .format(len(stars), self))
        for star in stars:
            try:
                m = Measurement.measure(self, star)
            except Exception as e:
                logging.warn(
                    "Run ({0.run}, {0.camcol}, {0.band}) failed with: "
                        .format(self)
                    + str(e))
                break
            m.save()


class Measurement(Model):
    cname = "photometry"
    fields = ["star", "position", "run", "band", "tai", "flux", "bg",
            "dx", "dy"]
    coords = "position"

    @classmethod
    def measure(cls, run, star, clobber=False):
        doc = {"star": star._id, "run": run._id, "band": run.band}

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
            doc["bg"] = {"value": 0, "ivar": 0}
            doc["dx"] = {"value": 0, "ivar": 0}
            doc["dy"] = {"value": 0, "ivar": 0}
            return cls(**doc)

        bg, flux, fx, fy = val
        bg_var, flux_var, fx_var, fy_var = var

        f2 = flux ** 2
        doc["dx"] = {"value": fx / flux,
                      "ivar": 1 / (fx_var / f2 + flux_var * (fx / f2) ** 2)}
        doc["dy"] = {"value": fy / flux,
                      "ivar": 1 / (fy_var / f2 + flux_var * (fy / f2) ** 2)}
        doc["flux"] = {"value": flux, "ivar": 1. / flux_var}
        doc["bg"] = {"value": bg, "ivar": 1. / bg_var}

        return cls(**doc)


class CalibPatch(Model):
    cname = "patches"
    fields = ["runs", "stars", "band", "position", "rng", "maxmag",
            "zero", "mean_flux", "delta2", "beta2", "eta2", "ramin",
            "ramax", "decmin", "decmax", "calibid"]
    coords = "position"

    def get_lightcurve(self, sid):
        """
        Get the calibrated lightcurve for a star with a given `_id`.

        ## Arguments

        * `sid` (int): The star `_id`.

        ## Returns

        * `tai` (numpy.ndarray): The timestamps of the observations in TAI.
        * `flux` (numpy.ndarray): The calibrated lightcurve.
        * `ferr` (numpy.ndarray): The standard deviation of the calibrated
          flux.

        """
        i = self.stars.index(sid)
        m = np.array(self.ivar)[:, i] > 0
        f0 = np.array(self.zero)[m]
        return np.array(self.tai)[m], np.array(self.flux)[m, i] / f0, \
                1 / np.sqrt(np.array(self.ivar)[m, i]) / f0, m

    def get_photometry(self, limit=None):
        if hasattr(self, "flux"):
            return self.runs, self.stars, self.flux, self.ivar, self.fp, \
                    self.ivp

        band = self.band
        ra, dec = self.position
        rng = self.rng
        maxmag = self.maxmag

        q = {"out_of_bounds": {"$exists": False}, "band": band}

        # Construct the box for the bounded search.
        box = [[ra - rng[0], dec - rng[1]], [ra + rng[0], dec + rng[1]]]
        q[Measurement.coords] = {"$within": {"$box": box}}

        # Find the measurements in the box.
        ms = Measurement.find(q, limit=limit)

        stars = set([])
        runs = set([])

        # Get the unique stars and runs.
        for m in ms:
            stars.add(m.star)
            runs.add(m.run)
        stars, runs = list(stars), list(runs)

        # Get the priors on the stellar flux.
        star_prior = Star.find({"_id": {"$in": stars}})
        stars, fp, ivp = [], [], []
        coords = []
        for s in star_prior:
            # Loop over the stars and only take the ones with mgnitudes
            # brighter than the given limit.
            mag = s[band]
            if maxmag < 0 or mag < maxmag:
                stars.append(s._id)
                coords.append([s.ra, s.dec])

                f0 = mag2nmgy(mag)
                fp.append(f0)

                # FIXME: This MAGIC shouldn't be here!
                # sig = 0.5 * (mag2nmgy(mag+0.01)-mag2nmgy(mag-0.01))
                ivp.append(1. / (0.05 * f0) ** 2)

        fp, ivp = np.array(fp), np.array(ivp)

        # Build the data arrays.
        flux = np.zeros((len(runs), len(stars)))
        ivar = np.zeros_like(flux)
        tai = np.zeros_like(flux)

        for m in ms:
            try:
                i = runs.index(m.run)
                j = stars.index(m.star)
            except ValueError:
                pass
            else:
                flux[i, j] = m.flux["value"]
                ivar[i, j] = m.flux["ivar"]
                tai[i, j] = m.tai

        # Deal with epochs with no measurements.
        m = np.sum(ivar, axis=1) > 0
        ivar, flux, tai = ivar[m], flux[m], tai[m]
        runs = [runs[i] for i in range(len(runs)) if m[i]]

        self.runs = runs
        self.stars = stars
        self.flux = flux
        self.ivar = ivar
        self.fp = fp
        self.ivp = ivp
        self.tai = np.array([np.median(tai[i, tai[i, :] > 0])
                                        for i in range(tai.shape[0])])

        return runs, stars, flux, ivar, fp, ivp

    @classmethod
    def calibrate(cls, band, ra, dec, rng, maxmag=22, limit=None, calibid=""):
        """
        Given a band and coordinates, calibrate a patch and return the patch.

        ## Arguments

        * `band` (str): The band to use.
        * `ra`, `dec` (float): The coordinates at the center of the patch.
        * `rng` (tuple): The range of the patch in degrees. This should have
          the for `(ra_range, dec_range)`.

        ## Keyword Arguments

        * `maxmag` (float): The limiting magnitude to use for calibration.
        * `limit` (int or None): Passed to `Measurement.find()`. This should
          probably always be `None` except for testing.

        """
        # Create a new empty `CalibPatch` object.
        doc = dict([(k, None) for k in cls.fields])

        doc["band"] = band
        doc["rng"] = rng
        doc["position"] = [ra, dec]
        doc["ramin"], doc["ramax"] = ra - rng[0], ra + rng[0]
        doc["decmin"], doc["decmax"] = dec - rng[1], dec + rng[1]
        doc["maxmag"] = maxmag
        doc["calibid"] = calibid

        self = cls(**doc)

        # Get the photometry.
        runs, stars, flux, ivar, fp, ivp = self.get_photometry(limit=limit)

        patch = Patch(flux, ivar)
        patch.optimize(fp, ivp)

        self.doc["runs"] = runs
        self.doc["stars"] = stars

        self.doc["zero"] = Binary(pickle.dumps(patch.f0, -1))
        self.doc["mean_flux"] = Binary(pickle.dumps(patch.fs, -1))
        self.doc["delta2"] = Binary(pickle.dumps(patch.d2, -1))
        self.doc["beta2"] = Binary(pickle.dumps(patch.b2, -1))
        self.doc["eta2"] = Binary(pickle.dumps(patch.e2, -1))

        return self

    def save(self):
        super(CalibPatch, self).save()

        # Update the stars.
        for i, _id in enumerate(self.stars):
            _db[Star.cname].update({"_id": _id},
                    {"$push": {"eta2": {"calibid": self.calibid,
                                        "value": self.eta2[i]},
                               "f": {"calibid": self.calibid,
                                     "value": self.mean_flux[i]}}})

        # Update the runs.
        for i, _id in enumerate(self.runs):
            _db[Run.cname].update({"_id": _id},
                    {"$push": {"beta2": {"calibid": self.calibid,
                                         "value": self.beta2[i]},
                               "zero": {"calbid": self.calibid,
                                        "value": self.zero[i],
                                        "patch": self._id,
                                        "ramin": self.ramin,
                                        "ramax": self.ramax,
                                        "decmin": self.decmin,
                                        "decmax": self.decmax}},
                     "$set": {"calibrated": True}})

        # Push the calibrated measurements.
        for i, sid in enumerate(self.stars):
            tai, flux, ferr, mask = self.get_lightcurve(sid)
            for j, r in enumerate(np.arange(len(self.runs))[mask]):
                rid = self.runs[r]
                _db[Measurement.cname].update({"star": sid, "run": rid},
                    {"$push": {"calibrated": {"calibid": self.calibid,
                                              "flux": flux[j],
                                              "ferr": ferr[j]}}})


def _do_photo(doc):
    """
    A simple wapper around the `do_photometry` method defined at the top
    level so that it can be pickled. `multiprocessing` y u so annoying?!?!

    ## Arguments

    * `doc` (dict): The document specifying the `Run` to use.

    """
    run = Run(**doc)
    run.do_photometry()


def _do_calib(doc):
    print "Calibrating:", doc
    s = time.time()
    try:
        p = CalibPatch.calibrate(doc["band"], doc["ra"], doc["dec"],
                doc["rng"], calibid=doc["calibid"])
        p.save()
    except Exception as e:
        print doc, " failed to calibrate\n{0}".format(str(e))
    print doc, " took ", time.time() - s, " seconds to calibrate with ", \
            len(p.stars), " stars in ", len(p.runs), " runs"


if __name__ == "__main__":
    import sys

    star = Star.find({"_id": 8647475120364651147})
    print star[0].get_lightcurve()

    sys.exit(0)

    import hashlib

    from multiprocessing import Pool
    pool = Pool(10)

    band = "g"
    calibid = hashlib.md5(str(time.time())).hexdigest()
    ra = -5.0
    dec = 0.0
    targets = [{"ra": ra, "dec": dec, "rng": [0.8, 0.1], "band": band,
        "calibid": calibid}]
    map(_do_calib, targets)
    sys.exit(0)

    if False:
        import lyrae
        import matplotlib.pyplot as pl
        import time

        # p = CalibPatch.calibrate("g", -2.925071, -0.022093, rng=[0.08, 0.1])
        s = time.time()
        p = CalibPatch.calibrate("g", -3.6, -0.1, rng=[0.08, 0.1])
        p.save()
        print "Calibration took ", time.time() - s, " seconds"
        # p = CalibPatch.calibrate("g", -2.145994, 0.437689, rng=[0.08, 0.1])

        b2, e2 = p.beta2, p.eta2

        cs = np.zeros((len(b2), 4))
        cs[:, -1] = 1 - 0.9 * np.array(b2) / np.max(b2)

        for order, i in enumerate(np.argsort(np.sqrt(e2) *
                                            np.array(p.mean_flux))[::-1]):
            pl.clf()

            star = Star(_id=p.stars[i])
            ra, dec = star.ra, star.dec

            tai, flux, ferr = p.get_lightcurve(p.stars[i])

            # Convert to MJD.
            mjd = tai / 24 / 3600

            T = lyrae.find_period({"g": mjd}, {"g": flux}, {"g": ferr})
            lcmodel = lyrae.get_model(T, mjd, flux, ferr)
            print "Period: ", T, "chi2: ", lcmodel[-1]

            # Plot the light curve with errorbars and alpha weight
            # based on badness of the run.
            pl.errorbar(mjd % T, flux, yerr=ferr, ls="None", capsize=0, lw=2,
                    marker="o", zorder=1, barsabove=False, color="k")
            pl.errorbar(mjd % T + T, flux, yerr=ferr, ls="None", capsize=0,
                    lw=2, marker="o", zorder=1, barsabove=False, color="k")

            # Plot the fit.
            t = np.linspace(0, 2 * T, 5000)
            pl.plot(t, lcmodel[0](t), "--k")

            # Plot the fit stellar flux.
            pl.gca().axhline(p.mean_flux[i], color="k")

            ymin = min(pl.gca().get_ylim()[0], 0)
            pl.ylim(ymin, 2 * p.mean_flux[i] - ymin)
            pl.xlim(0, 2 * T)

            pl.ylabel(r"$f \, [\mathrm{nMgy}]$")
            pl.xlabel(r"$t \, [\mathrm{days}]$")

            pl.title(r"$%d:\quad %s=%.4f$"
                % (p.stars[i], p.band, nmgy2mag(p.mean_flux[i])))
            pl.annotate(("eta = %.5f\nT = %.5f days\n"
                    + "(R.A., Dec.) = (%.4f, %.4f)")
                    % (np.sqrt(e2[i]), T, ra, dec), [1, 1],
                    xycoords="axes fraction", ha="right", va="top")
            pl.savefig("lc0/%03d.png" % order)
