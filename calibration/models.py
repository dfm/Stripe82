__all__ = ["Model", "ModelError"]

import os
import time
import logging
import cPickle as pickle

import numpy as np
import psycopg2

from data import SDSSRun
from patch import Patch
from conversions import mag2nmgy, nmgy2mag


class ModelError(Exception):
    """
    An exception thrown if a model is not *fully specified*.

    """
    pass


def get_db_connection():
    """
    Get the existing database connection or create a new one if needed.

    """
    return psycopg2.connect("dbname='{0}'".format(
                            os.environ.get("SDSS_DB_NAME", "sdss")))


class Model(object):
    """
    The base class for the database access model. The initializer takes a
    set of keyword arguments defining the object. The model needs to be
    fully specified (i.e. all of the columns need to be provided (except
    for the id) otherwise, a :class:`ModelError` is thrown.

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
            raise ModelError("{0} is not fully specified by {1}"
                    .format(type(self), kwargs.keys()))

    def save(self):
        """
        Upsert the object into the database.

        """
        connection = get_db_connection()
        cursor = connection.cursor()

        args = [self.get(k) for k in self.columns]

        # If the document has an id already, do the update... otherwise,
        # do an insert.
        if "id" in self.doc:
            update_cmd = ", ".join(["{0}=%s".format(k)
                                    for k in self.columns])
            cmd = "UPDATE {table_name} SET {update_cmd} WHERE id={doc_id}" \
                    .format(table_name=self.table_name,
                            update_cmd=update_cmd,
                            doc_id=self.get("id"))
        else:
            cmd = """INSERT INTO {table_name} ({columns}) VALUES ({values})
                     RETURNING id""".format(
                             table_name=self.table_name,
                             columns=", ".join(self.columns),
                             values=", ".join(["%s"] * len(self.columns)))

        cursor.execute(cmd, args)
        connection.commit()
        connection.close()

        # Save the id if we ran an insert.
        if not "id" in self.doc:
            self.doc["id"] = cursor.fetchone()[0]

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

        :returns:
            A list of :class:`Model` objects (or subclasses thereof) returned
            by the query. This will return an empty list if nothing was found.

        """
        connection = get_db_connection()
        cursor = connection.cursor()
        cmd = "SELECT id, {0} FROM {1}".format(",".join(cls.columns),
                                             cls.table_name)
        if q != "":
            cmd += " WHERE " + q
        if order is not None:
            cmd += " ORDER BY {0}".format(order)
        if limit is not None:
            cmd += " LIMIT {0}".format(limit)
        cursor.execute(cmd, args)
        results = [cls(**dict(zip(["id"] + cls.columns, d)))
                for d in cursor.fetchall()]
        connection.close()
        return results

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
    """
    Access objects in the ``stars`` table. These are objects from the SDSS
    co-add catalog. This object and the associated table have the following
    attributes:

    .. cssclass:: schema

    * ``id`` (integer primary key): The SDSS id of this star.
    * ``ra`` (real): The R.A. position of the star.
    * ``dec`` (real): The Dec. position of the star.
    * ``u`` (real): The CAS flux in the u-band.
    * ``g`` (real): The CAS flux in the g-band.
    * ``r`` (real): The CAS flux in the r-band.
    * ``i`` (real): The CAS flux in the i-band.
    * ``z`` (real): The CAS flux in the z-band.
    * ``has_lightcurve`` (bool): Does this star have a calibrated light curve?

    """
    table_name = "stars"
    columns = ["has_lightcurve", "ra", "dec"] + "u g r i z".split()

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
    """
    Access objects in the ``runs`` table. This class is essentially a pointer
    to an :class:`calibration.data.SDSSRun` object with a few database
    helpers. One entry corresponds to a ``(run, camcol)`` pair from the SDSS
    catalog. It has the following attributes:

    .. cssclass:: schema

    * ``id`` (integer primary key): The id of this run.
    * ``run`` (integer): The SDSS run number.
    * ``camcol`` (integer): The SDSS camcol number (1-6 inclusive).
    * ``field_min`` (integer): The first field contained in this run.
    * ``field_max`` (integer): The last field contained in this run
      (inclusive).
    * ``band`` (integer): The SDSS filter (0-4 inclusive).
    * ``ramin`` (real): The absolute minimum R.A. hit by the run.
    * ``ramax`` (real): The absolute maximum R.A. hit by the run.
    * ``decmin`` (real): The absolute minimum Dec. hit by the run.
    * ``decmax`` (real): The absolute maximum Dec. hit by the run.

    """
    table_name = "runs"
    columns = ["run", "camcol", "field_min", "field_max", "band",
              "ramin", "ramax", "decmin", "decmax"]

    @classmethod
    def point(cls, ra, dec, q="", **kwargs):
        """
        Find all of the runs that overlap with a given point.

        :param ra:
            The R.A. coordinate to search for.

        :param dec:
            The Dec. coordinate to search for.

        :param q: (optional)
            An additional SQL query that results must satisfy.

        :returns:
            A list of :class:`Run` objects returned by the query or
            ``None`` if nothing was found.

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

        :returns:
            A :class:`calibration.data.SDSSRun` object pointing to the right
            data file.

        """
        try:
            return self._data
        except AttributeError:
            self._data = SDSSRun(self.run, self.camcol, "ugriz"[self.band])
            return self._data

    def get_stars(self):
        """
        Get all the stars within the "bounds" of this run as defined by
        ``ramin``, ``ramax``, ``decmin`` and ``decmax``.

        :returns:
            The list of :class:`Star` objects within the bounds of the run.

        """
        q = """ra BETWEEN {ramin} AND {ramax}
            AND dec BETWEEN {decmin} AND {decmax}""".format(**self.doc)
        stars = Star.find(q)
        return stars

    def do_photometry(self):
        """
        Do the photometry for all of the stars in the run.

        """
        s0 = self.get_stars()
        sids = [s["id"] for s in s0]
        done = Measurement.find(q="starid IN ({0}) AND runid=%s".format(
            ",".join([str(s) for s in sids])), args=[self.get("id")])
        dids = [m.starid for m in done]
        stars = []
        for s in s0:
            if s["id"] not in dids:
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
    """
    Access objects in the ``raw`` table. These objects are raw photometric
    measurements of point sources. The object as the following attributes:

    .. cssclass:: schema

    * ``id`` (integer primary key): The id of this measurement.
    * ``runid`` (integer): The associated [run](/models/runs).
    * ``starid`` (integer): The associated [star](/models/stars).
    * ``ra`` (real): The R.A. position of the measurement (based on the
      associated :class:`Star` model.
    * ``dec`` (real): The Dec. position of the measurement.
    * ``tai`` (real): The time of the measurement (based on the
      interpolation of times associated with the :class:`Run` model). The
      units are seconds.
    * ``band`` (integer): The SDSS filter.
    * ``out_of_bounds`` (bool): Is the star *actually* out of the bounds of
      the run despite the estimated cut based on ``ramin``, ``ramax``,
      ``decmin`` and ``decmax``.
    * ``flux`` (real): The measured counts of the source.
    * ``fluxivar`` (real): The inverse variance in ``flux``.
    * ``sky`` (real): The measured background sky level.
    * ``skyivar`` (real): The inverse variance in ``sky``.
    * ``dx`` (real): The offset of the center of the star along the "x"
      coordinate measured in pixels.
    * ``dxivar`` (real): The inverse variance in ``dx``.
    * ``dy`` (real): The offset of the center of the star along the "y"
      coordinate measured in pixels.
    * ``dyivar`` (real): The inverse variance in ``dy``.

    """
    table_name = "raw"
    columns = ["starid", "runid", "ra", "dec", "tai", "band", "out_of_bounds",
               "flux", "fluxivar", "sky", "skyivar", "dx", "dxivar",
               "dy", "dyivar"]

    @classmethod
    def measure(cls, run, star):
        """
        Measure the photometry of a particular star in a particular run.

        :param run:
            A :class:`Run` object.

        :param star:
            A :class:`Star` object.

        :returns:
            The :class:`Measurement` object.

        """
        doc = {"starid": star.get("id"), "runid": run.get("id"),
                "band": run.band}

        ra, dec = star.ra, star.dec
        while ra < 0:
            ra += 360.

        doc["ra"] = star.ra
        doc["dec"] = star.dec
        doc["tai"] = run.data.get_tai(ra, dec)

        try:
            val, var = run.data.photometry(ra, dec)
        except IndexError:
            doc["out_of_bounds"] = True
            for k in ["flux", "fluxivar", "sky", "skyivar", "dx", "dxivar",
                    "dy", "dyivar"]:
                doc[k] = 0.0
            return cls(**doc)

        doc["out_of_bounds"] = False

        bg, flux, fx, fy = val
        bg_var, flux_var, fx_var, fy_var = var

        f2 = flux * flux
        doc["dx"] = fx / flux
        doc["dxivar"] = 1. / (fx_var / f2 + flux_var * (fx / f2) ** 2)
        doc["dy"] = fy / flux,
        doc["dyivar"] = 1. / (fy_var / f2 + flux_var * (fy / f2) ** 2)
        doc["flux"] = flux
        doc["fluxivar"] = 1. / flux_var
        doc["sky"] = bg
        doc["skyivar"] = 1. / bg_var

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


def _build_indices():
    print "Building indices:"
    connection = get_db_connection()
    cursor = connection.cursor()

    # Measurements need indexes on the run and star ids.
    print "... Raw photometry"
    cursor.execute("CREATE INDEX ON raw (starid);")
    cursor.execute("CREATE INDEX ON raw (runid);")

    connection.commit()
    connection.close()


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
