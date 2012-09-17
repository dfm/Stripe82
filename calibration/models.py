__all__ = ["Model", "ModelError"]

import time
import logging
import datetime

import numpy as np

from data import SDSSRun
from patch import Patch
from conversions import mag2nmgy, nmgy2mag
from db import DBConnection


class ModelError(Exception):
    """
    An exception thrown if a model is not *fully specified*.

    """
    pass


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

    @property
    def _save_cmd(self):
        args = [self.get(k) for k in self.columns]

        # If the document has an id already, do the update... otherwise,
        # do an insert.
        if "id" in self.doc:
            update_cmd = ", ".join(["{0}=%s".format(k)
                                    for k in self.columns])
            cmd = "UPDATE {table_name} SET {update_cmd} WHERE id={doc_id}"
            cmd = cmd.format(table_name=self.table_name,
                                update_cmd=update_cmd,
                                doc_id=self.get("id"))
        else:
            cmd = """INSERT INTO {table_name} ({columns}) VALUES ({values})
                    RETURNING id""".format(
                            table_name=self.table_name,
                            columns=", ".join(self.columns),
                            values=", ".join(["%s"] * len(self.columns)))
        return cmd, args

    def save(self):
        """
        Upsert the object into the database.

        """
        with DBConnection() as cursor:
            cmd, args = self._save_cmd
            cursor.execute(cmd, args)

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
        with DBConnection() as cursor:
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

        # Figure out which measurements have already been made.
        q = {
                "q": "starid IN ({0}) AND runid=%s".format(",".join([str(s)
                                                            for s in sids])),
                "args": [self.get("id")]
            }

        done = Measurement.find(**q) + OutOfBoundsMeasurement.find(**q)
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


class OutOfBoundsMeasurement(Model):
    """
    This is a "placeholder" model. When a :class:`Measurement` would be out
    of bounds, it is added to the ``out_of_bounds`` table.

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

    """
    table_name = "out_of_bounds"
    columns = ["starid", "runid", "ra", "dec", "tai", "band"]


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
    columns = ["starid", "runid", "ra", "dec", "tai", "band",
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
            return OutOfBoundsMeasurement(**doc)

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


class CalibRun(Model):
    """
    Access objects in the ``calibruns`` table. These objects provide metadata
    for a particular pass at calibration.

    .. cssclass:: schema

    * ``id`` (integer primary key): The id of this pass.
    * ``start_date`` (timestamp): The start date of the pass.
    * ``band`` (integer): The SDSS band calibrated in this pass.

    """
    table_name = "calibruns"
    columns = ["start_date", "band"]

    @classmethod
    def new(cls, band):
        doc = {"band": band, "start_date": datetime.datetime.now()}
        self = cls(**doc)
        return self


class CalibPatch(Model):
    """
    Access objects in the ``patches`` table. These are spatially defined
    patches on the sky that hit particular :class:`Run` and :class:`Star`
    objects.

    .. cssclass:: schema

    * ``id`` (integer primary key): The id of this patch.
    * ``calibid`` (integer): The id of the particular :class:`CalibRun`.
    * ``ramin`` (real): The lower bound in R.A. for the patch.
    * ``ramax`` (real): The upper bound in R.A. for the patch.
    * ``decmin`` (real): The lower bound in Dec. for the patch.
    * ``decmax`` (real): The upper bound in Dec. for the patch.
    * ``stars`` (integer[]): The list of :class:`Star` objects in the patch.
    * ``runs`` (integer[]): The list of :class:`Run` objects that hit the
      patch.

    """
    table_name = "patches"
    columns = ["runs", "stars", "ramin", "ramax", "decmin", "decmax",
               "calibid"]

    # def get_lightcurve(self, sid):
    #     """
    #     Get the calibrated lightcurve for a star with a given `_id`.

    #     ## Arguments

    #     * `sid` (int): The star `_id`.

    #     ## Returns

    #     * `tai` (numpy.ndarray): The timestamps of the observations in TAI.
    #     * `flux` (numpy.ndarray): The calibrated lightcurve.
    #     * `ferr` (numpy.ndarray): The standard deviation of the calibrated
    #       flux.

    #     """
    #     i = self.stars.index(sid)
    #     m = np.array(self.ivar)[:, i] > 0
    #     f0 = np.array(self.zero)[m]
    #     return np.array(self.tai)[m], np.array(self.flux)[m, i] / f0, \
    #             1 / np.sqrt(np.array(self.ivar)[m, i]) / f0, m

    def get_photometry(self, maxmag):
        """
        Get the photometric measurements within the bounds of this patch.

        :param maxmag:
            The maximum limit on the magnitudes of the stars.

        """
        band = self.band
        ramin, ramax = self.ramin, self.ramax
        decmin, decmax = self.decmin, self.decmax

        # Construct the box for the bounded search.
        q = "band = %s"
        q += " AND ra BETWEEN %s AND %s AND dec BETWEEN %s AND %s"

        # Also, limit the range of pixel offsets allowed (<2.5px).
        q += " AND dx * dx + dy * dy < 6.25"

        args = [band, ramin, ramax, decmin, decmax]

        # Find the measurements in the box.
        ms = Measurement.find(q=q, args=args)

        stars = set([])
        runs = set([])

        # Get the unique stars and runs.
        for m in ms:
            stars.add(m.starid)
            runs.add(m.runid)
        stars, runs = list(stars), list(runs)

        # Get the priors on the stellar flux.
        star_prior = Star.find(q="id IN ({0})".format(",".join([str(s)
                                                          for s in stars])))
        stars, fp, ivp = [], [], []
        coords = []
        for s in star_prior:
            # Loop over the stars and only take the ones with magnitudes
            # brighter than the given limit.
            mag = s["ugriz"[band]]
            if maxmag < 0 or mag < maxmag:
                stars.append(s["id"])
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
                i = runs.index(m.runid)
                j = stars.index(m.starid)
            except ValueError:
                pass
            else:
                flux[i, j] = m.flux
                ivar[i, j] = m.fluxivar
                tai[i, j] = m.tai

        # Deal with epochs with no measurements.
        m = np.sum(ivar, axis=1) > 0
        ivar, flux, tai = ivar[m], flux[m], tai[m]
        runs = [runs[i] for i in range(len(runs)) if m[i]]

        # Compute the array of times.
        tai = np.array([np.median(tai[i, tai[i, :] > 0])
                                        for i in range(tai.shape[0])])

        return tai, runs, stars, flux, ivar, fp, ivp

    @classmethod
    def calibrate(cls, band, ra, dec, rng, calibid, maxmag=22, limit=None):
        """
        Given a band and coordinates, calibrate a patch and return the patch.

        :param band:
            The band to use.

        :param ra:
            The R.A. of the center of the patch.

        :param dec:
            The Dec. of the center of the patch.

        :param rng:
            A tuple giving the range of the patch in degrees. This should
            have the form ``(ra_range, dec_range)``.

        :param calibid:
            The ``id`` of the :class:`CalibRun`.

        :param limit: (optional)
            Passed to :func:`Measurement.find()`. This should probably always
            be ``None`` except for testing.

        :param maxmag: (optional)
            The limiting magnitude to use for calibration.

        :param limit: (optional)
            Passed to :func:`Measurement.find()`. This should probably always
            be ``None`` except for testing.

        """
        # Create a new empty `CalibPatch` object.
        doc = dict([(k, None) for k in cls.columns])

        doc["band"] = band
        doc["ramin"], doc["ramax"] = ra - rng[0], ra + rng[0]
        doc["decmin"], doc["decmax"] = dec - rng[1], dec + rng[1]
        doc["calibid"] = calibid

        self = cls(**doc)

        # Get the photometry.
        tai, runs, stars, flux, ivar, fp, ivp = self.get_photometry(maxmag)

        self.doc["runs"] = runs
        self.doc["stars"] = stars
        self.save()
        _id = self["id"]

        patch = Patch(flux, ivar)
        patch.optimize(fp, ivp)

        with DBConnection() as cursor:
            # Build the photometrically calibrated models too.
            zero = patch.f0
            for i, run in enumerate(runs):
                for j, star in enumerate(stars):
                    if ivar[i, j] > 0 and zero[i] > 0.0:
                        doc = {"calibid": calibid, "patchid": _id,
                                "band": band, "runid": run, "starid": star}
                        doc["tai"] = tai[i]
                        doc["flux"] = flux[i, j] / zero[i]
                        doc["fluxivar"] = ivar[i, j] * zero[i] ** 2
                        p = Photometry(**doc)
                        cursor.execute(*p._save_cmd)

            # Save the zero points.
            beta2 = patch.b2
            delta2 = patch.d2
            for i, run in enumerate(runs):
                doc = {"calibid": calibid, "patchid": _id, "band": band,
                    "runid": run, "ramin": self.ramin, "ramax": self.ramax,
                    "decmin": self.decmin, "decmax": self.decmax}
                doc["zero"] = zero[i]
                doc["beta2"] = beta2[i]
                doc["delta2"] = delta2[i]
                z = Zero(**doc)
                cursor.execute(*z._save_cmd)

            # Save the mean fluxes.
            fs = patch.fs
            eta2 = patch.e2
            for j, star in enumerate(stars):
                doc = {"calibid": calibid, "patchid": _id, "band": band,
                        "starid": star, "mean_flux": fs[j], "eta2": eta2[j]}
                f = Flux(**doc)
                cursor.execute(*f._save_cmd)

        return self


class Photometry(Model):
    """
    Access objects in the ``photometry`` table. These are the calibrated
    versions of the objects in the :class:`Measurement`.

    .. cssclass:: schema

    * ``id`` (integer primary key): The id of this calibrated measurement.
    * ``calibid`` (integer): The id of the associated :class:`CalibRun`.
    * ``patchid`` (integer): The id of the associated :class:`CalibPatch`.
    * ``runid`` (integer): The associated :class:`Run`.
    * ``starid`` (integer): The associated :class:`Star`.
    * ``tai`` (real): The time of the measurement (based on the interpolation
      of times associated with the :class:`Run` model). The units are seconds.
    * `band` (integer): The SDSS filter.
    * `flux` (real): The calibrated photometry of the source.
    * `fluxivar` (real): The inverse variance in `flux`.

    """
    table_name = "photometry"
    columns = ["calibid", "patchid", "runid", "starid", "tai", "band",
               "flux", "fluxivar"]


class Zero(Model):
    """
    Access objects in the ``zeros`` table. These are the estimates of the
    photometric zero points of the runs.

    .. cssclass:: schema

    * ``id`` (integer primary key): The id of this zero point.
    * ``calibid`` (integer): The id of the associated :class:`CalibRun`.
    * ``patchid`` (integer): The id of the associated :class:`CalibPatch`.
    * ``runid`` (integer): The associated :class:`Run`.
    * ``ramin`` (real): The minimum R.A. bound of the patch.
    * ``ramax`` (real): The maximum R.A. bound of the patch.
    * ``decmin`` (real): The minimum Dec. bound of the patch.
    * ``decmax`` (real): The maximum Dec. bound of the patch.
    * ``band`` (integer): The SDSS filter.
    * ``zero`` (real): The actual value of the zero point.
    * ``beta2`` (real): The relative variability parameter found for the run.
    * ``delta2`` (real): The absolute variability parameter found for the run.

    """
    table_name = "zeros"
    columns = ["calibid", "patchid", "runid", "ramin", "ramax", "decmin",
               "decmax", "band", "zero", "beta2", "delta2"]


class Flux(Model):
    """
    Access objects in the ``fluxes`` table. These are the estimates of the
    mean fluxes and variability of the stars.

    .. cssclass:: schema

    * `id` (integer primary key): The id of this calibrated measurement.
    * `calibid` (integer): The id of the associated :class:`CalibRun`.
    * `patchid` (integer): The id of the associated :class:`CalibPatch`.
    * `starid` (integer): The associated :class:`Star`.
    * `band` (integer): The SDSS filter.
    * `mean_flux` (real): The inferred mean flux of the source.
    * `eta2` (real): The relative variability parameter found for the star.

    """
    table_name = "fluxes"
    columns = ["calibid", "patchid", "starid", "band", "mean_flux", "eta2"]


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
    print("Calibrating: {0}".format(doc))
    s = time.time()
    try:
        p = CalibPatch.calibrate(doc["band"], doc["ra"], doc["dec"],
                doc["rng"], calibid=doc["calibid"])
        p.save()
    except Exception as e:
        print("{0} failed to calibrate\n{1}".format(doc, str(e)))
    print("{0} took {1} seconds to calibrate with {2} stars in {3} runs"
            .format(doc, time.time() - s, len(p.stars), len(p.runs)))


def _build_indices():
    print "Building indices:"

    with DBConnection() as cursor:
        # Measurements need indexes on the run and star ids.
        print "... Raw photometry"
        cursor.execute("CREATE INDEX ON raw (starid);")
        cursor.execute("CREATE INDEX ON raw (runid);")

        cursor.execute("CREATE INDEX ON out_of_bounds (starid);")
        cursor.execute("CREATE INDEX ON out_of_bounds (runid);")


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
