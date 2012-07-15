#!/usr/bin/env python
# encoding: utf-8
"""
Utilities for preprocessing a batch of fields and combining them into one.

"""

import logging
logging.basicConfig(level=logging.INFO)

import os
import shutil
import subprocess
import multiprocessing
import pymongo

import numpy as np

import pyfits
import h5py

import requests

import pysdss


# Data access variables.
_http_root = os.environ.get("SDSS_HTTP_ROOT", None)
_use_http = os.environ.get("USE_HTTP", "False").upper() == "TRUE"
_remote_server = os.environ["SDSS_SERVER"]
_remote_data_dir = os.environ["SDSS_REMOTE"]
_local_tmp_dir = os.path.join(os.environ.get("SDSS_LOCAL", "."), ".sdss")
_local_data_dir = os.path.join(os.environ.get("SDSS_LOCAL", "."), "data")

# Make sure that the data directory exists.
try:
    os.makedirs(_local_data_dir)
except os.error as e:
    if not os.path.exists(_local_data_dir):
        raise e

# Field dimensions.
_f_width = 1489
_f_height = 2048
_f_overlap = 128

# HDF5 file tags.
IMG_TAG = "img"
INV_TAG = "inv"
TAI_TAG = "tai"
BOUNDS_TAG = "bounds"
CENTERS_TAG = "centers"
PSF_TAG = "psf"
EIGEN_TAG = "eigen"
INFO_TAG = "info"
TS_TAG = "ts"
AST_TAG = "ast"

# Connect to the database.
_db_server = os.environ.get("MONGO_SERVER", "localhost")
_db_port = int(os.environ.get("MONGO_PORT", 27017))
_db_name = os.environ.get("MONGO_DB", "sdss")
_db = pymongo.Connection(_db_server, _db_port)[_db_name]
_run_collection = _db.runs
_run_collection.ensure_index([("run",    pymongo.ASCENDING),
                              ("camcol", pymongo.ASCENDING),
                              ("band",   pymongo.ASCENDING)])


class SDSSFileError(Exception):
    pass


class _DR7(pysdss.DR7):
    def __init__(self, *args, **kwargs):
        super(_DR7, self).__init__(*args, **kwargs)
        self.filenames = {
            "fpObjc":  "fpObjc-%(run)06i-%(camcol)i-%(field)04i.fit",
            "fpM":     "fpM-%(run)06i-%(band)s%(camcol)i-%(field)04i.fit.gz",
            "fpC":     "fpC-%(run)06i-%(band)s%(camcol)i-%(field)04i.fit.gz",
            "fpAtlas": "fpAtlas-%(run)06i-%(camcol)i-%(field)04i.fit",
            "psField": "psField-%(run)06i-%(camcol)i-%(field)04i.fit",
            "tsObj":   "tsObj-%(run)06i-%(camcol)i-%(rerun)i-%(field)04i.fit",
            "tsField": \
                    "tsField-%(run)06i-%(camcol)i-%(rerun)i-%(field)04i.fit",
            }
        if _use_http:
            self.filenames["fpM"] = \
                    "fpM-%(run)06i-%(band)s%(camcol)i-%(field)04i.fit"

    def _fullpath(self, filetype, *args, **kwargs):
        for k, v in zip(["run", "camcol", "field", "rerun", "band"], args):
            kwargs[k] = v
        fn = self.filenames[filetype] % kwargs

        if _use_http:
            prefix = '%(run)i/%(rerun)i/' % kwargs
        else:
            prefix = '%(rerun)i/%(run)i/' % kwargs
        if filetype in ['fpC']:
            return prefix + 'corr/%(camcol)i/' % kwargs + fn
        elif filetype in ['psField', 'fpAtlas', 'fpObjc', 'fpM']:
            return prefix + 'objcs/%(camcol)i/' % kwargs + fn
        elif filetype in ['tsObj', 'tsField']:
            return prefix + 'calibChunks/%(camcol)i/' % kwargs + fn
        else:
            return None

    def _open(self, fn):
        local_path = os.path.join(_local_tmp_dir, fn)
        return pyfits.open(local_path)

    def fetch(self, targets):
        try:
            os.makedirs(_local_tmp_dir)
        except os.error as e:
            if not os.path.exists(_local_tmp_dir):
                raise e

        objs = []
        for i, o in enumerate(targets):
            k = dict(zip(["run", "camcol", "field", "rerun", "band"], o[1:]))
            f = os.path.join(_local_tmp_dir, self.filenames[o[0]] % k)
            if not os.path.exists(f):
                objs += [o]

        if len(objs) > 0:
            if _use_http:
                assert _http_root is not None
                for o in objs:
                    fn = self._fullpath(*o)
                    remote_path = os.path.join(_http_root, fn)
                    r = requests.get(remote_path)

                    if r.status_code != 200:
                        raise SDSSFileError("Couldn't copy %s." %
                                (remote_path))

                    f = open(os.path.join(_local_tmp_dir,
                        os.path.split(fn)[-1]), "w")
                    f.write(r.content)
                    f.close()
            else:
                # Copy the remote file to the local temp directory.
                remote_path = _remote_server + ":\""\
                        + " ".join([os.path.join(_remote_data_dir,
                            self._fullpath(*o)) for o in objs]) + "\""
                cmd = "scp %s %s" % (remote_path, _local_tmp_dir)
                proc = subprocess.Popen(cmd, shell=True,
                        stdout=subprocess.PIPE, close_fds=True)

                if proc.wait() != 0:
                    raise SDSSFileError("Couldn't copy %s." % (remote_path))


def get_filename(run, camcol, band):
    return os.path.join(_local_data_dir, "%d-%d-%s.hdf5" % (run, camcol, band))


def preprocess(run, camcol, fields, rerun, band, clobber=False):
    """
    Given a list of fields and specific (run, camcol, band), combine the
    fields into a single data stream and save it to an HDF5 file. In the
    overlapping regions, we just take the mean value of the two fields.
    This seems to work better than a weighted average but either should be
    fine because the data should be essentially _the same_ data in the
    overlap.

    ### Arguments

    * `run` (int): The run number.
    * `camcol` (int): The camera column number.
    * `rerun` (int): The data reduction pipeline rerun number.
    * `band` (str): One of "u", "g", "r", "i" or "z".

    ### Keyword Arguments

    * `clobber` (bool): Overwrite existing preprocessed data? (default: True)

    """
    band_id = "ugriz".index(band)

    # First, we have to check to make sure that the list of fields can be
    # coerced into a consecutive list. The way that we combine the fields
    # depends on this.
    fields.sort()
    if not (len(fields) == 1 or \
            np.all(np.array(fields) == np.arange(min(fields),
                                                 max(fields) + 1))):
        logging.warn("The fields must be consecutive. Skipping %d.%d." %
                (run, camcol))
        logging.warn(fields)
        return

    # Next, we check to see if a listing already exists in the database for
    # this particular (run, camcol, band).
    d = _run_collection.find_one({"run": run, "camcol": camcol, "band": band})

    # If there is, fail or overwrite it depending on the value of the
    # `clobber` option.
    if d is not None and clobber is False:
        logging.warn("An entry already exists for (run, camcol, band)"\
                + " = (%d, %d, %s). Skipping." % (run, camcol, band))
        return

    # Fetch all the needed FITS files from the server.
    files = [("tsField", run, camcol, f, rerun, band) for f in fields]
    files += [("fpC", run, camcol, f, rerun, band) for f in fields]
    files += [("fpM", run, camcol, f, rerun, band) for f in fields]
    files += [("psField", run, camcol, f, rerun, band) for f in fields]
    sdss = _DR7()
    try:
        sdss.fetch(files)
    except SDSSFileError:
        logging.warn("Couldn't fetch fields for %d.%d. Skipping." %
                (run, camcol))
        return

    # Allocate the output HDF5 file. Note: this will overwrite the existing
    # file if there is one.
    data = h5py.File(get_filename(run, camcol, band), "w")

    # What is the full shape of the imaging data?
    shape = [len(fields) * (_f_width - _f_overlap) + _f_overlap, _f_height]

    # Allocate space for image and inverse variance map/mask.
    data.create_dataset(IMG_TAG, shape, int, compression='gzip')
    data.create_dataset(INV_TAG, shape, np.float32, compression='gzip')

    # Get the astrometry.
    ast = []
    data.create_group(TS_TAG)
    for f in fields:
        try:
            ts = sdss.readTsField(run, camcol, f, rerun)
        except IOError:
            logging.warn("Couldn't read tsField for (%d, %d, %d)" %
                    (run, camcol, f))
            return
        g = data[TS_TAG].create_dataset(str(f), data=ts.hdus[1].data)
        hdu = ts.hdus[0]
        for k in hdu.header:
            g.attrs[k] = hdu.header[k]
        hdu = ts.hdus[1]
        for k in hdu.header:
            g.attrs[k] = hdu.header[k]
        ast.append(ts.getAsTrans(band))
        ts.hdus.close()

    # Get the PSF, etc.
    try:
        ps = [sdss.readPsField(run, camcol, f) for f in fields]
    except IOError:
        logging.warn("Couldn't read psField for (%d, %d)" % (run, camcol))
        return

    # Save the PSF information.
    data.create_group(PSF_TAG)
    for i, f in enumerate(fields):
        # Find the index of the PSF eigenimage HDU for this particular band.
        ind = ps[i].hdus[0].header["FILTERS"].split().index(band) + 1

        # Save the eigenimages & PSF.
        g = data[PSF_TAG].create_group(str(f))

        hdu = ps[i].hdus[ind]

        sg = g.create_group(EIGEN_TAG)
        for k in hdu.data.dtype.names:
            # This is a hack to deal with the fact that `h5py` can't handle
            # the `object` dtype. In this particular case, the object is
            # always a list of lists.
            if hdu.data[k].dtype == "object":
                o = [l for l in hdu.data[k]]
            else:
                o = hdu.data[k]
            sg.create_dataset(k, data=o)
        for k in hdu.header:
            sg.attrs[k] = hdu.header[k]

        # The general PSF information should be stored in HDU 6.
        hdu = ps[i].hdus[6]
        g.create_dataset(INFO_TAG, data=hdu.data)
        for k in hdu.header:
            g[INFO_TAG].attrs[k] = hdu.header[k]

    # Metadata.
    tai = []
    centers = np.zeros((len(fields), 2))
    bounds = np.zeros((len(fields), 4))  # ramin, ramax, decmin, decmax

    for ind, field in enumerate(fields):
        # Get the image and mask.
        try:
            fpC = sdss.readFpC(run, camcol, field, band)
            img = fpC.getImage()
        except IOError:
            logging.warn("Couldn't read fpC for (%d, %d, %d, %s)" %
                    (run, camcol, field, band))
            return

        try:
            fpM = sdss.readFpM(run, camcol, field, band)
            inv = sdss.getInvvar(img, fpM,
                             ps[ind].getGain(band_id),
                             ps[ind].getDarkVariance(band_id),
                             ps[ind].getSky(band_id),
                             ps[ind].getSkyErr(band_id))
        except IOError:
            logging.warn("Couldn't read fpM for (%d, %d, %d, %s)" %
                    (run, camcol, field, band))
            return

        # Position this field in the full image.
        dw = _f_width - _f_overlap
        if ind == 0:
            pos = np.arange(ind * dw, ind * dw + _f_width)
            data[IMG_TAG][pos, :] = img
            data[INV_TAG][pos, :] = inv
        else:
            pos = np.arange(ind * dw + _f_overlap, ind * dw + _f_width)
            overlap = ind * dw + np.arange(_f_overlap)

            # Fill in the non-overlapping region.
            data[IMG_TAG][pos, :] = img[_f_overlap:, :]
            data[INV_TAG][pos, :] = inv[_f_overlap:, :]

            # Average in the overlapping region.
            data[IMG_TAG][overlap, :] = 0.5 * (data[IMG_TAG][overlap, :] \
                            + img[:_f_overlap, :])
            data[INV_TAG][overlap, :] = 0.5 * (data[INV_TAG][overlap, :]\
                            + inv[:_f_overlap, :])

        # Get the time of the observation.
        tai.append(fpC.hdus[0].header["TAI"])

        # Get the approximate astrometric bounds of field.
        c = np.array([ast[ind].pixel_to_radec(0,         0),
                      ast[ind].pixel_to_radec(0,         _f_width),
                      ast[ind].pixel_to_radec(_f_height, _f_width),
                      ast[ind].pixel_to_radec(_f_height, 0)])
        bounds[ind] = [min(c[:, 0]), max(c[:, 0]), min(c[:, 1]), max(c[:, 1])]
        centers[ind] = ast[ind].pixel_to_radec(0.5 * _f_height, 0.5 * _f_width)

        fpC.hdus.close()
        fpM.hdus.close()

    data.create_dataset(TAI_TAG, data=tai)
    data.create_dataset(BOUNDS_TAG, data=bounds)
    data.create_dataset(CENTERS_TAG, data=centers)

    [f.hdus.close() for f in ps]

    data.close()

    # Write to the database.
    d = _run_collection.find_one({"run": run, "camcol": camcol, "band": band})
    if d is None:
        d = {"run": run, "camcol": camcol, "band": band}
    d["fields"] = fields
    d["raMin"], d["raMax"] = np.min(bounds[:, 0]), np.max(bounds[:, 1])
    d["decMin"], d["decMax"] = np.min(bounds[:, 2]), np.max(bounds[:, 3])
    _run_collection.save(d)


def _pp_wrapper(r):
    preprocess(*r)


def preprocess_multiple(runs, pool=None):
    if pool is None:
        pool = multiprocessing.Pool(8)
    pool.map(_pp_wrapper, runs)


def cleanup():
    shutil.rmtree(_local_tmp_dir)


def main(argv):
    import argparse
    import pymongo

    # Parse input files.
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bands",
        help="The SDSS bands to analyze.",
        type=str,
        default="ugriz")
    parser.add_argument("-r", "--ra",
        help="The range of R.A.s to analyze.",
        nargs=2,
        type=float)
    args = parser.parse_args(argv)

    # Connect to database.
    collection = _db.fields
    collection.ensure_index("raMin")
    collection.ensure_index("raMax")
    collection.ensure_index([("run",    pymongo.ASCENDING),
                             ("camcol", pymongo.ASCENDING)])

    # Find the fields.
    rng = sorted(args.ra)
    q = {"raMin": {"$lt": rng[1]}, "raMax": {"$gt": rng[0]}, "rerun": 40}
    f = {"_id": 0, "run": 1, "camcol": 1, "field": 1, "rerun": 1}

    logging.info("Finding fields in specified range.")
    docs = [d for d in collection.find(q, f)
                        .sort([("run", 1), ("camcol", 1)])]

    field_list = []
    for i, d in enumerate(docs):
        if i > 0:
            o = docs[i - 1]
        if i == 0 or o["run"] != d["run"] or o["camcol"] != d["camcol"]:
            field_list.append([d["run"], d["camcol"], [d["field"], ],
                d["rerun"]])
        else:
            field_list[-1][2].append(d["field"])

    for b in args.bands:
        runs = [f[:2] + [range(min(f[2]), max(f[2]) + 1)] + [f[3]] + [b]
                for f in field_list]
        logging.info("Processing %d %s-band runs." % (len(runs), b))
        preprocess_multiple(runs)


if __name__ == "__main__":
    import sys
    main(sys.argv)
