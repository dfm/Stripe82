#!/usr/bin/env python
# encoding: utf-8
"""
Utilities for preprocessing a batch of fields and combining them into one.

"""

__all__ = ["SDSSFileError", "preprocess", "cleanup"]

import logging
import os
import shutil
import subprocess

import numpy as np

import pyfits
import pysdss

import h5py

# Data access variables.
_remote_server   = os.environ["SDSS_SERVER"]
_remote_data_dir = os.environ["SDSS_REMOTE"]
_local_tmp_dir   = os.path.join(os.environ.get("SDSS_LOCAL", "."), ".sdss")
_local_data_dir  = os.path.join(os.environ.get("SDSS_LOCAL", "."), "data")

# Field dimensions.
_f_width   = 1489
_f_height  = 2048
_f_overlap = 128

# HDF5 file tags.
IMG_TAG         = "img"
INV_TAG         = "inv"
TAI_TAG         = "tai"
BOUNDS_TAG      = "bounds"
CENTERS_TAG     = "centers"
PSF_TAG         = "psf"
EIGEN_TAG       = "eigen"
INFO_TAG        = "info"

class SDSSFileError(Exception):
    pass

class _DR7(pysdss.DR7):
    def __init__(self, *args, **kwargs):
        super(_DR7, self).__init__(*args, **kwargs)
        self.filenames = {
            "fpObjc":  'fpObjc-%(run)06i-%(camcol)i-%(field)04i.fit',
            "fpM":     'fpM-%(run)06i-%(band)s%(camcol)i-%(field)04i.fit.gz',
            "fpC":     'fpC-%(run)06i-%(band)s%(camcol)i-%(field)04i.fit.gz',
            "fpAtlas": 'fpAtlas-%(run)06i-%(camcol)i-%(field)04i.fit',
            "psField": 'psField-%(run)06i-%(camcol)i-%(field)04i.fit',
            "tsObj":   'tsObj-%(run)06i-%(camcol)i-%(rerun)i-%(field)04i.fit',
            "tsField": \
                    'tsField-%(run)06i-%(camcol)i-%(rerun)i-%(field)04i.fit',
            }

    def _fullpath(self, filetype, *args, **kwargs):
        for k,v in zip(["run", "camcol", "field", "rerun", "band"], args):
            kwargs[k] = v
        fn = self.filenames[filetype] % kwargs

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
            # Copy the remote file to the local temp directory.
            remote_path = _remote_server+":\""\
                    +" ".join([os.path.join(_remote_data_dir, self._fullpath(*o))
                                for o in objs])+"\""
            logging.info("Running: scp %s %s"%(remote_path, _local_tmp_dir))
            proc = subprocess.Popen("scp %s %s"%(remote_path, _local_tmp_dir),
                        shell=True, stdout=subprocess.PIPE, close_fds=True)

            if proc.wait() != 0:
                logging.warn("Couldn't copy %s."%(remote_path))
                raise SDSSFileError()

def get_filename(run, camcol, band):
    return os.path.join(_local_data_dir, "%d-%d-%s.hdf5"%(run, camcol, band))

def preprocess(run, camcol, fields, rerun, band):
    # Make sure that the output directory exists.
    try:
        os.makedirs(_local_data_dir)
    except os.error as e:
        if not os.path.exists(_local_data_dir):
            raise e

    band_id = "ugriz".index(band)

    # Check to make sure that the fields are consecutive.
    fields.sort()
    assert len(fields) == 1 or \
            np.all(np.array(fields)==np.arange(min(fields), max(fields)+1)),\
                        "The fields must be consecutive."

    # Fetch all the needed files in one pass.
    files =  [("tsField", run, camcol, f, rerun, band) for f in fields]
    files += [("fpC", run, camcol, f, rerun, band) for f in fields]
    files += [("fpM", run, camcol, f, rerun, band) for f in fields]
    files += [("psField", run, camcol, f, rerun, band) for f in fields]
    sdss = _DR7()
    sdss.fetch(files)

    # Set up the output file.
    data = h5py.File(get_filename(run, camcol, band), "w")

    # What is the full shape of the imaging data?
    shape = [len(fields)*(_f_width-_f_overlap) + _f_overlap, _f_height]

    # Allocate space for image and inverse variance map/mask.
    data.create_dataset(IMG_TAG, shape, int, compression='gzip')
    data.create_dataset(INV_TAG, shape, np.float32, compression='gzip')

    # Get the astrometry.
    tsFields = [sdss.readTsField(run, camcol, f, rerun) for f in fields]
    ast      = [t.getAsTrans(band) for t in tsFields]

    # Get the PSF, etc.
    ps = [sdss.readPsField(run, camcol, f) for f in fields]

    # Save the PSF information.
    data.create_group(PSF_TAG)
    for i, f in enumerate(fields):
        # Find the index of the PSF eigenimage HDU for this particular band.
        ind = ps[i].hdus[0].header["FILTERS"].split().index(band) + 1

        # Save the eigenimages & PSF.
        g = data[PSF_TAG].create_group(str(f))

        hdu = ps[i].hdus[ind]
        g.create_dataset(EIGEN_TAG, data=hdu.data)
        for k in hdu.header:
            g[EIGEN_TAG].attrs[k] = hdu.header[k]

        # The general PSF information should be stored in HDU 6.
        hdu = ps[i].hdus[6]
        g.create_dataset(INFO_TAG, data=hdu.data)
        for k in hdu.header:
            g[INFO_TAG].attrs[k] = hdu.header[k]

    # Metadata.
    tai = []
    centers = np.zeros((len(fields), 2))
    bounds = np.zeros((len(fields), 4)) # ramin, ramax, decmin, decmax

    for ind, field in enumerate(fields):
        # Get the image and mask.
        fpC = sdss.readFpC(run, camcol, field, band)
        img = fpC.getImage()

        fpM = sdss.readFpM(run, camcol, field, band)
        inv = sdss.getInvvar(img, fpM,
                             ps[ind].getGain(band_id),
                             ps[ind].getDarkVariance(band_id),
                             ps[ind].getSky(band_id),
                             ps[ind].getSkyErr(band_id))

        # Position this field in the full image.
        dw = _f_width - _f_overlap
        if ind == 0:
            pos = np.arange(ind*dw, ind*dw+_f_width)
            data[IMG_TAG][pos,:] = img
            data[INV_TAG][pos,:] = inv
        else:
            pos     = np.arange(ind*dw+_f_overlap, ind*dw+_f_width)
            overlap = ind*dw + np.arange(_f_overlap)

            # Fill in the non-overlapping region.
            data[IMG_TAG][pos,:] = img[_f_overlap:,:]
            data[INV_TAG][pos,:] = inv[_f_overlap:,:]

            # Average in the overlapping region.
            data[IMG_TAG][overlap,:] = 0.5 * (data[IMG_TAG][overlap,:]\
                            + img[:_f_overlap,:])
            data[INV_TAG][overlap,:] = 0.5 * (data[INV_TAG][overlap,:]\
                            + inv[:_f_overlap,:])

        # Get the time of the observation.
        tai.append(fpC.hdus[0].header["TAI"])

        # Get the approximate astrometric bounds of field.
        c = np.array([ast[ind].pixel_to_radec(0,         0),
                      ast[ind].pixel_to_radec(0,         _f_width),
                      ast[ind].pixel_to_radec(_f_height, _f_width),
                      ast[ind].pixel_to_radec(_f_height, 0)])
        bounds[ind]  = [min(c[:,0]), max(c[:,0]), min(c[:,1]), max(c[:,1])]
        centers[ind] = ast[ind].pixel_to_radec(0.5*_f_height, 0.5*_f_width)

        fpC.hdus.close()
        fpM.hdus.close()

    data.create_dataset(TAI_TAG, data=tai)
    data.create_dataset(BOUNDS_TAG, data=bounds)
    data.create_dataset(CENTERS_TAG, data=centers)

    [f.hdus.close() for f in tsFields]
    [f.hdus.close() for f in ps]

    data.close()

def cleanup():
    shutil.rmtree(_local_tmp_dir)

if __name__ == '__main__':
    import argparse
    logging.basicConfig(level=logging.INFO)

    preprocess(4263, 4, [83, 84, 85, 86], 40, "g")

