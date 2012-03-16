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

_remote_server   = os.environ["SDSS_SERVER"]
_remote_data_dir = os.environ["SDSS_REMOTE"]
_local_tmp_dir   = os.path.join(os.environ.get("SDSS_LOCAL", "."), ".sdss")

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

def preprocess(run, camcol, fields, rerun, band):
    # Check to make sure that the fields are consecutive.
    fields.sort()
    assert len(fields) == 1 or \
            np.all(np.array(fields)==np.arange(min(fields), max(fields)+1)),\
                        "The fields must be consecutive."

    # Fetch all the needed files in one pass.
    files =  [("tsField", run, camcol, f, rerun, band) for f in fields]
    files += [("fpC", run, camcol, f, rerun, band) for f in fields]
    files += [("fpM", run, camcol, f, rerun, band) for f in fields]
    sdss = _DR7()
    sdss.fetch(files)

    # Get the astrometry.
    tsFields = [sdss.readTsField(run, camcol, f, rerun) for f in fields]
    ast      = [t.getAsTrans(band) for t in tsFields]

def cleanup():
    shutil.rmtree(_local_tmp_dir)

if __name__ == '__main__':
    import argparse
    logging.basicConfig(level=logging.INFO)

    preprocess(4263, 4, [83, 85], 40, "g")

