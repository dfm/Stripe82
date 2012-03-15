#!/usr/bin/env python
# encoding: utf-8
"""
Utilities for preprocessing a batch of fields and combining them into one.

"""

__all__ = [""]

import logging
import os
import shutil
import subprocess

import numpy as np

import pyfits
import pysdss

_remote_data_dir = os.environ["SDSS_REMOTE"]
_local_tmp_dir   = os.path.join(os.environ.get("SDSS_LOCAL", "."), ".sdss")

class _DR7(pysdss.DR7):
    def __init__(self, *args, **kwargs):
        super(_DR7, self).__init__(*args, **kwargs)
        self.filenames = {
            "fpObjc":  'fpObjc-%(run)06i-%(camcol)i-%(field)04i.fit',
            "fpM":     'fpM-%(run)06i-%(band)s%(camcol)i-%(field)04i.fit',
            "fpC":     'fpC-%(run)06i-%(band)s%(camcol)i-%(field)04i.fit',
            "fpAtlas": 'fpAtlas-%(run)06i-%(camcol)i-%(field)04i.fit',
            "psField": 'psField-%(run)06i-%(camcol)i-%(field)04i.fit',
            "tsObj":   'tsObj-%(run)06i-%(camcol)i-%(rerun)i-%(field)04i.fit',
            "tsField": \
                    'tsField-%(run)06i-%(camcol)i-%(rerun)i-%(field)04i.fit',
            }

    def getFilename(self, filetype, *args, **kwargs):
        kwargs = dict(zip(["run", "camcol", "field", "band"], args))
        kwargs["rerun"] = kwargs.pop("rerun", 40)
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
        # Set up the local file system to accept the file.
        local_path  = os.path.join(_local_tmp_dir, fn)
        local_dir = os.path.join(*os.path.split(local_path)[:-1])
        try:
            os.makedirs(local_dir)
        except os.error as e:
            if not os.path.exists(local_dir):
                raise e

        # Copy the remote file to the local temp directory.
        remote_path = os.path.join(_remote_data_dir, fn)
        logging.info("Running: scp %s %s"%(remote_path, local_path))
        proc = subprocess.Popen("scp %s %s"%(remote_path, local_path),
                    shell=True, stdout=subprocess.PIPE, close_fds=True)

        if proc.wait() != 0:
            logging.warn("Couldn't copy %s."%(remote_path))
            return None
        return pyfits.open(local_path)

def preprocess(run, camcol, fields, rerun, band):
    sdss = _DR7()
    tsField = sdss.readTsField(run, camcol, fields[0], rerun)

def cleanup():
    shutil.rmtree(_local_tmp_dir)

if __name__ == '__main__':
    import argparse
    logging.basicConfig(level=logging.INFO)

    preprocess(4263, 4, [84], 40, "g")

