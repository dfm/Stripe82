#!/usr/bin/env python
# encoding: utf-8
"""


History
-------
2011-07-13 - Created by Dan Foreman-Mackey

"""

import time

import numpy as np
import h5py

import calibration as calib

calib.database.photomodel.drop()
calib.database.obslist.drop()

np.random.seed()
N = 10
rs = [1.5,5,10,15]
outfile = h5py.File('timings.hdf5','w')
outfile['rs'] = rs
outfile['coords'] = np.zeros([len(rs),N,2])
outfile['data'] = np.zeros([len(rs),N,3])
for ri,radius in enumerate(rs):
    for i in range(N):
        res = None
        while res is None:
            ra = 20 + 2*np.random.rand()
            dec = -0.25 + 2*np.random.rand()
            strt = time.time()
            res = calib.calibrate(ra,dec,radius)
            strt = time.time()-strt
            if res is not None:
                outfile['coords'][ri,i,:] = [ra,dec]
                outfile['data'][ri,i,:] = [strt,res.data.nstars,res.data.nobs]

outfile.close()

