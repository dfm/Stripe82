#!/usr/bin/env python
# encoding: utf-8
"""
Command line utility for extracting a light-curve from Stripe 82 at a given RA/Dec

History
-------
2011-08-03 - Created by Dan Foreman-Mackey

"""

__all__ = ['extract_lightcurve']

import pylab as pl
import numpy as np
from calibration import calibrate

def extract_lightcurve(ra,dec,radius):
    model = calibrate(ra,dec,radius=radius)
    mjd = model.data.mjd()
    flux = model.data.flux/model.zero[:,np.newaxis]
    for i in range(np.shape(flux)[1]):
        pl.clf()
        pl.plot(mjd,flux[:,i],'+k')
        pl.savefig('lcs/%03d.png'%i)

if __name__ == '__main__':
    extract_lightcurve(29.47942,0.383557,10)


