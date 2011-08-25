#!/usr/bin/env python
# encoding: utf-8
"""
Extract the calibrated light-curve at a given position

History
-------
2011-06-22 - Created by Dan Foreman-Mackey

"""

__all__ = ['extract_lightcurves']

import numpy as np
import numpy.ma as ma

from calibrate import calibrate

def extract_lightcurves(*args,**kwargs):
    if 'units' in kwargs:
        units = kwargs['units']
    else:
        units = 1.0
    if 'model' in kwargs:
        model = kwargs['model']
    else:
        (ra,dec,radius) = args
        kws = {'radius': radius}
        if 'band' in kwargs:
            kws['band'] = kwargs['band']
        model = calibrate(ra,dec,**kws)
    mjd = model.data.mjd()
    flux = model.data.flux/model.zero[:,np.newaxis]*units
    mask = model.data.ivar <= 1e-8
    inv = ma.array(model.data.ivar,mask=mask)
    err = np.sqrt(1.0/inv+model.jitterabs2+\
            model.jitterrel2*np.outer(model.zero,model.flux)**2)/\
            (model.zero[:,np.newaxis])
    mjd_inds = mjd > 1
    return mjd[mjd_inds],flux[mjd_inds],err[mjd_inds],model

