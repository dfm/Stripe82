#!/usr/bin/env python
# encoding: utf-8
"""
Extract the calibrated light-curve at a given position

History
-------
2011-06-22 - Created by Dan Foreman-Mackey

"""

__all__ = ['extract_calibrated_lightcurve']

import numpy as np

from calibrate import calibrate
from photometry import force_photometry

def extract_calibrated_lightcurve(ra,dec,radius=3.0):
    """
    Extract the lightcurve at a given RA/Dec
    
    Parameters
    ----------
    ra : float
        In degrees

    dec : float
        In degrees

    Optional
    --------
    radius : float (default : 3)
        In arcmin. Search radius for finding calibration stars.
    
    Returns
    -------
    time : numpy.ndarray (N,)
        Times for each data point (in MJD)

    maggies : numpy.ndarray (N,2)
        (Maggies, error) at each time
    
    History
    -------
    2011-06-22 - Created by Dan Foreman-Mackey
    
    """
    model = calibrate(ra,dec,radius=radius)
    lc = []
    for oi,obs in enumerate(model.data.observations):
        df,dferr = tuple(force_photometry(ra,dec,obs['_id']))
        df /= model.zero[oi]
        dferr /= model.zero[oi]
        lc.append([df,dferr])
    return model.data.mjd(),np.array(lc)

if __name__ == '__main__':
    time,maggies = extract_calibrated_lightcurve(-23.431965,-0.227934)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as pl

    pl.errorbar(time,maggies[:,0],yerr=maggies[:,1],fmt='.k')
    pl.savefig('extracted.pdf')


