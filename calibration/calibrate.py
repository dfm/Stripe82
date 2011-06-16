#!/usr/bin/env python
# encoding: utf-8
"""
Fit the calibration model

History
-------
2011-06-16 - Created by Dan Foreman-Mackey

"""

__all__ = ['init_model','calibrate_grid']

import numpy as np
import numpy.ma as ma
import scipy.optimize as op

import database
from model import *
from opt import survey
from photometry import photometry

def init_model(data):
    """
    Initialize the model parameters to a reasonable guess given the data
    
    Parameters
    ----------
    data : PhotoData
        The data
    
    Returns
    -------
    p0 : numpy.ndarray
        A vector of model parameters
    
    History
    -------
    2011-06-16 - Created by Dan Foreman-Mackey
    
    """
    # we'll use the magnitude values from CAS
    # NOTE: I'm using masked array operations here because some of the 
    # measured fluxes will be be <= 0
    p0 = ma.mean(-2.5*ma.log10(data.flux)-data.magprior[:,0], axis=-1)
    p0 = ma.concatenate([p0,data.magprior[:,0]])

    zero = 10**(-p0[:data.nobs]/2.5)
    flux = 10**(-p0[data.nobs:]/2.5)

    # ln(jitterabs2), ln(jitterrel2), pvar
    p0 = np.append(p0,[-10.0,-10.0,0.5])

    # ln(sigvar2)
    ff = np.outer(zero,flux)
    delta = data.flux - ff
    sigvar2 = np.max(np.var(delta,axis=0))
    p0 = np.append(p0,np.log(sigvar2))
    # pbad, ln(sigmabad2)
    p0 = np.append(p0,[0.5,np.log(sigvar2)])

    return p0

def calibrate_grid(coords,radius=3.):
    """
    Fit calibration model on a grid of RA/Dec points
    
    Parameters
    ----------
    coords : numpy.ndarray (shape : [npoints,2])
        A list of coordinates (ra,dec) in degrees

    Optional
    --------
    radius : float (default : 3.0)
        Selection radius (passed to survey.find_stars) in arcmin
    
    Returns
    -------
    points : list
        A list of PhotoModel objects
    
    History
    -------
    2011-06-16 - Created by Dan Foreman-Mackey
    
    """
    points = []
    for pos in coords:
        print "calibrate_grid:",pos
        obs = survey.find_observations(pos[0],pos[1])
        stars = survey.find_stars(pos[0],pos[1],radius=radius)
        print "calibrate_grid: found %d stars in %d observations"%\
                (len(stars),len(obs))
        data = photometry(obs,stars)
        photo_data = PhotoData(data,obs,stars)
        p0 = init_model(photo_data)

        chi2 = lambda p: -lnprob(p,photo_data)
        print "calibrate_grid: optimizing model"
        p1 = op.fmin(chi2,p0)

        photo_model = PhotoModel(photo_data,p0)
        points.append(photo_model)

    return points

if __name__ == '__main__':
    print calibrate_grid([[-23.431965,-0.227934]])


