#!/usr/bin/env python
# encoding: utf-8
"""
Fit the calibration model

History
-------
2011-06-16 - Created by Dan Foreman-Mackey

"""

__all__ = ['init_model','calibrate_grid','calibrate','do_calibration']

import datetime

import numpy as np
import numpy.ma as ma
import scipy.optimize as op

import database
from photomodel import *
#from opt import survey
from photometry import *
from opt import *

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
    tmp = np.mean(data.flux/mag2nmgy(data.magprior[:,0]), axis=-1)
    p0 = []
    for i in range(data.nobs):
        p0.append(tmp[i])
    p0 = ma.concatenate([p0,data.magprior[:,0]])

    zero = p0[:data.nobs]
    flux = mag2nmgy(p0[data.nobs:])

    # jitterabs, jitterrel
    #p0 = np.append(p0,0.01)

    # ln(sigvar2)
    #ff = np.outer(zero,flux)
    #sigvar2 = np.max(np.var(data.flux,axis=0))

    # lnoddsvar, Q
    #p0 = np.append(p0,[prob2logodds(0.03),1.0])#np.log(sigvar2))
    #p0 = np.append(p0,1.0)#np.log(sigvar2))
    #p0 = np.append(p0,[prob2logodds(0.03)])#np.log(sigvar2))

    # pbad, ln(sigmabad2)
    #sigbad2 = np.median(np.var(data.flux,axis=0))
    #p0 = np.append(p0,1e4)#,np.log(sigbad2)])

    return p0

def calibrate_grid(coords,radius,meta=None):
    """
    Fit calibration model on a grid of RA/Dec points

    Parameters
    ----------
    coords : numpy.ndarray (shape : [npoints,2])
        A list of coordinates (ra,dec) in degrees

    radius : float
        Selection radius (passed to survey.find_stars) in arcmin

    Optional
    --------
    meta : dict (default : None)
        Meta data to include in model database

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
        points.append(calibrate(pos[0],pos[1],radius,meta=meta))
    return points

def calibrate(ra,dec,radius,meta=None):
    """
    Optimize calibration model using stars found in annulus at RA/Dec

    Parameters
    ----------
    ra : float
        In degrees

    dec : float
        In degrees

    radius : float
        Selection radius (passed to survey.find_stars) in arcmin

    Returns
    -------
    photo_model : model.PhotoModel
        The optimized model object

    History
    -------
    2011-06-22 - Created by Dan Foreman-Mackey

    """
    print datetime.datetime.today()
    print "calibrate:",ra,dec,radius,meta
    kwargs = {}
    if meta is not None:
        if 'mgroup' in meta:
            kwargs['mgroup'] = meta['mgroup']
        if 'resample' in meta:
            kwargs['resample'] = meta['resample']
    obs,stars = find_photometry(ra,dec,radius,**kwargs)
    print "calibrate: found %d stars in %d observations"%\
            (len(stars),len(obs))
    if len(stars) <= 1 or len(obs) <= 1:
        print "calibrate: Couldn't find any measurements!"
        return None
    data = get_photometry(obs,stars)
    photo_data = PhotoData(data,obs,stars)

    return do_calibration(photo_data,addtodb=True,meta=meta,ra=ra,dec=dec,radius=radius)

def do_calibration(photo_data,addtodb=False,meta=None,ra=0.0,dec=0.0,radius=0.0,
        p0=None,fix_probs=None):
    if p0 is None:
        p0 = init_model(photo_data)
    model0 = PhotoModel(photo_data,p0)

    chi2 = lambda p: -lnprob(p,photo_data,fix_probs=fix_probs)
    print "calibrate: optimizing model"
    #p1 = op.fmin_l_bfgs_b(chi2,p0,approx_grad=True)[0]
    import time
    strt = time.time()
    p1 = op.fmin_bfgs(chi2,p0)
    print "optimize took: %f seconds"%(time.time()-strt)
    photo_model = PhotoModel(photo_data,p1)

    if addtodb:
        doc = {'pos': {'ra':ra,'dec':dec},'radius':radius,
                'model': photo_model}
        if meta is not None: # append meta data
            for k in list(meta):
                if k not in doc and not k == '_id':
                    doc[k] = meta[k]
        modelid = database.photomodel.insert(doc)
        for oi in range(photo_data.nobs):
            doc = {'pos': {'ra':ra,'dec':dec},'radius':radius,
                    'modelid':modelid,
                    'zero':photo_model.zero[oi],
                    'measurements': []}
            #for i in np.arange(photo_data.nobs)\
            #        [np.array(photo_data.obsorder) == oi]:
            #    doc['measurements'].append(photo_data.observations[i])
            if meta is not None: # append meta data
                for k in list(meta):
                    if k not in doc and not k == '_id':
                        doc[k] = meta[k]
            database.obslist.insert(doc)

    return photo_model

if __name__ == '__main__':
    database.photomodel.drop()
    database.obslist.drop()

    # params
    grid_spacing = 5.0 # in arcmin
    delta = grid_spacing/60.0
    radius = 5.

    # run the grid
    #for grid_spacing in [60.,30.,20.,10.,5.]:
    for mgroup in [None]: #range(10):
        for dec in np.arange(-1.25+delta,0.75,delta):
            for ra in np.arange(20.0+delta,22.0,delta):
                for radius in [1.5,5.,15.]:
                    calibrate(ra,dec,radius,
                            meta={'grid_spacing': grid_spacing,
                                'mgroup': mgroup})

    # make sure that we have all the indexes set up
    import pymongo
    database.photomodel.ensure_index([('pos',pymongo.GEO2D)])
    database.photomodel.ensure_index('radius')
    database.photomodel.ensure_index('grid_spacing')

    database.obslist.ensure_index([('pos',pymongo.GEO2D)])
    database.obslist.ensure_index('radius')
    database.obslist.ensure_index('grid_spacing')
    database.obslist.ensure_index('obsid')
    database.obslist.ensure_index('modelid')
    database.obslist.ensure_index('measurements._id')
    database.obslist.ensure_index('measurements.mjd_g')
    database.obslist.ensure_index([('measurements.run',pymongo.ASCENDING),
        ('measurements.camcol',pymongo.ASCENDING),
        ('measurements.field',pymongo.ASCENDING)])



