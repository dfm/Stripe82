#!/usr/bin/env python
# encoding: utf-8

import timeit

import numpy as np
np.random.seed()

import markovpy

import calibration as calib

def crossvalidation(ra,dec):
    """
    Perform cross validation in a patch to test calibration model stability
    
    Parameters
    ----------
    ra : float
        In degrees

    dec : float
        In degrees
    
    Returns
    -------
    ret : type
        Description
    
    History
    -------
    2011-07-20 - Created by Dan Foreman-Mackey
    
    """
    models = []
    rs = [1.5,3,5]
    for radius in rs:
        models.append([])
        for mgroup in range(5):
            p = calib.calibrate(ra,dec,radius,
                meta={'mgroup': mgroup,
                      'cross_validation': 1})
            models[-1].append(p)

    likelihood = []
    for i0,row0 in enumerate(models):
        print 'Radius: ',rs[i0]
        for j0,p0 in enumerate(row0):
            obs,stars = calib.find_photometry(ra,dec,rs[i0],
                    mgroup={'$ne': j0})
            data = calib.get_photometry(obs,stars)
            photo_data = calib.PhotoData(data,obs,stars)
            
            like = calib.lnlikelihood(p0,photo_data)
            likelihood.append(like)
            print '\t',j0,like

def sample_posterior(ra,dec,radius):
    """
    Sample the posterior likelihood function
    
    Parameters
    ----------
    ra : float
        In degrees

    dec : float
        In degrees

    radius : float
        In arcmin
    
    History
    -------
    2011-07-21 - Created by Dan Foreman-Mackey
    
    """
    obs,stars = calib.find_photometry(ra,dec,radius)
    print "calibrate: found %d stars in %d observations"%\
            (len(stars),len(obs))
    data = calib.get_photometry(obs,stars)
    photo_data = calib.PhotoData(data,obs,stars)
    p0 = calib.init_model(photo_data)

    ndim     = len(p0)
    nwalkers = 2*ndim
    
    init_pos = [p0+0.1*np.random.randn(ndim) for i in range(nwalkers)]
    lnprob = lambda x: calib.lnprob(x,photo_data)
    sampler = markovpy.EnsembleSampler(nwalkers,ndim,lnprob)
    print "sampler: starting"
    pos,prob,state = sampler.run_mcmc(init_pos, None, 100)
    print "sampler: burnin finished"
    sampler.clear_chain()
    sampler.run_mcmc(pos, state, 500)

    return sampler
    

if __name__ == '__main__':
    crossvalidation(21,0)

