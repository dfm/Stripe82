#!/usr/bin/env python
# encoding: utf-8
"""
Do batch photometry for a set of stars in a set of fields

History
-------
2011-06-14 - Created by Dan Foreman-Mackey

"""

__all__ = ['force_photometry','do_photometry','get_photometry','find_photometry']

import time as timer
import datetime

import numpy as np

import pymongo

# options
import opt
from opt import survey

# databases
import database

def force_photometry(ra,dec,observation):
    """
    Force the photometry at a given R.A./Dec. in a given observation
    
    Parameters
    ----------
    ra : float
        In degrees

    dec : float
        In degrees

    observation : bson.ObjectID
        Pointer to specific observation
    
    Returns
    -------
    res : list (len : 2)
        (counts, error)

    Notes
    -----
    Returns [0,0] if nothing is found.

    Warnings
    --------
    If you change this, change do_photometry too!
    
    History
    -------
    2011-06-22 - Created by Dan Foreman-Mackey
    
    """
    info = survey.get_observation(observation)
    if info is None:
        print 'info is None'
        return [0,0]
    if info['ramin'] < ra < info['ramax'] and \
            info['decmin'] < dec < info['decmax']:
        try:
            obs = survey.Observation(observation)
        except survey.ObservationAccessError:
            print "couldn't access data"
            return [0,0]

        try:
            res = obs.photometry(ra,dec)
            return res
        except survey.PhotometryError:
            print "couldn't measure photometry"
            return [0,0]
    print "out of bounds"
    return [0,0]
    

def do_photometry():
    """
    Do the photometry for all of the stars in all of the observations
    
    Parameters
    ----------
    observations : list
        List of bson.ObjectID objects for observations

    stars : list
        List of bson.ObjectID objects for stars

    Warnings
    --------
    If you change this, change force_photometry too!

    History
    -------
    2011-06-14 - Created by Dan Foreman-Mackey
    
    """
    np.random.seed()
    database.photoraw.drop()
    nstarsperobs = []
    timeperobs = []
    observations = survey.find_all_observations()
    print len(observations)
    for oi,obsid in enumerate(observations):
        strt = timer.time()
        info = survey.get_observation(obsid)
        try:
            obs  = survey.Observation(obsid)
        except survey.ObservationAccessError:
            print "Couldn't access data"
            obs = None
        if obs is not None:
            stars = survey.find_stars_in_observation(obsid)
            nstarsperobs.append(len(stars))
            for starid in stars:
                star = survey.get_star(starid)
                try:
                    res,cov = obs.photometry(star['ra'],star['dec'])
                    doc = {'obsid':obsid,'starid':starid,
                        'model':res,'cov':cov,
                        'pos':{'ra':star['pos']['ra'],'dec':star['pos']['dec']},
                        'mgroup': np.random.randint(opt.nmgroups)}
                    for k in list(star):
                        if k not in doc and not k == '_id':
                            doc[k] = star[k]
                    for k in list(info):
                        if k not in doc and not k == '_id':
                            doc[k] = info[k]
                    database.photoraw.insert(doc)
                except survey.PhotometryError:
                    # couldn't do photometry
                    pass
        timeperobs.append(timer.time()-strt)
        dt = datetime.timedelta(seconds=(len(observations)-oi-1)*np.mean(timeperobs))
        print "do_photometry: %d/%d, approx. %.0f stars/obs and %s remaining"\
                %(oi+1,len(observations),
                np.mean(nstarsperobs),dt)

    # create the indexes
    print "do_photometry: generating indexes"
    database.photoraw.create_index([('pos',pymongo.GEO2D)])
    database.photoraw.create_index('obsid')
    database.photoraw.create_index('starid')
    database.photoraw.create_index('mjd_g')
    database.photoraw.create_index([('run',1), ('camcol',1), ('field',1)])

def find_photometry(ra,dec,radius):
    """
    Find all of the photometric measurements in radius around (ra,dec)
    
    Parameters
    ----------
    ra : float
        In degrees

    dec : float
        In degrees

    radius : float
        Search radius in arcmin
    
    Returns
    -------
    obsids : list
        List of bson.ObjectID objects for observations
    
    stars : list
        List of bson.ObjectID objects for stars

    History
    -------
    2011-07-04 - Created by Dan Foreman-Mackey
    
    """
    radius /= 60. # to degrees
    radius = np.radians(radius)

    res = database.photoraw.find({'pos':{'$within':{'$centerSphere': [[ra,dec],radius]}}})
    
    obsids = []
    stars  = []

    for doc in res:
        obsids.append(res['obsid'])
        stars.append(res['starid'])

    return obsids,stars

def get_photometry(observations,stars):
    """
    Get the photometry for given stars in given observations
    
    Parameters
    ----------
    observations : list
        List of bson.ObjectID objects for observations

    stars : list
        List of bson.ObjectID objects for stars
    
    Returns
    -------
    data : numpy.ndarray (shape: [len(observations),len(stars),2]
        The observation matrix the final column has the form (counts,error)

    History
    -------
    2011-07-01 - Created by Dan Foreman-Mackey
    
    """
    data = np.zeros([len(observations),len(stars),2])
    for oi,obsid in enumerate(observations):
        for si,starid in enumerate(stars):
            entry = database.photoraw.find_one({'obsid': obsid, 'starid': starid})
            if entry is not None:
                data[oi,si,:] = [entry['model'][1],entry['covar'][1][1]]
    return data

if __name__ == '__main__':
    # do_photometry()
    
    obs,stars = find_photometry(21,0,5)
    print get_photometry(obs,stars)

