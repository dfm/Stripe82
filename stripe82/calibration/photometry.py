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
import multiprocessing

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

_meta = {'band': 'g'}

def photometer_run(obsid):
    band = _meta['band']
    info = survey.get_observation(obsid)
    try:
        obs  = survey.Observation(obsid,band=band)
    except survey.ObservationAccessError:
        print "Couldn't access data"
        obs = None
    if obs is not None:
        doc = {}
        for k in list(info):
            if k not in doc and not k == '_id':
                doc[k] = info[k]
        for k in list(_meta):
            doc[k] = _meta[k]

        stars = survey.find_stars_in_observation(obsid)
        for starid in stars:
            star = survey.get_star(starid)
            for k in list(star):
                if not k == '_id':
                    doc[k] = star[k]
            try:
                res,cov = obs.photometry(star['ra'],star['dec'])
                doc['obsid']  = obsid
                doc['starid'] = starid
                doc['model']  = res
                doc['cov']    = cov
                doc['pos']    = [star['ra'],star['dec']]
                doc['mgroup'] = np.random.randint(opt.nmgroups)
                doc['mjd']    = obs.mjd(star['ra'],star['dec'])
                database.photoraw.insert(doc)
            except survey.PhotometryError:
                # couldn't do photometry
                pass


def do_photometry(threads=6,band='g'):
    """
    Do the photometry for all of the stars in all of the observations

    Parameters
    ----------
    threads : int (default: 6)
        Number of threads to run

    band : str (default: 'g')
        Which band to do photometry in

    Warnings
    --------
    If you change this, change force_photometry too!

    History
    -------
    2011-06-14 - Created by Dan Foreman-Mackey

    """
    global _meta
    _meta['band'] = band
    np.random.seed()
    #database.photoraw.drop()
    #observations = survey.find_all_observations()
    print "going?"
    observations = survey.find_observations(-23.431966,-0.227934)
    print len(observations),"runs"

    strt = timer.time()
    pool = multiprocessing.Pool(threads)
    pool.map(photometer_run,observations)
    print "took",timer.time()-strt

    # create the indexes
    print "do_photometry: generating indexes"
    database.photoraw.create_index([('pos',pymongo.GEO2D)])
    database.photoraw.create_index('obsid')
    database.photoraw.create_index('starid')
    database.photoraw.create_index('mjd_g')
    database.photoraw.create_index([('run',1), ('camcol',1), ('field',1)])
    database.photoraw.create_index('mgroup')
    database.photoraw.create_index('rank')

def find_photometry(ra,dec,radius,band='g',mgroup=None,resample=None,minnum=3):
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

    Optional
    --------
    mgroup : int
        Select only measurements in a particular measurement group

    resample : int
        Only select

    minnum : int
        What's the minimum number of observations that we need of a star?

    band : str (default: 'g')
        SDSS band

    Returns
    -------
    obsids : list
        List of (run,camcol) tuples

    stars : list
        List of bson.ObjectID objects for stars

    History
    -------
    2011-07-04 - Created by Dan Foreman-Mackey

    """
    radius /= 60. # to degrees
    radius = np.radians(radius)

    q = {'pos':{'$within':{'$centerSphere': [[ra,dec],radius]}}}
    q['band'] = band
    if mgroup is not None:
        q['mgroup'] = mgroup
    res = database.photoraw.find(q,{'obsid':1,'starid':1})
    if resample is not None:
        res = res.sort([('rank',pymongo.ASCENDING)])

    tries = 0
    while tries <= 1:
        tries += 1

        try:
            obsids = set([])
            stars  = []

            for doc in res:
                if resample is not None and len(stars) >= resample\
                        and doc['starid'] not in stars:
                    break
                print "Warning: FIXME: Line 232 of photometry.py"
                obsids.add((doc['run'],doc['camcol']))
                #obsids.add(tuple(doc['obsid']))
                stars.append(doc['starid'])
            break
        except pymongo.errors.OperationFailure as e:
            raise e
            print 'failed... too many results'
            q['rank'] = {'$lt': 0.5}
            res = database.photoraw.find(q,
                    {'obsid':1,'starid':1}).sort([('rank',pymongo.ASCENDING)])

    star_set = set(stars)
    for obj in list(star_set):
        if stars.count(obj) < minnum:
            print "Warning: less than %d observations of star %r"%(minnum,obj)
            star_set.remove(obj)
    return list(obsids),list(star_set)

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
    data : dict
        model : numpy.ndarray (shape: [len(observations),len(stars),4]
            The mean values for the photometric variables
        cov : numpy.ndarray (shape: [len(observations),len(stars),4,4]
            The covariance matrix

    History
    -------
    2011-07-01 - Created by Dan Foreman-Mackey

    """
    data = {}
    for oi,obsid in enumerate(observations):
        for si,starid in enumerate(stars):
            entry = database.photoraw.find_one({'obsid': obsid, 'starid': starid})
            if entry is not None and entry['cov'][1][1] > 0:
                if 'model' not in data:
                    dim = len(entry['model'])
                    data['model'] = np.zeros([len(observations),len(stars),
                                              dim])
                    data['cov']   = np.zeros([len(observations),len(stars),
                                              dim,dim])
                # data[oi,si,:] = [entry['model'][1],1/entry['cov'][1][1]]
                data['model'][oi,si,:] = entry['model']
                data['cov'][oi,si,:,:] = entry['cov']
    return data

if __name__ == '__main__':
    for b in ['r']:
        do_photometry(band = b)

#    obs,stars = find_photometry(21,0,5)
#    print get_photometry(obs,stars)

