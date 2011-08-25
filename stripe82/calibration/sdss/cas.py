#!/usr/bin/env python
# encoding: utf-8
"""
Interface to list of stars and fields from CAS

This uses a mongodb database generated by the populate_cas_db.py script

History
-------
2011-06-13 - Created by Dan Foreman-Mackey

"""

__all__ = ['find_stars','find_observations','get_star','get_observation',
        'find_all_stars','find_all_observations','find_stars_in_observation']

import numpy as np

import db
import pymongo

def _angular_distance(ra1,dec1,ra2,dec2):
    """
    Calculate the angular separation of 2 objects in RA/Dec space

    Parameters
    ----------
    ra1 : float
        In degrees

    dec1 : float
        In degrees

    ra2 : float
        In degrees

    dec2 : float
        In degrees

    Returns
    -------
    rho : float
        In degrees

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey

    """
    ra1,dec1 = np.radians(ra1),np.radians(dec1)
    ra2,dec2 = np.radians(ra2),np.radians(dec2)
    crho = np.cos(dec1)*np.cos(dec2)\
            *(np.cos(ra1)*np.cos(ra2)+np.sin(ra1)*np.sin(ra2))\
            +np.sin(dec1)*np.sin(dec2)
    return np.degrees(np.arccos(crho))

def find_stars(ra,dec,radius=3.0):
    """
    Find bright stars within a specific annulus around (ra/dec)

    Parameters
    ----------
    ra : float
        In degrees

    dec : float
        In degrees

    Optional
    --------
    radius : float (default = 3.0)
        In arcmin

    Returns
    -------
    stars : list
        List of ObjectId objects

    TODO
    ----
    - Check that ra == 360 works properly

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey

    """
    radius /= 60. # to degrees
    radius = np.radians(radius)
    while ra > 180.:
        ra -= 360.0

    res = db.stardb.find({'pos':{'$within':{'$centerSphere': [[ra,dec],radius]}}},
            {'_id': 1})
    stars = [star['_id'] for star in res]

    return stars

def find_observations(ra,dec):
    """
    Find all Stripe 82 fields that overlap with (RA/Dec)

    Parameters
    ----------
    ra : float
        In degrees

    dec : float
        In degrees

    Returns
    -------
    observations : list
        List of ObjectId objects

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey

    """
    res = db.obsdb.find({'ramin': {'$lt': ra}, 'ramax': {'$gt': ra},
                        'decmin': {'$lt': dec}, 'decmax': {'$gt': dec}},
                        {'run':1,'camcol':1})
    observations = [(obs['run'],obs['camcol']) for obs in res]
    return list(set(observations))

def get_star(objid):
    """
    Retrieve entry for star with given id

    Parameters
    ----------
    objid : bson.ObjectId
        The object ID

    Returns
    -------
    star : dict
        The star entry

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey

    """
    return db.stardb.find_one({'_id': objid})

def get_observation(obsid):
    """
    Retrieve entry for observation with given id

    Parameters
    ----------
    obsid : tuple
        (run,camcol)

    Returns
    -------
    observation : dict
        The observation entry

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey

    """
    return db.obsdb.find_one(dict(zip(('run','camcol'),obsid)))

def find_stars_in_observation(obsid):
    """
    Find all the stars within the bounds of a given observation

    Parameters
    ----------
    obsid : BSON.ObjectId
        The identifier for the field

    Returns
    -------
    stars : list
        List of BSON.ObjectIDs

    History
    -------
    2011-07-03 - Created by Dan Foreman-Mackey

    """
    q = dict(zip(('run','camcol'),obsid))
    keys = {'ramin':pymongo.ASCENDING,'ramax':pymongo.DESCENDING,
            'decmin':pymongo.ASCENDING,'decmax':pymongo.DESCENDING}
    rect = {}
    for k in list(keys):
        fields = db.obsdb.find(q,{k:1}).sort([(k,keys[k])]).limit(1)
        rect[k] = fields[0][k]
    return [star['_id'] for star in db.stardb.find({'pos': {'$within': {'$box':
        [[rect['ramin'],rect['decmin']],[rect['ramax'],rect['decmax']]]}}},
        {'_id':1})]


# RANDOM SHITE: FIXME
# 2 faint sesar lyrae, bright(ish) Sesar lyrae, random other faint source
#ap = [(29.47942,0.383557),(17.114624,1.05495),(43.511726,-0.420139),
#        (10.0018734334081,0.791580301596976)]
ap = [(29.47942,0.383557),(10.0018734334081,0.791580301596976)]

def find_all_stars():
    """
    Retrieve list of objids for all stars

    Returns
    -------
    stars : list
        List of BSON.ObjectIDs

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey

    """
    # FIXME FIXME FIXME
    ret = []
    for i in range(len(ap)):
        ret += find_stars(ap[i][0],ap[i][1],30)
    return ret

    return [star['_id'] for star in db.stardb.find(q)]

def find_all_observations():
    """
    Retrieve list of objids for all observations

    Returns
    -------
    observations : list
        List of BSON.ObjectIDs

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey

    """
    return list(set([(obs['run'],obs['camcol']) for obs in db.obsdb.find()]))


