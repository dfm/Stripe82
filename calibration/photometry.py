#!/usr/bin/env python
# encoding: utf-8
"""
Do batch photometry for a set of stars in a set of fields

History
-------
2011-06-14 - Created by Dan Foreman-Mackey

"""

__all__ = ['force_photometry','photometry']

import numpy as np

# options
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
    If you change this, change photometry too!
    
    History
    -------
    2011-06-22 - Created by Dan Foreman-Mackey
    
    """
    info = survey.get_observation(observation)
    
    if info['ramin'] < ra < info['ramax'] and \
            info['decmin'] < dec < info['decmax']:
        try:
            obs = survey.Observation(observation)
        except survey.ObservationAccessError:
            return [0,0]

        try:
            res = obs.photometry(ra,dec)
            return res
        except survey.PhotometryError:
            return [0,0]

    return [0,0]
    

def photometry(observations,stars):
    """
    Force photometry at given positions in given observations
    
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
    
    Warnings
    --------
    If you change this, change photometry too!

    History
    -------
    2011-06-14 - Created by Dan Foreman-Mackey
    
    """
    data = np.zeros([len(observations),len(stars),2])
    for oi,obsid in enumerate(observations):
        obs = -1
        for si,starid in enumerate(stars):
            entry = database.photoraw.find_one({'obsid': obsid, 'starid': starid})
            if entry is not None:
                data[oi,si,:] = [entry['counts'],entry['invvar']]
            else:
                star = survey.get_star(starid)
                info = survey.get_observation(obsid)
                if info['ramin'] < star['ra'] < info['ramax'] and \
                        info['decmin'] < star['dec'] < info['decmax']:
                    if obs is -1:
                        try:
                            obs = survey.Observation(obsid)
                        except survey.ObservationAccessError:
                            obs = None
                    if obs is not None:
                        try:
                            res = obs.photometry(star['ra'],star['dec'])
                            data[oi,si,:] = res
                            database.photoraw.insert({'obsid':obsid,'starid':starid,
                                    'counts':data[oi,si,0],'invvar':data[oi,si,1]})
                        except survey.PhotometryError:
                            res = None
                    if obs is None or res is None:
                        database.photoraw.insert({'obsid':obsid,'starid':starid,
                                'counts':0,'invvar':0})
    return data

if __name__ == '__main__':
    pass
    #database.photoraw.drop()
    #ra0,dec0 = -23.431965,-0.227934
    #for ra in np.arange(-49.5,58.5,1.):
    #    print "R.A. ==========> ",ra
    #    dec = -0.227934
    #    stars = survey.find_stars(ra,dec)
    #    observations = survey.find_observations(ra,dec)
    #    print len(stars)," stars and ",len(observations)," observations"
    #    force_photometry(observations,stars)

