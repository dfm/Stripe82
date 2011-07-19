#!/usr/bin/env python
# encoding: utf-8
"""
A general probabilistic calibration model

History
-------
2011-06-14 - Created by Dan Foreman-Mackey

"""

__all__ = ['CalibrationModel']

import numpy as np
import scipy.interpolate as inter

try:
    import progressbar
except:
    progressbar = None

from database import obslist,photomodel,calibcache

# ============================== #
# Full calibration model wrapper #
# ============================== #

class CalibrationModel:
    """
    Encapsulates and provides an interface to the general calibration model
    
    Parameters
    ----------
    query : dict
        Query for the photometry.obslist database
    
    History
    -------
    2011-07-05 - Created by Dan Foreman-Mackey
    
    """
    def __init__(self,query):
        self.pnames = ['jitterabs2','jitterrel2','pvar','sigvar2',
                    'pbad','sigbad2']
        
        doc = calibcache.find_one(query)
        if doc is not None:
            self.ra   = doc['ra']
            self.dec  = doc['dec']
            self.zero = doc['zero']
            self.pars = doc['pars']
        else:
            self.ra   = {}
            self.dec  = {}
            self.zero = {}
            self.pars = {}
            for p in self.pnames:
                self.pars[p] = {}

            if progressbar is not None:
                tot = obslist.find(query).count()
                widgets = ['Loading control points: ',
                        progressbar.Percentage(), ' ',
                        progressbar.Bar(marker=progressbar.RotatingMarker()),
                        ' ', progressbar.ETA()]
                pbar = progressbar.ProgressBar(widgets=widgets,
                        maxval=tot).start()
            
            for i,doc in enumerate(obslist.find(query,timeout=False)):
                if progressbar is not None:
                    pbar.update(i)
                runid = doc['obsid']
                if runid not in self.ra:
                    self.ra[runid]   = []
                    self.dec[runid]  = []
                    self.zero[runid] = []
                    for p in self.pnames:
                        self.pars[p][runid] = []

                self.ra[runid].append(doc['pos']['ra'])
                self.dec[runid].append(doc['pos']['dec'])
                self.zero[runid].append(doc['zero'])
                mod = photomodel.find_one({'_id': doc['modelid']})['model']
                for p in self.pnames:
                    self.pars[p][runid].append(getattr(mod,p))
            if progressbar is not None:
                pbar.finish()

            doc = query
            doc['ra'] = self.ra
            doc['dec'] = self.dec
            doc['zero'] = self.zero
            doc['pars'] = self.pars
            calibcache.insert(doc)

        self.splines = {}
        self.uniquedecs = {}
        for k in self.ra:
            # HACKHACKHACK FIXME:
            self.dec[k] = [d[0] for d in self.dec[k]]
            # END HACK
            kind = 'cubic'
            decs = list(set(self.dec[k]))
            self.uniquedecs[k] = []
            for dec in decs:
                inds = np.array(self.dec[k]) == dec
                x = np.array(self.ra[k])[inds]
                if np.shape(x)[0] >= 4:
                    if k not in self.splines:
                        self.splines[k] = []
                    self.uniquedecs[k].append(dec)
                    self.splines[k].append({})
                    sinds = np.argsort(x)
                    y = np.array(self.zero[k])[inds][sinds]
                    x = x[sinds]
                    self.splines[k][-1]['zero'] = \
                            inter.interp1d(x,y,kind=kind,
                                    fill_value=np.mean(y),
                                    bounds_error=False)
                    for p in self.pnames:
                        y = np.array(self.pars[p][k])[inds][sinds]
                        self.splines[k][-1][p] = \
                                inter.interp1d(x,y,kind=kind,
                                        fill_value=np.median(y),
                                        bounds_error=False)
                                        

