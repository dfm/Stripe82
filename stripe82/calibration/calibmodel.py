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

try:
    import progressbar
except:
    progressbar = None

from database import obslist,photomodel
from photomodel import PhotoModel

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
        self.ra   = {}
        self.dec  = {}
        self.zero = {}
        self.pnames = ['jitterabs2','jitterrel2','pvar','sigvar2',
                'pbad','sigbad2']
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
        
        for i,doc in enumerate(obslist.find(query)):
            if progressbar is not None:
                pbar.update(i)
            runid = "%05d%d"%(doc['run'],doc['camcol'])
            if runid not in self.ra:
                self.ra[runid]   = []
                self.dec[runid]  = []
                self.zero[runid] = []
                for p in self.pnames:
                    self.pars[p][runid] = []

            self.ra[runid].append(doc['pos']['ra'])
            self.dec[runid].append([doc['pos']['dec']])
            self.zero[runid].append(doc['zero'])
            mod = photomodel.find_one({'_id': doc['modelid']})['model']
            for p in self.pnames:
                self.pars[p][runid].append(getattr(mod,p))
        if progressbar is not None:
            pbar.finish()


