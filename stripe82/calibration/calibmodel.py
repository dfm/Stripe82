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

from database import obslist

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
        splineorder = 3

        self.ra   = dict([])
        self.dec  = dict([])
        self.zero = dict([])
        self.counts = dict([])

        for doc in obslist.find(query):
            runid = "%05d%d"%(doc['run'],doc['camcol'])
            if runid not in self.ra:
                self.ra[runid]   = []
                self.dec[runid]  = []
                self.zero[runid] = []
                self.counts[runid] = []

            if doc['pos']['ra'] not in self.ra[runid]:
                self.ra[runid].append(doc['pos']['ra'])
                self.dec[runid].append([doc['pos']['dec']])
                self.zero[runid].append(doc['zero'])
                self.counts[runid].append(1)
            else:
                ind = self.ra[runid].index(doc['pos']['ra'])
                self.zero[runid][ind] += doc['zero']
                self.counts[runid][ind] += 1
                self.dec[runid][ind].append(doc['pos']['dec'])

        self.splines = dict([])
        for k in list(self.ra):
            self.ra[k] = np.array(self.ra[k])
            self.zero[k] = np.array(self.zero[k])/np.array(self.counts[k])
            
            if np.shape(self.ra[k])[0] > splineorder+1:
                try:
                    self.splines[k] = inter.interp1d(
                        self.ra[k],self.zero[k],
                        kind=splineorder)
                except:
                    pass



