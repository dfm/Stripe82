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
        splineorder = 1

        self.runs = dict([])
        for doc in obslist.find(query):
            runid = "%05d%d"%(doc['run'],doc['camcol'])
            if runid not in self.runs:
                self.runs[runid] = []
            self.runs[runid].append([doc['pos']['ra'],doc['pos']['dec'],doc['zero']])
        self.splines = dict([])
        for k in list(self.runs):
            self.runs[k] = np.array(self.runs[k])
            if np.shape(self.runs[k])[0] >= (splineorder+1):
                self.splines[k] = inter.interp1d(
                        self.runs[k][:,0], self.runs[k][:,2],
                        kind='linear')



