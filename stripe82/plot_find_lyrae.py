#!/usr/bin/env python
# encoding: utf-8
"""
Plot the results of find_lyrae

NOTE: this should only be a temporary script... merge this with find lyrae eventually

History
-------
2011-08-16 - Created by Dan Foreman-Mackey

"""

__all__ = ['LyraeSearch']

import os
import os.path
import cPickle as pickle

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

from calib_results import CalibrationPatch
from calibration import PhotoModel

class LyraeSearch:
    def __init__(self, bp):
        self._bp = bp
        self._patches = []
        self._periods = []
        self._ra  = []
        self._dec = []
        self._params = {}
        for k in ["jitterabs2","jitterrel2","pbad","pvar","sigbad2","Q2"]:
            self._params[k] = []
        for fn in os.listdir(os.path.join(bp,'models')):
            number,ext = os.path.splitext(fn)
            if ext == '.pkl':
                print number
                data,vector,ra,dec,radius,period =\
                        pickle.load(open(os.path.join(bp,'models',fn),'rb'))
                model = PhotoModel(data,vector)
                self._periods.append(period)
                patch = CalibrationPatch(model,model.data,ra,dec,radius)
                self._patches.append(patch)
                self._ra.append(ra)
                self._dec.append(dec)
                for k in list(self._params):
                    self._params[k].append(getattr(model,k))

    def plot_model_params(self):
        for k in list(self._params):
            pl.figure()
            if k in ["pbad","pvar"]:
                pl.plot(self._ra,self._params[k],'.k')
                mu = np.median(self._params[k])
                pl.ylabel(k)
            else:
                pl.plot(self._ra,np.log(self._params[k]),'.k')
                mu = np.median(np.log(self._params[k]))
                pl.ylabel("$\ln$ "+k)
            pl.gca().axhline(mu,color="k",ls="--")
            pl.xlabel("R.A.")
            pl.title("%f"%mu)
            pl.savefig(os.path.join(self._bp,"%s.png"%k))

    def plot_lightcurves(self):
        constbp = os.path.join(self._bp,'plots','constant')
        varbp = os.path.join(self._bp,'plots','variable')
        if not os.path.exists(constbp):
            os.makedirs(constbp)
        if not os.path.exists(varbp):
            os.makedirs(varbp)
        pl.figure()
        for number,patch in enumerate(self._patches):
            print number,
            target_id,target_lc = patch.get_target()

            pl.clf()
            target_lc.plot(ax=pl.gca(),period=self._periods[number],nperiods=2,
                    calcspline=True,hyperparams=True)

            lnoddsvar = target_lc.get_lnoddsvar()
            if lnoddsvar > 0:
                print "lnoddsvar =",lnoddsvar,
                pl.savefig(os.path.join(varbp,'%s.png'%number))
            else:
                pl.savefig(os.path.join(constbp,'%s.png'%number))
            print

if __name__ == '__main__':
    import sys
    search = LyraeSearch(sys.argv[1])
    search.plot_model_params()
    search.plot_lightcurves()


