#!/usr/bin/env python
# encoding: utf-8
"""
Test the calibration

History
-------
2011-08-14 - Created by Dan Foreman-Mackey

"""

__all__ = []

import numpy as np
np.random.seed()

from calibration import PhotoData,do_calibration,lnprob,odds_bad,odds_variable
from calibration.opt import *


class SyntheticData(PhotoData):
    def __init__(self,nstars,nobs,nvar,pbad):
        self.nstars = nstars
        self.nobs   = nobs
        self.stars = []
        self.zeros = 3000+500*np.random.rand(nobs)
        self.flux = np.zeros([nobs,nstars])
        self.ivar = np.zeros([nobs,nstars])
        self.truth = []
        self.err = 0.1
        self.jitterabs2 = 1.0
        self.jitterrel2 = 0.0001 # 1e-5
        for i in range(nstars):
            star = {}
            star['g'] = 20+2*np.random.randn()
            star['Err_g'] = 0.1
            self.stars.append(star)

            self.truth.append(mag2nmgy(star['g']))
            self.flux[:,i] = self.truth[-1]*self.zeros
            self.flux[:,i] += \
                    np.sqrt(self.err**2+self.jitterrel2*self.truth[-1]**2+\
                        self.jitterabs2)*np.random.randn(nobs)
            self.ivar[:,i] = 1/1.0/(self.err)**2
        self.magprior = np.array([[s['g'],s['Err_g']**2] for s in self.stars])
        self.sigvar2 = (100.+5.*np.random.randn())**2
        print "sigvar_true...",self.sigvar2
        self.nvar = nvar
        self.pvar = float(self.nvar)/self.nstars
        for i in range(nvar):
            self.flux[:,i] += np.sqrt(self.sigvar2)*np.random.randn(nobs)
        self.pbad = pbad
        self.nbad = int(nobs*nstars*pbad)
        self.sigbad2 = (5+5.*np.random.rand())**2
        inds = set([])
        while len(inds) < self.nbad:
            inds.add((np.random.randint(nobs),np.random.randint(nstars)))
        for n in inds:
            self.flux[n[0],n[1]] += np.sqrt(self.sigbad2)*np.random.randn()

        self.truth = np.array(self.truth)

    def true_vector(self):
        p0 = np.append(self.zeros,nmgy2mag(self.truth))
        p0 = np.append(p0,[np.log(self.jitterabs2),np.log(self.jitterrel2),
            prob2logodds(self.pvar),np.log(self.sigvar2),
            prob2logodds(self.pbad),np.log(self.sigbad2)])
        return p0


if __name__ == '__main__':
    data = SyntheticData(50,60,1,0.01)
    p0 = None
    for i in range(1):
        if 0: #i == 0:
            fix_probs = [0.01,0.01]
        else:
            fix_probs = None
        model = do_calibration(data,addtodb=False,p0=p0,fix_probs=fix_probs)

        print "zeropoint % error = ",100*(data.zeros-model.zero)/model.zero
        print "odds_var = ",odds_variable(model,data)
        print "fluxes = ",data.truth,model.flux
        print "truth,fit"
        print "pvar = ",float(data.nvar)/data.nstars,model.pvar
        print "sigvar2 = ",data.sigvar2,model.sigvar2
        print "pbad = ",data.pbad,model.pbad
        print "sigbad2 = ",data.sigbad2,model.sigbad2
        print "jitterrel2 = ",data.jitterrel2,model.jitterrel2
        print "jitterabs2 = ",data.jitterabs2,model.jitterabs2

        print "lnprob = ",lnprob(data.true_vector(),data),lnprob(model.vector(),data)
        f_vec = np.array(model.vector())
        f_vec[:data.nstars+data.nobs] = data.true_vector()[:data.nstars+data.nobs]
        print "frankenstein = ",lnprob(f_vec,data)
        #print data.true_vector(),model.vector()
        p0 = model.vector()

