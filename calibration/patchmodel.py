#!/usr/bin/env python
# encoding: utf-8
"""
The probabilistic patch model

"""

__all__ = ['PatchProbModel', 'PatchMedianModel']

import numpy as np
import numpy.ma as ma

import scipy.optimize as op

import _likelihood

# options
from conversions import mag2nmgy, nmgy2mag

# ================================== #
#  Light curve model wrapper classes #
# ================================== #

class PatchProbModel(object):
    """
    Provide a probabilistic model for the calibration of lightcurves

    N stars and M runs.

    Parameters
    ----------
    time : numpy.ndarray (M,)
        List of timestamps for the observations

    obs_flux : numpy.ndarray (M, N)
        The observed lightcurves

    obs_ivar : numpy.ndarray (M, N)
        The inverse variances of obs_flux

    calib_mag : numpy.ndarray (N,)
        The cataloged values for the magnitude of the star

    """
    model_id = 1

    def __init__(self, time, obs_flux, obs_ivar, mag_prior):
        self._t = time
        self._f = obs_flux
        self._ivar_f = obs_ivar
        self._mag_prior = mag_prior

        self._nstars = self._f.shape[0]
        self._nruns  = self._t.size

        self.init_params()

    def init_params(self):
        tmp = np.mean(self._f/mag2nmgy(self._mag_prior), axis=-1)
        vector = ma.concatenate([tmp,self._mag_prior])

        self._f0    = vector[:self._nruns]
        self._mstar = vector[self._nruns:]
        self._fstar = mag2nmgy(self._mstar)

        # nuisance parameters
        self._pbad = 0.01
        self._pvar = 0.01
        self._sigbad2 = np.exp(17)
        self._Q2 = 0.3**2
        self._jitterabs2 = 0.0
        self._jitterrel2 = 1e-3

    def calibrate(self):
        p0 = self.vector
        chi2 = lambda p: -self(p)
        p1 = op.fmin_bfgs(chi2,p0)
        self.set_vector(p1)

    @property
    def vector(self):
        return np.concatenate((self._f0, self._mstar))

    def set_vector(self, vector):
        self._f0    = vector[:self._nruns]
        self._mstar = vector[self._nruns:]
        self._fstar = mag2nmgy(self._mstar)

    @property
    def zeros(self):
        return self._f0

    @property
    def star_flux(self):
        return self._fstar

    @property
    def star_mag(self):
        return self._mstar

    @property
    def nuisance(self):
        return [self._pbad, self._pvar, self._sigbad2, self._Q2,
                self._jitterabs2, self._jitterrel2]

    @nuisance.setter
    def set_nuisance(self, value):
        self._pbad, self._pvar, self._sigbad2, self._Q2, \
                self._jitterabs2, self._jitterrel2 = value

    def __unicode__(self):
        st = u"\ng-mags of stars\n"
        st +=   "---------------\n"
        st += repr(self._mstar)
        st += "\nzero points of runs [ADU/nMgy]\n"
        st +=   "------------------------------\n"
        st += repr(self._f0)
        st += "\n"
        for k in ['_jitterabs2','_jitterrel2','_pvar','_Q2',
                '_pbad','_sigbad2']:
            st += "%10s\t"%k
            st += "%e\n"%getattr(self,k)
        return st

    def __str__(self):
        return unicode(self)

    def __call__(self, p):
        self.set_vector(p)
        prior = self.lnprior()
        if np.isinf(prior):
            return -np.inf
        lnpost = prior + _likelihood.lnlikelihood(self)
        return lnpost

    def lnprior(self):
        if not (0 <= self._pbad <= 1 and 0 <= self._pvar <= 1):
            return -np.inf
        # g-band magnitude prior
        err = 0.03
        return -0.5*np.sum((self._mstar-self._mag_prior)**2/err+np.sum(np.log(err)))

    def odds_bad(self, p,data):
        oarr = np.zeros(np.shape(self._f))
        _likelihood.lnoddsbad(self, oarr)
        return oarr

    def odds_variable(self, p,data):
        oarr = np.zeros(data._nstars)
        _likelihood.lnoddsvar(self, oarr)
        return oarr

    def lnlikeratio_bad(self, p,data):
        oarr = np.zeros(np.shape(data._f))
        _likelihood.lnlikeratiobad(self, oarr)
        return oarr

class PatchMedianModel(PatchProbModel):
    model_id = 0
    
    def init_params(self):
        self._f0 = np.sum(self._ivar_f \
                * self._f/mag2nmgy(self._mag_prior)[None,:], axis=-1) \
                / np.sum(self._ivar_f, axis=-1)
        self._f0[self._f0 < 0] = np.abs(np.median(self._f0))
        self._fstar = np.sum(self._ivar_f * self._f/self._f0[:,None], axis=0) \
                / np.sum(self._ivar_f, axis=0)
        self._fstar[self._fstar == 0.0] = mag2nmgy(self._mag_prior[self._fstar == 0.0])
        self._mstar = nmgy2mag(self._fstar)

        # nuisance parameters
        self._pbad = 0.01
        self._pvar = 0.01
        self._sigbad2 = np.exp(17)
        self._Q2 = 0.3**2
        self._jitterabs2 = 0.0
        self._jitterrel2 = 1e-3

    def calibrate(self, maxiter=1000, tol=1.25e-8):
        iv0 = np.sum(self._ivar_f, axis=0)
        iv1 = np.sum(self._ivar_f, axis=1)
        for i in xrange(maxiter):
            f0 = np.sum(self._ivar_f * self._f/self._fstar[None,:], axis=-1)/iv1
            f0[f0 < 0] = np.abs(np.median(f0))
            self._fstar = np.sum(self._ivar_f * self._f/f0[:,None], axis=0)/iv0
            self._fstar[self._fstar == 0.0] = mag2nmgy(self._mag_prior[self._fstar == 0.0])
            if i > 5 and np.sum(np.fabs(f0-self._f0)/f0) < tol:
                self._f0 = f0
                break
            self._f0 = f0

        if i < maxiter-1:
            print "Median model converged in %d iterations"%(i)
        else:
            print "Warning: Median model didn't converge after %d iterations"%(i)
            
        self._mstar = nmgy2mag(self._fstar)

