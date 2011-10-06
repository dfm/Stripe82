#!/usr/bin/env python
# encoding: utf-8
"""
The probabilistic patch model

"""

__all__ = ['PatchData','PatchModel',
           'lnprob','lnlike','lnprior',
           'odds_bad','odds_variable','lnlikeratio_bad']

import numpy as np

import _likelihood

# options
from conversions import *

# ================================== #
#  Light curve model wrapper classes #
# ================================== #

class PatchData:
    """
    Wrapper class for the photometric data

    Parameters
    ----------
    data : numpy.ndarray (shape: [nobs,nstars,2])
        Matrix of observations of stars. The last axis is (counts,inverse variance)

    observations : list
        List of bson.ObjectID objects for the observations (same order as data)

    stars : list
        List of bson.ObjectID objects for the stars (same order as data)

    History
    -------
    2011-06-14 - Created by Dan Foreman-Mackey

    """
    def __init__(self,data,observations,stars,band='g'):
        self.band = band
        self.data = data
        self.stars = []
        # mean ra and dec
        self.ra  = 0.0
        self.dec = 0.0
        for sid in stars:
            self.stars.append(survey.get_star(sid))
            self.ra  += self.stars[-1]['ra']
            self.dec += self.stars[-1]['dec']
        self.nstars = len(stars)
        self.ra /= self.nstars
        self.dec /= self.nstars
        self.observations = []
        for oid in observations:
            doc = survey.get_observation(oid)
            self.observations.append(doc)
        self.magprior = np.array([[s[band],0.0] for s in self.stars])
        self.flux = data['model'][:,:,1]
        tmp = data['cov'][:,:,1,1]
        self.ivar = np.zeros(np.shape(tmp))
        self.ivar[tmp > 0] = 1.0/tmp[tmp > 0]
        self.nobs = len(observations) #np.shape(data)[0]

    def mjd(self):
        ra,dec = self.ra,self.dec
        mjds = []
        for obs in self.observations:
            try:
                mjds.append(survey.db.obsdb.find_one(
                    {'run':obs['run'],'camcol':obs['camcol'],
                    'ramin': {'$lt': ra}, 'ramax': {'$gt': ra},
                    'decmin': {'$lt': dec}, 'decmax': {'$gt': dec}},
                    {'mjd_g': 1})['mjd_g'])
            except TypeError:
                mjds.append(0.0)
        return np.array(mjds)


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

    obs_err : numpy.ndarray (M, N)
        The uncertainties on obs_flux

    calib_mag : numpy.ndarray (N,)
        The cataloged values for the magnitude of the star

    """
    model_id = 1

    def __init__(self, time, obs_flux, obs_err, mag_prior):
        self._t = time
        self._f = obs_flux
        self._sig_f = obs_err
        self._mag_prior = mag_prior

        self._nstars = self._f.shape[0]
        self._nruns  = self._t.size

        self._init_params()

    def init_params(self):
        tmp = np.mean(self._f/mag2nmgy(self._mag_prior), axis=-1)
        self.vector = ma.concatenate([tmp,data.magprior[:,0]])

        # nuisance parameters
        self._pbad = 0.01
        self._pvar = 0.01
        self._sigbad2 = np.exp(17)
        self._Q2 = 0.3**2
        self._jitterabs2 = 0.0
        self._jitterrel2 = 1e-3

    def calibrate(self):
        p0 = self.vector
        p1 = op.fmin_bfgs(chi2,p0)
        self.vector = p1

    @property
    def vector(self):
        return np.concatenate(self._f0, self._fstar)

    @vector.setter
    def set_vector(self, vector):
        self._f0    = vector[:self._nruns]
        self._mstar = vector[self._nruns:]
        self._fstar = mag2nmgy(self._mstar)

    @property
    def nuisance(self):
        return [self._pbad, self._pvar, self._sigbad2, self._Q2,
                self._jitterabs2, self._jitterrel2]

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
        self.vector = p
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
        return -0.5*np.sum((sefl._mstar-self._mag_prior)**2/err+np.sum(np.log(err)))

    def odds_bad(p,data):
        oarr = np.zeros(np.shape(self._f))
        _likelihood.lnoddsbad(self, oarr)
        return oarr

    def odds_variable(p,data):
        oarr = np.zeros(data._nstars)
        _likelihood.lnoddsvar(self, oarr)
        return oarr

    def lnlikeratio_bad(p,data):
        oarr = np.zeros(np.shape(data._f))
        _likelihood.lnlikeratiobad(self, oarr)
        return oarr

