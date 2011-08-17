#!/usr/bin/env python
# encoding: utf-8
"""
A general probabilistic calibration model

History
-------
2011-06-14 - Created by Dan Foreman-Mackey

"""

__all__ = ['PhotoData','PhotoModel',
           'lnprob','lnlike','lnprior',
           'odds_bad','odds_variable','lnlikeratio_bad']

import numpy as np

import _likelihood

# options
from opt import *


# ================================== #
#  Light curve model wrapper classes #
# ================================== #

class PhotoData:
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
    def __init__(self,data,observations,stars):
        print "nobs =",len(observations)
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
        self.magprior = np.array([[s['g'],s['Err_g']**2] for s in self.stars])
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

def unpickle_model(data,vector):
    __safe_for_unpickling__ = True
    return PhotoModel(data,vector)
class PhotoModel:
    """
    Wrapper class around calibration model.

    Note that the actual probability calculations are separate (so that they can
    be used with the multiprocessing module)

    Parameters
    ----------
    data : PhotoData
        The PhotoData object used to constrain the model
        FIXME: this is SERIOUSLY bad form

    vector : list
        A vector of model parameters

    History
    -------
    2011-06-14 - Created by Dan Foreman-Mackey

    """
    def __init__(self,data,vector):
        self.model = 2 # hardcoded to use best model
        self.data = data
        self.conv,self.npars = self.param_names()
        self.from_vector(vector)
        self.pbad = 0.01
        self.pvar = 0.01
        self.sigbad2 = np.exp(17)
        self.Q2 = 0.3**2
        self.jitterabs2 = 0.0
        self.jitterrel2 = 1e-3

    def __reduce__(self):
        return (unpickle_model,(self._data,self.vector()))

    def param_names(self):
        """
        Get the list of model parameters and their conversion functions

        Here, we have function pointers that tell us how to convert from the
        sampling space (often log-space) to the linear space that is probably
        more enlightening. There are also some convenient conversions (e.g.
        mag -> flux)

        Returns
        -------
        conv : dict
            A dictionary where the keys are strings giving the names of the
            parameters and the values are tuples (n,f,invf) where
                n : int
                    The index of this parameter in the sampling vector
                f : function
                    A function pointer to convert from sampling space
                    to linear space
                invf : function
                    The inverse of f

        npars : int
            The number of parameters

        History
        -------
        2011-06-14 - Created by Dan Foreman-Mackey

        """
        no = lambda x: x
        sq = lambda x: x**2
        conv = {'mag': (np.arange(self.data.nobs,self.data.nobs+self.data.nstars)
                    ,no,no)}

        # we sample in magnitudes (ie ~log(flux)) so we need to convert
        # to fluxes and back
        conv['zero'] = (np.arange(self.data.nobs),no,no)
        conv['flux'] = (np.arange(self.data.nobs,self.data.nobs+self.data.nstars),
                mag2nmgy, nmgy2mag)

        # in the more complicated models, we have more parameters
        zero = self.data.nobs+self.data.nstars
        #if self.model >= 1:
        #    #conv['jitterabs2'] = (zero,sq,np.sqrt)
        #    conv['jitterrel2'] = (zero,sq,np.sqrt)#sq,np.sqrt)
        #    #conv['pvar']      = (zero+2,logodds2prob,prob2logodds)
        #    conv['Q2']         = (zero+1,sq,np.sqrt)
        #    #conv['sigbad2']    = (zero+3,sq,np.sqrt)
        #    zero += 2
        #if self.model >= 2:
        #    conv['pbad']       = (zero,logodds2prob,prob2logodds)
        #    zero += 2
        return conv,zero

    def from_vector(self,p0):
        """
        Given a vector of parameter values populate the class attributes

        Parameters
        ----------
        p0 : list
            A vector of parameters given in the same order as self.conv

        History
        -------
        2011-06-15 - Created by Dan Foreman-Mackey

        """
        p0 = np.array(p0)
        for k in self.conv:
            setattr(self,k,self.conv[k][1](p0[self.conv[k][0]]))

    def vector(self):
        """
        The inverse of from_vector

        Returns
        -------
        vec : 1-D numpy.ndarray
            The vector of parameters in the same order as self.conv

        History
        -------
        2011-06-15 - Created by Dan Foreman-Mackey

        """
        vec = np.zeros(self.npars)
        for k in self.conv:
            vec[self.conv[k][0]] = self.conv[k][2](getattr(self,k))
        return vec

    def __unicode__(self):
        st = u"\ng-mags of stars\n"
        st +=   "---------------\n"
        st += repr(self.mag)
        st += "\nzero points of runs [ADU/nMgy]\n"
        st +=   "------------------------------\n"
        st += repr(self.zero)
        st += "\n"
        for k in ['jitterabs2','jitterrel2','pvar','Q2',
                'pbad','sigbad2']:
            st += "%10s\t"%k
            st += "%e\n"%getattr(self,k)
        return st

    def __str__(self):
        return unicode(self)

# ===================== #
#  Likelihood Function  #
# ===================== #

# precalculate this:
# log(1/sqrt(2 pi))
lisqrt2pi = - 0.5*np.log(2.0*np.pi)
def _lnnormal(x,mu,var):
    return -0.5*(x-mu)**2/var - 0.5*np.log(var) + lisqrt2pi

def lnprob(p,data,model=2,fix_probs=None):
    """
    NAME:
        lnprob
    PURPOSE:
        ln-probability of p given data
    INPUT:
        p  - vector for use with PhotoModel object
        data - a PhotoData object
        model - which model should we use
    OUTPUT:

    HISTORY:
        Created by Dan Foreman-Mackey on Jun 07, 2011
    """
    #print p
    params = PhotoModel(data,p)
    if fix_probs is not None:
        params.pvar = fix_probs[0]
        params.pbad = fix_probs[1]
    prior = lnprior(params)
    if np.isinf(prior):
        return -np.inf
    lnpost = prior + lnlike(params,data) #lnlike(params,data)
    if np.isnan(lnpost):
        print params
        print lnpost
        raise Exception()
    return lnpost

def lnprior(p):
    """
    NAME:
        lnprior
    PURPOSE:
        the prior on p
    INPUT:
        p - a PhotoParams object
    OUTPUT:
        the ln-prior
    HISTORY:
        Created by Dan Foreman-Mackey on Jun 07, 2011
    """
    if not (0 <= p.pbad <= 1 and 0 <= p.pvar <= 1):
        return -np.inf
    # g-band magnitude prior
    modelmag = p.mag
    mag = p.data.magprior[:,0]
    err = 0.03 # MAGIC p.data.magerr[:,1]
    lnprior = -0.5*np.sum((modelmag-mag)**2/err+np.sum(np.log(err)))

    return lnprior

def lnlike(p,data):
    return _likelihood.lnlikelihood(p,data)

def odds_bad(p,data):
    oarr = np.zeros(np.shape(data.flux))
    _likelihood.lnoddsbad(p,data,oarr)
    return oarr

def odds_variable(p,data):
    oarr = np.zeros(data.nstars)
    _likelihood.lnoddsvar(p,data,oarr)
    return oarr

def lnlikeratio_bad(p,data):
    oarr = np.zeros(np.shape(data.flux))
    _likelihood.lnlikeratiobad(p,data,oarr)
    return oarr

