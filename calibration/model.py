#!/usr/bin/env python
# encoding: utf-8
"""
A general probabilistic calibration model

History
-------
2011-06-14 - Created by Dan Foreman-Mackey

"""

__all__ = ['PhotoSONManipulator','PhotoData','PhotoModel',
           'lnprob','lnprior','lnlike',
           'lnprob_badobs','lnprob_variable']

import cPickle as pickle

import numpy as np

from bson.binary import Binary
from pymongo.son_manipulator import SONManipulator

# options
from opt import survey

# ================#
#  BSON Encodings #
# ================#

class PhotoSONManipulator(SONManipulator):
    """
    SON manipulator to deal with my custom data types
    
    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html#automatic-encoding-and-decoding
    [2] https://github.com/FlaPer87/django-mongodb-engine/blob/master/django_mongodb_engine/serializer.py
    
    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    def transform_incoming(self, value, collection):
        if isinstance(value, (list,tuple,set)):
            return [self.transform_incoming(item,collection) for item in value]
        if isinstance(value,dict):
            return dict((key,self.transform_incoming(item,collection))
                    for key,item in value.iteritems())
        if isinstance(value,PhotoModel):
            return encode_model(value)
        if isinstance(value,PhotoData):
            return encode_data(value)
        return value

    def transform_outgoing(self, son, collection):
        if isinstance(son,(list,tuple,set)):
            return [self.transform_outgoing(value,collection) for value in son]
        if isinstance(son,dict):
            if son.get('_type') == 'PhotoModel':
                return decode_model(son)
            if son.get('_type') == 'PhotoData':
                return decode_data(son)
            return dict((key,self.transform_outgoing(value,collection))
                    for key,value in son.iteritems())
        return son

def encode_model(model):
    """
    Encode the model parameters for BSON

    Parameters
    ----------
    model : PhotoModel
        The model object
    
    Returns
    -------
    encoded : dict
        Encoded version of the class

    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html
    
    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    return {'_type': 'PhotoModel', 'vector': list(model.vector()),
            'data': encode_data(model.data)}

def decode_model(document):
    """
    Decode a BSON PhotoModel document
    
    Parameters
    ----------
    document : dict
        BSON document
    
    Returns
    -------
    model : PhotoModel
        The model object
    
    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html

    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    assert document['_type'] == 'PhotoModel'
    return PhotoModel(decode_data(document['data']),document['vector'])

def encode_data(data):
    """
    Encode the PhotoData class for BSON
    
    Parameters
    ----------
    data : PhotoData
        The data object
    
    Returns
    -------
    encoded : dict
        Save-able by BSON
    
    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html

    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    return {'_type': 'PhotoData','data': Binary(pickle.dumps(data.data,-1)),
            'observations': [obs['_id'] for obs in data.observations],
            'stars': [star['_id'] for star in data.stars]}

def decode_data(document):
    """
    Decode a BSON PhotoData document
    
    Parameters
    ----------
    document : dict
        BSON document
    
    Returns
    -------
    data : PhotoData
        The data object
    
    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html

    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    assert document['_type'] == 'PhotoData'
    return PhotoData(pickle.loads(document['data']),document['observations'],
        document['stars'])

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
        self.data = data
        self.stars = []
        for sid in stars:
            self.stars.append(survey.get_star(sid))
        self.observations = []
        for oid in observations:
            self.observations.append(survey.get_observation(oid))
        self.magprior = np.array([[s['g'],s['Err_g']**2] for s in self.stars])
        self.flux = data[:,:,0]
        self.ivar = data[:,:,1]**2
        self.nobs = np.shape(data)[0]
        self.nstars = np.shape(data)[1]

    def mjd(self):
        return np.array([obs['mjd_g'] for obs in self.observations])

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
        conv = {'magzero': (np.arange(self.data.nobs),no,no),
                'mag': (np.arange(self.data.nobs,self.data.nobs+self.data.nstars)
                    ,no,no)}

        # we sample in magnitudes (ie ~log(flux)) so we need to convert
        # to fluxes and back
        fluxcon  = lambda x: 10**(-x/2.5)
        ifluxcon = lambda x: -2.5*np.log10(x)
        conv['zero'] = (np.arange(self.data.nobs),fluxcon,ifluxcon)
        conv['flux'] = (np.arange(self.data.nobs,self.data.nobs+self.data.nstars),
                fluxcon, ifluxcon)

        # in the more complicated models, we have more parameters
        zero = self.data.nobs+self.data.nstars
        if self.model >= 1:
            conv['jitterabs2'] = (zero,np.exp,np.log)
            conv['jitterrel2'] = (zero+1,np.exp,np.log)
            conv['pvar']       = (zero+2,no,no)
            conv['sigvar2']    = (zero+3,np.exp,np.log)
            zero += 4
        if self.model >= 2:
            conv['pbad']       = (zero,no,no)
            conv['sigbad2']    = (zero+1,np.exp,np.log)
            zero += 2
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

# ===================== #
#  Likelihood Function  #
# ===================== #

# precalculate this:
# log(1/sqrt(2 pi))
lisqrt2pi = - 0.5*np.log(2.0*np.pi)
def _lnnormal(x,mu,var):
    return -0.5*(x-mu)**2/var - 0.5*np.log(var) + lisqrt2pi

def lnprob(p,data,model=2):
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
    params = PhotoModel(data,p)
    prior = lnprior(params)
    if np.isinf(prior):
        return -np.inf
    return prior + lnlike(params,data)

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
    # pbad,pvar should be in [0,1]
    if p.model >= 1 and not 0 < p.pvar < 1:
        return -np.inf
    if p.model >= 2 and not 0 < p.pbad < 1:
        return -np.inf
    
    # g-band magnitude prior
    modelmag = p.mag
    mag = p.data.magprior[:,0]
    err = 0.5 # MAGIC p.data.magerr[:,1]
    lnprior = -0.5*np.sum((modelmag-mag)**2/err+np.sum(np.log(err)))

    return lnprior

def lnlike(p,data):
    """
    NAME:
        lnlike
    PURPOSE:
        return the ln likelihood of the data given model p
    INPUT:
        p    - PhotoModel object
        data - PhotoData object
    OUTPUT:
        the ln-likelihood
    HISTORY:
        Created by Dan Foreman-Mackey on Jun 07, 2011
    """
    # ff is the model flux matrix
    ff = np.outer(p.zero,p.flux)
    df = data.flux
    
    if p.model is 0:
        # the most basic possible model assuming Gaussian uncertainty
        # check for NaNs propagated from the inverse variance
        if np.any(np.isnan(data.ivar)):
            raise Exception("got some bad ivars")
        like = ((ff-df)**2*data.ivar).flatten()
        return -0.5*np.sum(like)
    if p.model >= 1:
        # allow for stars to vary
        # some data.ivar will be zero neglect those
        inds = data.ivar > 0.0
        sig2 = 1.0/(data.ivar[inds])
        
        ff = ff[inds]
        df = df[inds] 
        
        # calculate likelihood
        delta2 = p.jitterabs2 + p.jitterrel2*ff**2
        pow1 = np.log(1.0-p.pvar) \
                + _lnnormal(df,ff,sig2+delta2)
        pow2 = np.log(p.pvar) \
                + _lnnormal(df,ff,
                        sig2+delta2+p.sigvar2)

        lnlike1 = np.zeros(np.shape(data.flux))
        lnlike1[inds] = np.logaddexp(pow1,pow2)
        
        if p.model is 1:
            return np.sum(lnlike1)

        # some observations are bad too
        badlike = np.ones(np.shape(data.flux))
        badlike[inds] = _lnnormal(df,ff,sig2+delta2+p.sigbad2)
        pow1 = np.log(1-p.pbad) + np.sum(lnlike1,axis=-1)
        pow2 = np.log(p.pbad) + np.sum(badlike,axis=-1)
        lnlike2 = np.logaddexp(pow1,pow2)

        return np.sum(lnlike2)

def lnprob_badobs(p,data):
    """
    NAME:
        lnprob_badobs
    PURPOSE:
        returns the vector of probabilities that each observation is bad
    INPUT:
    OUTPUT:
        ln-probability that the observation is bad
    HISTORY:
        Created by Daniel Foreman-Mackey on 2011-04-03
    """
    if p.model < 2:
        print "lnprob_bad obs only works for model 2"
        return 0

    # choose only the data where inverse variance is non-zero
    inds = data.ivar > 0
    sig2 = 1.0/(data.ivar[inds])
    
    ff = np.outer(p.zero,p.flux)[inds]
    df = data.flux[inds]
    
    # calculate likelihood1
    delta2 = p.jitterabs2 + p.jitterrel2*ff**2
    pow1 = np.log(1.0-p.pvar) \
            + _lnnormal(df,ff,sig2+delta2)
    pow2 = np.log(p.pvar) \
            + _lnnormal(df,ff,
                    sig2+delta2+p.sigvar2)

    lnlike1 = np.zeros(np.shape(data.flux))
    lnlike1[inds] = np.logaddexp(pow1,pow2)

    badlike = np.ones(np.shape(data.flux))
    badlike[inds] = _lnnormal(df,ff,sig2+delta2+p.sigbad2)
    pow1 = np.log(1-p.pbad) + np.sum(lnlike1,axis=-1)
    pow2 = np.log(p.pbad) + np.sum(badlike,axis=-1)

    return pow2-pow1

def lnprob_variable(p,data,infield=False):
    """
    NAME:
        lnprob_variable
    PURPOSE:
        returns the probability that a given star is variable
    INPUT:
        starid  - integer id number of the star
    OPTIONAL:
        infield - if True this this will return a list (len nobs) with the
                  odds ratio that a star will be variable in a given field
    OUTPUT:
        ln-probability that the star is variable
    HISTORY:
        Created by Daniel Foreman-Mackey on 2011-04-03
    """
    if p.model < 2:
        print "lnprob_variable only works for model 2"
        return 0

    # choose only the data where inverse variance is non-zero
    inds = data.ivar > 0
    sig2 = 1.0/(data.ivar[inds])
    
    ff = np.outer(p.zero,p.flux)[inds]
    df = data.flux[inds]
    
    # calculate likelihood1
    delta2 = p.jitterabs2 + p.jitterrel2*ff**2
    pow1 = np.log(1.0-p.pvar) \
            + _lnnormal(df,ff,sig2+delta2)
    pow2 = np.log(p.pvar) \
            + _lnnormal(df,ff,
                    sig2+delta2+p.sigvar2)

    ret = np.zeros(np.shape(data.flux))
    ret[inds] = pow2-pow1
    
    # odds ratio
    if infield:
        return np.sum(ret,axis=-1)
    return np.sum(ret,axis=0)
    
