#!/usr/bin/env python
# encoding: utf-8
"""
A general probabilistic calibration model

History
-------
2011-06-14 - Created by Dan Foreman-Mackey

"""

__all__ = ['PhotoData','PhotoModel','lnprob','lnprior','lnlike',
           'lnprob_badobs','lnprob_variable']

import numpy as np

from pymongo.son_manipulator import SONManipulator

class PhotoDataManipulator:


class PhotoData:
    """
    NAME:
        PhotoData
    PURPOSE:
        wrap the data into a consistent class
    INPUT:
        data - array from batch_photometry (nobs,nstars,:)
    HISTORY:
        Created by Dan Foreman-Mackey on Jun 07, 2011
    """
    def __init__(self,data,observations,stars):
        self.data = data
        self.stars = []
        for sid in stars:
            self.stars.append(survey.get_star(sid))
        self.observations = []
        for oid in observations:
            self.observations.append(survey.get_observation(oid))
        self.dasmag = np.array([[s['g'],s['Err_g']**2] for s in self.stars])
        self.flux = data[:,:,0]
        self.ivar = data[:,:,1]**2
        self.nobs = np.shape(data)[0]
        self.nstars = np.shape(data)[1]

    def dump(self,fn):
        """
        NAME:
            dump
        PURPOSE:
            save the contents of the class to fn
        INPUT:
            fn - where to save
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 08, 2011
        """
        f = h5py.File(fn,'w')
        f['data']   = self.data
        f['stars']  = self.stars
        f['fields'] = self.fields
        f.close()

class PhotoModel:
    """
    NAME:
        PhotoModel
    PURPOSE:
        interface for model parameters
    INPUT:
        fields - observations
        stars  - stars
        model  - (0,1,2) how good is the model? (default=0)
    HISTORY:
        Created by Dan Foreman-Mackey on Jun 07, 2011
    """
    def __init__(self,*args,**kwargs):
        if 'vector' in kwargs:
            vector = kwargs['vector']
        else:
            vector = None
        if len(args) is 2:
            fields,stars = args
            if 'model' in kwargs:
                model = kwargs['model']
            else:
                model = 0
        else:
            f = h5py.File(args[0])
            fields = f['fields'][...]
            stars  = f['stars'][...]
            vector = f['vector'][...]
            model  = f['model'][...]
            f.close()
        if model not in [0,1,2]:
            raise Exception("%d is not a valid PhotoModel id"%model)
        self.dasmag = np.array([[s[2],s[3]**2] for s in stars])
        self.fields = fields
        self.stars  = stars
        self.nobs   = len(fields)
        self.nstars = len(stars)
        self.model  = model
        self.conv,self.npars = self.param_names()
        if vector is not None:
            self.from_vector(vector)
    
    def param_names(self):
        """
        NAME:
            param_names
        PURPOSE:
            return the conversion conventions for the model parameters
        OUTPUT:
            
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 07, 2011
        """
        no = lambda x: x
        conv = {'magzero': (np.arange(self.nobs),no,no),
                'mag': (np.arange(self.nobs,self.nobs+self.nstars),no,no)}

        # we sample in magnitudes (ie log(flux)) so we need to convert
        # to fluxes and back
        fluxcon  = lambda x: 10**(-x/2.5)
        ifluxcon = lambda x: -2.5*np.log10(x)
        conv['zero'] = (np.arange(self.nobs),fluxcon,ifluxcon)
        conv['flux'] = (np.arange(self.nobs,self.nobs+self.nstars),
                fluxcon, ifluxcon)

        # in the more complicated models, we have more parameters
        zero = self.nobs+self.nstars
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
        NAME:
            from_vector
        PURPOSE:
            return a PhotoModel class given a vector of values
        INPUT:
            
        OUTPUT:
            
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 07, 2011
        """
        p0 = np.array(p0)
        for k in self.conv:
            setattr(self,k,self.conv[k][1](p0[self.conv[k][0]]))

    def vector(self):
        """
        NAME:
            vector
        PURPOSE:
            return a PhotoModel class given a vector of values
        INPUT:
            
        OUTPUT:
            
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 07, 2011
        """
        vec = np.zeros(self.npars)
        for k in self.conv:
            vec[self.conv[k][0]] = self.conv[k][2](getattr(self,k))
        return vec

    def dump(self,fn):
        """
        NAME:
            dump
        PURPOSE:
            save the contents of the class to fn
        INPUT:
            fn - where to save
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 08, 2011
        """
        f = h5py.File(fn,'w')
        f['vector'] = self.vector()
        f['fields'] = self.fields
        f['stars']  = self.stars
        f['model']  = self.model
        f.close()

# ===================== #
#  Likelihood Function  #
# ===================== #

# precalculate this:
# log(1/sqrt(2 pi))
lisqrt2pi = - 0.5*np.log(2.0*np.pi)
def _lnnormal(x,mu,var):
    return -0.5*(x-mu)**2/var - 0.5*np.log(var) + lisqrt2pi

def lnprob(p,data,model=0):
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
    params = PhotoModel(data.fields,data.stars,model=model,vector=p)
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
    dasmag = p.dasmag[:,0]
    daserr = 0.5 # MAGIC p.dasmag[:,1]
    lnprior = -0.5*np.sum((modelmag-dasmag)**2/daserr+np.sum(np.log(daserr)))

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
    



if __name__ == '__main__':
    main()


