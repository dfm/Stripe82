#!/usr/bin/env python
# encoding: utf-8
"""
My Gaussian process object

"""

from __future__ import division

__all__ = ['GaussianProcess']

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
from scipy.linalg import det
import scipy.optimize as op

class GaussianProcess(object):
    """
    A Gaussian process object

    Parameters
    ----------
    kernel : function
        The kernel to use

    """
    def __init__(self,s2=1.0,a=1.0,l2=1.0):
        self._s2 = s2
        self._a  = a
        self._l2 = l2
        self._L,self._alpha = None,None

    def __repr__(self):
        return "GaussianProcess(s2=%f,a=%f,l2=%f)"%(self._s2,self._a,self._l2)

    def k(self,x1,x2,chi2max=25.0):
        """
        The default kernel function

        Parameters
        ----------
        x1,x2 : numpy.ndarray
            Vectors of positions.

        chi2max : float, optional
            Set clipping for sparseness.

        Returns
        -------
        k : numpy.ndarray
            Covariance matrix between x1 and x2.

        Note
        ----
        This works well for small matrices but it is poorly implemented for larger
        matrices --- especially if they are actually sparse!

        """
        d = (x1-x2)**2/self._l2
        k = sp.lil_matrix(d.shape)
        k = self._a*np.exp(-0.5*d)
        k[d > chi2max] = 0.0
        return sp.lil_matrix(k).tocsc()

    def K(self,x,y=None, ep2=0.0):
        if y is None:
            y = x
        b = self.k(x[:,np.newaxis],y[np.newaxis,:]) \
            + (self._s2 + ep2) *sp.identity(len(x),format="csc")
        return b

    def fit(self,x,y):
        self._x = x
        self._y = y
        self._Kxx = self.K(x)
        self._detK = None
        self._L = sp.linalg.splu(self._Kxx)
        self._alpha = self._L.solve(y)

    def __call__(self,x,cov=False):
        assert self._L is not None
        ks = self.k(x[:,np.newaxis],self._x[np.newaxis,:])
        # calculate the mean
        f = ks.dot(np.atleast_2d(self._alpha).T)[:,0]
        if not cov:
            return f
        kss = self.k(x[:,np.newaxis],x[np.newaxis,:])
        kik = np.zeros(ks.shape)
        for i in xrange(ks.shape[0]):
            kik[i,:] = self._L.solve(np.array(ks[i,:].todense())[0])
        v = kss - ks.dot(kik.T)
        return f,v

    def optimize(self,x,y):
        def chi2(p):
            self._a,self._l2 = p[:-1]**2
            self._s2 = np.exp(p[-1])
            self.fit(x,y)
            detK = np.log(det(self._Kxx.todense())) \
                     + len(y)*np.log(2*np.pi)
            c2 = np.dot(y,self._alpha) + detK
            return c2

        p0 = np.sqrt([self._a,self._l2])
        p0 = np.append(p0, np.log(self._s2))
        p1 = op.fmin(chi2,p0,disp=0)
        # print p1,(self._a,self._l2,self._s2)

    def sample_prior(self,x):
        """
        Return N samples from the prior

        Parameters
        ----------
        x : numpy.ndarray
            Positions in parameter space

        Returns
        -------
        ret : type
            Description

        """
        return np.random.multivariate_normal(np.zeros(len(x)),self.K(x).todense())

    def sample(self,x,N=1):
        """
        Sample given some data

        """
        mean,cov = self(x,cov=True)
        return np.random.multivariate_normal(mean,cov,N)

if __name__ == '__main__':
    import pylab as pl
    np.random.seed(5)

    N = 50
    s = 10
    x = s*np.random.rand(N)
    y = np.sin(x) + 0.1*np.random.randn(N)

    p = GaussianProcess(l2=0.5)
    p._s2 = 0.001
    p.optimize(x,y)

    # plot fit
    x0 = np.linspace(min(x),max(x),500)
    y0 = p.sample(x0,100)
    pl.plot(x0,y0.T,'.k',alpha=0.1)

    mean,cov = p(x0,cov=True)
    pl.plot(x0,mean)

    #plot data
    pl.plot(x,y,'.r')

    pl.savefig("%d.png"%s)

    pl.figure()
    pl.imshow(cov)

    pl.savefig("cov-%d.png"%s)

