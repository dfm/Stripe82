#!/usr/bin/env python
# encoding: utf-8
"""
My Gaussian process object

History
-------
2011-09-07 - Created by Dan Foreman-Mackey

"""

from __future__ import division

__all__ = ['GaussianProcess']

import numpy as np

import scipy.sparse as sp
import scipy.sparse.linalg
from scipy.linalg import det
import scipy.optimize as op

from _gp import _sparse_k

class GaussianProcess(object):
    """
    A Gaussian process object

    Parameters
    ----------
    kernel : function
        The kernel to use

    History
    -------
    2011-09-07 - Created by Dan Foreman-Mackey

    """
    def __init__(self, s2=1.0, a2=1.0, la2=1.0, b2=1.0, lb2=1.0):
        self._s2  = s2
        self._a2  = a2
        self._b2  = b2
        self._la2 = la2
        self._lb2 = lb2
        self._L,self._alpha = None,None

    def __repr__(self):
        return "GaussianProcess(s2=%f,a2=%f,la2=%f,b2=%f,lb2=%f)"\
                %(self._s2, self._a2, self._la2, self._b2, self._lb2)

    def K(self,x,y=None):
        if y is None:
            y = x
        b = _sparse_k(self._a2, self._la2, self._b2, self._lb2, x, y) \
            + self._s2 * sp.identity(len(x),format="csc")
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
        ks = _sparse_k(self._a2, self._la2, self._b2, self._lb2, x, self._x)
        # calculate the mean
        f = ks.dot(np.atleast_2d(self._alpha).T)[:,0]
        if not cov:
            return f
        kss = _sparse_k(self._a2, self._la2, self._b2, self._lb2, x, x)
        kik = np.zeros(ks.shape)
        for i in xrange(ks.shape[0]):
            kik[i,:] = self._L.solve(np.array(ks[i,:].todense())[0])
        v = kss - ks.dot(kik.T)
        return f,v

    def optimize(self,x,y):
        def chi2(p):
            self._a2 = p[:1]**2
            self._b2, self._lb2 = p[1:-1]**2
            if self._lb2 < self._la2:
                return np.inf
            self._s2 = np.exp(p[-1])
            self.fit(x,y)
            s,lndet = np.linalg.slogdet(self._Kxx.todense())
            detK = lndet + len(y)*np.log(2*np.pi)
            c2 = np.dot(y,self._alpha) + detK
            return c2

        p0 = np.sqrt([self._a2])
        p0 = np.append(p0, np.sqrt([self._b2, self._lb2]))
        p0 = np.append(p0, np.log(self._s2))
        op.fmin(chi2,p0,disp=0)

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

    N = 100
    s = 10
    x = np.linspace(0,s,N)#s*np.random.rand(N)
    y = np.sin(x) + 0.01*np.random.randn(N)

    p = GaussianProcess(la2=0.25)
    p._s2 = 0.001
    p.optimize(x,y)
    p.fit(x,y)
    print p

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

