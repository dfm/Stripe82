#!/usr/bin/env python
# encoding: utf-8
"""
My Gaussian process object

TODO: use the gradient for optimization

"""

from __future__ import division

__all__ = ['GaussianProcess']

import numpy as np

import scipy.sparse as sp
# import scipy.sparse.linalg
import scipy.linalg as linalg

try:
    from numpy.linalg import slogdet
except ImportError:
    print "Warning: upgrade numpy to get slogdet!"
    def slogdet(A):
        n = A.shape[0]
        lu,piv = linalg.lu_factor(A)
        s = 1. - 2. * (np.sum(piv != np.arange(0, n )) % 2)
        d = lu.diagonal()
        absd = np.abs(d)
        s *= np.prod(np.sign(d))

        return s, np.sum(np.log(absd))

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
        K = _sparse_k(self._a2, self._la2, self._b2, self._lb2, x, y,
                s2=self._s2, grad=0)
        b = K + self._s2 * sp.identity(len(x),format="csc")
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
        x -= np.mean(x)
        x /= np.var(x)
        y -= np.mean(y)
        y /= np.var(y)
        def chi2(p):
            self._a2, self._b2, self._lb2, self._s2 = np.exp(p)
            # if self._lb2 < self._la2:
            #     return np.inf
            # self._s2 = np.exp(p[-1])
            self.fit(x,y)
            s,lndet = slogdet(self._Kxx.todense())
            detK = lndet + len(y)*np.log(2*np.pi)
            c2 = np.dot(y,self._alpha) + detK
            print c2
            return c2

        # def gradchi2(p):
        #     print "grad"
        #     dKdtheta = self._gradK
        #     gradc2 = []
        #     for i in range(4):
        #         dK = dKdtheta[i]
        #         gradc2.append(-0.5*np.trace(\
        #                 np.dot(np.outer(self._alpha, self._alpha), dK)\
        #                 -self._L.solve(dK.todense())))
        #     print gradc2
        #     return gradc2

        p0 = np.sqrt([self._a2, self._b2, self._lb2, self._s2])
        # p0 = np.append(p0, np.sqrt([self._b2, self._lb2]))
        # p0 = np.append(p0, np.log(self._s2))
        self.fit(x,y)
        op.fmin_bfgs(chi2,p0,disp=1)

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
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as pl
    np.random.seed(5)

    N = 15
    s = 10
    x = np.linspace(0,s,N)#s*np.random.rand(N)
    y = np.sin(x) + 0.05*np.random.randn(N)

    p = GaussianProcess(a2=0.0, la2=0.001)
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

    # pl.figure()
    # pl.imshow(cov)

    # pl.savefig("cov-%d.png"%s)

