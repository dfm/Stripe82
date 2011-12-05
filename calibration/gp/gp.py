#!/usr/bin/env python
# encoding: utf-8
"""
My Gaussian process object

TODO: use the gradient for optimization

"""

from __future__ import division

__all__ = ['GaussianProcess']

import numpy as np

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
import scipy.linalg as linalg

class GaussianProcess(object):
    def __init__(self, s2=1.0, a2=1.0, la2=1.0):
        self.s2  = s2
        self.a2  = a2
        self.la2 = la2
        self.L,self.alpha = None,None

    def __repr__(self):
        return "GaussianProcess(s2=%f,a2=%f,la2=%f,b2=%f,lb2=%f)"\
                %(self.s2, self.a2, self.la2, self.b2, self.lb2)

    def k(self, x, y):
        p = (x[:, None] - y[None, :])**2/self.la2
        k = self.a2*np.exp(-0.5*p)
        return k, p

    def K(self, x, y):
        return self.s2*np.identity(len(y))+self.k(x,y)[0]

    def dK(self, x):
        k, p = self.k(x, x)
        dS = 2*np.sqrt(self.s2)*np.identity(len(self.y))
        da = 2*k/np.sqrt(self.a2)
        dl = k*p/np.sqrt(self.la2)
        return self.s2*np.identity(len(self.y))+k, [dS, da, dl]

    def condition(self, x, y):
        self.mu_x, self.std_x = np.mean(x, axis=0), np.std(x, axis=0)
        self.mu_y, self.std_y = np.mean(y, axis=0), np.std(y, axis=0)
        self.x, self.y = (x-self.mu_x)/self.std_x, (y-self.mu_y)/self.std_y

        self.Kxx, self.dKxx = self.dK(self.x)
        self.L = linalg.cho_factor(self.Kxx)
        self.alpha = linalg.cho_solve(self.L, self.y)
        self.denom = self.y.size*np.log(2*np.pi)

    def predict(self, x0):
        x0 = (x0-self.mu_x)/self.std_x

        # mean
        k0 = self.k(self.x,x0)[0]
        mean = np.dot(self.alpha, k0)

        # cov
        v = np.dot(linalg.cho_solve(self.L, k0).T, k0)
        cov = self.k(x0, x0)[0] - v

        return mean, cov

    def sample(self, x, N=1):
        mean,cov = self.predict(x)
        return mean*self.std_y+self.mu_y, \
                np.random.multivariate_normal(mean,cov,N)*self.std_y+self.mu_y

    def marginal_loglike(self):
        logp = -0.5*np.dot(self.y, self.alpha) - np.trace(self.L[0]) - self.denom
        return logp

    def grad_loglike(self):
        # gradient
        grad = [np.trace(np.dot(np.outer(self.alpha, self.alpha), self.dKxx[i])
            -linalg.cho_solve(self.L, self.dKxx[i])) for i in range(3)]
        return 0.5*np.array(grad)

    def train(self, x, y):
        def chi2(p):
            self.s2, self.a2, self.la2 = p**2
            self.condition(x,y)
            lnlike = self.marginal_loglike()
            print p**2, lnlike
            return -lnlike

        def dchi2(p):
            self.s2, self.a2, self.la2 = p**2
            self.condition(x,y)
            dlnlike = self.grad_loglike()
            # print dp, np.sum(np.abs(dp))
            return -dlnlike

        p0 = np.sqrt([self.s2, self.a2, self.la2])
        print op.fmin_bfgs(chi2,p0,fprime=dchi2,maxiter=10)**2

if __name__ == '__main__':
    # import matplotlib
    # matplotlib.use('Agg')
    import matplotlib.pyplot as pl
    np.random.seed(5)

    N = 400
    s = 10
    x = np.linspace(0,s,N)#s*np.random.rand(N)
    y = np.sin(x) + 0.1*np.random.randn(N)

    p = GaussianProcess(a2=0.05, la2=0.05)
    p.s2 = 0.1
    p.train(x,y)

    # plot fit
    x0 = np.linspace(min(x),max(x),100)
    y0, samples = p.sample(x0, N=500)
    pl.plot(x0,samples.T,'k',alpha=0.1)
    pl.plot(x0,y0,'b')

    #plot data
    pl.plot(x,y,'.r')

    pl.show()

    # pl.figure()
    # pl.imshow(cov)

    # pl.savefig("cov-%d.png"%s)

