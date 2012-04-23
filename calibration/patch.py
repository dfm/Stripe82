"""
The model of a single calibration patch.

"""

__all__ = ["Patch"]

import numpy as np

import scipy.optimize as op

class Patch(object):
    """
    The patch model.

    ## Arguments

    * `fobs` (numpy.ndarray): The list of observations. This object should
      have the shape `(nruns, nstars)`.
    * `ivar` (numpy.ndarray): The observational inverse variance. This should
      be the same shape as `fobs`.

    """
    def __init__(self, fobs, ivar):
        self.fobs = np.atleast_2d(fobs)
        self.ivar = np.atleast_2d(ivar)
        self._mask = ivar > 0.0
        self.var = np.zeros_like(self.ivar)
        self.var[self._mask] = 1./self.ivar[self._mask]

        self.M, self.N = self.fobs.shape

    @property
    def nruns(self):
        return self.M

    @property
    def nstars(self):
        return self.N

    def get_initial_params(self, fs, maxiter=20, tol=1.24e-4):
        """
        Given an initial guess for the stellar fluxes, compute a reasonable
        first guess for all the other parameters.

        ## Arguments

        * `fs` (numpy.ndarray): A initial guess for the fluxes of the stars.

        ## Keyword Arguments

        * `maxiter` (int): The maximum number of iterations to run when
          estimating the parameters.
        * `tol` (float): The convergence tolerance constraint on the stellar
          flux estimate.

        """
        # Iteratively determine a first guess for the fluxes and zero points.
        f0 = np.mean(self.fobs/fs[None, :], axis=1)
        for i in xrange(maxiter):
            fs2 = np.mean(self.fobs/f0[:, None], axis=0)
            f0 = np.mean(self.fobs/fs2[None, :], axis=1)
            d = np.sum(np.abs(fs-fs2))
            fs = fs2
            if d <= tol:
                break
        w = np.log(np.std(self.fobs/f0[:, None], axis=0))
        s = np.log(np.std(self.fobs/fs[None, :], axis=1))
        return np.concatenate([np.log([2, 0.01]), np.log(fs), w, np.log(f0), s])

    def _preprocess(self, p):
        """
        Given a vector of input parameters, compute the quantities needed
        by `nll` and `grad_nll`.

        ## Arguments

        * `p` (list): The input parameters.

        """
        M, N = self.M, self.N

        # First, parse the input parameters. There's a lot of magic here.
        d2, e2 = np.exp(2*p[:2])
        if not np.isfinite(d2):
            d2 = 0.0
        if not np.isfinite(e2):
            e2 = 0.0

        fs = np.exp(p[2:2+N])
        w2 = np.exp(2*p[2+N:2*(1+N)])
        w2[~np.isfinite(w2)] = 0
        f0 = np.exp(p[2*(1+N):2*(1+N)+M])
        s2 = np.exp(2*p[2*(1+N)+M:2*(1+N+M)])
        s2[~np.isfinite(s2)] = 0

        # `Cs.shape = (nruns, nstars)`.
        Cs = np.outer(f0, fs)

        sig2 = np.array(self.var)
        sig2 += d2 + e2 * Cs**2
        sig2 += s2[:, None] + w2[None, :]
        ivar = np.zeros_like(Cs)
        ivar[self._mask] = 1./sig2[self._mask]

        return d2, e2, fs, w2, f0, s2, Cs, ivar

    def nll(self, p, fp, ivfp):
        """
        Calculate the negative log-likelihood of the parameters `p`.

        ## Arguments

        * `p` (list): The list of parameters for which the likelihood should
          be calculated. The layout of this object should be:
              [log(delta), log(eta), f_star, log(w), f_0, log(sigma)]
        * `fp` (numpy.ndarray): The mean value for the Gaussian prior on
          stellar flux. There should be one entry for each star.
        * `ivfp` (numpy.ndarray): The inverse variance of the stellar flux
          prior.

        """
        d2, e2, fs, w2, f0, s2, Cs, ivar  = self._preprocess(p)
        nll = 0.5*np.sum(((self.fobs - Cs)**2 * ivar)[self._mask]
                - np.log(ivar[self._mask]))
        return nll + 0.5*np.sum((fs - fp)**2 * ivfp)

    def grad_nll(self, p, fp, ivfp):
        """
        Calculate the analytic expression for the gradient of the negative
        log-likelihood function. The arguments are the same as for the `nll`
        function.

        """
        d2, e2, fs, w2, f0, s2, Cs, ivar  = self._preprocess(p)

        delta = self.fobs - Cs

        dlds2 = ivar-(delta*ivar)**2 # 2 * dl/dsig2
        dldC  = -delta*ivar + dlds2 * e2 * Cs

        dldf0 = np.sum(dldC * fs[None, :], axis=1)
        dldfs = np.sum(dldC * f0[:, None], axis=0) + (fs-fp) * ivfp

        dldw = np.sum(dlds2, axis=0) * w2
        dlds = np.sum(dlds2, axis=1) * s2

        dldd = np.sum(dlds2) * d2
        dlde = np.sum(dlds2 * Cs**2) * e2

        grad = np.concatenate([[dldd, dlde], dldfs, dldw, dldf0, dlds])
        return grad

    def optimize(self, fp, ivfp, **kwargs):
        """


        """
        p0 = self.get_initial_params(fp)
        p1 = op.fmin_bfgs(self.nll, p0,
                args=(fp, ivfp), disp=0, **kwargs)
        # fprime=self.grad_nll,
        self.d2, self.e2, self.fs, self.w2, self.f0, self.s2, Cs, ivar\
                = self._preprocess(p1)
        return p1

    @classmethod
    def test_grad(cls, d=1e-6):
        """
        A unit test to make sure that the gradient calculation is right.

        ## Keyword Arguments

        * `d` (float): The step size to use to calculate the numerical
          gradient vector. This also sets the tolerance of the test.

        """
        # Generate some fake data.
        nstars, nobs = 10, 50

        # Synthetic stellar fluxes...
        fp = 100+5*np.random.randn(nstars)
        ivfp = 25./np.ones_like(fp)

        # ...and zero points.
        f0 = 50+100*np.random.rand(nobs)
        fo = np.outer(f0, fp)

        # The observational variance.
        var = 10*np.random.rand(nobs*nstars).reshape(nobs, nstars)
        fo += np.sqrt(var)*np.random.randn(nobs*nstars).reshape(nobs, nstars)
        ivar = 1.0/var

        # Set up the patch object.
        patch = cls(fo, ivar)
        p0 = patch.get_initial_params(fp)

        # What is the analytic gradient?
        grad0 = patch.grad_nll(p0, fp, ivfp)

        # Iterate through the components and use the centered finite
        # difference to calculate the numerical gradient.
        grad1 = np.zeros_like(grad0)
        for i in range(len(p0)):
            p0[i] += d
            nllp = patch.nll(p0, fp, ivfp)
            p0[i] -= 2*d
            nllm = patch.nll(p0, fp, ivfp)
            p0[i] += d
            grad1[i] = (nllp-nllm)/(2.*d)

        assert np.all(np.abs(grad0-grad1)/np.abs(patch.nll(p0,fp,ivfp)) < d)

class SyntheticPatchData(object):
    """
    Generate synthetic data from a more complex model of what the data should
    look like in a patch.

    ## Arguments

    * `nstars` and `nobs` (int): The number of stars/runs.

    ## Keyword Arguments

    * `Qvar` (float): The probability that a star is variable.
    * `Qbad` (float): The probability that a run is bad.

    """
    def __init__(self, nstars, nobs, Qvar=0.1, Qbad=0.1):
        self.Svar = 0.5
        self.Sbad = 600.
        self.delta = 1.
        self.eta = 0.001

        # Synthetic stellar fluxes...
        self.fp = 100+5*np.random.randn(nstars)
        self.ivfp = 5**-2 * np.ones(nstars)
        fs = self.fp + np.sqrt(1.0/self.ivfp) * np.random.randn(nstars)

        # ...and zero points.
        self.f0 = 50+10*np.random.rand(nobs)

        # Which stars are variable?
        self.isvar = np.arange(nstars)[np.random.rand(nstars) < Qvar]
        self.isbad = np.arange(nobs)[np.random.rand(nobs) < Qbad]

        # Construct the data.
        fo = np.outer(self.f0, fs)

        # Add jitter.
        fo += self.eta*fo*np.random.randn(nobs*nstars).reshape(nobs, nstars)
        fo += self.delta*np.random.randn(nobs*nstars).reshape(nobs, nstars)

        # Add variable stars and bad observations.
        for alpha in self.isvar:
            fo[:, alpha] += self.fp[alpha] * self.Svar * np.random.randn(nobs)
        for i in self.isbad:
            fo[i, :] += self.Sbad * np.random.randn(nstars)

        # The observational variance.
        var = np.random.rand(nobs*nstars).reshape(nobs, nstars)
        fo += np.sqrt(var)*np.random.randn(nobs*nstars).reshape(nobs, nstars)

        # Save the dataset.
        self.ivar = 1.0/var
        self.fobs = fo

        # Randomly censor some data.
        self.ivar[np.random.randint(nobs, size=5),
                np.random.randint(nstars, size=5)] = 0

if __name__ == "__main__":
    Patch.test_grad()

    data = SyntheticPatchData(10, 50)
    patch = Patch(data.fobs, data.ivar)
    p = patch.optimize(data.fp, data.ivfp)

    print np.abs(data.fp - patch.fs)/data.fp

    print patch.w2[data.isvar]
    print patch.w2

    print patch.s2[data.isbad]
    print patch.s2

    import matplotlib.pyplot as pl

    [pl.plot(data.fobs[:, i], ".") for i in range(10)]

    pl.figure()
    [pl.plot(data.fobs[:, i]/p[2*(1+10):2*(1+10)+50], ".") for i in range(10)]
    [pl.gca().axhline(p[2+i]) for i in range(10)]
    pl.show()

