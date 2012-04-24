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
        ivar = np.atleast_2d(ivar)
        self._mask = ivar > 0.0
        self.var = np.zeros_like(ivar)
        self.var[self._mask] = 1./ivar[self._mask]

        self.N, self.M = self.fobs.shape

    @property
    def nruns(self):
        return self.N

    @property
    def nstars(self):
        return self.M

    def get_initial_params(self, fs, maxiter=100, tol=1.24e-3):
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
        # f0 = np.mean(self.fobs/fs[None, :], axis=1)
        ivar = np.zeros_like(self.var)
        ivar[self._mask] = 1./self.var[self._mask]
        iv0 = np.sum(ivar, axis=0)
        iv1 = np.sum(ivar, axis=1)

        f0 = np.abs(np.sum(ivar*self.fobs/fs[None,:], axis=1)/iv1)
        # for i in xrange(maxiter):
        #     fs2 = np.sum(ivar*self.fobs/f0[:, None], axis=0)/iv0
        #     f0 = np.sum(ivar*self.fobs/fs2[None, :], axis=1)/iv1
        #     d = np.sum(np.abs(fs-fs2))
        #     fs = fs2
        #     if d <= tol:
        #         break
        p = np.concatenate([fs, f0])
        return p

    def _preprocess(self, p, hyper):
        """
        Given a vector of input parameters, compute the quantities needed
        by `nll` and `grad_nll`.

        ## Arguments

        * `p` (list): The input parameters.

        """
        p = np.atleast_1d(p)

        M, N = self.M, self.N

        # First, parse the input parameters. There's a lot of magic here.
        self.fs = p[:M]
        self.f0 = p[M:M+N]
        if hyper:
            self.d2 = p[M+N:M+2*N]
            self.b2 = p[M+2*N:M+3*N]
            self.e2 = p[M+3*N:]
        else:
            self.d2 = np.zeros(N)
            self.b2 = np.zeros(N)
            self.e2 = np.zeros(M)

        # `Cs.shape = (nruns, nstars)`.
        self.Cs   = np.outer(self.f0, self.fs)
        self.rel2 = self.b2[:,None] + self.e2[None,:]

        sig2 = np.array(self.var)
        sig2 += self.d2[:,None] + self.rel2 * self.Cs**2
        self.ivar = np.zeros_like(self.Cs)
        self.ivar[self._mask] = 1./sig2[self._mask]

    def nll(self, p, fp, ivfp, hyper):
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
        self._preprocess(p, hyper)
        nll = 0.5*np.sum(((self.fobs - self.Cs)**2 * self.ivar)[self._mask]
                - np.log(self.ivar[self._mask]))
        nll += 0.5*np.sum((self.fs - fp)**2 * ivfp)
        return nll

    def grad_nll(self, p, fp, ivfp, hyper):
        """
        Calculate the analytic expression for the gradient of the negative
        log-likelihood function. The arguments are the same as for the `nll`
        function.

        """
        self._preprocess(p, hyper)
        Cs2 = self.Cs**2

        delta = self.fobs - self.Cs

        # Cache some gradients for speed.
        dlds2 = 0.5*(self.ivar-(delta*self.ivar)**2)
        dldC  = -delta*self.ivar + 2 * dlds2 * self.rel2 * self.Cs

        # Gradients with respect to the parameters of interest.
        dldf0 = np.sum(dldC * self.fs[None, :], axis=1)
        dldfs = np.sum(dldC * self.f0[:, None], axis=0) + (self.fs-fp) * ivfp

        if not hyper:
            return np.concatenate([dldfs, dldf0])

        # Gradients with respect to the variances.
        dldd = np.sum(dlds2, axis=1)
        Cs2dlds2 = Cs2 * dlds2
        dldb = np.sum(Cs2dlds2, axis=1)
        dlde = np.sum(Cs2dlds2, axis=0)

        grad = np.concatenate([dldfs, dldf0, dldd, dldb, dlde])
        return grad

    def optimize(self, fp, ivfp, **kwargs):
        """


        """
        p0 = self.get_initial_params(fp)
        res = op.fmin_l_bfgs_b(self.nll, p0, fprime=self.grad_nll,
                args=(fp, ivfp, False), disp=0,
                bounds=[(0,None) for i in range(len(p0))], **kwargs)
        p1 = np.concatenate([res[0], 0.01*np.ones(2*self.N+self.M)])
        res = op.fmin_l_bfgs_b(self.nll, p1, fprime=self.grad_nll,
                args=(fp, ivfp, True), disp=0,
                bounds=[(0,None) for i in range(len(p1))], **kwargs)
        p2 = res[0]
        self._preprocess(p2, True)
        return p2

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
        p0 = np.concatenate([patch.get_initial_params(fp),
                                        0.01*np.ones(2*nobs+nstars)])

        # What is the analytic gradient?
        grad0 = patch.grad_nll(p0, fp, ivfp, True)

        # Iterate through the components and use the centered finite
        # difference to calculate the numerical gradient.
        grad1 = np.zeros_like(grad0)
        for i in range(len(p0)):
            p0[i] += d
            nllp = patch.nll(p0, fp, ivfp, True)
            p0[i] -= 2*d
            nllm = patch.nll(p0, fp, ivfp, True)
            p0[i] += d
            grad1[i] = (nllp-nllm)/(2.*d)

        assert np.all(np.abs(grad0-grad1)/np.abs(patch.nll(p0,fp,ivfp,True))
                            < d)

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
    # Patch.test_grad()

    # M, N = 5, 10

    # data = SyntheticPatchData(M, N)
    # patch = Patch(data.fobs, data.ivar)
    # p = patch.optimize(data.fp, data.ivfp)

    import h5py
    from conversions import *

    f = h5py.File("test_data.h5")
    flux = f["flux"][...]
    ivar = f["ivar"][...]
    prior = f["prior"][...]
    ivp = np.ones_like(prior)

    patch = Patch(flux, ivar)
    # patch._preprocess(patch.get_initial_params(prior))
    patch.optimize(prior, ivp)

    import matplotlib.pyplot as pl

    cs = np.zeros((patch.nruns, 4))
    cs[:, -1] = 1-0.9*patch.b2/np.max(patch.b2)

    for i in np.argsort(patch.fs):
        pl.clf()
        pl.errorbar(np.arange(patch.nruns), flux[:,i]/patch.f0,
                yerr=np.sqrt(patch.var[:,i])/patch.f0, ls="None",
                marker="None", zorder=1, barsabove=False, color="k")
        pl.scatter(np.arange(patch.nruns), flux[:, i]/patch.f0,
                c=cs, zorder=2, s=40)
        pl.gca().axhline(patch.fs[i], color="k")
        ymin = min(pl.gca().get_ylim()[0], 0)
        pl.ylim(ymin, 2*patch.fs[i]-ymin)
        pl.xlim(0, patch.nruns)
        pl.savefig("lc/%d.png")

