# encoding: utf-8
"""
Some RR Lyrae lightcurve analysis tools.

"""

__all__ = ["lc_model", "fit", "chi2", "find_period", "get_model"]

from multiprocessing import Pool
import numpy as np
from numpy.linalg import lstsq
import scipy.optimize as op

def lc_model(omega, amplitudes, order):
    a = amplitudes
    return lambda t: a[0] + np.sum([a[2*i+1] * np.sin((i+1)*omega*t)+\
                                    a[2*i+2] * np.cos((i+1)*omega*t)
                                            for i in range(order)],axis=0)

def fit(omega, time, flux, order=3):
    a = np.ones((len(time), 1+2*order))
    ip1 = np.arange(1, order+1)
    a[:,1::2] = np.sin(ip1[None,:] * time[:,None] * omega)
    a[:,2::2] = np.cos(ip1[None,:] * time[:,None] * omega)
    amplitudes = lstsq(a, flux)[0]
    return lc_model(omega, amplitudes, order), amplitudes

def get_model(period, time, flux, order=3):
    r = fit(2*np.pi/period, time, flux, order=order)
    return r

def chi2(model, time, flux, ferr=None):
    if ferr is None:
        ferr = np.ones_like(flux)
    return np.sum( (flux - model(time))**2 / ferr**2)

class _fit_wrapper(object):
    def __init__(self, time, flux, ferr, order):
        self.time  = time
        self.flux  = flux
        self.ferr  = ferr
        self.order = order

    def __call__(self, omega):
        # Check to make sure that the period isn't _very_ close to a day.
        if np.abs(omega - 2*np.pi) < 0.05:
            return 1e10

        # Do the fit.
        r = 0.0
        for k in self.time:
            model, a = fit(omega, self.time[k], self.flux[k],
                    order=self.order)

            # Calculate the amplitudes and make sure that the 1st order
            # dominates.
            a = a[1::2]**2 + a[2::2]**2
            if np.any(a[0] < a[1:]):
                return 1e10

            r += chi2(model, self.time[k], self.flux[k], ferr=self.ferr[k])

        return r

class _op_wrapper(object):
    def __init__(self, time, flux, ferr, order):
        self.time  = time
        self.flux  = flux
        self.ferr  = ferr
        self.order = order
        self.wrapped = _fit_wrapper(time, flux, ferr, order)

    def __call__(self, omega):
        res = op.fmin(self.wrapped, omega, disp=False, full_output=True)
        return res[:2]

def find_period(time, flux, ferr=None, order=3, N=30, Ts=[0.2, 1.3],
        pool=None):
    """
    Find the best fit period of an RR Lyrae lightcurve by doing a grid
    search then a non-linear refinement step.

    """
    # Set up a grid to do a grid search in frequency space.
    domega = 0.2/min([time[k].max()-time[k].min() for k in time])
    omegas = 2 * np.pi * np.arange(1./max(Ts), 1./min(Ts), domega)

    # Do a parallel grid search.
    if pool is None:
        pool = Pool()
    chi2 = pool.map(_fit_wrapper(time, flux, ferr, order), omegas)

    # Sort the results by chi2.
    inds = np.argsort(chi2)

    # Refine the top `N` best fits.
    ref = pool.map(_op_wrapper(time, flux, ferr, order), omegas[inds[:N]])

    # Clean up... otherwise we get the error: _Too many open files_.
    pool.close()
    pool.join()
    del pool

    # Sort the refinements and return the best.
    omega = min(ref, key=lambda x: x[1])
    return 2*np.pi/float(omega[0])

