"""
Unit conversions for mags/fluxes/probabilities

"""

__all__ = ['nmgy2mag', 'mag2nmgy', 'odds2prob', 'prob2odds',
                'logodds2prob', 'prob2logodds']


import numpy as np


def nmgy2mag(nmgy):
    return 22.5 - 2.5 * np.log10(nmgy)


def mag2nmgy(mag):
    return 10 ** (-0.4 * (mag - 22.5))


def odds2prob(odds):
    return odds / (1. + odds)


def prob2odds(prob):
    return prob / (1. - prob)


def logodds2prob(logodds):
    logoddsflat = logodds.flatten()
    prob = odds2prob(np.exp(logoddsflat))
    # fix overflow
    isinfornan = np.isinf(prob) + np.isnan(prob)
    if np.any(isinfornan):
        prob[isinfornan] = 1.0
    return prob.reshape(np.shape(logodds))


def prob2logodds(prob):
    return np.log(prob2odds(prob))


# From:
#    trac.astrometry.net/browser/trunk/src/astrometry/util/starutil_numpy.py

def ra_normalize(ra):
    return np.mod(ra, 360.)


def ra2hms(ra):
    ra = ra_normalize(ra)
    h = ra * 24. / 360.
    hh = int(np.floor(h))
    m = (h - hh) * 60.
    mm = int(np.floor(m))
    s = (m - mm) * 60.
    return (hh, mm, s)


def dec2dms(dec):
    sgn = (dec >= 0) and 1. or -1.
    d = dec * sgn
    dd = int(np.floor(d))
    m = (d - dd) * 60.
    mm = int(np.floor(m))
    s = (m - mm) * 60.
    if s >= 60.:
        m += 1.
        s -= 60.
    # Don't just return sgn * d because values between 0 and 1 deg will get
    # you!
    return (sgn, d, m, s)
