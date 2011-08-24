#!/usr/bin/env python
# encoding: utf-8
"""
Options for the stripe82 module

History
-------
2011-08-16 - Created by Dan Foreman-Mackey

"""

__all__ = ['iau_name']

import numpy as np

# RA in degrees
def ra2hms(ra):
    while ra < 0.:
        ra += 360.
    while ra > 360.:
        ra -= 360.
    h = ra * 24. / 360.
    hh = int(np.floor(h))
    m = (h - hh) * 60.
    mm = int(np.floor(m))
    s = (m - mm) * 60.
    return (hh, mm, s)

# Dec in degrees
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
    # don't just return sgn*d because values between 0 and 1 deg will get you!
    return (sgn, d, m, s)

def iau_name(prefix,ra, dec):
    (rh,rm,rs) = ra2hms(ra)
    (sgn,dd,dm,ds) = dec2dms(dec)
    # According to http://www.sdss.org/dr3/coverage/IAU.html
    # the coordinates are truncated, not rounded.
    rcs = int(rs * 100.)
    dds = int(ds * 10.)

    return '%s~J%02i%02i%02i.%02i%s%02i%02i%02i.%01i' % (
        prefix,
        rh, rm, rcs / 100, rcs % 100,
        '+' if sgn >= 0 else '-',
        dd, dm, dds / 10, dds % 10)


