#!/usr/bin/env python
# encoding: utf-8
"""


"""

__all__ = ['']

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

from models import Star, CalibRun
from conversions import nmgy2mag

def plot_colorcolor(run, camcol, ax=None):
    try:
        runs = CalibRun.find({'run': run, 'camcol': camcol})
    except:
        return None
    [r.fit_gp() for r in runs]
    bands = [r._band for r in runs]
    stars = Star.find_with_photometry_in_run(runs[0])

    rminusi = []
    gminusr = []
    rminusi0 = []
    gminusr0 = []

    for star in stars:
        photo = [getattr(star, "_photo_%s"%(r._band)) for r in runs]
        if np.all([str(runs[ri]._id) in photo[ri] for ri in range(len(runs))]):
            f0 = [r.sample_gp([star._ra], N=0)[0] for r in runs]
            fobs = [photo[ri][str(r._id)][0][1] for ri, r in enumerate(runs)]
            fstar = np.array([fobs[i]/f0[i] for i in range(len(f0))])
            if np.all(fstar > 0):
                g = nmgy2mag(fstar[bands.index('g')])
                r = nmgy2mag(fstar[bands.index('r')])
                i = nmgy2mag(fstar[bands.index('i')])
                if g < 20:
                    rminusi.append(r - i)
                    gminusr.append(g - r)
                    rminusi0.append(star.r - star.i)
                    gminusr0.append(star.g - star.r)

    print len(gminusr)
    pl.plot(gminusr0, rminusi0, '.g', ms=6, alpha=0.3)
    pl.plot(gminusr, rminusi, '.k', ms=5)
    if ax is not None:
        ax.plot(gminusr, rminusi, '.k', ms=4, alpha=0.15)
    return gminusr, rminusi

if __name__ == '__main__':
    ax = pl.figure().add_subplot(111)
    pl.figure()
    runs = [r for r in CalibRun._collection.find({"band": "g", "failed": {"$exists": False}})]
    for i in range(len(runs)):
        if "failed" not in runs[i]:
            pl.clf()
            if plot_colorcolor(runs[i]["run"],runs[i]["camcol"], ax=ax) is not None:
                pl.xlim(-1, 2)
                pl.ylim(-1, 2)
                pl.xlabel(r"$g-r$", fontsize=16)
                pl.ylabel(r"$r-i$", fontsize=16)
                pl.savefig("color-color2/%d-%d.png"%(runs[i]["run"],runs[i]["camcol"]))

                ax.set_xlim(-1, 2)
                ax.set_ylim(-1, 2)
                ax.set_xlabel(r"$g-r$", fontsize=16)
                ax.set_ylabel(r"$r-i$", fontsize=16)

                ax.figure.savefig("color-color2/all.png")

