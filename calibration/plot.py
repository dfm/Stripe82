#!/usr/bin/env python
# encoding: utf-8
"""
Plot results of calibration

History
-------
2011-06-16 - Created by Dan Foreman-Mackey

"""

__all__ = ['']

import os
import os.path
import shutil

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

import numpy as np

from calibrate import *
from model import *

def plot_lightcurves(ra,dec,radius=3,basepath='.',period=None):
    """
    Plot the lightcurves for the stars around RA/Dec
    
    Parameters
    ----------
    ra : float
        RA in degrees

    dec : float
        Dec. in degrees
    
    Optional
    --------
    radius : float (default : 3.0)
        Search radius (passed to survey.find_stars) in arcmin

    basepath : string (default : '.')
        Directory to save figures

    period : float (default : None)
        If not None, all of the plots will be folded with this period
        (in MJD)

    History
    -------
    2011-06-16 - Created by Dan Foreman-Mackey
    
    """
    model = calibrate_grid([[ra,dec]])[0]
    data = model.data
    bp = basepath

    lcdir = os.path.join(bp,'lcs')
    if os.path.exists(lcdir):
        shutil.rmtree(lcdir)
    os.makedirs(lcdir)

    # array of probability that each observation is bad
    lnprobbad = lnprob_badobs(model,data)
    lnprobvar = lnprob_variable(model,data)

    for si in range(data.nstars):
        pl.clf()

        # select where inverse variance is non-zero
        inds = data.ivar[:,si] > 1e-8
        if np.any(inds):
            t = data.mjd()[inds]
            if period is not None:
                t = (t%period)/period
            df = data.flux[inds,si]
            dferr = 1/np.sqrt(data.ivar[inds,si])

            # learned errorbars
            delta2 = model.jitterabs2 + model.jitterrel2*\
                    (model.zero[inds]*model.flux[si])**2
            lrnerr = np.sqrt(dferr**2+delta2)

            # use calibration
            df /= model.zero[inds]
            dferr /= model.zero[inds]
            lrnerr /= model.zero[inds]
            dmag = -2.5*np.log10(df)
            dmagerr = 2.5/df/np.log(10)*dferr
            lrnmagerr = 2.5/df/np.log(10)*lrnerr

            # calculate colors and sizes MAGIC MAGIC MAGIC
            # don't ask...
            size = -lnprobbad[inds,:]
            size -= min(size)
            size /= max(size)*2
            size += 0.75
            clrs = np.log(2-size)
            clrs /= max(clrs)
            size = 10*size**2

            pl.scatter(t,dmag,c=clrs,zorder=10)
            pl.errorbar(t,dmag,yerr=lrnmagerr,fmt='.',color="#666666")
            pl.errorbar(t,dmag,yerr=dmagerr,fmt='.k')
            if period is not None:
                pl.errorbar(t+1,dmag,yerr=lrnmagerr,fmt='.',color="#666666")
                pl.errorbar(t+1,dmag,yerr=dmagerr,fmt='.k')
                pl.scatter(t+1,dmag,c=clrs,zorder=10)
                pl.xlim([0,2])

            # mean flux
            pl.axhline(-2.5*np.log10(model.flux[si]),color='k',ls='-')

            # flux from DAS
            pl.axhline(data.magprior[si,0],color='r',ls=':')
            xlim = pl.gca().get_xlim()
            p,m = data.magprior[si,0]+data.magprior[si,1],\
                    data.magprior[si,0]-data.magprior[si,1]
            pl.fill(np.append(xlim,xlim[::-1]),[p,p,m,m],color='r',alpha=0.3,zorder=-3)
            
            pl.ylim(pl.gca().get_ylim()[::-1])
            pl.ylabel(r'$m_g$',fontsize=16.)
            pl.xlabel(r'$t/T$',fontsize=16.)
            pl.title(r'$\ln\,r^\mathrm{var} = %.0f$'%(lnprobvar[si]),fontsize=16.)

            pl.savefig(os.path.join(lcdir,'%04d.pdf'%si))

if __name__ == '__main__':
    plot_lightcurves(-23.431965,-0.227934,period=0.48399)


