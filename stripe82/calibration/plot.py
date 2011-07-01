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
import database
from opt import survey

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
    print data.stars
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

def plot_grid(basepath='.'):
    """
    Plot the calibration model parameters as a function of RA/Dec
    
    Parameters
    ----------
    basepath : string
        Path to save files in
    
    History
    -------
    2011-06-16 - Created by Dan Foreman-Mackey
    
    """
    plotdir = os.path.join(basepath,'grid')
    if os.path.exists(plotdir):
        shutil.rmtree(plotdir)
    os.makedirs(plotdir)
    
    runs = []
    ras = []
    models = []
    zeros = []
    for doc in database.obslist.find():
        run = survey.get_observation(doc['obsid'])['run']
        runs.append(run)
        model = database.photomodel.find_one({'_id': doc['modelid']})
        models.append(model)
        ras.append(model['ra'])
        zeros.append(doc['zero'])
    uniqueruns = set(runs)
    runs = np.array(runs)
    zeros = np.array(zeros)
    ras = np.array(ras)
    
    pl.figure(figsize=(7,14))
    for run in uniqueruns:
        pl.clf()
        pl.subplot(411)
        inds = np.arange(len(runs))[runs == run]
        pl.plot(ras[inds],10**(-zeros[inds]/2.5-9),'.k')
        pl.ylim([1e2,4e2])
        pl.ylabel(r'$\mathrm{ADU\,nMgy}^{-1}$',fontsize=16.)
        pl.gca().set_xticklabels([])
        pl.title(r'$\mathrm{run} = %d$'%run,fontsize=16.)

        ax1 = pl.gcf().add_subplot(412)
        ax2 = pl.gcf().add_subplot(413)
        ax3 = pl.gcf().add_subplot(414)
        for i in inds:
            model = models[i]['model']
            ax3.plot(ras[i],model.jitterabs2,'.k')
            ax2.plot(ras[i],model.jitterrel2,'.k')
            ax1.plot(ras[i],model.data.nstars,'.k')
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])

        ax3.set_ylabel(r'$\delta^2$',fontsize=16.)
        ax2.set_ylabel(r'$\eta^2$',fontsize=16.)
        ax1.set_ylabel(r'$N_\mathrm{stars}$',fontsize=16.)

        ax3.set_xlabel(r'$\mathrm{R.A.}\,[\mathrm{deg}]$',fontsize=16.)
        pl.savefig(os.path.join(plotdir,'%04d.png'%run))
        
if __name__ == '__main__':
    #plot_lightcurves(-23.431965,-0.227934,period=0.48399)
    plot_grid()


