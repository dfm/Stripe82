#!/usr/bin/env python
# encoding: utf-8
"""
Command line utility for extracting a light-curve from Stripe 82 at a given RA/Dec

History
-------
2011-08-03 - Created by Dan Foreman-Mackey

"""

__all__ = ['extract_lightcurve']

import os
import cPickle as pickle

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

import numpy as np

from calibration import calibrate,PhotoModel,odds_variable,odds_bad
from lyrae import sesar,fit,find_period

def extract_lightcurve(ra,dec,radius):
    return calibrate(ra,dec,radius=radius)

def _angular_distance(ra1,dec1,ra2,dec2):
    ra1,dec1 = np.radians(ra1),np.radians(dec1)
    ra2,dec2 = np.radians(ra2),np.radians(dec2)
    crho = np.cos(dec1)*np.cos(dec2)\
            *(np.cos(ra1)*np.cos(ra2)+np.sin(ra1)*np.sin(ra2))\
            +np.sin(dec1)*np.sin(dec2)
    return np.degrees(np.arccos(crho))

if __name__ == '__main__':
    tmpfn = 'tmp.pkl'
    ra,dec = 29.47942,0.383557 #10.0018734334081,0.791580301596976
    if os.path.exists(tmpfn):
        model = PhotoModel(*pickle.load(open(tmpfn,'rb')))
    else:
        model = calibrate(ra,dec,radius=10)
        pickle.dump((model.data,model.vector()),open(tmpfn,'wb'),-1)
    mjd = model.data.mjd()
    flux = model.data.flux/model.zero[:,np.newaxis]*1e9
    varodds = odds_variable(model,model.data)
    badodds = odds_bad(model,model.data)

    # plot 1
    pl.figure(figsize=(8.,15.))

    ax1 = pl.subplot(311)
    ax1.plot(mjd,model.zero/1e9,'.k')
    ax1.set_xticklabels([])
    ax1.set_ylabel(r'$(\mathrm{nMgy/ADU})_i$',fontsize=16.)

    Nstars = np.sum(model.data.ivar>0,axis=1)
    ax2 = pl.subplot(312)
    ax2.plot(mjd,np.sum(badodds,axis=1)/Nstars,'.k')
    ax2.set_xticklabels([])
    ax2.set_ylabel(r'$(1/N_\mathrm{stars})\sum_\alpha \ln \, r^\mathrm{bad}_{i\alpha}$',fontsize=16.)

    ax3 = pl.subplot(313)
    ax3.plot(mjd,Nstars,'.k')
    ax3.set_ylabel(r'$N_\mathrm{stars}$',fontsize=16.)
    ax3.set_xlabel(r'$\mathrm{MJD}$',fontsize=16.)

    pl.savefig('lcs/plot1.png')

    # plot 2
    pl.figure(figsize=(8.,15.))
    inds = np.argsort(model.flux)
    sids = np.arange(len(model.flux))
    dist = lambda i: _angular_distance(ra,dec,
            model.data.stars[inds[i]]['ra'],model.data.stars[inds[i]]['dec'])
    target_id = inds[sorted(range(len(inds)),key = dist)[0]]
    target = model.data.stars[target_id]

    ax1 = pl.subplot(311)
    ax1.plot(sids,model.flux[inds]*1e9,'.k')
    ax1.axvline(target_id,color='r')
    ax1.set_xticklabels([])
    ax1.set_ylabel(r'$f^*_\alpha\,[\mathrm{nMgy}]$',fontsize=16.)

    ax2 = pl.subplot(312)
    ax2.plot(sids,varodds[inds],'.k')
    ax2.axvline(target_id,color='r')
    ax2.set_xticklabels([])
    ax2.set_ylabel(r'$\ln \, r^\mathrm{var}_\alpha$',fontsize=16.)

    ax3 = pl.subplot(313)
    ax3.plot(sids,np.mean(badodds,axis=0)[inds],'.k')
    ax3.axvline(target_id,color='r')
    ax3.set_ylabel(r'$(1/N_\mathrm{obs})\sum_i \ln \, r^\mathrm{bad}_{i\alpha}$',fontsize=16.)
    ax3.set_xlabel(r'$\mathrm{Star\,ID}$',fontsize=16.)

    pl.savefig('lcs/plot2.png')
    


    #if np.any(sesar.coords['ra'] == ra):
    #    ind = sesar.coords[sesar.coords['ra'] == ra]['Num']
    #    period = sesar.coords[sesar.coords['ra'] == ra]['Per']
    #    sesardata = sesar.table1['%d'%(ind)][...]
    #    inds = sesardata['g'] > 0
    #    s_time = sesardata['gmjd'][inds]
    #    s_data = 1e9*10**(-sesardata['g']/2.5)[inds]
    #else:
    #    s_data = None
    #    period = None
    #
    #if os.path.exists(tmpfn):
    #    model = PhotoModel(*pickle.load(open(tmpfn,'rb')))
    #else:
    #    model = calibrate(ra,dec,radius=10)
    #    pickle.dump((model.data,model.vector()),open(tmpfn,'wb'),-1)
    #mjd = model.data.mjd()
    #if period is None:
    #    period = max(mjd)+1
    #flux = model.data.flux/model.zero[:,np.newaxis]*1e9
    #varodds = odds_variable(model,model.data)

    #if False: # find period
    #    i = 24
    #    inv = model.data.ivar
    #    inds = inv[:,i] > 0
    #    err0_i = 1e9/np.sqrt(inv[inds,i]*model.zero[inds]**2)
    #    flux_i = flux[inds,i]
    #    mjd_i = mjd[inds]
    #    
    #    data = np.zeros([len(flux_i),2])
    #    data[:,0] = flux_i
    #    data[:,1] = err0_i
    #    my_period = find_period(mjd_i,data)
    #    print my_period,period
    #    period = my_period
    #

    #for i in range(np.shape(flux)[1]):
    #    inv = model.data.ivar
    #    inds = inv[:,i] > 0
    #    err0_i = 1e9/np.sqrt(inv[inds,i]*model.zero[inds]**2)
    #    # learning error bars
    #    #err_i  = 1/np.sqrt(inv[inds,i]*model.zero[inds]**2+model.jitterabs2 \
    #    #        +model.jitterrel2*model.flux[i]**2)*(1e9)
    #    flux_i = flux[inds,i]
    #    mjd_i = mjd[inds]

    #    pl.clf()
    #    if s_data is not None:
    #        pl.plot(s_time%period,s_data,'og',alpha=0.3)
    #        pl.plot(s_time%period+period,s_data,'og',alpha=0.3)
    #    
    #    #pl.errorbar(mjd_i%period,flux_i,yerr=err_i,fmt='.k',alpha=0.5)
    #    pl.errorbar(mjd_i%period,flux_i,yerr=err0_i,fmt='.k')
    #    if s_data is not None:
    #        pl.errorbar(mjd_i%period+period,flux_i,yerr=err0_i,fmt='.k')

    #    pl.gca().axhline(model.flux[i]*1e9,color='r',ls='--')

    #    if s_data is not None:
    #        data = np.zeros([len(flux_i),2])
    #        data[:,0] = flux_i
    #        lyrae_fit,chi2 = fit(2*np.pi/period,mjd_i,data)
    #        x = np.linspace(0,period*2,500)
    #        pl.plot(x,lyrae_fit(x),'-k')

    #    star = model.data.stars[i]
    #    pl.title('(%f,%f) and %f'%(star['ra'],star['dec'],varodds[i]))
    #    print i,star['ra'],star['dec'],varodds[i]
    #    pl.savefig('lcs/%03d.png'%i)

