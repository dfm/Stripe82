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

from calibration import calibrate,PhotoModel,odds_variable
from lyrae import sesar,fit,find_period

def extract_lightcurve(ra,dec,radius):
    return calibrate(ra,dec,radius=radius)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fit M31 dynamics')
    parser.add_argument('basepath',help='Directory for results')
    parser.add_argument('f','tmpfile',help='Temp pickle file',
                        default=None)
    parser.add_argument('--all',
                        help='all figures',
                        action='store_true')
    parser.add_argument('-r','-a','--ra',default=29.47942)
    parser.add_argument('-d','--dec',default=0.383557)

    args = parser.parse_args()
    
    ra,dec = args.ra,args.dec#29.47942,0.383557 #10.0018734334081,0.791580301596976

    if np.any(sesar.coords['ra'] == ra):
        ind = sesar.coords[sesar.coords['ra'] == ra]['Num']
        period = sesar.coords[sesar.coords['ra'] == ra]['Per']
        sesardata = sesar.table1['%d'%(ind)][...]
        inds = sesardata['g'] > 0
        s_time = sesardata['gmjd'][inds]
        s_data = 1e9*10**(-sesardata['g']/2.5)[inds]
    else:
        s_data = None
        period = None
    
    if args.tmpfile is not None and os.path.exists(args.tmpfile):
        model = PhotoModel(*pickle.load(open(args.tmpfile,'rb')))
    else:
        model = calibrate(ra,dec,radius=10)
        if args.tmpfile is not None:
            pickle.dump((model.data,model.vector()),open(args.tmpfile,'wb'),-1)
    mjd = model.data.mjd()
    if period is None:
        period = max(mjd)+1
    flux = model.data.flux/model.zero[:,np.newaxis]*1e9
    varodds = odds_variable(model,model.data)

    if True: # find period
        i = 24
        inv = model.data.ivar
        inds = inv[:,i] > 0
        err0_i = 1e9/np.sqrt(inv[inds,i]*model.zero[inds]**2)
        flux_i = flux[inds,i]
        mjd_i = mjd[inds]
        
        data = np.zeros([len(flux_i),2])
        data[:,0] = flux_i
        data[:,1] = err0_i
        my_period = find_period(mjd_i,data)
        print my_period,period
        period = my_period
    

    for i in range(np.shape(flux)[1]):
        inv = model.data.ivar
        inds = inv[:,i] > 0
        err0_i = 1e9/np.sqrt(inv[inds,i]*model.zero[inds]**2)
        # learning error bars
        #err_i  = 1/np.sqrt(inv[inds,i]*model.zero[inds]**2+model.jitterabs2 \
        #        +model.jitterrel2*model.flux[i]**2)*(1e9)
        flux_i = flux[inds,i]
        mjd_i = mjd[inds]

        pl.clf()
        if s_data is not None:
            pl.plot(s_time%period,s_data,'og',alpha=0.3)
            pl.plot(s_time%period+period,s_data,'og',alpha=0.3)
        
        #pl.errorbar(mjd_i%period,flux_i,yerr=err_i,fmt='.k',alpha=0.5)
        pl.errorbar(mjd_i%period,flux_i,yerr=err0_i,fmt='.k')
        if s_data is not None:
            pl.errorbar(mjd_i%period+period,flux_i,yerr=err0_i,fmt='.k')

        pl.gca().axhline(model.flux[i]*1e9,color='r',ls='--')

        if s_data is not None:
            data = np.zeros([len(flux_i),2])
            data[:,0] = flux_i
            lyrae_fit,chi2 = fit(2*np.pi/period,mjd_i,data)
            x = np.linspace(0,period*2,500)
            pl.plot(x,lyrae_fit(x),'-k')

        star = model.data.stars[i]
        pl.title('(%f,%f) and %f'%(star['ra'],star['dec'],varodds[i]))
        print i,star['ra'],star['dec'],varodds[i]
        pl.savefig('lcs/%03d.png'%i)

