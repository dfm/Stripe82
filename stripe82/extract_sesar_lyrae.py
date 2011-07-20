#!/usr/bin/env python
# encoding: utf-8
"""


History
-------
2011-06-23 - Created by Dan Foreman-Mackey

"""

import os
import os.path

import time as timer

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl

from calibration import extract_calibrated_lightcurve
import lyrae.sesar as sesar
from lyrae import fit



for obj in sesar.coords:
    sid = obj['Num']
    print 'Lyrae number: ',sid

    bp = os.path.join('lcs','%d'%sid)
    if not os.path.exists(bp):
        os.makedirs(bp)

    sesardata = sesar.table1['%d'%sid][...]
    inds = sesardata['g'] > 0
    s_time = sesardata['gmjd'][inds]
    s_data = 1e9*10**(-sesardata['g']/2.5)[inds]

    time,data = extract_calibrated_lightcurve(obj['ra'],obj['dec'],radius=5)

    inds = data[:,0] > 0
    data = data[inds,:]
    time = time[inds]
    data[:,0] *= 1e9
    data[:,1] *= 1e9
    inds = data[:,0] > data[:,1]
    time = time[inds]
    data = data[inds,:]

    # grid in frequency
    domega = 1.0/(time.max()-time.min())
    omegas = 2*np.pi*np.arange(1.0/1.3,1.0/0.2,domega)
    print "omega_min,omega_max,d_omega = ",omegas.min(),omegas.max(),domega
    print "nomega = ",len(omegas)
    chi2 = []
    strt = timer.time()
    for omega in omegas:
        model,diff = fit(omega,time,data)
        chi2.append(diff)
    print "time:",timer.time()-strt
    T = 2*np.pi/omegas
    inds = np.argsort(chi2)
    sortedT = T[inds]
    sortedchi2 = np.array(chi2)[inds]

    # plotting
    pl.figure(figsize=(8.,8.))
    pl.subplot(211)
    pl.plot(T,chi2,'k')
    pl.xlim([T.min(),T.max()])
    pl.xlabel(r'$T\,[\mathrm{days}]$',fontsize=16.)
    pl.ylabel(r'$\chi^2$',fontsize=16.)
    truth = sesar.coords['Per'][sesar.coords['Num'] == sid][0]
    print truth
    T0 = float(T[np.array(chi2) == min(chi2)])
    pl.title('"truth" = %f days --- best-fit = %f days'%(truth,T0))
    pl.gca().axvline(truth,color='r',ls='--')

    pl.subplot(212)
    pl.plot(T,chi2,'k')
    pl.gca().axvline(truth,color='r',ls='--')
    pl.xlim([truth-0.005,truth+0.005])
    pl.xlabel(r'$T\,[\mathrm{days}]$',fontsize=16.)
    pl.ylabel(r'$\chi^2$',fontsize=16.)

    pl.savefig(os.path.join(bp,'chi2.png'))

    pl.figure()
    omegas = 2*np.pi/sortedT[:10]
    omegas = np.append(2*np.pi/truth,omegas)
    for i,omega in enumerate(omegas):
        pl.clf()
        T = 2*np.pi/omega

        # sesar
        pl.plot(s_time%T,s_data,'.r')
        pl.plot(s_time%T+T,s_data,'.r')

        pl.errorbar(time%T,data[:,0],yerr=data[:,1],fmt='.k')
        pl.errorbar(time%T+T,data[:,0],yerr=data[:,1],fmt='.k')

        for order in np.arange(5)+1:
            model,diff = fit(omega,time,data,order=order)
            t = np.linspace(0,2*T,2000)
            f = model(t)
            pl.plot(t,f,'--',label=r'$%d,\,\chi^2 = %.2e$'%(order,diff))

        pl.legend()
        pl.xlim([0,2*T])
        pl.xlabel(r'$t\%T\,[\mathrm{days}]$',fontsize=16.)
        pl.ylabel(r'$\mathrm{nMgy}$',fontsize=16.)
        pl.title(r'$T=%f\,\mathrm{days}$'%(T),fontsize=16.)
        pl.savefig(os.path.join(bp,'%04d.png'%i))


