#!/usr/bin/env python
# encoding: utf-8
"""
Lyrae module

History
-------
2011-06-22 - Created by Dan Foreman-Mackey

"""

__all__ = ['fit','find_period']

import numpy as np
from numpy.linalg import lstsq

def fit(omega,time,data,order=3,full_output=False):
    """
    Fit a general RR Lyrae model to a time series

    Solve A.x = F

    Parameters
    ----------
    time : numpy.ndarray (N,)
        The list of time steps

    data : numpy.ndarray (N,2)
        (flux,counts) for each time step

    Optional
    --------
    order : int (default : 3)
        What order of model should we use

    Returns
    -------
    model : function
        A function that returns the model flux at a given time

    History
    -------
    2011-06-22 - Created by Dan Foreman-Mackey

    """
    # MAGIC: construct A matrix
    # MAGIC: 7 parameters in model
    A = np.ones((len(time),1+2*order))
    for i in range(order):
        A[:,2*i+1] = np.sin((i+1)*omega*time)
        A[:,2*i+2] = np.cos((i+1)*omega*time)

    res = lstsq(A,data[:,0])[0]
    model = lambda t: res[0]+\
            np.sum([res[2*i+1]*np.sin((i+1)*omega*t)+\
                    res[2*i+2]*np.cos((i+1)*omega*t) for i in range(order)],axis=0)

    if full_output:
        return res

    chi2 = np.sum((data[:,0]-model(time))**2/data[:,1]**2)

    return model, chi2

def find_period(time,data,order=3,nmin=10,N=50,full_output=False):
    """
    Find the period of some data

    time : numpy.ndarray (N,)
        The list of time steps

    data : numpy.ndarray (N,2)
        (flux,counts) for each time step

    Optional
    --------
    order : int (default : 3)
        What order of model should we use

    nmin : int (default : 10)
        How many chi2 minima should we investigate

    N : int (default : 100)
        How many samples in iteration

    Returns
    -------
    period : float
        In the same units as time

    History
    -------
    2011-08-05 - Created by Dan Foreman-Mackey

    """
    # initial grid in frequency
    domega = 0.3/(time.max()-time.min()) # times 2pi to be precise
    omegas = 2*np.pi*np.arange(1.0/1.3,1.0/0.2,domega)
    print "omega_min,omega_max,d_omega = ",omegas.min(),omegas.max(),2*np.pi*domega
    print "nomega = ",len(omegas)
    chi2 = []
    for omega in omegas:
        model,diff = fit(omega,time,data)
        chi2.append(diff)
    inds = np.argsort(chi2)

    minchi2 = chi2[inds[0]]
    omega_f = omegas[inds[0]]
    for i in inds[:nmin]:
        print "T_%d = "%i,2*np.pi/omegas[i]
        ws = np.linspace(omegas[i-1],omegas[i+1],N)
        c2 = []
        for w in ws:
            model,diff = fit(w,time,data)
            c2.append(diff)
        mc2 = min(c2)
        if mc2 < minchi2:
            minchi2 = mc2
            omega_f = ws[c2.index(mc2)]

    print "Final T = ",2.0*np.pi/omega_f
    if full_output:
        return 2.0*np.pi/omega_f, fit(omega_f,time,data,full_output=True)

    return 2.0*np.pi/omega_f

if __name__ == '__main__':
    import h5py
    import time as timer
    import sesar

    sid = 1052471
    sesardata = sesar.table1['%d'%sid][...]
    inds = sesardata['g'] > 0
    s_time = sesardata['gmjd'][inds]
    s_data = 1e9*10**(-sesardata['g']/2.5)[inds]

    f = h5py.File('testlc.hdf5')
    data = f['data'][...]
    time = f['time'][...]
    f.close()

    inds = data[:,0] > 0
    data = data[inds,:]
    time = time[inds]
    data[:,0] *= 1e9
    data[:,1] *= 1e9
    #data[:,0] = -2.5*np.log10(data[:,0])
    #data[:,1] = 2.5/data[:,0]/np.log(10)*data[:,1]

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
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as pl

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

    pl.savefig('chi2.png')

    pl.figure()
    omegas = 2*np.pi/sortedT[:10]
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
        pl.savefig('grid/%04d.png'%i)


