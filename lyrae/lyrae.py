#!/usr/bin/env python
# encoding: utf-8
"""
Lyrae module

History
-------
2011-06-22 - Created by Dan Foreman-Mackey

"""

__all__ = ['fit']

import numpy as np
from numpy.linalg import lstsq

def fit(omega,time,data):
    """
    Fit a general RR Lyrae model to a time series

    Solve A.x = F
    
    Parameters
    ----------
    time : numpy.ndarray (N,)
        The list of time steps
    
    data : numpy.ndarray (N,2)
        (flux,counts) for each time step

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
    A = np.ones((len(time),7))
    A[:,1] = np.sin(omega*time)
    A[:,2] = np.cos(omega*time)
    A[:,3] = np.sin(2*omega*time)
    A[:,4] = np.cos(2*omega*time)
    A[:,5] = np.sin(3*omega*time)
    A[:,6] = np.cos(3*omega*time)

    res = lstsq(A,data[:,0])[0]
    model = lambda t: res[0]+res[1]*np.sin(omega*t)+res[2]*np.cos(omega*t)+ \
        res[3]*np.sin(2*omega*t)+res[4]*np.cos(2*omega*t)+ \
        res[5]*np.sin(3*omega*t)+res[6]*np.cos(3*omega*t)
    
    chi2 = np.sum((data[:,0]-model(time))**2/data[:,1]**2)

    return model, chi2

if __name__ == '__main__':
    import h5py
    import time as timer
    f = h5py.File('testlc.hdf5')
    data = f['data'][...]
    time = f['time'][...]
    f.close()

    inds = data[:,0] > 0
    data = data[inds,:]
    time = time[inds]
    data[:,0] = -2.5*np.log10(data[:,0])
    data[:,1] = 2.5/data[:,0]/np.log(10)*data[:,1]

    domega = 1.0/(time.max()-time.min())
    omegas = np.arange(11.0,15.0,domega)
    print "omega_min,omega_max,d_omega = ",omegas.min(),omegas.max(),domega
    print "nomega = ",len(omegas)
    chi2 = []
    strt = timer.time()
    for omega in omegas:
        model,diff = fit(omega,time,data)
        chi2.append(diff)
    print "time:",timer.time()-strt
    import pylab as pl
    T = 2*np.pi/omegas
    pl.plot(T,chi2,'k')
    pl.xlim([T.min(),T.max()])
    pl.xlabel(r'$T\,[\mathrm{days}]$',fontsize=16.)
    pl.ylabel(r'$\chi^2$',fontsize=16.)

    truth = 0.48399
    T0 = float(T[np.array(chi2) == min(chi2)])
    pl.title('"truth" = %f days --- best-fit = %f days'%(truth,T0))
    pl.gca().axvline(truth,color='r',ls='--')
    
    pl.savefig('chi2.pdf')
