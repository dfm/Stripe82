#!/usr/bin/env python
# encoding: utf-8
"""
A couple of wrappers around the results of a calibration optimization

History
-------
2011-08-16 - Created by Dan Foreman-Mackey

"""

__all__ = ['CalibrationPatch','Lightcurve']

import numpy as np
import matplotlib
import matplotlib.pyplot as pl

from extract_lightcurve import extract_lightcurves
import calibration
from calibration.opt import *
import lyrae

class CalibrationPatch:
    """
    A patch where the calibration model has been optimized

    Parameters
    ----------
    model : PhotoModel
        The optimized model

    data : PhotoData
        The data in the patch

    ra : float
        The patch center

    dec : float
        The patch center

    radius : float
        The radius in which stars were selected [arcmin]

    History
    -------
    2011-08-16 - Created by Dan Foreman-Mackey

    """
    def __init__(self,model,data,ra,dec,radius):
        self._model  = model
        self._data   = data
        self._ra     = ra
        self._dec    = dec
        self._radius = radius
        mjd,flux,err,m = extract_lightcurves(model=model)
        self._lnoddsbad = calibration.odds_bad(model,data)
        self._lnoddsvar = calibration.odds_variable(model,data)
        self._lightcurves = [Lightcurve(mjd,flux[:,i],err[:,i],self,
                                        self._lnoddsbad[:,i],self._lnoddsvar[i],
                                        model.flux[i])
                                for i in range(np.shape(flux)[-1])]

class Lightcurve:
    """
    An individual lightcurve

    Parameters
    ----------
    starnumber : int
        Index in model.data.stars

    mjd : list
        The time stamps

    flux : list
        The list of fluxes

    err : list
        The uncertainties

    calibpatch : CalibrationPatch
        The associated fit

    lnoddsbad : list
        The list of odds that an observations is bad

    lnoddsvar : float
        Odds that this star is variable

    meanflux : float
        The fitted man flux

    Optional
    --------
    period : float
        If the period is known you can supply it

    History
    -------
    2011-08-16 - Created by Dan Foreman-Mackey

    """
    def __init__(self,starnumber,mjd,flux,err,calibpatch,
            lnoddsbad,lnoddsvar,meanflux,period=None):
        inds = ~err.mask[:,i]
        self._starnumber = starnumber
        self._mjd  = mjd[inds]
        self._err  = err[inds]
        self._flux = flux[inds]
        self._meanflux = meanflux
        self._lnoddsbad = lnoddsbad[inds]
        self._lnoddsvar = lnoddsvar
        self._calibpatch = calibpatch
        self._period = period
        self._spline = None

    def __unicode__(self):
        star = self._calibpatch.data.stars[self._starnumber]
        return u"Lightcurve for star at (%.2f,%.2f) in %r"%\
                (star['ra'],star['dec'],self._calibpatch)

    def __str__(self):
        return unicode(self)

    def get_period(self):
        """
        Return the period (fit for it if it doesn't exist)

        Returns
        -------
        period : float
            The fit period

        History
        -------
        2011-08-16 - Created by Dan Foreman-Mackey

        """
        if self._period is not None:
            return self._period
        data = np.zeros([len(self._mjd),2])
        data[:,0] = self._flux
        data[:,1] = self._err
        self._period = lyrae.find_period(self._mjd,data)
        self._spline,chi2 = lyrae.fit(2*np.pi/self._period,self._mjd,data)
        return self._period

    def plot(self,ax=None,period=None,nperiods=None):
        """
        Plot the lightcurve

        Parameters
        ----------
        ax : matplotlib.Axes
            Axis to plot in

        period : float (defualt : None)
            The period to fold on

        nperiod : int (default : None)
            The number of periods to plot. If None and period is None, don't fold.

        History
        -------
        2011-08-16 - Created by Dan Foreman-Mackey

        """
        if ax is None:
            ax = pl.gca()
        if nperiods is not None and period is None:
            period = self.get_period()
        if nperiods is None:
            nperiods = 1
        if period is None:
            period = max(mjd) + 1

        clrs = logodds2prob(self._lnoddsbad)
        for n in range(nperiods):
            ax.errorbar(self._mjd%period+n*period,self._flux,yerr=self._err,
                    fmt='.k',zorder=-1,capsize=0)
            ax.scatter(mjd[inds]%period+n*period,flux[inds,i],
                    c=clrs,edgecolor='k',zorder=100,cmap='gray',
                    vmin=0.0,vmax=1.0)

        if self._spline is not None:
            ts = np.linspace(0,nperiods*period,nperiods*500)
            ax.plot(ts,self._spline(ts),'k')

        ax.axhline(self.meanflux,color="k",ls="--")

        ax.set_xlim([0.0,nperiods*period])
        ax.set_ylim([0,2*self._meanflux])

        ax.set_ylabel(r'$\mathrm{nMgy}$',fontsize=16)
        ax.set_xlabel(r'$\mathrm{days}$',fontsize=16)

        ax.set_title(r'$\ln\,p_\mathrm{var} = %.3f$'%\
                (np.log(logodds2prob(self._lnoddsvar))),fontsize=16)

if __name__ == '__main__':
    pass

