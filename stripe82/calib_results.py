#!/usr/bin/env python
# encoding: utf-8
"""
A couple of wrappers around the results of a calibration optimization

History
-------
2011-08-16 - Created by Dan Foreman-Mackey

"""

__all__ = ['CalibrationPatch','Lightcurve']

import cPickle as pickle

import numpy as np
import matplotlib
from matplotlib import rc
rc('font',**{'family':'serif','serif':'Computer Modern Roman'})
rc('text', usetex=True)
import matplotlib.pyplot as pl

from calibration.extract import extract_lightcurves
import calibration
from calibration.opt import *
import lyrae
from opt import *

def _angular_distance(ra1,dec1,ra2,dec2):
    ra1,dec1 = np.radians(ra1),np.radians(dec1)
    ra2,dec2 = np.radians(ra2),np.radians(dec2)
    crho = np.cos(dec1)*np.cos(dec2)\
            *(np.cos(ra1)*np.cos(ra2)+np.sin(ra1)*np.sin(ra2))\
            +np.sin(dec1)*np.sin(dec2)
    return np.degrees(np.arccos(crho))

#def unpickle_calibpatch(fn):
#    (vector,self._data,self._ra,self._dec,
#        self._radius,self.lnoddsbad,self.lnoddsvar,self._lightcurves,
#        self._target) = load(open(fn,"rb"))

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
        self.lnoddsbad = calibration.odds_bad(model,data)
        self.lnoddsvar = calibration.odds_variable(model,data)
        self._lightcurves = [Lightcurve(i,mjd,flux[:,i],err[:,i],self,
                                        self.lnoddsbad[:,i],
                                        self.lnoddsvar[i],
                                        model.flux[i])
                                for i in range(np.shape(flux)[-1])]
        dist = lambda i: _angular_distance(ra,dec,
                model.data.stars[i]['ra'],model.data.stars[i]['dec'])
        self._target = sorted(np.arange(len(model.flux)),key = dist)[0]

    def dump(self,fn):
        # synced with unpickle_calibpatch
        data = (self._model.vector(),self._data,self._ra,self._dec,
                self._radius,self.lnoddsbad,self.lnoddsvar,self._lightcurves,
                self._target)
        pickle.dump(data,open(fn,"wb"),-1)

    def get_target(self):
        """
        Get the target lightcurve

        Returns
        -------
        target_id : int
            The index of the target star

        lightcurve : Lightcurve
            The lightcurve object of the target star

        History
        -------
        2011-08-16 - Created by Dan Foreman-Mackey

        """
        return self._target,self._lightcurves[self._target]

    def get_radec(self,starindex):
        star = self._data.stars[starindex]
        return star['ra'],star['dec']

    def get_nstars(self):
        """
        Return the number of stars

        Returns
        -------
        nstars : int
            Number of stars used

        History
        -------
        2011-08-16 - Created by Dan Foreman-Mackey

        """
        return len(self._lightcurves)

    def get_model(self):
        """
        Get the model

        Returns
        -------
        model : PhotoModel
            The optimized model

        History
        -------
        2011-08-16 - Created by Dan Foreman-Mackey

        """
        return self._model

    def plot(self,period=None,nperiods=None):
        """
        Plot the lightcurves and return the list of figures

        Parameters
        ----------
        period : float
            Define a period to fold on

        nperiods : int
            Number of periods to plot

        Returns
        -------
        figs : list
            List of matplotlib figures

        History
        -------
        2011-08-16 - Created by Dan Foreman-Mackey

        """
        figs = []
        for lc in self._lightcurves:
            fig = pl.figure()
            figs.append(fig)
            lc.plot(ax=fig.add_subplot(111),period=period,nperiods=nperiods)
        return figs


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
        inds = ~err.mask
        self._starnumber = starnumber
        self._mjd  = mjd[inds]
        self._err  = err[inds]#np.sqrt(err[inds]**2+model.jitterabs2+\
        #model.jitterrel2*meanflux**2)
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

    @property
    def mjd(self):
        return self._mjd

    def get_period(self,subscript=None):
        """
        Return the period (fit for it if it doesn't exist)

        Parameters
        ----------
        subscript : slice
            Slice of the data to use for the fit. Defaults to all data.

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
        if subscript is None:
            subscript = np.arange(len(self._mjd))
        self._period,self._spline,chi2 = self.get_period_in_data(self._mjd[subscript],
                self._flux[subscript],self._err[subscript],full_output=True)
        return self._period

    def get_period_in_data(self,mjd,flux,err,full_output=False):
        data = np.zeros([len(mjd),2])
        data[:,0] = flux
        data[:,1] = err
        period = lyrae.find_period(mjd,data)
        if full_output:
            spline,chi2 = lyrae.fit(2*np.pi/period,mjd,data)
            return period,spline,chi2
        return period

    def get_lnoddsvar(self):
        """
        Get the odds that a lightcurve is variable

        Returns
        -------
        lnoddsvar : float
            The odds that the LC is variable

        History
        -------
        2011-08-16 - Created by Dan Foreman-Mackey

        """
        return self._lnoddsvar

    def get_lnpvar(self):
        """
        Get the probability that a lightcurve is variable

        Returns
        -------
        lnpvar : float
            The probability that a star is variable

        History
        -------
        2011-08-16 - Created by Dan Foreman-Mackey

        """
        return np.log(logodds2prob(self._lnoddsvar))

    def plot(self,*args,**kwargs):
        """
        Plot the lightcurve

        Parameters
        ----------
        ax : matplotlib.Axes
            Axis to plot in

        period : float (defualt : None)
            The period to fold on

        nperiods : int (default : None)
            The number of periods to plot. If None and period is None, don't fold.

        calcspline : bool (default : False)
            Calculate the spline fit if we don't already have one?

        hyperparams : bool (default : False)
            List the hyper-parameters.

        History
        -------
        2011-08-16 - Created by Dan Foreman-Mackey

        """
        self.plot_slice(np.arange(len(self._mjd)),*args,**kwargs)

    def plot_slice(self,subscript,ax=None,period=None,nperiods=None,calcspline=False,
            hyperparams=False,show_title=True):
        if ax is None:
            ax = pl.gca()
        if nperiods is not None and period is None:
            period = self.get_period(subscript=subscript)
            print period
        if nperiods is None:
            nperiods = 1
        #print "variance=",np.var(self._flux),self._meanflux
        clrs = logodds2prob(self._lnoddsbad[subscript])
        if period is not None:
            for n in range(nperiods):
                ax.errorbar(self._mjd[subscript]%period+n*period,
                        self._flux[subscript],yerr=self._err[subscript],
                        fmt='.k',zorder=-1,capsize=0)
                ax.scatter(self._mjd[subscript]%period+n*period,
                        self._flux[subscript],
                        c=clrs,edgecolor='k',zorder=100,cmap='gray',
                        vmin=0.0,vmax=1.0)
        else:
            ax.errorbar(self._mjd[subscript],self._flux[subscript],
                        yerr=self._err[subscript],
                        fmt='.k',zorder=-1,capsize=0)
            ax.scatter(self._mjd[subscript],self._flux[subscript],
                        c=clrs,edgecolor='k',zorder=100,cmap='gray',
                        vmin=0.0,vmax=1.0)

        mx,mn = max(self._flux[subscript]),min(self._flux[subscript])
        a = (mx-mn)/(mx+mn) # asymmetry

        if period is not None and calcspline and self._spline is None:
            data = np.zeros([len(self._mjd[subscript]),2])
            data[:,0] = self._flux[subscript]
            data[:,1] = self._err[subscript]
            self._spline = lyrae.fit(2*np.pi/period,self._mjd[subscript],data)[0]

        if period is not None and self._spline is not None:
            ts = np.linspace(0,nperiods*period,nperiods*500)
            ax.plot(ts,self._spline(ts),'k')

        model = self._calibpatch.get_model()
        ax.axhline(self._meanflux,color="k",ls="-")

        if period is not None:
            ax.set_xlim([0.0,nperiods*period])
        ax.set_ylim([0,2*self._meanflux])

        ax.set_ylabel(r'$\mathrm{nMgy}$',fontsize=16)
        ax.set_xlabel(r'$\mathrm{days}$',fontsize=16)

        if show_title:
            ax.set_title(r'%s / $\ln\,r_\mathrm{var} = %.3f$ / $a = %.3f$'%\
                (iau_name('Stripe82',*(self._calibpatch.get_radec(self._starnumber))),
                    self._lnoddsvar,a),fontsize=16)
        if hyperparams:
            # an = "\n".join(
            #         ["%10s = %.3e"%(k,getattr(model,k))
            #             for k in ["jitterrel2","jitterabs2","pvar","Q2",
            #                 "pbad","sigbad2"]])
            # an += "\n%10s = %.3f"%("lnrvar",self._lnoddsvar)
            # an += "\n%10s = %d"%("nstars",self._calibpatch.get_nstars())
            # if period is not None:
            #     an += "\n%10s = %.3f"%("period",period)
            an = "$%10s = %.3f$"%("a",a)
            if period is not None:
                an += "\n%10s = %.5f"%("period",period)
            ax.annotate(an, xy=(0.,0.),  xycoords='axes fraction',
                xytext=(0, 0), textcoords='offset points',
                size=9,family="monospace",
                bbox=dict(fc="w",alpha=0.25),alpha=0.5)
        return a


if __name__ == '__main__':
    pass

