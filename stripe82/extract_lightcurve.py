#!/usr/bin/env python
# encoding: utf-8
"""
Command line utility for extracting a light-curve from Stripe 82 at a given RA/Dec

History
-------
2011-08-03 - Created by Dan Foreman-Mackey

"""

__all__ = ['extract_lightcurve','extract_lightcurves']

import os
import cPickle as pickle
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

import numpy as np
import numpy.ma as ma

from calibration import calibrate,PhotoModel,odds_variable,odds_bad
from calibration.opt import *
from lyrae import sesar,fit,find_period

def hogg_expand(y,scale=1.1):
    """
    Expand plot limits by SOME percent

    Parameters
    ----------
    y : list
        The data

    Optional
    --------
    scale : float (default : 1.1)
        How much to scale by

    Returns
    -------
    (xmin,xmax) : tuple
        The limits

    History
    -------
    2011-08-15 - Created by Dan Foreman-Mackey and David W. Hogg

    """
    mx,mn = max(y),min(y)
    a = 0.5*(mx+mn)
    b = 0.5*(mx-mn)
    return (a-scale*b,a+scale*b)

def _angular_distance(ra1,dec1,ra2,dec2):
    ra1,dec1 = np.radians(ra1),np.radians(dec1)
    ra2,dec2 = np.radians(ra2),np.radians(dec2)
    crho = np.cos(dec1)*np.cos(dec2)\
            *(np.cos(ra1)*np.cos(ra2)+np.sin(ra1)*np.sin(ra2))\
            +np.sin(dec1)*np.sin(dec2)
    return np.degrees(np.arccos(crho))

def extract_lightcurves(*args,**kwargs):
    if 'units' in kwargs:
        units = kwargs['units']
    else:
        units = 1.0
    if 'model' in kwargs:
        model = kwargs['model']
    else:
        (ra,dec,radius) = args
        model = calibrate(ra,dec,radius=radius)
    mjd = model.data.mjd()
    flux = model.data.flux/model.zero[:,np.newaxis]*units
    mask = model.data.ivar <= 1e-8
    inv = ma.array(model.data.ivar,mask=mask)
    err = units/np.sqrt(inv*model.zero[:,np.newaxis]**2)
    return mjd,flux,err,model

def extract_lightcurve(ra,dec,radius):
    mjd,flux,err,model = extract_lightcurves(ra,dec,radius)
    dist = lambda i: _angular_distance(ra,dec,
            model.data.stars[i]['ra'],model.data.stars[i]['dec'])
    target_id = sorted(range(len(model.data.stars)),key = dist)[0]
    inds = ~err.mask[:,target_id]
    return mjd[inds],flux[inds,target_id],err[inds,target_id],model

def plot_lightcurve(i,mjd,flux,err,model,ax=None,badodds=None,period=None,
        nperiods=1,fit_period=False,pvar=None):
    if ax is None:
        ax = pl.gca()
    if badodds is None:
        badodds = odds_bad(model,model.data)
    if period is None and not fit_period:
        period = max(mjd)+1

    inds = ~err.mask[:,i]

    if fit_period:
        data = np.zeros([np.sum(inds),2])
        data[:,0] = flux[inds,i]
        data[:,1] = err[inds,i]
        period = find_period(mjd[inds],data)
        m,chi2 = fit(2*np.pi/period,mjd[inds],data)

        ts = np.linspace(0,nperiods*period,nperiods*500)
        ax.plot(ts,m(ts),'k')

    err[inds,i] = np.sqrt(err[inds,i]**2+model.jitterabs2+model.flux[i]**2*model.jitterrel2)

    for n in range(nperiods):
        ax.errorbar(mjd[inds]%period+n*period,flux[inds,i],yerr=err[inds,i],
                fmt='.k',zorder=-1,capsize=0)

    # colors based on r_bad
    clrs = logodds2prob(badodds[inds,i])

    for n in range(nperiods):
        ax.scatter(mjd[inds]%period+n*period,flux[inds,i],
                c=clrs,edgecolor='k',zorder=100,cmap='gray',
                vmin=0.0,vmax=1.0)

    ax.axhline(model.flux[i],color='k',ls='--')

    ax.set_ylim([0,2*model.flux[i]])
    ax.set_xlim([0,nperiods*period])

    ax.set_ylabel(r'$\mathrm{nMgy}$',fontsize=16)
    ax.set_xlabel(r'$\mathrm{days}$',fontsize=16)

    #if i == target_id:
    #    title = r"Target: "
    #else:
    #    title = r""
    title = r""
    title += r"$N_\mathrm{obs} = %d"%(np.sum(inds))
    if pvar is not None:
        title += ",\,\ln\,p_\mathrm{var} = %.3f$"%(np.log(pvar[i]))
    ax.set_title(title,fontsize=16)

    return period

def plot_lightcurves(basepath,mjd,flux,err,model,
        badodds=None,period=None,s_data=None,pvar=None):
    if badodds is None:
        badodds = odds_bad(model,model.data)

    pl.figure()
    for i in range(np.shape(flux)[1]):
        pl.clf()
        ax = pl.subplot(111)

        if 0: #s_data is not None:
            plot_lightcurve(i,mjd,flux,err,model,ax=ax,badodds=badodds,period=period,
                    nperiods=2,pvar=pvar)
            ax.plot(s_time%period,s_data,'og',alpha=0.3)
            ax.plot(s_time%period+period,s_data,'og',alpha=0.3)
            ax.set_xlim([0,2*period])
        else:
            plot_lightcurve(i,mjd,flux,err,model,ax=ax,badodds=badodds,period=period,
                    nperiods=2,pvar=pvar)


        #ax = pl.subplot(212)
        #ax.plot(mjd_i%period,badodds[inds,i],'.k')
        #if s_data is not None:
        #    ax.plot(mjd_i%period+period,badodds[inds,i],'.k')

        #ax.set_ylabel(r'$\ln\,r^\mathrm{bad}_{i\alpha}$',fontsize=16)
        #ax.set_xlabel(r'$t\%T$',fontsize=16)

        pl.savefig(os.path.join(basepath,'%04d.png'%i))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract, calibrate and plot a lightcurve')
    parser.add_argument('basepath',help='Directory for results')
    parser.add_argument('-f','--tmpfile',help='Temp pickle file',
                        default=None)
    parser.add_argument('--all',
                        help='all figures',
                        action='store_true')
    parser.add_argument('--debug',
                        help='only plot debug plots',
                        action='store_true')
    parser.add_argument('--period',
                        help='plot target and fit period',
                        action='store_true')
    parser.add_argument('--ra',default=29.47942)
    parser.add_argument('--dec',default=0.383557)
    parser.add_argument('--radius',default=5.0)

    args = parser.parse_args()
    bp = str(args.basepath)
    if not os.path.exists(bp):
        os.makedirs(bp)
    if args.tmpfile is not None:
        tmpfile = os.path.join(bp,args.tmpfile)
    else:
        tmpfile = None

    ra,dec = float(args.ra),float(args.dec)#29.47942,0.383557 #10.0018734334081,0.791580301596976
    radius = float(args.radius)

    if np.any(sesar.coords['ra'] == ra):
        ind = sesar.coords[sesar.coords['ra'] == ra]['Num']
        period = sesar.coords[sesar.coords['ra'] == ra]['Per']
        print "Sesar's period: ",period
        sesardata = sesar.table1['%d'%(ind)][...]
        inds = sesardata['g'] > 0
        s_time = sesardata['gmjd'][inds]
        s_data = mag2nmgy(sesardata['g'][inds])
        #s_err  = mag2nmgy(sesardata['g'][inds]+sesardata['gerr'][inds]) - s_data #+\
        #s_err  = s_data-mag2nmgy(sesardata['g'][inds]-sesardata['gerr'][inds])
    else:
        s_data = None
        period = None

    if tmpfile is not None and os.path.exists(tmpfile):
        model = PhotoModel(*pickle.load(open(tmpfile,'rb')))
        (mjd,flux,err,model) = extract_lightcurves(model=model)
    else:
        mjd,flux,err,model = extract_lightcurves(ra,dec,radius)
    if args.tmpfile is not None:
            pickle.dump((model.data,model.vector()),open(tmpfile,'wb'),-1)
    print "final model= ",model

    if period is None:
        period = max(mjd)+1
    varodds = odds_variable(model,model.data)
    pvar = logodds2prob(varodds)

    badodds = odds_bad(model,model.data)
    sorted_inds = np.argsort(model.flux)
    sids = np.arange(len(model.flux))
    dist = lambda i: _angular_distance(ra,dec,
            model.data.stars[i]['ra'],model.data.stars[i]['dec'])
    target_id = sorted(range(len(sorted_inds)),key = dist)[0]
    target = model.data.stars[target_id]

    if args.all or args.period:
        ax = pl.figure().add_subplot(111)
        period = plot_lightcurve(target_id,mjd,flux,err,model,ax=ax,badodds=badodds,
                nperiods=2,fit_period=True,pvar=pvar)
        if s_data is not None:
            for i in range(2):
                ax.plot(s_time%period+i*period,s_data,'s',
                        color=[0.5,0.5,0.5],alpha=0.5)
        pl.savefig(os.path.join(bp,'target.png'))

    if args.all or args.debug:
        # plot 1
        pl.figure(figsize=(8.,8.))

        ax1 = pl.subplot(311)
        ax1.plot(model.zero,'.k')
        ax1.set_xticklabels([])
        ax1.set_ylabel(r'$(\mathrm{nMgy/ADU})_i$',fontsize=16.)

        Nstars = np.sum(model.data.ivar>0,axis=1)
        ax2 = pl.subplot(312)
        nbad = np.sum(badodds > 0,axis=1)
        ax2.plot(nbad,'.k')
        ax2.set_xticklabels([])
        ax2.set_ylim(hogg_expand(nbad))
        ax2.set_ylabel(r'$N_\mathrm{bad}$',fontsize=16.)

        ax3 = pl.subplot(313)
        ax3.plot(Nstars,'.k')
        ax3.set_ylim(hogg_expand(Nstars))
        ax3.set_ylabel(r'$N_\mathrm{stars}$',fontsize=16.)
        ax3.set_xlabel(r'$\mathrm{Obs.ID}$',fontsize=16.)

        pl.savefig(os.path.join(bp,'plot1.png'))


        # plot 2
        pl.figure(figsize=(8.,8.))

        taget_star = sids[sorted_inds == target_id]
        ax1 = pl.subplot(411)
        ax1.plot(sids,model.mag[sorted_inds],'.k')
        ax1.axvline(taget_star,color='r')
        ax1.set_xticklabels([])
        ax1.set_ylim(hogg_expand(model.mag)[::-1])
        ax1.set_ylabel(r'$g-\mathrm{mag}$',fontsize=16.)

        ax2 = pl.subplot(412)
        Nobs = np.sum(model.data.ivar>0,axis=0)
        ax2.plot(sids,Nobs[sorted_inds],'.k')
        ax2.axvline(taget_star,color='r')
        ax2.set_xticklabels([])
        ax2.set_ylim(hogg_expand(Nobs))
        ax2.set_ylabel(r'$N_\mathrm{obs}$',fontsize=16.)

        ax3 = pl.subplot(413)
        ax3.plot(sids,pvar[sorted_inds],'.k')
        ax3.axhline(0., color='k', alpha=0.5)
        ax3.axhline(1., color='k', alpha=0.5)
        ax3.set_ylim(hogg_expand([0., 1.]))
        ax3.axvline(taget_star,color='r')
        ax3.set_xticklabels([])
        ax3.set_ylabel(r'$p_\mathrm{var}$',fontsize=16.)

        ax4 = pl.subplot(414)
        nbad = np.sum(badodds > 0,axis=0)
        ax4.plot(sids,nbad[sorted_inds],'.k')
        ax4.axvline(taget_star,color='r')
        ax4.set_xlim(ax1.get_xlim())
        ax4.set_ylim(hogg_expand(nbad))
        ax4.set_ylabel(r'$N_\mathrm{bad}$',fontsize=16.)
        ax4.set_xlabel(r'$\mathrm{Star\,ID}$',fontsize=16.)

        pl.savefig(os.path.join(bp,'plot2.png'))

        posbp = os.path.join(bp,'offsets')
        if not os.path.exists(posbp):
            os.makedirs(posbp)
        pl.figure(figsize=(8.,8.))
        for j in range(np.shape(model.data.flux)[0]):
            pl.clf()
            ax = pl.subplot(211,aspect='equal')

            # plot 3
            ras,decs = [],[]
            dr,dd = [],[]
            targetradec = None
            for i,doc in enumerate(model.data.stars):
                ras.append(doc['ra'])
                decs.append(doc['dec'])
                photo = model.data.data['model'][j,i]
                if photo[1] > 0:
                    dr.append(photo[2]/photo[1])
                    dd.append(photo[3]/photo[1])
                else:
                    dr.append(0)
                    dd.append(0)
                if i == target_id:
                    targetradec = (doc['ra'],doc['dec'])
            dr,dd = np.array(dr),np.array(dd)
            ax.quiver(ras,decs,dr,dd)
            ax.plot(targetradec[0],targetradec[1],'or')
            ax.set_xlabel('R.A.')
            ax.set_ylabel('Dec.')

            ax = pl.subplot(212)

            length = np.sqrt(dr**2+dd**2)
            ax.plot(model.flux,length,'.k')
            ax.set_xlabel('nMgy')
            ax.set_ylabel('length of offset (pixels)')
            pl.savefig(os.path.join(posbp,'%04d.png'%j))

        starbp = os.path.join(bp,'stars')
        if not os.path.exists(starbp):
            os.makedirs(starbp)

        plot_lightcurves(starbp,mjd,flux,err,model,
                badodds=badodds,period=period,s_data=s_data,pvar=pvar)

#        pl.figure(figsize=(8.,8.))
#        for sid,i in enumerate(sorted_inds):
#            inv = model.data.ivar
#            inds = inv[:,i] > 0
#            err0_i = 1e9/np.sqrt(inv[inds,i]*model.zero[inds]**2)
#            flux_i = flux[inds,i]
#            mjd_i = mjd[inds]
#
#            if len(mjd_i)>2:
#                pl.clf()
#
#                # colors based on r_bad
#                clrs = badodds[inds,i]
#                clrs -= min(clrs)
#                clrs /= max(clrs)/256.0
#
#                #pl.errorbar(mjd_i%period,flux_i,yerr=err_i,fmt='.k',alpha=0.5)
#                ax = pl.subplot(211)
#                ax.errorbar(mjd_i%period,flux_i,yerr=err0_i,fmt='.k',zorder=-1)
#                ax.scatter(mjd_i%period,flux_i,s=40,c=clrs,edgecolor='none',zorder=100)
#                #ax.set_xlim([0,period])
#
#                if s_data is not None:
#                    ax.errorbar(mjd_i%period+period,flux_i,yerr=err0_i,fmt='.k')
#                    ax.scatter(mjd_i%period+period,flux_i,s=40,c=clrs,edgecolor='none',zorder=100)
#                    pl.plot(s_time%period,s_data,'og',alpha=0.3)
#                    pl.plot(s_time%period+period,s_data,'og',alpha=0.3)
#                    ax.set_xlim([0,2*period])
#
#
#                ax.axhline(model.flux[i]*1e9,color='r',ls='--')
#
#                #if model.flux[i] < 2*np.median(model.flux):
#                #    ax.set_ylim([0,2*np.median(model.flux)*1e9])
#                #else:
#                ax.set_ylim([0,2*model.flux[i]*1e9])
#
#                ax.set_ylabel(r'$\mathrm{nMgy}$',fontsize=16)
#
#                if i == target_id:
#                    ax.set_title(
#                            r"Target: $N_\mathrm{obs} = %d,\,\ln\,r^\mathrm{var}_{\alpha} = %.3f$"%\
#                                (np.sum(inds),varodds[i]),fontsize=16)
#                else:
#                    ax.set_title(
#                        r"$N_\mathrm{obs} = %d,\,\ln\,r^\mathrm{var}_{\alpha} = %.3f$"%\
#                                (np.sum(inds),varodds[i]),fontsize=16)
#
#                ax = pl.subplot(212)
#                ax.plot(mjd_i%period,badodds[inds,i],'.k')
#                if s_data is not None:
#                    ax.plot(mjd_i%period+period,badodds[inds,i],'.k')
#
#                ax.set_ylabel(r'$\ln\,r^\mathrm{bad}_{i\alpha}$',fontsize=16)
#                ax.set_xlabel(r'$t\%T$',fontsize=16)
#
#                pl.savefig(os.path.join(starbp,'%04d.png'%sid))


#        obsbp = os.path.join(bp,'obs')
#        if not os.path.exists(obsbp):
#            os.makedirs(obsbp)
#
#        for oid in range(np.shape(flux)[0]):
#            pl.clf()
#            ax2 = pl.subplot(312)
#            ax2.plot(sids,(badodds[:,])[sorted_inds],'.k')
#
#            ax2.axvline(sorted_inds[target_id],color='r')
#            ax2.set_xticklabels([])
#            ax2.set_ylabel(r'$f^*_\alpha\,[\mathrm{nMgy}]$',fontsize=16.)
#
            #pl.clf()
            ## if s_data is not None:
            ##     pl.plot(s_time%period,s_data,'og',alpha=0.3)
            ##     pl.plot(s_time%period+period,s_data,'og',alpha=0.3)
            #
            ## colors based on r_bad
            #clrs = badodds[inds,i]
            #clrs -= min(clrs)
            #clrs /= max(clrs)/256.0

            ##pl.errorbar(mjd_i%period,flux_i,yerr=err_i,fmt='.k',alpha=0.5)
            #ax = pl.subplot(211)
            #ax.errorbar(mjd_i%period,flux_i,yerr=err0_i,fmt='.k',zorder=-1)
            #ax.scatter(mjd_i%period,flux_i,s=40,c=clrs,edgecolor='none',zorder=100)
            #ax.set_xlim([0,period])
            #
            #if s_data is not None:
            #    ax.errorbar(mjd_i%period+period,flux_i,yerr=err0_i,fmt='.k')
            #    ax.scatter(mjd_i%period+period,flux_i,s=40,c=clrs,edgecolor='none',zorder=100)
            #    ax.set_xlim([0,2*period])
            #

            #ax.axhline(model.flux[i]*1e9,color='r',ls='--')

            #if model.flux[i] < 2*np.median(model.flux):
            #    ax.set_ylim([0,2*np.median(model.flux)*1e9])
            #else:
            #    ax.set_ylim([0,2*model.flux[i]*1e9])

            #ax.set_ylabel(r'$\mathrm{nMgy}$',fontsize=16)

            #if i == target_id:
            #    ax.set_title(
            #            r"Target: $N_\mathrm{obs} = %d,\,\ln\,r^\mathrm{var}_{\alpha} = %.3f$"%\
            #                (np.sum(inds),varodds[i]),fontsize=16)
            #else:
            #    ax.set_title(
            #        r"$N_\mathrm{obs} = %d,\,\ln\,r^\mathrm{var}_{\alpha} = %.3f$"%\
            #                (np.sum(inds),varodds[i]),fontsize=16)

            #ax = pl.subplot(212)
            #ax.plot(mjd_i%period,badodds[inds,i],'.k')
            #if s_data is not None:
            #    ax.plot(mjd_i%period+period,badodds[inds,i],'.k')

            #ax.set_ylabel(r'$\ln\,r^\mathrm{bad}_{i\alpha}$',fontsize=16)
            #ax.set_xlabel(r'$t\%T$',fontsize=16)

            #pl.savefig(os.path.join(obsbp,'%04d.png'%oid))


