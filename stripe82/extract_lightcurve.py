#!/usr/bin/env python
# encoding: utf-8
"""
Command line utility for extracting a light-curve from Stripe 82 at a given RA/Dec

History
-------
2011-08-03 - Created by Dan Foreman-Mackey

"""

__all__ = ['hogg_expand']

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
from calib_results import *
from lyrae import sesar,fit,find_period
from calibration.extract import extract_lightcurves

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
    mx,mn = np.max(y),np.min(y)
    a = 0.5*(mx+mn)
    b = 0.5*(mx-mn)
    return (a-scale*b,a+scale*b)

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

    patch = CalibrationPatch(model,model.data,ra,dec,radius)

    varodds = patch.lnoddsvar
    pvar = logodds2prob(varodds)
    badodds = patch.lnoddsbad
    target_id,target_lc = patch.get_target()
    period = target_lc.get_period()

    #dist = lambda i: _angular_distance(ra,dec,
    #        model.data.stars[i]['ra'],model.data.stars[i]['dec'])
    #target_id = sorted(range(len(sorted_inds)),key = dist)[0]
    #target = model.data.stars[target_id]

    if args.all or args.period:
        ax = pl.figure().add_subplot(111)
        if s_data is not None:
            for i in range(2):
                ax.plot(s_time%period+i*period,s_data,'s',
                        color=[0.5,0.5,0.5],alpha=0.5)

        target_lc.plot(ax=ax,nperiods=2,hyperparams=True)
        pl.savefig(os.path.join(bp,'target.png'))

    if args.all or args.debug:
        # plot 1
        pl.figure(figsize=(8.,8.))
        Nstars = np.sum(model.data.ivar>0,axis=1)
        nbad = np.sum(badodds > 0,axis=1)

        ax1 = pl.subplot(311)
        ax1.plot(model.zero,'.k')
        ax1.set_xticklabels([])
        ax1.set_ylabel(r'$(\mathrm{nMgy/ADU})_i$',fontsize=16.)

        ax2 = pl.subplot(312)
        ax2.plot(nbad,'.k')
        ax2.set_xticklabels([])
        ax2.set_ylim(hogg_expand(nbad))
        ax2.set_xlim(ax1.get_xlim())
        ax2.set_ylabel(r'$N_\mathrm{bad}$',fontsize=16.)

        ax3 = pl.subplot(313)
        ax3.plot(Nstars,'.k')
        ax3.set_ylim(hogg_expand(Nstars))
        ax3.set_xlim(ax1.get_xlim())
        ax3.set_ylabel(r'$N_\mathrm{stars}$',fontsize=16.)
        ax3.set_xlabel(r'$\mathrm{Run\,ID}$',fontsize=16.)

        pl.savefig(os.path.join(bp,'plot1.png'))

        # plot 2
        pl.figure(figsize=(8.,8.))

        sorted_inds = np.argsort(model.flux)
        sids = np.arange(len(model.flux))
        taget_star = sids[sorted_inds == target_id]
        nobs = np.sum(model.data.ivar>0,axis=0)[sorted_inds]
        nbad = np.sum(badodds > 0,axis=0)[sorted_inds]
        ngood = nobs-nbad

        ax1 = pl.subplot(311)
        ax1.plot(sids,model.mag[sorted_inds],'.k')
        ax1.set_xticklabels([])
        ax1.set_ylim(hogg_expand(model.mag)[::-1])
        ax1.set_ylabel(r'$g-\mathrm{mag}$',fontsize=16.)

        ax2 = pl.subplot(312)
        ax2.plot(sids,nobs,'.k',alpha=0.6)
        ax2.plot(sids,ngood,'.k')
        ax2.set_xticklabels([])
        ax2.set_ylim(hogg_expand(np.append(nobs,ngood)))
        ax2.set_ylabel(r'$N_\mathrm{obs}$',fontsize=16.)

        ax3 = pl.subplot(313)
        ax3.plot(sids,pvar[sorted_inds],'.k')
        ax3.axhline(0., color='k', alpha=0.5)
        ax3.axhline(1., color='k', alpha=0.5)
        ax3.set_ylim(hogg_expand([0., 1.]))
        ax3.set_ylabel(r'$p_\mathrm{var}$',fontsize=16.)
        ax3.set_xlabel(r'$\mathrm{Star\,ID}$',fontsize=16.)

        for ax in [ax1,ax2,ax3]:
            ax.axvline(taget_star,color='r')

        pl.savefig(os.path.join(bp,'plot2.png'))

        # offsets
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

        # stars
        starbp = os.path.join(bp,'stars')
        if not os.path.exists(starbp):
            os.makedirs(starbp)

        figs = patch.plot(period=period,nperiods=2)
        for i,fig in enumerate(figs):
            fig.savefig(os.path.join(starbp,'%04d.png'%i))

        #plot_lightcurves(starbp,mjd,flux,err,model,
        #        badodds=badodds,period=period,s_data=s_data,pvar=pvar)

