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
import argparse

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
    parser = argparse.ArgumentParser(description='Fit M31 dynamics')
    parser.add_argument('basepath',help='Directory for results')
    parser.add_argument('-f','--tmpfile',help='Temp pickle file',
                        default=None)
    parser.add_argument('--all',
                        help='all figures',
                        action='store_true')
    parser.add_argument('--debug',
                        help='only plot debug plots',
                        action='store_true')
    parser.add_argument('-r','-a','--ra',default=29.47942)
    parser.add_argument('-d','--dec',default=0.383557)

    args = parser.parse_args()
    bp = str(args.basepath)
    if not os.path.exists(bp):
        os.makedirs(bp)

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
    flux = model.data.flux/model.zero[:,np.newaxis]*1e9
    varodds = odds_variable(model,model.data)
    badodds = odds_bad(model,model.data)
    sorted_inds = np.argsort(model.flux)
    sids = np.arange(len(model.flux))
    dist = lambda i: _angular_distance(ra,dec,
            model.data.stars[i]['ra'],model.data.stars[i]['dec'])
    target_id = sorted(range(len(sorted_inds)),key = dist)[0]
    target = model.data.stars[target_id]

    if args.all or args.debug:
        # plot 1
        pl.figure(figsize=(8.,8.))

        ax1 = pl.subplot(311)
        ax1.plot(model.zero/1e9,'.k')
        ax1.set_xticklabels([])
        ax1.set_ylabel(r'$(\mathrm{nMgy/ADU})_i$',fontsize=16.)

        Nstars = np.sum(model.data.ivar>0,axis=1)
        ax2 = pl.subplot(312)
        ax2.plot(np.sum(badodds,axis=1),'.k')
        ax2.set_xticklabels([])
        ax2.set_ylabel(r'$\sum_\alpha \ln \, r^\mathrm{bad}_{i\alpha}$',fontsize=16.)

        ax3 = pl.subplot(313)
        ax3.plot(Nstars,'.k')
        ax3.set_ylabel(r'$N_\mathrm{stars}$',fontsize=16.)
        ax3.set_xlabel(r'$\mathrm{Obs.ID}$',fontsize=16.)

        pl.savefig(os.path.join(bp,'plot1.png'))


        # plot 2
        pl.figure(figsize=(8.,8.))

        ax1 = pl.subplot(411)
        ax1.plot(sids,model.flux[sorted_inds]*1e9,'.k')
        ax1.axvline(sorted_inds[target_id],color='r')
        ax1.set_xticklabels([])
        ax1.set_ylabel(r'$f^*_\alpha\,[\mathrm{nMgy}]$',fontsize=16.)

        ax2 = pl.subplot(412)
        Nobs = np.sum(model.data.ivar>0,axis=0)
        ax2.plot(sids,Nobs[sorted_inds],'.k')
        ax2.axvline(sorted_inds[target_id],color='r')
        ax2.set_xticklabels([])
        ax2.set_ylabel(r'$N_\mathrm{obs}$',fontsize=16.)

        ax3 = pl.subplot(413)
        ax3.plot(sids,varodds[sorted_inds],'.k')
        ax3.axvline(sorted_inds[target_id],color='r')
        ax3.set_xticklabels([])
        ax3.set_ylabel(r'$\ln \, r^\mathrm{var}_\alpha$',fontsize=16.)

        ax4 = pl.subplot(414)
        ax4.plot(sids,(np.sum(badodds,axis=0))[sorted_inds],'.k')
        ax4.axvline(sorted_inds[target_id],color='r')
        ax4.set_ylabel(r'$\sum_i \ln \, r^\mathrm{bad}_{i\alpha}$',fontsize=16.)
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
            ax.plot(model.flux*1e9,length,'.k')
            ax.set_xlabel('nMgy')
            ax.set_ylabel('length of offset (pixels)')
            pl.savefig(os.path.join(posbp,'%04d.png'%j))

        starbp = os.path.join(bp,'stars')
        if not os.path.exists(starbp):
            os.makedirs(starbp)

        pl.figure(figsize=(8.,8.))
        for sid,i in enumerate(sorted_inds):
            inv = model.data.ivar
            inds = inv[:,i] > 0
            err0_i = 1e9/np.sqrt(inv[inds,i]*model.zero[inds]**2)
            flux_i = flux[inds,i]
            mjd_i = mjd[inds]

            pl.clf()
            # if s_data is not None:
            #     pl.plot(s_time%period,s_data,'og',alpha=0.3)
            #     pl.plot(s_time%period+period,s_data,'og',alpha=0.3)
            
            # colors based on r_bad
            clrs = badodds[inds,i]
            clrs -= min(clrs)
            clrs /= max(clrs)/256.0

            #pl.errorbar(mjd_i%period,flux_i,yerr=err_i,fmt='.k',alpha=0.5)
            ax = pl.subplot(211)
            ax.errorbar(mjd_i%period,flux_i,yerr=err0_i,fmt='.k',zorder=-1)
            ax.scatter(mjd_i%period,flux_i,s=40,c=clrs,edgecolor='none',zorder=100)
            ax.set_xlim([0,period])
            
            if s_data is not None:
                ax.errorbar(mjd_i%period+period,flux_i,yerr=err0_i,fmt='.k')
                ax.scatter(mjd_i%period+period,flux_i,s=40,c=clrs,edgecolor='none',zorder=100)
                ax.set_xlim([0,2*period])
                

            ax.axhline(model.flux[i]*1e9,color='r',ls='--')

            if model.flux[i] < 2*np.median(model.flux):
                ax.set_ylim([0,2*np.median(model.flux)*1e9])
            else:
                ax.set_ylim([0,2*model.flux[i]*1e9])

            ax.set_ylabel(r'$\mathrm{nMgy}$',fontsize=16)

            if i == target_id:
                ax.set_title(
                        r"Target: $N_\mathrm{obs} = %d,\,\ln\,r^\mathrm{var}_{\alpha} = %.3f$"%\
                            (np.sum(inds),varodds[i]),fontsize=16)
            else:
                ax.set_title(
                    r"$N_\mathrm{obs} = %d,\,\ln\,r^\mathrm{var}_{\alpha} = %.3f$"%\
                            (np.sum(inds),varodds[i]),fontsize=16)

            ax = pl.subplot(212)
            ax.plot(mjd_i%period,badodds[inds,i],'.k')
            if s_data is not None:
                ax.plot(mjd_i%period+period,badodds[inds,i],'.k')

            ax.set_ylabel(r'$\ln\,r^\mathrm{bad}_{i\alpha}$',fontsize=16)
            ax.set_xlabel(r'$t\%T$',fontsize=16)

            pl.savefig(os.path.join(starbp,'%04d.png'%sid))


        obsbp = os.path.join(bp,'obs')
        if not os.path.exists(obsbp):
            os.makedirs(obsbp)

        for oid in range(np.shape(flux)[0]):
            inv = model.data.ivar
            inds = inv[oid,:] > 0
            err0_i = 1e9/np.sqrt(inv[oid,inds]*model.zero[inds]**2)
            flux_i = flux[inds,i]
            mjd_i = mjd[inds]
            
            pl.clf()
            ax2 = pl.subplot(312)
            ax2.plot(sids,(badodds[:,])[sorted_inds],'.k')

            ax2.axvline(sorted_inds[target_id],color='r')
            ax2.set_xticklabels([])
            ax2.set_ylabel(r'$f^*_\alpha\,[\mathrm{nMgy}]$',fontsize=16.)
            
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
        

