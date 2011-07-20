#!/usr/bin/env python
# encoding: utf-8
"""


History
-------
2011-07-12 - Created by Dan Foreman-Mackey

"""


import numpy as np
import pylab as pl

import calibration.calibmodel as mod

cmodels = {}
rs = [1.5,5,15]
for r in rs:
    cmodels[str(r)] = mod.CalibrationModel({'radius': r})

cs = ['g','b','r']
for k in list(cmodels['15'].ra):
    pl.clf()
    ax1 = pl.subplot(111)
    #ax2 = pl.subplot(312)
    #ax3 = pl.subplot(313)
    for i,r in enumerate(rs):
        if k in cmodels[str(r)].ra:
            y = 10**(-np.array(cmodels[str(r)].zero[k])/2.5-9)
            ax1.plot(cmodels[str(r)].ra[k],y,'+%s'%(cs[i]),
                    label=str(r))
            # ax2.plot(cmodels[str(r)].ra[k],cmodels[str(r)].pars['jitterabs2'][k],
            #         '+%s'%(cs[i]))
            if k in cmodels[str(r)].splines:
                x = np.linspace(min(cmodels[str(r)].ra[k]),
                        max(cmodels[str(r)].ra[k]),500)
                decs = cmodels[str(r)].uniquedecs[k]
                dmean = np.median(decs)
                for j,spl in enumerate(cmodels[str(r)].splines[k]):
                    wdth = np.exp(-(decs[j]-dmean)**2/2/(7./60.0)**2)
                    y = spl['zero'](x)
                    ax1.plot(x,10**(-y/2.5-9),'-%s'%(cs[i]),lw=wdth)
                # ax2.plot(x,cmodels[str(r)].splines[k]['jitterabs2'](x),'-%s'%(cs[i]),lw=0.5)

    ax1.legend()
    ax1.set_ylim([0,600])
    ax1.set_ylabel(r'$\mathrm{ADU\,nMgy}^{-1}$',fontsize=16.)
    ax1.set_xlabel(r'$\mathrm{R.A.}$',fontsize=16)
    ax1.set_title(r'$\mathrm{run} = %d \,\, \mathrm{camcol} = %d$'%\
            (int(k[:-1]),int(k[-1])),fontsize=16)
    pl.savefig('../calib_grid/%s.png'%k)


