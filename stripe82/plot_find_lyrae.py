#!/usr/bin/env python
# encoding: utf-8
"""
Plot the results of find_lyrae

NOTE: this should only be a temporary script... merge this with find lyrae eventually

History
-------
2011-08-16 - Created by Dan Foreman-Mackey

"""

__all__ = ['LyraeSearch']

import os
import os.path
import cPickle as pickle

import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

from calib_results import CalibrationPatch
from calibration import PhotoModel

class LyraeSearch:
    def __init__(self, bp):
        lsd = os.listdir(os.path.join(bp,'models'))
        nra = len(lsd)
        self._bp = bp
        self._patches = []
        self._periods = []
        self._ra  = np.zeros(nra)
        self._dec = np.zeros(nra)
        self._zeros = {}
        self._params = {}
        for k in ["jitterabs2","jitterrel2","pbad","pvar","sigbad2","Q2"]:
            self._params[k] = []
        count = 0
        for fn in lsd:
            number,ext = os.path.splitext(fn)
            if ext == '.pkl':
                print number
                data,vector,ra,dec,radius,period =\
                        pickle.load(open(os.path.join(bp,'models',fn),'rb'))
                model = PhotoModel(data,vector)
                for i in range(data.nobs):
                    obs = data.observations[i]
                    k = "%d%d"%(obs['run'],obs['camcol'])
                    if k not in self._zeros:
                        self._zeros[k] = np.zeros(nra)
                    self._zeros[k][count] = model.zero[i]
                self._periods.append(period)
                patch = CalibrationPatch(model,model.data,ra,dec,radius)
                self._patches.append(patch)
                self._ra[count]  = ra
                self._dec[count] = dec
                for k in list(self._params):
                    self._params[k].append(getattr(model,k))

                count += 1
        pickle.dump((self._ra,self._dec,self._zeros,self._params),
                open(os.path.join(self._bp,'cache.pkl'),'wb'),-1)

    def plot_model_params(self):
        for k in list(self._params):
            pl.figure()
            if k in ["pbad","pvar"]:
                pl.plot(self._ra,self._params[k],'.k')
                mu = np.median(self._params[k])
                pl.ylabel(k)
            else:
                pl.plot(self._ra,np.log(self._params[k]),'.k')
                mu = np.median(np.log(self._params[k]))
                pl.ylabel("$\ln$ "+k)
            pl.gca().axhline(mu,color="k",ls="--")
            pl.xlabel("R.A.")
            pl.title("%f"%mu)
            pl.savefig(os.path.join(self._bp,"%s.png"%k))

    def plot_lightcurves(self):
        constbp = os.path.join(self._bp,'plots','constant')
        varbp = os.path.join(self._bp,'plots','variable')
        nowrapbp = os.path.join(self._bp,'plots','nowrap')
        if not os.path.exists(constbp):
            os.makedirs(constbp)
        if not os.path.exists(varbp):
            os.makedirs(varbp)
        if not os.path.exists(nowrapbp):
            os.makedirs(nowrapbp)
        pl.figure()
        for number,patch in enumerate(self._patches):
            target_id,target_lc = patch.get_target()

            pl.clf()
            a = target_lc.plot(ax=pl.gca(),period=self._periods[number],nperiods=2,
                    calcspline=True,hyperparams=True)

            lnoddsvar = target_lc.get_lnoddsvar()
            if lnoddsvar > 0 and a > 0.25:
                pl.savefig(os.path.join(varbp,'%s.png'%number))
            else:
                pl.savefig(os.path.join(constbp,'%s.png'%number))

            pl.clf()
            a = target_lc.plot(ax=pl.gca(),period=None,nperiods=None,
                    calcspline=False,hyperparams=True)
            pl.savefig(os.path.join(nowrapbp,'%snowrap.png'%number))

    def plot_chunks(self):
        chunkbp = os.path.join(self._bp,'plots','chunks')
        if not os.path.exists(chunkbp):
            os.makedirs(chunkbp)

        # MAGIC: clusters in SDSS data by inspection
        chunks = [53500,53833,54250,54500]
        N = len(chunks)
        fig = pl.figure(figsize=(6,10))
        for number,patch in enumerate(self._patches):
            pl.clf()
            axes = [fig.add_subplot(N,1,i+1) for i in range(N)]
            target_id,target_lc = patch.get_target()
            for i in range(N):
                show_title = True
                if i > 0:
                    show_title = False
                    chunk = (target_lc.mjd >= chunks[i-1]) * \
                            (target_lc.mjd < chunks[i])
                else:
                    chunk = np.arange(len(target_lc.mjd))
                if np.sum(chunk) > 0:
                    target_lc._period = None
                    try:
                        target_lc.plot_slice(chunk,ax=axes[i],
                            nperiods=2,calcspline=True,show_title=show_title,
                            hyperparams=True)
                    except Exception as e:
                        print e
            pl.savefig(os.path.join(chunkbp,'%s.png'%number))

if __name__ == '__main__':
    import sys
    search = LyraeSearch(sys.argv[1])
    search.plot_chunks()
    sys.exit()
    search.plot_lightcurves()
    ra,dec,zero,params =\
            pickle.load(open(os.path.join(sys.argv[1],'cache.pkl'),'rb'))
    zerobp = os.path.join(sys.argv[1],'zeros')
    if not os.path.exists(zerobp):
        os.makedirs(zerobp)
    for k in zero:
        pl.clf()
        inds = zero[k] > 0
        pl.plot(ra[inds],zero[k][inds],'.k')
        mu = np.median(zero[k][inds])
        five = mu*0.05
        pl.gca().axhline(mu,color='k')
        pl.gca().axhline(mu+five,color='k',ls=':')
        pl.gca().axhline(mu-five,color='k',ls=':')
        pl.xlim(pl.gca().get_xlim()[::-1])
        pl.title("run: %s / camcol: %s"%(k[:-1],k[-1]))
        pl.xlabel("R.A.")
        pl.ylabel("ADU/nMgy")
        pl.savefig(os.path.join(zerobp,'%s.png'%(k)))

