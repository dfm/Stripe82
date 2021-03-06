import os

import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as pl

import lyrae
import sesar

pl.figure(figsize=(8, 12))

fn = "sesar.out"
f = open(fn, "w")
f.write("# i, id, T, order, chi2, A\n")
f.close()

for ind in [144]: #range(len(sesar.table2)):
    for order in [12]: #[3, 6, 9, 12]:
        _id = str(sesar.table2[ind]["Num"])
        print "%d:"%ind, _id, order
        d = sesar.table1[_id][...]

        bands = "gri"
        inds = dict([(b, (d[b] > 1) * (d[b] < 99)) for b in bands])
        time = dict([(b, d[b+"mjd"][inds[b]]) for b in bands])
        flux = dict([(b, 10**(9-0.4*d[b][inds[b]])) for b in bands])
        ferr = dict([(b, flux[b] * np.log(10)*0.4 *d[b+"err"][inds[b]])
            for b in bands])

        p = [sesar.table2[ind]["Per"], lyrae.find_period(time, flux, ferr,
            order=order)]

        m = dict([(b, [lyrae.get_model(p0, time[b], flux[b], ferr[b],
            order=order) for p0 in p]) for b in time])

        f = open(fn, "a")
        s = ["%4d"%ind, _id]
        s += ["%.8f"%p[1]]
        s += [str(order)]
        s += ["%.0f"%m[b][1][2] for b in time]
        for b in m:
            a = m[b][1][1]
            s += ["%e"%a0 for a0 in a]
        f.write(" ".join(s)+"\n")
        f.close()

        def plot_lc(b, i):
            N = 5413
            t = np.linspace(min(time[b]), max(time[b]), N)
            pl.plot(t%(2*p[i]), m[b][(i+1)%2][0](t), "k.", alpha=0.1)

            t = np.linspace(0, 2*p[i], N)
            pl.plot(t, m[b][i][0](t), 'r')

            pl.errorbar(time[b]%p[i], flux[b], yerr=ferr[b], fmt=".k")
            pl.errorbar(time[b]%p[i]+p[i], flux[b], yerr=ferr[b], fmt=".k")
            pl.annotate(r"$\chi^2 = %.0f$"%m[b][i][2], [1,1], xytext=[-3, -3],
                    xycoords="axes fraction", textcoords="offset points",
                    horizontalalignment="right", verticalalignment="top",
                    backgroundcolor=[1,1,1,0.8])
            pl.annotate(r"$T = %.6f$"%p[i], [1,1], xytext=[-3, -3-14-3],
                    xycoords="axes fraction", textcoords="offset points",
                    horizontalalignment="right", verticalalignment="top",
                    backgroundcolor=[1,1,1,0.8])

        pl.clf()
        t = r"$\mathrm{Sesar\,%s}:\,T_0=%.6f,\,T_1=%.6f,"%(_id,p[0],p[1])
        t += "\,\log_{10}\Delta T=%.2f$"%(np.log10(np.abs(p[0]-p[1])))
        pl.title(t, fontsize=14.)

        for i,b in enumerate(bands):
            # Plot Sesar's results
            ax1 = pl.subplot(len(bands), 2, 2*i+1)

            if i == 0:
                t = r"$\mathrm{Sesar\,%s}:\,T_0=%.6f$"%(_id,p[0])
                pl.title(t, fontsize=14.)

            plot_lc(b, 0)
            pl.xlim(0, 2*p[0])
            pl.ylabel(r"$f_%s$"%b, fontsize=14.)

            if i < len(bands)-1:
                pl.gca().set_xticklabels([])
            else:
                pl.xlabel(r"$t/T$", fontsize=14.)

            # Plot my results
            ax2 = pl.subplot(len(bands), 2, 2*i+2)

            if i == 0:
                t = r"$\mathrm{DFM}:\,T_1=%.6f,"%(p[1])
                t += "\,\log_{10}\Delta T=%.2f$"%(np.log10(np.abs(p[0]-p[1])))
                pl.title(t, fontsize=14.)

            plot_lc(b, 1)
            pl.xlim(0, 2*p[1])
            pl.gca().set_yticklabels([])

            r1, r2 = ax1.get_ylim(), ax2.get_ylim()
            ax1.set_ylim([min(r1[0], r2[0]), max(r1[1], r2[1])])
            ax2.set_ylim([min(r1[0], r2[0]), max(r1[1], r2[1])])

            if i < len(bands)-1:
                pl.gca().set_xticklabels([])
            else:
                pl.xlabel(r"$t/T$", fontsize=14.)

        try:
            os.makedirs("sesar_test/%d"%order)
        except os.error:
            pass
        pl.savefig("sesar_test/%d/sesar-%04d.png"%(order, ind))

