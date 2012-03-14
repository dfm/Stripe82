import numpy as np
import matplotlib.pyplot as pl

import lyrae
import sesar

pl.figure(figsize=(8,8))

for ind in [207]: #range(91, len(sesar.table2)):
    _id = str(sesar.table2[ind]["Num"])
    print "%d:"%ind, _id
    d = sesar.table1[_id][...]

    bands = "gri"
    inds = dict([(b, (d[b] > 1) * (d[b] < 99)) for b in bands])
    time = dict([(b, d[b+"mjd"][inds[b]]) for b in bands])
    flux = dict([(b, 10**(9-0.4*d[b][inds[b]])) for b in bands])
    ferr = dict([(b, flux[b] * np.log(10)*0.4 *d[b+"err"][inds[b]])
        for b in bands])

    p = [sesar.table2[ind]["Per"], lyrae.find_period(time, flux, ferr)]

    print p

    # m = [lyrae.get_model(p0, time, flux) for p0 in p]
    # c = [lyrae.chi2(m0, time, flux, ferr) for m0 in m]

    # def plot_lc(i, b="g"):
    #     N = 5413
    #     t = np.linspace(min(time), max(time), N)
    #     pl.plot(t%(2*p[i]), m[(i+1)%2](t), "k.", alpha=0.1)

    #     t = np.linspace(0, 2*p[i], N)
    #     pl.plot(t, m[i](t), 'r')

    #     pl.errorbar(time[b]%p[i], flux, yerr=ferr, fmt=".k")
    #     pl.errorbar(time%p[i]+p[i], flux, yerr=ferr, fmt=".k")
    #     pl.annotate(r"$\chi^2 = %.0f$"%c[i], [1,1], xytext=[-3, -3],
    #             xycoords="axes fraction", textcoords="offset points",
    #             horizontalalignment="right", verticalalignment="top")
    #     pl.annotate(r"$T = %.6f$"%p[i], [1,1], xytext=[-3, -3-14-3],
    #             xycoords="axes fraction", textcoords="offset points",
    #             horizontalalignment="right", verticalalignment="top")

    # pl.clf()
    # pl.subplot(211)
    # t = r"$\mathrm{Sesar\,%s}:\,T_0=%.6f,\,T_1=%.6f,"%(_id,p[0],p[1])
    # t += "\,\log_{10}\Delta T=%.2f$"%(np.log10(np.abs(p[0]-p[1])))

    # pl.title(t, fontsize=14.)
    # plot_lc(0)
    # pl.xlim(0, 2*max(p))

    # pl.subplot(212)
    # plot_lc(1)
    # pl.xlim(0, 2*max(p))

    # pl.xlabel(r"$t/T_i$", fontsize=14.)

    # pl.savefig("sesar_test/sesar-%04d.png"%ind)

