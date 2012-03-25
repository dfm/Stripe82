#!/usr/bin/env python

import sys
import itertools

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator

f = open("sesar_all.out")
f.readline()

order = 12

data = []
for line in f:
    r = [float(n) for n in line.split()]
    if r[3] == order:
        data.append(r)
data = np.array(data)
A = data[:, -(2*order+1):]

mask = (A[:,0] > 0) * (A[:,0] < 500) \
        * np.all(np.abs(A[:,1:]/A[:,0][:,None]) < 1, axis=1)

ind = np.argmin(np.sqrt(A[mask,1]**2+A[mask,2]**2)/A[mask,0])
print ind
print (np.arange(len(A))[mask, :])[ind]
print (np.sqrt(A[mask,1]**2+A[mask,2]**2)/A[mask,0])[ind]

pl.figure()

# First calculate the covariance matrix and plot the diagonal.
C = np.cov((A[mask,1:]).T) # np.cov((A[mask,1:]/A[mask,0][:,None]).T)
d = C.diagonal()

# The order of each entry in C.diagonal()
n = np.array([[i]*2 for i in range(1, len(C)/2+1)]).flatten()

pl.plot(n, np.log10(C.diagonal()), ".k")
pl.xlabel(r"$n$")
pl.ylabel(r"$\log_{10}C_{ii}$")
pl.savefig("scatter/diagonal_linear.png")

# Also plot it in log-log.
pl.clf()
pl.plot(np.log10(n), np.log10(C.diagonal()), ".k")

m,b = np.polyfit(np.log10(n)[:8], np.log10(C.diagonal())[:8], 1)
print m, b
pl.plot(np.log10(n), b+m*np.log10(n))
pl.plot(np.log10(n), b-2*np.log10(n))

pl.xlabel(r"$\log_{10}n$")
pl.ylabel(r"$\log_{10}C_{ii}$")
pl.title(r"$C_{ii} \propto n^{%.3f}$"%m)
pl.savefig("scatter/diagonal_log.png")

# Then plot the covariance matrix itself.
pl.clf()
pl.pcolor(C)
pl.colorbar()
pl.ylim(pl.gca().get_ylim()[::-1])
pl.savefig("scatter/cov.png")

# pl.clf()
# m = np.mean(A[mask,1::2]**2+A[mask,2::2]**2, axis=0)
# pl.plot(np.log10(np.arange(len(C)/2)+1), np.log10(m), ".k")
# pl.savefig("scatter_mean.png")

sys.exit(0)

nvars = order+1 #A.shape[-1]

fig = pl.figure(figsize=(30, 30))

for i, (yi,xi) in enumerate(itertools.product(np.arange(nvars), np.arange(nvars))):
    if xi <= yi:
        ax = fig.add_subplot(nvars, nvars, i+1)

        if xi == yi:
            ax.hist(A[mask, xi]/A[mask, 0], 50, histtype="step", color="k")
            ax.set_yticklabels([])
        else:
            ax.plot(A[mask, xi]/A[mask, 0], A[mask, yi]/A[mask, 0], ".k")

        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_locator(MaxNLocator(3))

        if xi > 0:
            ax.set_yticklabels([])
        else:
            if yi > 0:
                if yi % 2:
                    ax.set_ylabel(r"$A_{%d}$"%(np.floor(yi/2)+1))
                else:
                    ax.set_ylabel(r"$B_{%d}$"%(yi/2))

        if yi < nvars-1:
            ax.set_xticklabels([])
        else:
            if xi > 0:
                if xi % 2:
                    ax.set_xlabel(r"$A_{%d}$"%(np.floor(xi/2)+1))
                else:
                    ax.set_xlabel(r"$B_{%d}$"%(xi/2))
            for t in ax.get_xticklabels():
                t.set_rotation(60)

pl.savefig("scatter/scatter.png")

