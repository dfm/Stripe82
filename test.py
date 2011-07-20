#!/usr/bin/env python
# encoding: utf-8

import timeit

import numpy as np

import stripe82.calibration as calib
from stripe82.calibration._likelihood import lnlikelihood as like_C
from stripe82.calibration.photomodel import lnlike as like_py

np.random.seed()

obs,stars = calib.find_photometry(21,0,5)
data = calib.get_photometry(obs,stars)
photo_data = calib.PhotoData(data,obs,stars)
p0 = calib.init_model(photo_data)
photo_model = calib.PhotoModel(photo_data,p0)

N = 5000
print "Time in C (us):", 1e6*np.mean(
        timeit.repeat("like_C(photo_model,photo_data)",
        "from __main__ import like_C,photo_model,photo_data",
        number=N))/N
print "Time in Python (us):", 1e6*np.mean(
        timeit.repeat("like_py(photo_model,photo_data)",
        "from __main__ import like_py,photo_model,photo_data",
        number=N))/N

res_py = like_py(photo_model,photo_data)
res_C  = like_C(photo_model,photo_data)
print res_C,res_py,res_C-res_py

