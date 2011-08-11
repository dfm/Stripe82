#!/usr/bin/env python
# encoding: utf-8
"""
Test the likelihood function

History
-------
2011-08-05 - Created by Dan Foreman-Mackey

"""

import pickle
import numpy as np
from calibration import PhotoModel,odds_variable
from calibration._likelihood import lnoddsvar,lnoddsbad

model = PhotoModel(*pickle.load(open('tmp.pkl','rb')))

#print lnlikelihood(model,model.data)+7367.27456661
#oarr = np.zeros(model.data.nstars)
#lnoddsvar(model,model.data,oarr)
#print oarr-odds_variable(model,model.data)

oarr = np.zeros(np.shape(model.data.flux))
lnoddsbad(model,model.data,oarr)
print np.shape(oarr)
print oarr


