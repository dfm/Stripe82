#!/usr/bin/env python
# encoding: utf-8
"""


History
-------
2011-08-14 - Created by Dan Foreman-Mackey

"""

import sys
import os
import time
import cPickle as pickle

import numpy as np
import pymongo

from calibration import calibrate
from calib_results import CalibrationPatch

def main():
    bp = sys.argv[1]
    modbp = os.path.join(bp,"models")
    resbp = os.path.join(bp,"results")
    if not os.path.exists(resbp):
        os.makedirs(resbp)
    if not os.path.exists(modbp):
        os.makedirs(modbp)

    starsdb = pymongo.Connection().cas.stars
    candidates = [obj for obj in starsdb.find(
        {"pos": {"$within": {"$box": [[-29.,-0.4],[-20.,-0.2]]}},
            "lyrae_candidate": True})]
    ind = np.random.randint(len(candidates))
    for ind in np.arange(len(candidates)):
        print "candidate:",ind,"of",len(candidates)
        radius = 5.
        ra,dec = candidates[ind]['ra'],candidates[ind]['dec']
        strt = time.time()

        # extracting and calibrating
        try:
            model = calibrate(ra,dec,radius=radius)
            patch = CalibrationPatch(model,model.data,ra,dec,radius)
            period = patch.get_target()[1].get_period()
            pickle.dump((model.data,model.vector(),ra,dec,radius,period),
                    open(os.path.join(modbp,"%03d.pkl"%(ind)),"wb"),-1)
        except Exception as e:
            print e
            print "Failure"

        #    #fitting lightcurve
        #    data = np.zeros([len(mjd),2])
        #    data[:,0] = flux
        #    data[:,1] = err
        #    period, A = find_period(mjd,data,full_output=True)
        #    print period
        #    pickle.dump((period,A),open(os.path.join(resbp,"%03d.pkl"%(ind)),"wb"),-1)

        print "it took: ",time.time()-strt


if __name__ == '__main__':
    main()


