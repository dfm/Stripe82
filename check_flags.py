#!/usr/bin/env python
# encoding: utf-8
"""
Find the flags for all of the lyrae

History
-------
2011-07-19 - Created by Dan Foreman-Mackey

"""

import os
import os.path
import time

import pyfits

import astrometry.util.casjobs as casjobs
import stripe82.lyrae.sesar as sesar
from stripe82.calibration.sdss.flags import get_flags

casjobs.setup_cookies()
# load username and password from environment
casuser = os.environ['CAS_USER']
caspass = os.environ['CAS_PASS']
# make scratch folder if it doesn't exist
casscratch = os.path.join(os.environ['SDSS_SCRATCH'],'castmp')
if os.path.exists(casscratch) is False:
    os.mkdir(casscratch)

#resultsfile = open('check_flags.dat','w')
#resultsfile.close()

for c in sesar.coords[450:]:
    tries = 0
    while tries < 5:
        try:
            cas = casjobs.get_known_servers()['dr7']
            cas.login(casuser,caspass)
            cas.drop_table('tmp')
            query = '''SELECT p.flags,p.g
INTO mydb.tmp 
FROM fGetNearbyObjEq(%f,%f,0.1) n, PhotoPrimary p
WHERE n.objID=p.objID
'''%(c[2]%360.,c[3])

            jobid = cas.submit_query(query)

            # wait until the job finishes
            while True:
                jobstatus = cas.get_job_status(jobid)
                if jobstatus is 'Finished':
                    break
                elif jobstatus in ['Failed','Cancelled']:
                    raise Exception(jobstatus)
                time.sleep(10)

            outputfn = 'tmp.fits'
            cas.output_and_download('tmp', outputfn, True)
            hdu = pyfits.open(outputfn)[1]
            tries = 100
        except:
            hdu = None
            print "casutils failed!"
            tries += 1

    if hdu is not None and hdu.data is not None:
        flags = get_flags(int(hdu.data[0][0]))
        resultsfile = open('check_flags.dat','a')
        resultsfile.write("%10d\t%16.8e"%(c[0],hdu.data[0][1]))
        for flag in flags:
            resultsfile.write("\t%s"%(flag))
        resultsfile.write("\n")
        resultsfile.close()

