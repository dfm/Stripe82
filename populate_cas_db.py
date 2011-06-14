#!/usr/bin/env python
# encoding: utf-8
"""
Load the needed tables from CAS into local mongodb

History
-------
2011-06-13 - Created by Dan Foreman-Mackey

"""

import time

import pyfits

from util import add_fits_table_to_db

if False: # copy the fields list into database
    hdu = pyfits.open('calibration/sdss/large_cas_queries/fieldlist.fit')[1]
    add_fits_table_to_db('cas','fields',hdu,clobber=True)

if True: # get the info for all the stars
    import os
    import os.path
    import numpy as np
    import astrometry.util.casjobs as casjobs

    casjobs.setup_cookies()
    
    # load username and password from environment
    casuser = os.environ['CAS_USER']
    caspass = os.environ['CAS_PASS']
    
    # make scratch folder if it doesn't exist
    casscratch = os.path.join(os.environ['SDSS_SCRATCH'],'castmp')
    if os.path.exists(casscratch) is False:
        os.mkdir(casscratch)
    
    grange = (0,20)
    delta = 1 # 1x1 degree patches
    ras = np.arange(0,59,delta)
    decs = np.arange(-1.25,1.25,delta)
    for ra in ras:
        for dec in decs:
            tries = 0
            while tries < 5:
                try:
                    cas = casjobs.get_known_servers()['dr7']
                    cas.login(casuser,caspass)

                    ramin = ra%360.
                    ramax = (ra+delta)%360.
                    if ramax == 0:
                        ramax = 360.
                    cas.drop_table('stars')
                    query = '''SELECT
    p.ra,p.dec,p.g,p.Err_g
INTO mydb.stars
FROM PhotoPrimary p
WHERE p.ra BETWEEN %f AND %f
AND p.dec BETWEEN %f AND %f
AND (p.flags &
    (dbo.fPhotoFlags('BRIGHT')+dbo.fPhotoFlags('EDGE')+dbo.fPhotoFlags('BLENDED')
    +dbo.fPhotoFlags('SATURATED')+dbo.fPhotoFlags('NOTCHECKED')
    +dbo.fPhotoFlags('NODEBLEND')+dbo.fPhotoFlags('INTERP_CENTER')
    +dbo.fPhotoFlags('DEBLEND_NOPEAK')+dbo.fPhotoFlags('PEAKCENTER'))) = 0
    AND p.type = 6 AND p.g BETWEEN %f AND %f
'''%(ramin,ramax,dec,dec+delta,grange[0],grange[1])
                    # NOTE: (above) R.A.s in CAS need to be mod 360

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
                    cas.output_and_download('stars', outputfn, True)
                    hdu = pyfits.open(outputfn)[1]
                    tries = 100
                except:
                    hdu = None
                    print "casutils failed!"
                    tries += 1

            if hdu is not None:
                add_fits_table_to_db('cas','stars',hdu)

