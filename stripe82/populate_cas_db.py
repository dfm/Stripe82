#!/usr/bin/env python
# encoding: utf-8
"""
Load the needed tables from CAS into local mongodb

History
-------
2011-06-13 - Created by Dan Foreman-Mackey

"""

import time
import os
import os.path
import threading

import numpy as np

import pyfits
import pymongo

from util import add_fits_table_to_db
import astrometry.util.casjobs as casjobs

casjobs.setup_cookies()
# load username and password from environment
casuser = os.environ['CAS_USER']
caspass = os.environ['CAS_PASS']
# make scratch folder if it doesn't exist
casscratch = os.path.join(os.environ['SDSS_SCRATCH'],'castmp')
if os.path.exists(casscratch) is False:
    os.mkdir(casscratch)

class MyThread(threading.Thread):
    def __init__(self,hdu,col):
        threading.Thread.__init__(self)
        self.hdu = hdu
        self.col = col is True
    def run(self):
        add_fits_table_to_db('cas','stars',self.hdu,opts={'objID':'_id'},
                            meta={'lyrae_candidate': self.col})
        

def get_fields():
    """
    Get the list of fields in Stripe 82 from CAS
    
    History
    -------
    2011-07-01 - Created by Dan Foreman-Mackey
    
    """
    #hdu = pyfits.open('calibration/sdss/large_cas_queries/fieldlist.fit')[1]
    tries = 0
    while tries < 5:
        try:
            cas = casjobs.get_known_servers()['dr7']
            cas.login(casuser,caspass)
            cas.drop_table('fields')
            query = get_fields_query()

            jobid = cas.submit_query(query)

            # wait until the job finishes
            while True:
                jobstatus = cas.get_job_status(jobid)
                if jobstatus is 'Finished':
                    break
                elif jobstatus in ['Failed','Cancelled']:
                    raise Exception(jobstatus)
                time.sleep(10)

            outputfn = 'tmp_fields.fits'
            cas.output_and_download('fields', outputfn, True)
            hdu = pyfits.open(outputfn)[1]
            tries = 100
        except:
            hdu = None
            print "casutils failed!"
            tries += 1

    if hdu is not None:
        add_fits_table_to_db('cas','fields',hdu,clobber=True)
        db = pymongo.Connection().cas
        db.eval("""function() {
            db.fields.find({ramax: {'$gt':180.}}).forEach(function(obj) {
                obj.ramax -= 360.0;
                if (obj.ramin > 180.0) {
                    obj.ramin -= 360.0;
                } else {
                    ramax = obj.ramin;
                    obj.ramin = obj.ramax;
                    obj.ramax = ramax;
                }
                db.fields.save(obj);
            })}""")
        db.fields.ensure_index('ramin')
        db.fields.ensure_index('ramax')
        db.fields.ensure_index('decmin')
        db.fields.ensure_index('decmax')


def get_calibstars():
    """
    Get a list of calibstars from CAS
    
    Parameters
    ----------
    ras : list
        [ramin,ramax] in degrees

    decs : list
        [decmin,decmax] in degrees
    
    History
    -------
    2011-07-01 - Created by Dan Foreman-Mackey
    
    """
    delta_ra = 3 #60
    delta_dec = 3
    ras = range(303.0,360.0,delta_ra) #range(55,60,delta_ra)+range(300,360,delta_ra)
    decs = [-1.5]
    print ras,decs
    for col in [False]:#[True, False]:
        for ra in ras:
            for dec in decs:
                tries = 0
                while tries < 5:
                    try:
                        cas = casjobs.get_known_servers()['dr7']
                        cas.login(casuser,caspass)
    
                        ramin = ra
                        ramax = ra+delta_ra
                        decmin = dec
                        decmax = dec+delta_dec
                        cas.drop_table('stars')
                        if col is False:
                            query = get_calibstar_query(ramin,ramax,decmin,decmax)
                        else:
                            query = get_lyrae_query(ramin,ramax,decmin,decmax)

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
                        f = pyfits.open(outputfn)
                        hdu = f[1]
                        f.close()
                        tries = 100
                    except:
                        hdu = None
                        print "casutils failed!"
                        tries += 1
    
                if hdu is not None:
                    new_thread = MyThread(hdu,col)
                    new_thread.start()
                    #add_fits_table_to_db('cas','stars',hdu,opts={'objID':'_id'},
                    #        meta={'lyrae_candidate': col})

    db = pymongo.Connection().cas
    db.eval("""function() {
            db.stars.find().forEach(function(obj) {
            if (obj.ra > 180) {
                obj.ra -= 360.;
            }
            obj.pos = [obj.ra,obj.dec];
            obj.rank = Math.random();
            db.stars.save(obj);
        })}""")
    db.stars.create_index([('pos',pymongo.GEO2D)])

def get_star_query(ramin,ramax,decmin,decmax,grange):
    return '''SELECT
    p.objID,p.ra,p.dec,p.g,p.Err_g,p.u,p.r,p.i,p.z
INTO mydb.stars
FROM Stripe82..PhotoPrimary p
WHERE p.ra BETWEEN %f AND %f
AND p.dec BETWEEN %f AND %f
AND p.type = 6 AND p.g BETWEEN %f AND %f
'''%(ramin,ramax,decmin,decmax,grange[0],grange[1])

def get_calibstar_query(ramin,ramax,decmin,decmax):
    grange = (15,20)
    return get_star_query(ramin,ramax,decmin,decmax,grange)+'''AND (p.flags &
    (dbo.fPhotoFlags('BRIGHT')+dbo.fPhotoFlags('EDGE')+dbo.fPhotoFlags('BLENDED')
    +dbo.fPhotoFlags('SATURATED')+dbo.fPhotoFlags('NOTCHECKED')
    +dbo.fPhotoFlags('NODEBLEND')+dbo.fPhotoFlags('INTERP_CENTER')
    +dbo.fPhotoFlags('DEBLEND_NOPEAK')+dbo.fPhotoFlags('PEAKCENTER'))) = 0

'''

def get_lyrae_query(ramin,ramax,decmin,decmax):
    grange = (18,22)
    return get_star_query(ramin,ramax,decmin,decmax,grange)+\
            '''AND p.u - p.g BETWEEN 0.7 AND 1.35
AND p.g - p.r BETWEEN -0.15 AND 0.40
AND p.r - p.i BETWEEN -0.15 AND 0.22
AND p.i - p.z BETWEEN -0.21 AND 0.25

'''

def get_fields_query():
    return '''SELECT
    p.run,p.camcol,p.field,p.rerun,p.mjd_g,p.ramin,p.ramax,p.decmin,p.decmax
INTO mydb.fields
FROM Stripe82..Field AS p
WHERE p.mjd_g > 0

'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Scrape data from CAS')
    parser.add_argument('-s','--stars',
                        help='Stars',
                        action='store_true')
    parser.add_argument('-f','--fields',
                        help='Fields',
                        action='store_true')
    parser.add_argument('-q','--queries',
                        help='List queries',
                        action='store_true')
    args = parser.parse_args()
    
    if args.stars:
        get_calibstars()
    if args.fields:
        get_fields()
    if args.queries:
        f = open('fields.sql','w')
        f.write(get_fields_query())
        f.close()
        f = open('stars.sql','w')
        f.write(get_calibstar_query(0,10,-2,2))
        f.close()
        f = open('lyrae.sql','w')
        f.write(get_lyrae_query(0,10,-2,2))
        f.close()

