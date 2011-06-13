#!/usr/bin/env python
# encoding: utf-8
"""
Load the needed tables from CAS into local mongodb

History
-------
2011-06-13 - Created by Dan Foreman-Mackey

"""

import pyfits

from util import add_fits_table_to_db

if False: # copy the fields list into database
    hdu = pyfits.open('calibration/sdss/large_cas_queries/fieldlist.fit')[1]
    add_fits_table_to_db('cas','fields',hdu,clobber=True)

if True: # get the info for all the stars
    

