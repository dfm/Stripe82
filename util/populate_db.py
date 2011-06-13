#!/usr/bin/env python
# encoding: utf-8
"""
Load a FITS table and store it in a local mongodb

History
-------
2011-06-13 - Created by Dan Foreman-Mackey

"""

__all__ = ['populate_db']

import pymongo

def populate_db(database,collection,hdu,host='localhost',port=27017):
    """
    Use a list of fields from CAS to populate a local mongodb
    
    Parameters
    ----------
    database : string
        Name of database to use

    collection : string
        Name of collection to use (will be dropped first)

    hdu : pyfits.HDU
        The pyfits HDU object containing the data

    Optional
    --------
    host : string (default : 'localhost')
        Host to connect to for database

    port : int (default : 27017)
        If you don't want to use the default mongodb port

    Examples
    --------
    ```python
    >>> hdu = pyfits.open('path/to/table.fits')[1]
    >>> populate_db('cas','fields',hdu)
    ```

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey
    
    """
    # pymongo connection
    con = pymongo.Connection(host=host,port=port)
    db  = con[database]
    db[collection].drop()
    coll = db[collection]
    
    # read in column names
    colname = []
    colconv = [] 
    for c in hdu.columns:
        colname.append(c.name)
        # conversions to types that mongodb understands
        if c.format[-1] in ['I','J','K','B','L']:
            colconv.append(int)
        elif c.format[-1] in ['E','D']:
            colconv.append(float)
        else:
            print "Warning: treating FITS data type %s as str"%(c.format)
            colconv.append(str)

    # loop over rows in FITS table
    for ri in xrange(len(hdu.data)):
        row = hdu.data[ri]
        if ri%(len(hdu.data)/10) is 0:
            print "%.2f percent - %d"%(100*float(ri)/len(hdu.data),ri)
        entry = {}
        for i in xrange(len(colconv)):
            entry[colname[i]] = colconv[i](row[i])

        coll.insert(entry)

if __name__ == '__main__':
    import pyfits
    hdu = pyfits.open('calibration/sdss/large_cas_queries/fieldlist.fit')[1]
    populate_db('cas','fields',hdu)

