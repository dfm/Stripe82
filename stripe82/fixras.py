#!/usr/bin/env python
# encoding: utf-8
"""
Fix the R.A.s

History
-------
2011-06-14 - Created by Dan Foreman-Mackey

"""

import pymongo

db = pymongo.Connection().cas
# for fields
if True:
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

#for stars
if True:
    db.eval("""function() {
            db.stars.find({ra: {'$gt':180.}}).forEach(function(obj) {
                obj.ra -= 360.0;
                db.stars.save(obj);
            })}""")

