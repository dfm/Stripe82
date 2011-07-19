#!/usr/bin/env python
# encoding: utf-8
"""
Databases for calibration model

History
-------
2011-06-14 - Created by Dan Foreman-Mackey

"""

__all__ = ['photoraw','photomodel','obslist']

import pymongo
from customBSON import PhotoSONManipulator

# raw photometry database
_db  = pymongo.Connection().photometry
_db.add_son_manipulator(PhotoSONManipulator())
photoraw = _db.raw
photoforced = _db.forced
photomodel = _db.model
obslist = _db.obslist
calibcache = _db.calib

