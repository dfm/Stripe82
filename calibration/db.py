#!/usr/bin/env python
# encoding: utf-8
"""
Databases for calibration model

History
-------
2011-06-14 - Created by Dan Foreman-Mackey

"""

__all__ = ['photoraw','photomodel']

import pymongo

# raw photometry database
_db  = pymongo.Connection().photometry
photoraw = _db.raw
photomodel = _db.model

