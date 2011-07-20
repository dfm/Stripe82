#!/usr/bin/env python
# encoding: utf-8
"""
Convenience functions for accessing the MongoDB instance

History
-------
2011-06-13 - Created by Dan Foreman-Mackey

"""

__all__ = ['stardb','obsdb']

import pymongo

# MongoDB databaase names
casdb = 'cas'
casstars = 'stars'
casfields = 'fields'

# connect to mongodb
_db = pymongo.Connection()[casdb]
stardb  = _db[casstars]
obsdb = _db[casfields]


