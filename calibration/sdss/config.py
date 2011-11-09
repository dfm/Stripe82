#!/usr/bin/env python
# encoding: utf-8
"""
Configuration file for SDSS imaging interface

"""

__all__ = ['das_base','nyu_base','scratch_base',
    'field_width','field_height','field_overlap','field_selection_query']

import os

field_selection_query = {'ramax': {'$gt': -29.0}, 'ramin': {'$lt': -20.0}}

# URL for local NYU data
nyu_base = "bootes:/mount/coma1/bw55/sdss3/mirror/bosswork/groups/boss/photo/redux/"

# URL for the imaging root for DR7 DAS server
das_base = "http://das.sdss.org/imaging"

# local scratch directory
scratch_base = None
if 'SDSS_SCRATCH' in os.environ.keys():
    scratch_base = os.environ['SDSS_SCRATCH']

field_width   = 1489
field_height  = 2048
field_overlap = 128

