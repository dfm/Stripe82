#!/usr/bin/env python
# encoding: utf-8
"""
Options file for SDSS imaging interface

History
-------
Jun 12, 2011 - Created by Dan Foreman-Mackey

"""

__all__ = ['das_base','scratch_base']

import os

# URL for the imaging root for DR7 DAS server
das_base = "http://das.sdss.org/imaging"

# local scratch directory
scratch_base = None
if 'SDSS_SCRATCH' in os.environ.keys():
    scratch_base = os.environ['SDSS_SCRATCH']

