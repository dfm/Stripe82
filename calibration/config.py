#!/usr/bin/env python
# encoding: utf-8
"""
Configuration options for the calibration module

"""

__all__ = ['survey','nmgroups']

import numpy as np

#
# Default options
#

import sdss as survey

# the number of measurement groups for cross-validation
nmgroups = 10

