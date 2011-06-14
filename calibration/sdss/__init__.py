#!/usr/bin/env python
# encoding: utf-8
"""
An interface used to access SDSS imaging for use in calibration project

Classes
-------
SDSSField - A wrapper class for an SDSS imaging field

Exceptions
----------
SDSSDASFileError - Raised when a file is not available on DAS

History
-------
Jun 12, 2011 - Created by Dan Foreman-Mackey

"""

from sdss import *
from cas import *


