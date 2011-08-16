#!/usr/bin/env python
# encoding: utf-8

import os
from distutils.core import setup
from distutils.extension import Extension
import numpy.distutils.misc_util

if 'HOME' in os.environ and os.environ['HOME'] == '/Users/dfm':
    extra_compile_args = []
    extra_link_args = []
else:
    extra_compile_args=['-fopenmp','-DUSEOPENMP']
    extra_link_args=['-lgomp']

ext = Extension('stripe82.calibration._likelihood',
                ['stripe82/calibration/_likelihood.c'],
                extra_compile_args=extra_compile_args,
                extra_link_args=extra_link_args)

setup(packages=['stripe82'],
        ext_modules = [ext],
        include_dirs = numpy.distutils.misc_util.get_numpy_include_dirs()
        )

