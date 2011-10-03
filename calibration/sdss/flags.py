#!/usr/bin/env python
# encoding: utf-8
"""
Get the flags for a given bitmask

History
-------
2011-07-19 - Created by Dan Foreman-Mackey

"""

__all__ = ['get_flags']

import os
path = os.path.abspath(__file__)
dir_path = os.path.dirname(path)
flags = [line.split() for line in open(os.path.join(dir_path,'flags.dat'))]

def get_flags(mask):
    result = []
    for flag in flags:
        if int(flag[1],16) & mask:
            result.append(flag[0])
    return result

