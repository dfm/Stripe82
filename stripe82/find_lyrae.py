#!/usr/bin/env python
# encoding: utf-8
"""


History
-------
2011-08-14 - Created by Dan Foreman-Mackey

"""

import numpy as np
import pymongo


def main():
    starsdb = pymongo.Connection().cas.stars
    candidates = [obj for obj in starsdb.find(
        {"pos": {"$within": {"$box": [[-29.,-0.4],[-20.,-0.2]]}},
            "lyrae_candidate": True})]
    print len(candidates)

if __name__ == '__main__':
    main()


