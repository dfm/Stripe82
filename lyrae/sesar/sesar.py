#!/usr/bin/env python
# encoding: utf-8
"""
Load the catalog of RR Lyrae from Sesar et al. (2010)

History
-------
2011-06-23 - Created by Dan Foreman-Mackey

"""

__all__ = ['coords','table2','table1']

import sys
import os
import os.path

import h5py
import numpy as np

path = os.path.dirname(__file__)
hdf5_fn = os.path.join(path,'sesar.hdf5')

# load tables if not already cached
if not os.path.exists(hdf5_fn):
    hdf5_f = h5py.File(hdf5_fn,'w')
    f = open(os.path.join(path,'Sesar2010_table2.dat'))

    # ignore header
    while 1:
        line = f.readline()
        if not line[0] == '#':
            break
    
    # parse dtype
    dtype = []
    while 1:
        fmt = line[9:15]
        lbl = line[22:26]
        dtype.append((lbl.strip(),fmt))
        line = f.readline()
        if line[0] == '#':
            break
    dtype = np.dtype(dtype)

    # ignore notes
    while 1:
        line = f.readline()
        if not line[0] == '#':
            break

    # read the data
    table2 = []
    while 1:
        table2.append(tuple(line.split()))
        line = f.readline()
        if line == '':
            break
    table2 = np.array(table2,dtype=dtype)
    hdf5_f['table2'] = table2
    f.close()

    dtype = [('ra', np.float64),('dec',np.float64),
             ('umjd',np.float64),('u',np.float64),('uerr',np.float64),
             ('gmjd',np.float64),('g',np.float64),('gerr',np.float64),
             ('rmjd',np.float64),('r',np.float64),('rerr',np.float64),
             ('imjd',np.float64),('i',np.float64),('ierr',np.float64),
             ('zmjd',np.float64),('z',np.float64),('zerr',np.float64)]
    hdf5_f.create_group('lyrae')
    coords = []
    for i,num in enumerate(table2['Num']):
        f = open(os.path.join(path,'Sesar2010_table1','%d.dat'%num))
        data = np.array([tuple(line.split()) for line in f],dtype=dtype)
        f.close()
        hdf5_f['lyrae'][str(num)] = data
        coords.append((num,table2['Per'][i],data['ra'][0],data['dec'][0]))
    coords = np.array(coords,
            dtype=[('Num',np.int),('Per',np.float64),('ra',np.float64),('dec',np.float64)])
    hdf5_f['coords'] = coords

    hdf5_f.close()

hdf5_f = h5py.File(hdf5_fn)
table2 = hdf5_f['table2'][...]
coords = hdf5_f['coords'][...]
table1 = hdf5_f['lyrae']

