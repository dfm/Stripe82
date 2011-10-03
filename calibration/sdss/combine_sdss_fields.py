# encoding: utf-8
"""
A function for combining adjacent SDSS fields into one image.

History
-------
2011-08-10 - Created by Dan Foreman-Mackey

"""

__all__ = ['combine_sdss_fields']

import numpy as np

def combine_sdss_fields(images,invar_maps,overlap=128):
    """
    Combine adjacent SDSS fields into one image.

    SDSS fields are (for our purposes) images with dimension (1489,2048) and they
    overlap by 128 pixels along the zeroth dimension. That (128,2048) region should
    contain almost _exactly_ the same data. The fields should be arranged with
    increasing field number.

    (1) We take the inverse variance map and ramp one down to zero and the other one
        up from zero in the overlapping range.

    (2) Then, we take the mean of the images in the overlapping region weighted by
        the inverse variance.

    (3) Finally, the inverse variance map in the overlapping region is the _average_
        of the two maps.
    
    Parameters
    ----------
    images : list
        A list of images to combine. Note: we assume that this list is sorted
        by increasing field number.

    invar_maps : list
        The inverse variance maps for each image in images. Should be the same
        length and in the same order as images

    overlap : int
        The number of overlapping pixels in the zeroth dimension.

    Returns
    -------
    img : numpy.ndarray
        The combined image.

    inv : numpy.ndarray
        The combined inverse variance map.

    History
    -------
    2011-08-10 - Created by Dan Foreman-Mackey

    """
    assert(len(images) == len(invar_maps) and len(images) > 1)
    shape = np.shape(images[0])
    assert(len(shape) == 2)
    nfields = len(images)
    combinedShape = [nfields*shape[0]-overlap*(nfields-1),shape[1]]

    img = np.zeros(combinedShape)
    inv = np.zeros(combinedShape)

    img[:shape[0],:] = images[0]
    inv[:shape[0],:] = invar_maps[0]

    for i in np.arange(nfields-1)+1:
        dw = (shape[0]-overlap)
        no_overlap_slice = np.arange(i*dw+overlap,i*dw+shape[0])
        overlap_slice = i*dw + np.arange(overlap)

        # weighted average
        inv0 = inv[overlap_slice,:]*np.linspace(1.,0.,overlap)[:,np.newaxis]
        inv1 = invar_maps[i][:overlap,:]*np.linspace(0.,1.,overlap)[:,np.newaxis]
        totinv = inv0+inv1

        # this line avoids dividing by zero... these pixels don't matter anyways
        # though because the inverse variance is zero
        inds = totinv == 0.0
        totinv[inds] = np.ones(np.sum(inds))
        
        img[overlap_slice,:] = (img[overlap_slice,:]*inv0
                + images[i][:overlap,:]*inv1)/totinv

        # average inverse variance in overlap region
        inv[overlap_slice,:] = (inv[overlap_slice,:]+invar_maps[i][:overlap,:])/2.0

        # new field in rest of frame
        img[no_overlap_slice,:] = images[i][overlap:,:]
        inv[no_overlap_slice,:] = invar_maps[i][overlap:,:]

    return img,inv

if __name__ == '__main__':
    import h5py
    import pylab as pl
    imgs = []
    invs = []
    f = h5py.File('combine_test.hdf5')
    for k in list(f):
        print k[4:]
        if k[:3] == 'img':
            imgs.append(f[k])
        else:
            invs.append(f[k])
    img,inv = combine_sdss_fields(imgs,invs)
    pl.imshow(img,vmax=2000,vmin=1000)
    pl.colorbar()
    pl.show()

