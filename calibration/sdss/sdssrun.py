#!/usr/bin/env python
# encoding: utf-8
"""
Interface for accessing SDSS imaging data

Notes
-----
Requires utilities from the astrometry.net codebase (trac.astrometry.net) in
Python path.

"""

__all__ = ['SDSSDASFileError','SDSSRunMissingFieldsError','SDSSRun']

import os
import os.path
import subprocess
import hashlib
import cPickle as pickle

import numpy as np
import pyfits
import h5py

import astrometry.util.sdss_psf as sdss_psf
from astrometry.sdss.common import *
from astrometry.sdss import DR7

from config import *

class SDSSDASFileError(Exception):
    """
    Custom exception raised when an image cannot be loaded from DAS

    Parameters
    ----------
    url : string
        URL that was requested

    History
    -------
    Created by Dan Foreman-Mackey on Jun 12, 2011

    """
    def __init__(self,url):
        self.url = url

    def __str__(self):
        return str(self.url)

class SDSSRunMissingFieldsError(Exception):
    pass

class SDSSOutOfBounds(Exception):
    """
    Raised when a star is not in the requested image

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey

    """
    pass

class DFMDR7(DR7):
    """
    A Hack of Lang and Hogg's DR7 class

    Less verbose. Automatically retrieves image from DAS if it's not
    available locally

    Examples
    --------
    ```python
    >>> sdss = DFMDR7()
    >>> # get the g-band image of (run,camcol,field)=(1056,3,84)
    >>> fn = sdss.getFilename('fpC', 1056, 3, 84, 'g')
    ```

    TODO
    ----
    - Hack all the file opening functions so that they are less verbose

    History
    -------
    2011-06-12 - Created by Dan Foreman-Mackey

    """
    def __init__(self):
        super(DFMDR7, self).__init__()
        self.filenames = {
                'fpObjc': 'fpObjc-%(run)06i-%(camcol)i-%(field)04i.fit',
                'fpM': 'fpM-%(run)06i-%(band)s%(camcol)i-%(field)04i.fit',
                'fpC': 'fpC-%(run)06i-%(band)s%(camcol)i-%(field)04i.fit',
                'fpAtlas': 'fpAtlas-%(run)06i-%(camcol)i-%(field)04i.fit',
                'psField': 'psField-%(run)06i-%(camcol)i-%(field)04i.fit',
                'tsObj': 'tsObj-%(run)06i-%(camcol)i-%(rerun)i-%(field)04i.fit',
                'tsField': 'tsField-%(run)06i-%(camcol)i-%(rerun)i-%(field)04i.fit',
            }

    def getFilename(self, filetype, *args, **kwargs):
        """
        Get the relative path to an image from SDSS

        Parameters
        ----------
        filetype : string
            'fpC', 'psField', 'fpAtlas', 'fpObjc', 'fpM', 'tsObj' or 'tsField'

        run : int
            SDSS run number

        camcol : int
            Camera column

        field : int
            Field number

        Optional
        --------
        rerun : int (default: 40)
            Processing rerun id

        Returns
        -------
        path : string
            Path to image

        History
        -------
        2011-06-12 - Created by Dan Foreman-Mackey

        """
        x = zip(['run', 'camcol', 'field', 'band'], args)
        for k,v in x:
            kwargs[k] = v
        if not filetype in self.filenames:
            return None
        pat = self.filenames[filetype]
        fn = pat % kwargs

        # hack!
        # maybe we should always request the rerun number?
        if 'rerun' not in kwargs.keys():
            kwargs['rerun'] = 40

        prefix = '%(rerun)i/%(run)i/' % kwargs

        if filetype in ['fpC']:
            return prefix + 'corr/%(camcol)i/' % kwargs + fn
        elif filetype in ['psField', 'fpAtlas', 'fpObjc', 'fpM']:
            return prefix + 'objcs/%(camcol)i/' % kwargs + fn
        elif filetype in ['tsObj', 'tsField']:
            return prefix + 'calibChunks/%(camcol)i/' % kwargs + fn
        else:
            return None

    def readTsField(self, run, camcol, field, rerun):
        """
        Get the astrometry data

        Parameters
        ----------
        run : int
            SDSS run number

        camcol : int
            Camera column

        field : int
            Field number

        Returns
        -------
        tsField : TsField
            The TsField object

        References
        ----------
        http://www.sdss.org/dr7/dm/flatFiles/tsField.html

        History
        -------
        2011-06-14 - Created by Dan Foreman-Mackey

        """
        f = TsField(run, camcol, field, rerun=rerun)
        fn = self.getFilename('tsField', run, camcol, field, rerun=rerun)
        p = self._open(fn)
        f.setHdus(p)
        return f

    def readFpM(self, run, camcol, field, band):
        """
        Get the mask image

        Parameters
        ----------
        run : int
            SDSS run number

        camcol : int
            Camera column

        field : int
            Field number

        Returns
        -------
        fpM : FpM
            The mask object

        References
        ----------
        http://data.sdss3.org/datamodel/files/PHOTO_REDUX/RERUN/RUN/objcs/CAMCOL/fpM.html

        History
        -------
        2011-06-14 - Created by Dan Foreman-Mackey

        """
        f = FpM(run, camcol, field, band)
        fn = self.getFilename('fpM', run, camcol, field, band)
        p = self._open(fn+'.gz')
        f.setHdus(p)
        return f

    def readPsField(self, run, camcol, field):
        """
        Get the data for the psField image

        Parameters
        ----------
        run : int
            SDSS run number

        camcol : int
            Camera column

        field : int
            Field number

        Returns
        -------
        f : PsField
            The resulting PsField object

        p : ??
            The result of pyfits.open

        References
        ----------
        http://data.sdss3.org/datamodel/files/PHOTO_REDUX/RERUN/RUN/objcs/CAMCOL/psField.html

        History
        -------
        2011-06-13 - Created by Dan Foreman-Mackey

        """
        f = PsField(run, camcol, field)
        # ...
        fn = self.getFilename('psField', run, camcol, field)
        p = self._open(fn)
        f.setHdus(p)
        return (f,p)

    def _open(self, fn):
        """
        Load the data for a given filename locally or from DAS server

        Parameters
        ----------
        fn : string
            The relative path to the image (the result of self.getFilename)

        Returns
        -------
        hdus : ???
            The results of a pyfits.open operation

        Raises
        ------
        SDSSDASFileError :
            If the file is not available locally or on DAS

        History
        -------
        2011-06-13 - Created by Dan Foreman-Mackey

        """
        if scratch_base is not None:
            local_path = os.path.join(scratch_base, fn)
            try: # do we have it locally? This should also fail if it's corrupted???
                # FIXME: It actually DOESN'T FAIL IF THE FILE IS CORRUPTED
                return pyfits.open(local_path)
            except:
                file_dir = os.path.join(os.path.split(local_path)[:-1])[0]
                if os.path.exists(file_dir) is False:
                    os.makedirs(file_dir)

                nyu_path = os.path.join(nyu_base, fn)

                # scp nyu_path local_path
                print "scp %s %s"%(nyu_path,local_path)
                ret = subprocess.Popen("scp %s %s"%(nyu_path,local_path),
                        shell=True,stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE).wait()

                if ret is not 0:
                    os.remove(local_path)
                    raise SDSSDASFileError(das_path)

                return pyfits.open(local_path)
        else:
            path = os.path.join(das_base, fn)
            return pyfits.open(path)

# hdf5 tags
metaTag = 'meta'
imgTag = 'img'
invTag = 'inv'
astTag = 'ast'
psfTag = 'psf'
approxAstTag = 'approxAst'
fieldsTag = 'fields'
mjdTag = 'mjd'

class SDSSRun:
    """
    Wrapper class for SDSS image file with various convenience functions

    Parameters
    ----------
    run : int
        SDSS run number

    camcol : int
        Camera column

    Optional
    --------
    filename : str, optional
        Path to a pre-calculated saved HDF5 file to load from

    fields : list, optional
        List of Field objects

    band : str, optional
        SDSS band pass (default: g)

    Raises
    ------
    SDSSDASFileError :
        If any of the requested files are not available locally or on DAS

    TODO
    ----
    - Decide what to do about multiple bands

    """
    def __init__(self, arg, band='g'):
        if isinstance(arg, str) or isinstance(arg, unicode):
            filename = arg
            fields   = None
        elif isinstance(arg, list):
            filename = None
            fields   = arg
        else:
            raise TypeError()

        if filename is not None:
            assert os.path.exists(filename), "The file %s doesn't exist"%filename
            self.filename = filename

            self.data = h5py.File(self.filename)

            self._run, self._camcol, self._band = \
                    self.data[metaTag]

            self._ast  = {}
            for k in list(self.data[astTag]):
                self._ast[k] = pickle.loads(str(self.data[astTag][k][...]))

            self._fields    = self.data[fieldsTag][...]
            self._minField  = self._fields[0]
            self._maxField  = self._fields[-1]
            self._approxAst = self.data[approxAstTag][...]
            self._mjds      = self.data[mjdTag][...]
            # print "Finished loading"

        else:
            f = fields[0]
            self._run    = f['run']
            self._camcol = f['camcol']
            self._band   = band
            band_id = 'ugriz'.index(self._band)

            assert(scratch_base is not None)
            self.filename = os.path.join(scratch_base,
                    '%d-%d-%s.hdf5'%(self._run,self._camcol,self._band))

            self._psfData,self._psField = {},{}
            self._fields = []

            self.data   = h5py.File(self.filename,'w')
            self.data.create_group(astTag)

            self.data[metaTag] = [self._run, self._camcol, self._band]

            self._fieldlist   = []
            self._fieldcoords = [] # mean RA/Dec points in fields
            self._mjds        = []
            self._ast         = {}

            # copy documents first (avoid timeout!)
            # q = {'run': run,'camcol': camcol}
            # for k in list(field_selection_query):
            #     q[k] = field_selection_query[k]
            # docs = [doc for doc in obsdb.find(q,{'field':1,'rerun':1,'mjd_%s'%(band): 1}).sort('field')]
            # print 'done query'

            self._minField = fields[0]['field']
            self._maxField = fields[-1]['field']
            deltaFields = self._maxField - self._minField
            if len(fields) < deltaFields-1:
                raise SDSSRunMissingFieldsError("Warning: missing fields")

            # MAGIC: shape of SDSS fields
            fullShape = [(deltaFields+1)*field_width-field_overlap*deltaFields,
                            field_height]
            self.data.create_dataset(imgTag,fullShape,float,compression='gzip')
            self.data.create_dataset(invTag,fullShape,float,compression='gzip')

            print "SDSSRun is loading all ",len(fields)," fields for run: ",\
                    self._run," and camcol: ",self._camcol
            for f in fields:
                assert(f['run'] == self._run and f['camcol'] == self._camcol)
                field = f['field']
                self._fields.append(field)
                rerun = f['rerun']
                self._mjds.append(f['mjd_%s'%band])
                sdss    = DFMDR7()
                tsField = sdss.readTsField(self._run, self._camcol, field, rerun)
                self._ast[str(field)] = tsField.getAsTrans(self._band)

                (psField,psfData) = sdss.readPsField(self._run, self._camcol, field)
                self._psfData["%d"%field] = psfData
                self._psField["%d"%field] = psField

                fpCPath = sdss.getFilename('fpC',
                        self._run, self._camcol, field, self._band)
                fpCPath += '.gz' # always zipped?
                fpCfits = sdss._open(fpCPath)
                fpC = FpC(self._run, self._camcol, field, self._band)
                fpC.image = fpCfits[0].data
                img = fpC.getImage()

                # mask file
                fpM = sdss.readFpM(self._run, self._camcol, field, self._band)

                # inverse variance
                inv = sdss.getInvvar(fpC.getImage(), fpM,
                                        psField.getGain(band_id),
                                        psField.getDarkVariance(band_id),
                                        psField.getSky(band_id),
                                        psField.getSkyErr(band_id))

                # position of this field in the full image
                field_id = field-self._minField
                dw = (field_width-field_overlap)
                pos = np.arange(field_id*dw+field_overlap,field_id*dw+field_width)
                overlap = field_id*dw + np.arange(field_overlap)

                if np.any(self.data[invTag][overlap,:] > 0.0):
                    # weighted average
                    inv0 = self.data[invTag][overlap,:]\
                            *np.linspace(1.,0.,field_overlap)[:,np.newaxis]
                    inv1 = inv[:field_overlap,:]\
                            *np.linspace(0.,1.,field_overlap)[:,np.newaxis]
                    totinv = inv0+inv1
                    totinv[totinv == 0.0] = \
                            np.ones(np.shape(totinv[totinv == 0.0]))

                    self.data[imgTag][overlap,:] = \
                            (self.data[imgTag][overlap,:]*inv0
                            + img[:field_overlap,:]*inv1)/totinv

                    # average inverse variance in overlap region
                    self.data[invTag][overlap,:] = (self.data[invTag][overlap,:]\
                            + inv[:field_overlap,:])/2.0

                    # new field in rest of frame
                    self.data[imgTag][pos,:] = img[field_overlap:,:]
                    self.data[invTag][pos,:] = inv[field_overlap:,:]
                else:
                    pos = np.arange(field_id*dw,field_id*dw+field_width)
                    self.data[imgTag][pos,:] = img
                    self.data[invTag][pos,:] = inv

                self.data[astTag][str(field)] = \
                        pickle.dumps(self._ast[str(field)],-1)

                # approximate astrometry
                # get the RA/Dec coordinate at the centre of the frame
                # check order!
                ra,dec = self._ast[str(field)].pixel_to_radec(field_height/2,
                        field_width/2,
                        color=0.0)
                if ra >= 180.:
                    ra -= 360.
                self._fieldcoords.append([ra,dec])

                # close all the open fits files?
                del sdss

            self.data[approxAstTag] = np.array(self._fieldcoords)
            self.data[fieldsTag]    = np.array(self._fields)
            self.data[mjdTag]       = np.array(self._mjds)
            self.data.close()

    def find_closest_field(self, ra, dec):
        """
        NAME:
            find_closest_field
        PURPOSE:
            find the field center closest to ra/dec position
        INPUT:
            ra,dec - (float,float) in degrees
        OUTPUT:
            the field number (int)
        HISTORY:
            Created by Dan Foreman-Mackey on Aug 10, 2011
        """
        ras = self.data[approxAstTag][:,0]
        inds = np.arange(len(ras))
        f0 = inds[ras <= ra][-1]
        if f0 == len(ras)-1:
            return self._fields[f0]

        p0 = self._approxAst[f0,:]
        p1 = self._approxAst[f0+1,:]
        d0 = (p0[0]-ra)**2+(p0[1]-dec)**2
        d1 = (p1[0]-ra)**2+(p1[1]-dec)**2

        if d0 < d1:
            return self._fields[f0]

        return self._fields[f0+1]

    def radec_to_pixel(self, ra, dec, color=None,infield=None):
        """
        NAME:
            radec_to_pixel
        PURPOSE:
            return the pixel position of a given RA/Dec within the field
        INPUT:
            ra,dec - (float,float) in degrees
        OPTIONAL
            color  - (float; default = 0.0)
        OUTPUT:
            (x0,y0) - in field coordinates
            (x,y)   - in combined image coordinates
            infield - which field is it in
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 06, 2011
        """
        if color == None:
            color = 0.0
        if infield is None:
            infield = self.find_closest_field(ra,dec)
        x0,y0 = self._ast["%d"%infield].radec_to_pixel(ra, dec, color=color)

        field_id = infield-self._minField
        dw = (field_width-field_overlap)

        return (x0,y0),(x0,y0+field_id*dw),infield

    def mjd_at_radec(self,ra,dec):
        """
        Get the time at a given RA/Dec interpolated to the pixel location

        Parameters
        ----------
        ra : float
            In degrees

        dec : float
            In degrees

        Returns
        -------
        mjd : float
            Time at that position

        History
        -------
        2011-08-24 - Created by Dan Foreman-Mackey

        """
        dw = (field_width-field_overlap)
        xy0,xy,field_id = self.radec_to_pixel(ra,dec)
        f_id = field_id - self._minField
        y = xy[1]
        if f_id+1 >= np.shape(self.mjds):
            y0 = (len(self._mjds)-2)*dw
            y1 = (len(self._mjds)-1)*dw
            t0 = self._mjds[-2]
            t1 = self._mjds[-1]
        else:
            y0 = f_id*dw
            y1 = (f_id+1)*dw
            t0 = self._mjds[f_id]
            t1 = self._mjds[f_id+1]
        return (t1-t0)/(y1-y0)*(y-y0)+t0

    def psf_at_radec(self, ra, dec):
        """
        NAME:
            psf_at_radec
        PURPOSE:
            return the set psf images centered at (ra,dec)
        INPUT:
            ra,dec - (float,float) in degrees
        OUTPUT:
            a dictionary with one entry for each self.bands: each
            entry is the psf image in that band
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 06, 2011
        """
        xy0,xy,field_id = self.radec_to_pixel(ra,dec)
        return self.psf_at_point(xy0[0],yy0[1],field_id)

    def psf_at_point(self, x, y, field, radius=25, dblgauss=True):
        """
        NAME:
            psf_at_point
        PURPOSE:
            return the psf image centered at (x,y)
        INPUT:
            x,y  - (float,float) pixel positions in _field_ coordinates
            field  - (int) which field is it in?
        OPTIONAL:
            dblgauss - (bool; default=True) use the double Gaussian model
                       otherwise, use the KL basis functions
            radius   - (int, default=25) size of returned image (only implemented
                       for double Gaussian)
        OUTPUT:
            the psf image
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 06, 2011
        """
        field_id = "%s"%(field)
        if field_id not in self._psfData:
            (psField,psfData) = DFMDR7()\
                    .readPsField(self._run, self._camcol, field)
            self._psField[field_id] = psField
            self._psfData[field_id] = psfData

        if dblgauss: # Double Gaussian
            a,s1,b,s2 = self._psField[field_id].getDoubleGaussian(
                    'ugriz'.index(self._band))
            a /= (s1**2 + b*s2**2)
            psf = np.zeros((2*radius+1,2*radius+1))
            xy = np.arange(-radius,radius+1)
            rad = np.sqrt((xy[np.newaxis,:])**2+(xy[:,np.newaxis])**2)
            psf += a*(np.exp(-rad**2/2.0/s1**2) \
                    + b*np.exp(-rad**2/2.0/s2**2))/(2*np.pi)
            return psf

        # KL...
        psf = sdss_psf.sdss_psf_at_points(
                self._psfData[field_id]['ugriz'.index(self._band)], x, y)

        return psf

    def photo_at_radec(self, ra, dec, delta=25):
        """
        NAME:
            photo_at_radec
        PURPOSE:
            fit photometry of a point source at coordinates (RA/Dec)
        INPUT:
            ra,dec - (float) in degrees
        OPTIONAL:
            delta  - (int; default = 25) the output images will have shape (2*delta,2*delta)
        OUTPUT:
            (res,img,inv,psf):
                res - (dictionary) one entry for each self.bands: each tuple is
                      in the output from fit_psf_leastsq
                img - (dictionary) one entry for each self.bands: each entry is
                      a thumbnail image of the raw pixel counts (in arbitrary,
                      uncalibrated flux units)
                inv - (dictionary) one entry for each self.bands: same as img but
                      for the inverse variance map
                psf - (dictionary) one entry for each self.bands: same as img but
                      for the psf image instead of the raw pixels
        RAISES:
            SDSSOutOfBounds - if star is not in field
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 06, 2011
        """
        res = []
        psf_out = []
        img_out = []
        inv_out = []

        xy0,xy,field_id = self.radec_to_pixel(ra,dec)

        psf = self.psf_at_point(xy0[0], xy0[1], field_id)

        x,y = int(xy[0]),int(xy[1])
        xmin,xmax = x-delta,x+delta+1
        ymin,ymax = y-delta,y+delta+1

        fullShape = self.data[imgTag].shape

        if not (0 <= xmin and xmax < fullShape[1] \
                and 0 <= ymin and ymax < fullShape[0]):
            raise SDSSOutOfBounds('Not in bounds')

        #if xmin < 0:
        #    xmin = 0
        #if ymin < 0:
        #    ymin = 0
        #if xmax >= fullShape[1]:
        #    xmax = fullShape[1]-1
        #if ymax >= fullShape[0]:
        #    ymax = fullShape[0]-1
        #if xmin >= xmax or ymin >= ymax:
        #    raise SDSSOutOfBounds('Not in bounds')

        # trim image and inv. variance
        img = self.data[imgTag][ymin:ymax,xmin:xmax]
        inv = self.data[invTag][ymin:ymax,xmin:xmax]

        # trim psf
        #psf = psf[ymin-(y-delta):ymax-(y-delta), \
        #        xmin-(x-delta):xmax-(x-delta)]
        #psf_out.append((b,psf))

        return self.fit_psf_leastsq(img,inv,psf),img,inv,psf

    def fit_psf_leastsq(self,img,inv,psf):
        """
        NAME:
            fit_psf_leastsq
        PURPOSE:
            return stellar flux based on our photometry model using linear
            regression
        INPUT:
            img  -  image data in counts
            inv  -  the inverse variance model of the pixels
            psf  -  the psf image at this point
        OUTPUT:
            return (f0, fstar, dx, dy),cov_matrix
            f0    - the background flux level
            fstar - the total stellar flux
            dx,dy - offsets of psf center (to get in units of
                    pixels, you have to divide by fstar)
        HISTORY:
            Created by Dan Foreman-Mackey on Apr 30, 2011
        """
        # calculate derivatives using finite difference
        dpdx,dpdy = np.zeros(np.shape(psf)),np.zeros(np.shape(psf))
        dpdx[:-1,:] = psf[1:,:]-psf[:-1,:]
        dpdy[:,:-1] = psf[:,1:]-psf[:,:-1]
        X = np.ones([len(psf.flatten()),4])
        X[:,1] = psf.flatten()
        X[:,2] = dpdx.flatten()
        X[:,3] = dpdy.flatten()
        Y = img.flatten().T
        Cinv = np.diag(inv.flatten())
        XC     = np.dot(X.T,Cinv)
        XCY    = np.dot(XC,Y)
        XCXinv = np.linalg.inv(np.dot(XC,X))

        # returns [[b,m],[b_err,m_err]]
        return np.dot(XCXinv,XCY), \
           XCXinv


if __name__ == '__main__':
    import time
    strt = time.time()
    obs = SDSSObservation(run=4858, camcol=6)#run=4253,camcol=3,usecache=True)
    print time.time()-strt
    import sys
    sys.exit()

    radecs = [[30.69153051,-4.779E-4],
        [30.69639824,-0.01486342],
        [30.6891575,0.01390988  ],
        [30.70598525,0.03174024 ],
        [30.65279247,-0.0192601 ],
        [30.73283836,-0.03315539],
        [30.6463008,0.02177736  ],
        [30.67756003,0.0466423  ],
        [30.67582076,0.04486459 ],
        [30.71956806,0.03948578 ],
        [30.68072529,0.03519326 ],
        [30.68523937,0.04608791 ],
        [30.72507914,0.03467313 ],
        [30.73878796,0.01064782 ],
        [30.67761626,-0.06545483],
        [30.64620022,-0.06272369],
        [30.71772923,-0.05313227],
        [30.7329595,-0.04568421 ],
        [30.75866333,-0.01027296],
        [30.63542478,0.01097651 ],
        [30.64070452,0.03897599 ],
        [30.6232155,0.03546206  ],
        [30.73897639,0.0582784  ],
        [30.73377848,0.05628999 ],
        [30.67619731,0.06146714 ],
        [30.59960924,-5.89344E-3],
        [30.65667631,-0.08545852],
        [30.6792161,-0.08202459 ],
        [30.70129416,-0.08481   ],
        [30.79013993,-0.04793415],
        [30.59124629,0.03654471 ],
        [30.63070591,0.06456958 ],
        [30.70059653,0.08581666 ],
        [30.7649538,0.06131464  ],
        [30.67272733,0.09855124 ],
        [30.70295564,0.08398645 ],
        [30.53498664,-0.03416116],
        [30.5342788,-0.03241634 ],
        [30.54378633,-0.01388675],
        [30.563629,-0.09885316  ],
        [30.57872783,-0.02066558],
        [30.57057097,-0.04849638],
        [30.63945783,-0.14167086],
        [30.66095008,-0.16018594],
        [30.68429081,-0.12436866],
        [30.71254481,-0.15475524],
        [30.75191108,-0.11395759],
        [30.7193307,-0.14650902 ],
        [30.79652836,-0.0456693 ],
        [30.82361762,-0.04303993]]

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as pl

    for i,radec in enumerate(radecs):
        try:
            res,img,inv,psf = obs.photo_at_radec(*(tuple(radec)))

            pl.figure()
            pl.subplot(221)
            pl.imshow(img)
            pl.xlabel('img')

            pl.subplot(222)
            pl.imshow(inv)
            pl.xlabel('inv')

            pl.subplot(223)
            pl.imshow(psf)
            pl.xlabel('psf')

            # residuals
            dpdx,dpdy = np.zeros(np.shape(psf)),np.zeros(np.shape(psf))
            dpdx[:-1,:] = psf[1:,:]-psf[:-1,:]
            dpdy[:,:-1] = psf[:,1:]-psf[:,:-1]

            pl.subplot(224)
            pl.imshow(res[0][0]+res[0][1]*psf+res[0][2]*dpdx+res[0][3]*dpdy - img)
            pl.xlabel('diff')

            pl.savefig('blurgle/%d.pdf'%i)
        except:
            print "failed"

