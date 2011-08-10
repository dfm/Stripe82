#!/usr/bin/env python
# encoding: utf-8
"""
Interface for accessing SDSS imaging data

Notes
-----
Requires utilities from the astrometry.net codebase (trac.astrometry.net) in
Python path.

TODO
----
- Incorporate astrometry.net tools into distribution

History
-------
2011-02-10 - Created by Dan Foreman-Mackey
2011-06-12 - Updated to be used as interface file for calibration project

"""


__all__ = ['SDSSDASFileError','SDSSField']

import os
import os.path
import subprocess
import hashlib

import numpy as np
import pyfits
import h5py

import astrometry.util.sdss_psf as sdss_psf
from astrometry.sdss.common import *
from astrometry.sdss import DR7

from opt import das_base,scratch_base,casdb

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
        
        if filetype in ['fpC']:
            return '%(run)i/%(rerun)i/corr/%(camcol)i/' % kwargs + fn
        elif filetype in ['psField', 'fpAtlas', 'fpObjc', 'fpM']:
            return '%(run)i/%(rerun)i/objcs/%(camcol)i/' % kwargs + fn
        elif filetype in ['tsObj', 'tsField']:
            return '%(run)i/%(rerun)i/calibChunks/%(camcol)i/' % kwargs + fn
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
        p = self._open(fn)
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
                return pyfits.open(local_path)
            except:
                das_path = os.path.join(das_base, fn)
                
                file_dir = os.path.join(os.path.split(local_path)[:-1])[0]
                if os.path.exists(file_dir) is False:
                    os.makedirs(file_dir)
                
                # wget -nH das_path -O local_path
                ret = subprocess.Popen("wget -nH %s -O %s"%(das_path,local_path),
                        shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).wait()
                if ret is not 0:
                    os.remove(local_path)
                    raise SDSSDASFileError(das_path)
                
                return pyfits.open(local_path)
        else:
            path = os.path.join(das_base, fn)
            return pyfits.open(path)

class SDSSField:
    """
    Wrapper class for SDSS image file with various convenience functions
    
    Parameters
    ----------
    run : int
        SDSS run number

    camcol : int
        Camera column

    Raises
    ------
    SDSSDASFileError :
        If any of the requested files are not available locally or on DAS
    
    TODO
    ----
    - Decide what to do about multiple bands

    History
    -------
    2011-06-13 - Created by Dan Foreman-Mackey
    
    """
    def __init__(self,run,camcol,band='g'):
        self.run = run
        self.camcol = camcol

        self.band = band
        band_id = 'ugriz'.index(b)

        self.fields = []
        self.means  = [] # mean RA/Dec points in fields
        self.decmin,self.decmax = None,None
        self.ramin,self.ramax   = None,None
        for doc in casdb.find({'run': run,'camcol': camcol}):
            field = doc['field']
            rerun = doc['rerun']
            sdss    = DFMDR7()
            tsField = sdss.readTsField(run, camcol, field, rerun)
            ast     = tsField.getAsTrans(band)
            (psField,psfData) = sdss.readPsField(run, camcol, field)
            fpCPath = sdss.getFilename('fpC', run, camcol, field, band)
            fpCPath += '.gz' # always zipped?
            fpCfits = sdss._open(fpCPath)
            fpC = FpC(run, camcol, field, band)
            fpC.image = fpCfits[0].data
            img = fpC
            
            # mask file
            fpM = sdss.readFpM(run, camcol, field, band)
            
            # inverse variance
            inv = self.sdss.getInvvar(fpC.getImage(), fpM,
                                    psField.getGain(band_id),
                                    psField.getDarkVariance(band_id),
                                    psField.getSky(band_id),
                                    psField.getSkyErr(band_id))

            self.fields.append({'ast': ast,
                                'img': img,
                                'inv': inv,
                            'psfData': psfData})

            # approximate astrometry



    def hash(self):
        """
        NAME:
            hash
        PURPOSE:
            return an unique identifier for the field
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 06, 2011
        """
        rep = [self.run,self.camcol]
        rep.append(self.band)
        print "-".join([str(r) for r in rep])
        return hashlib.md5("-".join([str(r) for r in rep])).hexdigest()

    def radec_to_pixel(self, ra, dec, color=None):
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
            returns a dictionary of (x,y) tuples with keys from each member of
            self.bands
        HISTORY:
            Created by Dan Foreman-Mackey on Jun 06, 2011
        """
        if color == None:
            color = 0.0
        
        # color will probably always pass 0.0...
        # this is bad for blue objects
        return self.fields[i]['ast'].radec_to_pixel(ra, dec, color=color)
    
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
        psf = []
        xy = self.radec_to_pixel(ra,dec)
        for b in self.bands:
            (x,y) = xy[b]
            psf.append((b,self.psf_at_point(x,y,b)))
        
        return dict(psf)
    
    def psf_at_point(self, x, y, band, radius=25, dblgauss=True):
        """
        NAME:
            psf_at_point
        PURPOSE:
            return the psf image centered at (x,y)
        INPUT:
            x,y  - (float,float) pixel positions
            band - (str) SDSS observation band
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
        if dblgauss: # Double Gaussian
            a,s1,b,s2 = self.psField.getDoubleGaussian('ugriz'.index(band))
            a /= (s1**2 + b*s2**2)
            psf = np.zeros((2*radius+1,2*radius+1))
            xy = np.arange(-radius,radius+1)
            rad = np.sqrt((xy[np.newaxis,:])**2+(xy[:,np.newaxis])**2)
            psf += a*(np.exp(-rad**2/2.0/s1**2) + b*np.exp(-rad**2/2.0/s2**2))/(2*np.pi)
            return psf
        
        # KL...
        psf = sdss_psf.sdss_psf_at_points(self.psfData['ugriz'.index(band)+1], x, y)
        
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
        
        xy = self.radec_to_pixel(ra,dec)
        for b in self.bands:
            (x,y) = xy[b]
            x0,y0 = int(x),int(y)
            
            psf = self.psf_at_point(x, y, b)
            
            xmin,xmax = x0-delta,x0+delta
            ymin,ymax = y0-delta,y0+delta
            
            if xmin < 0:
                xmin = 0
            if ymin < 0:
                ymin = 0
            if xmax >= np.shape(self.img[b].getImage())[1]:
                xmax = np.shape(self.img[b].getImage())[1]-1
            if ymax >= np.shape(self.img[b].getImage())[0]:
                ymax = np.shape(self.img[b].getImage())[0]-1
            if xmin >= xmax or ymin >= ymax:
                raise SDSSOutOfBounds('Not in bounds')
            
            # trim image and inv. variance
            img = self.img[b].getImage()[ymin:ymax,xmin:xmax]
            inv = self.inv[b][ymin:ymax,xmin:xmax]
            img_out.append((b,img))
            inv_out.append((b,inv))
            
            # trim psf
            psf = psf[ymin-(y0-delta):ymax-(y0-delta), \
                    xmin-(x0-delta):xmax-(x0-delta)]
            psf_out.append((b,psf))
            
            # do fit
            res.append((b,self.fit_psf_leastsq(img,inv,psf)))
        
        return dict(res), dict(img_out), dict(inv_out), dict(psf_out)
    
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


