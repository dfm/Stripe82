"""
Interact with SDSS data saved in an HDF5 file generated by the preprocess
script.

"""

__all__ = ['']

import os

import numpy as np
import scipy.optimize as op
import h5py

# Where are the data files stored locally?
_local_data_dir  = os.path.join(os.environ.get("SDSS_LOCAL", "."), "data")

# HDF5 file tags.
IMG_TAG         = "img"
INV_TAG         = "inv"
TAI_TAG         = "tai"
BOUNDS_TAG      = "bounds"
CENTERS_TAG     = "centers"
PSF_TAG         = "psf"
EIGEN_TAG       = "eigen"
INFO_TAG        = "info"
TS_TAG          = "ts"
AST_TAG         = "ast"

# Field dimensions.
_f_width   = 1489
_f_height  = 2048
_f_overlap = 128

def get_filename(run, camcol, band):
    """
    Note: make sure that this stays in sync with the function in the
    preprocess script!

    """
    return os.path.join(_local_data_dir, "%d-%d-%s.hdf5"%(run, camcol, band))

class _Astrometry(object):
    """
    Compute the astronometric transformations for an SDSS field. Ported
    from code written by Dustin Lang (Princeton). See the
    [original code](https://github.com/astrometry/pysdss).

    """
    def __init__(self, band_id, node, incl, data):
        self.node = node
        self.incl = incl
        self.trans = {}
        for f in ['a','b','c','d','e','f', 'riCut',
                    'dRow0', 'dRow1', 'dRow2', 'dRow3',
                    'dCol0', 'dCol1', 'dCol2', 'dCol3',
                    'csRow', 'csCol', 'ccRow', 'ccCol']:
            self.trans[f] = data[f][0][band_id]

    def _get_abcdef(self):
        return [self.trans[x] for x in "abcdef"]

    def _get_drow(self):
        return [self.trans[x] for x in ["dRow0", "dRow1", "dRow2", "dRow3"]]

    def _get_dcol(self):
        return [self.trans[x] for x in ["dCol0", "dCol1", "dCol2", "dCol3"]]

    def _get_cscc(self):
        return [self.trans[x] for x in ["csRow", "csCol", "ccRow", "ccCol"]]

    def _get_ricut(self):
        return self.trans['riCut']

    def _get_ghpq(self):
        # FIXME: I'm not sure if this is right... --DFM
        g0, g1, g2, g3 = self._get_drow()
        h0, h1, h2, h3 = self._get_dcol()
        px, py, qx, qy = self._get_cscc()

        # #$(%*&^(%$%*& bad documentation. --Lang
        (px,py) = (py,px)
        (qx,qy) = (qy,qx)

        return g0, g1, g2, g3, h0, h1, h2, h3, px, py, qx, qy

    def pixel_to_radec(self, x, y, color=0):
        mu, nu = self.pixel_to_munu(x, y, color)
        return self.munu_to_radec(mu, nu)

    def radec_to_pixel(self, ra, dec, color=0):
        mu, nu = self.radec_to_munu(ra, dec)
        return self.munu_to_pixel(mu, nu, color)

    def munu_to_pixel(self, mu, nu, color=0):
        xprime, yprime = self.munu_to_prime(mu, nu, color)
        return self.prime_to_pixel(xprime, yprime)

    def munu_to_prime(self, mu, nu, color=0):
        a, b, c, d, e, f = self._get_abcdef()
        determinant = b * f - c * e
        B, C, E, F =  np.array([f, -c, -e, b])/determinant
        mua = mu - a
        while mua < -180.:
            mua += 360.
        yprime = B * mua + C * (nu - d)
        xprime = E * mua + F * (nu - d)
        return xprime, yprime

    def pixel_to_munu(self, x, y, color=0):
        xprime, yprime = self.pixel_to_prime(x, y, color)
        a, b, c, d, e, f = self._get_abcdef()
        mu = a + b * yprime + c * xprime
        nu = d + e * yprime + f * xprime
        return mu, nu

    def pixel_to_prime(self, x, y, color=0):
        """
        (Secret decoder ring)
        [http://www.sdss.org/dr7/products/general/astrometry.html]

        * color0 is called `riCut`
        * g0, g1, g2, and g3 are called `dRow0`, `dRow1`, `dRow2` and `dRow3`
        * h0, h1, h2, and h3 are called `dCol0`, `dCol1`, `dCol2` and `dCol3`
        * px and py are called `csRow` and `csCol`
        * qx and qy are called `ccRow` and `ccCol`

        """
        color0 = self._get_ricut()
        g0, g1, g2, g3, h0, h1, h2, h3, px, py, qx, qy = \
                self._get_ghpq()

        yprime = y + g0 + g1 * x + g2 * x**2 + g3 * x**3
        xprime = x + h0 + h1 * x + h2 * x**2 + h3 * x**3

        qx = qx * np.ones_like(x)
        qy = qy * np.ones_like(y)
        xprime += np.where(color < color0, px * color, qx)
        yprime += np.where(color < color0, py * color, qy)

        return xprime, yprime

    def prime_to_pixel(self, xprime, yprime,  color=0):
        color0 = self._get_ricut()
        g0, g1, g2, g3, h0, h1, h2, h3, px, py, qx, qy = \
                self._get_ghpq()

        qx = qx * np.ones_like(xprime)
        qy = qy * np.ones_like(yprime)
        xprime -= np.where(color < color0, px * color, qx)
        yprime -= np.where(color < color0, py * color, qy)

        f  = lambda x: x + h0 + h1 * x + h2 * x**2 + h3 * x**3 - xprime
        fp = lambda x: 1 + h1 + 2 * h2 * x + 3 * h3 * x**2
        x0 = xprime - h0

        x = op.newton(f, x0, fprime=fp)
        y = yprime - (g0 + g1 * x + g2 * x**2 + g3 * x**3)
        return x, y

    def radec_to_munu(self, ra, dec):
        node, incl = self.node, self.incl
        ra, dec = np.radians(ra), np.radians(dec)
        mu = node + np.arctan2(np.sin(ra - node) * np.cos(dec) * np.cos(incl)
                               + np.sin(dec) * np.sin(incl),
                               np.cos(ra - node) * np.cos(dec))
        nu = np.arcsin(-np.sin(ra - node) * np.cos(dec) * np.sin(incl) +
                       np.sin(dec) * np.cos(incl))
        mu, nu = np.degrees(mu), np.degrees(nu)
        mu += (360. * (mu < 0))
        return mu, nu

    def munu_to_radec(self, mu, nu):
        node,incl = self.node, self.incl
        mu, nu = np.radians(mu), np.radians(nu)
        ra = node + np.arctan2(np.sin(mu - node) * np.cos(nu) * np.cos(incl)
                               - np.sin(nu) * np.sin(incl),
                               np.cos(mu - node) * np.cos(nu))
        dec = np.arcsin(np.sin(mu - node) * np.cos(nu) * np.sin(incl)
                        + np.sin(nu) * np.cos(incl))
        ra, dec = np.degrees(ra), np.degrees(dec)
        ra += (360. * (ra < 0))
        return ra, dec

class _Photometry(object):
    """
    Compute the point-spread function for an SDSS field. Ported
    from code written by Dustin Lang (Princeton). See the
    [original code](https://github.com/astrometry/pysdss).

    """
    def __init__(self, info, eigen):
        # The eigenimages.
        self.eigen = eigen

        # Save the metadata.
        t = info[0]

        self.gain           = t["gain"]
        self.dark_variance  = t["dark_variance"]
        self.sky            = t["sky"]
        self.skyerr         = t["skyerr"]
        self.psp_status     = t["status"]

        self.psf_fwhm       = t["psf_width"] * (2.*np.sqrt(2.*np.log(2.)))

        # Parameters for the double Gaussian PSF.
        self.dgpsf_s1 = t["psf_sigma1_2G"]
        self.dgpsf_s2 = t["psf_sigma2_2G"]
        self.dgpsf_b  = t["psf_b_2G"]

    def psf_at_points(self, x, y):
        """
        Reconstruct the SDSS model PSF from KL basis functions. `x` and `y`
        can be scalars or 1D `numpy.ndarray`s.

        """
        rtnscalar = np.isscalar(x) and np.isscalar(y)
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        assert len(x.shape) == 1 and len(y.shape) == 1

        psf = self.eigen
        psfimgs = None
        outh, outw = psf["RNROW"][0], psf["RNCOL"][0]

        # From the IDL docs:
        # http://photo.astro.princeton.edu/photoop_doc.html#SDSS_PSF_RECON
        #   acoeff_k = SUM_i{ SUM_j{ (0.001*ROWC)^i * (0.001*COLC)^j * C_k_ij } }
        #   psfimage = SUM_k{ acoeff_k * RROWS_k }
        for k in range(len(psf)):
            nrb = psf["nrow_b"][k]
            ncb = psf["ncol_b"][k]
            c = psf["c"][k].reshape(5, 5)
            c = c[:nrb,:ncb]
            (gridi,gridj) = np.meshgrid(range(nrb), range(ncb))
            print "RROWS", psf["RROWS"][k]

            if psfimgs is None:
                psfimgs = [np.zeros_like(psf["RROWS"][k])
                                        for xy in np.broadcast(x,y)]
            assert(psf["RNROW"][k] == outh)
            assert(psf["RNCOL"][k] == outw)

            for i,(xi,yi) in enumerate(np.broadcast(x,y)):
                acoeff_k = np.sum(((0.001 * xi)**gridi * (0.001 * yi)**gridj * c))
                psfimgs[i] += acoeff_k * psf["RROWS"][k]
                print acoeff_k
                print psf["RROWS"][k]

        print outh,outw
        print [img for img in psfimgs]

        psfimgs = [img.reshape((outh,outw)) for img in psfimgs]
        if rtnscalar:
            return psfimgs[0]
        return psfimgs

class SDSSRun(object):
    def __init__(self, run, camcol, band):
        self.fn = get_filename(run, camcol, band)
        if not os.path.exists(self.fn):
            raise Exception("The file %s does not exist."%self.fn)

        # Link to the file.
        self.f = h5py.File(self.fn)
        self.img     = self.f[IMG_TAG]
        self.inv     = self.f[INV_TAG]
        self.centers = self.f[CENTERS_TAG]
        self.bounds  = self.f[BOUNDS_TAG]

        # Some properties of the run.
        self.run     = run
        self.camcol  = camcol
        self.band    = band
        self.band_id = "ugriz".index(band)

        self.ast = []
        self.psf = []
        for field in self.f[TS_TAG]:
            info = self.f[TS_TAG][field].attrs
            node = np.radians(info["NODE"])
            incl = np.radians(info["INCL"])
            self.ast.append(_Astrometry(self.band_id, node, incl,
                                                self.f[TS_TAG][field][...]))

            self.psf.append(_Photometry(
                                    self.f[PSF_TAG][field][INFO_TAG][...],
                                    self.f[PSF_TAG][field][EIGEN_TAG][...]))

    def _find_closest_field(self, ra, dec):
        delta = np.sum((np.array([ra, dec]) - self.centers)**2, axis=1)
        return np.argmin(delta)

    def radec_to_pixel(self, ra, dec, fid=None):
        if fid is None:
            fid = self._find_closest_field(ra, dec)

        px, py = self.ast[fid].radec_to_pixel(ra, dec)
        if not (0 < px < _f_height and 0 < py < _f_width):
            raise IndexError("(%.2f, %.2f) is out of bounds."%(ra, dec))
        py += fid * (_f_width - _f_overlap)
        return px, py

    def get_image_patch(self, ra, dec, dim=25, fid=None):
        px, py = self.radec_to_pixel(ra, dec, fid=fid)
        px, py = int(px), int(py)
        shape = self.img.shape

        # Dimensions within the image.
        xmin, xmax = np.max([px-dim, 0]), np.min([px+dim+1, shape[1]])
        ymin, ymax = np.max([py-dim, 0]), np.min([py+dim+1, shape[0]])

        # Dimensions in the target image.
        txmin, tymin = xmin-(px-dim), ymin-(py-dim)
        txmax, tymax = txmin + (xmax-xmin), tymin + (ymax-ymin)

        im = np.zeros((2*dim+1, 2*dim+1))
        im[tymin:tymax, txmin:txmax] = self.img[ymin:ymax, xmin:xmax]

        iv = np.zeros((2*dim+1, 2*dim+1))
        iv[tymin:tymax, txmin:txmax] = self.inv[ymin:ymax, xmin:xmax]

        return im, iv

    def get_psf(self, ra, dec, dim=25, fid=None):
        if fid is None:
            fid = self._find_closest_field(ra, dec)
        px, py = self.psf[fid].psf_at_points(ra, dec)

if __name__ == "__main__":
    run = SDSSRun(4263, 4, "g")

    import matplotlib.pyplot as pl
    img = run.get_image_patch(355.31597205, 0.08818654)
    psf = run.get_psf(355.31597205, 0.08818654)
    pl.imshow(img[0])
    pl.colorbar()
    pl.figure()
    pl.imshow(img[1])
    pl.colorbar()
    pl.show()
