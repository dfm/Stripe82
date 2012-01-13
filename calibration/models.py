#!/usr/bin/env python
# encoding: utf-8
"""
The data model with MongoDB as a backend

"""

__all__ = ['Model', 'CalibRun', 'CalibPatch', 'Star', 'Field']

from datetime import datetime
import cPickle as pickle

from multiprocessing import Pool

import numpy as np
import numpy.ma as ma

import pymongo
from pymongo.code import Code
from pymongo.objectid import ObjectId
from pymongo.son_manipulator import SONManipulator
from bson.binary import Binary

from sdss import SDSSRun, SDSSOutOfBounds
from patchmodel import PatchProbModel, PatchMedianModel
from conversions import *
from config import TESTING
import gp


_connection = pymongo.Connection()

class NumpySONManipulator(SONManipulator):
    """
    SON manipulator for dealing with NumPy arrays

    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html#automatic-encoding-and-decoding

    """
    def transform_incoming(self, value, collection):
        if isinstance(value, (list,tuple,set)):
            return [self.transform_incoming(item,collection) for item in value]
        if isinstance(value,dict):
            return dict((key,self.transform_incoming(item,collection))
                                         for key,item in value.iteritems())
        if isinstance(value,np.ndarray):
            # return {'_type': 'np.ndarray', 'data': value.tolist()}
            return {'_type': 'np.ndarray',
                     'data': Binary(pickle.dumps(value,-1))}
        return value

    def transform_outgoing(self, son, collection):
        if isinstance(son,(list,tuple,set)):
            return [self.transform_outgoing(value,collection) for value in son]
        if isinstance(son,dict):
            if son.get('_type') == 'np.ndarray':
                # return np.array(son.get('data'))
                return pickle.loads(son.get('data'))
            return dict((key,self.transform_outgoing(value,collection))
                                         for key,value in son.iteritems())
        return son

class Model(object):
    """
    The generic data model subclass

    Parameters
    ----------
    _id : str or pymongo.ObjectId, optional
        The _id used to select an existing object

    doc : dict
        A pymongo document to use to construct the object

    band : str, optional
        The SDSS band pass used

    """
    _collection = None

    def __init__(self, _id=None, doc=None, band='g'):
        if doc is not None:
            self.doc = doc
        elif _id is not None:
            if isinstance(_id, str):
                _id = ObjectId(_id)

            doc = self._collection.find_one({'_id': _id})
            assert(doc is not None)
            self.doc = doc
        else:
            self._id = _id
            self.date_created = datetime.now()
            self._band = band

        self._doc = self.doc

        self._collection.ensure_index('band')

    def __getitem__(self, k):
        if k in self._doc:
            return self._doc[k]
        return None

    @classmethod
    def find(cls, q={}, sort=None, limit=None):
        """
        Construct a list of Model objects based on a particular query

        Parameters
        ----------
        q : dict, optional
            pymongo query

        sort : str or list, optional
            The pymongo sort query

        limit : int, optional
            Limit the number of objects returned

        Returns
        -------
        obj : list of Model objects
            The constructed objects or None if no documents matched q

        """
        docs = cls._collection.find(q)
        if sort is not None:
            docs = docs.sort(sort)
        if limit is not None:
            docs = docs.limit(limit)
        r = [cls(doc=doc) for doc in docs]
        if len(r) == 0:
            return None
        return r

    @classmethod
    def find_one(cls, q={}):
        """
        Construct a Model object with a particular query

        Parameters
        ----------
        q : dict
            pymongo query

        Returns
        -------
        obj : Model
            The constructed object or None if no documents matched q

        """
        doc = cls._collection.find_one(q)
        if doc is None:
            return None
        return cls(doc=doc)

    @classmethod
    def find_sphere(cls, coords, radius, c_label=None):
        """
        Construct a list of Model objects based on a particular query

        Parameters
        ----------
        coords : tuple
            A tuple with shape (RA, Dec) in degrees for the center of the search

        radius : float
            Search radius in arcmins

        Returns
        -------
        obj : list of Model objects
            The constructed objects or None if no documents matched the query

        """
        if c_label is None:
            c_label = cls._coord_label

        radius = np.radians(radius/60.0)
        while coords[0] > 180.:
            coords[0] -= 360.0

        q = {c_label: {'$within': {'$centerSphere': [coords,radius]}}}
        docs = cls._collection.find(q)
        if docs is None:
            return None
        return [cls(doc=doc) for doc in docs]

    def __repr__(self):
        return "%s(_id=%s)" % ( type(self).__name__, str(self._id) )

    def __str__(self):
        s = ",\n    ".join(["%s: %s"%(k,str(v)) for k,v in self.doc.iteritems()])
        return """{
    %s
}"""%s

    def dump(self):
        """
        Return a dictionary of values that should be added to th document

        Overridden by subclasses to add values

        Returns
        -------
        properties : dict
            Dictionary to append to the document

        """
        return {}

    def load(self, doc):
        """
        Load a document into the class

        Overridden by subclasses to provide custom rules

        Parameters
        ----------
        doc : dict
            The document

        """
        pass

    @property
    def doc(self):
        doc = {'date_created': self.date_created, 'band': self._band}
        if self._id is not None:
            doc['_id'] = self._id
        for k,v in self.dump().iteritems():
            doc[k] = v
        return doc

    @doc.setter
    def doc(self, doc):
        self._id = doc.pop('_id', None)
        self._band = doc.pop('band', None)
        self.date_created = doc.pop('date_created', datetime.now())
        self.load(doc)

    def save(self):
        doc = self.doc
        doc['date_modified'] = datetime.now()
        self._id = self._collection.save(doc)
        return self._id

#
# Calibration Wrapping Objects
#

class CalibObject(Model):
    _db = _connection.calibration
    _db.add_son_manipulator(NumpySONManipulator())

class CalibRun(CalibObject):
    """
    Object wrapping a 'run' document

    Parameters
    ----------
    run : int
        The run number

    camcol : int
        The camera column

    band : str, optional
        The bandpass to use (default: 'g')

    using_file : str, optional
        If provided, we will try to use the precomputed file

    """
    if TESTING:
        _collection  = CalibObject._db.runs_test
    else:
        _collection  = CalibObject._db.runs

    def __init__(self, *args, **kwargs):
        if len(args) > 1:
            assert('_id' not in kwargs)
            self._run, self._camcol = args
            args = ()

            # make a new run document
            q = {'run': self._run, 'camcol': self._camcol}
            fields = Field.find(q, sort='field')
            assert fields is not None, "There were no fields matching: %s"%(str(q))
            self._fields = [f['_id'] for f in fields]
            band = 'g'
            if 'band' in kwargs:
                band = kwargs['band']
            fn = kwargs.pop('using_file', None)
            if fn is not None:
                self._sdssrun = SDSSRun(fn, band=band)
            else:
                self._sdssrun = SDSSRun(fields, band=band)
            self._filename = self._sdssrun.filename
            self._patches = []

        super(CalibRun, self).__init__(*args, **kwargs)

    def __getattr__(self, name):
        return getattr(self._sdssrun, name)

    def __repr__(self):
        if self._id is not None:
            return "%s(_id=%s)" % ( type(self).__name__, repr(self._id) )
        return "%s(%s, %s)" % \
                (type(self).__name__, repr(self._run), repr(self._camcol))

    def __str__(self):
        return "%s(run = %f, camcol = %f)" % \
                (type(self).__name__, self._run, self._camcol)

    @classmethod
    def find_coords(cls, ra, dec, band='g'):
        """
        Find all the runs that overlap a given coordinate

        Parameters
        ----------
        ra : float
            R.A. in degrees

        dec : flaot
            Dec. in degrees

        Returns
        -------
        runs : list of CalibRun objects
            A list of the matching CalibRun objects or None if there are no
            matches

        """
        runs = []
        runcamcols = Field.find_coords(ra, dec)
        for runcamcol in runcamcols:
            if runcamcol is None:
                return None
            runcamcol['band'] = band
            run = cls._collection.find_one(runcamcol)
            if run is None:
                # run object hasn't been created yet...
                runs.append(cls(runcamcol['run'], runcamcol['camcol']))
            else:
                if 'failed' not in run:
                    runs.append(cls(doc=run))

        return runs

    def dump(self):
        doc = {'run': self._run, 'camcol': self._camcol}
        doc['filename'] = self._filename
        doc['fields'] = self._fields
        doc['patches'] = self._patches
        return doc

    def load(self, doc):
        self._run = doc['run']
        self._camcol = doc['camcol']
        self._filename = doc['filename']
        self._fields = doc['fields']
        self._patches = doc.pop('patches', [])
        self._ras = doc.pop('ras', [])
        self._zeros = doc.pop('zeros', [])
        self._zeros0 = doc.pop('zeros0', [])
        self._sdssrun = SDSSRun(self._filename, band=self._band)

    @classmethod
    def generate(cls, band='g'):
        # cls._collection.drop()
        unique_runs = [r['_id'] for r in UniqueRun._collection.find(\
                fields={'_id': 1})]
        for i in range(len(unique_runs)):
            unique_runs[i]['band'] = band

        pool = Pool(4)
        pool.map(_build_run, unique_runs)

    def fit_gp(self):
        x0, y0 = np.array(self._ras), np.array(self._zeros)
        # MAGIC MAGIC MAGIC
        inds = (y0 > 10) * (~np.isnan(y0))
        x0, y0 = x0[inds], y0[inds]
        self._gp = gp.GaussianProcess()
        self._gp.train(x0, y0)

    def sample_gp(self, x, N=500):
        try:
            return self._gp.sample(x, N)
        except AttributeError as e:
            print "You need to run fit_gp first!"
            raise(e)

def _build_run(info):
    try:
        run = CalibRun(info['run'], info['camcol'], band=info['band'])
    except Exception as e:
        print "Couldn't generate run..."
        print e
        CalibRun._collection.insert( \
                {'run': info['run'], 'camcol': info['camcol'], 'band': info['band'],
                    'failed': True})
    else:
        run.save()

class CalibPatch(CalibObject):
    """
    Describes a patch on the sky where the calibration will be performed

    Parameters
    ----------
    ra : float
        In degrees

    dec : float
        In degrees

    radius : float
        Search radius in arcmin

    model_id : int, optional
        1 - probablistic model, 0 - median model

    """
    if TESTING:
        _collection  = CalibObject._db.patches_test
    else:
        _collection  = CalibObject._db.patches

    _coord_label = 'coords'
    _collection.ensure_index([(_coord_label, pymongo.GEO2D)])
    _collection.ensure_index('runs')
    _collection.ensure_index('stars')

    def __init__(self, *args, **kwargs):
        assert(len(args) in (1,3) or (len(args) == 0 and
            ('_id' in kwargs or 'doc' in kwargs)))
        if len(args) > 1:
            assert('_id' not in kwargs)
            self._ra, self._dec, self._radius = args
            args = ()

            self._model_id = kwargs.pop('model_id', 1)
            band = kwargs.pop('band', 'g')
            kwargs['band'] = band

            self._stars = Star.find_sphere([self._ra, self._dec], self._radius)
            self._star_ids = [s._id for s in self._stars]
            self._runs  = CalibRun.find_coords(self._ra, self._dec, band=band)
            self._run_ids = [r._id for r in self._runs]

            self._t = np.array([r.mjd_at_radec(self._ra, self._dec) for r in self._runs])

            # some of the observations have a messed up timestamp... don't use them!
            self._runs = [self._runs[i] for i in np.where(self._t > 1)[0]]
            self._t = self._t[self._t > 1]

            self._f = np.zeros([len(self._runs), len(self._stars)])
            self._ivar_f = np.zeros(self._f.shape)
            self._calib_mag = np.zeros(len(self._stars))

            # do the photometry
            for ri,run in enumerate(self._runs):
                for si,star in enumerate(self._stars):
                    flux, ivar = star.do_photometry_in_run(run)
                    self._f[ri, si] = flux
                    self._ivar_f[ri, si] = ivar
                    self._calib_mag[si] = getattr(star, run._band)

            # check for stars with < 5?? observations
            crit = np.sum(self._f > 0, axis=0) > 5
            self._f = self._f[:,crit]
            self._ivar_f = self._ivar_f[:, crit]
            self._calib_mag = self._calib_mag[crit]
            for i in np.arange(len(self._stars))[~crit][::-1]:
                self._stars.pop(i)
                self._star_ids.pop(i)

            # also, negatives are bad for us too!
            crit = np.sum(self._f < 0, axis=-1) == 0
            self._f = self._f[crit,:]
            self._ivar_f = self._ivar_f[crit, :]
            self._t = self._t[crit]
            for i in np.arange(len(self._runs))[~crit][::-1]:
                self._runs.pop(i)
                self._run_ids.pop(i)

            self.calibrate()

        super(CalibPatch, self).__init__(*args, **kwargs)

    def __repr__(self):
        if self._id is not None:
            return "%s(_id=%s)" % ( type(self).__name__, str(self._id) )
        return "%s(%s, %s, %s)" % ( type(self).__name__,
                repr(self._ra), repr(self._dec), repr(self._radius) )

    def __str__(self):
        return "Calibration patch @ (ra, dec) = (%f, %f) w/ radius = %f" % \
                (self._ra, self._dec, self._radius)

    def dump(self):
        doc = {self._coord_label: [self._ra,self._dec], 'radius': self._radius}
        doc['stars'] = [star._id for star in self._stars]
        doc['runs']  = [run._id for run in self._runs]
        doc['model_id'] = self._model_id

        doc['t']         = Binary(pickle.dumps(self._t, -1))
        doc['f']         = Binary(pickle.dumps(self._f, -1))
        doc['ivar_f']    = Binary(pickle.dumps(self._ivar_f, -1))
        doc['calib_mag'] = Binary(pickle.dumps(self._calib_mag, -1))

        doc['zeros'] = self._zeros
        doc['star_mag'] = self._mstar

        doc['zeros0'] = self._zeros0
        doc['star_mag0'] = self._mstar0

        doc['nuisance'] = self._nuisance

        return doc

    def load(self, doc):
        self._model_id = doc.pop('model_id', 1)
        self._ra,self._dec = doc[self._coord_label]
        self._radius = doc['radius']

        self._stars = [Star(s) for s in doc['stars']]
        self._star_ids = doc['stars']
        self._runs  = [CalibRun(s) for s in doc['runs']]
        self._run_ids = doc['runs']

        self._t = pickle.loads(doc['t'])
        self._f = pickle.loads(doc['f'])
        self._ivar_f = pickle.loads(doc['ivar_f'])
        self._calib_mag = pickle.loads(doc['calib_mag'])

        self._zeros  = doc['zeros']
        self._mstar  = doc['star_mag']
        self._zeros0 = doc['zeros0']
        self._mstar0 = doc['star_mag0']

        self._nuisance = doc['nuisance']

    def calibrate(self):
        model0 = PatchMedianModel(self._t, self._f, self._ivar_f, self._calib_mag)
        model0.calibrate()
        self._zeros0 = model0.zeros
        self._mstar0 = model0.star_mag

        model = PatchProbModel(self._t, self._f, self._ivar_f, self._calib_mag)
        model._f0 = self._zeros0
        model._mstar = self._mstar0
        model.calibrate()
        self._zeros = model.zeros
        self._mstar = model.star_mag
        self._nuisance = model.nuisance

    def update_runs(self):
        """
        Update the zero point measurements in the run documents hit by this patch

        """
        [CalibRun._collection.update({'_id': run._id},
                {'$push': {'patches': self._id, 'ras': self._ra,
                    'zeros': self.zero_for_run(run),
                    'zeros0': self.zero_for_run(run,model=0)}})
            for run in self._runs]

    def zero_for_run(self, run, model=1):
        try:
            if model == 1:
                return self._zeros[self._run_ids.index(run._id)]
            return self._zeros0[self._run_ids.index(run._id)]
        except ValueError:
            return None

    def lightcurve(self, star):
        try:
            sid = self._star_ids.index(star._id)
        except ValueError:
            return None
        else:
            flux = self._f[:,sid]/self._zeros
            inv  = self._ivar_f[:,sid]

            mask = inv <= 1e-10
            inv = ma.array(inv, mask=mask)
            err = np.sqrt(1.0/inv+self._nuisance[4]+ \
                self._nuisance[5]*(self._zeros*mag2nmgy(self._mstar[sid]))**2)/\
                self._zeros
            err[mask] = np.inf

            return self._t, flux, err

#
# Data Wrapping Objects
#

class SDSSObject(Model):
    _bands = [b for b in 'ugriz']
    _db = _connection.cas
    _db.add_son_manipulator(NumpySONManipulator())

class Star(SDSSObject):
    if TESTING:
        _collection = SDSSObject._db.stars_test
    else:
        _collection = SDSSObject._db.stars
    _coord_label = 'pos'
    _collection.ensure_index([(_coord_label, pymongo.GEO2D)])
    [_collection.ensure_index('photo_%s'%(b)) for b in SDSSObject._bands]

    def dump(self):
        doc = {self._coord_label: [self._ra,self._dec],
                'lyrae_candidate': self._lyrae_candidate,
                'rank': self._rank}
        for b in self._bands:
            doc[b] = getattr(self, b)
            doc['photo_%s'%b] = getattr(self, '_photo_%s'%b)
        return doc

    def load(self, doc):
        self._ra,self._dec = doc.pop(self._coord_label)
        self._lyrae_candidate = doc.pop('lyrae_candidate', False)
        self._rank = doc.pop('rank')
        for b in self._bands:
            setattr(self, b, doc[b])
            setattr(self, '_photo_%s'%b, doc.pop('photo_%s'%b, {}))
            d = getattr(self, '_photo_%s'%b)
            for k in list(d):
                d[k] = pickle.loads(d[k])

    def do_photometry_in_run(self, run):
        self_photo = getattr(self, '_photo_%s'%(run._band))
        rid = str(run._id)
        if rid not in self_photo:
            # measurement is not already in the database
            try:
                photo, cov = run.photo_at_radec(self._ra, self._dec)
                inv = 1.0/cov.diagonal()
            except SDSSOutOfBounds:
                ndim = 4
                photo, inv = np.zeros(ndim), np.zeros(ndim)
                cov = None

            # update self
            self_photo[rid] = (photo, cov, inv)
            # update database
            self._collection.update({'_id': self._id},
                    {'$set': {'photo_%s.%s'%(run._band, rid):
                        Binary(pickle.dumps((photo, cov, inv),-1))}})
        else:
            photo, cov, inv = self_photo[rid]

        return photo[1], inv[1]

    @classmethod
    def find_with_photometry_in_run(cls, run):
        cursor = cls._collection.find({"photo_%s.%s"%(run._band, run._id): {"$exists": True}}, {"_id": 1})
        return [Star(d["_id"]) for d in cursor]

class Field(SDSSObject):
    if TESTING:
        _collection = SDSSObject._db.fields_test
    else:
        _collection = SDSSObject._db.fields
    _collection.ensure_index([('ramin', pymongo.ASCENDING),
        ('ramax', pymongo.ASCENDING), ('decmin', pymongo.ASCENDING),
        ('decmax', pymongo.ASCENDING)])
    _collection.ensure_index([('run', pymongo.ASCENDING),
        ('camcol', pymongo.ASCENDING)])

    _mjd_labels = ['mjd_%s'%b for b in SDSSObject._bands]

    def dump(self):
        doc = {'decmin': self._decmin, 'decmax': self._decmax,
                'ramin': self._ramin, 'ramax': self._ramax,
                'run': self._run, 'camcol': self._camcol,
                'field': self._field, 'rerun': self._rerun}
        for b in self._mjd_labels:
            doc[b] = getattr(self, b)
        return doc

    def load(self, doc):
        self._decmin, self._decmax = doc['decmin'], doc['decmax']
        self._ramin, self._ramax = doc['ramin'], doc['ramax']
        self._run, self._camcol, self._field, self._rerun = \
                doc['run'], doc['camcol'], doc['field'], doc['rerun']
        for b in self._mjd_labels:
            setattr(self, b, doc.pop(b, None))

    @classmethod
    def find_coords(cls, ra, dec):
        """
        Find all the fields that overlap a given coordinate (in Dec only!)

        Parameters
        ----------
        ra : float
            R.A. in degrees

        dec : flaot
            Dec. in degrees

        Returns
        -------
        runcamcols : list
            A list of (unique) matching run/camcol dicts

        """
        cursor = cls._collection.find(
                {'decmin': {'$lt': dec}, 'decmax': {'$gt': dec}},
                {'_id': 0, 'run':1, 'camcol':1})
        results = []
        for doc in cursor:
            if doc not in results:
                # only include the document if there isn't already a matching one
                # stupid unhashable dictionaries
                results.append(doc)
        return results

class UniqueRun(SDSSObject):
    """
    Data model that encapsulates the set of unique runs in the SDSS fields database

    """
    _coll_name = "unique_runs"
    _collection = SDSSObject._db[_coll_name]

    @classmethod
    def generate(cls):
        """
        Generate the database based on the Field dataset

        Returns
        -------
        count : int
            The number of unique run/camcol pairs

        """
        cls._collection.drop()
        m = Code("function () {"
                 "    emit({run: this.run, camcol: this.camcol}, 1);"
                 "}")
        r = Code("function (key, values) {"
                 "    var total = 0;"
                 "    for (var i = 0; i < values.length; i++) {"
                 "        total += values[i];"
                 "    }"
                 "    return total;"
                 "}")

        result = Field._collection.map_reduce(m,r, cls._coll_name)
        return result.count()

if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rc('text', usetex=True)
    import matplotlib.pyplot as pl

    import os
    import os.path

    band = 'i'

    runs = CalibRun.find({'failed': {'$exists': False}, 'band': band})
    # runs = [CalibRun.find_one({'run': 6537, 'camcol': 4, 'band': 'g'})]
    for i, run in enumerate(runs):
        pl.clf()

        try:
            run.fit_gp()
            x = np.linspace(-1,1,500)
            y = run.sample_gp(x)
        except Exception as e:
            print "Exception"
            print repr(e)
            print e
            pl.title('run: %d, camcol: %d, Singular'%(run._run, run._camcol))
        else:
            print i, run._gp._la2, run._gp._a2, \
                    run._gp._lb2, run._gp._b2, run._gp._s2
            mu = run.mean_gp(x)
            pl.plot(x, y,'k',alpha=0.05)
            pl.plot(x,mu,'r')
            pl.title('run: %d, camcol: %d, $L$: %.4f'%(run._run, run._camcol,
                np.sqrt(run._gp._lb2)))

            an = "\n".join(["$%s = %.5f$"%(k[1:-1].upper(), np.sqrt(getattr(run._gp,k)))\
                    for k in ['_s2', '_a2', '_la2', '_b2', '_lb2']])
            pl.gca().annotate(an, xy=(0.05,0.05),  xycoords='axes fraction',
                xytext=(0, 0), textcoords='offset points',
                bbox=dict(fc="w",alpha=0.25),
                size=9,alpha=0.5)

        x0, y0 = np.array(run._ras), np.array(run._zeros)
        y2 = np.array(run._zeros0)
        inds = y0 > 10.0
        x0, y0 = x0[inds], y0[inds]
        y2 = y2[inds]

        pl.plot(x0,y2,'og', alpha=0.2)
        pl.plot(x0,y0,'.b')

        mu = np.median(y0)
        five = mu*0.05
        pl.gca().axhline(mu,color='b', alpha=0.5)
        pl.gca().axhline(mu+five,color='b',ls='--', alpha=0.5)
        pl.gca().axhline(mu-five,color='b',ls='--', alpha=0.5)

        try:
            pl.plot([-0.9, -0.9],
                    np.median(y0) + np.sqrt(run._gp._s2)*np.array([1,-1]),
                    'r', lw=2.0)
        except Exception as e:
            print "Couldn't show error bar"
            print e

        pl.xlabel('R.A.')
        pl.ylabel('ADU nMgy$^{-1}$')

        bfn = '/home/dfm265/public_html/zero_test/%s'%band
        ffn = os.path.join(bfn, "%d"%(run._run))
        print ffn
        if not os.path.exists(ffn):
            os.makedirs(ffn)
        pl.savefig(os.path.join(ffn, '%d.png'%(run._camcol)))

