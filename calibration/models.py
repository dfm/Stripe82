#!/usr/bin/env python
# encoding: utf-8
"""
The data model with MongoDB as a backend

"""

__all__ = ['Model', 'CalibRun', 'CalibPatch']

from datetime import datetime
import cPickle as pickle

from multiprocessing import Pool

import numpy as np

import pymongo
from pymongo.code import Code
from pymongo.objectid import ObjectId
from pymongo.son_manipulator import SONManipulator
from bson.binary import Binary

from sdss import SDSSRun
from config import TESTING


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
            assert(isinstance(_id, ObjectId))

            doc = self._collection.find_one({'_id': _id})
            assert(doc is not None)
            self.doc = doc
        else:
            self._id = _id
            self.date_created = datetime.now()
            self._band = band

        self._doc = self.doc

    def __getitem__(self, k):
        if k in self._doc:
            return self._doc[k]
        return None

    @classmethod
    def find(cls, q={}, sort=None):
        """
        Construct a list of Model objects based on a particular query

        Parameters
        ----------
        q : dict, optional
            pymongo query

        sort : str or list, optional
            The pymongo sort query

        Returns
        -------
        obj : list of Model objects
            The constructed objects or None if no documents matched q

        """
        docs = cls._collection.find(q)
        if sort is not None:
            docs = docs.sort(sort)
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
    def find_coords(cls, ra, dec):
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
            run = cls._collection.find_one(runcamcol)
            if run is None:
                # run object hasn't been created yet...
                runs.append(cls(runcamcol['run'], runcamcol['camcol']))
            else:
                runs.append(cls(doc=run))

        return runs

    def dump(self):
        doc = {'run': self._run, 'camcol': self._camcol}
        doc['filename'] = self._filename
        doc['fields'] = self._fields
        doc['patches'] = self._patches
        print self._patches
        return doc

    def load(self, doc):
        self._run = doc['run']
        self._camcol = doc['camcol']
        self._filename = doc['filename']
        self._fields = doc['fields']
        self._patches = doc.pop('patches', [])
        self._sdssrun = SDSSRun(self._filename, band=self._band)

    @classmethod
    def generate(cls):
        cls._collection.drop()
        unique_runs = [r['_id'] for r in UniqueRun._collection.find(\
                fields={'_id': 1})]
        pool = Pool(4)
        pool.map(_build_run, unique_runs)

def _build_run(info):
    try:
        run = CalibRun(info['run'], info['camcol'])
    except Exception as e:
        print "Couldn't generate run..."
        print e
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

    """
    if TESTING:
        _collection  = CalibObject._db.patches_test
    else:
        _collection  = CalibObject._db.patches

    _coord_label = 'coords'
    _collection.ensure_index([(_coord_label, pymongo.GEO2D)])

    def __init__(self, *args, **kwargs):
        assert(len(args) in (1,3) or (len(args) == 0 and '_id' in kwargs))
        if len(args) > 1:
            assert('_id' not in kwargs)
            self._ra, self._dec, self._radius = args
            args = ()
            self._stars = Star.find_sphere([self._ra, self._dec], self._radius)
            self._runs  = CalibRun.find_coords(self._ra, self._dec)

            t = np.array([r.mjd_at_radec(self._ra, self._dec) for r in self._runs])

            print t

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
        doc['stars'] = self._stars
        doc['runs']  = self._runs
        return doc

    def load(self, doc):
        self._ra,self._dec = doc[self._coord_label]
        self._radius = doc['radius']
        self._stars = doc['stars']
        self._runs  = doc['runs']


#
# Data Wrapping Objects
#

class SDSSObject(Model):
    _bands = [b for b in 'ugriz']
    _db = _connection.cas

class Star(SDSSObject):
    if TESTING:
        _collection = SDSSObject._db.stars_test
    else:
        _collection = SDSSObject._db.stars
    _coord_label = 'pos'
    _collection.ensure_index([(_coord_label, pymongo.GEO2D)])

    def dump(self):
        doc = {self._coord_label: [self._ra,self._dec],
                'lyrae_candidate': self._lyrae_candidate,
                'rank': self._rank}
        for b in self._bands:
            doc[b] = getattr(self, b)
        return doc

    def load(self, doc):
        self._ra,self._dec = doc[self._coord_label]
        self._lyrae_candidate = doc['lyrae_candidate']
        self._rank = doc['rank']
        for b in self._bands:
            setattr(self, b, doc[b])

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
        Find all the fields that overlap a given coordinate

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
                {'ramin': {'$lt': ra}, 'ramax': {'$gt': ra},
                        'decmin': {'$lt': dec}, 'decmax': {'$gt': dec}},
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

