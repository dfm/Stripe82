#!/usr/bin/env python
# encoding: utf-8
"""
The data model with MongoDB as a backend

"""

__all__ = ['Model', 'CalibRun', 'CalibPatch']

from datetime import datetime
import cPickle as pickle

import numpy as np

import pymongo
from pymongo.objectid import ObjectId
from pymongo.son_manipulator import SONManipulator
from bson.binary import Binary

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
            return {'_type': 'np.ndarray', 'data': value.tolist()}
            # return {'_type': 'np.ndarray',
            #         'data': Binary(pickle.dumps(value,-1))}
        return value

    def transform_outgoing(self, son, collection):
        if isinstance(son,(list,tuple,set)):
            return [self.transform_outgoing(value,collection) for value in son]
        if isinstance(son,dict):
            if son.get('_type') == 'np.ndarray':
                return np.array(son.get('data'))
                # return pickle.loads(son.get('data'))
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


    @classmethod
    def find(cls, q={}):
        """
        Construct a list of Model objects based on a particular query

        Parameters
        ----------
        q : dict
            pymongo query

        Returns
        -------
        obj : list of Model objects
            The constructed objects or None if no documents matched q

        """
        docs = cls._collection.find(q)
        if docs is None:
            return None
        return [cls(doc=doc) for doc in docs]

    @classmethod
    def find_one(cls, q):
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
        doc = {'date_created': self.date_created}
        if self._id is not None:
            doc['_id'] = self._id
        for k,v in self.dump().iteritems():
            doc[k] = v
        return doc

    @doc.setter
    def doc(self, doc):
        for k in ['_id', 'date_created']:
            setattr(self, k, doc.pop(k))
        self.load(doc)

    def save(self):
        doc = self.doc
        doc['date_modified'] = datetime.now()
        self._id = self._collection.insert(doc)

_calib_db = _connection.calibration
_calib_db.add_son_manipulator(NumpySONManipulator())

class CalibRun(Model):
    _collection = _calib_db.runs

    def __init__(self, *args, **kwargs):
        super(CalibRun, self).__init__(*args, **kwargs)

class CalibPatch(Model):
    _collection = _calib_db.patches

    def __init__(self, *args, **kwargs):
        assert(len(args) in (1,3) or (len(args) == 0 and '_id' in kwargs))
        if len(args) > 1:
            assert('_id' not in kwargs)
            self._ra, self._dec, self._radius = args
            args = ()
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
        return {'coords': [self._ra,self._dec], 'radius': self._radius}

    def load(self, doc):
        self._ra,self._dec = doc['coords']
        self._radius = doc['radius']

