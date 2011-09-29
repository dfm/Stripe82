#!/usr/bin/env python
# encoding: utf-8
"""
The data model with MongoDB as a backend

"""

__all__ = ['Model']

from datetime import datetime

import pymongo
from pymongo.objectid import ObjectId

_connection = pymongo.Connection()

class Model(object):
    """
    The generic data model subclass

    Parameters
    ----------
    _id : str or pymongo.ObjectId, optional
        The _id used to select an existing object

    """
    _collection = None

    def __init__(self, _id=None):
        if _id is None:
            self._id = _id
            self.date_created = datetime.now()
        else:
            if isinstance(_id, str):
                _id = ObjectId(_id)
            assert(isinstance(_id, ObjectId))

            doc = self._collection.find_one({'_id': _id})
            assert(doc is not None)
            self.doc = doc

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

class CalibRun(Model):
    _collection = _connection.calibration.runs

    def __init__(self, *args, **kwargs):
        super(CalibRun, self).__init__(*args, **kwargs)

class CalibPatch(Model):
    _collection = _connection.calibration.patches

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

