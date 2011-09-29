#!/usr/bin/env python
# encoding: utf-8
"""
The data model with MongoDB as a backend

"""

__all__ = ['Model']

from datetime import datetime

import pymongo
from pymongo.objectid import ObjectId

class Model(object):
    """
    The generic data model subclass

    Parameters
    ----------
    query : dict or pymongo.ObjectId, optional
        The query used to select an existing object

    """
    _collection = pymongo.Connection().data.model

    def __init__(self, _id=None):
        if _id is not None:
            if isinstance(_id, str):
                _id = ObjectId(_id)
            assert(isinstance(_id, ObjectId))

            doc = self._collection.find_one({'_id': _id})
            assert(doc is not None)
            self.doc = doc

        self._id = _id
        self.date_created = datetime.now()

    def __repr__(self):
        return "%s(_id=%s)" % ( "Model", str(self._id) )

    @property
    def doc(self):
        doc = {'date_created': self.date_created}
        if self._id is not None:
            doc['_id'] = self._id
        return doc

    @doc.setter
    def doc(self, doc):
        for k, v in doc.iteritems():
            setattr(self, k, doc[k])

    def save(self):
        doc = self.doc
        doc['date_modified'] = datetime.now()
        self._id = self._collection.insert(doc)

if __name__ == '__main__':
    Model()


