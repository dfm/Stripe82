#!/usr/bin/env python
# encoding: utf-8
"""


History
-------
2011-07-06 - Created by Dan Foreman-Mackey

"""

__all__ = ['PhotoSONManipulator']

import cPickle as pickle

import numpy as np

from bson.binary import Binary
from pymongo.son_manipulator import SONManipulator

from photomodel import PhotoData,PhotoModel

# ================#
#  BSON Encodings #
# ================#

class PhotoSONManipulator(SONManipulator):
    """
    SON manipulator to deal with my custom data types AND multi-D NumPy arrays
    
    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html#automatic-encoding-and-decoding
    [2] https://github.com/FlaPer87/django-mongodb-engine/blob/master/django_mongodb_engine/serializer.py
    
    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    def transform_incoming(self, value, collection):
        if isinstance(value, (list,tuple,set)):
            return [self.transform_incoming(item,collection) for item in value]
        if isinstance(value,dict):
            return dict((key,self.transform_incoming(item,collection))
                    for key,item in value.iteritems())
        if isinstance(value,PhotoModel):
            return encode_model(value)
        if isinstance(value,PhotoData):
            return encode_data(value)
        if isinstance(value,np.ndarray):
            return {'_type': 'np.ndarray', 
                    'data': Binary(pickle.dumps(value,-1))}
        return value

    def transform_outgoing(self, son, collection):
        if isinstance(son,(list,tuple,set)):
            return [self.transform_outgoing(value,collection) for value in son]
        if isinstance(son,dict):
            if son.get('_type') == 'PhotoModel':
                return decode_model(son)
            if son.get('_type') == 'PhotoData':
                return decode_data(son)
            if son.get('_type') == 'np.ndarray':
                return pickle.loads(son.get('data'))
            return dict((key,self.transform_outgoing(value,collection))
                    for key,value in son.iteritems())
        return son

# ============================ #
#  Light curve model encodings #
# ============================ #

def encode_model(model):
    """
    Encode the model parameters for BSON

    Parameters
    ----------
    model : PhotoModel
        The model object
    
    Returns
    -------
    encoded : dict
        Encoded version of the class

    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html
    
    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    return {'_type': 'PhotoModel', 'vector': list(model.vector()),
            'data': encode_data(model.data)}

def decode_model(document):
    """
    Decode a BSON PhotoModel document
    
    Parameters
    ----------
    document : dict
        BSON document
    
    Returns
    -------
    model : PhotoModel
        The model object
    
    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html

    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    assert document['_type'] == 'PhotoModel'
    return PhotoModel(decode_data(document['data']),document['vector'])

def encode_data(data):
    """
    Encode the PhotoData class for BSON
    
    Parameters
    ----------
    data : PhotoData
        The data object
    
    Returns
    -------
    encoded : dict
        Save-able by BSON
    
    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html

    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    return {'_type': 'PhotoData','data': Binary(pickle.dumps(data.data,-1)),
            'observations': [obs['_id'] for obs in data.observations],
            'stars': [star['_id'] for star in data.stars]}

def decode_data(document):
    """
    Decode a BSON PhotoData document
    
    Parameters
    ----------
    document : dict
        BSON document
    
    Returns
    -------
    data : PhotoData
        The data object
    
    References
    ----------
    [1] http://api.mongodb.org/python/current/examples/custom_type.html

    History
    -------
    2011-06-15 - Created by Dan Foreman-Mackey
    
    """
    assert document['_type'] == 'PhotoData'
    return PhotoData(pickle.loads(document['data']),document['observations'],
        document['stars'])


