"""


"""

__all__ = ["Model", "SDSSModel"]

import os

from db import Database

_db = Database()
_sdss_db = Database(name=os.environ.get("MONGO_DB", "sdss"))

class Model(object):
    cname  = None
    fields = []
    coords = None

    def __init__(self, **kwargs):
        if "_id" in kwargs:
            self.doc = kwargs
        else:
            f = dict([(k, 1) for k in self.fields])
            self.doc = self.collection.find_one(kwargs, f)

    def __getitem__(self, k):
        return self._lazy_get_one(k)

    def __getattr__(self, k):
        try:
            return self._lazy_get_one(k)
        except KeyError:
            raise AttributeError("'%s' object has not attribute '%s'"
                    % (type(self), k))

    def _lazy_get_one(self, k):
        try:
            return self.doc[k]
        except KeyError:
            d = self.collection.find_one({"_id": self._id}, {"_id":0, "k": 1})
            self.doc[k] = d[k]
            return d[k]

    def _lazy_get(self, f):
        result = {}
        get_f  = {"_id": 0}
        for k in f:
            try:
                result[k] = self.doc[k]
            except KeyError:
                get_f[k]  = 1
        if len(get_f) > 1:
            d = self.collection.find_one({"_id": self._id}, get_f)
            for k in get_f:
                self.doc[k] = d[k]
                result[k] = d[k]
        return result

    @classmethod
    def find(cls, q):
        f = dict([(k, 1) for k in cls.fields])
        c = cls.collection.find(q, f)
        return [cls(**d) for d in c]

    @classmethod
    def collection(cls):
        return _db[cls.cname]

class SDSSModel(Model):
    @property
    def collection(self):
        return _sdss_db[self.cname]

class Star(SDSSModel):
    cname  = "stars"
    fields = ["ra", "dec"] + "u g r i z".split()

class Field(SDSSModel):
    cname  = "fields"
    fields = [""]

