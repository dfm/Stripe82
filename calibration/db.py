"""


"""

__all__ = ["Database"]

import os
import pymongo


class Database(object):
    def __init__(self, **kwargs):
        self.server = kwargs.pop("server",
                os.environ.get("MONGO_SERVER", "localhost"))
        self.port = int(kwargs.pop("port",
            os.environ.get("MONGO_PORT", 27017)))
        self.name = kwargs.pop("name", "sdss")
        self.db = pymongo.Connection(self.server, self.port)[self.name]

    def __getitem__(self, i):
        return self.db[i]

    def __getattr__(self, attr):
        return getattr(self.db, attr)
