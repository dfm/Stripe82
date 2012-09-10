"""


"""

__all__ = ["Database", "DBConnection"]

import os
import pymongo
import psycopg2


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


class DBConnection(object):
    def __init__(self, dbname="sdss"):
        self._connection = psycopg2.connect("dbname='{0}'".format(
                                os.environ.get("SDSS_DB_NAME", "sdss")))
        self._cursor = self._connection.cursor()

    def __enter__(self):
        return self._cursor

    def __exit__(self, exc_type, exc_value, traceback):
        if traceback is None:
            self._connection.commit()
        self._cursor.close()
        self._connection.close()
