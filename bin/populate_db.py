#!/usr/bin/env python
# encoding: utf-8
"""
Scrape CAS for all of the needed tables and populate the local MongoDB.

"""

import os

import pymongo

import casjobs

_local_tmp_dir   = os.path.join(os.environ.get("SDSS_LOCAL", "."), ".sdss")

# Connect to database.
_db_server = os.environ.get("MONGO_SERVER", "localhost")
_db_port   = int(os.environ.get("MONGO_PORT", 27017))
_db_name   = os.environ.get("MONGO_DB", "sdss")
_db = pymongo.Connection(_db_server, _db_port)[_db_name]

def populate_db():
    # Get the fields.
    field_table = "fields"
    jobs = casjobs.CasJobs()

    # How many fields are there?
    field_count = jobs.count("FROM Stripe82..Field")

    # MAGIC: An approximation of the maximum number of fields that can be
    # collected without filling MyDB.
    max_rows = 75000

    q = """SELECT * INTO mydb.{0} FROM
(SELECT ROW_NUMBER() OVER (ORDER BY fieldID ASC) AS ROW_NUMBER, *
    FROM Stripe82..Field) foo
WHERE ROW_NUMBER BETWEEN {1} AND {2}
""".format(field_table, 0, max_rows)
    try:
        jobs.drop_table(field_table)
    except:
        pass
    job_id = jobs.submit(q)
    status = jobs.monitor(job_id)
    if status[0] != 5:
        raise Exception("Couldn't complete field list request.")

    # Download the output file.
    outfn = os.path.join(_local_tmp_dir, "fields.fits")
    jobs.request_and_get_output("fields", "FITS", outfn)

if __name__ == '__main__':
    populate_db()

