#!/usr/bin/env python

"""
Tools for one off migration from MongoDB to PostgreSQL

"""


import pymongo
import psycopg2


def to_bandid(band):
    return "ugriz".index(band)


def get_dbs():
    mongodb = pymongo.Connection().sdss
    connection = psycopg2.connect("dbname='sdss'")
    return mongodb, connection


def migrate_stars():
    db, connection = get_dbs()
    cursor = connection.cursor()
    c = db.stars.find()
    for star in c:
        cursor.execute("""INSERT INTO stars
                (ra, dec, u, g, r, i, z, has_lightcurve)
                VALUES
                (%s, %s, %s, %s, %s, %s, %s, FALSE)""",
                [star["ra"], star["dec"], star["u"], star["g"], star["r"],
                star["i"], star["z"]])
    connection.commit()
    connection.close()


def migrate_runs():
    db, connection = get_dbs()
    cursor = connection.cursor()
    c = db.runs.find()
    for run in c:
        ramin, ramax = run["raMin"], run["raMax"]
        decmin, decmax = run["decMin"], run["decMax"]

        print ramin, ramax,

        # Fix the wrapping problems with R.A.
        # If the min is less than 180 and max is greater than 180 then they
        # have been swapped. Let's fix that.
        if ramin < 180.0 and ramax > 180.0:
            ramin, ramax = ramax, ramin

        # Now, make sure that everything is in the range ``[-180, 180]``.
        while ramin > 180.0:
            ramin -= 360.0
        while ramax > 180.0:
            ramax -= 360.0

        print ramin, ramax

        cursor.execute("""INSERT INTO runs
                (run, camcol, field_min, field_max, band, ramin, ramax,
                    decmin, decmax)
                VALUES
                (%s, %s, %s, %s, %s, %s, %s, %s, %s)""",
                [run["run"], run["camcol"], min(run["fields"]),
                    max(run["fields"]), to_bandid(run["band"]),
                    ramin, ramax, decmin, decmax])
    connection.commit()
    connection.close()


if __name__ == "__main__":
    print "Migrating stars table..."
    migrate_stars()

    print "Migrating runs table..."
    migrate_runs()
