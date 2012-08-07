#!/usr/bin/env python

__all__ = ["init_schema"]


import psycopg2


def init_schema():
    connection = psycopg2.connect("dbname='sdss' host='localhost'")
    cursor = connection.cursor()

    # Create the stars table.
    cursor.execute("""CREATE TABLE IF NOT EXISTS stars
            (id SERIAL,
             ra REAL, dec REAL,
             u REAL, g REAL, r REAL, i REAL, z REAL,
             has_lightcurve BOOLEAN)
            """)

    # Create the runs table.
    cursor.execute("""CREATE TABLE IF NOT EXISTS runs
            (id SERIAL,
             run INTEGER, camcol INTEGER,
             field_min INTEGER, field_max INTEGER,
             band INTEGER,
             ramin REAL, ramax REAL, decmin REAL, decmax REAL)
            """)

    # Create the raw photometry table.
    cursor.execute("""CREATE TABLE IF NOT EXISTS raw
            (id SERIAL,
             runid INTEGER, starid INTEGER,
             ra REAL, dec REAL, tai REAL,
             band INTEGER,
             out_of_bounds BOOLEAN,
             flux REAL, fluxerr REAL,
             sky REAL, skyerr REAL,
             dx REAL, dxerr REAL,
             dy REAL, dyerr REAL)
            """)

    # Create the calibration run table.
    cursor.execute("""CREATE TABLE IF NOT EXISTS calibruns
            (id SERIAL,
             start_data TIMESTAMP,
             band INTEGER)
            """)

    # Create the calibration patch table.
    cursor.execute("""CREATE TABLE IF NOT EXISTS patches
            (id SERIAL,
             calibid INTEGER,
             ramin REAL, ramax REAL, decmin REAL, decmax REAL,
             stars INTEGER[], runs INTEGER[])
            """)

    # Create the calibrated photometry table.
    cursor.execute("""CREATE TABLE IF NOT EXISTS photometry
            (id SERIAL,
             rawid INTEGER, calibid INTEGER, patchid INTEGER,
             runid INTEGER, starid INTEGER,
             ra REAL, dec REAL, tai REAL,
             band INTEGER,
             flux REAL, fluxerr REAL)
            """)

    # Create the zero points table.
    cursor.execute("""CREATE TABLE IF NOT EXISTS zeros
            (id SERIAL,
             calibid INTEGER, patchid INTEGER, runid INTEGER,
             ramin REAL, ramax REAL, decmin REAL, decmax REAL,
             band INTEGER,
             zero REAL, beta2 REAL, delta2 REAL)
            """)

    # Create the mean fluxes table.
    cursor.execute("""CREATE TABLE IF NOT EXISTS fluxes
            (id SERIAL,
             calibid INTEGER, patchid INTEGER, starid INTEGER,
             ra REAL, dec REAL,
             band INTEGER,
             mean_flux REAL, eta2 REAL)
            """)

    connection.commit()
    connection.close()


if __name__ == "__main__":
    init_schema()