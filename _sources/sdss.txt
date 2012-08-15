.. _sdss:

SDSS Run Interface
==================

.. module:: calibration.data


I decided to combine all of the SDSS fields for one ``(run, camcol)`` pair
into a single object and store it as an HDF5 file. This has the advantage
that you never have to worry about edge effects. This is also really the
natural way to store the data because that's how it was taken! I've wrapped
these files with an object that handles all of the data access details.
Here's how it works:

.. autoclass:: SDSSRun
   :inherited-members:
