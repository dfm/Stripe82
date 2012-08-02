---
    layout: default
    title: Calibrated Photometry Model
---

Members
-------

* `photoid` (integer primary key): The id of this calibrated measurement.
* `rawid` (integer): The id of the associated [raw measurement](/models/raw).
* `calibid` (integer): The id of the associated [calibration
  run](/models/calibruns).
* `patchid` (integer): The id of the associated [calibration
  patch](/models/patches).
* `runid` (integer): The associated [run](/models/runs).
* `starid` (integer): The associated [star](/models/stars).
* `ra` (real): The R.A. position of the measurement (based on the associated
  [star](/models/stars) model).
* `dec` (real): The Dec. position of the measurement.
* `tai` (real): The time of the measurement (based on the interpolation of
  times associated with the [run](/models/runs) model). The units are seconds.
* `band` (integer): The SDSS filter. See [bands](/bands.html) for definitions.
* `flux` (real): The measured counts of the source.
* `ferr` (real): The uncertainty in `flux`.
* `eta2` (real): The relative variability parameter found for the star.
* `beta2` (real): The relative variability parameter found for the run.
