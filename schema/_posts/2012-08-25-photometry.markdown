---
    layout: default
    title: Calibrated Photometry Model
---

Members
-------

* `id` (integer primary key): The id of this calibrated measurement.
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
* `flux` (real): The calibrated photometry of the source.
* `fluxerr` (real): The uncertainty in `flux`.
