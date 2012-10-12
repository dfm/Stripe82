from __future__ import print_function

import os
import sys
import time
from functools import wraps

import numpy as np
import flask
import george

import config

try:
    import calibration
    import calibration.db
    calibration.db
except ImportError:
    sys.path.append(os.path.abspath(os.path.dirname(
                                                os.path.dirname(__file__))))
    import calibration
    import calibration.db
    calibration.db

# import lightcurve

app = flask.Flask(__name__)
app.config.from_object(config)
app.debug = True


def dfmtime(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        s = time.time()
        result = f(*args, **kwargs)
        print("'{0}' took {1} seconds".format(f.__name__, time.time() - s))
        return result
    return wrapper


@app.before_request
def before_request():
    flask.g.db = calibration.db.DBConnection()


# Basic routes.

@app.route("/")
def index():
    return flask.redirect(flask.url_for(".stars"))


@app.route("/stars/")
@app.route("/stars/<int:starid>")
def stars(starid=None):
    if starid is None:
        return flask.render_template("stars.html")

    columns = ["id", "ra", "dec"] + [b for b in "ugriz"]

    with flask.g.db as c:
        c.execute("SELECT {0} FROM stars WHERE id = %s LIMIT 1"
                    .format(",".join(columns)), (starid,))
        star = c.fetchone()

    if star is None:
        flask.abort(404)

    star = dict(zip(columns, star))

    return flask.render_template("starinfo.html", star=star)


@app.route("/runs/<int:runid>")
def runs(runid=None):
    # columns = ["id"]

    # with flask.g.db as c:
    #     c.execute("SELECT {0} FROM runs WHERE id = %s LIMIT 1"
    #                 .format(",".join(columns)), (runid,))
    #     run = c.fetchone()

    # if run is None:
    #     flask.abort(404)

    # run = dict(zip(columns, run))

    return flask.render_template("run.html", run={"id": runid})


@app.route("/starlist/")
@app.route("/starlist/<int:page>")
def starlist(page=0):
    if page < 0:
        flask.abort(404)

    delta = 50
    columns = ["id", "ra", "dec", "g"]

    with flask.g.db as c:
        c.execute("SELECT {0} FROM stars ORDER BY ra LIMIT %s OFFSET %s"
                    .format(",".join(columns)),
                (delta, page * delta))
        docs = c.fetchall()

    if len(docs) is 0:
        flask.abort(404)

    docs = [dict(zip(columns, d)) for d in docs]

    next_page = page + 1
    prev_page = page - 1
    if prev_page < 0:
        prev_page = None

    return flask.render_template("starlist.html", stars=docs,
            prev_page=prev_page, next_page=next_page)


# API

@app.route("/api/sample/stars")
@app.route("/api/sample/stars/<int:number>")
@dfmtime
def sample_stars(number=1000):
    columns = ["id", "ra", "dec", "eta2"] + [b for b in "ugriz"]
    q = ",".join(columns)

    with flask.g.db as c:
        c.execute(("SELECT {0} FROM starview WHERE eta2 != 'NaN' ORDER BY "
                + "random() LIMIT %s")
                .format(q), (number,))
        docs = c.fetchall()

    if len(docs) is 0:
        flask.abort(404)

    docs = [dict(zip(columns, d)) for d in docs]

    return flask.jsonify(data=docs)


@app.route("/api/<method>/<int:starid>/<int:band>")
@dfmtime
def lightcurve(method=None, starid=None, band=None):
    if method not in ["raw", "calib"]:
        flask.abort(404)

    table = "raw"
    columns = ["tai", "flux", "fluxivar"]

    if method == "calib":
        table = "photometry"
        columns += ["calibid", "patchid"]

    with flask.g.db as c:
        c.execute("SELECT {0} FROM {1} WHERE starid = %s AND band = %s"
                .format(",".join(columns), table) + " ORDER BY tai",
                (starid, band))
        docs = c.fetchall()

    if len(docs) is 0:
        flask.abort(404)

    docs = [dict(zip(columns, d)) for d in docs]

    return flask.jsonify(data=docs)


@app.route("/api/run/<int:runid>")
@dfmtime
def run_zeropoint(runid=None):
    table = "zeros"
    columns = ["zero", "beta2", "delta2", "ramin", "ramax",
               "decmin", "decmax"]

    with flask.g.db as c:
        c.execute("SELECT {0} FROM {1} WHERE runid = %s"
                .format(",".join(columns), table),
                (runid, ))
        docs = c.fetchall()

    if len(docs) is 0:
        flask.abort(404)

    docs = [dict(zip(columns, d)) for d in docs]

    return flask.jsonify(data=docs)


@app.route("/api/run/<int:runid>/fit")
@dfmtime
def fit_run(runid=None):
    table = "zeros"
    columns = ["zero", "beta2", "delta2", "ramin", "ramax",
               "decmin", "decmax"]

    with flask.g.db as c:
        c.execute("SELECT {0} FROM {1} WHERE runid = %s"
                .format(",".join(columns), table),
                (runid, ))
        docs = c.fetchall()

    if len(docs) is 0:
        flask.abort(404)

    docs = [dict(zip(columns, d)) for d in docs]
    x = np.empty((len(docs), 2))
    y = np.empty(len(docs))
    for i, d in enumerate(docs):
        x[i, :] = [0.5 * (d["ramin"] + d["ramax"]),
                   0.5 * (d["decmin"] + d["decmax"])]
        y[i] = d["zero"]

    m, v = np.median(y), 1.0  # np.var(y)

    y -= m
    y /= v

    gp = george.GaussianProcess([0.1, 0.001, 0.005])

    # Iteratively reweighted fit.
    yerr = 0.01 * np.ones_like(y)
    gp.fit(x, y, yerr=yerr, normalize=False)

    for i in range(2):
        s = time.time()
        mu, var = gp.predict(x, full_cov=False)
        print(time.time() - s)
        var = var + yerr * yerr
        chi2 = (mu - y) ** 2 / var + np.log(2 * np.pi * var)
        yerr *= np.sqrt(1 + chi2 / 9.0)
        gp.fit(x, y, yerr=yerr, normalize=False)

    # Output prediction.
    nra = int(flask.request.args.get("nra", 25))
    ndec = int(flask.request.args.get("ndec", 25))
    X, Y = np.meshgrid(np.linspace(np.min(x[:, 0]), np.max(x[:, 0]), nra),
                       np.linspace(np.min(x[:, 1]), np.max(x[:, 1]), ndec))
    X = np.vstack([X.flatten(), Y.flatten()]).T
    mu, var = gp.predict(X, full_cov=False)

    dec = X[:, 1].reshape([ndec, nra])

    ra = X[:, 0].reshape([ndec, nra])
    mu = (v * mu + m).reshape([ndec, nra])
    std = np.sqrt(v * var).reshape([ndec, nra])
    cols = ["ra", "mu", "std"]

    result = [{"dec": dec[j, 0],
               "data": [dict(zip(cols, [ra[j, i], mu[j, i], std[j, i]]))
                                for i in range(len(ra[j]))]}
                                    for j in range(len(dec))]

    return flask.jsonify(data=result)


if __name__ == "__main__":
    app.run()
