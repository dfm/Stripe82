from __future__ import print_function

import os
import sys
import time
from functools import wraps

import flask

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
    columns = ["id", "ra", "dec"] + [b for b in "ugriz"]
    q = ",".join(["stars.{0}".format(c) for c in columns])

    # Get the average eta2 too.
    columns += ["eta2"]
    q += ",(SELECT AVG(fluxes.eta2) FROM fluxes " \
            + "WHERE fluxes.starid = stars.id)"

    with flask.g.db as c:
        c.execute("SELECT {0} FROM stars ORDER BY random() LIMIT %s"
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


if __name__ == "__main__":
    app.run()
