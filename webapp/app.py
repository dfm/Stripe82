import os
import sys
import json
import cPickle as pickle

import flask
from bson import ObjectId

import config

sys.path.append(os.path.abspath(os.path.join(__file__, '..', '..')))
import calibration
import calibration.db
import calibration.models


app = flask.Flask(__name__)
app.config.from_object(config)


@app.before_request
def before_request():
    flask.g.db = calibration.db.Database()


# Basic routes.

@app.route("/")
def index():
    return flask.render_template("index.html")


@app.route("/patches")
def patches():
    docs = flask.g.db.patches.find({}, {"ramin": 1, "ramax": 1,
        "decmin": 1, "decmax": 1})

    if docs is None:
        flask.abort(404)

    return flask.render_template("patches.html", patches=docs)


@app.route("/stars")
@app.route("/stars/<int:page>")
def stars(page=0):
    if page < 0:
        flask.abort(404)

    delta = 50
    docs = flask.g.db.stars.find({}, {"ra": 1, "dec": 1}) \
                           .sort("ra") \
                           .skip(page * delta) \
                           .limit(delta)
    if docs is None:
        flask.abort(404)

    next_page = page + 1
    prev_page = page - 1
    if prev_page < 0:
        prev_page = None

    return flask.render_template("stars.html", stars=docs, prev_page=prev_page,
            next_page=next_page)


@app.route("/lightcurve/<int:sid>")
def patch(sid):
    return flask.render_template("lightcurve.html", star=sid)


# API spec.

@app.route("/api/lightcurve/<int:sid>")
def api_lightcurve(sid):
    patch = calibration.models.Star.find_one({"_id": sid})
    if patch is None:
        flask.abort(404)

    t, f, ferr = patch.get_lightcurve(sid)

    doc = [{"t": t[i], "f": f[i], "ferr": ferr[i]} for i in range(len(t))]
    return json.dumps(doc, default=str)


if __name__ == "__main__":
    app.run()
