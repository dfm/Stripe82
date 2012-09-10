import os
import sys
# import json

import flask

import config

sys.path.append(os.path.abspath(os.path.join(__file__, '..', '..')))
import calibration
import calibration.db
import calibration.models

# import lightcurve


app = flask.Flask(__name__)
app.config.from_object(config)


@app.before_request
def before_request():
    flask.g.db = calibration.db.DBConnection()


# Basic routes.

@app.route("/")
def index():
    return flask.redirect(flask.url_for(".stars"))
    # return flask.render_template("index.html")


# @app.route("/patches")
# def patches():
#     docs = flask.g.db.patches.find({}, {"ramin": 1, "ramax": 1,
#         "decmin": 1, "decmax": 1})

#     if docs is None:
#         flask.abort(404)

#     return flask.render_template("patches.html", patches=docs)


@app.route("/stars")
@app.route("/stars/<int:page>")
def stars(page=0):
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

    return flask.render_template("stars.html", stars=docs, prev_page=prev_page,
            next_page=next_page)


# @app.route("/lightcurve/<int:sid>")
# def patch(sid):
#     return flask.render_template("lightcurve.html", star=sid)


# API spec.

# def get_lightcurve(sid):
#     star = calibration.models.Star.find_one({"_id": sid})
#     if star is None:
#         flask.abort(404)

#     t, f, ferr = star.get_lightcurve(sid)
#     t /= 86400.0

#     return t, f, ferr


# @app.route("/api/lightcurve/<int:sid>")
# def api_lightcurve(sid):
#     t, f, ferr = get_lightcurve(sid)
#     doc = [{"t": t[i], "f": f[i], "ferr": ferr[i]} for i in range(len(t))]
#     return json.dumps(doc, default=str)


# @app.route("/api/period/<int:sid>")
# def api_period(sid):
#     t, f, ferr = get_lightcurve(sid)

#     period = lightcurve.find_period({"g": t}, {"g": f}, {"g": ferr}, order=5)
#     t = t % period

#     doc = [{"t": t[i], "f": f[i], "ferr": ferr[i]} for i in range(len(t))]
#     return json.dumps({"period": period, "lc": doc}, default=str)


if __name__ == "__main__":
    app.run()
