(function () {

    var root = this;

    var astro = {};

    // From:
    // trac.astrometry.net/browser/trunk/src/astrometry/util/starutil_numpy.py
    var ra_normalize = astro.ra_normalize = function (ra) {
        var result = ra;
        while (result < 0.0) result += 360;
        return result % 360.;
    };

    var ra2hms = astro.ra2hms = function (ra0) {
        var ra = ra_normalize(ra0);
        var h = ra * 24. / 360.;
        var hh = Math.floor(h);
        var m = (h - hh) * 60.;
        var mm = Math.floor(m);
        var s = (m - mm) * 60.
        return [hh, mm, s];
    };

    var dec2dms = astro.dec2dms = function (dec) {
        var sgn = (dec >= 0) ? 1 : -1;
        var d = dec * sgn;
        var dd = Math.floor(d);
        var m = (d - dd) * 60;
        var mm = Math.floor(m);
        var s = (m - mm) * 60;

        if (s >= 60.) {
            m += 1.
            s -= 60.
        }

        // Don't just return sgn * d because values between 0 and 1 deg
        // will get you!
        return [sgn, dd, mm, s];
    };

    var zero_pad = function (num, size) {
        var result = "" + num;
        var decimal = size + 1;

        while (decimal > size) {
            decimal = result.indexOf(".");
            if (decimal < 0) decimal = result.length;
            if (decimal < size) result = "0" + result;
        }
        return result;
    };

    astro.get_sdss_name = function (ra, dec) {
        var n = "SDSS J";
        var r = ra2hms(ra), d = dec2dms(dec);

        // R.A.
        n += zero_pad(r[0].toFixed(0), 2);
        n += zero_pad(r[1].toFixed(0), 2);
        n += zero_pad(r[2].toFixed(2), 2);

        // Dec.
        if (d[0] < 0) n += "-";
        else n += "+";
        n += zero_pad(d[1].toFixed(0), 2);
        n += zero_pad(d[2].toFixed(0), 2);
        n += zero_pad(d[3].toFixed(2), 2);

        return n;
    };

    root.Astro = astro;

})();
