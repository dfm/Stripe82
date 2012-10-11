var dfmplot = function () {
    // Defaults.
    var width = 800,
        height = 500,
        margin = 75,
        xticks = 4,
        yticks = 5,
        xlabel = "",
        ylabel = "",
        xlim = null,
        ylim = null,
        padding = 0.01;

    // Plot elements.
    var element,
        xscale = d3.scale.linear(),
        yscale = d3.scale.linear(),
        xaxis = d3.svg.axis().scale(xscale),
        yaxis = d3.svg.axis().scale(yscale).orient("left"),
        xaxisel, yaxisel,
        elements = [];

    var formatAxes = function (el) {
        el.selectAll("line")
            .attr("stroke-width", "1px")
            .attr("stroke", "#000")
            .attr("fill", "none");

        el.selectAll("path")
            .attr("stroke-width", "1px")
            .attr("stroke", "#000")
            .attr("fill", "none");
    };

    var plot = function (el) {
        // Scales.
        xscale.range([0, width]);
        yscale.range([height, 0]);

        // Base elements.
        element = d3.select(el).append("svg")
                        .attr("width", width + 2 * margin)
                        .attr("height", height + 2 * margin)
                    .append("g")
                        .attr("transform", "translate(" + margin + ", "
                                                        + margin + ")");
        element.append("rect")
                        .attr("class", "background")
                        .attr("x", 0)
                        .attr("y", 0)
                        .attr("width", width)
                        .attr("height", height)
                        .attr("fill", "#fff")
                        .attr("stroke", "none");

        // Axes elements.
        xaxisel = element.append("g").attr("class", "x axis")
                           .attr("transform", "translate(0, " + height + ")")
                           .call(yaxis);
        yaxisel = element.append("g").attr("class", "y axis")
                           .call(xaxis);
        formatAxes(xaxisel);
        formatAxes(yaxisel);

        // Axes labels.
        // MAGIC number 6 pixel offset for labels—should depend on font size.
        element.append("text").attr("class", "x label")
                              .attr("text-anchor", "end")
                              .attr("x", width)
                              .attr("y", height - 6)
                              .text(xlabel);
        element.append("text").attr("class", "y label")
                              .attr("text-anchor", "end")
                              .attr("y", 6)
                              .attr("dy", ".75em")
                              .attr("transform", "rotate(-90)")
                              .text(ylabel);

        return plot;
    };

    var _scale_axis = function (data, fmin, fmax, current) {
        if (typeof(fmax) === "undefined" || fmax === null)
            fmax = fmin;

        var mn = d3.min(data, fmin),
            mx = d3.max(data, fmax),
            delta;

        if (typeof(current) !== "undefined" && current !== null
                && current.length == 2) {
            mn = d3.min([current[0], mn]);
            mx = d3.min([current[1], mx]);
        }

        delta = padding * (mx - mn);

        return [mn - delta, mx + delta];
    };

    var _pop = function (d, k, v) {
        if (typeof(d) !== "undefined")
            if (k in d) return d[k];
        if (typeof(v) !== "undefined") return v;
        throw "Key Error: " + k + " not in object.";
    };

    plot.plot = function (data, cx, cy, args) {
        var cssclass = _pop(args, "cssclass", "datapoint"),
            fill = _pop(args, "fill", "#000"),
            radius = _pop(args, "radius", 2.0);

        if (typeof(cx) !== "function") {
            var _cx = cx;
            cx = function (d) { return d[_cx]; };
        }
        if (typeof(cy) !== "function") {
            var _cy = cy;
            cy = function (d) { return d[_cy]; };
        }

        // Auto-scale the axes.
        xlim = _scale_axis(data, cx, xlim);
        ylim = _scale_axis(data, cy, ylim);

        // Push the plotting function to the stack.
        plot.append(function () {
            var selection = element.selectAll("circle." + cssclass)
                            .data(data)
            selection.enter().append("circle")
                            .attr("class", cssclass)
                            .attr("cx", function (d) { return xscale(cx(d)); })
                            .attr("cy", function (d) { return yscale(cy(d)); })
                            .attr("r", radius);
            selection.exit().remove();
        });

        return plot;
    };

    plot.append = function (el) {
        elements.push(el);
        return plot;
    };

    plot.draw = function () {
        xscale.domain(xlim);
        yscale.domain(ylim);
        xaxis.ticks(xticks);
        yaxis.ticks(yticks);
        xaxisel.call(xaxis);
        yaxisel.call(yaxis);

        for (var k in elements) {
            elements[k](plot);
        }

        return plot;
    };

    // Getters and setters.
    plot.width = function (value) {
        if (!arguments.length) return width;
        width = value;
        return plot;
    };

    plot.height = function (value) {
        if (!arguments.length) return height;
        height = value;
        return plot;
    };

    plot.margin = function (value) {
        if (!arguments.length) return margin;
        margin = value;
        return plot;
    };

    plot.xticks = function (value) {
        if (!arguments.length) return xticks;
        xticks = value;
        return plot;
    };

    plot.yticks = function (value) {
        if (!arguments.length) return yticks;
        yticks = value;
        return plot;
    };

    plot.xlabel = function (value) {
        if (!arguments.length) return xlabel;
        xlabel = value;
        return plot;
    };

    plot.ylabel = function (value) {
        if (!arguments.length) return ylabel;
        ylabel = value;
        return plot;
    };

    plot.xlim = function (value) {
        if (!arguments.length) return xlim;
        xlim = value;
        return plot;
    };

    plot.ylim = function (value) {
        if (!arguments.length) return ylim;
        ylim = value;
        return plot;
    };

    plot.padding = function (value) {
        if (!arguments.length) return padding;
        padding = value;
        return plot;
    };

    // Getters only.
    plot.element = function () {
        return element;
    };

    plot.xscale = function (value) {
        if (!arguments.length) return xscale;
        return xscale(value);
    };

    plot.yscale = function (value) {
        if (!arguments.length) return yscale;
        return yscale(value);
    };

    return plot;
};
