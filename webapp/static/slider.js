(function () {

    var root = this;

    root.Slider = function (target, min, max, on_finish, on_change) {
        var _this = this;
        var w = 260, h = 20, margin = 10, radius = 5;

        // On change event.
        this.on_change = function () {};
        if (typeof(on_change) !== "undefined") this.on_change = on_change;

        // And the on finish function.
        this.on_finish = function () {};
        if (typeof(on_finish) !== "undefined") this.on_finish = on_finish;

        // Set up the scale.
        this.min = min;
        this.max = max;
        this.current_min = min;
        this.current_max = max;
        this.scale = d3.scale.linear().range([radius,
                                              w - radius])
                                      .domain([min, max]);

        // Select the target element.
        this.el = d3.select(target)
                    .append("svg:svg")
                        .attr("width", w + 2 * margin)
                        .attr("height", h)
                    .append("svg:g")
                        .attr("transform", "translate(" + margin + ", "
                                                        + 0 + ")");

        // Draw the bar.
        this.bar = this.el.append("line")
                        .attr("y1", 0.5 * (h - 1))
                        .attr("y2", 0.5 * (h - 1))
                        .attr("x1", this.scale(min) - radius)
                        .attr("x2", this.scale(max) + radius)
                        .attr("stroke-width", 1.5)
                        .attr("stroke", "#111");

        // Set up the buttons.
        this.selected = undefined;
        this.buttons = this.el.selectAll("circle.slider-button")
                        .data([{order: -1, pos: min}, {order: 1, pos: max}])
                        .enter()
                            .append("circle")
                                .attr("class", "slider-button")
                                .attr("cx", function (d) {
                                    return _this.scale(d.pos) + d.order * radius;
                                })
                                .attr("cy", 0.5 * (h - 1))
                                .attr("r", 5)
                                .attr("fill", "#d88")
                                .attr("stroke", "#111")
                                .attr("stroke-width", 1.5)
                                .attr("opacity", 1)
                                .style("cursor", "pointer")
                                .attr("order", function (d) { return d.order; })
                                .on("mousedown", function () {
                                    _this.selected = d3.select(this);
                                });

        // Body-wide event listeners.
        d3.select("body").on("mousemove", function () {
            if (typeof(_this.selected) !== "undefined") {
                // Find the value implied by the mouse position.
                var value = _this.scale.invert(d3.mouse(_this.el[0][0])[0]);

                // Bound the motion of the buttons.
                var x, order = _this.selected.attr("order");
                if (order == 1) {
                    x = Math.min(_this.max, Math.max(_this.current_min, value));
                    _this.current_max = x;
                } else {
                    x = Math.min(_this.current_max, Math.max(_this.min, value));
                    _this.current_min = x;
                }

                // Move the button sprite.
                _this.selected.attr("cx", _this.scale(x)
                                    + _this.selected.attr("order") * radius);

                // Fire the change event.
                _this.on_change([_this.current_min, _this.current_max]);
            }
        }).on("mouseup", function () {
            // Clear the selection.
            _this.selected = undefined;

            _this.on_finish([_this.current_min, _this.current_max]);
        });
    };

})();
