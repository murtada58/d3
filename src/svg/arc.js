import "../core/functor";
import "../core/zero";
import "../math/trigonometry";
import "../geom/polygon";
import "svg";

d3.svg.arc = function() {
  var innerRadius = d3_svg_arcInnerRadius,
      outerRadius = d3_svg_arcOuterRadius,
      cornerRadius = d3_zero,
      padRadius = d3_svg_arcAuto,
      startAngle = d3_svg_arcStartAngle,
      endAngle = d3_svg_arcEndAngle,
      padAngle = d3_svg_arcPadAngle;

  function arc() {
    var r0 = Math.max(0, +innerRadius.apply(this, arguments)),
        r1 = Math.max(0, +outerRadius.apply(this, arguments)),
        a0 = startAngle.apply(this, arguments) - halfπ,
        a1 = endAngle.apply(this, arguments) - halfπ,
        da = Math.abs(a1 - a0),
        cw = a0 > a1 ? 0 : 1;

    // Ensure that the outer radius is always larger than the inner radius.
    if (r1 < r0) rc = r1, r1 = r0, r0 = rc;

    // Special case for an arc that spans the full circle.
    if (da >= τε) return circleSegment(r1, cw) + (r0 ? circleSegment(r0, 1 - cw) : "") + "Z";

    var rc,
        cr,
        rp,
        ap,
        p0 = 0,
        p1 = 0,
        x0,
        y0,
        x1,
        y1,
        x2,
        y2,
        x3,
        y3,
        path = [];

    // The recommended minimum inner radius when using padding is outerRadius *
    // padAngle / sin(θ), where θ is the angle of the smallest arc (without
    // padding). For example, if the outerRadius is 200 pixels and the padAngle
    // is 0.02 radians, a reasonable θ is 0.04 radians, and a reasonable
    // innerRadius is 100 pixels.

    if (ap = (+padAngle.apply(this, arguments) || 0) / 2) {
      rp = padRadius === d3_svg_arcAuto ? Math.sqrt(r0 * r0 + r1 * r1) : +padRadius.apply(this, arguments);
      if (!cw) p1 *= -1;
      if (r1) p1 = d3_asin(rp / r1 * Math.sin(ap));
      if (r0) p0 = d3_asin(rp / r0 * Math.sin(ap));
    }

    // Compute the two outer corners.
    if (r1) {
      x0 = r1 * Math.cos(a0 + p1);
      y0 = r1 * Math.sin(a0 + p1);
      x1 = r1 * Math.cos(a1 - p1);
      y1 = r1 * Math.sin(a1 - p1);

      // Detect whether the outer corners are collapsed.
      var l1 = Math.abs(a1 - a0 - 2 * p1) <= π ? 0 : 1;
      if (p1 && d3_svg_arcSweep(x0, y0, x1, y1) === cw ^ l1) {
        var h1 = (a0 + a1) / 2;
        x0 = r1 * Math.cos(h1);
        y0 = r1 * Math.sin(h1);
        x1 = y1 = null;
      }
    } else {
      x0 = y0 = 0;
    }

    // Compute the two inner corners.
    if (r0) {
      x2 = r0 * Math.cos(a1 - p0);
      y2 = r0 * Math.sin(a1 - p0);
      x3 = r0 * Math.cos(a0 + p0);
      y3 = r0 * Math.sin(a0 + p0);

      // Detect whether the inner corners are collapsed.
      var l0 = Math.abs(a0 - a1 + 2 * p0) <= π ? 0 : 1;
      if (p0 && d3_svg_arcSweep(x2, y2, x3, y3) === (1 - cw) ^ l0) {
        var h0 = (a0 + a1) / 2;
        x2 = r0 * Math.cos(h0);
        y2 = r0 * Math.sin(h0);
        x3 = y3 = null;
      }
    } else {
      x2 = y2 = 0;
    }

    // Compute the rounded corners.
    var cornerRadii = cornerRadius.apply(this, arguments);
    if (typeof(cornerRadii) === 'number') { cornerRadii = [cornerRadii, cornerRadii, cornerRadii, cornerRadii]; }
    else if (cornerRadii.length === 1) { cornerRadii = [cornerRadii[0], cornerRadii[0], cornerRadii[0], cornerRadii[0]]; }
    else if (cornerRadii.length === 2) { cornerRadii = [cornerRadii[0], cornerRadii[1], cornerRadii[1], cornerRadii[0]]; }
    else if (cornerRadii.length === 3) { cornerRadii = [cornerRadii[0], cornerRadii[1], cornerRadii[1], cornerRadii[2]]; }
    for (var i = 0; i < 4; i++) { cornerRadii[i] = Math.min(Math.abs(r1 - r0) / 2, cornerRadii[i]); }
    
    if (Math.max(cornerRadii[0], cornerRadii[1], cornerRadii[2], cornerRadii[3]) > 1e-3) {
      cr = r0 < r1 ^ cw ? 0 : 1;

      // Compute the angle of the sector formed by the two sides of the arc.
      var oc = x3 == null ? [x2, y2] : x1 == null ? [x0, y0] : d3_geom_polygonIntersect([x0, y0], [x3, y3], [x1, y1], [x2, y2]),
          ax = x0 - oc[0],
          ay = y0 - oc[1],
          bx = x1 - oc[0],
          by = y1 - oc[1],
          kc = 1 / Math.sin(Math.acos((ax * bx + ay * by) / (Math.sqrt(ax * ax + ay * ay) * Math.sqrt(bx * bx + by * by))) / 2),
          lc = Math.sqrt(oc[0] * oc[0] + oc[1] * oc[1]);

      // Compute the outer corners.
      if (x1 != null) {
        var crtl = Math.min(cornerRadii[0], (r1 - lc) / (kc + 1)),
            crtr = Math.min(cornerRadii[1], (r1 - lc) / (kc + 1)),
            t30 = d3_svg_arcCornerTangents(x3 == null ? [x2, y2] : [x3, y3], [x0, y0], r1, crtl, cw), // rc1 here controlls top left corner
            t12 = d3_svg_arcCornerTangents([x1, y1], [x2, y2], r1, crtr, cw); // rc1 here controlls top right corner

        path.push(
          "M", t30[0],
          "A", crtl, ",", crtl, " 0 0,", cr, " ", t30[1],
          "A", r1, ",", r1, " 0 ", (1 - cw) ^ d3_svg_arcSweep(t30[1][0], t30[1][1], t12[1][0], t12[1][1]), ",", cw, " ", t12[1],
          "A", crtr, ",", crtr, " 0 0,", cr, " ", t12[0]);

      } else {
        path.push("M", x0, ",", y0);
      }

      // Compute the inner corners.
      if (x3 != null) {
        var crbl = Math.min(cornerRadii[2], (r0 - lc) / (kc - 1)),
            crbr = Math.min(cornerRadii[3], (r0 - lc) / (kc - 1)),
            t03 = d3_svg_arcCornerTangents([x0, y0], [x3, y3], r0, -crbl, cw), // rc0 here controlls bottom left corner
            t21 = d3_svg_arcCornerTangents([x2, y2], x1 == null ? [x0, y0] : [x1, y1], r0, -crbr, cw); // rc0 here controlls bottom right corner

        path.push(
          "L", t21[0],
          "A", crbr, ",", crbr, " 0 0,", cr, " ", t21[1],
          "A", r0, ",", r0, " 0 ", cw ^ d3_svg_arcSweep(t21[1][0], t21[1][1], t03[1][0], t03[1][1]), ",", 1 - cw, " ", t03[1],
          "A", crbl, ",", crbl, " 0 0,", cr, " ", t03[0]);
      } else {
        path.push("L", x2, ",", y2);
      }
    }

    // Compute straight corners.
    else {
      path.push("M", x0, ",", y0);
      if (x1 != null) path.push("A", r1, ",", r1, " 0 ", l1, ",", cw, " ", x1, ",", y1);
      path.push("L", x2, ",", y2);
      if (x3 != null) path.push("A", r0, ",", r0, " 0 ", l0, ",", 1 - cw, " ", x3, ",", y3);
    }

    path.push("Z");
    return path.join("");
  }

  function circleSegment(r1, cw) {
    return "M0," + r1
        + "A" + r1 + "," + r1 + " 0 1," + cw + " 0," + -r1
        + "A" + r1 + "," + r1 + " 0 1," + cw + " 0," + r1;
  }

  arc.innerRadius = function(v) {
    if (!arguments.length) return innerRadius;
    innerRadius = d3_functor(v);
    return arc;
  };

  arc.outerRadius = function(v) {
    if (!arguments.length) return outerRadius;
    outerRadius = d3_functor(v);
    return arc;
  };

  arc.cornerRadius = function(v) {
    if (!arguments.length) return cornerRadius;
    cornerRadius = d3_functor(v);
    return arc;
  };

  arc.padRadius = function(v) {
    if (!arguments.length) return padRadius;
    padRadius = v == d3_svg_arcAuto ? d3_svg_arcAuto : d3_functor(v);
    return arc;
  };

  arc.startAngle = function(v) {
    if (!arguments.length) return startAngle;
    startAngle = d3_functor(v);
    return arc;
  };

  arc.endAngle = function(v) {
    if (!arguments.length) return endAngle;
    endAngle = d3_functor(v);
    return arc;
  };

  arc.padAngle = function(v) {
    if (!arguments.length) return padAngle;
    padAngle = d3_functor(v);
    return arc;
  };

  arc.centroid = function() {
    var r = (+innerRadius.apply(this, arguments) + +outerRadius.apply(this, arguments)) / 2,
        a = (+startAngle.apply(this, arguments) + +endAngle.apply(this, arguments)) / 2 - halfπ;
    return [Math.cos(a) * r, Math.sin(a) * r];
  };

  return arc;
};

var d3_svg_arcAuto = "auto";

function d3_svg_arcInnerRadius(d) {
  return d.innerRadius;
}

function d3_svg_arcOuterRadius(d) {
  return d.outerRadius;
}

function d3_svg_arcStartAngle(d) {
  return d.startAngle;
}

function d3_svg_arcEndAngle(d) {
  return d.endAngle;
}

function d3_svg_arcPadAngle(d) {
  return d && d.padAngle;
}

// Note: similar to d3_cross2d, d3_geom_polygonInside
function d3_svg_arcSweep(x0, y0, x1, y1) {
  return (x0 - x1) * y0 - (y0 - y1) * x0 > 0 ? 0 : 1;
}

// Compute perpendicular offset line of length rc.
// http://mathworld.wolfram.com/Circle-LineIntersection.html
function d3_svg_arcCornerTangents(p0, p1, r1, rc, cw) {
  var x01 = p0[0] - p1[0],
      y01 = p0[1] - p1[1],
      lo = (cw ? rc : -rc) / Math.sqrt(x01 * x01 + y01 * y01),
      ox = lo * y01,
      oy = -lo * x01,
      x1 = p0[0] + ox,
      y1 = p0[1] + oy,
      x2 = p1[0] + ox,
      y2 = p1[1] + oy,
      x3 = (x1 + x2) / 2,
      y3 = (y1 + y2) / 2,
      dx = x2 - x1,
      dy = y2 - y1,
      d2 = dx * dx + dy * dy,
      r = r1 - rc,
      D = x1 * y2 - x2 * y1,
      d = (dy < 0 ? -1 : 1) * Math.sqrt(r * r * d2 - D * D),
      cx0 = (D * dy - dx * d) / d2,
      cy0 = (-D * dx - dy * d) / d2,
      cx1 = (D * dy + dx * d) / d2,
      cy1 = (-D * dx + dy * d) / d2,
      dx0 = cx0 - x3,
      dy0 = cy0 - y3,
      dx1 = cx1 - x3,
      dy1 = cy1 - y3;

  // Pick the closer of the two intersection points.
  // TODO Is there a faster way to determine which intersection to use?
  if (dx0 * dx0 + dy0 * dy0 > dx1 * dx1 + dy1 * dy1) cx0 = cx1, cy0 = cy1;

  return [
    [cx0 - ox, cy0 - oy],
    [cx0 * r1 / r, cy0 * r1 / r]
  ];
}
