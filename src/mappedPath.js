//------------------------------------------------------------------------------
//
// Eschersketch - A drawing program for exploring symmetrical designs
//
//
// Copyright (c) 2026 Anselm Levskaya (http://anselmlevskaya.com)
// Licensed under the MIT (http://www.opensource.org/licenses/mit-license.php)
// license.
//
//------------------------------------------------------------------------------

// Drawing a shape through one symmetry transform
//------------------------------------------------------------------------------
// Tools describe shapes in drawing coordinates (before any symmetry is applied) and
// mappedPath(ctx, af) draws their image under the transform af.
//
// Affine transforms map lines to lines and ellipses to ellipses, so mapping the
// control points is exact and AffinePath issues the same canvas calls as before.
//
// Mobius transforms (hyperbolic symmetries) map lines and circles to circles, but
// do not preserve centers, Bezier curves or ellipses. MobiusPath draws each line
// segment as the exact arc through its mapped endpoints and midpoint, each circle
// or circular arc exactly, and Bezier curves and true ellipses by mapping the
// control points of pieces that are split in half until they are accurate to
// BEZIER_TOLERANCE_PX.

import {l2dist, pointToAngle} from './math_utils';

const {abs, sin, cos, tan, atan2, ceil, max, sqrt, PI} = Math;
const TAU = 2 * PI;

// mapped bulge (px) below which an arc is drawn as a straight segment
const ARC_SAGITTA_PX = 0.1;
// allowed error (px) at the midpoint of each mapped Bezier piece, and the maximum number of
// times a curve is halved (2^6 = 64 pieces)
const BEZIER_TOLERANCE_PX = 0.1;
const MAX_BEZIER_DEPTH = 6;

const lerp = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t];
// (Math.hypot is much slower than sqrt in V8, and these run for every copy)
const norm = (x, y) => sqrt(x * x + y * y);
const dist = (a, b) => norm(b[0] - a[0], b[1] - a[1]);

// point on the cubic with control points P, evaluated at parameters (u, v, w) of
// its blossom; the piece of the curve over [t0, t1] has control points
// blossom(t0,t0,t0), blossom(t0,t0,t1), blossom(t0,t1,t1), blossom(t1,t1,t1)
const blossom = function(P, u, v, w) {
  const a0 = lerp(P[0], P[1], u), a1 = lerp(P[1], P[2], u), a2 = lerp(P[2], P[3], u);
  return lerp(lerp(a0, a1, v), lerp(a1, a2, v), w);
};

// point at parameter t on the ellipse ctx.ellipse(c, rx, ry, rot) would draw,
// and its derivative with respect to t
const ellipsePoint = (c, rx, ry, rot, t) =>
  [c[0] + rx * cos(t) * cos(rot) - ry * sin(t) * sin(rot),
   c[1] + rx * cos(t) * sin(rot) + ry * sin(t) * cos(rot)];
const ellipseTangent = (rx, ry, rot, t) =>
  [-rx * sin(t) * cos(rot) - ry * cos(t) * sin(rot),
   -rx * sin(t) * sin(rot) + ry * cos(t) * cos(rot)];

// circle through three points, or null if they are (nearly) collinear
const circumcircle = function(a, b, c) {
  const bx = b[0] - a[0], by = b[1] - a[1], cx = c[0] - a[0], cy = c[1] - a[1];
  const d = 2 * (bx * cy - by * cx);
  if (d === 0) { return null; }
  const b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
  const ux = (cy * b2 - by * c2) / d, uy = (bx * c2 - cx * b2) / d;
  return {x: a[0] + ux, y: a[1] + uy, r: norm(ux, uy)};
};

// continues the current path (which ends at a) along the circle through a, b, c,
// from a via b to c
const arcThrough = function(ctx, a, b, c) {
  const chord = dist(a, c);
  // distance of b from the chord a-c
  const sagitta = chord > 0 ?
    abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / chord : 0;
  const circ = sagitta < ARC_SAGITTA_PX ? null : circumcircle(a, b, c);
  if (!circ) {
    ctx.lineTo(c[0], c[1]);
    return;
  }
  const angA = atan2(a[1] - circ.y, a[0] - circ.x);
  const angB = atan2(b[1] - circ.y, b[0] - circ.x);
  const angC = atan2(c[1] - circ.y, c[0] - circ.x);
  const sweep = ang => (((ang - angA) % TAU) + TAU) % TAU; // clockwise angle from a
  ctx.arc(circ.x, circ.y, circ.r, angA, angC, sweep(angB) > sweep(angC));
};

class AffinePath {
  constructor(ctx, af) {
    this.ctx = ctx;
    this.af = af;
  }

  moveTo(p) {
    const Tp = this.af.onVec(p);
    this.ctx.moveTo(Tp[0], Tp[1]);
  }

  lineTo(p) {
    const Tp = this.af.onVec(p);
    this.ctx.lineTo(Tp[0], Tp[1]);
  }

  bezierCurveTo(c0, c1, p) {
    const Tc0 = this.af.onVec(c0), Tc1 = this.af.onVec(c1), Tp = this.af.onVec(p);
    this.ctx.bezierCurveTo(Tc0[0], Tc0[1], Tc1[0], Tc1[1], Tp[0], Tp[1]);
  }

  // the ellipse arc ctx.ellipse(center, |p1 - center|, |p2 - center|, angle of p1, 0, arcAngle)
  // would draw; starts a new subpath
  ellipse(center, p1, p2, arcAngle) {
    const Tp0 = this.af.onVec(center), Tp1 = this.af.onVec(p1), Tp2 = this.af.onVec(p2);
    this.ctx.ellipse(Tp0[0], Tp0[1], l2dist(Tp0, Tp1), l2dist(Tp0, Tp2), pointToAngle(Tp0, Tp1), 0, arcAngle, 0);
  }

  closePath() {
    this.ctx.closePath();
  }
}

class MobiusPath {
  constructor(ctx, af) {
    this.ctx = ctx;
    this.af = af;
    this.last = null;   // current point, in drawing coordinates
    this.Tlast = null;  // and mapped
    this.start = null;  // start of the current subpath, in drawing coordinates
  }

  moveTo(p) {
    this.last = this.start = p;
    this.Tlast = this.af.onVec(p);
    this.ctx.moveTo(this.Tlast[0], this.Tlast[1]);
  }

  lineTo(p) {
    const Tp = this.af.onVec(p);
    if (this.af.imageIsStraight && this.af.imageIsStraight(this.Tlast, Tp, ARC_SAGITTA_PX)) {
      this.ctx.lineTo(Tp[0], Tp[1]); // short enough that the arc can't bulge visibly
    } else {
      arcThrough(this.ctx, this.Tlast, this.af.onVec(lerp(this.last, p, 0.5)), Tp);
    }
    this.last = p;
    this.Tlast = Tp;
  }

  bezierCurveTo(c0, c1, p) {
    const P = [this.last, c0, c1, p];
    const T = q => this.af.onVec(q);
    // draws the piece of the curve over [t0, t1] (whose mapped endpoints are Tq0, Tq3) as the
    // Bezier through its mapped control points, splitting it in half until that Bezier is within
    // BEZIER_TOLERANCE_PX of the true mapped curve at its quarter points and midpoint (checking
    // only the midpoint misses S-shaped errors, which cancel there)
    const piece = (t0, t1, Tq0, Tq3, depth) => {
      const Tq1 = T(blossom(P, t0, t0, t1)), Tq2 = T(blossom(P, t0, t1, t1));
      const tm = (t0 + t1) / 2, Tm = T(blossom(P, tm, tm, tm));
      const Q = [Tq0, Tq1, Tq2, Tq3];
      const withinTolerance = [0.25, 0.5, 0.75].every(s => {
        const t = t0 + s * (t1 - t0);
        const approx = blossom(Q, s, s, s);
        const exact = s === 0.5 ? Tm : T(blossom(P, t, t, t));
        return norm(approx[0] - exact[0], approx[1] - exact[1]) < BEZIER_TOLERANCE_PX;
      });
      if (depth >= MAX_BEZIER_DEPTH || withinTolerance) {
        this.ctx.bezierCurveTo(Tq1[0], Tq1[1], Tq2[0], Tq2[1], Tq3[0], Tq3[1]);
      } else {
        piece(t0, tm, Tq0, Tm, depth + 1);
        piece(tm, t1, Tm, Tq3, depth + 1);
      }
    };
    const Tp = T(p);
    piece(0, 1, this.Tlast, Tp, 0);
    this.last = p;
    this.Tlast = Tp;
  }

  ellipse(center, p1, p2, arcAngle) {
    const rx = dist(center, p1), ry = dist(center, p2);
    const rot = atan2(p1[1] - center[1], p1[0] - center[0]);
    const at = t => ellipsePoint(center, rx, ry, rot, t);
    const full = arcAngle >= TAU - 1e-9;
    if (abs(rx - ry) <= 1e-9 * max(rx, ry)) {
      // a circle maps to a circle, which three mapped points determine exactly
      if (full) {
        const circ = circumcircle(this.af.onVec(at(0)), this.af.onVec(at(TAU / 3)), this.af.onVec(at(2 * TAU / 3)));
        if (circ) {
          this.last = this.start = at(0);
          this.Tlast = this.af.onVec(this.last);
          this.ctx.moveTo(circ.x + circ.r, circ.y);
          this.ctx.arc(circ.x, circ.y, circ.r, 0, TAU);
          return;
        }
      } else {
        this.moveTo(at(0));
        arcThrough(this.ctx, this.Tlast, this.af.onVec(at(arcAngle / 2)), this.af.onVec(at(arcAngle)));
        this.last = at(arcAngle);
        this.Tlast = this.af.onVec(this.last);
        return;
      }
    }
    // true ellipse: approximate by cubic Beziers spanning at most a quarter turn each
    const pieces = max(1, ceil(arcAngle / (PI / 2)));
    const dt = arcAngle / pieces;
    const k = 4 / 3 * tan(dt / 4); // standard cubic approximation of an arc
    this.moveTo(at(0));
    for (let i = 0; i < pieces; i++) {
      const t0 = i * dt, t1 = (i + 1) * dt;
      const p0 = at(t0), p3 = at(t1);
      const d0 = ellipseTangent(rx, ry, rot, t0), d1 = ellipseTangent(rx, ry, rot, t1);
      this.bezierCurveTo([p0[0] + k * d0[0], p0[1] + k * d0[1]],
                         [p3[0] - k * d1[0], p3[1] - k * d1[1]], p3);
    }
  }

  closePath() {
    if (this.start && this.last && dist(this.start, this.last) > 0) {
      this.lineTo(this.start);
    }
    this.ctx.closePath();
  }
}

export const mappedPath = (ctx, af) =>
  af.preservesLines === false ? new MobiusPath(ctx, af) : new AffinePath(ctx, af);
