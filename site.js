/* ==========================================================================
   site.js — the small live parts of dayvalpraskihakim.github.io
   - Draws the model illustration at the top of each project page
     (elements with a data-plot attribute) as SVG, from fixed seeds.
   - Runs the homepage map: pins show the work linked to each place.
   - Shows the local time in Surabaya and the date the page was last updated.
   Colours come from style.css, so the charts follow light and dark mode.
   ========================================================================== */
(function () {
  "use strict";

  var NS = "http://www.w3.org/2000/svg";
  var uidCount = 0;

  /* ---------- random numbers ---------- */
  function rng(seed) { // mulberry32
    var a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) | 0;
      var t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function gauss(r) {
    var u = 0, v = 0;
    while (u === 0) u = r();
    while (v === 0) v = r();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }
  function studentT(r, nu) {
    var chi = 0;
    for (var i = 0; i < nu; i++) { var z = gauss(r); chi += z * z; }
    return gauss(r) / Math.sqrt(chi / nu);
  }
  function quantile(arr, q) {
    var s = arr.slice().sort(function (a, b) { return a - b; });
    var i = (s.length - 1) * q, lo = Math.floor(i), hi = Math.ceil(i);
    return s[lo] + (s[hi] - s[lo]) * (i - lo);
  }

  /* ---------- svg helpers ---------- */
  function el(name, attrs, parent) {
    var e = document.createElementNS(NS, name);
    for (var k in attrs) if (Object.prototype.hasOwnProperty.call(attrs, k)) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function txt(parent, x, y, s, anchor, cls) {
    var t = el("text", { x: x.toFixed(1), y: y.toFixed(1), "text-anchor": anchor || "start", "class": "p-label" + (cls ? " " + cls : "") }, parent);
    t.textContent = s;
    return t;
  }
  function line(parent, x1, y1, x2, y2, cls) {
    return el("line", { x1: x1.toFixed(1), y1: y1.toFixed(1), x2: x2.toFixed(1), y2: y2.toFixed(1), "class": cls }, parent);
  }
  function pathD(pts) {
    var s = "";
    for (var i = 0; i < pts.length; i++) s += (i ? "L" : "M") + pts[i][0].toFixed(1) + " " + pts[i][1].toFixed(1);
    return s;
  }
  function areaD(pts, yBase) {
    if (!pts.length) return "";
    return "M" + pts[0][0].toFixed(1) + " " + yBase.toFixed(1) + pathD(pts).replace(/^M/, "L") +
      "L" + pts[pts.length - 1][0].toFixed(1) + " " + yBase.toFixed(1) + "Z";
  }
  function scale(d0, d1, r0, r1) {
    var k = (r1 - r0) / ((d1 - d0) || 1);
    return function (v) { return r0 + (v - d0) * k; };
  }
  function frame(box) {
    while (box.firstChild) box.removeChild(box.firstChild);
    var W = Math.max(60, Math.round(box.clientWidth));
    var H = Math.max(28, Math.round(box.clientHeight));
    var svg = el("svg", { width: W, height: H, viewBox: "0 0 " + W + " " + H, "aria-hidden": "true", focusable: "false" }, box);
    return { svg: svg, W: W, H: H };
  }
  function clip(svg, x, y, w, h) {
    var id = "clip" + (++uidCount);
    var defs = el("defs", {}, svg);
    var cp = el("clipPath", { id: id }, defs);
    var rect = el("rect", { x: x, y: y, width: w, height: h }, cp);
    return { id: id, rect: rect };
  }

  /* ---------- Hawkes process ---------- */
  // multivariate version: A[i][j] is how much an event in lane j lifts lane i
  function simMultiHawkes(r, mu, A, beta, T) {
    var m = mu.length, ev = [], t = 0;
    var S = mu.map(function () { return 0; }), last = 0;
    while (t < T) {
      var d0 = Math.exp(-beta * (t - last)), lamBar = 0, i;
      for (i = 0; i < m; i++) lamBar += mu[i] + S[i] * d0;
      t += -Math.log(1 - r()) / lamBar;
      if (t >= T) break;
      var d = Math.exp(-beta * (t - last)), lam = [], tot = 0;
      for (i = 0; i < m; i++) { S[i] *= d; lam.push(mu[i] + S[i]); tot += lam[i]; }
      last = t;
      var u = r() * lamBar;
      if (u <= tot) {
        var acc = 0, j = 0;
        for (j = 0; j < m; j++) { acc += lam[j]; if (u <= acc) break; }
        j = Math.min(j, m - 1);
        ev.push({ t: t, j: j });
        for (i = 0; i < m; i++) S[i] += A[i][j];
      }
    }
    return ev;
  }
  function laneCurve(ev, mu, A, beta, T, lane, n) {
    var out = [], S = 0, last = 0, k = 0;
    for (var i = 0; i <= n; i++) {
      var t = (i / n) * T;
      while (k < ev.length && ev[k].t <= t) {
        var e = ev[k];
        S *= Math.exp(-beta * (e.t - last)); last = e.t;
        out.push([e.t, mu[lane] + S]);
        S += A[lane][e.j];
        out.push([e.t, mu[lane] + S]);
        k++;
      }
      S *= Math.exp(-beta * (t - last)); last = t;
      out.push([t, mu[lane] + S]);
    }
    return out;
  }
  function peakOf(curve) {
    var best = curve[0];
    for (var i = 1; i < curve.length; i++) if (curve[i][1] > best[1]) best = curve[i];
    return best;
  }

  /* ---------- Model illustrations (fixed seeds) ---------- */
  var DRAW = {};

  // Three currencies with cross-excitation: a shock in FX 1 spills into FX 2 and FX 3
  DRAW.hawkes3 = function (box, seed, mini) {
    var r = rng(seed), T = 252, beta = 0.5;
    var mu = [0.05, 0.03, 0.03];
    var A = [[0.3, 0.02, 0.02], [0.17, 0.22, 0.03], [0.17, 0.04, 0.22]];
    var ev = simMultiHawkes(r, mu, A, beta, T);
    var curves = [0, 1, 2].map(function (i) { return laneCurve(ev, mu, A, beta, T, i, 520); });
    var top = 0;
    curves.forEach(function (cv) { top = Math.max(top, peakOf(cv)[1]); });
    var f = frame(box), svg = f.svg, W = f.W, H = f.H;
    var padL = mini ? 0 : 46, padT = mini ? 2 : 20, padB = mini ? 2 : 24;
    var laneH = (H - padT - padB) / 3;
    var x = scale(0, T, padL, W);
    var shockT = peakOf(curves[0])[0];
    if (!mini) {
      line(svg, x(shockT), padT - 6, x(shockT), H - padB, "p-shock-dash");
      var rightSide = x(shockT) < W - 160;
      txt(svg, x(shockT) + (rightSide ? 8 : -8), padT - 8, "shock in FX 1, then contagion", rightSide ? "start" : "end", "shock");
    }
    curves.forEach(function (cv, i) {
      var t0 = padT + i * laneH, base = t0 + laneH - 5;
      var y = scale(0, top * 1.05, base, t0 + 3);
      var pts = cv.map(function (p) { return [x(p[0]), y(p[1])]; });
      el("path", { d: areaD(pts, base), "class": "p-area" }, svg);
      el("path", { d: pathD(pts), "class": "p-line" }, svg);
      line(svg, padL, base, W, base, "p-axis");
      ev.forEach(function (e) { if (e.j === i) line(svg, x(e.t), base - 6, x(e.t), base, "p-shock"); });
      if (!mini) txt(svg, 0, base - 4, "FX " + (i + 1), "start", "strong");
    });
    if (!mini) {
      txt(svg, padL, H - 6, "day 0", "start");
      txt(svg, W, H - 6, "day 252", "end");
    }
  };

  // Ornstein–Uhlenbeck paths with a shift in the mean (the US-yield shock)
  DRAW.ou = function (box, seed, mini) {
    var r = rng(seed), T = 60, dt = 0.25, n = Math.round(T / dt), K = mini ? 7 : 28;
    var theta = 0.12, sigma = 0.3, tShock = 18, jump = 1.1;
    var paths = [];
    for (var k = 0; k < K; k++) {
      var v = 0, pts = [[0, 0]];
      for (var i = 1; i <= n; i++) {
        var t = i * dt, m = t >= tShock ? jump : 0;
        v += theta * (m - v) * dt + sigma * Math.sqrt(dt) * gauss(r);
        pts.push([t, v]);
      }
      paths.push(pts);
    }
    var med = [], lo = [], hi = [];
    for (var j = 0; j <= n; j++) {
      var col = paths.map(function (p) { return p[j][1]; });
      med.push([j * dt, quantile(col, 0.5)]);
      lo.push([j * dt, quantile(col, 0.1)]);
      hi.push([j * dt, quantile(col, 0.9)]);
    }
    var all = [];
    paths.forEach(function (p) { p.forEach(function (q) { all.push(q[1]); }); });
    var yMin = quantile(all, 0.01), yMax = quantile(all, 0.995);
    var f = frame(box), svg = f.svg, W = f.W, H = f.H;
    var padT = mini ? 3 : 20, padB = mini ? 3 : 24;
    var x = scale(0, T, mini ? 1 : 0, W - (mini ? 1 : 0)), y = scale(yMin, yMax, H - padB, padT);
    var c = clip(svg, 0, padT - 2, W, H - padT - padB + 4);
    var g = el("g", { "clip-path": "url(#" + c.id + ")" }, svg);
    if (!mini) {
      var band = hi.map(function (p) { return [x(p[0]), y(p[1])]; })
        .concat(lo.slice().reverse().map(function (p) { return [x(p[0]), y(p[1])]; }));
      el("path", { d: pathD(band) + "Z", "class": "p-area" }, g);
    }
    paths.forEach(function (p) { el("path", { d: pathD(p.map(function (q) { return [x(q[0]), y(q[1])]; })), "class": "p-faint" }, g); });
    var meanPts = [[x(0), y(0)], [x(tShock), y(0)], [x(tShock), y(jump)], [x(T), y(jump)]];
    el("path", { d: pathD(meanPts), "class": "p-base" }, g);
    el("path", { d: pathD(med.map(function (q) { return [x(q[0]), y(q[1])]; })), "class": "p-line" }, g);
    if (!mini) {
      line(svg, x(tShock), padT - 8, x(tShock), H - padB, "p-shock-dash");
      txt(svg, x(tShock) + 8, padT - 8, W < 640 ? "US yield shock" : "US 10-year yield shock enters the drift", "start", "shock");
      txt(svg, W, y(jump) - 8, "new mean", "end");
      txt(svg, 0, H - 6, "day 0", "start");
      txt(svg, W, H - 6, "day 60", "end");
    }
  };

  // GARCH(1,1) returns with one large shock and slowly decaying volatility
  DRAW.garch = function (box, seed, mini) {
    var r = rng(seed), n = mini ? 64 : 128, omega = 0.06, a = 0.11, b = 0.83;
    var tS = Math.round(n * 0.45), zS = -4.4;
    var s2 = omega / (1 - a - b), eps = [], sig = [];
    for (var t = 0; t < n; t++) {
      var z = t === tS ? zS : gauss(r), s = Math.sqrt(s2), e = s * z;
      eps.push(e); sig.push(s);
      s2 = omega + a * e * e + b * s2;
    }
    var m = 0;
    for (var i = 0; i < n; i++) m = Math.max(m, Math.abs(eps[i]), 2 * sig[i]);
    var f = frame(box), svg = f.svg, W = f.W, H = f.H;
    var padT = mini ? 2 : 20, padB = mini ? 2 : 24;
    var step = (W - 2) / n;
    var x = function (k) { return 1 + (k + 0.5) * step; };
    var y = scale(-m * 1.05, m * 1.05, H - padB, padT), y0 = y(0);
    var up = [], dn = [];
    for (var k = 0; k < n; k++) { up.push([x(k), y(2 * sig[k])]); dn.push([x(k), y(-2 * sig[k])]); }
    el("path", { d: pathD(up.concat(dn.reverse())) + "Z", "class": "p-area" }, svg);
    line(svg, 0, y0, W, y0, "p-axis");
    for (k = 0; k < n; k++) line(svg, x(k), y0, x(k), y(eps[k]), k === tS ? "p-shock" : "p-needle");
    el("path", { d: pathD(up), "class": "p-line-2" }, svg);
    if (!mini) {
      txt(svg, x(tS) + 8, y(eps[tS]) + 4, "uncertainty shock", "start", "shock");
      txt(svg, x(Math.min(n - 1, tS + 16)), y(2 * sig[Math.min(n - 1, tS + 16)]) - 8, "±2σ, α + β = 0.94", "start");
      txt(svg, 0, H - 6, "month 1", "start");
      txt(svg, W, H - 6, "month 128", "end");
    }
  };

  // Merton / distance to default: asset paths, a debt barrier, and the horizon density
  DRAW.merton = function (box, seed, mini) {
    var r = rng(seed), V0 = 100, mu = 0.06, sg = 0.25, D = 62, T = 1, n = 120, K = mini ? 9 : 30;
    var dt = T / n, paths = [];
    for (var k = 0; k < K; k++) {
      var v = V0, pts = [[0, v]];
      for (var i = 1; i <= n; i++) { v *= Math.exp((mu - sg * sg / 2) * dt + sg * Math.sqrt(dt) * gauss(r)); pts.push([i * dt, v]); }
      paths.push(pts);
    }
    var DD = (Math.log(V0 / D) + (mu - sg * sg / 2) * T) / (sg * Math.sqrt(T));
    var EV = V0 * Math.exp(mu * T);
    var f = frame(box), svg = f.svg, W = f.W, H = f.H;
    var padT = mini ? 3 : 20, padB = mini ? 3 : 24;
    var split = mini ? W : Math.round(W * 0.74);
    var vmax = 0;
    paths.forEach(function (p) { p.forEach(function (q) { vmax = Math.max(vmax, q[1]); }); });
    var yLo = mini ? D - 6 : 40, yHi = mini ? vmax + 4 : 178;
    var x = scale(0, T, mini ? 1 : 0, split - (mini ? 1 : 0)), y = scale(yLo, yHi, H - padB, padT);
    var c = clip(svg, 0, padT - 2, W, H - padT - padB + 4);
    var g = el("g", { "clip-path": "url(#" + c.id + ")" }, svg);
    paths.forEach(function (p) { el("path", { d: pathD(p.map(function (q) { return [x(q[0]), y(q[1])]; })), "class": "p-faint" }, g); });
    var exp = [];
    for (var s = 0; s <= 40; s++) exp.push([x(s / 40 * T), y(V0 * Math.exp(mu * s / 40 * T))]);
    el("path", { d: pathD(exp), "class": "p-line" }, g);
    line(svg, 0, y(D), W, y(D), "p-shock-dash");
    if (mini) return;
    // horizon density of V_T, drawn sideways
    var m = Math.log(V0) + (mu - sg * sg / 2) * T, sd = sg * Math.sqrt(T);
    var dens = function (v) { return Math.exp(-Math.pow(Math.log(v) - m, 2) / (2 * sd * sd)) / (v * sd * Math.sqrt(2 * Math.PI)); };
    var x0 = split + 18, maxW = W - x0 - 4, fmax = 0, vs = [];
    for (var vv = 41; vv <= 177; vv += 1) { vs.push(vv); fmax = Math.max(fmax, dens(vv)); }
    var curve = vs.map(function (v) { return [x0 + dens(v) / fmax * maxW, y(v)]; });
    var full = [[x0, y(vs[0])]].concat(curve).concat([[x0, y(vs[vs.length - 1])]]);
    el("path", { d: pathD(full) + "Z", "class": "p-area" }, svg);
    var tail = vs.filter(function (v) { return v <= D; }).map(function (v) { return [x0 + dens(v) / fmax * maxW, y(v)]; });
    if (tail.length) {
      var tp = [[x0, y(vs[0])]].concat(tail).concat([[x0, y(D)]]);
      el("path", { d: pathD(tp) + "Z", "class": "p-wash" }, svg);
    }
    el("path", { d: pathD(curve), "class": "p-line-2" }, svg);
    line(svg, x0, padT, x0, H - padB, "p-axis");
    // distance-to-default bracket
    var bx = split + 6;
    line(svg, bx, y(EV), bx, y(D), "p-line");
    line(svg, bx - 4, y(EV), bx + 4, y(EV), "p-line");
    line(svg, bx - 4, y(D), bx + 4, y(D), "p-line");
    txt(svg, bx - 9, (y(EV) + y(D)) / 2 + 4, "DD ≈ " + DD.toFixed(1) + "σ", "end", "strong");
    txt(svg, 0, y(D) - 7, "default barrier: debt D", "start", "shock");
    txt(svg, 0, padT - 8, "firm asset value V", "start", "strong");
    txt(svg, W, padT - 8, "V in one year", "end");
    txt(svg, 0, H - 6, "today", "start");
    txt(svg, split, H - 6, "1 year", "end");
  };

  // Peaks over threshold for daily returns, above a storm-exposure lane
  DRAW.pot = function (box, seed, mini) {
    var r = rng(seed), n = mini ? 90 : 260, rets = [];
    for (var t = 0; t < n; t++) rets.push(0.9 * studentT(r, 4));
    var u = quantile(rets, 0.05);
    var storms = [];
    for (var k = 0; k < 4; k++) storms.push({ c: n * (0.55 + 0.33 * r()), w: n * (0.012 + 0.012 * r()), h: 0.45 + 0.55 * r() });
    var expo = [];
    for (t = 0; t < n; t++) {
      var e = 0;
      storms.forEach(function (s) { e += s.h * Math.exp(-Math.pow((t - s.c) / s.w, 2)); });
      expo.push(e);
    }
    var f = frame(box), svg = f.svg, W = f.W, H = f.H;
    var padT = mini ? 2 : 20, padB = mini ? 2 : 24;
    var laneH = mini ? 0 : Math.round((H - padT - padB) * 0.24);
    var mainB = H - padB - laneH - (mini ? 0 : 10);
    var mR = Math.max(Math.abs(quantile(rets, 0.005)), Math.abs(quantile(rets, 0.995))) * 1.08;
    var step = (W - 2) / n;
    var x = function (i) { return 1 + (i + 0.5) * step; };
    var y = scale(-mR * 1.05, mR * 1.05, mainB, padT), y0 = y(0);
    line(svg, 0, y0, W, y0, "p-axis");
    var cl = function (v) { return Math.max(-mR, Math.min(mR, v)); };
    for (var i = 0; i < n; i++) line(svg, x(i), y0, x(i), y(cl(rets[i])), rets[i] < u ? "p-shock" : "p-needle");
    line(svg, 0, y(u), W, y(u), "p-shock-dash");
    for (i = 0; i < n; i++) if (rets[i] < u) el("circle", { cx: x(i).toFixed(1), cy: y(cl(rets[i])).toFixed(1), r: mini ? 1.8 : 2.4, "class": "p-dot-shock" }, svg);
    if (mini) return;
    var lb = H - padB, lt = lb - laneH;
    var ey = scale(0, 1.05, lb, lt + 2);
    var ep = expo.map(function (v, k) { return [x(k), ey(v)]; });
    el("path", { d: areaD(ep, lb), "class": "p-area" }, svg);
    el("path", { d: pathD(ep), "class": "p-line-2" }, svg);
    line(svg, 0, lb, W, lb, "p-axis");
    txt(svg, 0, lt + 10, "storm exposure", "start", "strong");
    txt(svg, W, padT - 8, "threshold u: worst 5% of days", "end", "shock");
    txt(svg, 0, padT - 8, "daily return", "start", "strong");
    txt(svg, 0, H - 6, "day 1", "start");
    txt(svg, W, H - 6, "day " + n, "end");
  };

  // Gravity: log trade against log distance, divergent-regime pairs sitting lower
  DRAW.gravity = function (box, seed, mini) {
    var r = rng(seed), N = mini ? 34 : 150, pts = [];
    for (var i = 0; i < N; i++) {
      var xv = r(), div = r() < 0.3;
      pts.push({ x: xv, y: 0.9 - 0.62 * xv + 0.085 * gauss(r) - (div ? 0.12 : 0), d: div });
    }
    var same = pts.filter(function (p) { return !p.d; });
    var mx = 0, my = 0;
    same.forEach(function (p) { mx += p.x; my += p.y; });
    mx /= same.length; my /= same.length;
    var sxy = 0, sxx = 0;
    same.forEach(function (p) { sxy += (p.x - mx) * (p.y - my); sxx += (p.x - mx) * (p.x - mx); });
    var slope = sxy / sxx, icpt = my - slope * mx, shift = 0, nd = 0;
    pts.forEach(function (p) { if (p.d) { shift += p.y - (icpt + slope * p.x); nd++; } });
    shift = nd ? shift / nd : 0;
    var ys = pts.map(function (p) { return p.y; });
    var f = frame(box), svg = f.svg, W = f.W, H = f.H;
    var padT = mini ? 4 : 22, padB = mini ? 4 : 24, padL = mini ? 3 : 4;
    var x = scale(-0.02, 1.02, padL, W - 3), y = scale(Math.min.apply(null, ys) - 0.04, Math.max.apply(null, ys) + 0.04, H - padB, padT);
    line(svg, x(-0.02), y(icpt + slope * -0.02), x(1.02), y(icpt + slope * 1.02), "p-line");
    if (!mini) line(svg, x(-0.02), y(icpt + shift + slope * -0.02), x(1.02), y(icpt + shift + slope * 1.02), "p-shock-dash");
    pts.forEach(function (p) {
      if (p.d) el("circle", { cx: x(p.x).toFixed(1), cy: y(p.y).toFixed(1), r: mini ? 2 : 3, "class": "p-ring-shock" }, svg);
      else el("circle", { cx: x(p.x).toFixed(1), cy: y(p.y).toFixed(1), r: mini ? 1.8 : 2.6, "class": "p-dot-soft" }, svg);
    });
    if (mini) return;
    line(svg, 0, H - padB, W, H - padB, "p-axis");
    txt(svg, 0, padT - 9, W < 520 ? "log trade" : "log bilateral trade", "start", "strong");
    txt(svg, W, H - 6, "log distance →", "end");
    var lx = W - 218;
    el("circle", { cx: lx, cy: padT - 13, r: 3, "class": "p-dot-soft" }, svg);
    txt(svg, lx + 9, padT - 9, "same regime", "start");
    el("circle", { cx: lx + 100, cy: padT - 13, r: 3, "class": "p-ring-shock" }, svg);
    txt(svg, lx + 109, padT - 9, "divergent regime", "start", "shock");
  };

  /* ---------- mounting ---------- */
  function render(box) {
    var fn = DRAW[box.getAttribute("data-plot")];
    if (!fn) return;
    var seed = parseInt(box.getAttribute("data-seed") || "1", 10);
    fn(box, seed, box.hasAttribute("data-mini"));
  }

  function mountAll() {
    var boxes = Array.prototype.slice.call(document.querySelectorAll("[data-plot]"));
    var ro = "ResizeObserver" in window ? new ResizeObserver(function (entries) {
      entries.forEach(function (en) {
        var b = en.target, w = Math.round(en.contentRect.width), h = Math.round(en.contentRect.height);
        if (b._w === w && b._h === h) return;
        b._w = w; b._h = h;
        render(b);
      });
    }) : null;
    boxes.forEach(function (b) {
      b._w = Math.round(b.clientWidth); b._h = Math.round(b.clientHeight);
      try { render(b); } catch (err) { if (window.console) console.error(err); }
      if (ro) ro.observe(b);
    });
  }

  /* ---------- local time and last update ---------- */
  function clocks() {
    var els = document.querySelectorAll("[data-clock]");
    if (!els.length) return;
    var fmt;
    try { fmt = new Intl.DateTimeFormat("en-GB", { hour: "2-digit", minute: "2-digit", hour12: false, timeZone: "Asia/Jakarta" }); }
    catch (e) { return; }
    function tick() { var s = fmt.format(new Date()); for (var i = 0; i < els.length; i++) els[i].textContent = s; }
    tick();
    setInterval(tick, 30000);
  }
  function updated() {
    var els = document.querySelectorAll("[data-updated]");
    var d = new Date(document.lastModified);
    if (!els.length || isNaN(d.getTime())) return;
    var s = d.toLocaleDateString("en-GB", { month: "long", year: "numeric" });
    for (var i = 0; i < els.length; i++) els[i].textContent = s;
  }

  /* ---------- homepage map ---------- */
  function atlas() {
    var map = document.querySelector("[data-atlas]");
    var src = document.getElementById("atlas-data");
    var card = document.getElementById("atlas-card");
    if (!map || !src || !card) return;
    var data;
    try { data = JSON.parse(src.textContent); } catch (e) { return; }
    var pins = Array.prototype.slice.call(map.querySelectorAll(".pin"));
    var byKey = {};
    pins.forEach(function (p) { byKey[p.getAttribute("data-place")] = p; });
    var active = null, locked = false, timer = 0;

    function el(tag, cls, text) {
      var e = document.createElement(tag);
      if (cls) e.className = cls;
      if (text) e.textContent = text;
      return e;
    }
    function fill(key) {
      var pl = data.places[key];
      while (card.firstChild) card.removeChild(card.firstChild);
      var head = el("div", "card-head");
      head.appendChild(el("p", "card-place", pl[0]));
      head.appendChild(el("p", "card-note", pl[1]));
      var ul = el("ul", "card-work");
      pl[2].forEach(function (file) {
        var w = data.work[file], li = el("li"), a = el("a", "", w[0]);
        a.href = file;
        if (w[1]) a.appendChild(el("span", "k", w[1]));
        li.appendChild(a);
        ul.appendChild(li);
      });
      card.appendChild(head);
      card.appendChild(ul);
    }
    function show(pin) {
      if (!pin) return;
      clearTimeout(timer);
      if (active === pin) return;
      if (active) { active.classList.remove("on"); active.setAttribute("aria-expanded", "false"); }
      active = pin;
      pin.classList.add("on");
      pin.setAttribute("aria-expanded", "true");
      fill(pin.getAttribute("data-place"));
    }
    function rest() { locked = false; show(byKey[data.start]); }

    pins.forEach(function (pin) {
      pin.addEventListener("mouseenter", function () { if (!locked) show(pin); });
      pin.addEventListener("focus", function () { show(pin); });
      pin.addEventListener("click", function (e) {
        e.stopPropagation();
        locked = !(active === pin && locked);
        if (locked) show(pin); else rest();
      });
    });
    map.addEventListener("mouseleave", function () { if (!locked) timer = setTimeout(rest, 600); });
    card.addEventListener("mouseenter", function () { clearTimeout(timer); });
    document.addEventListener("click", function (e) {
      if (locked && !map.contains(e.target) && !card.contains(e.target)) rest();
    });
    document.addEventListener("keydown", function (e) { if (e.key === "Escape" && locked) rest(); });

    rest();
    // on phones the map scrolls sideways: start with home in view
    var scroller = map.parentNode;
    if (scroller.scrollWidth > scroller.clientWidth + 4) {
      var hx = parseFloat(byKey[data.start].style.left) / 100 * map.clientWidth;
      scroller.scrollLeft = Math.max(0, hx - scroller.clientWidth * 0.6);
    }
  }

  function init() { mountAll(); atlas(); clocks(); updated(); }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
