// Ki-Engir planet model: an energy-balance climate on the canon relief, a moisture march
// toward the substellar point, Ki-Engir Climate Zones classes, and ecological succession.
// Pure JS, no DOM. Ki frame: lat N+, lon E+, SSP at (0,0). Row 0 = southernmost.
function createKiModel(D) {
  const NX = D.NX, NY = D.NY, N = NX * NY, RES = D.RES, h = RES * Math.PI / 180;
  const dec = (s, T) => { const bin = atob(s), u8 = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i); return new T(u8.buffer); };
  const land8 = dec(D.land, Uint8Array), zl16 = dec(D.zland, Int16Array), zm16 = dec(D.zmean, Int16Array), zr16 = dec(D.zrel, Int16Array);
  const g67 = dec(D.gcm067, Uint16Array), g70 = dec(D.gcm070, Uint16Array), gPr = dec(D.gcmPr, Uint16Array), gTas = dec(D.gcmTas, Uint16Array);
  const lf = new Float32Array(N), zl = new Float32Array(N), zmean = new Float32Array(N), zrel = new Float32Array(N);
  const gcm = { rc067: new Float32Array(N), rc070: new Float32Array(N), pr: new Float32Array(N), tas: new Float32Array(N) };
  for (let k = 0; k < N; k++) {
    lf[k] = land8[k] / 255; zl[k] = Math.max(0, zl16[k]); zmean[k] = zm16[k]; zrel[k] = zr16[k];
    gcm.rc067[k] = g67[k] / 100 + 150; gcm.rc070[k] = g70[k] / 100 + 150; gcm.pr[k] = gPr[k]; gcm.tas[k] = gTas[k] / 100 + 150;
  }
  const lonc = new Float32Array(NX), latc = new Float32Array(NY), cosl = new Float64Array(NY), cosN = new Float64Array(NY);
  for (let i = 0; i < NX; i++) lonc[i] = -180 + RES / 2 + RES * i;
  for (let j = 0; j < NY; j++) { latc[j] = -90 + RES / 2 + RES * j; cosl[j] = Math.cos(latc[j] * Math.PI / 180); cosN[j] = j === NY - 1 ? 0 : Math.cos((latc[j] + RES / 2) * Math.PI / 180); }
  const mu = new Float32Array(N), sunEl = new Float32Array(N), wgt = new Float32Array(N);
  for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
    const k = j * NX + i, m = Math.cos(latc[j] * Math.PI / 180) * Math.cos(lonc[i] * Math.PI / 180);
    mu[k] = m; sunEl[k] = Math.asin(Math.max(-1, Math.min(1, m))) * 180 / Math.PI; wgt[k] = cosl[j];
  }
  const idx = (lat, lon) => { const j = Math.max(0, Math.min(NY - 1, Math.floor((lat + 90) / RES))); const i = ((Math.floor((lon + 180) / RES) % NX) + NX) % NX; return j * NX + i; };

  // ---- distance to the sea (km), Dijkstra on the grid, R = 5,861 km
  const RKM = 5861, coastKm = new Float64Array(N).fill(1e9);
  {
    const heap = [], push = (d, k) => { heap.push([d, k]); let i = heap.length - 1; while (i > 0) { const p = (i - 1) >> 1; if (heap[p][0] <= heap[i][0]) break; [heap[p], heap[i]] = [heap[i], heap[p]]; i = p; } };
    const pop = () => { const top = heap[0], last = heap.pop(); if (heap.length) { heap[0] = last; let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i; if (l < heap.length && heap[l][0] < heap[m][0]) m = l; if (r < heap.length && heap[r][0] < heap[m][0]) m = r; if (m === i) break; [heap[m], heap[i]] = [heap[i], heap[m]]; i = m; } } return top; };
    for (let k = 0; k < N; k++) if (lf[k] < 0.5) { coastKm[k] = 0; push(0, k); }
    while (heap.length) {
      const [d, k] = pop(); if (d > coastKm[k]) continue;
      const j = (k / NX) | 0, i = k % NX, dx = RKM * h * cosl[j], dy = RKM * h;
      const nb = [[j, (i + 1) % NX, dx], [j, (i + NX - 1) % NX, dx]];
      if (j < NY - 1) nb.push([j + 1, i, dy]); if (j > 0) nb.push([j - 1, i, dy]);
      for (const [jj, ii, s] of nb) { const kk = jj * NX + ii, nd = d + s; if (nd < coastKm[kk]) { coastKm[kk] = nd; push(nd, kk); } }
    }
    // a land cell's own half-width counts: the coast sits inside a mixed cell
    for (let k = 0; k < N; k++) if (lf[k] >= 0.5) coastKm[k] = Math.max(0, coastKm[k] - RKM * h * 0.5);
  }

  // ---- polar (SSP-centred) frame for the moisture march: bearing b from Ki-north, distance d from SSP
  const NB = 180, NR = 90, dB = 360 / NB, dR = 180 / NR;
  const p2k = new Int32Array(NB * NR), k2p = new Int32Array(N), sinD = new Float64Array(NR);
  for (let r = 0; r < NR; r++) sinD[r] = Math.max(Math.sin((r + 0.5) * dR * Math.PI / 180), 1e-3);
  for (let b = 0; b < NB; b++) for (let r = 0; r < NR; r++) {
    const d = (r + 0.5) * dR * Math.PI / 180, be = (-180 + (b + 0.5) * dB) * Math.PI / 180;
    const pE = -(Math.PI / 2 - d); // Earth latitude
    const lat = Math.asin(Math.cos(pE) * Math.cos(be)) * 180 / Math.PI;
    const lon = Math.atan2(Math.cos(pE) * Math.sin(be), -Math.sin(pE)) * 180 / Math.PI;
    p2k[b * NR + r] = idx(lat, lon);
  }
  for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
    const k = j * NX + i, p = latc[j] * Math.PI / 180, l = lonc[i] * Math.PI / 180;
    const d = Math.acos(Math.max(-1, Math.min(1, Math.cos(p) * Math.cos(l)))) * 180 / Math.PI;
    const be = Math.atan2(Math.cos(p) * Math.sin(l), Math.sin(p)) * 180 / Math.PI;
    const r = Math.min(NR - 1, Math.floor(d / dR)), b = ((Math.floor((be + 180) / dB) % NB) + NB) % NB;
    k2p[k] = b * NR + r;
  }
  // relief along the march, for orographic lift: polar-cell surface height
  const pz = new Float32Array(NB * NR);
  for (let q = 0; q < NB * NR; q++) { const k = p2k[q]; pz[q] = lf[k] >= 0.5 ? zl[k] : 0; }

  const P = Object.assign({
    S0: 1600, co2: 1000, cloudBase: 0.10, cloudMax: 0.852, albCloud: 0.615, cloudLW: 0.30,
    // star and orbit: with useOrbit the flux and the rotation follow the period by Kepler's third law
    useOrbit: 1, starT: 2566, starL: 5.53e-4, starM: 0.0898, period: 5.0, spec: { fvis: 0.015, ratio: [1.046, 3.192, 4.331] },
    // spectrum: 1 = water vapour, cloud and ice respond to the star's colour; 0 = the old gray shortwave
    spectral: 1, cloudNir: 0.7, iceNir: 0.33, wvScale: 1.49, wvTop: 0.25, cloudSens: 0.000, rhCol: 0.8,
    Datm: 7.40, Docean: 1.33, lapse: 3.89, invRate: 25,
    hUnstable: 18, hStable: 5, swSurf: 0.65, iceK: 0.73,
    mld: 50, albIce: 0.45, albSea: 0.07, albBare: 0.28, albVeg: 0.10,
    tauDry: 0.242, tauWv: 0.363, evap: 1.0, recycle: 0.9, petScale: 0.42,
    tauRes: 15, rhCrit: 0.8, Hw: 2000, hBL: 500, fricDays: 1, dragCd: 1.5e-3, stabK: 30, inflowH: 1500,
    succYears: 35, disturb: 0.004,
  }, D.params || {});

  // ---- star and orbit ------------------------------------------------------------------------
  // Near-infrared water-vapour absorption that responds to the star, the same construction as
  // ExoPlaSim's nirabsorption option: the Lacis & Hansen (1974) solar curve mapped onto the
  // three-band k-distribution of Chou & Lee (1996) / Chou & Suarez (1999), bands re-weighted by
  // the star's flux relative to a 5,772 K blackbody.
  const XK = [0.0010, 0.0133, 0.0422, 0.1334, 0.4217, 1.334, 5.623, 31.62, 177.8, 1000.0];
  const HK = [[.20673, .03497, .03011, .02260, .01336, .00696, .00441, .00115, .00026, 0],
              [.08236, .01157, .01133, .01143, .01240, .01258, .01381, .00650, .00244, .00094],
              [.01074, .00360, .00411, .00421, .00389, .00326, .00499, .00465, .00245, .00145]];
  const EDGE = [0.70e-6, 1.22e-6, 2.27e-6, 10e-6], NTAB = 201, LY0 = -7, DLY = 0.05;
  const wvTab = new Float64Array(NTAB), star = { T: 0, fvis: 0.517, ratio: [1, 1, 1], key: '' };
  let OMEGA = 2 * Math.PI / (6.097 * 86400);
  function bandFractions(Tk) {
    const n = 4000, l0 = Math.log(316.036e-9), l1 = Math.log(1e-4), f = [0, 0, 0]; let tot = 0, vis = 0;
    const B = w => Math.pow(w * 1e6, -5) / Math.expm1(0.0143877735383 / (w * Tk));
    for (let i = 0; i < n; i++) {
      const wa = Math.exp(l0 + (l1 - l0) * i / n), wb = Math.exp(l0 + (l1 - l0) * (i + 1) / n), wm = 0.5 * (wa + wb);
      const seg = 0.5 * (B(wa) + B(wb)) * (wb - wa); tot += seg; if (wm < 0.75e-6) vis += seg;
      for (let b = 0; b < 3; b++) if (wm >= EDGE[b] && wm < EDGE[b + 1]) f[b] += seg;
    }
    return { f: f.map(v => v / tot), vis: vis / tot };
  }
  const SUN = bandFractions(5772);
  const lh = y => 2.9 * y / (Math.pow(1 + 141.5 * y, 0.635) + 5.925 * y);
  function refreshStar() {
    const key = [P.starT, P.starL, P.starM, P.period, P.useOrbit, P.spectral, P.spec ? P.spec.ratio.join() + P.spec.fvis : ''].join();
    if (key === star.key) return;
    star.key = key;
    if (P.useOrbit) {
      const a = Math.cbrt(P.starM * Math.pow(P.period / 365.25, 2)); // AU
      star.a = a; P.S0 = 1361 * P.starL / (a * a);
      OMEGA = 2 * Math.PI / (P.period * 86400);
    }
    const bf = bandFractions(P.spectral ? P.starT : 5772);
    star.T = P.starT; star.fvis = bf.vis; star.ratio = bf.f.map((v, b) => v / SUN.f[b]);
    // a measured spectrum's band weights, where the star has one, replace the blackbody's
    if (P.spectral && P.spec) { star.fvis = P.spec.fvis; star.ratio = P.spec.ratio.slice(); bf.vis = P.spec.fvis; }
    const chou = (u, r) => { let a = 0; for (let b = 0; b < 3; b++) for (let k = 0; k < 10; k++) a += r[b] * HK[b][k] * (1 - Math.exp(-Math.min(XK[k] * u, 80))); return a; };
    const one = [1, 1, 1], capA = 0.98 * (1 - bf.vis);
    for (let t = 0; t < NTAB; t++) {
      const y = Math.pow(10, LY0 + t * DLY), target = lh(y); let lo = -12, hi = 6;
      for (let it = 0; it < 50; it++) { const m = 0.5 * (lo + hi); if (chou(Math.pow(10, m), one) < target) lo = m; else hi = m; }
      wvTab[t] = Math.min(chou(Math.pow(10, 0.5 * (lo + hi)), star.ratio), capA);
    }
  }
  // fraction of the starlight absorbed by an effective water amount y (cm)
  function absWv(y) {
    if (!(y > 0)) return 0;
    const x = (Math.log10(y) - LY0) / DLY;
    if (x <= 0) return wvTab[0] * y / Math.pow(10, LY0);
    if (x >= NTAB - 1) return wvTab[NTAB - 1];
    const j = Math.floor(x), w = x - j; return (1 - w) * wvTab[j] + w * wvTab[j + 1];
  }
  // how much darker a surface that is bright in the visible and dull in the near-infrared looks under this star
  const colour = rel => (star.fvis + (1 - star.fvis) * rel) / (SUN.vis + (1 - SUN.vis) * rel);

  const T = new Float32Array(N).fill(250), Ts = new Float32Array(N), alb = new Float32Array(N), asr = new Float32Array(N);
  const swSfc = new Float32Array(N), colW = new Float32Array(N);
  const prec = new Float32Array(N), pet = new Float32Array(N), ai = new Float32Array(N), evapF = new Float32Array(N);
  const stage = new Float32Array(N), cap = new Float32Array(N), cls = new Uint8Array(N), iceF = new Float32Array(N);
  const Tair = new Float32Array(N), landE = new Float32Array(N), latent = new Float64Array(N);
  let years = 0, rng = 12345;
  const rand = () => { rng = (rng * 1664525 + 1013904223) >>> 0; return rng / 4294967296; };
  const SIG = 5.670e-8;
  const esat = t => 6.112 * Math.exp(17.67 * (t - 273.15) / (t - 29.65)); // hPa
  const ptf = t => { const tc = t - 273.15, es = esat(t), dl = 4098 * es / Math.pow(tc + 237.3, 2); return dl / (dl + 0.66); };
  const olr = (t, c) => { const tau = P.tauDry + 0.048 * Math.log(P.co2 / 280) + P.tauWv * Math.min(Math.exp(0.03 * (t - 288)), 2.5) + P.cloudLW * c; return SIG * t * t * t * t / (1 + 0.75 * Math.max(tau, 0)); };
  const rhoVsat = t => esat(t) * 100 / (461.5 * t);
  const smooth = (e0, e1, x) => { const t = Math.max(0, Math.min(1, (x - e0) / (e1 - e0))); return t * t * (3 - 2 * t); };

  // Downward longwave at the surface: Brutsaert clear-sky emissivity at 70% humidity, cloud filling the gap.
  function backRad(ta, c) {
    const e = 0.7 * esat(Math.max(ta, 180));
    const ec = Math.min(1, 1.24 * Math.pow(e / ta, 1 / 7));
    return (ec + (1 - ec) * 0.9 * c) * SIG * ta * ta * ta * ta;
  }
  // Skin temperature from the surface energy balance:
  // sigma Ts^4 + h (Ts - Ta) + LE + k_c (Ts - T_c) = down-LW + absorbed SW, with h unstable when the skin is warmer.
  function seb(ta, sw, le, c, kc, tc) {
    const down = backRad(ta, c) + sw - le;
    const f = t => SIG * t * t * t * t + (P.hStable + (P.hUnstable - P.hStable) * smooth(-1, 1, t - ta)) * (t - ta) + kc * (t - tc) - down;
    let t = ta;
    for (let i = 0; i < 6; i++) { const v = f(t), dv = (f(t + 0.05) - v) / 0.05; t -= v / dv; }
    return t;
  }
  function surfaceTemps() {
    for (let k = 0; k < N; k++) {
      const c = cloud(k), sw = P.spectral ? swSfc[k] : P.S0 * Math.max(mu[k], 0) * (1 - alb[k]) * P.swSurf;
      const ta = T[k] - P.lapse * zl[k] / 1000;
      Tair[k] = lf[k] >= 0.5 ? ta : T[k];
      const li = 1 - smooth(271.15, 275.15, ta), si = 1 - smooth(269.35, 273.35, T[k]);
      let tl = ta;
      if (lf[k] > 0) {
        tl = seb(ta, sw, landE[k] / 12.6, c, 0, 0);
        // slopes above the cold-pool floor: the canon 25 K/km, no warmer than the air
        if (tl < ta) tl += Math.min(P.invRate * zrel[k] / 1000, ta - tl);
        if (li > 0.5) tl = Math.min(tl, 273.15);
      }
      let tsea = T[k];
      if (si > 0) { const tice = Math.min(273.15, seb(T[k], sw, 0, c, P.iceK, 271.35)); tsea = si * tice + (1 - si) * T[k]; }
      Ts[k] = lf[k] >= 0.5 ? tl : tsea;
      iceF[k] = lf[k] * li + (1 - lf[k]) * si;
    }
  }
  // cloud cover: thickest under the star; with the spectral scheme it thins as the column warms
  // (the GCM's dayside cover falls about 0.1 per 12 K), which is the feedback behind the hot branch
  const cloud = k => {
    if (mu[k] <= 0) return P.cloudBase * 0.6;
    let c = P.cloudBase + (P.cloudMax - P.cloudBase) * Math.pow(mu[k], 0.8);
    if (P.spectral) c -= P.cloudSens * (T[k] - 270) / 10;
    return Math.max(0.05, Math.min(0.98, c));
  };
  function radiation() {
    for (let k = 0; k < N; k++) {
      const c = cloud(k);
      const ground = P.albBare + (P.albVeg - P.albBare) * (stage[k] / 5);
      const tl = T[k] - P.lapse * zl[k] / 1000;
      const li = 1 - smooth(271.15, 275.15, tl), si = 1 - smooth(269.35, 273.35, T[k]);
      const aI = P.spectral ? P.albIce * colour(P.iceNir) : P.albIce;
      const aL = li * aI + (1 - li) * ground, aS = si * aI * 0.85 + (1 - si) * P.albSea;
      const as_ = lf[k] * aL + (1 - lf[k]) * aS;
      if (!P.spectral) {
        alb[k] = c * P.albCloud + (1 - c) * as_;
        asr[k] = P.S0 * Math.max(mu[k], 0) * (1 - alb[k]); swSfc[k] = asr[k] * P.swSurf;
        continue;
      }
      // column water vapour (kg/m2) from the surface air temperature, and its slant path
      const w = P.rhCol * rhoVsat(Math.max(T[k], 180)) * P.Hw; colW[k] = w;
      const m0 = Math.max(mu[k], 0);
      if (m0 <= 0) { alb[k] = as_; asr[k] = 0; swSfc[k] = 0; continue; }
      const slant = 35 / Math.sqrt(1 + 1224 * m0 * m0), y = P.wvScale * slant * w / 10;
      const aTop = absWv(P.wvTop * y), aAll = absWv(y);
      // 1. vapour above the cloud deck absorbs; 2. cloud reflects what is left; 3. vapour below absorbs;
      // 4. the surface reflects part of what arrives, and that leaves through the gaps in the cloud
      const cA = c * P.albCloud * colour(P.cloudNir);
      const i1 = 1 - aTop, i2 = i1 * (1 - cA), i3 = i2 * (1 - (aAll - aTop) / Math.max(1 - aTop, 1e-6));
      const refl = i1 * cA + i3 * as_ * (1 - c);
      alb[k] = refl; asr[k] = P.S0 * m0 * (1 - refl); swSfc[k] = P.S0 * m0 * i3 * (1 - as_);
    }
  }

  // cyclic tridiagonal solve (a: sub, b: diag, c: super), Sherman-Morrison
  const ta = new Float64Array(NX), tb = new Float64Array(NX), tc = new Float64Array(NX), tr = new Float64Array(NX);
  const cp = new Float64Array(NX), dp = new Float64Array(NX), xs = new Float64Array(NX), zz = new Float64Array(NX), uu = new Float64Array(NX);
  function triSolve(a, b, c, d, x, n) {
    cp[0] = c[0] / b[0]; dp[0] = d[0] / b[0];
    for (let i = 1; i < n; i++) { const m = b[i] - a[i] * cp[i - 1]; cp[i] = c[i] / m; dp[i] = (d[i] - a[i] * dp[i - 1]) / m; }
    x[n - 1] = dp[n - 1]; for (let i = n - 2; i >= 0; i--) x[i] = dp[i] - cp[i] * x[i + 1];
  }
  const bb = new Float64Array(NX);
  function cyclicSolve(a, b, c, d, x) {
    const n = NX, alpha = c[n - 1], beta = a[0], gamma = -b[0];
    for (let i = 0; i < n; i++) bb[i] = b[i];
    bb[0] = b[0] - gamma; bb[n - 1] = b[n - 1] - alpha * beta / gamma;
    triSolve(a, bb, c, d, x, n);
    for (let i = 0; i < n; i++) uu[i] = 0; uu[0] = gamma; uu[n - 1] = alpha;
    triSolve(a, bb, c, uu, zz, n);
    const fact = (x[0] + beta * x[n - 1] / gamma) / (1 + zz[0] + beta * zz[n - 1] / gamma);
    for (let i = 0; i < n; i++) x[i] -= fact * zz[i];
  }

  const Cdt = new Float64Array(N), Bk = new Float64Array(N), Rk = new Float64Array(N);
  const Dx = new Float64Array(N), Dn = new Float64Array(N); // face D: east face of k, north face of k
  function energyStep(dtYears) {
    const dt = dtYears * 3.156e7;
    const Cs = 4.18e6 * P.mld, Cl = 2.0e7;
    for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
      const k = j * NX + i, ke = j * NX + (i + 1) % NX, kn = (j + 1) * NX + i;
      const oc = 1 - lf[k];
      Dx[k] = P.Datm + P.Docean * Math.min(oc, 1 - lf[ke]) * (1 - 0.7 * Math.max(iceF[k], iceF[ke]) * (oc > 0.5 ? 1 : 0));
      Dn[k] = j < NY - 1 ? P.Datm + P.Docean * Math.min(oc, 1 - lf[kn]) * (1 - 0.7 * Math.max(iceF[k], iceF[kn]) * (oc > 0.5 ? 1 : 0)) : 0;
      Cdt[k] = (lf[k] * Cl + oc * Cs) / dt;
      const c = cloud(k), o0 = olr(T[k], c), b = (olr(T[k] + 0.5, c) - olr(T[k] - 0.5, c));
      Bk[k] = b; Rk[k] = Cdt[k] * T[k] + asr[k] - o0 + b * T[k] + latent[k];
    }
    // Symmetric (cos-weighted) implicit system A x = b, solved by PCG with the zonal rows as preconditioner.
    // A_kk = c_j (Cdt + B) + w + e + n + s ;  zonal coupling D/(c_j h^2), meridional cos(phi_face) D / h^2
    for (let j = 0; j < NY; j++) {
      const cj = cosl[j], zx = 1 / (cj * h * h), mn = cosN[j] / (h * h), ms = j > 0 ? cosN[j - 1] / (h * h) : 0;
      for (let i = 0; i < NX; i++) {
        const k = j * NX + i, kw = j * NX + (i + NX - 1) % NX;
        cW[k] = zx * Dx[kw]; cE[k] = zx * Dx[k]; cN[k] = mn * Dn[k]; cS[k] = j > 0 ? ms * Dn[k - NX] : 0;
        cD[k] = cj * (Cdt[k] + Bk[k]) + cW[k] + cE[k] + cN[k] + cS[k];
        bv[k] = cj * Rk[k];
      }
    }
    const Aop = (x, y) => {
      for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
        const k = j * NX + i;
        let v = cD[k] * x[k] - cW[k] * x[j * NX + (i + NX - 1) % NX] - cE[k] * x[j * NX + (i + 1) % NX];
        if (j < NY - 1) v -= cN[k] * x[k + NX];
        if (j > 0) v -= cS[k] * x[k - NX];
        y[k] = v;
      }
    };
    const prec = (r, z) => {
      for (let j = 0; j < NY; j++) {
        for (let i = 0; i < NX; i++) { const k = j * NX + i; ta[i] = -cW[k]; tc[i] = -cE[k]; tb[i] = cD[k]; tr[i] = r[k]; }
        cyclicSolve(ta, tb, tc, tr, xs);
        for (let i = 0; i < NX; i++) z[j * NX + i] = xs[i];
      }
    };
    for (let k = 0; k < N; k++) xv[k] = T[k];
    Aop(xv, qv); let bn = 0;
    for (let k = 0; k < N; k++) { rv[k] = bv[k] - qv[k]; bn += bv[k] * bv[k]; }
    prec(rv, zv); let rz = 0; for (let k = 0; k < N; k++) { pv[k] = zv[k]; rz += rv[k] * zv[k]; }
    const tol = 1e-9 * bn;
    for (let it = 0; it < 40; it++) {
      Aop(pv, qv); let pq = 0; for (let k = 0; k < N; k++) pq += pv[k] * qv[k];
      const al = rz / pq; let rr = 0;
      for (let k = 0; k < N; k++) { xv[k] += al * pv[k]; rv[k] -= al * qv[k]; rr += rv[k] * rv[k]; }
      lastIters = it + 1;
      if (rr < tol) break;
      prec(rv, zv); let rz2 = 0; for (let k = 0; k < N; k++) rz2 += rv[k] * zv[k];
      const be = rz2 / rz; rz = rz2;
      for (let k = 0; k < N; k++) pv[k] = zv[k] + be * pv[k];
    }
    for (let k = 0; k < N; k++) T[k] = xv[k];
  }
  let lastIters = 0;
  const cW = new Float64Array(N), cE = new Float64Array(N), cN = new Float64Array(N), cS = new Float64Array(N), cD = new Float64Array(N), bv = new Float64Array(N);
  const xv = new Float64Array(N), rv = new Float64Array(N), zv = new Float64Array(N), pv = new Float64Array(N), qv = new Float64Array(N);

  // ---- low-level wind, diagnosed from the model's own temperature field.
  // Air is pulled toward warm columns (the thermal low over the dayside), turned by Ki-Engir's rotation
  // (6.097-day period) and slowed by surface friction: r v + f k x v = G grad T, G = g h_BL / T0.
  // Air chilled under an inversion drains downslope: the katabatic term g (Ta - Ts) / T0 grad z.
  const GRAV = 9.81 * 0.815, RM = 5.861e6, YR = 3.156e7;
  const uW = new Float32Array(N), vW = new Float32Array(N), divW = new Float32Array(N);
  const Wv = new Float32Array(N), Wn = new Float32Array(N), Pacc = new Float64Array(N), convQ = new Float32Array(N);
  const dI = new Int32Array(N), dJ = new Int32Array(N), dX = new Float32Array(N), dY = new Float32Array(N);
  const zs = k => lf[k] >= 0.5 ? zl[k] : 0;
  // Velocity potential: L chi = R^2 * target divergence, solved by PCG with the zonal rows as preconditioner.
  const chi = new Float64Array(N), pB = new Float64Array(N), pR = new Float64Array(N), pZ = new Float64Array(N), pP2 = new Float64Array(N), pQ = new Float64Array(N);
  function solvePoisson(target) {
    const eps = 1e-6;
    const Aop = (x, y) => {
      for (let j = 0; j < NY; j++) {
        const cj = cosl[j], zx = 1 / (cj * h * h), mn = cosN[j] / (h * h), ms = j > 0 ? cosN[j - 1] / (h * h) : 0;
        for (let i = 0; i < NX; i++) {
          const k = j * NX + i;
          let v = (cj * eps + 2 * zx + mn + ms) * x[k] - zx * (x[j * NX + (i + NX - 1) % NX] + x[j * NX + (i + 1) % NX]);
          if (j < NY - 1) v -= mn * x[k + NX]; if (j > 0) v -= ms * x[k - NX];
          y[k] = v;
        }
      }
    };
    const precond = (r, z) => {
      for (let j = 0; j < NY; j++) {
        const cj = cosl[j], zx = 1 / (cj * h * h), mn = cosN[j] / (h * h), ms = j > 0 ? cosN[j - 1] / (h * h) : 0;
        for (let i = 0; i < NX; i++) { ta[i] = -zx; tc[i] = -zx; tb[i] = cj * eps + 2 * zx + mn + ms; tr[i] = r[j * NX + i]; }
        cyclicSolve(ta, tb, tc, tr, xs);
        for (let i = 0; i < NX; i++) z[j * NX + i] = xs[i];
      }
    };
    let mean = 0, mw = 0; for (let k = 0; k < N; k++) { mean += target[k] * wgt[k]; mw += wgt[k]; } mean /= mw;
    let bn = 0; for (let k = 0; k < N; k++) { pB[k] = -wgt[k] * RM * RM * (target[k] - mean); bn += pB[k] * pB[k]; }
    Aop(chi, pQ); for (let k = 0; k < N; k++) pR[k] = pB[k] - pQ[k];
    precond(pR, pZ); let rz = 0; for (let k = 0; k < N; k++) { pP2[k] = pZ[k]; rz += pR[k] * pZ[k]; }
    for (let it = 0; it < 80; it++) {
      Aop(pP2, pQ); let pq = 0; for (let k = 0; k < N; k++) pq += pP2[k] * pQ[k];
      const al = rz / pq; let rr = 0;
      for (let k = 0; k < N; k++) { chi[k] += al * pP2[k]; pR[k] -= al * pQ[k]; rr += pR[k] * pR[k]; }
      if (rr < 1e-10 * bn) break;
      precond(pR, pZ); let rz2 = 0; for (let k = 0; k < N; k++) rz2 += pR[k] * pZ[k];
      const be = rz2 / rz; rz = rz2; for (let k = 0; k < N; k++) pP2[k] = pZ[k] + be * pP2[k];
    }
  }
  const tgtDiv = new Float64Array(N);
  function winds() {
    // 1. divergent flow from the weak-temperature-gradient balance: radiative heating of a column is
    //    carried off by ascent (low-level convergence), radiative cooling by subsidence (divergence)
    for (let k = 0; k < N; k++) {
      const net = asr[k] - olr(T[k], cloud(k));
      tgtDiv[k] = -net / (1004 * P.stabK) / (1.2 * P.inflowH);
    }
    solvePoisson(tgtDiv);
    // 2. turned by rotation and slowed by drag that grows with speed; katabatic drainage added on slopes
    const r0 = 1 / (P.fricDays * 86400);
    for (let j = 0; j < NY; j++) {
      const f = 2 * OMEGA * Math.sin(latc[j] * Math.PI / 180);
      const dx = 2 * RM * Math.max(cosl[j], 0.02) * h;
      const jn = Math.min(NY - 1, j + 1), js = Math.max(0, j - 1), dy = (jn - js) * RM * h;
      for (let i = 0; i < NX; i++) {
        const k = j * NX + i, kE = j * NX + (i + 1) % NX, kW = j * NX + (i + NX - 1) % NX, kN = jn * NX + i, kS = js * NX + i;
        const ud = (chi[kE] - chi[kW]) / dx, vd = (chi[kN] - chi[kS]) / dy;
        let Fx = r0 * ud, Fy = r0 * vd;
        const inv = lf[k] >= 0.5 ? Math.max(0, Tair[k] - Ts[k]) : 0;
        if (inv > 0) { const a = GRAV * inv / 280; Fx -= a * (zs(kE) - zs(kW)) / dx; Fy -= a * (zs(kN) - zs(kS)) / dy; }
        let u = ud, v = vd;
        for (let it = 0; it < 5; it++) {
          const re = r0 + P.dragCd * Math.hypot(u, v) / P.hBL, de = re * re + f * f;
          u = 0.5 * u + 0.5 * (re * Fx + f * Fy) / de; v = 0.5 * v + 0.5 * (re * Fy - f * Fx) / de;
        }
        uW[k] = u; vW[k] = v;
      }
    }
    for (let j = 0; j < NY; j++) {
      const jn = Math.min(NY - 1, j + 1), js = Math.max(0, j - 1);
      const cn = Math.cos(latc[jn] * Math.PI / 180), cs = Math.cos(latc[js] * Math.PI / 180), cj = Math.max(cosl[j], 0.02);
      for (let i = 0; i < NX; i++) {
        const k = j * NX + i;
        divW[k] = (uW[j * NX + (i + 1) % NX] - uW[j * NX + (i + NX - 1) % NX]) / (2 * RM * cj * h)
                + (vW[jn * NX + i] * cn - vW[js * NX + i] * cs) / ((jn - js) * RM * cj * h);
      }
    }
  }
  // semi-Lagrangian departure points on the sphere for one moisture sub-step
  function departures(dt) {
    for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
      const k = j * NX + i, la = latc[j] * Math.PI / 180, lo = lonc[i] * Math.PI / 180;
      const cl = Math.cos(la), sl = Math.sin(la), co = Math.cos(lo), so = Math.sin(lo);
      const X = [cl * co, cl * so, sl], E = [-so, co, 0], Nn = [-sl * co, -sl * so, cl];
      const V = [uW[k] * E[0] + vW[k] * Nn[0], uW[k] * E[1] + vW[k] * Nn[1], uW[k] * E[2] + vW[k] * Nn[2]];
      const sp = Math.hypot(V[0], V[1], V[2]);
      let latd = latc[j], lond = lonc[i];
      if (sp > 1e-6) {
        const al = sp * dt / RM, ca = Math.cos(al), sa = Math.sin(al);
        const Xd = [X[0] * ca - V[0] / sp * sa, X[1] * ca - V[1] / sp * sa, X[2] * ca - V[2] / sp * sa];
        latd = Math.asin(Math.max(-1, Math.min(1, Xd[2]))) * 180 / Math.PI; lond = Math.atan2(Xd[1], Xd[0]) * 180 / Math.PI;
      }
      const gx = (lond + 180) / RES - 0.5, gy = Math.max(0, Math.min(NY - 1.001, (latd + 90) / RES - 0.5));
      const i0 = Math.floor(gx), j0 = Math.floor(gy);
      dI[k] = ((i0 % NX) + NX) % NX; dJ[k] = Math.min(NY - 2, j0); dX[k] = gx - i0; dY[k] = gy - dJ[k];
    }
  }
  // ---- steady-state vapour transport, flux form, first-order upwind on the finite-volume grid:
  //   div(v W) + (1/tau + ascent + sat) W = E      (exact water conservation: total P = total E)
  // ascent = low-level convergence carries vapour up where it condenses; sat = an extra sink that
  // holds the column at saturation. Solved by BiCGSTAB with the zonal rows as preconditioner.
  const mD = new Float64Array(N), mW = new Float64Array(N), mE = new Float64Array(N), mN = new Float64Array(N), mS = new Float64Array(N);
  const satK = new Float64Array(N), sinkK = new Float64Array(N), Ev = new Float64Array(N), Wsol = new Float64Array(N);
  const bR = new Float64Array(N), bR0 = new Float64Array(N), bP = new Float64Array(N), bV = new Float64Array(N), bS = new Float64Array(N), bT = new Float64Array(N), bY = new Float64Array(N), bZ = new Float64Array(N);
  function buildTransport() {
    for (let j = 0; j < NY; j++) {
      const cj = Math.max(cosl[j], 1e-3), A = cj * h * h * RM; // area / RM (all faces scaled by RM too)
      const lz = h, lnN = cosN[j] * h, lnS = j > 0 ? cosN[j - 1] * h : 0;
      for (let i = 0; i < NX; i++) {
        const k = j * NX + i, ke = j * NX + (i + 1) % NX, kw = j * NX + (i + NX - 1) % NX;
        const ue = 0.5 * (uW[k] + uW[ke]), uw = 0.5 * (uW[k] + uW[kw]);
        const vn = j < NY - 1 ? 0.5 * (vW[k] + vW[k + NX]) : 0, vs = j > 0 ? 0.5 * (vW[k] + vW[k - NX]) : 0;
        // outward normal velocities: east +ue, west -uw, north +vn, south -vs
        let d = 0; mE[k] = 0; mW[k] = 0; mN[k] = 0; mS[k] = 0;
        const face = (un, len, set) => { const q = un * len / A; if (q > 0) d += q; else set(q); };
        face(ue, lz, q => { mE[k] = q; }); face(-uw, lz, q => { mW[k] = q; });
        face(vn, lnN, q => { mN[k] = q; }); face(-vs, lnS, q => { mS[k] = q; });
        mD[k] = d;
      }
    }
  }
  function solveTransport(maxit) {
    const Aop = (x, y) => {
      for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
        const k = j * NX + i;
        let v = (mD[k] + sinkK[k] + satK[k]) * x[k] + mE[k] * x[j * NX + (i + 1) % NX] + mW[k] * x[j * NX + (i + NX - 1) % NX];
        if (j < NY - 1) v += mN[k] * x[k + NX]; if (j > 0) v += mS[k] * x[k - NX];
        y[k] = v;
      }
    };
    const pre = (r, z) => {
      for (let j = 0; j < NY; j++) {
        for (let i = 0; i < NX; i++) { const k = j * NX + i; ta[i] = mW[k]; tc[i] = mE[k]; tb[i] = mD[k] + sinkK[k] + satK[k]; tr[i] = r[k]; }
        cyclicSolve(ta, tb, tc, tr, xs);
        for (let i = 0; i < NX; i++) z[j * NX + i] = xs[i];
      }
    };
    Aop(Wsol, bT); let bn = 0;
    for (let k = 0; k < N; k++) { bR[k] = Ev[k] - bT[k]; bR0[k] = bR[k]; bP[k] = 0; bV[k] = 0; bn += Ev[k] * Ev[k]; }
    let rho = 1, alpha = 1, omega = 1;
    for (let it = 0; it < maxit; it++) {
      let rho1 = 0; for (let k = 0; k < N; k++) rho1 += bR0[k] * bR[k];
      const beta = (rho1 / rho) * (alpha / omega); rho = rho1;
      for (let k = 0; k < N; k++) bP[k] = bR[k] + beta * (bP[k] - omega * bV[k]);
      pre(bP, bY); Aop(bY, bV);
      let r0v = 0; for (let k = 0; k < N; k++) r0v += bR0[k] * bV[k];
      alpha = rho / r0v;
      for (let k = 0; k < N; k++) bS[k] = bR[k] - alpha * bV[k];
      pre(bS, bZ); Aop(bZ, bT);
      let tt = 0, ts = 0; for (let k = 0; k < N; k++) { tt += bT[k] * bT[k]; ts += bT[k] * bS[k]; }
      omega = tt > 0 ? ts / tt : 0;
      let rr = 0;
      for (let k = 0; k < N; k++) { Wsol[k] += alpha * bY[k] + omega * bZ[k]; bR[k] = bS[k] - omega * bT[k]; rr += bR[k] * bR[k]; }
      if (rr < 1e-12 * bn || !isFinite(rr)) break;
    }
    for (let k = 0; k < N; k++) if (!(Wsol[k] >= 0)) Wsol[k] = 0;
  }
  let lastScale = 1;
  function moisture() {
    for (let k = 0; k < N; k++) {
      const tsea = T[k], si = 1 - smooth(269.35, 273.35, tsea);
      const rn = asr[k] * 0.62;
      const eo = (1 - si) * P.evap * (0.8 * rn * ptf(Math.max(tsea, 250)) * 12.9 + 8 * Math.max(0, tsea - 262));
      const petk = P.petScale * 1.26 * ptf(Math.max(Tair[k], 250)) * Math.max(rn, 0) * 12.9 + 4 * Math.max(0, Tair[k] - 263);
      pet[k] = petk;
      const veg = 0.25 + 0.75 * stage[k] / 5;
      landE[k] = P.recycle * veg * Math.min(prec[k], petk);
      evapF[k] = (1 - lf[k]) * eo + lf[k] * landE[k];
      Ev[k] = evapF[k] / YR;
    }
    winds(); buildTransport();
    const tau = P.tauRes * 86400;
    for (let k = 0; k < N; k++) sinkK[k] = 1 / tau + Math.max(-divW[k], 0);
    // saturation: a few passes, each raising the sink where the solved column exceeds its capacity
    for (let pass = 0; pass < 3; pass++) {
      solveTransport(pass === 0 ? 60 : 30);
      for (let k = 0; k < N; k++) {
        const ta = T[k] - P.lapse * zs(k) / 1000, cap = P.rhCrit * rhoVsat(Math.max(ta, 180)) * P.Hw;
        if (Wsol[k] > cap) satK[k] = (satK[k] + mD[k] + sinkK[k]) * (Wsol[k] / cap) - mD[k] - sinkK[k];
        else satK[k] *= 0.7;
      }
    }
    let sp = 0, se = 0;
    for (let k = 0; k < N; k++) {
      const p = Wsol[k] * (sinkK[k] + satK[k]);
      prec[k] = p * YR; sp += p * wgt[k]; se += Ev[k] * wgt[k];
    }
    lastScale = sp > 0 ? se / sp : 1;
    for (let k = 0; k < N; k++) { ai[k] = pet[k] > 1 ? prec[k] / pet[k] : (prec[k] > 0 ? 9 : 0); latent[k] = (prec[k] - evapF[k]) / 12.6; }
    Wv.set(Wsol);
  }
  // classes: 0 ocean, 1 sea ice, 2 ice, 3 tundra, then desert/steppe/subhumid/humid x cool/temperate/hot
  const CLS = ['open sea', 'sea ice', 'ice', 'tundra', 'cool desert', 'temperate desert', 'hot desert', 'cool steppe', 'temperate steppe', 'hot steppe', 'cool sub-humid', 'temperate sub-humid', 'savanna', 'cool humid', 'temperate humid', 'tropical humid'];
  const CAP = [0, 0, 0, 2, 1, 1, 1, 2, 2, 2, 3, 3.5, 3.5, 5, 5, 5];
  function classify() {
    const tPhoto = P.co2 <= 500 ? 0 : 11.7 * Math.log(P.co2 / 500) / Math.log(2); // CO2-dependent cool/temperate cut (vault: ~0 at 500, 11.7 at 1000, 24.6 at 2000)
    for (let k = 0; k < N; k++) {
      if (lf[k] < 0.5) { cls[k] = T[k] < 271.35 ? 1 : 0; cap[k] = 0; continue; }
      const tc = Tair[k] - 273.15;
      let c;
      if (tc < 0) c = 2; else if (tc < 6) c = 3; else {
        const band = tc < tPhoto ? 0 : tc < 27 ? 1 : 2;
        const a = ai[k], m = a < 0.2 ? 0 : a < 0.5 ? 1 : a < 0.65 ? 2 : 3;
        c = 4 + m * 3 + band;
      }
      cls[k] = c; cap[k] = CAP[c];
    }
  }
  function succession(dtYears) {
    for (let k = 0; k < N; k++) {
      if (lf[k] < 0.5) { stage[k] = 0; continue; }
      const tc = Tair[k] - 273.15;
      const g = Math.max(0, Math.min(1, tc / 18)) * Math.max(0.15, Math.min(1, ai[k] / 0.65));
      if (stage[k] < cap[k]) stage[k] = Math.min(cap[k], stage[k] + dtYears * g / P.succYears * (1 + 0.4 * (stage[k] < 1 ? 0 : 1)));
      else stage[k] = cap[k] + (stage[k] - cap[k]) * Math.exp(-dtYears / 6);
      if (P.disturb > 0 && stage[k] > 0.5 && rand() < 1 - Math.exp(-P.disturb * dtYears)) stage[k] = 0;
    }
  }

  // ---- river routing on a 0.5-degree grid: priority-flood fill (with a small gradient across flats),
  // D8 steepest descent on the filled surface, then water routed upstream-first. A depression is open
  // water that evaporates at the local PET: it spills when its inflow beats that, and otherwise ends
  // as a terminal lake that keeps its water.
  const RX = D.RX, RY = D.RY, RN = RX * RY, RH = D.RRES * Math.PI / 180;
  const rz = dec(D.rz, Int16Array), rl8 = dec(D.rland, Uint8Array);
  const rLand = new Uint8Array(RN), rFill = new Float64Array(RN).fill(Infinity), rRecv = new Int32Array(RN).fill(-1);
  const rOrder = new Int32Array(RN), rArea = new Float64Array(RN), rClim = new Int32Array(RN), rDep = new Uint8Array(RN);
  const rCatch = new Float64Array(RN), rRoot = new Int32Array(RN).fill(-1), rLat = new Float32Array(RY), rLon = new Float32Array(RX);
  let rCount = 0;
  {
    for (let j = 0; j < RY; j++) rLat[j] = -90 + D.RRES * (j + 0.5);
    for (let i = 0; i < RX; i++) rLon[i] = -180 + D.RRES * (i + 0.5);
    for (let j = 0; j < RY; j++) for (let i = 0; i < RX; i++) {
      const k = j * RX + i;
      rLand[k] = rl8[k] >= 128 ? 1 : 0;
      rArea[k] = (RKM * 1000) * (RKM * 1000) * Math.cos(rLat[j] * Math.PI / 180) * RH * RH;
      rClim[k] = idx(rLat[j], rLon[i]);
    }
    const nbs = (k, out) => {
      const j = (k / RX) | 0, i = k % RX; let n = 0;
      for (let dj = -1; dj <= 1; dj++) { const jj = j + dj; if (jj < 0 || jj >= RY) continue;
        for (let di = -1; di <= 1; di++) { if (!dj && !di) continue; out[n++] = jj * RX + (i + di + RX) % RX; } }
      return n;
    };
    const dist = (a, b) => { // metres between neighbour centres
      const ja = (a / RX) | 0, jb = (b / RX) | 0, ia = a % RX, ib = b % RX;
      const di = Math.min(Math.abs(ia - ib), RX - Math.abs(ia - ib)), cl = Math.cos((rLat[ja] + rLat[jb]) / 2 * Math.PI / 180);
      return RKM * 1000 * RH * Math.hypot(di * cl, ja - jb);
    };
    // binary heap keyed on filled height
    const hk = [], hv = [];
    const push = (z, k) => { hk.push(z); hv.push(k); let i = hk.length - 1; while (i > 0) { const p = (i - 1) >> 1; if (hk[p] <= hk[i]) break; [hk[p], hk[i]] = [hk[i], hk[p]]; [hv[p], hv[i]] = [hv[i], hv[p]]; i = p; } };
    const pop = () => { const k = hv[0], z = hk[0], lz = hk.pop(), lv = hv.pop(); if (hk.length) { hk[0] = lz; hv[0] = lv; let i = 0; for (;;) { const l = 2 * i + 1, r = l + 1; let m = i; if (l < hk.length && hk[l] < hk[m]) m = l; if (r < hk.length && hk[r] < hk[m]) m = r; if (m === i) break; [hk[m], hk[i]] = [hk[i], hk[m]]; [hv[m], hv[i]] = [hv[i], hv[m]]; i = m; } } return [z, k]; };
    const nb = new Int32Array(8);
    for (let k = 0; k < RN; k++) if (!rLand[k]) {
      rFill[k] = Math.min(rz[k], 0);
      const n = nbs(k, nb); let coast = false; for (let q = 0; q < n; q++) if (rLand[nb[q]]) coast = true;
      if (coast) push(rFill[k], k);
    }
    while (hk.length) {
      const [z, k] = pop();
      if (rLand[k]) rOrder[rCount++] = k;
      const n = nbs(k, nb);
      for (let q = 0; q < n; q++) { const m = nb[q]; if (!rLand[m] || rFill[m] !== Infinity) continue;
        rFill[m] = Math.max(rz[m], z + 1e-3); push(rFill[m], m); }
    }
    // D8 receivers on the filled surface
    for (let k = 0; k < RN; k++) {
      if (!rLand[k] || rFill[k] === Infinity) continue;
      const n = nbs(k, nb); let best = -1, bs = 0;
      for (let q = 0; q < n; q++) { const m = nb[q]; const sl = (rFill[k] - rFill[m]) / dist(k, m); if (sl > bs) { bs = sl; best = m; } }
      rRecv[k] = best;
      rDep[k] = rFill[k] > rz[k] + 10 ? 1 : 0;
    }
    // catchment areas and the outlet each cell drains to
    for (let q = rCount - 1; q >= 0; q--) { const k = rOrder[q]; rCatch[k] += rArea[k]; const r = rRecv[k]; if (r >= 0 && rLand[r]) rCatch[r] += rCatch[k]; }
    for (let q = 0; q < rCount; q++) { const k = rOrder[q]; const r = rRecv[k]; rRoot[k] = r >= 0 && rLand[r] ? rRoot[r] : k; }
  }
  const rQ = new Float64Array(RN), rLakeF = new Float32Array(RN);
  function route() {
    rQ.fill(0); rLakeF.fill(0);
    for (let q = 0; q < rCount; q++) {
      const k = rOrder[q], c = rClim[k];
      const run = Tair[c] > 273.15 ? Math.max(0, prec[c] - landE[c]) : 0; // frozen ground stores its water as ice
      rQ[k] = run / 1000 / YR * rArea[k];
    }
    for (let q = rCount - 1; q >= 0; q--) {
      const k = rOrder[q];
      if (rDep[k]) {
        const demand = pet[rClim[k]] / 1000 / YR * rArea[k];
        if (rQ[k] > 0) rLakeF[k] = demand > 0 ? Math.min(1, rQ[k] / demand) : 1;
        rQ[k] = Math.max(0, rQ[k] - demand);
      }
      const r = rRecv[k]; if (r >= 0 && rLand[r]) rQ[r] += rQ[k];
    }
  }
  function riverStats() {
    let qmax = 0, terminal = 0, through = 0, area = 0, rootA = new Map();
    for (let q = 0; q < rCount; q++) {
      const k = rOrder[q], r = rRecv[k];
      if (r < 0 || !rLand[r]) qmax = Math.max(qmax, rQ[k]);
      if (rLakeF[k] > 0) { if (rLakeF[k] < 1 || rQ[k] === 0) terminal += rArea[k] * rLakeF[k]; else through += rArea[k]; }
      const c = rClim[k]; if (mu[c] > 0) { area += rArea[k]; rootA.set(rRoot[k], (rootA.get(rRoot[k]) || 0) + rArea[k]); }
    }
    let big = 0; for (const v of rootA.values()) big = Math.max(big, v);
    return { qmax, terminalKm2: terminal / 1e6, throughKm2: through / 1e6, daysideLandKm2: area / 1e6, biggestCatchKm2: big / 1e6 };
  }

  function step(dtYears, fast) {
    refreshStar();
    surfaceTemps(); radiation(); energyStep(dtYears); surfaceTemps(); moisture(); classify();
    if (!fast) { succession(dtYears); route(); }
    years += dtYears;
  }
  function reset(opts = {}) {
    T.fill(opts.T0 || 250); stage.fill(0); prec.fill(0); landE.fill(0); latent.fill(0); alb.fill(0.4); Wv.fill(0); chi.fill(0); Wsol.fill(0); satK.fill(0); years = 0; rng = 12345;
    if (opts.mature) { /* spin climate first, then set stages to climax */ }
  }

  // ---- diagnostics against canon gates and the GCM
  const SITES = {
    quilatria: idx(39, -50), quedros: idx(55, -25), silver: idx(23, -21),
  };
  const highlands = [], nightLow = [];
  for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
    const k = j * NX + i;
    if (latc[j] > 60 && Math.abs(lonc[i]) <= 5 && lf[k] >= 0.5) highlands.push(k);
    if (mu[k] < 0 && lf[k] >= 0.5 && zl[k] < 1000) nightLow.push(k);
  }
  // sun-elevation profile: 10-degree bands from -90 to 90
  const NBAND = 18;
  function bandProfile(field) {
    const sum = new Float64Array(NBAND), w = new Float64Array(NBAND);
    for (let k = 0; k < N; k++) { const b = Math.min(NBAND - 1, Math.floor((sunEl[k] + 90) / 10)); sum[b] += field[k] * wgt[k]; w[b] += wgt[k]; }
    return Array.from(sum, (v, i) => v / w[i]);
  }
  const gcmBands = { rc067: bandProfile(gcm.rc067), rc070: bandProfile(gcm.rc070) };
  function rainStats() {
    let lp = 0, lg = 0, lw = 0, ea = 0, eb = 0, ab = 0, aa = 0, bb = 0, ww = 0, ta = 0;
    for (let k = 0; k < N; k++) {
      const w = wgt[k];
      if (lf[k] >= 0.5) { lp += prec[k] * w; lg += gcm.pr[k] * w; lw += w; }
      const a = Math.log10(prec[k] + 30), b = Math.log10(gcm.pr[k] + 30);
      ea += a * w; eb += b * w; ww += w; ta += (Tair[k] - gcm.tas[k]) ** 2 * w;
    }
    ea /= ww; eb /= ww;
    for (let k = 0; k < N; k++) { const w = wgt[k], a = Math.log10(prec[k] + 30) - ea, b = Math.log10(gcm.pr[k] + 30) - eb; ab += a * b * w; aa += a * a * w; bb += b * b * w; }
    return { landRain: lp / lw, gcmLandRain: lg / lw, rainCorr: ab / Math.sqrt(aa * bb), airRms: Math.sqrt(ta / ww), waterScale: lastScale };
  }
  function stats() {
    let asrSum = 0, g = 0, gw = 0, d = 0, dw = 0, n = 0, nw = 0, mx = -1e9, mn = 1e9, e67 = 0, e70 = 0, dl = 0, dlw = 0, des = 0, hum = 0, ste = 0;
    for (let k = 0; k < N; k++) {
      const w = wgt[k], t = Ts[k];
      g += t * w; gw += w; asrSum += asr[k] * w;
      if (mu[k] > 0) { d += t * w; dw += w; } else { n += t * w; nw += w; }
      if (Ts[k] > mx) mx = Ts[k]; if (Ts[k] < mn) mn = Ts[k];
      e67 += (t - gcm.rc067[k]) ** 2 * w; e70 += (t - gcm.rc070[k]) ** 2 * w;
      if (mu[k] > 0 && lf[k] >= 0.5) {
        dlw += w; const c = cls[k];
        if (c >= 4 && c <= 6) des += w; else if (c >= 7 && c <= 9) ste += w; else if (c >= 13) hum += w;
      }
    }
    let hmax = -1e9; for (const k of highlands) hmax = Math.max(hmax, Tair[k]);
    let nf = 0; for (const k of nightLow) if (Ts[k] < 273.15) nf++;
    const frost = 3031 / (3031 / 194.7 - Math.log(P.co2 * 1e-6 * 101325 / 101325)) ; // CO2 frost point, K
    // extra diagnostics for the star/orbit sweep
    let oc = 0, ocOpen = 0, b1 = 0, b1w = 0, b2 = 0, b2w = 0, cw = 0;
    for (let k = 0; k < N; k++) {
      const w = wgt[k]; cw += colW[k] * w;
      if (lf[k] < 0.5) { oc += w; if (T[k] >= 271.35) ocOpen += w; }
      else if (sunEl[k] >= 15 && sunEl[k] < 30) { b1 += Ts[k] * w; b1w += w; }
      else if (sunEl[k] >= 30 && sunEl[k] < 60) { b2 += Ts[k] * w; b2w += w; }
    }
    return {
      S0: P.S0, period: P.period, orbitAU: star.a, albedo: 1 - (asrSum / gw) / (P.S0 / 4), iceFree: ocOpen / oc,
      land1530: b1 / b1w, land3060: b2 / b2w, colWater: cw / gw, starRatio: star.ratio.slice(), fvis: star.fvis,
      years, asrGlobal: asrSum / gw, global: g / gw, day: d / dw, night: n / nw, max: mx, min: mn,
      rms067: Math.sqrt(e67 / gw), rms070: Math.sqrt(e70 / gw),
      quilatria: T[SITES.quilatria], quedros: T[SITES.quedros], silver: T[SITES.silver] - P.lapse * 3.215, highlands: hmax,
      nightFrozen: nf / Math.max(1, nightLow.length), frost,
      desert: des / dlw, steppe: ste / dlw, humid: hum / dlw, bands: bandProfile(Ts), ...rainStats(),
    };
  }

  return { uW, vW, divW, Wv, river: { RX, RY, RRES: D.RRES, rLand, rRecv, rQ, rLakeF, rDep, rCatch, rOrder, count: () => rCount, rLat, rLon }, riverStats, cosLat: j => cosl[j], coastKm, gcmBands, bandProfile, iters: () => lastIters, NX, NY, RES, N, P, T, Ts, Tair, latent, evapF, alb, asr, prec, pet, ai, stage, cls, iceF, lf, zl, zmean, mu, sunEl, gcm, lonc, latc, CLS, idx, step, reset, stats, surfaceTemps, moisture, classify };
}
if (typeof module !== 'undefined') module.exports = { createKiModel };
