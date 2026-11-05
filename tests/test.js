// Tests du simulateur d'IRM (node tests/test.js)
'use strict';
const I = require('../irm.js');
const T = I.TISSUES;
let n = 0, ko = 0;
const fails = [];
function ok(cond, msg) { n++; if (!cond) { ko++; fails.push(msg); } }
const near = (a, b, tol) => Math.abs(a - b) <= tol;
const rel = (a, b, tol) => Math.abs(a - b) <= tol * Math.max(Math.abs(a), Math.abs(b), 1e-12);
const R = I.rng(2024);
const uni = (a, b) => a + (b - a) * R();
const BRAIN = I.PHANTOMS.brain, QA = I.PHANTOMS.qa;
const single = (x, y, r, t) => ({ fov: 240, ellipses: [{ x, y, a: r, b: r, ang: 0, t: t || 'wm', p: null }], metalAt: null, extent: Math.max(Math.abs(x), Math.abs(y)) + r + 2 });
function image(ph, seq, acq, opts) { const K = I.acquire(ph, seq, Object.assign({ noise: false }, acq), opts); return { K, img: I.reconstruct(K), mag: null }; }
function mag(ph, seq, acq, opts) { const r = image(ph, seq, acq, opts); r.mag = I.magnitude(r.img); return r; }
function centroid(m, N, mask) { let sx = 0, sy = 0, s = 0; for (let row = 0; row < N; row++) for (let j = 0; j < N; j++) { const i = row * N + j; if (mask && !mask[i]) continue; const v = m[i]; sx += v * j; sy += v * row; s += v; } return [sx / s, sy / s]; }

// ------------------------------------------------------------ 1. constantes
ok(near(I.larmor(), 63.8655, 1e-3), 'fréquence de Larmor à 1,5 T ≈ 63,87 MHz');
ok(near(I.fatHz(), -223.53, 0.01), 'graisse à −3,5 ppm ≈ −224 Hz à 1,5 T');
ok(near(1000 / (2 * Math.abs(I.fatHz())), 2.237, 0.001), 'opposition de phase eau/graisse à 2,24 ms');

// ------------------------------------------------------------ 2. Bloch exact contre les formules
const tissues = ['wm', 'gm', 'csf', 'fat', 'muscle', 'blood', 'edema', 'hematoma'];
for (let k = 0; k < 60; k++) {
  const id = tissues[k % tissues.length], TR = uni(200, 10000), TE = uni(5, Math.min(250, TR * 0.8));
  const b = I.steadyState(id, { kind: 'SE', TR, TE, ETL: 1 }).echoes[0];
  ok(rel(b, I.formula.se(T[id], TR, TE), 1e-9), 'écho de spin ' + id + ' TR ' + TR.toFixed(0) + ' TE ' + TE.toFixed(0));
}
for (let k = 0; k < 60; k++) {
  const id = tissues[k % tissues.length], TR = uni(1500, 12000), TE = uni(5, 150), TI = uni(20, TR - TE - 50);
  const b = I.steadyState(id, { kind: 'IR', TR, TE, TI, ETL: 1 }).echoes[0];
  ok(near(b, I.formula.ir(T[id], TR, TI, TE), 1e-9), 'inversion-récupération ' + id + ' TI ' + TI.toFixed(0));
}
for (let k = 0; k < 60; k++) {
  const id = tissues[k % tissues.length], TR = uni(5, 1000), TE = uni(1, Math.min(40, TR * 0.9)), flip = uni(2, 90);
  const b = I.steadyState(id, { kind: 'GRE', TR, TE, flip }).echoes[0];
  ok(rel(b, I.formula.gre(T[id], TR, TE, flip), 1e-9), 'écho de gradient ' + id + ' α ' + flip.toFixed(0));
}
// TSE : en CPMG idéal, l'écho k vaut le premier × e^(−(k−1)·ESP/T2)
for (const id of tissues) {
  const s = I.steadyState(id, { kind: 'SE', TR: 4000, TE: 100, ETL: 16, ESP: 10 }).echoes;
  ok(s.every((v, k) => rel(v, s[0] * Math.exp((-k * 10) / T[id].T2), 1e-9)), 'train d\'échos ' + id + ' : décroissance T2 exacte');
}
// angle d'Ernst : maximum du signal en écho de gradient
for (const id of tissues) {
  const TR = 50, ae = I.formula.ernst(T[id], TR), f = (a) => I.formula.gre(T[id], TR, 0, a);
  ok(f(ae) >= f(ae - 1) && f(ae) >= f(ae + 1), 'angle d\'Ernst ' + id + ' (' + ae.toFixed(1) + '°)');
}
// TI d'annulation : formule et Bloch (bisection) d'accord
for (let k = 0; k < 20; k++) {
  const id = tissues[k % tissues.length], TR = uni(3000, 12000), TE = uni(10, 100);
  const a = I.formula.tiNull(T[id], TR, TE), b = I.tiNullBloch(id, { kind: 'IR', TR, TE, ETL: 1 });
  ok(b !== null && near(a, b, 1e-6), 'TI d\'annulation ' + id + ' : ' + a.toFixed(1) + ' ms');
  ok(near(I.formula.ir(T[id], TR, a, TE), 0, 1e-12), 'signal nul au TI d\'annulation ' + id);
}
ok(near(I.formula.tiNull(T.fat, 1e9, 0), T.fat.T1 * Math.LN2, 1e-6), 'STIR : TI = T1·ln 2 quand TR ≫ T1');
// régime établi atteint
ok(I.steadyState('csf', { kind: 'GRE', TR: 10, TE: 5, flip: 10 }).iters < 2000, 'régime établi atteint (LCR, TR court)');

// ------------------------------------------------------------ 3. ensemble d'isochromats
{
  const id = 'wm', seq = { kind: 'GRE', TR: 400, TE: 20, flip: 90 }, r2p = 30;
  const en = I.ensemble(id, seq, { n: 4000, r2p, samples: 400 });
  let worst = 0;
  en.times.forEach((t, i) => { if (t > 0 && t < 2000 / r2p) { const m = Math.hypot(en.mean[i][0], en.mean[i][1]), th = en.mz0 * Math.exp(-t / T[id].T2) * Math.exp((-r2p * t) / 1000); worst = Math.max(worst, Math.abs(m - th)); } });
  ok(worst < 0.01 * en.mz0, 'T2* : la somme des isochromats décroît en e^(−t/T2)·e^(−R2\'t) (écart ' + (worst / en.mz0).toExponential(1) + ')');
}
for (const TE of [20, 60, 120]) {
  const id = 'gm', seq = { kind: 'SE', TR: 2000, TE, ETL: 1 }, en = I.ensemble(id, seq, { n: 500, r2p: 40 });
  const i = en.times.findIndex((t) => Math.abs(t - TE) < 1e-9), m = Math.hypot(en.mean[i][0], en.mean[i][1]);
  ok(rel(m, en.mz0 * Math.exp(-TE / T[id].T2), 1e-9), 'écho de spin : tous les isochromats se rephasent à TE = ' + TE + ' ms');
}
{
  const offs = I.lorentzOffsets(1001, 20);
  ok(near(offs[500], 0, 1e-12) && offs.every((v, i) => i === 0 || v > offs[i - 1]), 'quantiles de Lorentz centrés et croissants');
}

// ------------------------------------------------------------ 4. FFT, Bessel, transformée de l'ellipse
for (const N of [8, 16, 64, 256]) {
  const re = Float64Array.from({ length: N }, () => uni(-1, 1)), im = Float64Array.from({ length: N }, () => uni(-1, 1));
  const [Rr, Ri] = I.dft(re, im, false), a = Float64Array.from(re), b = Float64Array.from(im); I.fft(a, b, false);
  ok(a.every((v, k) => near(v, Rr[k], 1e-9)) && b.every((v, k) => near(v, Ri[k], 1e-9)), 'FFT = transformée directe (N = ' + N + ')');
  I.fft(a, b, true); ok(a.every((v, k) => near(v / N, re[k], 1e-12)), 'FFT aller-retour (N = ' + N + ')');
}
{
  const J1int = (x) => { let s = 0; const M = 2000; for (let k = 0; k <= M; k++) { const t = (Math.PI * k) / M, w = k === 0 || k === M ? 1 : k % 2 ? 4 : 2; s += w * Math.cos(t - x * Math.sin(t)); } return (s * Math.PI) / (3 * M) / Math.PI; };
  for (let k = 0; k < 40; k++) { const x = uni(-30, 30); ok(near(I.besselJ1(x), J1int(x), 2e-7), 'J1(' + x.toFixed(2) + ') = représentation intégrale'); }
  ok(near(I.besselJ1(3.8317059702), 0, 1e-7) && near(I.besselJ1(7.0155866698), 0, 1e-7), 'zéros de J1 : 3,8317 et 7,0156');
  ok(near(I.jinc(0), 1, 1e-15), 'jinc(0) = 1');
}
{
  // transformée de Fourier d'une ellipse tournée et décalée contre une intégration numérique fine
  const e = { x: 12, y: -7, a: 20, b: 11, ang: 33 }, h = 0.1;
  for (const [kx, ky] of [[0, 0], [0.01, 0], [0.02, 0.015], [-0.03, 0.04], [0.05, -0.02]]) {
    let sr = 0, si = 0;
    for (let x = e.x - 21; x <= e.x + 21; x += h) for (let y = e.y - 21; y <= e.y + 21; y += h) if (I.inEllipse(e, x + h / 2, y + h / 2)) { const a = -2 * Math.PI * (kx * (x + h / 2) + ky * (y + h / 2)); sr += Math.cos(a) * h * h; si += Math.sin(a) * h * h; }
    const [fr, fi] = I.ellipseFT(e, kx, ky);
    ok(Math.hypot(fr - sr, fi - si) < 2e-3 * Math.PI * e.a * e.b, 'transformée de l\'ellipse exacte en k = (' + kx + ', ' + ky + ')');
  }
}

// ------------------------------------------------------------ 5. patients virtuels
for (const [name, ph] of [['cerveau', BRAIN], ['contrôle', QA]]) {
  const ells = I.withMetal(ph, true);
  let nest = true, disjoint = true, additive = true;
  const vals = Object.fromEntries(I.TISSUE_IDS.map((t, i) => [t, i + 1]));
  for (let k = 0; k < 20000; k++) {
    const x = uni(-ph.extent, ph.extent), y = uni(-ph.extent, ph.extent);
    const inside = ells.map((e) => I.inEllipse(e, x, y));
    ells.forEach((e, i) => { if (inside[i] && e.p !== null && !inside[e.p]) nest = false; });
    ells.forEach((e, i) => ells.forEach((f, j) => { if (j > i && inside[i] && inside[j] && e.p === f.p) disjoint = false; }));
    let v = 0; ells.forEach((e, i) => { if (inside[i]) v += vals[e.t] - (e.p === null ? 0 : vals[ells[e.p].t]); });
    if (v !== (I.tissueAt(ells, x, y) === 'air' ? 0 : vals[I.tissueAt(ells, x, y)])) additive = false;
  }
  ok(nest, name + ' : chaque ellipse est dans son parent');
  ok(disjoint, name + ' : les ellipses sœurs ne se chevauchent pas');
  ok(additive, name + ' : somme des contributions = tissu le plus profond');
}

// ------------------------------------------------------------ 6. image : valeurs mesurées = théorie
for (const p of ['t1', 't1ir', 'gre']) {
  const seq = I.PRESETS[p], acq = { N: 128, fovx: 200, fovy: 200 }, { img, mag: m } = mag(QA, seq, acq);
  I.QA_ORDER.forEach((id, i) => {
    if (id === 'bone') return;
    const mask = I.roiMask(QA.ellipses, img, i, 6), mean = I.stats(m, mask).mean, th = I.theory(id, seq, acq);
    ok(rel(mean, th, 0.02), p + ' : tube ' + id + ' mesuré ' + mean.toFixed(4) + ' / théorie ' + th.toFixed(4));
  });
}
{
  // écho de spin rapide : le contraste suit l'écho effectif (à quelques % près, modulation du train d'échos)
  for (const p of ['t2', 'pd']) {
  const seq = I.PRESETS[p], acq = { N: 128, fovx: 200, fovy: 200 }, { img, mag: m } = mag(QA, seq, acq);
  I.QA_ORDER.forEach((id, i) => { if (id === 'bone') return; const mean = I.stats(m, I.roiMask(QA.ellipses, img, i, 6)).mean, th = I.theory(id, seq, acq); ok(rel(mean, th, 0.07), p + ' turbo : tube ' + id + ' ' + mean.toFixed(4) + ' / ' + th.toFixed(4)); });
  }
}
{
  // phénomène de Gibbs : dépassement d'environ 9 % au bord d'un disque
  const ph = single(0, 0, 60), seq = I.PRESETS.t1, acq = { N: 256, r2p: 0 }, { mag: m } = mag(ph, seq, acq);
  let mx = 0; for (let j = 0; j < 256; j++) mx = Math.max(mx, m[128 * 256 + j]);
  const ov = mx / I.theory('wm', seq, acq) - 1;
  ok(ov > 0.085 && ov < 0.095, 'Gibbs : dépassement ' + (ov * 100).toFixed(2) + ' % (≈ 9 %)');
  const h = mag(ph, seq, Object.assign({ filter: 'hamming' }, acq)).mag; let mh = 0; for (let j = 0; j < 256; j++) mh = Math.max(mh, h[128 * 256 + j]);
  ok(mh / I.theory('wm', seq, acq) - 1 < 0.01, 'filtre de Hamming : plus d\'oscillations');
}

// ------------------------------------------------------------ 7. contrastes connus en clinique (cerveau)
{
  const vals = {};
  for (const p of ['t1', 't1ir', 't2', 'flair', 'stir', 'gre']) vals[p] = Object.fromEntries(['wm', 'gm', 'csf', 'fat', 'edema', 'hematoma'].map((id) => [id, I.theory(id, I.PRESETS[p], {})]));
  ok(vals.t1.csf < 0.6 * vals.t1.wm, 'T1 : liquide céphalo-rachidien sombre');
  ok(vals.t1.fat > 1.5 * vals.t1.wm, 'T1 : graisse claire');
  ok(vals.t1ir.wm > 1.2 * vals.t1ir.gm, 'T1 par inversion : substance blanche plus claire que la grise');
  ok(vals.t2.csf > 2 * vals.t2.wm, 'T2 : liquide céphalo-rachidien clair');
  ok(vals.t2.edema > 1.5 * vals.t2.wm, 'T2 : œdème clair');
  ok(vals.flair.csf < 0.05 * vals.flair.wm, 'FLAIR : liquide céphalo-rachidien éteint');
  ok(vals.flair.edema > 1.3 * vals.flair.wm, 'FLAIR : œdème clair');
  ok(vals.stir.fat < 0.05 * vals.stir.wm, 'STIR : graisse éteinte');
  ok(vals.gre.hematoma < 0.5 * vals.gre.wm, 'écho de gradient T2* : hématome sombre');
  const fl = mag(BRAIN, I.PRESETS.flair, { N: 128 }), mWM = I.stats(fl.mag, I.roiMask(BRAIN.ellipses, fl.img, 6, 4)).mean, mV = I.stats(fl.mag, I.roiMask(BRAIN.ellipses, fl.img, 7, 2.5)).mean;
  ok(mV < 0.15 * mWM, 'FLAIR, image : ventricules éteints (' + (mV / mWM * 100).toFixed(1) + ' % de la substance blanche)');
}

// ------------------------------------------------------------ 8. bruit et rapport signal sur bruit
{
  const seq = I.PRESETS.t1;
  const sdNoise = (acq) => { const K = I.acquire(QA, seq, Object.assign({ fovx: 200, fovy: 200 }, acq)); const a = I.reconstruct(K), b = I.reconstruct(K, 'clean'); const d = a.re.map((v, i) => v - b.re[i]); return [I.stats(d).sd, K.sigma]; };
  const base = { N: 128, bw: 130, thk: 5, nex: 1 };
  const [s0, p0] = sdNoise(base);
  ok(rel(s0, p0, 0.03), 'bruit mesuré = bruit prévu (' + s0.toExponential(3) + ' / ' + p0.toExponential(3) + ')');
  for (const [lab, ch, fac] of [['moyennages ×4', { nex: 4 }, 1 / 2], ['épaisseur ×2', { thk: 10 }, 1 / 2], ['bande passante ×4', { bw: 520 }, 2], ['matrice 64', { N: 64 }, Math.pow(0.5, 1.5)]]) {
    const [s] = sdNoise(Object.assign({}, base, ch, { seed: 99 }));
    ok(rel(s / s0, fac, 0.05), 'loi du rapport signal sur bruit : ' + lab + ' → bruit × ' + fac.toFixed(3) + ' (mesuré × ' + (s / s0).toFixed(3) + ')');
  }
  const snr = I.predictedSNR('wm', seq, { N: 256, fovx: 240, fovy: 240 });
  ok(snr > 15 && snr < 60, 'T1 256², 5 mm : rapport signal sur bruit réaliste dans la substance blanche (' + snr.toFixed(1) + ')');
}

// ------------------------------------------------------------ 9. déplacement chimique
for (const bw of [40, 80, 160, 320]) {
  for (const engine of ['exact', 'grille']) {
    const seq = { kind: 'SE', TR: 500, TE: 15, ETL: 1 }, acq = { N: 128, fovx: 200, fovy: 200, bw, r2p: 0 };
    const { img, mag: m } = mag(QA, seq, acq, engine === 'grille' ? { grid: true } : null);
    const box = (cx, cy) => { const mk = new Uint8Array(128 * 128); for (let row = 0; row < 128; row++) for (let j = 0; j < 128; j++) { const [x, y] = I.pixelXY(img, row, j); if (Math.abs(x - cx) < 31 && Math.abs(y - cy) < 31) mk[row * 128 + j] = 1; } return mk; };
    // tube de graisse en (−52, 0), tube de substance blanche (référence) en (−52, 52)
    const [fx, fy] = centroid(m, 128, box(-52, 0)), [wx, wy] = centroid(m, 128, box(-52, 52));
    const shift = fx - wx;
    ok(near(shift, -(I.fatHz() / bw), 0.12), engine + ' : graisse déplacée de ' + shift.toFixed(2) + ' px en lecture à ' + bw + ' Hz/px (prévu ' + (-(I.fatHz() / bw)).toFixed(2) + ')');
    ok(near(fy - wy, 52 / (200 / 128), 0.05), engine + ' : aucun déplacement en phase à ' + bw + ' Hz/px');
  }
}
{
  // phase et opposition de phase eau / graisse en écho de gradient
  const ph = (te) => (2 * Math.PI * I.fatHz() * te) / 1000;
  ok(near(Math.cos(ph(I.PRESETS.opposed.TE)), -1, 2e-3), 'TE 2,24 ms : eau et graisse en opposition');
  ok(near(Math.cos(ph(I.PRESETS.inphase.TE)), 1, 2e-3), 'TE 4,47 ms : eau et graisse en phase');
  // un voxel moitié eau moitié graisse s'éteint en opposition de phase
  const mix = { fov: 200, extent: 60, metalAt: null, ellipses: [{ x: 0, y: 0, a: 50, b: 50, ang: 0, t: 'muscle', p: null }, { x: 0, y: 0, a: 30, b: 30, ang: 0, t: 'fat', p: 0 }] };
  const edge = (p) => { const r = mag(mix, I.PRESETS[p], { N: 128, fovx: 200, fovy: 200, bw: 400 }); let lo = Infinity; for (let j = 64; j < 128; j++) { const [x] = I.pixelXY(r.img, 64, j); if (x > 25 && x < 35) lo = Math.min(lo, r.mag[64 * 128 + j]); } return lo; };
  ok(edge('opposed') < 0.5 * edge('inphase'), 'opposition de phase : liseré noir à l\'interface eau / graisse');
}

// ------------------------------------------------------------ 10. repliement
{
  const seq = I.PRESETS.t1;
  const r1 = mag(single(0, 90, 8), seq, { N: 128, fovx: 240, fovy: 120 });
  const [, cy] = centroid(r1.mag, 128), yApp = I.pixelXY(r1.img, cy, 0)[1];
  ok(near(yApp, 90 - 120, 1.5), 'phase : objet à y = 90 mm replié à y = ' + yApp.toFixed(1) + ' mm (prévu −30)');
  const r2 = mag(single(90, 0, 8), seq, { N: 128, fovx: 120, fovy: 240 });
  let mx = 0; for (const v of r2.mag) mx = Math.max(mx, v);
  ok(mx < 0.02 * I.theory('wm', seq, {}), 'lecture : sur-échantillonnage, pas de repliement (' + (mx / I.theory('wm', seq, {}) * 100).toFixed(2) + ' %)');
  const r3 = mag(single(90, 0, 8), seq, { N: 128, fovx: 120, fovy: 240 }, { grid: true }); let m3 = 0; for (const v of r3.mag) m3 = Math.max(m3, v);
  ok(m3 < 0.02 * I.theory('wm', seq, {}), 'lecture, moteur sur grille : pas de repliement');
}

// ------------------------------------------------------------ 11. mouvement
{
  const seq = I.PRESETS.t1, N = 128;
  for (const P of [2, 4, 8]) {
    const acq = { N, fovx: 240, fovy: 240, motion: { type: 'breath', amp: 1.5, period: P } };
    const r = mag(single(0, 0, 5), seq, acq);
    const g = Math.round((N * seq.TR) / (P * 1000)); // décalage prévu des fantômes, en pixels
    const col = N / 2, prof = Array.from({ length: N }, (_, row) => r.mag[row * N + col]);
    // fantôme repéré au maximum du profil lissé sur 7 pixels (le fantôme d'un petit disque a deux lobes autour de sa position)
    const c = N / 2 - 1; let best = -1, bv = -1;
    for (let row = 3; row < N - 3; row++) { if (Math.abs((N / 2 - 1 - row) * (240 / N)) <= 12) continue; let sm = 0; for (let q = -3; q <= 3; q++) sm += prof[row + q]; if (sm > bv) { bv = sm; best = row; } }
    const d = Math.abs(best - c);
    ok(Math.abs(d - g) <= 1, 'respiration ' + P + ' s : fantôme à ' + d + ' px (prévu ' + g + ' px = N·TR/période)');
  }
  const still = mag(single(0, 0, 25), seq, { N }), jump = mag(single(0, 0, 25), seq, { N, motion: { type: 'jump', amp: 5, period: 4 } });
  let e0 = 0, e1 = 0; for (let row = 0; row < N; row++) { const y = (N / 2 - 1 - row) * (240 / N); if (Math.abs(y) > 40) { e0 += still.mag[row * N + 64]; e1 += jump.mag[row * N + 64]; } }
  ok(e1 > 5 * e0 + 1e-6, 'mouvement brusque : signal étalé en dehors de l\'objet, dans le sens de la phase');
}

// ------------------------------------------------------------ 12. Fourier partiel, temps d'examen, plan des lignes
{
  ok(near(I.scanTime(I.PRESETS.t1, { N: 256 }), 128, 1e-9), 'T1 256 lignes, TR 500 ms : 2 min 08 s');
  ok(near(I.scanTime(I.PRESETS.t2, { N: 256 }), 64, 1e-9), 'T2 turbo ×16 : 256 × 4 s / 16 = 64 s');
  ok(near(I.scanTime(I.PRESETS.t1, { N: 256, nex: 2 }), 256, 1e-9), 'deux moyennages : temps doublé');
  ok(near(I.scanTime(I.PRESETS.t1, { N: 256, pf: 0.625 }), 80, 1e-9), 'Fourier partiel 5/8 : 160 lignes');
  for (const p of ['t2', 'pd', 'flair', 'stir']) {
    const seq = I.PRESETS[p], pl = I.lineplan(seq, I.defaults({ N: 128 }));
    ok(pl.echo[64] === I.effEcho(seq), p + ' : la ligne centrale est acquise à l\'écho effectif');
    const cnt = {}; pl.lines.forEach((m) => (cnt[pl.echo[m]] = (cnt[pl.echo[m]] || 0) + 1));
    ok(Object.keys(cnt).length === seq.ETL && Object.values(cnt).every((c) => c === 128 / seq.ETL), p + ' : chaque écho remplit le même nombre de lignes');
  }
  const seq = I.PRESETS.t1, ph = QA, base = { N: 128, fovx: 200, fovy: 200 };
  const full = mag(ph, seq, base).mag, zero = mag(ph, seq, Object.assign({ pf: 0.625, pfMode: 'zero' }, base)).mag, conj = mag(ph, seq, Object.assign({ pf: 0.625, pfMode: 'conj' }, base)).mag;
  const rms = (a) => Math.sqrt(a.reduce((s, v, i) => s + (v - full[i]) ** 2, 0) / a.length);
  ok(rms(conj) < 0.5 * rms(zero), 'Fourier partiel : la symétrie hermitienne fait mieux que les zéros (écho de spin)');
  const g = { kind: 'GRE', TR: 300, TE: 10, flip: 30, ETL: 1 }, gfull = mag(ph, g, base).mag, gconj = mag(ph, g, Object.assign({ pf: 0.625, pfMode: 'conj' }, base)).mag;
  const fatMask = I.roiMask(ph.ellipses, mag(ph, g, base).img, 3, 4), errFat = I.stats(gconj.map((v, i) => Math.abs(v - gfull[i])), fatMask).mean;
  const se = mag(ph, seq, base), seconj = mag(ph, seq, Object.assign({ pf: 0.625, pfMode: 'conj' }, base)).mag, errSe = I.stats(seconj.map((v, i) => Math.abs(v - se.mag[i])), I.roiMask(ph.ellipses, se.img, 3, 4)).mean;
  ok(errFat > 3 * errSe, 'Fourier partiel : la symétrie se trompe sur la graisse quand sa phase tourne (écho de gradient)');
}

// ------------------------------------------------------------ 13. métal
{
  const f = I.metalField(I.withMetal(BRAIN, true), 'ti');
  const c = I.METAL_CLIP;
  ok(near(f(c.x + c.a + 1e-9, c.y), (-I.GAMMA_BAR * 1e6 * 1.5 * 190e-6) / 3, 1e-3), 'sphère en titane : Δf = −γ̄·B0·Δχ/3 au bord (≈ −4 kHz)');
  ok(rel(f(c.x + 2 * c.a, c.y) / f(c.x + c.a + 1e-9, c.y), 1 / 8, 1e-6), 'champ de la sphère en 1/r³');
  // moteur sur grille contre moteur exact, sans métal
  for (const p of ['t1', 't2', 'gre']) {
    const seq = I.PRESETS[p], acq = { N: 128, r2p: p === 'gre' ? 10 : 0 };
    const a = mag(BRAIN, seq, acq).mag, b = mag(BRAIN, seq, acq, { grid: true }).mag;
    let num = 0, den = 0; for (let i = 0; i < a.length; i++) { num += (a[i] - b[i]) ** 2; den += a[i] ** 2; }
    ok(Math.sqrt(num / den) < 0.06, p + ' : grille fine = calcul exact à ' + (Math.sqrt(num / den) * 100).toFixed(1) + ' % près');
  }
  const voidArea = (p, metal) => {
    const seq = I.PRESETS[p], acq = { N: 128, metal }, r = mag(BRAIN, seq, acq), ref = mag(BRAIN, seq, { N: 128 }, { grid: true });
    let k = 0; for (let i = 0; i < r.mag.length; i++) { const [x, y] = I.pixelXY(r.img, Math.floor(i / 128), i % 128); if (Math.hypot(x - c.x, y - c.y) < 25 && r.mag[i] < 0.4 * ref.mag[i]) k++; }
    return k;
  };
  const se = voidArea('t1', 'ti'), gre = voidArea('gre', 'ti'), steel = voidArea('t1', 'steel');
  ok(gre > 1.5 * se, 'métal : trou plus grand en écho de gradient (' + gre + ' px) qu\'en écho de spin (' + se + ' px)');
  ok(steel > 2 * se, 'acier inoxydable : artefact bien plus grand que le titane (' + steel + ' / ' + se + ' px)');
  const distort = (bw) => { const seq = I.PRESETS.t1, r = mag(BRAIN, seq, { N: 128, metal: 'ti', bw }), ref = mag(BRAIN, seq, { N: 128, bw }, { grid: true }); let k = 0; for (let i = 0; i < r.mag.length; i++) { const [x, y] = I.pixelXY(r.img, Math.floor(i / 128), i % 128); if (Math.hypot(x - c.x, y - c.y) < 30 && Math.abs(r.mag[i] - ref.mag[i]) > 0.3 * ref.mag[i]) k++; } return k; };
  const d130 = distort(130), d520 = distort(520);
  ok(d520 < d130, 'métal : une bande passante plus large réduit la distorsion (' + d130 + ' → ' + d520 + ' px)');
}

// ------------------------------------------------------------ 14. séquences types valides
for (const [p, seq] of Object.entries(I.PRESETS)) ok(I.seqValid(seq), p + ' : tous les échos tiennent dans le TR');

console.log(n - ko + ' / ' + n + ' tests réussis');
if (ko) { fails.forEach((f) => console.log('ÉCHEC : ' + f)); process.exit(1); }
