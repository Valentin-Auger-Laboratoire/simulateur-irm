/* Simulateur d'IRM : moteur physique, sans dépendance.
   - Contraste : équations de Bloch résolues exactement (rotations des impulsions, précession, relaxation T1/T2)
     jusqu'au régime établi, pour l'écho de spin, l'écho de spin rapide (turbo), l'écho de gradient « spoilé »
     et l'inversion-récupération (FLAIR, STIR).
   - Image : le patient virtuel est fait d'ellipses dont la transformée de Fourier est exacte (fonction de Bessel J1) :
     l'espace k est rempli avec le signal exact, ligne par ligne, puis l'image est reconstruite par FFT.
     Déplacement chimique de la graisse, décroissance T2* pendant la lecture, repliement, sur-échantillonnage
     en lecture, Fourier partiel, mouvement, bruit thermique, moyennages.
   - Métal : champ d'une sphère aimantée calculé exactement dans le plan ; le signal est alors intégré sur une grille fine
     (distorsion en lecture, déphasage intravoxel en écho de gradient).
   Outil pédagogique : ce n'est pas un logiciel d'imagerie médicale.
   Utilisable dans le navigateur (window.IRM) et avec Node (module.exports). */
(function (root, factory) {
  const M = factory();
  if (typeof module === 'object' && module.exports) module.exports = M;
  else root.IRM = M;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  const DEG = Math.PI / 180, TAU = 2 * Math.PI;
  const GAMMA_BAR = 42.577; // MHz/T, rapport gyromagnétique du proton / 2π
  const B0 = 1.5; // T
  const FAT_PPM = -3.5; // pic principal (CH2) de la graisse par rapport à l'eau
  const larmor = (b0) => GAMMA_BAR * (b0 || B0); // MHz
  const fatHz = (b0) => FAT_PPM * GAMMA_BAR * (b0 || B0); // Hz (≈ −224 Hz à 1,5 T)

  // ---------------------------------------------------------------- tissus (1,5 T)
  // T1, T2 en ms ; PD = densité de protons relative ; r2p = R2' propre au tissu (s⁻¹, en plus de l'inhomogénéité du champ) ;
  // ppm = déplacement chimique ; hu = densité au scanner. src : origine des valeurs.
  const TISSUES = {
    air: { T1: 1, T2: 1, PD: 0, r2p: 0, ppm: 0, hu: -1000, src: '—' },
    wm: { T1: 884, T2: 72, PD: 0.69, r2p: 0, ppm: 0, hu: 28, src: 'Stanisz 2005' },
    gm: { T1: 1124, T2: 95, PD: 0.8, r2p: 0, ppm: 0, hu: 38, src: 'Stanisz 2005' },
    csf: { T1: 4000, T2: 2000, PD: 1, r2p: 0, ppm: 0, hu: 5, src: 'usuelle' },
    fat: { T1: 343, T2: 58, PD: 0.9, r2p: 0, ppm: FAT_PPM, hu: -100, src: 'Gold 2004' },
    muscle: { T1: 1008, T2: 44, PD: 0.72, r2p: 0, ppm: 0, hu: 45, src: 'Stanisz 2005' },
    blood: { T1: 1441, T2: 290, PD: 0.85, r2p: 0, ppm: 0, hu: 40, src: 'Stanisz 2005' },
    bone: { T1: 250, T2: 0.5, PD: 0.05, r2p: 0, ppm: 0, hu: 1000, src: 'indicative' },
    edema: { T1: 1300, T2: 130, PD: 0.84, r2p: 0, ppm: 0, hu: 22, src: 'illustrative' },
    hematoma: { T1: 950, T2: 50, PD: 0.8, r2p: 60, ppm: 0, hu: 70, src: 'illustrative' },
    metal: { T1: 1, T2: 1, PD: 0, r2p: 0, ppm: 0, hu: 30000, src: '—' },
  };
  const TISSUE_IDS = Object.keys(TISSUES);

  // ---------------------------------------------------------------- patients virtuels (mm)
  // ellipse : { x, y, a, b, ang (degrés), t : tissu, p : indice de l'ellipse qui la contient (null = air) }
  // Règle : chaque ellipse est entièrement dans son parent, deux ellipses de même parent ne se touchent pas.
  function E(x, y, a, b, ang, t, p) { return { x, y, a, b, ang, t, p }; }
  const BRAIN = [
    E(0, 0, 92, 114, 0, 'fat', null), // 0 cuir chevelu (graisse sous-cutanée)
    E(0, 0, 88, 110, 0, 'bone', 0), // 1 table externe du crâne
    E(0, 0, 85.5, 107.5, 0, 'fat', 1), // 2 diploé (moelle grasse)
    E(0, 0, 83.5, 105.5, 0, 'bone', 2), // 3 table interne
    E(0, 0, 81.5, 103.5, 0, 'csf', 3), // 4 liquide céphalo-rachidien des sillons
    E(0, -1, 79, 101, 0, 'gm', 4), // 5 cortex (substance grise)
    E(0, -3, 66, 87, 0, 'wm', 5), // 6 substance blanche
    E(-9, 9, 6, 23, 12, 'csf', 6), // 7 ventricule latéral gauche
    E(9, 9, 6, 23, -12, 'csf', 6), // 8 ventricule latéral droit
    E(-24, -24, 8, 10, 0, 'gm', 6), // 9 noyaux gris
    E(24, -24, 8, 10, 0, 'gm', 6), // 10 noyaux gris
    E(32, 44, 12, 9, 20, 'edema', 6), // 11 œdème (AVC de quelques heures)
    E(-34, 44, 10, 7, -25, 'hematoma', 6), // 12 hématome
  ];
  const METAL_CLIP = E(22, -54, 3, 3, 0, 'metal', 6); // clip en métal (sphère de 3 mm de rayon)
  const QA_ORDER = ['wm', 'gm', 'csf', 'fat', 'muscle', 'blood', 'edema', 'hematoma', 'bone'];
  const QA = QA_ORDER.map((t, i) => E(((i % 3) - 1) * 52, (1 - Math.floor(i / 3)) * 52, 17, 17, 0, t, null));
  const PHANTOMS = {
    brain: { fov: 240, ellipses: BRAIN, metalAt: 6, extent: 115 },
    qa: { fov: 200, ellipses: QA, metalAt: null, extent: 72 },
  };
  function withMetal(ph, on) {
    if (!on || ph.metalAt === null) return ph.ellipses;
    return ph.ellipses.concat([Object.assign({}, METAL_CLIP, { p: ph.metalAt })]);
  }
  function inEllipse(e, x, y) {
    const c = Math.cos(e.ang * DEG), s = Math.sin(e.ang * DEG), dx = x - e.x, dy = y - e.y;
    const u = (dx * c + dy * s) / e.a, w = (-dx * s + dy * c) / e.b;
    return u * u + w * w <= 1;
  }
  // tissu au point (x, y) : l'ellipse la plus profonde qui contient le point (les enfants sont après leur parent)
  function tissueAt(ells, x, y) {
    let t = 'air';
    for (let i = 0; i < ells.length; i++) if (inEllipse(ells[i], x, y)) t = ells[i].t;
    return t;
  }
  // indice de la région (ellipse la plus profonde) au point
  function regionAt(ells, x, y) {
    let r = -1;
    for (let i = 0; i < ells.length; i++) if (inEllipse(ells[i], x, y)) r = i;
    return r;
  }

  // ---------------------------------------------------------------- Bloch exact (un isochromat à la résonance)
  function rot(M, angDeg, axis) {
    const a = angDeg * DEG, c = Math.cos(a), s = Math.sin(a);
    if (axis === 'x') { const y = M[1] * c + M[2] * s, z = -M[1] * s + M[2] * c; M[1] = y; M[2] = z; }
    else { const x = M[0] * c - M[2] * s, z = M[0] * s + M[2] * c; M[0] = x; M[2] = z; }
  }
  function relax(M, dt, T1, T2, M0, df) {
    if (dt <= 0) return;
    const e2 = Math.exp(-dt / T2), e1 = Math.exp(-dt / T1);
    if (df) { const ph = (TAU * df * dt) / 1000, c = Math.cos(ph), s = Math.sin(ph), x = M[0] * c - M[1] * s, y = M[0] * s + M[1] * c; M[0] = x; M[1] = y; }
    M[0] *= e2; M[1] *= e2; M[2] = M0 + (M[2] - M0) * e1;
  }
  // séquence : { kind: 'SE' | 'GRE' | 'IR', TR, TE, TI, flip, ETL, ESP }
  //   SE  : 90° puis 180° ; ETL > 1 = écho de spin rapide (train de 180°, espacement ESP, écho effectif proche de TE)
  //   GRE : impulsion α, écho de gradient à TE, aimantation transversale détruite (« spoilée ») en fin de TR
  //   IR  : inversion 180°, puis lecture en écho de spin (ou rapide) au temps TI
  function echoTimes(seq) {
    if (seq.kind === 'GRE') return [seq.TE];
    const etl = Math.max(1, seq.ETL | 0);
    if (etl === 1) return [seq.TE];
    return Array.from({ length: etl }, (_, k) => (k + 1) * seq.ESP);
  }
  function effEcho(seq) { // indice (0…) de l'écho qui remplit le centre de l'espace k
    const te = echoTimes(seq);
    if (te.length === 1) return 0;
    let best = 0; te.forEach((t, k) => { if (Math.abs(t - seq.TE) < Math.abs(te[best] - seq.TE)) best = k; });
    return best;
  }
  function events(seq) {
    const ev = [], t0 = seq.kind === 'IR' ? seq.TI : 0;
    if (seq.kind === 'IR') ev.push({ t: 0, kind: 'rf', ang: 180, axis: 'x' });
    if (seq.kind === 'GRE') { ev.push({ t: 0, kind: 'rf', ang: seq.flip, axis: 'x' }); ev.push({ t: seq.TE, kind: 'echo', e: 0 }); }
    else {
      ev.push({ t: t0, kind: 'rf', ang: 90, axis: 'x' });
      const etl = Math.max(1, seq.ETL | 0);
      if (etl === 1) { ev.push({ t: t0 + seq.TE / 2, kind: 'rf', ang: 180, axis: 'y' }); ev.push({ t: t0 + seq.TE, kind: 'echo', e: 0 }); }
      else for (let k = 0; k < etl; k++) { ev.push({ t: t0 + (k + 0.5) * seq.ESP, kind: 'rf', ang: 180, axis: 'y' }); ev.push({ t: t0 + (k + 1) * seq.ESP, kind: 'echo', e: k }); }
    }
    return ev.sort((a, b) => a.t - b.t || (a.kind === 'rf' ? -1 : 1));
  }
  function seqValid(seq) {
    const ev = events(seq), last = ev[ev.length - 1].t;
    return last < seq.TR;
  }
  // régime établi : on répète le TR jusqu'à ce que Mz avant la première impulsion ne bouge plus
  // renvoie { mz : Mz avant la séquence, echoes : signal signé (composante y) à chaque écho, sans R2' }
  function steadyState(tis, seq) {
    const T = typeof tis === 'string' ? TISSUES[tis] : tis;
    if (T.PD === 0) return { mz: 0, echoes: echoTimes(seq).map(() => 0), iters: 0 };
    const ev = events(seq), M = [0, 0, T.PD], echoes = [];
    let mzPrev = Infinity, it = 0;
    for (it = 0; it < 2000; it++) {
      M[0] = 0; M[1] = 0; // aimantation transversale détruite en fin de TR
      const mzStart = M[2];
      if (Math.abs(mzStart - mzPrev) < 1e-13 * T.PD) break;
      mzPrev = mzStart;
      let t = 0;
      for (const e of ev) {
        relax(M, e.t - t, T.T1, T.T2, T.PD); t = e.t;
        if (e.kind === 'rf') rot(M, e.ang, e.axis); else echoes[e.e] = M[1];
      }
      relax(M, seq.TR - t, T.T1, T.T2, T.PD);
    }
    return { mz: mzPrev, echoes: echoes.slice(), iters: it };
  }
  // formules analytiques classiques (pour vérifier Bloch et pour l'affichage)
  const formula = {
    se: (T, TR, TE) => T.PD * (1 - 2 * Math.exp(-(TR - TE / 2) / T.T1) + Math.exp(-TR / T.T1)) * Math.exp(-TE / T.T2),
    gre: (T, TR, TE, flip) => { const E1 = Math.exp(-TR / T.T1), a = flip * DEG; return (T.PD * Math.sin(a) * (1 - E1)) / (1 - Math.cos(a) * E1) * Math.exp(-TE / T.T2); },
    ernst: (T, TR) => Math.acos(Math.exp(-TR / T.T1)) / DEG,
    // TI qui annule un tissu en inversion-récupération (lecture en écho de spin)
    tiNull: (T, TR, TE) => {
      // Mz(TI) = PD·[1 − 2e^(−TI/T1) + 2e^(−(TR−TE/2)/T1) − e^(−TR/T1)] = 0  →  e^(−TI/T1) = (1 + 2e^(−(TR−TE/2)/T1) − e^(−TR/T1)) / 2
      const q = (1 + 2 * Math.exp(-(TR - TE / 2) / T.T1) - Math.exp(-TR / T.T1)) / 2;
      return -T.T1 * Math.log(q);
    },
  };
 // inversion-récupération, lecture en écho de spin : Mz(TI) = PD·[1 − 2e^(−TI/T1) + 2e^(−(TR−TE/2)/T1) − e^(−TR/T1)]
  formula.ir = (T, TR, TI, TE) => T.PD * (1 - 2 * Math.exp(-TI / T.T1) + 2 * Math.exp(-(TR - TE / 2) / T.T1) - Math.exp(-TR / T.T1)) * Math.exp(-TE / T.T2);
  // TI d'annulation trouvé numériquement avec Bloch (vaut aussi pour la lecture en écho de spin rapide)
  function tiNullBloch(tis, seq) {
    const T = typeof tis === 'string' ? TISSUES[tis] : tis;
    const f = (ti) => steadyState(T, Object.assign({}, seq, { kind: 'IR', TI: ti })).echoes[effEcho(seq)];
    let lo = 1, hi = seq.TR - echoTimes(seq)[echoTimes(seq).length - 1] - 1;
    if (f(lo) * f(hi) > 0) return null;
    for (let k = 0; k < 80; k++) { const mid = (lo + hi) / 2; if (f(lo) * f(mid) <= 0) hi = mid; else lo = mid; }
    return (lo + hi) / 2;
  }

  // ---------------------------------------------------------------- ensemble d'isochromats (pour voir T2* et l'écho)
  // désaccords de fréquence répartis selon une loi de Lorentz (quantiles de Cauchy) : leur somme décroît en e^(−R2'·t)
  function lorentzOffsets(n, r2p) {
    const w = r2p / TAU; // demi-largeur en Hz
    return Float64Array.from({ length: n }, (_, i) => w * Math.tan(Math.PI * ((i + 0.5) / n - 0.5)));
  }
  // simule un TR (après le régime établi) : renvoie les instants, l'aimantation moyenne et celle de quelques isochromats
  function ensemble(tis, seq, opts) {
    const T = typeof tis === 'string' ? TISSUES[tis] : tis, o = Object.assign({ n: 400, show: 12, samples: 900, r2p: 10, tmax: seq.TR }, opts || {});
    const r2p = o.r2p + (T.r2p || 0), offs = lorentzOffsets(o.n, r2p);
    const ss = steadyState(T, seq), ev = events(seq);
    const spins = Array.from({ length: o.n }, () => [0, 0, ss.mz]);
    const times = [], mean = [], some = [];
    const pick = Array.from({ length: o.show }, (_, i) => Math.round(((i + 0.5) / o.show) * o.n - 0.5));
    // instants : réguliers + chaque événement (les échos tombent exactement sur un échantillon)
    const ts = new Set(); for (let k = 0; k <= o.samples; k++) ts.add((k / o.samples) * o.tmax);
    ev.forEach((e) => { if (e.t <= o.tmax) { ts.add(e.t); ts.add(Math.max(0, e.t - 1e-9)); } });
    const grid = [...ts].sort((a, b) => a - b);
    let t = 0, ei = 0;
    const record = () => {
      let mx = 0, my = 0, mz = 0; for (const s of spins) { mx += s[0]; my += s[1]; mz += s[2]; }
      times.push(t); mean.push([mx / o.n, my / o.n, mz / o.n]); some.push(pick.map((i) => spins[i].slice()));
    };
    for (const tg of grid) {
      while (ei < ev.length && ev[ei].t <= tg) {
        const e = ev[ei]; for (let i = 0; i < o.n; i++) relax(spins[i], e.t - t, T.T1, T.T2, T.PD, offs[i]); t = e.t;
        if (e.kind === 'rf') for (const s of spins) rot(s, e.ang, e.axis);
        ei++;
      }
      for (let i = 0; i < o.n; i++) relax(spins[i], tg - t, T.T1, T.T2, T.PD, offs[i]); t = tg;
      record();
    }
    return { times, mean, some, events: ev, mz0: ss.mz };
  }

  // ---------------------------------------------------------------- FFT (radix 2, en place)
  function fft(re, im, inv) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1; for (; j & bit; bit >>= 1) j ^= bit; j ^= bit;
      if (i < j) { let t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = ((inv ? 1 : -1) * TAU) / len, wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < len / 2; k++) {
          const a = i + k, b = a + len / 2, xr = re[b] * cr - im[b] * ci, xi = re[b] * ci + im[b] * cr;
          re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
          const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
        }
      }
    }
  }
  // transformée discrète directe (référence pour les tests)
  function dft(re, im, inv) {
    const n = re.length, R = new Float64Array(n), I = new Float64Array(n), sg = inv ? 1 : -1;
    for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) { const a = (sg * TAU * j * k) / n; R[k] += re[j] * Math.cos(a) - im[j] * Math.sin(a); I[k] += re[j] * Math.sin(a) + im[j] * Math.cos(a); }
    return [R, I];
  }

  // ---------------------------------------------------------------- Bessel J1 et transformée de Fourier d'une ellipse
  function besselJ1(x) { // approximations rationnelles classiques (Numerical Recipes), erreur ~1e-8
    const ax = Math.abs(x);
    if (ax < 8) {
      const y = x * x;
      const a1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.4826 + y * -30.16036606)))));
      const a2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))));
      return a1 / a2;
    }
    const z = 8 / ax, y = z * z, xx = ax - 2.356194491;
    const a1 = 1 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * -0.240337019e-6)));
    const a2 = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
    const r = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * a1 - z * Math.sin(xx) * a2);
    return x < 0 ? -r : r;
  }
  const jinc = (q) => { const u = TAU * q; return Math.abs(u) < 1e-8 ? 1 : (2 * besselJ1(u)) / u; };
  // ∫∫_ellipse e^(−i2π(kx·x + ky·y)) dx dy  (densité 1, k en mm⁻¹) → [re, im]
  function ellipseFT(e, kx, ky) {
    const c = Math.cos(e.ang * DEG), s = Math.sin(e.ang * DEG), ku = kx * c + ky * s, kv = -kx * s + ky * c;
    const amp = Math.PI * e.a * e.b * jinc(Math.hypot(e.a * ku, e.b * kv)), ph = -TAU * (kx * e.x + ky * e.y);
    return [amp * Math.cos(ph), amp * Math.sin(ph)];
  }

  // ---------------------------------------------------------------- acquisition
  // acq : { N, fovx, fovy, bw (Hz/pixel), thk (mm), nex, pf (fraction de lignes), pfMode ('zero'|'conj'), filter ('none'|'hamming'),
  //         r2p (s⁻¹, inhomogénéité), motion: { type: 'none'|'breath'|'jump', amp (mm), period (s) }, metal: null|'ti'|'steel',
  //         order: 'seq'|'centric', noise: bool, seed, sigmaRef }
  const OS = 2; // sur-échantillonnage en lecture (le récepteur filtre : pas de repliement gauche-droite)
  const METALS = { ti: { dchi: 190e-6 }, steel: { dchi: 3000e-6 } };
  const SIGMA_REF = 0.012;
  function defaults(acq) {
    return Object.assign({ N: 128, fovx: 240, fovy: 240, bw: 130, thk: 5, nex: 1, pf: 1, pfMode: 'zero', filter: 'none', r2p: 10, motion: { type: 'none', amp: 0, period: 4 }, metal: null, order: 'seq', noise: true, seed: 7, sigmaRef: SIGMA_REF }, acq || {});
  }
  // lignes acquises, écho de chaque ligne, instant d'acquisition de chaque ligne
  function lineplan(seq, acq) {
    const N = acq.N, nAcq = Math.min(N, Math.max(Math.ceil(N / 2) + 4, Math.round(acq.pf * N))), m0 = N - nAcq;
    const lines = []; for (let m = m0; m < N; m++) lines.push(m);
    const ky = (m) => m - N / 2;
    let ordered = lines.slice();
    if (acq.order === 'centric') ordered.sort((a, b) => Math.abs(ky(a)) - Math.abs(ky(b)) || ky(a) - ky(b));
    const te = echoTimes(seq), etl = te.length, eEff = effEcho(seq);
    const echo = new Int32Array(N).fill(-1), time = new Float64Array(N).fill(NaN);
    const t0 = seq.kind === 'IR' ? seq.TI : 0;
    if (etl === 1) ordered.forEach((m, i) => { echo[m] = 0; time[m] = i * seq.TR + t0 + te[0]; });
    else {
      // k-space découpé en ETL segments contigus ; le segment le plus proche du centre reçoit l'écho effectif
      const L = Math.ceil(nAcq / etl), segs = [];
      for (let s = 0; s < etl; s++) { const ls = lines.slice(s * L, (s + 1) * L); if (ls.length) segs.push(ls); }
      const segOrder = segs.map((ls, i) => [i, Math.abs(ls.reduce((a, m) => a + ky(m), 0) / ls.length)]).sort((a, b) => a[1] - b[1]).map((x) => x[0]);
      const echoOrder = Array.from({ length: etl }, (_, k) => k).sort((a, b) => Math.abs(a - eEff) - Math.abs(b - eEff) || b - a);
      segOrder.forEach((si, r) => segs[si].forEach((m, j) => { echo[m] = echoOrder[r]; time[m] = j * seq.TR + t0 + te[echoOrder[r]]; }));
    }
    const shots = Math.ceil(nAcq / etl);
    return { lines, nAcq, m0, echo, time, shots, etl, tAcq: shots * seq.TR, eEff };
  }
  function scanTime(seq, acq) { const p = lineplan(seq, defaults(acq)); return (p.shots * seq.TR * defaults(acq).nex) / 1000; } // s
  // signal complexe d'un tissu à l'instant t (ms depuis l'impulsion de lecture), pour l'écho e
  function tissueTables(seq, acq, tissues) {
    const te = echoTimes(seq), tab = {};
    for (const id of tissues) {
      const T = TISSUES[id], ss = steadyState(T, seq);
      tab[id] = { B: ss.echoes, f: T.ppm ? (T.ppm / FAT_PPM) * fatHz() : 0, r2p: acq.r2p + T.r2p, T2: T.T2, te };
    }
    return tab;
  }
  // valeur du signal (complexe) : B·e^(−(t−TE)/T2)·e^(−R2'·|t − tr|)·e^(i2πf(t − tr)) ; tr = TE (écho de spin) ou 0 (écho de gradient)
  function sigAt(tt, e, t, gre) {
    const TE = tt.te[e], tr = gre ? 0 : TE, B = tt.B[e];
    if (!B) return [0, 0];
    const amp = B * Math.exp(-(t - TE) / tt.T2) * Math.exp((-tt.r2p * Math.abs(t - tr)) / 1000), ph = (TAU * tt.f * (t - tr)) / 1000;
    return [amp * Math.cos(ph), amp * Math.sin(ph)];
  }
  function rng(seed) { let s = seed >>> 0 || 1; return () => ((s = (Math.imul(s, 1664525) + 1013904223) >>> 0) / 4294967296); }
  function gauss(r) { let u = r(); if (u < 1e-300) u = 1e-300; return Math.sqrt(-2 * Math.log(u)) * Math.cos(TAU * r()); }
  function noiseSigma(acq) { // écart-type du bruit (par voie, réelle ou imaginaire) dans l'image, acquisition complète
    const bwTot = acq.bw * acq.N;
    return (acq.sigmaRef * Math.sqrt((acq.N * acq.N * bwTot) / acq.nex)) / (acq.fovx * acq.fovy) * (5 / acq.thk);
  }
  // déplacement du patient (mm, selon y) à l'instant t (ms) de l'examen
  function motionAt(mo, t, tTot) {
    if (!mo || mo.type === 'none' || !mo.amp) return 0;
    if (mo.type === 'breath') return mo.amp * Math.sin((TAU * t) / (mo.period * 1000));
    if (mo.type === 'jump') return t > tTot / 2 ? mo.amp : 0;
    return 0;
  }
  function metalField(ells, metal) { // Δf (Hz) dans le plan de la sphère : −γ̄·B0·Δχ/3·(a/r)³, B0 perpendiculaire à la coupe
    if (!metal) return null;
    const clip = ells.find((e) => e.t === 'metal'); if (!clip) return null;
    const k = (-GAMMA_BAR * 1e6 * B0 * METALS[metal].dchi) / 3, a = clip.a;
    return (x, y) => { const r = Math.hypot(x - clip.x, y - clip.y); return r <= a ? 0 : k * (a / r) ** 3; };
  }
  // cache des transformées de Fourier des ellipses (ne dépendent que de la géométrie)
  const ftCache = new Map();
  function ellipseFTs(ells, acq) {
    const key = JSON.stringify([ells.map((e) => [e.x, e.y, e.a, e.b, e.ang]), acq.N, acq.fovx, acq.fovy]);
    if (ftCache.has(key)) return ftCache.get(key);
    const N = acq.N, Nro = OS * N, out = ells.map(() => ({ re: new Float32Array(Nro * N), im: new Float32Array(Nro * N) }));
    for (let m = 0; m < N; m++) {
      const ky = (m - N / 2) / acq.fovy;
      for (let n = 0; n < Nro; n++) {
        const kx = (n - Nro / 2) / (OS * acq.fovx), idx = m * Nro + n;
        ells.forEach((e, i) => { const [r, im] = ellipseFT(e, kx, ky); out[i].re[idx] = r; out[i].im[idx] = im; });
      }
    }
    if (ftCache.size > 6) ftCache.delete(ftCache.keys().next().value);
    ftCache.set(key, out);
    return out;
  }
  // moteur exact (ellipses) : renvoie l'espace k { re, im, Nro, N, plan, ... } sans bruit ni mouvement
  function kspaceAnalytic(ph, seq, acq) {
    const ells = withMetal(ph, false), N = acq.N, Nro = OS * N, plan = lineplan(seq, acq), gre = seq.kind === 'GRE';
    const tabs = tissueTables(seq, acq, [...new Set(ells.map((e) => e.t))]);
    const fts = ellipseFTs(ells, acq), dt = 1000 / (OS * acq.bw * N); // ms entre deux échantillons
    const re = new Float64Array(Nro * N), im = new Float64Array(Nro * N);
    // table des poids complexes par ellipse, écho et échantillon
    const ne = echoTimes(seq).length, W = ells.map(() => Array.from({ length: ne }, () => null));
    const weight = (i, e) => {
      if (W[i][e]) return W[i][e];
      const wr = new Float64Array(Nro), wi = new Float64Array(Nro), el = ells[i], par = el.p === null ? null : ells[el.p];
      for (let n = 0; n < Nro; n++) {
        const t = echoTimes(seq)[e] + (n - Nro / 2) * dt;
        const [a, b] = sigAt(tabs[el.t], e, t, gre), [c, d] = par ? sigAt(tabs[par.t], e, t, gre) : [0, 0];
        wr[n] = a - c; wi[n] = b - d;
      }
      return (W[i][e] = [wr, wi]);
    };
    for (const m of plan.lines) {
      const e = plan.echo[m];
      for (let i = 0; i < ells.length; i++) {
        const [wr, wi] = weight(i, e), F = fts[i];
        for (let n = 0; n < Nro; n++) { const idx = m * Nro + n, fr = F.re[idx], fi = F.im[idx]; re[idx] += wr[n] * fr - wi[n] * fi; im[idx] += wr[n] * fi + wi[n] * fr; }
      }
    }
    return { re, im, Nro, N, plan, engine: 'exact' };
  }
  // moteur sur grille fine (métal) : chaque petit élément est déplacé en lecture de Δf/(bande passante par pixel) pixels
  // et garde la phase 2πΔf·(TE − tr) ; la direction de phase est calculée par somme directe (repliement exact)
  const labelCache = new Map();
  function kspaceGrid(ph, seq, acq, opts) {
    const o = Object.assign({ subx: 4, suby: 3 }, opts || {});
    const ells = withMetal(ph, !!acq.metal), N = acq.N, Nro = OS * N, plan = lineplan(seq, acq), gre = seq.kind === 'GRE';
    const tabs = tissueTables(seq, acq, [...new Set(ells.map((e) => e.t))]), field = metalField(ells, acq.metal);
    const Lx = o.subx * Nro, hx = (OS * acq.fovx) / Lx, hy = acq.fovy / N / o.suby, Y = ph.extent + 2;
    const R = Math.ceil((2 * Y) / hy), ys = Float64Array.from({ length: R }, (_, r) => -Y + (r + 0.5) * hy);
    const xs = Float64Array.from({ length: Lx }, (_, j) => (j - Lx / 2 + 0.5) * hx);
    const key = JSON.stringify([ells.map((e) => [e.x, e.y, e.a, e.b, e.ang, e.t]), Lx, hx, R, hy]);
    let lab = labelCache.get(key);
    if (!lab) {
      lab = new Int8Array(R * Lx).fill(-1);
      for (let r = 0; r < R; r++) for (let j = 0; j < Lx; j++) { if (Math.abs(xs[j]) > ph.extent + 2) continue; const t = tissueAt(ells, xs[j], ys[r]); if (t !== 'air' && TISSUES[t].PD > 0) lab[r * Lx + j] = TISSUE_IDS.indexOf(t); }
      if (labelCache.size > 4) labelCache.delete(labelCache.keys().next().value);
      labelCache.set(key, lab);
    }
    const bwTot = acq.bw * N, area = hx * hy, ne = echoTimes(seq).length;
    const used = [...new Set(plan.lines.map((m) => plan.echo[m]))];
    // G[e][r] : transformée en lecture de la ligne r du patient, pour l'écho e (seulement les fréquences utiles)
    const G = {}; const rr = new Float64Array(Lx), ri = new Float64Array(Lx);
    const corr = Float64Array.from({ length: Nro }, (_, n) => { const q = (n - Nro / 2) / Lx; const s = Math.abs(q) < 1e-12 ? 1 : Math.sin(Math.PI * q) / (Math.PI * q); return 1 / (s * s); });
    for (const e of used) {
      const TE = echoTimes(seq)[e], tr = gre ? 0 : TE, gr = new Float32Array(R * Nro), gi = new Float32Array(R * Nro), rowUsed = new Uint8Array(R);
      for (let r = 0; r < R; r++) {
        rr.fill(0); ri.fill(0); let any = false;
        for (let j = 0; j < Lx; j++) {
          const li = lab[r * Lx + j]; if (li < 0) continue;
          const tt = tabs[TISSUE_IDS[li]], B = tt.B[e]; if (!B) continue;
          const df = tt.f + (field ? field(xs[j], ys[r]) : 0);
          const amp = B * Math.exp((-tt.r2p * Math.abs(TE - tr)) / 1000) * area, ph0 = (TAU * df * (TE - tr)) / 1000;
          const xp = xs[j] - (df * acq.fovx) / bwTot; // position apparente en lecture
          let u = xp / hx + Lx / 2 - 0.5; const j0 = Math.floor(u), f = u - j0;
          const c = Math.cos(ph0) * amp, s = Math.sin(ph0) * amp;
          const a0 = ((j0 % Lx) + Lx) % Lx, a1 = (((j0 + 1) % Lx) + Lx) % Lx;
          rr[a0] += c * (1 - f); ri[a0] += s * (1 - f); rr[a1] += c * f; ri[a1] += s * f; any = true;
        }
        if (!any) continue;
        rowUsed[r] = 1;
        fft(rr, ri, false);
        for (let n = 0; n < Nro; n++) {
          const q = n - Nro / 2, qi = ((q % Lx) + Lx) % Lx;
          // centre des éléments en x_j = (j − Lx/2 + 0,5)·hx : facteur e^(−i2π q (−Lx/2 + 0,5)/Lx)
          const ph1 = -TAU * q * (-Lx / 2 + 0.5) / Lx, c1 = Math.cos(ph1), s1 = Math.sin(ph1);
          gr[r * Nro + n] = (rr[qi] * c1 - ri[qi] * s1) * corr[n]; gi[r * Nro + n] = (rr[qi] * s1 + ri[qi] * c1) * corr[n];
        }
      }
      G[e] = { gr, gi, rowUsed };
    }
    const re = new Float64Array(Nro * N), im = new Float64Array(Nro * N);
    for (const m of plan.lines) {
      const ky = (m - N / 2) / acq.fovy, g = G[plan.echo[m]];
      for (let r = 0; r < R; r++) {
        if (!g.rowUsed[r]) continue;
        const a = -TAU * ky * ys[r], c = Math.cos(a), s = Math.sin(a), base = r * Nro;
        for (let n = 0; n < Nro; n++) { const xr = g.gr[base + n], xi = g.gi[base + n], idx = m * Nro + n; re[idx] += xr * c - xi * s; im[idx] += xr * s + xi * c; }
      }
    }
    return { re, im, Nro, N, plan, engine: 'grille' };
  }
  // acquisition complète : signal (exact ou grille), mouvement, moyennages, bruit
  function acquire(ph, seq, acqIn, opts) {
    const acq = defaults(acqIn), useGrid = !!acq.metal || (opts && opts.grid);
    const k = useGrid ? kspaceGrid(ph, seq, acq, opts) : kspaceAnalytic(ph, seq, acq);
    const { re, im, Nro, N, plan } = k;
    // mouvement : translation rigide selon y → phase e^(−i2π·ky·dy(t)) sur chaque ligne, moyennée sur les répétitions
    if (acq.motion && acq.motion.type !== 'none' && acq.motion.amp) {
      const tTot = plan.tAcq * acq.nex;
      for (const m of plan.lines) {
        const ky = (m - N / 2) / acq.fovy; let cr = 0, ci = 0;
        for (let rep = 0; rep < acq.nex; rep++) { const dy = motionAt(acq.motion, plan.time[m] + rep * plan.tAcq, tTot), a = -TAU * ky * dy; cr += Math.cos(a); ci += Math.sin(a); }
        cr /= acq.nex; ci /= acq.nex;
        for (let n = 0; n < Nro; n++) { const idx = m * Nro + n, a = re[idx], b = im[idx]; re[idx] = a * cr - b * ci; im[idx] = a * ci + b * cr; }
      }
    }
    const clean = { re: Float64Array.from(re), im: Float64Array.from(im) };
    let sigma = 0;
    if (acq.noise) {
      sigma = noiseSigma(acq);
      const dk = 1 / (OS * acq.fovx * acq.fovy), sk = sigma / (dk * Math.sqrt(Nro * N)), r = rng(acq.seed);
      for (const m of plan.lines) for (let n = 0; n < Nro; n++) { const idx = m * Nro + n; re[idx] += sk * gauss(r); im[idx] += sk * gauss(r); }
    }
    return { re, im, clean, Nro, N, plan, acq, seq, sigma, engine: k.engine };
  }
  // reconstruction : Fourier partiel, filtre, FFT 2D inverse, retrait du sur-échantillonnage
  function reconstruct(K, which) {
    const { Nro, N, plan, acq } = K, src = which === 'clean' ? K.clean : K;
    const re = Float64Array.from(src.re), im = Float64Array.from(src.im), acquired = new Uint8Array(N);
    plan.lines.forEach((m) => (acquired[m] = 1));
    if (acq.pf < 1 && acq.pfMode === 'conj') {
      for (let m = 0; m < N; m++) {
        if (acquired[m]) continue; const m2 = (N - m) % N; if (!acquired[m2]) continue;
        for (let n = 0; n < Nro; n++) { const n2 = (Nro - n) % Nro; re[m * Nro + n] = re[m2 * Nro + n2]; im[m * Nro + n] = -im[m2 * Nro + n2]; }
      }
    }
    if (acq.filter === 'hamming') {
      for (let m = 0; m < N; m++) { const wy = 0.54 + 0.46 * Math.cos((Math.PI * (m - N / 2)) / (N / 2)); for (let n = 0; n < Nro; n++) { const wx = 0.54 + 0.46 * Math.cos((Math.PI * (n - Nro / 2)) / (Nro / 2)), w = wx * wy; re[m * Nro + n] *= w; im[m * Nro + n] *= w; } }
    }
    const scale = 1 / (OS * acq.fovx * acq.fovy);
    // lecture : FFT inverse de chaque ligne, puis on garde les N pixels du centre
    const mid = new Float64Array(N * N), midI = new Float64Array(N * N), lr = new Float64Array(Nro), li = new Float64Array(Nro);
    for (let m = 0; m < N; m++) {
      for (let n = 0; n < Nro; n++) { const sg = n & 1 ? -1 : 1; lr[n] = re[m * Nro + n] * sg; li[n] = im[m * Nro + n] * sg; }
      fft(lr, li, true);
      for (let j = 0; j < N; j++) { const jj = j + Nro / 2 - N / 2, sg = jj & 1 ? -1 : 1; mid[m * N + j] = lr[jj] * sg; midI[m * N + j] = li[jj] * sg; }
    }
    // phase : FFT inverse de chaque colonne ; ligne 0 de l'image = y le plus grand (haut de l'écran)
    const outR = new Float64Array(N * N), outI = new Float64Array(N * N), cr = new Float64Array(N), ci = new Float64Array(N);
    for (let j = 0; j < N; j++) {
      for (let m = 0; m < N; m++) { const sg = m & 1 ? -1 : 1; cr[m] = mid[m * N + j] * sg; ci[m] = midI[m * N + j] * sg; }
      fft(cr, ci, true);
      for (let i = 0; i < N; i++) { const sg = i & 1 ? -1 : 1, row = N - 1 - i; outR[row * N + j] = cr[i] * sg * scale; outI[row * N + j] = ci[i] * sg * scale; }
    }
    // pixel (row, j) ↔ x = (j − N/2)·Δx, y = (N − 1 − row − N/2)·Δy
    return { re: outR, im: outI, N, dx: acq.fovx / N, dy: acq.fovy / N };
  }
  const magnitude = (img) => Float64Array.from(img.re, (r, i) => Math.hypot(r, img.im[i]));
  const pixelXY = (img, row, j) => [(j - img.N / 2) * img.dx, (img.N / 2 - 1 - row) * img.dy];
  // ROI : pixels dont le centre est dans la région, à plus de `margin` mm de tout bord
  function roiMask(ells, img, region, margin) {
    const N = img.N, mask = new Uint8Array(N * N), d = margin, dirs = [[0, 0], [d, 0], [-d, 0], [0, d], [0, -d], [d * 0.7, d * 0.7], [-d * 0.7, d * 0.7], [d * 0.7, -d * 0.7], [-d * 0.7, -d * 0.7]];
    for (let row = 0; row < N; row++) for (let j = 0; j < N; j++) {
      const [x, y] = pixelXY(img, row, j);
      if (dirs.every(([a, b]) => regionAt(ells, x + a, y + b) === region)) mask[row * N + j] = 1;
    }
    return mask;
  }
  function stats(arr, mask) {
    let s = 0, s2 = 0, n = 0;
    for (let i = 0; i < arr.length; i++) if (!mask || mask[i]) { s += arr[i]; s2 += arr[i] * arr[i]; n++; }
    const mean = n ? s / n : 0; return { mean, sd: n > 1 ? Math.sqrt(Math.max(0, s2 / n - mean * mean) * (n / (n - 1))) : 0, n };
  }
  // signal théorique d'un tissu au centre de l'espace k (écho effectif, R2' compris en écho de gradient)
  function theory(id, seq, acq) {
    const a = defaults(acq), T = TISSUES[id], ss = steadyState(T, seq), e = effEcho(seq), TE = echoTimes(seq)[e];
    const B = ss.echoes[e] || 0, r2p = a.r2p + T.r2p;
    return Math.abs(B) * (seq.kind === 'GRE' ? Math.exp((-r2p * TE) / 1000) : 1);
  }
  function predictedSNR(id, seq, acq) { const a = defaults(acq); return theory(id, seq, a) / noiseSigma(a); }

  // séquences types (paramètres de départ réalistes à 1,5 T)
  const PRESETS = {
    t1: { kind: 'SE', TR: 500, TE: 15, ETL: 1, ESP: 10, TI: 0, flip: 90 },
    t1ir: { kind: 'IR', TR: 2000, TE: 12, ETL: 1, ESP: 10, TI: 800, flip: 90 },
    t2: { kind: 'SE', TR: 4000, TE: 100, ETL: 16, ESP: 10, TI: 0, flip: 90 },
    pd: { kind: 'SE', TR: 3000, TE: 20, ETL: 8, ESP: 10, TI: 0, flip: 90 },
    flair: { kind: 'IR', TR: 9000, TE: 120, ETL: 16, ESP: 10, TI: 2356, flip: 90 },
    stir: { kind: 'IR', TR: 4000, TE: 40, ETL: 8, ESP: 10, TI: 238, flip: 90 },
    gre: { kind: 'GRE', TR: 600, TE: 20, ETL: 1, ESP: 10, TI: 0, flip: 20 },
    inphase: { kind: 'GRE', TR: 150, TE: 4.47, ETL: 1, ESP: 10, TI: 0, flip: 70 },
    opposed: { kind: 'GRE', TR: 150, TE: 2.24, ETL: 1, ESP: 10, TI: 0, flip: 70 },
  };

  return {
    GAMMA_BAR, B0, FAT_PPM, larmor, fatHz, TISSUES, TISSUE_IDS, PHANTOMS, METAL_CLIP, METALS, QA_ORDER, withMetal, inEllipse, tissueAt, regionAt,
    rot, relax, events, echoTimes, effEcho, seqValid, steadyState, formula, tiNullBloch, lorentzOffsets, ensemble,
    fft, dft, besselJ1, jinc, ellipseFT, OS, defaults, lineplan, scanTime, noiseSigma, motionAt, metalField,
    kspaceAnalytic, kspaceGrid, acquire, reconstruct, magnitude, pixelXY, roiMask, stats, theory, predictedSNR, PRESETS, rng, gauss,
  };
});
