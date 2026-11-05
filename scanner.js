/* Scanner (tomodensitométrie) : moteur de simulation et de reconstruction, sans dépendance.
   Fantômes faits d'ellipses (projections calculées exactement), acquisition en géométrie parallèle avec bruit
   de comptage des photons (loi de Poisson), rétroprojection simple, rétroprojection filtrée (filtres rampe,
   Shepp-Logan, Hann, appliqués par transformée de Fourier rapide), reconstruction itérative SART, unités Hounsfield.
   Outil pédagogique : ce n'est pas un logiciel d'imagerie médicale.
   Utilisable dans le navigateur (window.Scanner) et avec Node (module.exports). */
(function (root, factory) {
  const M = factory();
  if (typeof module === 'object' && module.exports) module.exports = M;
  else root.Scanner = M;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const MU_WATER = 0.19; // coefficient d'atténuation de l'eau vers 70 keV (cm⁻¹)
  const hu2mu = (hu) => MU_WATER * (1 + hu / 1000);
  const mu2hu = (mu) => 1000 * (mu / MU_WATER - 1);
  const DEG = Math.PI / 180;

  // ---------- Fantômes ----------
  // ellipse : [x0, y0, a, b (demi-axes, cm), angle (degrés), Δμ (cm⁻¹, ajouté à ce qui est dessous)]
  // Les valeurs s'additionnent : une ellipse incluse dans une autre porte la DIFFÉRENCE de μ.
  function nested(list) {
    // list : [x0, y0, a, b, ang, HU, HU du milieu qui l'entoure]
    return list.map(([x, y, a, b, ang, hu, around]) => [x, y, a, b, ang, hu2mu(hu) - (around === null ? 0 : hu2mu(around))]);
  }
  // Shepp-Logan modifié (Toft), mis à l'échelle d'une tête de 18 cm de haut ; valeurs en unités arbitraires × μ eau
  const SL = [[0, 0, 0.69, 0.92, 0, 1], [0, -0.0184, 0.6624, 0.874, 0, -0.8], [0.22, 0, 0.11, 0.31, -18, -0.2], [-0.22, 0, 0.16, 0.41, 18, -0.2], [0, 0.35, 0.21, 0.25, 0, 0.1], [0, 0.1, 0.046, 0.046, 0, 0.1], [0, -0.1, 0.046, 0.046, 0, 0.1], [-0.08, -0.605, 0.046, 0.023, 0, 0.1], [0, -0.605, 0.023, 0.023, 0, 0.1], [0.06, -0.605, 0.023, 0.046, 0, 0.1]];
  const PHANTOMS = {
    // la tête de Shepp-Logan : référence classique pour comparer les algorithmes (crâne → μ = 1 × échelle)
    shepp: { fov: 24, ellipses: SL.map(([x, y, a, b, ang, v]) => [x * 10, y * 10, a * 10, b * 10, ang, v * MU_WATER * 2]), units: 'mu' },
    // crâne avec un hématome (sang coagulé, 70 HU, à peine plus dense que le cerveau à 35 HU)
    head: {
      fov: 24, units: 'hu', ellipses: nested([
        [0, 0, 7.6, 9.4, 0, 1000, null], // os du crâne
        [0, -0.1, 7.0, 8.8, 0, 35, 1000], // cerveau
        [-1.2, 0.6, 0.55, 1.8, 15, 5, 35], [1.2, 0.6, 0.55, 1.8, -15, 5, 35], // ventricules (liquide céphalo-rachidien)
        [3.6, 2.8, 1.6, 1.1, 30, 70, 35], // hématome
        [0, -5.2, 1.6, 1.0, 0, 20, 35], // cervelet (un peu moins dense)
      ]),
    },
    // thorax : corps, poumons, cœur, aorte, colonne, côtes
    thorax: {
      fov: 36, units: 'hu', ellipses: nested([
        [0, 0, 16, 11, 0, 40, null],
        [-7, 1, 5.2, 7.5, 8, -850, 40], [7, 1, 5.2, 7.5, -8, -850, 40],
        [1.2, 0.5, 3.4, 4.2, -25, 45, 40], [-1.8, -3.8, 1.2, 1.2, 0, 200, 40],
        [0, -7.6, 1.8, 1.6, 0, 700, 40], [0, -7.6, 0.9, 0.8, 0, 250, 700],
        [-12.5, -2, 0.6, 1.8, 30, 900, 40], [12.5, -2, 0.6, 1.8, -30, 900, 40], [-11, 5.5, 0.6, 1.6, -40, 900, 40], [11, 5.5, 0.6, 1.6, 40, 900, 40],
      ]),
    },
  };
  // thorax avec deux implants en acier (μ ≈ 5,9 cm⁻¹) : source d'artefacts en étoile
  PHANTOMS.metal = { fov: 36, units: 'hu', ellipses: PHANTOMS.thorax.ellipses.concat(nested([[-1.8, -6.0, 0.45, 0.45, 0, 30000, 40], [1.9, -6.0, 0.45, 0.45, 0, 30000, 40]])) };

  // valeur du fantôme au point (x, y) en cm⁻¹
  function phantomAt(ph, x, y) {
    let v = 0;
    for (const [x0, y0, a, b, ang, d] of ph.ellipses) {
      const c = Math.cos(ang * DEG), s = Math.sin(ang * DEG), dx = x - x0, dy = y - y0;
      const u = (dx * c + dy * s) / a, w = (-dx * s + dy * c) / b;
      if (u * u + w * w <= 1) v += d;
    }
    return v;
  }
  // image du fantôme (N × N pixels couvrant le champ de vue), échantillonnée avec sur-échantillonnage ss × ss
  function phantomImage(ph, N, ss = 2) {
    const img = new Float64Array(N * N), px = ph.fov / N;
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      let v = 0;
      for (let a = 0; a < ss; a++) for (let b = 0; b < ss; b++) {
        const x = -ph.fov / 2 + (i + (a + 0.5) / ss) * px, y = ph.fov / 2 - (j + (b + 0.5) / ss) * px;
        v += phantomAt(ph, x, y);
      }
      img[j * N + i] = v / (ss * ss);
    }
    return img;
  }
  // projection exacte d'une ellipse : intégrale de ligne le long du rayon d'angle θ et de position s
  function ellipseProjection([x0, y0, a, b, ang, d], theta, s) {
    const t = theta - ang * DEG, a2 = a * a * Math.cos(t) ** 2 + b * b * Math.sin(t) ** 2;
    const sp = s - (x0 * Math.cos(theta) + y0 * Math.sin(theta));
    return sp * sp <= a2 ? ((2 * d * a * b) / a2) * Math.sqrt(a2 - sp * sp) : 0;
  }
  // ---------- Acquisition ----------
  // géométrie parallèle : nA angles répartis sur [0, range), nD cellules de détecteur couvrant le champ de vue
  function acquire(ph, opts) {
    const nA = opts.nA, nD = opts.nD, range = (opts.range || 180) * DEG, fov = ph.fov, ds = fov / nD;
    const thetas = Float64Array.from({ length: nA }, (_, k) => (k * range) / nA);
    const sino = new Float64Array(nA * nD);
    for (let k = 0; k < nA; k++) for (let i = 0; i < nD; i++) {
      const s = -fov / 2 + (i + 0.5) * ds;
      let p = 0; for (const e of ph.ellipses) p += ellipseProjection(e, thetas[k], s);
      sino[k * nD + i] = p;
    }
    let counts = null;
    if (opts.I0) {
      // bruit quantique : nombre de photons détectés ~ Poisson(I0·e^(−p)), puis p mesuré = −ln(N/I0)
      let seed = opts.seed || 1;
      const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
      const poisson = (lam) => {
        if (lam > 50) { const g = Math.sqrt(-2 * Math.log(rnd() + 1e-300)) * Math.cos(2 * Math.PI * rnd()); return Math.max(0, Math.round(lam + Math.sqrt(lam) * g)); }
        const L = Math.exp(-lam); let k = 0, q = 1; do { k++; q *= rnd(); } while (q > L); return k - 1;
      };
      counts = new Float64Array(nA * nD);
      for (let i = 0; i < sino.length; i++) {
        const nph = poisson(opts.I0 * Math.exp(-sino[i]));
        counts[i] = nph;
        sino[i] = -Math.log(Math.max(nph, 0.5) / opts.I0); // aucun photon : « famine de photons », valeur bornée
      }
    }
    return { sino, nA, nD, thetas, ds, fov, counts };
  }

  // ---------- Transformée de Fourier rapide (radix 2, en place) ----------
  function fft(re, im, inv) {
    const n = re.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      for (; j & bit; bit >>= 1) j ^= bit;
      j ^= bit;
      if (i < j) { [re[i], re[j]] = [re[j], re[i]]; [im[i], im[j]] = [im[j], im[i]]; }
    }
    for (let len = 2; len <= n; len <<= 1) {
      const ang = ((inv ? 2 : -2) * Math.PI) / len, wr = Math.cos(ang), wi = Math.sin(ang);
      for (let i = 0; i < n; i += len) {
        let cr = 1, ci = 0;
        for (let k = 0; k < len / 2; k++) {
          const ar = re[i + k + len / 2] * cr - im[i + k + len / 2] * ci, ai = re[i + k + len / 2] * ci + im[i + k + len / 2] * cr;
          re[i + k + len / 2] = re[i + k] - ar; im[i + k + len / 2] = im[i + k] - ai;
          re[i + k] += ar; im[i + k] += ai;
          const t = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = t;
        }
      }
    }
    if (inv) for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; }
  }

  // ---------- Filtres de la rétroprojection filtrée ----------
  // noyau rampe discret exact (Ram-Lak, Kak et Slaney) : h(0) = 1/(4τ²), h(n impair) = −1/(n²π²τ²), 0 sinon
  function filterResponse(nD, ds, type) {
    let M = 1; while (M < 2 * nD) M <<= 1;
    const re = new Float64Array(M), im = new Float64Array(M);
    for (let k = 0; k < M; k++) {
      const n = k <= M / 2 ? k : k - M;
      re[k] = n === 0 ? 1 / (4 * ds * ds) : n % 2 ? -1 / (n * n * Math.PI * Math.PI * ds * ds) : 0;
    }
    fft(re, im, false);
    // fenêtre d'apodisation (réduit le bruit aux hautes fréquences)
    for (let k = 0; k < M; k++) {
      const f = (k <= M / 2 ? k : M - k) / M; // fréquence en cycles par échantillon, de 0 à 0,5
      let w = 1;
      if (type === 'shepp') w = f === 0 ? 1 : Math.sin(Math.PI * f) / (Math.PI * f);
      else if (type === 'hann') w = 0.5 + 0.5 * Math.cos(2 * Math.PI * f);
      else if (type === 'cosine') w = Math.cos(Math.PI * f);
      re[k] *= w * ds; im[k] *= w * ds; // ds : facteur τ de la convolution discrète
    }
    return { re, im, M };
  }
  function filterSinogram(acq, type) {
    const { sino, nA, nD, ds } = acq, H = filterResponse(nD, ds, type), M = H.M;
    const out = new Float64Array(nA * nD), re = new Float64Array(M), im = new Float64Array(M);
    for (let k = 0; k < nA; k++) {
      re.fill(0); im.fill(0);
      for (let i = 0; i < nD; i++) re[i] = sino[k * nD + i];
      fft(re, im, false);
      for (let f = 0; f < M; f++) { const a = re[f] * H.re[f] - im[f] * H.im[f], b = re[f] * H.im[f] + im[f] * H.re[f]; re[f] = a; im[f] = b; }
      fft(re, im, true);
      for (let i = 0; i < nD; i++) out[k * nD + i] = re[i];
    }
    return out;
  }

  // ---------- Rétroprojection ----------
  // pour chaque pixel et chaque angle, on lit la projection (filtrée ou non) en s = x·cos θ + y·sin θ
  function backproject(acq, proj, N, upTo) {
    const { nA, nD, thetas, ds, fov } = acq, img = new Float64Array(N * N), px = fov / N;
    const K = upTo === undefined ? nA : Math.min(nA, upTo);
    const cs = thetas.map(Math.cos), sn = thetas.map(Math.sin);
    for (let j = 0; j < N; j++) {
      const y = fov / 2 - (j + 0.5) * px;
      for (let i = 0; i < N; i++) {
        const x = -fov / 2 + (i + 0.5) * px;
        let v = 0;
        for (let k = 0; k < K; k++) {
          const u = (x * cs[k] + y * sn[k] + fov / 2) / ds - 0.5, i0 = Math.floor(u), fr = u - i0;
          if (i0 >= 0 && i0 + 1 < nD) v += proj[k * nD + i0] * (1 - fr) + proj[k * nD + i0 + 1] * fr;
          else if (i0 === nD - 1) v += proj[k * nD + i0] * (1 - fr);
          else if (i0 === -1) v += proj[k * nD] * fr;
        }
        img[j * N + i] = v;
      }
    }
    return img;
  }
  // ajoute à img la contribution des angles k0 à k1 − 1 (pour l'animation, angle par angle)
  function backprojectAdd(acq, proj, N, k0, k1, img, scale) {
    const { nD, thetas, ds, fov } = acq, px = fov / N;
    for (let k = k0; k < k1; k++) {
      const c = Math.cos(thetas[k]), s = Math.sin(thetas[k]);
      for (let j = 0; j < N; j++) {
        const y = fov / 2 - (j + 0.5) * px;
        for (let i = 0; i < N; i++) {
          const u = ((-fov / 2 + (i + 0.5) * px) * c + y * s + fov / 2) / ds - 0.5, i0 = Math.floor(u), fr = u - i0;
          if (i0 >= 0 && i0 + 1 < nD) img[j * N + i] += scale * (proj[k * nD + i0] * (1 - fr) + proj[k * nD + i0 + 1] * fr);
        }
      }
    }
    return img;
  }
  // rétroprojection filtrée : f = (π / nA) · Σ q_θ(x cos θ + y sin θ), pour des angles répartis sur 180°
  function fbp(acq, N, type = 'ramlak', upTo) {
    const q = filterSinogram(acq, type), img = backproject(acq, q, N, upTo);
    const span = acq.thetas.length > 1 ? acq.thetas[1] - acq.thetas[0] : Math.PI;
    for (let i = 0; i < img.length; i++) img[i] *= span;
    return img;
  }
  // rétroprojection simple (sans filtre) : image floue, normalisée pour garder la même moyenne
  function simpleBP(acq, N, upTo) {
    const img = backproject(acq, acq.sino, N, upTo);
    const span = acq.thetas.length > 1 ? acq.thetas[1] - acq.thetas[0] : Math.PI;
    // la rétroprojection simple vaut f * 1/r : on la ramène à une échelle comparable (division par le diamètre)
    for (let i = 0; i < img.length; i++) img[i] *= span / (acq.fov * 1.0);
    return img;
  }

  // ---------- Projection d'une image (pour la reconstruction itérative) ----------
  // lancer de rayons avec interpolation bilinéaire, pas d'échantillonnage de la moitié d'un pixel
  function forwardProjectAngle(img, N, fov, theta, nD, ds, out) {
    const px = fov / N, c = Math.cos(theta), s = Math.sin(theta), step = px / 2, R = fov / 2;
    const at = (x, y) => {
      const fi = (x + R) / px - 0.5, fj = (R - y) / px - 0.5, i0 = Math.floor(fi), j0 = Math.floor(fj), a = fi - i0, b = fj - j0;
      let v = 0;
      if (i0 >= 0 && j0 >= 0 && i0 + 1 < N && j0 + 1 < N) v = img[j0 * N + i0] * (1 - a) * (1 - b) + img[j0 * N + i0 + 1] * a * (1 - b) + img[(j0 + 1) * N + i0] * (1 - a) * b + img[(j0 + 1) * N + i0 + 1] * a * b;
      return v;
    };
    for (let d = 0; d < nD; d++) {
      const sd = -R + (d + 0.5) * ds, half = Math.sqrt(Math.max(0, R * R - sd * sd));
      let sum = 0;
      for (let t = -half + step / 2; t < half; t += step) sum += at(sd * c - t * s, sd * s + t * c);
      out[d] = sum * step;
    }
  }
  // SART (Andersen et Kak) : on corrige l'image angle par angle, avec la contrainte μ ≥ 0
  function sart(acq, N, iters = 5, lambda = 0.25, init) {
    const { nA, nD, thetas, ds, fov, sino } = acq, px = fov / N, R = fov / 2;
    const img = init ? Float64Array.from(init) : new Float64Array(N * N);
    const proj = new Float64Array(nD), len = new Float64Array(nD);
    for (let d = 0; d < nD; d++) { const sd = -R + (d + 0.5) * ds; len[d] = 2 * Math.sqrt(Math.max(0, R * R - sd * sd)); }
    // ordre des angles entrelacé (accélère la convergence)
    const order = []; const used = new Uint8Array(nA); let k = 0; const stride = Math.max(1, Math.round(nA * 0.618));
    for (let m = 0; m < nA; m++) { while (used[k]) k = (k + 1) % nA; order.push(k); used[k] = 1; k = (k + stride) % nA; }
    const hist = [];
    for (let it = 0; it < iters; it++) {
      for (const a of order) {
        forwardProjectAngle(img, N, fov, thetas[a], nD, ds, proj);
        const c = Math.cos(thetas[a]), s = Math.sin(thetas[a]);
        const corr = new Float64Array(nD);
        for (let d = 0; d < nD; d++) corr[d] = len[d] > px ? (sino[a * nD + d] - proj[d]) / len[d] : 0;
        for (let j = 0; j < N; j++) {
          const y = R - (j + 0.5) * px;
          for (let i = 0; i < N; i++) {
            const x = -R + (i + 0.5) * px, u = (x * c + y * s + R) / ds - 0.5, i0 = Math.floor(u), fr = u - i0;
            if (i0 < 0 || i0 + 1 >= nD) continue;
            const v = img[j * N + i] + lambda * (corr[i0] * (1 - fr) + corr[i0 + 1] * fr);
            img[j * N + i] = v > 0 ? v : 0;
          }
        }
      }
      hist.push(Float64Array.from(img));
    }
    return { img, hist };
  }

  // ---------- Mesures de qualité (dans le cercle du champ de vue) ----------
  function compare(rec, truth, N, units) {
    let se = 0, n = 0, mx = -Infinity, mn = Infinity;
    for (let j = 0; j < N; j++) for (let i = 0; i < N; i++) {
      const x = (i + 0.5) / N - 0.5, y = (j + 0.5) / N - 0.5;
      if (x * x + y * y > 0.24) continue;
      const k = j * N + i, a = units === 'hu' ? mu2hu(rec[k]) : rec[k], b = units === 'hu' ? mu2hu(truth[k]) : truth[k];
      se += (a - b) ** 2; n++; mx = Math.max(mx, b); mn = Math.min(mn, b);
    }
    const rmse = Math.sqrt(se / n);
    return { rmse, psnr: 20 * Math.log10((mx - mn) / rmse) };
  }
  // profil le long d'une ligne horizontale
  const profile = (img, N, row) => Array.from({ length: N }, (_, i) => img[row * N + i]);

  return { MU_WATER, hu2mu, mu2hu, PHANTOMS, phantomAt, phantomImage, ellipseProjection, acquire, fft, filterResponse, filterSinogram, backproject, backprojectAdd, fbp, simpleBP, forwardProjectAngle, sart, compare, profile };
});
