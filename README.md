# Simulateur d'IRM : de l'aimantation à l'image

🇫🇷 [Français](#français) · 🇬🇧 [English](#english)

**Essayer l'outil :** https://valentin-auger-laboratoire.github.io/portfolio/outils/irm/

> Outil pédagogique : coupe 2D, antenne homogène, impulsions parfaites, pas de diffusion ni d'imagerie parallèle. Valeurs de l'œdème et de l'hématome illustratives. **Pas un logiciel d'imagerie médicale.**

---

## Français

### L'idée
Une IRM ne prend pas de photo : elle **mesure un signal radio** et **calcule** l'image. Ce simulateur refait toute la chaîne, en JavaScript pur et sans bibliothèque, avec une règle : **chaque chose affichée est calculée par la physique, et chaque calcul est vérifié contre un résultat connu.** Aucun artefact n'est dessiné : ils apparaissent tout seuls quand on change un réglage, parce que leur cause est simulée.

### La chaîne simulée
1. **Aimantation** (équations de Bloch, résolues exactement)
   - Une impulsion est une rotation ; entre deux impulsions, Mz revient vers M0 en T1 et Mxy décroît en T2.
   - Le TR est répété jusqu'au **régime établi**, puis le signal est lu à chaque écho.
   - Séquences : **écho de spin**, **écho de spin rapide** (turbo, train de 180° CPMG), **écho de gradient** « spoilé », **inversion-récupération** (FLAIR, STIR, T1 par inversion).
   - T2* : paquets de protons répartis selon une loi de Lorentz (inhomogénéité R2') ; l'écho de spin les rephase au temps TE, l'écho de gradient non.
2. **Espace k** (signal exact)
   - Le patient virtuel est fait d'ellipses emboîtées. La transformée de Fourier d'une ellipse est connue exactement (fonction de Bessel J1) : chaque échantillon est le signal que mesurerait l'antenne.
   - Pendant la lecture : décroissance T2*, fréquence propre de la graisse (−3,5 ppm, soit −224 Hz à 1,5 T).
   - En écho de spin rapide, chaque ligne est acquise à **son** écho (TE effectif au centre).
   - Sur-échantillonnage ×2 en lecture, comme le récepteur d'une vraie machine.
3. **Reconstruction**
   - FFT 2D écrite à la main.
   - Fourier partiel : remplissage par des zéros ou par **symétrie hermitienne**.
   - Filtre de Hamming possible.
4. **Bruit et temps d'examen**
   - Bruit thermique ajouté dans l'espace k.
   - Rapport signal sur bruit ∝ volume du voxel × √(lignes × colonnes × moyennages / bande passante).
   - Durée = TR × lignes × moyennages / facteur turbo.

### Les artefacts, et ce qui les cause
| Artefact | Cause simulée |
|---|---|
| Repliement | champ de vue trop petit en phase (en lecture, le sur-échantillonnage l'évite) |
| Déplacement chimique | graisse décalée de 224 / (bande passante par pixel) pixels en lecture |
| Opposition de phase | eau et graisse en opposition à TE = 2,24 ms : liseré noir aux interfaces |
| Mouvement | translation du patient entre deux lignes → erreur de phase → fantômes, à N·TR/période pixels pour la respiration |
| Phénomène de Gibbs | espace k tronqué → oscillations d'environ 9 % aux bords |
| Métal | champ d'une sphère aimantée en 1/r³ (titane, acier inoxydable) → distorsion en lecture et déphasage intravoxel, trou plus grand en écho de gradient |

Avec du métal, le signal est intégré sur une **grille fine** : chaque petit élément est déplacé de Δf / (bande passante par pixel) pixels et garde sa phase 2π·Δf·(TE − tr). Ce moteur est vérifié contre le moteur exact.

### Les cinq onglets
- **L'image** : image, espace k, patient virtuel, acquisition animée (les lignes arrivent dans l'ordre réel), tableau « mesuré sur l'image / prévu par Bloch ».
- **L'aimantation** : les protons d'un voxel en 3D, le chronogramme de la séquence, Mz et |Mxy| sur un TR.
- **Le contraste** : signal selon TR, TE, TI ou l'angle ; un clic change le réglage.
- **Scanner ou IRM ?** : le même patient vu par le scanner (calculé par le [simulateur de scanner](https://github.com/Valentin-Auger-Laboratoire/scanner-tomographie)) et par l'IRM en T1, T2 et FLAIR.
- **Vérifications** : 14 contrôles recalculés en direct dans le navigateur.

### Valeurs des tissus à 1,5 T
| Tissu | T1 (ms) | T2 (ms) | Densité de protons | Source |
|---|---|---|---|---|
| Substance blanche | 884 | 72 | 0,69 | Stanisz et al., 2005 |
| Substance grise | 1 124 | 95 | 0,80 | Stanisz et al., 2005 |
| Muscle | 1 008 | 44 | 0,72 | Stanisz et al., 2005 |
| Sang | 1 441 | 290 | 0,85 | Stanisz et al., 2005 |
| Graisse sous-cutanée | 343 | 58 | 0,90 | Gold et al., 2004 |
| Liquide céphalo-rachidien | 4 000 | 2 000 | 1,00 | valeurs usuelles |
| Os cortical | 250 | 0,5 | 0,05 | indicatif |
| Œdème (AVC récent) | 1 300 | 130 | 0,84 | illustratif |
| Hématome | 950 | 50 | 0,80 | illustratif (R2' + 60 s⁻¹) |

Les densités de protons sont des valeurs relatives usuelles.

### Fichiers
| Fichier | Rôle |
|---|---|
| `irm.js` | Moteur : tissus, patients virtuels, Bloch, ensemble d'isochromats, FFT, Bessel, transformée des ellipses, acquisition (moteur exact et moteur sur grille), reconstruction, mesures (testable avec Node) |
| `index.html` | Interface bilingue FR / EN ; les calculs tournent dans un Web Worker |
| `scanner.js` | Copie du moteur du [simulateur de scanner](https://github.com/Valentin-Auger-Laboratoire/scanner-tomographie), pour l'onglet « Scanner ou IRM ? » |
| `tests/test.js` | 422 tests |

### Lancer les tests
```bash
node tests/test.js
```
- **Bloch contre les formules analytiques**, à 10⁻⁹ près : écho de spin, inversion-récupération, écho de gradient (180 tirages aléatoires), train d'échos, angle d'Ernst.
- **TI d'annulation** : formule = Bloch ; STIR → T1·ln 2 quand TR ≫ T1.
- **Isochromats** : somme en e^(−t/T2)·e^(−R2'·t) ; rephasage exact à l'écho de spin.
- **FFT** contre la transformée directe ; **J1** contre sa représentation intégrale et ses zéros ; **transformée de l'ellipse** contre une intégration numérique.
- **Patients virtuels** : ellipses bien emboîtées, sans chevauchement, somme des contributions = tissu.
- **Image** : moyennes mesurées dans les tubes du fantôme de contrôle = Bloch à 2 % près (T1, T1 par inversion, écho de gradient), 7 % en turbo.
- **Gibbs** : dépassement mesuré entre 8,5 et 9,5 % ; supprimé par Hamming.
- **Contrastes connus en clinique** : LCR sombre en T1 et clair en T2, éteint en FLAIR ; graisse éteinte en STIR ; hématome sombre en T2*…
- **Bruit** : mesuré = prévu à 3 % ; lois en √moyennages, épaisseur, √bande passante, matrice.
- **Déplacement chimique** : mesuré = 224 / bande passante à 0,12 pixel près, dans les deux moteurs ; aucun en phase.
- **Repliement** en phase à la position prévue, aucun en lecture ; **fantômes de respiration** à N·TR/période pixels (±1).
- **Fourier partiel**, temps d'examen, plan des lignes du turbo.
- **Métal** : champ en 1/r³ ; moteur sur grille = moteur exact sans métal ; trou plus grand en écho de gradient ; acier > titane ; bande passante plus large = moins de distorsion.

### Limites
- **Une coupe 2D** : pas de profil de coupe, pas de déplacement hors du plan près du métal.
- Antenne homogène, impulsions parfaites, pas de diffusion, pas d'imagerie parallèle, pas de flux.
- La graisse est un seul pic à −3,5 ppm.
- Dans le moteur sur grille (métal), la décroissance T2* pendant la lecture est négligée.
- Valeurs de l'œdème et de l'hématome illustratives : un vrai hématome change d'aspect avec son âge.

### Références
- Bloch F. *Nuclear induction.* Phys Rev 70:460-474 (1946).
- Hahn EL. *Spin echoes.* Phys Rev 80:580-594 (1950).
- Stanisz GJ et al. *T1, T2 relaxation and magnetization transfer in tissue at 3T.* Magn Reson Med 54:507-512 (2005), valeurs à 1,5 T.
- Gold GE et al. *Musculoskeletal MRI at 3.0 T: relaxation times and image contrast.* AJR 183:343-351 (2004), valeurs à 1,5 T.
- Schenck JF. *The role of magnetic susceptibility in magnetic resonance imaging.* Med Phys 23:815-850 (1996).
- Haacke EM et al. *Magnetic Resonance Imaging: Physical Principles and Sequence Design.* Wiley, 1999.
- Press WH et al. *Numerical Recipes*, approximation de J1.

### Utiliser en local
Ouvrir `index.html` dans un navigateur (de préférence via un petit serveur local, par exemple `python3 -m http.server`, pour que le Web Worker se lance ; sinon les calculs se font sur la page).

---

## English

### The idea
An MRI scanner does not take a picture: it **measures a radio signal** and **computes** the image. This simulator rebuilds the whole chain in plain JavaScript, with one rule: **everything displayed is computed from physics, and every computation is checked against a known result.** No artefact is drawn: they appear by themselves when a setting changes, because their cause is simulated.

### What is simulated
- **Bloch equations solved exactly** to steady state: spin echo, turbo spin echo (CPMG), spoiled gradient echo, inversion recovery (FLAIR, STIR, T1 inversion). T2* from a Lorentzian spread of spin packets.
- **Exact k-space**: ellipse phantoms with an analytical Fourier transform (Bessel J1); T2* decay and fat frequency (−3.5 ppm) during readout; each turbo line at its own echo; ×2 readout oversampling.
- **Reconstruction** by hand-written 2D FFT, partial Fourier (zeros or Hermitian symmetry), optional Hamming filter.
- **Noise and scan time** following the standard laws.
- **Artefacts from their causes**: wrap-around, chemical shift, opposed phase, motion ghosts, Gibbs ringing, metal (field of a magnetised sphere, fine-grid engine).
- **CT or MRI?** tab: the same patient through the [CT simulator](https://github.com/Valentin-Auger-Laboratoire/scanner-tomographie) and through MRI in T1, T2 and FLAIR.

### Tests
```bash
node tests/test.js
```
422 tests, including:
- Bloch against analytical formulas (10⁻⁹);
- FFT, J1 and the ellipse transform;
- measured tube values against Bloch (2 %);
- Gibbs ≈ 9 %;
- noise laws;
- chemical shift (0.12 px);
- wrap-around;
- motion ghosts at N·TR/period;
- partial Fourier;
- metal, fine grid against the exact engine.

### Limits
- 2D slice, uniform coil, perfect pulses.
- No diffusion, parallel imaging or flow; single fat peak.
- Oedema and haematoma values are illustrative.
- Teaching tool, not medical imaging software.

---

← [Portfolio](https://valentin-auger-laboratoire.github.io/portfolio/#projets)
