# Evidence Artifact: SEQ-HYDRO-001

**Test Unit ID:** SEQ-HYDRO-001
**Algorithm:** Hydrophobicity Analysis (Kyte-Doolittle GRAVY index and sliding-window hydropathy profile)
**Date Collected:** 2026-06-13

---

## Online Sources

### Biopython `Bio/SeqUtils/ProtParamData.py` (kd scale, master branch)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParamData.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation in an established library)

Retrieved via WebFetch of the raw GitHub file, prompting for the complete `kd` dictionary.

**Key Extracted Points:**

1. **Kyte-Doolittle scale (`kd` dict), all 20 residues (verbatim values):**
   A 1.8, R −4.5, N −3.5, D −3.5, C 2.5, Q −3.5, E −3.5, G −0.4, H −3.2, I 4.5,
   L 3.8, K −3.9, M 1.9, F 2.8, P −1.6, S −0.8, T −0.7, W −0.9, Y −1.3, V 4.2.
2. **Provenance:** the file attributes the scale to Kyte & Doolittle, J Mol Biol (1982).

### Biopython `Bio/SeqUtils/ProtParam.py` (gravy / protein_scale, master branch)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParam.py
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation)

Retrieved via WebFetch of the raw GitHub file, prompting for the `gravy()` and `protein_scale()` method bodies.

**Key Extracted Points:**

1. **GRAVY formula:** `gravy()` returns `total_gravy / self.length`, where `total_gravy` is the sum of the scale value of each residue. GRAVY = (sum of hydropathy values) / (sequence length).
2. **Sliding window (`protein_scale`):** loops `for i in range(self.length - window + 1)`, so for a sequence of length N and window W it returns exactly **N − W + 1** values.
3. **Window weighting:** `protein_scale` default `edge=1.0` gives every position in the window equal weight (no edge down-weighting), i.e. an unweighted mean over the window.

### Expasy ProtParam documentation (GRAVY definition)

**URL:** https://web.expasy.org/protparam/protparam-doc.html
**Accessed:** 2026-06-13
**Authority rank:** 2 (official tool specification, EBI/SIB)

Retrieved via WebFetch, prompting for the verbatim GRAVY definition and its reference.

**Key Extracted Points:**

1. **GRAVY definition (verbatim):** "The GRAVY value for a peptide or protein is calculated as the sum of hydropathy values of all the amino acids, divided by the number of residues in the sequence."
2. **Hydropathy source:** the values come from Kyte, J. and Doolittle, R.F. (1982). Positive GRAVY = hydrophobic, negative = hydrophilic.

### GCAT Davidson — Kyte-Doolittle background

**URL:** https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm
**Accessed:** 2026-06-13
**Authority rank:** 4 (educational page citing the Kyte-Doolittle 1982 primary)

Retrieved via WebFetch, prompting for window-size recommendations and the transmembrane threshold.

**Key Extracted Points:**

1. **Window size 9:** "a window size of 9 was found to give the best results" for surface regions of globular proteins.
2. **Window size 19:** "a window size of 19 is needed" and "Transmembrane regions are identified by peaks with scores greater than 1.6 using a window size of 19."

### alakazam (CRAN) `gravy` documentation

**URL:** https://rdrr.io/cran/alakazam/man/gravy.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation in a maintained package)

Retrieved via WebFetch, prompting for default scale and citation.

**Key Extracted Points:**

1. **Default scale:** Kyte & Doolittle scale, citing "Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol. 157, 105-32 (1982)."

---

## Documented Corner Cases and Failure Modes

### From Biopython ProtParam

1. **Window larger than sequence:** the loop `range(N - W + 1)` produces zero iterations when W > N, so the profile is empty.
2. **Unknown residues:** Biopython `gravy()` indexes the scale dictionary directly and raises `KeyError` on a residue absent from the chosen scale. The standard scale defines only the 20 canonical amino acids; ambiguity codes (B, Z, X), gaps, and stop are undefined.

### From GCAT Davidson / Kyte-Doolittle 1982

1. **Window-size dependence:** the chosen window changes the profile values; W=9 for surface exposure, W=19 (peaks > 1.6) for transmembrane segment detection.

---

## Test Datasets

### Dataset: Single-residue and short-peptide GRAVY (explicit derivation from the kd scale)

**Source:** kd values from Biopython ProtParamData.py; GRAVY = sum/length per Expasy ProtParam doc.

| Sequence | Derivation | Expected GRAVY |
|----------|-----------|----------------|
| `A` | 1.8 / 1 | 1.8 |
| `AG` | (1.8 + (−0.4)) / 2 = 1.4/2 | 0.7 |
| `FLIV` | (2.8 + 3.8 + 4.5 + 4.2) / 4 = 15.3/4 | 3.825 |
| `RKDE` | (−4.5 + −3.9 + −3.5 + −3.5) / 4 = −15.4/4 | −3.85 |

### Dataset: Sliding-window hydropathy profile (window = 3)

**Source:** kd values (Biopython); window mean per Biopython `protein_scale` (edge=1.0).

| Sequence | Window | Windows (N−W+1) | Expected profile |
|----------|--------|-----------------|------------------|
| `FLIV` | 3 | 2 | [(2.8+3.8+4.5)/3, (3.8+4.5+4.2)/3] = [3.7, 4.1666666667] |
| `AG` | 3 | 0 (W > N) | empty |

### Dataset: Biopython 1.88 `ProteinAnalysis` reference values (review 2026-09, computed locally)

| Sequence | Call | Biopython 1.88 output |
|----------|------|-----------------------|
| ubiquitin (76 aa, `MQIFVKTLTG…RLRGG`) | `gravy()` | −0.48947368421052634 |
| ubiquitin | `protein_scale(kd, 9, 1.0)` | 68 values; [0.93333…, 0.67777…, 0.63333…, 0.05555…, 0.24444…, −0.3, …] |
| ubiquitin | `protein_scale(kd, 5, 0.4)` | 72 values; [2.015625, 1.853125, 1.546875, 0.371875, 0.196875, 0.34375, …] |
| ubiquitin | `protein_scale(kd, 9, 0.4)` | [1.095, 0.70416…, 0.5375, 0.10083…, −0.01333…, −0.315, …] |
| ubiquitin | `protein_scale(kd, 19, 0.1)` | 58 values; [0.28, 0.178, 0.138, 0.07, 0.002, −0.059, …] |
| `MKWVTFISLLLLFSSAYS` | `protein_scale(kd, 9, 0.4)` | [1.24417, 1.66083, 2.32167, 2.63, 2.77167, 2.8625, …] |
| `FLIVAG` | `protein_scale(kd, 5, 0.4)` / `(kd, 3, 0.0)` | [3.73125, 3.115625] / [3.8, 4.5, 4.2, 1.8] |
| `FXIVA` | `protein_scale(kd, 3)` | first value 2.43333 (unknown centre → 0, divisor 3) |
| `AV` | `protein_scale(kd, 2)` | 3.4 (even-window artefact: V counted twice, /3) — library returns 3.0 by design |

All reproduced by the library to ≤ 1e-10.

### Dataset: FIN-B03 WP6/WP7 — Biopython 1.88 `protein_scale` options (executed 2026-10-14)

Sequence `S = MKTAYIAKQRQISFVKSHFSRQ`; script: `ProteinAnalysis(seq).protein_scale(table, W, edge)`; normalised = the same call on `{k: (v − min)/(max − min)}` of the 20-entry table. Full arrays are locked in `SequenceStatistics_HydrophobicityProfileOptions_Tests` (tolerance 1e-12).

| Case | Biopython 1.88 output (first values) |
|------|--------------------------------------|
| `S`, kd, W=4, edge 1 | −0.32, −0.46, 0.6, 2.26, … (19 values; residue W/2 twice, /5) |
| `S`, kd, W=2, edge 1 | −1.96667, −1.76667, 0.96667, … (21 values) |
| `S`, kd, W=6, edge 0.5 | 0.451923…, 0.146154…, 1.461538…, … (even weighted window accepted) |
| `ACDXEFGHIK`, kd, W=5, edge 1 / 0.4 | −1.04, −0.34, −1.48, −0.22, 0.04, −0.04 / −1.30625, −0.86875, −1.58125, 0.021875, −0.0875, −0.240625 (lower-case input identical) |
| `AXCDE` / `ACXDE`, kd, W=5 | 0.16 (pair X,D dropped) / −0.54 |
| `AC-DEB`, kd, W=4 | −0.34, −0.9, −2.1 |
| `ACDE`, kd, W=1 | ZeroDivisionError |
| `ProteinAnalysis("ACDX").gravy()` | KeyError 'X' |
| `S`, hw / em / ja / Flex, W=9 | −0.0666…/1.0731111…/−0.3122222…/1.0015555… (14 values each) |
| `S`, hw / em / ja / Flex, W=5, edge 0.4 | −0.028125/1.111625/−0.358125/1.003 (18 values each) |
| `S`, kd normalised, W=9 / W=5 edge 0.4 | 0.459259…, 0.380247… / 0.432986…, 0.522222… |
| `S`, hw normalised, W=5 edge 0.4 | 0.526855…, 0.437988… |

ProtScale "exponential" weight-variation model: formula not obtainable in FIN-B03 (B03.md F30); **superseded by "SRC RO2" below — derived from the live ProtScale tool and implemented.**

---

## Assumptions

1. **ASSUMPTION: Unknown-residue handling diverges from Biopython.** The repository implementation *skips* residues not in the kd scale (GRAVY divides by the count of recognized residues; profile treats them as 0), whereas Biopython raises `KeyError`. Kyte-Doolittle 1982 and the Expasy doc define values only for the 20 standard residues and are silent on ambiguity codes/gaps, so neither rule is mandated by an authoritative source. This is an API-shape/robustness choice, not a scoring constant: every value that *is* produced for in-alphabet residues remains exactly source-conformant. Documented as a deviation in the algorithm doc §5.4; it does not affect any GRAVY/profile value over the 20 canonical residues.

---

## Recommendations for Test Coverage

1. **MUST Test:** GRAVY of single residue, short hydrophobic and hydrophilic peptides equals the exact sum/length derivation. — Evidence: kd scale (Biopython) + GRAVY def (Expasy).
2. **MUST Test:** GRAVY is case-insensitive (lowercase input gives the same value). — Evidence: implementation uppercases; scale defined on uppercase letters.
3. **MUST Test:** sliding-window profile returns exactly N−W+1 unweighted window means. — Evidence: Biopython `protein_scale` loop and `edge=1.0`.
4. **MUST Test:** window > length yields an empty profile; empty/null input yields GRAVY 0 / empty profile. — Evidence: Biopython `range(N−W+1)` produces 0 iterations.
5. **SHOULD Test:** unknown residues are skipped in GRAVY (divide by recognized count). — Rationale: documents the deviation from Biopython KeyError behavior.
6. **COULD Test:** transmembrane-style window (W=19) over a hydrophobic stretch exceeds 1.6. — Rationale: confirms the documented biological threshold.

---

## References

1. Kyte, J., Doolittle, R.F. (1982). A simple method for displaying the hydropathic character of a protein. J Mol Biol 157(1):105–132. https://doi.org/10.1016/0022-2836(82)90515-0
2. Biopython. Bio/SeqUtils/ProtParamData.py (kd scale), master. https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParamData.py
3. Biopython. Bio/SeqUtils/ProtParam.py (gravy, protein_scale), master. https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/ProtParam.py
4. Expasy. ProtParam documentation (GRAVY). https://web.expasy.org/protparam/protparam-doc.html
5. GCAT (Davidson College). Kyte-Doolittle background. https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm
6. alakazam (CRAN). gravy {alakazam}. https://rdrr.io/cran/alakazam/man/gravy.html

---

## Addendum 2026-10-14 (FIN-B03 F41, F52) — shipped scales, MCP options

- F52: Biopython 1.88 `Bio/SeqUtils/ProtParamData.py` tables shipped as `ProteinScales.HoppWoods / Emini / Janin / Flexibility` (verbatim; KD forwards to `ProteinPhysicochemistry.KyteDoolittle`). Label correction: `em` is "Emini Surface fractional probability" (Biopython comment), not Eisenberg.
- Reference arrays (Biopython `protein_scale` on MKTAYIAKQRQISFVKSHFSRQ, from F31): hw W9 −0.06666666666666662 …, em W9 1.0731111111111111 …, ja W9 −0.31222222222222223 …, Flex W5 e0.4 1.003 …, hw W5 e0.4 −0.028125 ….
- F41 MCP `hydrophobicity_profile(edgeWeight, convention, strict, normalizeScale)`, `scale` (F52): `protein_scale(kd)` ubiquitin W5 e0.4 2.015625, 1.853125, 1.546875, 0.371875, 0.196875, 0.34375 …; FLIVAG W5 e0.4 3.7312499999999997, 3.1156249999999996, W3 e0 3.8, 4.5, 4.2, 1.8.
- Tests: `ProteinScales_EqualBiopythonProtParamData`; F31 profile tests now use the shipped tables; MCP `Seqeron.Mcp.Analysis.Tests/B03SequenceStatisticsMcpOptionsTests` (F41, +5 F52).
- Exponential weight model: BLOCKED at FIN-B03 (B03 F30, round-2 B11); resolved by SRC RO2 (below).

## Change History

- **2026-06-13**: Initial documentation (SEQ-HYDRO-001).
- **2026-09-28**: Review 2026-09 — added Biopython 1.88 reference dataset incl. edge-weighted profiles (linear ProtScale model now implemented), window < 1 rejection.
- **2026-10-14**: FIN-B03 WP6/WP7 — Biopython window convention, strict mode, ProtScale 0..1 normalisation, caller-supplied scales (hw/em/ja/Flex) datasets; exponential model BLOCKED.
- **2026-10-14**: FIN-B03 R2-WP5 — addendum F41/F52.

## SRC RO2 (2026-10-14) — ProtScale exponential weight-variation model (SRC for B03)

Re-opened because `web.expasy.org` now answers 200. Sources fetched 2026-10-14:

- `https://web.expasy.org/protscale/protscale-doc.html` (200, 11 648 B). It defines the model only in words:
  "Exponential weight variation model — This option makes the weights decrease exponentially from the central
  position to the window edge. This parameter has an effect only if you set the window edge relative weight to a
  value other than 100%." Its worked example (window 7, edge 10 %) prints `10% 12% 39% 100% 39% 12% 10%`.
  **The live tool disagrees with that example** (it prints and uses `0.10 0.13 0.29 1.00 0.29 0.13 0.10`, below),
  so the doc example is stale and the equation was taken from the tool itself.
- `https://web.expasy.org/cgi-bin/protscale/protscale.pl?1` (POST form fields `sequence`, `scale`, `window`,
  `weight_edges`, `weight_var=exponential|linear`, `norm`), results page + `/tmp/scores<N>.txt` (verbose numeric
  output, 3 decimals). The results page prints the weight table ("Weights for window positions 1,..,W, using
  exponential weight variation model"), 2 decimals.

**Method.** A single W in a poly-G run on ProtScale's "Molecular weight" scale (G = 75.000, W = 204.000) gives
score(d) = 75 + 129·w(d)/Σw around the spike, so score(d) − 75 over score(0) − 75 is w(d)/w(0) to ~1e-5.
Collected for W ∈ {3, 5, 7, 9, 11, 15, 21} × edge ∈ {1, 10, 50, 90} % (28 runs) plus edge 0, 12.5, 99, 99.5,
100, −5 %. Observations: w(R) = e exactly; for edge 1 % the weights are exactly geometric (e^(d/R)); in general the
differences w(d) − w(d+1), d ≥ 1, form a geometric series whose ratio r satisfies (1 − e)·r^R = 0.01, with w(0) = 1.

**Equation (R = (W − 1)/2, e = edge fraction, 1 ≤ d ≤ R):**
`w(0) = 1`, `w(d) = (e − 0.01) + (1 − e)·(0.01/(1 − e))^(d/R)`; e = 1 → all weights 1.
Max |measured − equation| over all 28 runs = 6.9e-5 (3-decimal rounding of the score); every printed 2-decimal
weight equals the equation rounded to 2 decimals. Out-of-range inputs the tool accepts (edge −5 %, 99.5 %) also
follow the equation; edge 150 % returns no profile. The library accepts edge ∈ [0, 1] (unchanged contract).

**Profile check.** score = Σ w·v / Σ w over the window, value k at the window centre (ProtScale numbers positions
from 0 in the verbose file: first position = (W − 1)/2). Kyte & Doolittle, sequence
`MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK` (60 aa); max |ProtScale − equation| ≤ 0.000496 in every run (= 3-decimal rounding):

| W | edge | model | n | ProtScale values (3 dp) |
|---|---|---|---|---|
| 7 | 10 % | exponential | 54 | 0.814, -0.048, 1.922, 0.515, -2.376, -2.570, -3.098, -2.042, 1.506, 0.281, 1.677, 1.850, -1.348, -0.839, -1.600, 0.089, -0.925, -2.696, -2.122, 0.076, -2.224, -2.531, -2.037, 1.226, 0.541, 2.439, 2.501, -0.385, 1.576, -0.752, 0.741, 0.299, 2.181, 2.269, -0.015, -1.439, 1.236, -0.491, -1.980, -1.013, -1.544, -2.543, -2.553, -2.000, 0.904, -0.361, -0.394, 0.289, -1.746, -1.938, 0.752, 1.465, -0.639, 1.448 |
| 9 | 25 % | exponential | 52 | -0.088, 0.225, -0.339, -1.211, -1.714, -1.541, -1.005, 0.015, -0.081, 0.433, 0.899, -0.240, -0.764, -1.087, -0.246, -1.323, -1.997, -2.000, -0.558, -1.753, -1.685, -0.904, 0.494, 0.448, 1.113, 1.438, 0.497, 1.198, 0.578, 0.848, 0.144, 1.592, 1.163, 0.199, -0.526, 0.608, -0.539, -1.625, -1.511, -1.077, -1.921, -1.969, -1.277, -0.301, -0.926, -0.605, 0.005, -0.840, -0.810, 0.055, 0.805, -0.488 |
| 21 | 0 % | exponential | 40 | -1.108, 0.088, 0.229, 0.621, 0.495, -0.518, -0.737, -0.951, -0.634, -1.124, -1.678, -1.565, -0.967, -1.487, -1.509, -1.061, 0.204, 0.617, 1.272, 1.323, 0.613, 0.868, 0.309, 0.652, 0.700, 1.281, 1.087, 0.251, -0.343, 0.050, -0.535, -1.147, -1.194, -1.429, -1.821, -1.801, -1.374, -0.424, -0.474, -0.439 |
| 5 | 50 % | exponential | 56 | -0.505, 0.314, 0.891, 1.195, -0.084, -1.555, -2.903, -2.603, -1.810, 0.397, 1.167, 1.604, 0.892, -0.703, -0.410, -1.419, -0.670, -1.635, -1.157, -1.936, -1.322, -2.349, -1.341, -2.012, 0.385, 1.237, 2.002, 2.104, 0.490, 1.097, -0.860, 1.054, 0.667, 1.997, 0.905, 1.012, -0.286, -0.223, -0.787, -0.704, -1.507, -2.046, -2.500, -1.882, -1.756, -0.180, 0.082, 0.100, -0.877, -1.274, -0.594, -0.491, 1.088, 0.051, 1.362, -0.926 |
| 9 | 25 % | linear | 52 | 0.132, -0.065, -0.561, -1.106, -1.519, -1.438, -0.985, -0.324, 0.156, 0.406, 0.514, 0.079, -0.504, -0.918, -0.952, -1.315, -1.500, -1.712, -1.439, -1.688, -1.358, -0.594, 0.099, 0.550, 1.036, 1.252, 0.937, 0.845, 0.767, 0.790, 0.608, 1.050, 0.885, 0.527, 0.062, -0.006, -0.662, -1.245, -1.513, -1.514, -1.813, -1.667, -1.157, -0.748, -0.764, -0.583, -0.320, -0.527, -0.455, -0.149, 0.242, 0.056 |

Spike runs (Molecular-weight scale, poly-G + one W): ProtScale score at distance d = 0 … R from the spike and the
weight table it printed.

| W | edge | scores d = 0…R (3 dp) | printed weights (edge … centre) |
|---|---|---|---|
| 3 | 10 % | 182.500, 85.750 | 0.10, 1.00 |
| 3 | 50 % | 139.500, 107.250 | 0.50, 1.00 |
| 3 | 1 % | 201.471, 76.265 | 0.01, 1.00 |
| 3 | 90 % | 121.071, 116.464 | 0.90, 1.00 |
| 5 | 10 % | 157.179, 90.192, 83.218 | 0.10, 0.18, 1.00 |
| 5 | 50 % | 116.327, 98.173, 95.664 | 0.50, 0.56, 1.00 |
| 5 | 1 % | 180.825, 85.529, 76.058 | 0.01, 0.10, 1.00 |
| 5 | 90 % | 102.782, 100.605, 100.004 | 0.90, 0.92, 1.00 |
| 7 | 10 % | 137.887, 93.290, 83.478, 81.289 | 0.10, 0.13, 0.29, 1.00 |
| 7 | 50 % | 104.964, 93.749, 90.786, 89.982 | 0.50, 0.53, 0.63, 1.00 |
| 7 | 1 % | 158.737, 92.920, 78.874, 75.837 | 0.01, 0.05, 0.21, 1.00 |
| 7 | 90 % | 94.859, 93.596, 93.102, 92.873 | 0.90, 0.91, 0.94, 1.00 |
| 9 | 10 % | 125.083, 94.142, 84.259, 81.050, 80.008 | 0.10, 0.12, 0.18, 0.38, 1.00 |
| 9 | 50 % | 98.409, 90.872, 88.126, 87.093, 86.705 | 0.50, 0.52, 0.56, 0.68, 1.00 |
| 9 | 1 % | 142.547, 96.200, 81.721, 77.131, 75.675 | 0.01, 0.03, 0.10, 0.31, 1.00 |
| 9 | 90 % | 90.447, 89.616, 89.236, 89.022, 88.902 | 0.90, 0.91, 0.92, 0.95, 1.00 |
| 11 | 10 % | 116.352, 93.853, 84.874, 81.223, 79.739, 79.135 | 0.10, 0.11, 0.15, 0.24, 0.46, 1.00 |
| 11 | 50 % | 94.180, 88.784, 86.404, 85.315, 84.818, 84.590 | 0.50, 0.51, 0.54, 0.59, 0.72, 1.00 |
| 11 | 1 % | 131.069, 97.143, 83.833, 78.524, 76.406, 75.561 | 0.01, 0.03, 0.06, 0.16, 0.39, 1.00 |
| 11 | 90 % | 87.637, 87.045, 86.750, 86.565, 86.448, 86.374 | 0.90, 0.91, 0.92, 0.93, 0.95, 1.00 |
| 15 | 10 % | 105.470, 92.161, 85.324, 81.729, 79.838, 78.844, 78.322, 78.047 | 0.10, 0.11, 0.13, 0.16, 0.22, 0.34, 0.56, 1.00 |
| 15 | 50 % | 89.069, 85.917, 84.194, 83.209, 82.646, 82.324, 82.140, 82.035 | 0.50, 0.51, 0.52, 0.54, 0.58, 0.65, 0.78, 1.00 |
| 15 | 1 % | 116.449, 96.284, 86.040, 80.726, 77.970, 76.541, 75.799, 75.414 | 0.01, 0.02, 0.04, 0.07, 0.14, 0.27, 0.51, 1.00 |
| 15 | 90 % | 84.266, 83.913, 83.727, 83.592, 83.495, 83.425, 83.375, 83.339 | 0.90, 0.90, 0.91, 0.92, 0.93, 0.94, 0.96, 1.00 |
| 21 | 10 % | 96.746, 89.436, 84.914, 82.031, 80.192, 79.020, 78.273, 77.796, 77.492, 77.298, 77.175 | 0.10, 0.11, 0.11, 0.13, 0.15, 0.18, 0.24, 0.32, 0.46, 0.66, 1.00 |
| 21 | 50 % | 85.041, 83.315, 82.216, 81.473, 80.970, 80.630, 80.400, 80.245, 80.140, 80.069, 80.021 | 0.50, 0.50, 0.51, 0.52, 0.54, 0.56, 0.59, 0.64, 0.72, 0.83, 1.00 |
| 21 | 1 % | 104.586, 93.499, 86.684, 82.380, 79.661, 77.944, 76.859, 76.174, 75.742, 75.468, 75.296 | 0.01, 0.02, 0.03, 0.04, 0.06, 0.10, 0.16, 0.25, 0.39, 0.63, 1.00 |
| 21 | 90 % | 81.617, 81.415, 81.307, 81.221, 81.153, 81.098, 81.055, 81.021, 80.994, 80.972, 80.955 | 0.90, 0.90, 0.91, 0.91, 0.92, 0.92, 0.93, 0.94, 0.95, 0.97, 1.00 |

Edge cases (W = 9, ratios w(d)/w(0), d = 0…4): edge 0 % → 1, 0.306234, 0.090006, 0.021616, 0.0 (eq.: 0.306228,
0.09, 0.021623, 0); 99 % → all 0.989974 (eq. 0.99); 99.5 % → 0.990919, 0.992097, 0.993414, 0.995009 (eq. 0.990946,
0.992071, 0.993409, 0.995); 100 % → all 1; 12.5 % → 0.401096, 0.208551, 0.145575, 0.124989.

Locked by `SequenceStatistics_HydrophobicityProfileOptions_Tests` region "ProtScale exponential weight variation
model" (4 exponential + 1 linear KD profiles, 6 spike weight tables, all within 0.0005 of the ProtScale output) and
MCP `HydrophobicityProfile_ExponentialWeightModel_MatchesProtScaleWebTool`. Implementation:
`SequenceStatistics.HydropathyWeightModel.Exponential` via `HydrophobicityProfileOptions.WeightModel`
(default `Linear`, unchanged); MCP `hydrophobicity_profile` optional `weightModel`.
