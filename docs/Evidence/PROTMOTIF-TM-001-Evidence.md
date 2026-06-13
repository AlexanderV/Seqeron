# Evidence Artifact: PROTMOTIF-TM-001

**Test Unit ID:** PROTMOTIF-TM-001
**Algorithm:** Transmembrane Helix Prediction (Kyte-Doolittle hydropathy sliding window)
**Date Collected:** 2026-06-14

---

## Online Sources

### Kyte-Doolittle Background (Davidson College, Genomics / DGPB)

**URL:** https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm
**Accessed:** 2026-06-14
**Authority rank:** 4 (educational page citing the primary Kyte & Doolittle 1982 paper)
**Retrieved via:** WebSearch "Kyte Doolittle hydropathy scale 1982 window 19 transmembrane threshold paper", then WebFetch of the page above.

**Key Extracted Points:**

1. **Primary citation:** "Kyte, J. and Doolittle, R. 1982. A simple method for displaying the hydropathic character of a protein. *J. Mol. Biol.* 157: 105-132." (extracted verbatim).
2. **Window size for TM detection:** "a window size of 19 is needed" to locate a transmembrane region (extracted verbatim).
3. **Threshold:** "Transmembrane regions are identified by peaks with scores greater than 1.6" (extracted verbatim).
4. **Windowing rule:** "The computer program starts with the first window of amino acids and calculates the average of all the hydrophobicity scores in that window. Then the computer program moves down one amino acid and calculates the average of all the hydrophobicity scores in the second window." (extracted verbatim) — i.e. the profile point is the arithmetic mean of the window's per-residue scores, sliding one residue at a time.

### QIAGEN CLC Genomics Workbench — Hydrophobicity scales

**URL:** https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation / vendor reference table, citing Kyte & Doolittle 1982)
**Retrieved via:** WebSearch (same query), then WebFetch of the page above.

**Key Extracted Points:**

1. **Kyte-Doolittle scale (one-letter):** A 1.80, C 2.50, D −3.50, E −3.50, F 2.80, G −0.40, H −3.20, I 4.50, K −3.90, L 3.80, M 1.90, N −3.50, P −1.60, Q −3.50, R −4.50, S −0.80, T −0.70, V 4.20, W −0.90, Y −1.30 (extracted verbatim from the table). These are the exact values used in the implementation's `HydropathyScale`.
2. **Citation:** "Kyte and Doolittle, 1982" referenced as the source of the scale.

### Davidson College — per-amino-acid scores page

**URL:** https://gcat.davidson.edu/DGPB/kd/aminoacidscores.htm
**Accessed:** 2026-06-14
**Authority rank:** 4
**Retrieved via:** WebFetch.

**Key Extracted Points:**

1. **Hydropathy values (independent confirmation of QIAGEN):** I 4.5, V 4.2, L 3.8, F 2.8, C 2.5, M 1.9, A 1.8, G −0.4, T −0.7, S −0.8, W −0.9, Y −1.3, P −1.6, H −3.2, E −3.5, Q −3.5, D −3.5, N −3.5, K −3.9, R −4.5. Matches the QIAGEN table exactly.

### Biopython Bio.SeqUtils.ProtParam (reference implementation source)

**URL:** https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParam.py
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation, Biopython 1.x)
**Retrieved via:** WebSearch "Biopython ProtParam protein_scale Kyte-Doolittle window edge weight hydropathy implementation", then WebFetch of the GitHub source.

**Key Extracted Points:**

1. **`protein_scale(self, param_dict, window, edge=1.0)`** computes, for each window start `i`, a weighted sum normalized by `sum_of_weights = sum(weights) * 2 + 1`. With `edge = 1.0` every weight is 1.0, so the score is the **arithmetic mean** of the window's per-residue scale values — identical to the windowing rule in the Davidson page and to this unit's implementation.
2. **`gravy()`** confirms the Kyte-Doolittle hydropathy convention ("Grand Average of Hydropathy according to Kyte and Doolittle, 1982"): `total_gravy = sum(scale[aa] for aa in sequence) / length`.

### Transmembrane α-helix length (bilayer-spanning requirement)

**URL:** (WebSearch result summary) — "Transmembrane Domains" / biochemistry textbook problem, surfaced by WebSearch "transmembrane alpha helix typical length amino acid residues span lipid bilayer ~20".
**Accessed:** 2026-06-14
**Authority rank:** 4
**Retrieved via:** WebSearch.

**Key Extracted Points:**

1. **Helix length to span the bilayer:** "A single alpha-helix requires about 18 to 21 amino acid residues to span the width of a cell's lipid bilayer," each residue contributing ≈0.15 nm; a 20-residue helix is ≈3.0 nm, matching the ≈3–4 nm bilayer thickness. This justifies the 19-residue scanning window and the minimum-span filter equal to the window width.

---

## Documented Corner Cases and Failure Modes

### From Davidson background page

1. **Window longer than sequence:** A window average is undefined when fewer than `windowSize` residues exist; no profile point can be computed (implementation returns no segments).

### From Biopython ProtParam

1. **Non-standard residues:** scales are defined only for the 20 standard amino acids; characters outside the scale (X, B, Z, *) have no value. The implementation excludes them from the window mean (mean of the residues that do have values).

---

## Test Datasets

### Dataset: Synthetic single-segment hydrophobic stretch

**Source:** Direct application of the Kyte-Doolittle window=19 / threshold=1.6 rule (Davidson background page) to a constructed sequence; expected values computed by hand from the scale.

| Parameter | Value |
|-----------|-------|
| Sequence | `D`×10 + `L`×20 + `D`×10 (40 residues) |
| Window / threshold | 19 / 1.6 |
| Profile length | 40 − 19 + 1 = 22 |
| First window mean ≥ 1.6 at profile index | 5 |
| Last window mean ≥ 1.6 at profile index | 16 |
| Reported segment (Start, End) | (5, 35) |
| Peak score | 3.8 (any all-L window: mean of 19 × 3.8 = 3.8) |

### Dataset: All-hydrophilic sequence (negative control)

**Source:** Same rule; Asp (D) = −3.5, so every window mean = −3.5 < 1.6.

| Parameter | Value |
|-----------|-------|
| Sequence | `D`×40 |
| Reported segments | none (empty) |

### Dataset: Exactly one window of poly-Leu

**Source:** Same rule.

| Parameter | Value |
|-----------|-------|
| Sequence | `L`×19 |
| Profile length | 1 |
| Profile value | 3.8 (≥ 1.6) |
| Reported segment (Start, End, Score) | (0, 18, 3.8) |

---

## Assumptions

1. **ASSUMPTION: Segment `End` reporting convention** — The published rule defines the threshold-crossing window run; it does not prescribe how the window-indexed run maps back to exact residue boundaries. The implementation reports `End = lastProfileIndex + windowSize` (clamped to the last residue index). This is an output-coordinate convention, not a change to the detection rule (which residues qualify), and is documented in the algorithm doc §5.4. Expected boundary values in tests are derived from this stated convention.

---

## Recommendations for Test Coverage

1. **MUST Test:** Single hydrophobic stretch yields exactly one segment at the expected (Start, End) with peak score 3.8 — Evidence: Davidson window/threshold rule + Kyte-Doolittle scale (QIAGEN).
2. **MUST Test:** All-hydrophilic sequence yields no segments — Evidence: scale value D = −3.5 < 1.6.
3. **MUST Test:** Kyte-Doolittle scale value lookup via a 19-residue uniform window (mean = the residue's own score) reproduces the published per-residue values — Evidence: QIAGEN / Davidson scale tables.
4. **MUST Test:** Null / empty / shorter-than-window input returns empty — Evidence: window undefined when < windowSize residues.
5. **SHOULD Test:** Non-standard residues are excluded from the window mean — Rationale: Biopython scale-coverage behavior.
6. **SHOULD Test:** Custom threshold raises/lowers segment count deterministically — Rationale: threshold is a stated parameter.
7. **COULD Test:** Property — every reported `Score` ≥ `threshold` and `Start ≤ End` — Rationale: detection invariant.

---

## References

1. Kyte J, Doolittle RF. 1982. A simple method for displaying the hydropathic character of a protein. *Journal of Molecular Biology* 157(1):105-132. https://doi.org/10.1016/0022-2836(82)90515-0
2. Davidson College, Department of Biology — Genomics Project. Kyte-Doolittle background. https://gcat.davidson.edu/DGPB/kd/kyte-doolittle-background.htm (accessed 2026-06-14).
3. QIAGEN CLC Genomics Workbench manual — Hydrophobicity scales. https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Hydrophobicity_scales.html (accessed 2026-06-14).
4. Biopython — Bio.SeqUtils.ProtParam (`protein_scale`, `gravy`). https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/ProtParam.py (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
