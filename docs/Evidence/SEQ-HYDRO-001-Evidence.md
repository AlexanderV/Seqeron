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

ProtScale "exponential" weight-variation model: formula not obtainable (see B03.md F30 for the sources tried); not implemented.

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

## Change History

- **2026-06-13**: Initial documentation (SEQ-HYDRO-001).
- **2026-09-28**: Review 2026-09 — added Biopython 1.88 reference dataset incl. edge-weighted profiles (linear ProtScale model now implemented), window < 1 rejection.
- **2026-10-14**: FIN-B03 WP6/WP7 — Biopython window convention, strict mode, ProtScale 0..1 normalisation, caller-supplied scales (hw/em/ja/Flex) datasets; exponential model BLOCKED.
