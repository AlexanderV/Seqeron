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
