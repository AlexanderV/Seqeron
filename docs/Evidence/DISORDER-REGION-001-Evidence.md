# Evidence Artifact: DISORDER-REGION-001

**Test Unit ID:** DISORDER-REGION-001
**Algorithm:** Disordered Region Detection (IDR Identification and Classification)
**Date Collected:** 2026-02-12

---

## Online Sources

### Campen et al. (2008) — TOP-IDP Scale

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/
**Citation:** Campen A, Williams RM, Brown CJ, Meng J, Uversky VN, Dunker AK (2008). "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder." Protein Pept Lett. 15(9):956-963. PMC2676888, PMID 18991772.
**Accessed:** 2026-02-12
**Authority rank:** 1 (Peer-reviewed paper, primary research)

**Key Extracted Points:**

1. **TOP-IDP prediction cutoff:** 0.542 based on maximum-likelihood methods.
2. **Window size:** 21 residues used for evaluation in the web-based prediction server.
3. **Prediction method:** Average normalized TOP-IDP values over a sliding window. Positive prediction values indicate ordered, negative indicate disordered.
4. **Scale values (Table 2):** W(-0.884), F(-0.697), Y(-0.510), I(-0.486), M(-0.397), L(-0.326), V(-0.121), N(0.007), C(0.020), T(0.059), A(0.060), G(0.166), R(0.180), D(0.192), H(0.303), Q(0.318), K(0.586), S(0.341), E(0.736), P(0.987).
5. **Region identification:** The TOP-IDP web server performs window-by-window prediction; contiguous windows predicted as disordered constitute a disordered region.

### Dunker et al. (2001) — Intrinsically Disordered Protein

**URL:** https://doi.org/10.1016/s1093-3263(00)00138-8
**Citation:** Dunker AK, et al. (2001). "Intrinsically disordered protein." J Mol Graph Model. 19(1):26-59. PMID 11381529.
**Accessed:** 2026-02-12
**Authority rank:** 1 (Peer-reviewed, foundational paper)

**Key Extracted Points:**

1. **Disorder-promoting amino acids:** {A, R, G, Q, S, P, E, K} — polar, charged, small residues.
2. **Order-promoting amino acids:** {W, C, F, I, Y, V, L, N} — bulky hydrophobic and aromatic residues.
3. **Ambiguous amino acids:** {D, H, M, T} — found in both ordered and disordered regions.
4. **Long disordered regions:** >30 residues are functionally significant, found in ~33% of eukaryotic proteins.
5. **Region types by function:** Molecular recognition domains, flexible linkers, entropic springs, entropic bristles.

### van der Lee et al. (2014) — Classification of IDRs and IDPs

**URL:** https://doi.org/10.1021/cr400525m
**Citation:** van der Lee R, et al. (2014). "Classification of intrinsically disordered regions and proteins." Chemical Reviews. 114(13):6589-6631. PMC4095912, PMID 24773235.
**Accessed:** 2026-02-12
**Authority rank:** 1 (Peer-reviewed review, Chemical Reviews)

**Key Extracted Points:**

1. **IDR functional classification:** IDRs can be classified by amino acid composition biases — proline-rich, acidic, basic, and Ser/Thr-rich regions are recognized subtypes.
2. **Compositional bias:** Low-complexity regions within IDRs often have biased amino acid compositions (homo-repeats, polyQ, polyE, etc.).
3. **Length-based classification:** Short IDRs (<30 residues) vs long IDRs (≥30 residues) differ in functional properties. Long IDRs are more likely to be functionally autonomous.
4. **Region boundaries:** The boundaries of IDRs are defined by the transition from residues predicted as disordered to residues predicted as ordered.

### Wikipedia — Intrinsically Disordered Proteins

**URL:** https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins
**Accessed:** 2026-02-12
**Authority rank:** 4 (Wikipedia citing Dunker 2001 [ref 2], Campen 2008 [ref 71], van der Lee 2014 [ref 11])

**Key Extracted Points:**

1. **Prevalence of long disorder:** "long (>30 residue) disordered segments occur in 2.0% of archaean, 4.2% of eubacterial and 33.0% of eukaryotic proteins" — based on DISOPRED2 prediction (Ward et al. 2004).
2. **Low complexity and disorder:** "Many disordered proteins also reveal low complexity sequences, i.e. sequences with over-representation of a few residues."
3. **Disorder-promoting classification:** Confirms Dunker (2001) classification with {A, R, G, Q, S, P, E, K} as disorder-promoting.

---

## Documented Corner Cases and Failure Modes

### From Campen et al. (2008)

1. **Short sequences:** Sequences shorter than the window size (21 residues) have reduced prediction accuracy due to boundary effects.
2. **Border effects:** At region transitions (order→disorder), window averaging may blur boundaries by ±half-window residues.

### From Algorithm Design

1. **Empty predictions list:** An empty list of per-residue predictions should yield no regions.
2. **All ordered:** If no residue exceeds the disorder threshold, no regions should be identified.
3. **All disordered:** If all residues exceed the threshold and the sequence length ≥ minLength, exactly one region should span the entire sequence.
4. **Discontinuous disorder:** Isolated disordered residues (fewer than minLength) should not form a region.
5. **Trailing region:** A disordered region that extends to the end of the sequence must be captured (no off-by-one).
6. **Region classification with equal composition:** When multiple amino acid fractions tie or none exceeds 0.25, the fallback classification applies.

---

## Test Datasets

### Dataset: Pure Proline Sequence (Proline-rich Classification)

**Source:** Campen et al. (2008), P has highest TOP-IDP propensity (0.987). Proline-rich regions are a known IDR subtype per van der Lee et al. (2014).

| Parameter | Value |
|-----------|-------|
| Sequence | `PPPPPPPPPPPPPPPPPPPPPPPPPPPPPP` (30×P) |
| Expected normalized TOP-IDP score | (0.987 - (-0.884)) / 1.871 = 1.0 |
| Expected prediction | All residues disordered (score >> 0.542) |
| Expected region count | 1 (with minLen ≤ 30) |
| Expected classification | "Proline-rich" (P fraction = 1.0 > 0.25) |

### Dataset: Pure Glutamate Sequence (Acidic Classification)

**Source:** Campen et al. (2008), E has TOP-IDP propensity 0.736. Acidic IDRs are recognized per van der Lee et al. (2014).

| Parameter | Value |
|-----------|-------|
| Sequence | `EEEEEEEEEEEEEEEEEEEEEEEEEEEEEE` (30×E) |
| Expected normalized TOP-IDP score | (0.736 - (-0.884)) / 1.871 ≈ 0.866 |
| Expected prediction | All residues disordered (0.866 > 0.542) |
| Expected classification | "Acidic" (E fraction = 1.0 > 0.25) |

### Dataset: Lysine/Arginine Sequence (Basic Classification)

**Source:** K has TOP-IDP propensity 0.586, R has 0.180. Basic IDRs per van der Lee (2014).

| Parameter | Value |
|-----------|-------|
| Sequence | `KKKKKKKKKKRRRRRRRRRRKKKKKKKKKKRRRRRRRRRR` (40aa) |
| K fraction | 0.5 > 0.25 → basic dominates |
| Expected classification | "Basic" |

### Dataset: Serine/Threonine Sequence (Ser/Thr-rich Classification)

**Source:** S has TOP-IDP propensity 0.341. Ser/Thr-rich IDRs per van der Lee (2014). Note: T=0.059, normalized = (0.059+0.884)/1.871 ≈ 0.504, which is below 0.542 cutoff.

| Parameter | Value |
|-----------|-------|
| Sequence | `SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS` (31×S) |
| S normalized TOP-IDP | (0.341 - (-0.884)) / 1.871 ≈ 0.655 |
| Expected prediction | All residues disordered (0.655 > 0.542) |
| Expected classification | "Ser/Thr-rich" (S fraction = 1.0 > 0.25) |

### Dataset: Ordered Sequence (No Regions)

**Source:** W has lowest TOP-IDP propensity (-0.884), normalized = 0.0.

| Parameter | Value |
|-----------|-------|
| Sequence | `WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW` (30×W) |
| Expected prediction | All residues ordered (score = 0.0 << 0.542) |
| Expected region count | 0 |

### Dataset: Mixed Ordered-Disordered Sequence

**Source:** Derived from Campen (2008) scale values.

| Parameter | Value |
|-----------|-------|
| Sequence | `WWWWWWWWWW` + `PPPPPPPPPPPPPPPPPPPPPPPPPPPPPP` + `WWWWWWWWWW` |
| Expected regions | At least 1 disordered region in the central P-rich segment |
| Expected disordered start | Near position 10 (boundary may shift by window/2) |

---

## Recommendations for Test Coverage

1. **MUST Test:** Empty predictions list produces no regions — trivially correct
2. **MUST Test:** All-ordered predictions produce no regions — trivially correct
3. **MUST Test:** Contiguous disordered predictions produce exactly one region — Algorithm definition
4. **MUST Test:** Region boundaries (Start, End) are correct — Algorithm definition
5. **MUST Test:** Region MeanScore equals average of constituent residue scores — Algorithm definition
6. **MUST Test:** Regions shorter than minLength are excluded — Algorithm definition
7. **MUST Test:** Trailing region (disorder at end of sequence) is captured — Algorithm definition
8. **MUST Test:** Proline-rich classification (P fraction > 0.25 → "Proline-rich") — van der Lee (2014), Campen (2008)
9. **MUST Test:** Acidic classification (E/D fraction > 0.25 → "Acidic") — Das & Pappu (2013), van der Lee (2014)
10. **MUST Test:** Basic classification (K/R fraction > 0.25 → "Basic") — Das & Pappu (2013), van der Lee (2014)
11. **MUST Test:** Ser/Thr-rich classification (S/T fraction > 0.25 → "Ser/Thr-rich") — van der Lee (2014) Table 1
12. **MUST Test:** Long IDR classification (length > 30, no composition bias → "Long IDR") — Ward (2004)
13. **MUST Test:** Standard IDR classification (short, no composition bias → "Standard IDR") — fallback
14. **MUST Test:** Confidence values in [0, 1] — invariant
15. **SHOULD Test:** Multiple disordered regions separated by ordered segment — Algorithm definition
16. **SHOULD Test:** Region at very start of sequence — boundary case
17. **SHOULD Test:** MinLength edge case (region exactly == minLength) — boundary case
18. **COULD Test:** Classification priority when multiple biases present — most-specific-bias-first per Campen (2008)

---

## References

1. Campen A, et al. (2008). "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder." Protein Pept Lett. 15(9):956-963. https://doi.org/10.2174/092986608785849164
2. Dunker AK, et al. (2001). "Intrinsically disordered protein." J Mol Graph Model. 19(1):26-59. https://doi.org/10.1016/s1093-3263(00)00138-8
3. van der Lee R, et al. (2014). "Classification of intrinsically disordered regions and proteins." Chemical Reviews. 114(13):6589-6631. https://doi.org/10.1021/cr400525m
4. Ward JJ, et al. (2004). "Prediction and functional analysis of native disorder in proteins from the three kingdoms of life." J Mol Biol. 337(3):635-645. https://doi.org/10.1016/j.jmb.2004.02.002
5. Kyte J, Doolittle RF (1982). "A simple method for displaying the hydropathic character of a protein." J Mol Biol. 157(1):105-132. https://doi.org/10.1016/0022-2836(82)90515-0
6. Das RK, Pappu RV (2013). "Conformations of intrinsically disordered proteins are influenced by linear sequence distributions of oppositely charged residues." Proc Natl Acad Sci USA. 110(33):13392-13397. https://doi.org/10.1073/pnas.1304749110

---

## Change History

- **2026-02-12**: Initial documentation.
