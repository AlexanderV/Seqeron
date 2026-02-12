# Evidence Artifact: DISORDER-PRED-001

**Test Unit ID:** DISORDER-PRED-001
**Algorithm:** Disorder Prediction (Intrinsically Disordered Protein Prediction)
**Date Collected:** 2026-02-10

---

## Online Sources

### Campen et al. (2008) — TOP-IDP Scale (Primary Source)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/
**Citation:** Campen A, Williams RM, Brown CJ, Meng J, Uversky VN, Dunker AK (2008). "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder." Protein Pept Lett. 15(9):956-963. PMC2676888, PMID 18991772.
**Accessed:** 2026-02-11
**Authority rank:** 1 (Peer-reviewed, primary research, Table 2 provides exact numerical values)

**Key Extracted Points:**

1. **TOP-IDP scale:** An optimized amino acid scale for discriminating between order and disorder. Developed by surveying 517 amino acid scales and applying simulated annealing to maximize discrimination (ARV = 0.761, 11% improvement over the best existing scale).
2. **Table 2 values:** The 20 amino acid scale values from order-promoting to disorder-promoting: W(-0.884), F(-0.697), Y(-0.510), I(-0.486), M(-0.397), L(-0.326), V(-0.121), N(0.007), C(0.020), T(0.059), A(0.060), G(0.166), R(0.180), D(0.192), H(0.303), Q(0.318), K(0.586), S(0.341), E(0.736), P(0.987).
3. **Prediction cutoff:** 0.542 based on maximum-likelihood methods.
4. **Window size:** 21 residues used for evaluation.
5. **Confirms Dunker (2001) classification:** Order-promoting = {W, C, F, I, Y, V, L, N}; Disorder-promoting = {A, R, G, Q, S, P, E, K}.

### Wikipedia — Intrinsically Disordered Proteins

**URL:** https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins
**Accessed:** 2026-02-10
**Authority rank:** 3

**Key Extracted Points:**

1. **Definition:** Intrinsically disordered proteins (IDPs) are proteins that lack fixed or ordered 3D structure under physiological conditions. They exist as dynamic ensembles of interconverting structures.
2. **Disorder-promoting amino acids:** IDPs are characterized by high proportion of polar and charged amino acids (Gln, Ser, Pro, Glu, Lys, and on occasion Arg and Gly), and of Ala — these are disorder-promoting.
3. **Order-promoting amino acids:** IDPs have low content of bulky hydrophobic amino acids (Val, Leu, Ile, Met, Phe, Trp, Tyr) and of Cys as well as Asn.
4. **TOP-IDP scale:** Campen et al. (2008) derived the TOP-IDP scale ranking amino acids from most disorder-promoting to most order-promoting based on fractional difference between disordered and ordered proteins.
5. **Charge-hydropathy model:** Uversky et al. (2000) showed IDPs can be distinguished from structured proteins by a combination of low mean hydropathy and high mean net charge — the charge-hydropathy (CH) plot.
6. **MoRFs:** Molecular Recognition Features (MoRFs) are short segments within disordered regions that undergo disorder-to-order transition upon binding to a specific partner (Mohan et al. 2006).
7. **Prevalence:** ~33% of eukaryotic proteins contain long (>30 residues) disordered regions.

### Wikipedia — Hydrophilicity Plot (Kyte-Doolittle Scale)

**URL:** https://en.wikipedia.org/wiki/Hydrophilicity_plot
**Accessed:** 2026-02-10
**Authority rank:** 3

**Key Extracted Points:**

1. **Source:** Kyte & Doolittle (1982). "A simple method for displaying the hydropathic character of a protein." J Mol Biol. 157(1):105-132.
2. **Scale values (all 20 standard amino acids):**

| Amino Acid | Code | Hydropathy |
|------------|------|------------|
| Isoleucine | I | 4.5 |
| Valine | V | 4.2 |
| Leucine | L | 3.8 |
| Phenylalanine | F | 2.8 |
| Cysteine | C | 2.5 |
| Methionine | M | 1.9 |
| Alanine | A | 1.8 |
| Glycine | G | -0.4 |
| Threonine | T | -0.7 |
| Serine | S | -0.8 |
| Tryptophan | W | -0.9 |
| Tyrosine | Y | -1.3 |
| Proline | P | -1.6 |
| Histidine | H | -3.2 |
| Glutamic acid | E | -3.5 |
| Glutamine | Q | -3.5 |
| Aspartic acid | D | -3.5 |
| Asparagine | N | -3.5 |
| Lysine | K | -3.9 |
| Arginine | R | -4.5 |

3. **Range:** Values span from -4.5 (most hydrophilic, Arg) to 4.5 (most hydrophobic, Ile).
4. **Window method:** Sliding window approach averages hydropathy over a window of consecutive residues.

### Wikipedia — Amino Acid (Properties Table)

**URL:** https://en.wikipedia.org/wiki/Amino_acid
**Accessed:** 2026-02-10
**Authority rank:** 3

**Key Extracted Points:**

1. **Confirming Kyte-Doolittle values:** The amino acid properties table (ref 66 = Kyte & Doolittle 1982) lists the same hydropathy scores as the Hydrophilicity plot article, cross-validating all 20 values.
2. **Charge at pH 7:** Three amino acids with side chains that are cations at neutral pH: Arginine (R, +1), Lysine (K, +1), and Histidine (H, pKa=6.0, ~10% protonated). Two amino acids are anions: Aspartate (D, -1), Glutamate (E, -1).
3. **Proline flexibility:** Proline's cyclic side chain makes it particularly inflexible when incorporated into proteins, disrupting alpha-helices — this is why it promotes disorder.
4. **Glycine flexibility:** Glycine's lack of side chain provides unique flexibility among amino acids — it provides conformational flexibility that favors disorder.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia (IDP article)

1. **Short sequences below window size:** When the input protein is shorter than the sliding window, boundary effects dominate and predictions are less reliable.
2. **Flanking ordered regions:** Disorder scores at the boundaries of ordered/disordered transitions may be unreliable due to window averaging effects.
3. **Unknown amino acids:** Non-standard amino acid codes (X, B, Z, etc.) are not part of the Kyte-Doolittle or disorder propensity scales.

### From Implementation Analysis

1. **Single-residue input:** Window of size 1 has no averaging context; score is based solely on that residue's properties.
2. **Homopolymeric sequences:** All-same-residue proteins (e.g., poly-P, poly-I) produce uniform per-residue TOP-IDP scores with no window-averaging benefit.

---

## Test Datasets

### Dataset: Kyte-Doolittle Hydropathy Scale

**Source:** Kyte & Doolittle (1982), confirmed by Wikipedia "Hydrophilicity plot" and "Amino acid" articles

| Amino Acid | Code | Hydropathy Value |
|------------|------|-----------------|
| A | A | 1.8 |
| R | R | -4.5 |
| N | N | -3.5 |
| D | D | -3.5 |
| C | C | 2.5 |
| Q | Q | -3.5 |
| E | E | -3.5 |
| G | G | -0.4 |
| H | H | -3.2 |
| I | I | 4.5 |
| L | L | 3.8 |
| K | K | -3.9 |
| M | M | 1.9 |
| F | F | 2.8 |
| P | P | -1.6 |
| S | S | -0.8 |
| T | T | -0.7 |
| W | W | -0.9 |
| Y | Y | -1.3 |
| V | V | 4.2 |

### Dataset: Charge at pH 7

**Source:** Wikipedia "Amino acid" article (side chain ionization properties)

| Amino Acid | Code | Charge |
|------------|------|--------|
| R | R | +1.0 |
| K | K | +1.0 |
| H | H | +0.1 (pKa 6.0, ~10% protonated) |
| D | D | -1.0 |
| E | E | -1.0 |
| All others | - | 0.0 |

### Dataset: Disorder-Promoting Classification

**Source:** Wikipedia IDP article, Dunker et al. (2001)

| Classification | Amino Acids |
|---------------|-------------|
| Disorder-promoting | A, R, G, Q, S, P, E, K |
| Order-promoting | W, C, F, I, Y, V, L, N |
| Neutral/borderline | D, T, H, M |

### Dataset: TOP-IDP Disorder Propensity Scale

**Source:** Campen et al. (2008) "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder." Protein Pept Lett 15(9):956-963. PMC2676888, PMID 18991772. Table 2.

| Amino Acid | Code | TOP-IDP Value | Classification |
|------------|------|--------------|----------------|
| W | W | -0.884 | Order-promoting |
| F | F | -0.697 | Order-promoting |
| Y | Y | -0.510 | Order-promoting |
| I | I | -0.486 | Order-promoting |
| M | M | -0.397 | Ambiguous |
| L | L | -0.326 | Order-promoting |
| V | V | -0.121 | Order-promoting |
| N | N | 0.007 | Order-promoting |
| C | C | 0.020 | Order-promoting |
| T | T | 0.059 | Ambiguous |
| A | A | 0.060 | Disorder-promoting |
| G | G | 0.166 | Disorder-promoting |
| R | R | 0.180 | Disorder-promoting |
| D | D | 0.192 | Ambiguous |
| H | H | 0.303 | Ambiguous |
| Q | Q | 0.318 | Disorder-promoting |
| S | S | 0.341 | Disorder-promoting |
| K | K | 0.586 | Disorder-promoting |
| E | E | 0.736 | Disorder-promoting |
| P | P | 0.987 | Disorder-promoting |

**Ranking (order→disorder):** W, F, Y, I, M, L, V, N, C, T, A, G, R, D, H, Q, K, S, E, P

---

## Assumptions

None. All parameters traceable to published peer-reviewed sources:
- Disorder propensity: TOP-IDP scale — Campen et al. (2008) Table 2.
- Prediction cutoff (0.542): Campen et al. (2008) maximum-likelihood.
- Classification: Dunker et al. (2001) disorder/order/ambiguous sets.
- Hydropathy: Kyte & Doolittle (1982).

---

## Recommendations for Test Coverage

1. **MUST Test:** TOP-IDP propensity values for all 20 amino acids match Campen et al. (2008) Table 2 — M8
2. **MUST Test:** Disorder-promoting classification matches Dunker et al. (2001): {A, R, G, Q, S, P, E, K} — M9
3. **MUST Test:** Order-promoting classification matches Dunker et al. (2001): {W, C, F, I, Y, V, L, N} — M10
4. **MUST Test:** Ambiguous amino acids {D, H, M, T} are NOT disorder-promoting — M10b
5. **MUST Test:** Hydrophobic sequences (poly-I) produce low disorder scores — Uversky et al. (2000) — M4
6. **MUST Test:** Charged/polar sequences (poly-E, poly-P) produce high disorder scores — Uversky et al. (2000) — M5, M6
7. **MUST Test:** Empty sequence returns zero-initialized result — M1
8. **MUST Test:** Residue predictions count equals sequence length — M2
9. **MUST Test:** All disorder scores are in [0, 1] range — M3
10. **MUST Test:** Case insensitivity — M7
11. **MUST Test:** Residue predictions have correct positions — M13
12. **SHOULD Test:** DisorderPromotingAminoAcids property contains all 8 Dunker disorder-promoting AA — M11
13. **SHOULD Test:** OrderPromotingAminoAcids property contains all 8 Dunker order-promoting AA — M12
14. **SHOULD Test:** CalculateHydropathy returns mean Kyte-Doolittle value — C4
15. **SHOULD Test:** AmbiguousAminoAcids property contains {D, H, M, T} — C3
16. **SHOULD Test:** Three classification sets are disjoint and cover all 20 AA — C5
17. **COULD Test:** Performance on long sequences — O(n) complexity claim

---

## References

1. Kyte J, Doolittle RF (1982). "A simple method for displaying the hydropathic character of a protein." J Mol Biol. 157(1):105-132. PMID 7108955.
2. Dunker AK, et al. (2001). "Intrinsically disordered protein." J Mol Graph Model. 19(1):26-59.
3. Uversky VN, Gillespie JR, Fink AL (2000). "Why are "natively unfolded" proteins unstructured under physiologic conditions?" Proteins. 41(3):415-427.
4. Campen A, et al. (2008). "TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder." Protein Pept Lett. 15(9):956-963.
5. Mohan A, et al. (2006). "Analysis of molecular recognition features (MoRFs)." J Mol Biol. 362(5):1043-1059.
6. Shannon CE (1948). "A Mathematical Theory of Communication." Bell System Technical Journal. 27(3):379-423.
7. Wikipedia. "Intrinsically disordered proteins." https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins
8. Wikipedia. "Hydrophilicity plot." https://en.wikipedia.org/wiki/Hydrophilicity_plot
9. Wikipedia. "Amino acid." https://en.wikipedia.org/wiki/Amino_acid

---

## Change History

- **2026-02-10**: Initial documentation.
- **2026-02-12**: Replaced 5-component scoring model with pure TOP-IDP averaging (Campen 2008). Changed threshold to published 0.542. Removed Charge dictionary (dead code). All assumptions resolved.
