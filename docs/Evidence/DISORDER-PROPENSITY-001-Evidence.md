# Evidence Artifact: DISORDER-PROPENSITY-001

**Test Unit ID:** DISORDER-PROPENSITY-001
**Algorithm:** Disorder Propensity (TOP-IDP scale lookup + Dunker order/disorder amino-acid classification)
**Date Collected:** 2026-06-14

---

## Online Sources

### Campen et al. (2008) — TOP-IDP-Scale (PMC full text)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *Protein & Peptide Letters* 15(9):956-963)

**Retrieval:** WebSearch query `Campen 2008 TOP-IDP scale amino acid intrinsic disorder propensity Table values` → opened the PMC2676888 full text with WebFetch, prompted for the exact per-amino-acid TOP-IDP values from Table 2, the ranking order, and the prediction cutoff.

**Key Extracted Points:**

1. **TOP-IDP per-residue values (Table 2), quoted verbatim from the fetched text:**
   W: -0.884, F: -0.697, Y: -0.510, I: -0.486, M: -0.397, L: -0.326, V: -0.121,
   N: 0.007, C: 0.02, T: 0.059, A: 0.06, G: 0.166, R: 0.180, D: 0.192,
   H: 0.303, Q: 0.318, S: 0.341, K: 0.586, E: 0.736, P: 0.987.
2. **Ranking (order-promoting → disorder-promoting), quoted:** "W, F, Y, I, M, L, V, N, C, T, A, G, R, D, H, Q, K, S, E, P". (Note: the text-rendered ranking places K before S, but the numeric values give S=0.341 < K=0.586; the per-residue *values* are authoritative and used for testing — see Assumptions.)
3. **Prediction cutoff, quoted:** "a prediction cut-off of 0.542 was calculated"; formula `I_Top-IDP = -(<Top-IDP> - 0.542)`. (Cutoff governs `PredictDisorder`, not the four methods in this unit's scope; recorded for completeness.)
4. **Most order/disorder anchors:** W (-0.884) is the most order-promoting residue; P (+0.987) is the most disorder-promoting residue.

### Wikipedia — Intrinsically disordered proteins (cites Dunker et al. 2001 primary)

**URL:** https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia article; used only for the Dunker et al. 2001 primary classification it cites)

**Retrieval:** WebSearch query `Dunker 2001 "Intrinsically disordered protein" disorder-promoting order-promoting amino acids ...` → opened the Wikipedia article with WebFetch, prompted for the amino-acid classification with the cited primary (Dunker et al.).

**Key Extracted Points:**

1. **Disorder-promoting amino acids, quoted:** "A, R, G, Q, S, P, E and K" (8 residues).
2. **Order-promoting amino acids, quoted:** "W, C, F, I, Y, V, L, and N" (8 residues).
3. **Ambiguous amino acids, quoted:** "H, M, T and D" (4 residues).
4. **Continuum ranking, quoted:** "W, F, Y, I, M, L, V, N, C, T, A, G, R, D, H, Q, K, S, E, P" (agrees with Campen 2008 ordering, confirming the cross-source consistency of the scale).
5. **Property note, quoted/paraphrased:** disorder-promoting residues are hydrophilic and charged; order-promoting residues are hydrophobic and uncharged; the ambiguous group occurs in both ordered and unstructured regions.

### Dunker et al. (2001) — PubMed record (primary citation locator)

**URL:** https://pubmed.ncbi.nlm.nih.gov/11381529/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed paper, *J Mol Graph Model* 19(1):26-59)

**Retrieval:** Returned by the WebSearch query above; used to confirm the primary citation (author/year/journal/PMID) that the Wikipedia classification cites. The verbatim residue sets are extracted from the Wikipedia source (cited to this primary), not recalled.

**Key Extracted Points:**

1. **Citation:** Dunker AK, Lawson JD, Brown CJ, et al. (2001) "Intrinsically disordered protein." *J Mol Graph Model* 19(1):26-59. PMID 11381529.

---

## Documented Corner Cases and Failure Modes

### From Campen et al. (2008)

1. **Scale defined for the 20 standard amino acids only:** Table 2 lists values for the 20 standard residues. No value is defined for non-standard / ambiguity codes (B, J, O, U, X, Z) or gap characters.

### From the implementation contract (GetValueOrDefault)

1. **Unknown residue → 0.0:** `GetDisorderPropensity` returns 0.0 for any character not in the 20-residue scale (lookup default). This is an implementation contract, not a value defined by the source — recorded as an assumption.
2. **Case handling:** input is upper-cased before lookup, so 'p' and 'P' return the same value.

---

## Test Datasets

### Dataset: TOP-IDP scale (Campen et al. 2008, Table 2) — all 20 residues

**Source:** Campen et al. (2008) PMC2676888, Table 2 (fetched 2026-06-14).

| Amino acid | TOP-IDP value |
|-----------|---------------|
| A | 0.060 |
| R | 0.180 |
| N | 0.007 |
| D | 0.192 |
| C | 0.020 |
| Q | 0.318 |
| E | 0.736 |
| G | 0.166 |
| H | 0.303 |
| I | -0.486 |
| L | -0.326 |
| K | 0.586 |
| M | -0.397 |
| F | -0.697 |
| P | 0.987 |
| S | 0.341 |
| T | 0.059 |
| W | -0.884 |
| Y | -0.510 |
| V | -0.121 |

### Dataset: Dunker et al. (2001) order/disorder classification

**Source:** Wikipedia "Intrinsically disordered proteins" citing Dunker et al. (2001); fetched 2026-06-14.

| Class | Amino acids | Count |
|-------|-------------|-------|
| Disorder-promoting | A, R, G, Q, S, P, E, K | 8 |
| Order-promoting | W, C, F, I, Y, V, L, N | 8 |
| Ambiguous | H, M, T, D | 4 |

---

## Assumptions

1. **ASSUMPTION: Unknown-residue propensity = 0.0** — Campen (2008) defines values only for the 20 standard residues; returning 0.0 for any out-of-scale character is an implementation contract (`GetValueOrDefault(..., 0)`), not a source-defined value. Tested as a documented contract, not against the source.
2. **ASSUMPTION: Ranking-vs-value discrepancy for S/K** — the Campen (2008) and Wikipedia rendered ranking strings place "...Q, K, S, E, P", whereas the per-residue Table 2 values give S=0.341 < K=0.586 (so by value the order is Q, S, K, E, P). The numeric Table 2 values are authoritative and are what the implementation and tests use; the ranking string is a presentation-order artifact. No correctness impact on the four scope methods (which use values and set membership, not the rank string).

---

## Recommendations for Test Coverage

1. **MUST Test:** `GetDisorderPropensity` returns the exact Table 2 value for all 20 standard residues — Evidence: Campen et al. (2008) Table 2.
2. **MUST Test:** `IsDisorderPromoting` is true for each of {A, R, G, Q, S, P, E, K} and false for each order-promoting {W, C, F, I, Y, V, L, N} and ambiguous {H, M, T, D} residue — Evidence: Dunker et al. (2001).
3. **MUST Test:** `DisorderPromotingAminoAcids` = {A, E, G, K, P, Q, R, S} and `OrderPromotingAminoAcids` = {C, F, I, L, N, V, W, Y}, each exactly 8 members — Evidence: Dunker et al. (2001).
4. **MUST Test:** the three classification sets are pairwise disjoint and cover all 20 standard residues (8+8+4) — Evidence: Dunker et al. (2001) (derived consistency check).
5. **SHOULD Test:** unknown residue → 0.0; lowercase input equals uppercase — Rationale: documented implementation contract / case-insensitivity.
6. **COULD Test:** W is the global minimum (-0.884) and P the global maximum (0.987) of `GetDisorderPropensity` over the 20 residues — Rationale: confirms anchor residues from Campen Table 2.

---

## References

1. Campen A, Williams RM, Brown CJ, Meng J, Uversky VN, Dunker AK (2008). TOP-IDP-Scale: A New Amino Acid Scale Measuring Propensity for Intrinsic Disorder. Protein Pept Lett 15(9):956-963. https://pmc.ncbi.nlm.nih.gov/articles/PMC2676888/ (PMID 18991772)
2. Dunker AK, Lawson JD, Brown CJ, et al. (2001). Intrinsically disordered protein. J Mol Graph Model 19(1):26-59. https://pubmed.ncbi.nlm.nih.gov/11381529/ (PMID 11381529)
3. Wikipedia contributors. Intrinsically disordered proteins. https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins (accessed 2026-06-14; used for the Dunker 2001 classification it cites)

---

## Change History

- **2026-06-14**: Initial documentation.
