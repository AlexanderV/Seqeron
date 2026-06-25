# Evidence Artifact: RNA-PKPREDICT-001

**Test Unit ID:** RNA-PKPREDICT-001
**Algorithm:** Pseudoknot Structure Prediction (canonical H-type, pknotsRG class)
**Date Collected:** 2026-06-23

---

## Online Sources

### Reeder & Giegerich 2004 — pknotsRG (BMC Bioinformatics 5:104)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC514697/ (open-access full text; canonical DOI 10.1186/1471-2105-5-104)
**Retrieved via:** WebSearch "Reeder Giegerich 2004 pknotsRG canonical simple recursive pseudoknots BMC Bioinformatics recurrence" → WebFetch of the PMC full text.
**Accessed:** 2026-06-23
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Class:** "canonical simple recursive pseudoknots" — a simple recursive pseudoknot is "two crossing helices with three intervening loops." The motif (their grammar) is
   `knot = knt <<< a ~~~ u ~~~ b ~~~ v ~~~ a' ~~~ w ~~~ b'`, with boundaries i, e, k, g, f, l, h, j; segment *a* pairs with *a'*, segment *b* pairs with *b'*; *u*, *v*, *w* are the three loops that separate the helices and can fold internally.
2. **Complexity:** O(n⁴) time, O(n²) space to predict the energetically optimal structure possibly containing such pseudoknots.
3. **Canonization rule 1:** "(a) Both strands in a helix must have the same length, i.e. |a| = |a'| and |b| = |b'|. (b) Both helices must not have bulges."
4. **Canonization rule 2:** "The helices a, a' and b, b' facing each other must have maximal extent" (loops between facing strands as short as base-pairing allows).
5. **Canonization rule 3:** "If two maximal helices would overlap, their boundary is fixed at an arbitrary point between them."
6. **Energy model:** Turner nearest-neighbour stacking for both pseudoknot helices, plus dangling/coaxial terms; pseudoknot-specific penalties below.
7. **Pseudoknot initiation penalty:** "the pseudoknot initiation parameter … we found out, that setting this value to **9 kcal/mole** performs better."
8. **Unpaired loop penalty:** "we penalize each unpaired nucleotide inside a pseudoknot loop with **0.3 kcal/mole**."

### pknotsRG reference source (Energy.lhs)

**URL:** https://github.com/jensreeder/pknotsRG
**Retrieved via:** WebSearch → WebFetch of the repository page (Energy.lhs / Foldingspace.lhs).
**Accessed:** 2026-06-23
**Authority rank:** 3 (reference implementation by the paper's authors)

**Key Extracted Points:**

1. Pseudoknot destabilizing penalties, verbatim from the documentation: **"creating a new pseudoknot: 9.0"**, **"not paired base in pk: 0.3"**, **"basepair inside pseudoknot: 0.0"** (kcal/mol).
2. Stacking interactions use mfold-3.1 / Turner values; dangling and coaxial stacking use the standard nested-structure energies. This confirms the two helices are scored with the SAME nearest-neighbour model as nested structures, and there is **no** extra per-base-pair penalty inside the pseudoknot.

### Wikipedia — Pseudoknot (H-type geometry; cites Rivas & Eddy 1999)

**URL:** https://en.wikipedia.org/wiki/Pseudoknot
**Retrieved via:** WebFetch.
**Accessed:** 2026-06-23
**Authority rank:** 4 (uses cited primaries)

**Key Extracted Points:**

1. **H-type 5'→3' order:** stem1-5' → loop1 → stem2-5' → loop2 → stem1-3' → loop3 → stem2-3'; "half of one stem is intercalated between the two halves of another stem," producing the crossing (pseudoknotted) arrangement.
2. Telomerase P2b-P3 is given as an H-type pseudoknot example with a two-layer `()`/`[]` annotation.

### RCSB PDB 437D — BWYV frameshifting pseudoknot (Su et al. 1999)

**URL:** https://www.rcsb.org/structure/437D ; structure paper https://pmc.ncbi.nlm.nih.gov/articles/PMC7097825/
**Retrieved via:** WebSearch "BWYV pseudoknot 437D … sequence" → WebFetch of the RCSB entry and PMC paper.
**Accessed:** 2026-06-23
**Authority rank:** 1–5 (structural database + peer-reviewed paper)

**Key Extracted Points:**

1. The beet western yellows virus (BWYV) frameshifting pseudoknot is a documented 28-nt **H-type** pseudoknot; PDB 437D chain sequence is **`GGCGCGGCACCGUCCGCGGAACAAACGG`** (28 nt).
2. It has two stems S1/S2 and is stabilized in the crystal by tertiary interactions (a minor-groove triplex of loop 2 with stem 1, a quadruple interaction of loop 1 with stem 2, Mg²⁺/Na⁺ coordination) — interactions OUTSIDE the secondary-structure nearest-neighbour model.

---

## Documented Corner Cases and Failure Modes

### From Reeder & Giegerich 2004

1. **Spurious pseudoknots:** the 9 kcal/mol initiation penalty exists precisely so the optimal fold contains a pseudoknot only when the two crossing helices' stabilizing energy outweighs both the penalty and the best pseudoknot-free alternative. A pure thermodynamic predictor must therefore NOT report a pseudoknot when a competing nested structure is more stable.
2. **Tertiary-stabilized knots:** knots like BWYV are held together by triplexes / ion coordination not captured by the secondary-structure energy model; a thermodynamic NN predictor is not expected to recover them as the MFE structure (a documented limit of all NN-only pseudoknot predictors).

### From the canonization rules

3. **Equal-length, bulge-free helices (rule 1):** within the canonical class a helix is a contiguous run of pairs of equal strand length; bulged helices are not in the class.
4. **Maximal extent (rule 2):** facing helices are extended to maximal Watson–Crick/GU length, so the predicted helix length is determined by the sequence, not searched independently.

---

## Test Datasets

### Dataset: Designed canonical H-type pseudoknot (fully derivable)

**Source:** Constructed from the H-type geometry (Wikipedia / Reeder & Giegerich 2004) so that both crossing helices are strong and unambiguous.

| Parameter | Value |
|-----------|-------|
| Sequence | `GGGGAACCCCAACCCCAAGGGG` (22 nt) |
| Layout | S1a=GGGG[0–3] · L1=AA · S2a=CCCC[6–9] · L2=AA · S1b=CCCC[12–15] · L3=AA · S2b=GGGG[18–21] |
| Stem 1 (a·a') | (0,15)(1,14)(2,13)(3,12) — four G·C pairs |
| Stem 2 (b·b') | (6,21)(7,20)(8,19)(9,18) — four C·G pairs |
| Crossing check | S2a (6–9) lies inside S1's span (0–15); S2b (18–21) lies outside → i<k<j<l crossing |
| Two-layer dot-bracket | `((((..[[[[..))))..]]]]` |
| Pseudoknot ΔG | strictly below the plain-MFE ΔG of the same sequence (knot accepted) |

### Dataset: BWYV pseudoknot (real H-type; thermodynamic non-recovery documented)

**Source:** PDB 437D / Su et al. 1999.

| Parameter | Value |
|-----------|-------|
| Sequence | `GGCGCGGCACCGUCCGCGGAACAAACGG` (28 nt) |
| Observed | The NN-thermodynamic optimum is the pseudoknot-free stem-1 hairpin; the crystallographic knot is tertiary-stabilized and not the MFE structure (expected; documents the limit). |

### Dataset: Plain hairpin (no spurious pseudoknot)

**Source:** Trivial well-formed hairpin.

| Parameter | Value |
|-----------|-------|
| Sequence | `GGGGAAAACCCC` |
| Expected | `HasPseudoknot == false`; returned structure and ΔG equal the plain MFE (`((((....))))`). |

---

## Assumptions

1. **ASSUMPTION: PARTIAL coverage of the pknotsRG class.** The implementation realizes the canonical *single* H-type pseudoknot (two crossing helices + three internally-foldable loops) with the sourced energy model and penalties. The full pknotsRG O(n⁴) grammar additionally composes recursively-nested pseudoknots and over-arching/multiple knots within one structure; these are NOT implemented. Loops u/v/w fold with the existing pseudoknot-free MFE (`CalculateMinimumFreeEnergy`), which is consistent with "loops can fold internally" but does not re-search a second knot inside a loop. This is documented (PARTIAL), not an invented parameter.

---

## Recommendations for Test Coverage

1. **MUST Test:** Designed canonical H-type sequence → `HasPseudoknot==true`, the two crossing 4-bp helices recovered exactly, two-layer dot-bracket `((((..[[[[..))))..]]]]`, ΔG strictly below MFE — Evidence: H-type geometry + pknotsRG penalties.
2. **MUST Test:** Plain hairpin / non-knot sequence → `HasPseudoknot==false`, structure and ΔG equal MFE (no spurious pseudoknot) — Evidence: 9 kcal/mol penalty rationale, Reeder & Giegerich 2004.
3. **MUST Test (invariant):** for any sequence, `FreeEnergy ≤ CalculateMfeStructure(seq).FreeEnergy` (the predictor never returns a structure worse than the plain MFE) — Evidence: the MFE is the always-available fallback baseline.
4. **MUST Test (validity):** when a pseudoknot is returned, every index is in range and each position is paired at most once, and `DetectPseudoknots` finds ≥1 genuine crossing — Evidence: crossing condition i<k<j<l (Antczak 2018), canonization rule 1.
5. **SHOULD Test:** null / empty / too-short input → empty pseudoknot-free structure (no pairs, all dots, ΔG 0) — Rationale: contract parity with `CalculateMfeStructure`.
6. **SHOULD Test:** BWYV real knot is NOT recovered as MFE (documents the NN-thermodynamic limit) — Rationale: prevents over-claiming; matches the literature.
7. **COULD Test:** DNA input (T read as U) folds identically to the RNA spelling — Rationale: parity with the rest of the module.

---

## References

1. Reeder J, Giegerich R. (2004). Design, implementation and evaluation of a practical pseudoknot folding algorithm based on thermodynamics. BMC Bioinformatics 5:104. https://doi.org/10.1186/1471-2105-5-104 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC514697/)
2. Reeder J. pknotsRG source (Energy.lhs — pseudoknot penalties). https://github.com/jensreeder/pknotsRG
3. Pseudoknot. Wikipedia (cites Rivas & Eddy 1999, Staple & Butcher 2005). https://en.wikipedia.org/wiki/Pseudoknot
4. Su L, Chen L, Egli M, Berger JM, Rich A. (1999). Minor groove RNA triplex in the crystal structure of a ribosomal frameshifting viral pseudoknot. Nat Struct Biol 6(3):285–292. https://pmc.ncbi.nlm.nih.gov/articles/PMC7097825/ ; PDB 437D https://www.rcsb.org/structure/437D
5. Antczak M, et al. (2018). New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation. Bioinformatics 34(8):1304–1312 (crossing condition i<k<j<l). https://doi.org/10.1093/bioinformatics/btx783

---

## Change History

- **2026-06-23**: Initial documentation (RNA-STRUCT-001 limitation fix: add canonical H-type pseudoknot prediction).
