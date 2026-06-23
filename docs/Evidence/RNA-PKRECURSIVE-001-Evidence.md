# Evidence Artifact: RNA-PKRECURSIVE-001

**Test Unit ID:** RNA-PKRECURSIVE-001 (recursive-grammar extension of RNA-STRUCT-001 / RNA-PKPREDICT-001)
**Algorithm:** Recursive pknotsRG pseudoknot prediction (nested / multiple H-type pseudoknots)
**Date Collected:** 2026-06-23

---

## Online Sources

### Reeder & Giegerich 2004 — pknotsRG (BMC Bioinformatics 5:104)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC514697/ (open-access full text; canonical DOI 10.1186/1471-2105-5-104)
**Retrieved via:** WebSearch "Reeder Giegerich 2004 pknotsRG pseudoknot folding algorithm canonical recursive grammar BMC Bioinformatics" → WebFetch of the PMC full text (two passes, focused on the ADP grammar and recursion).
**Accessed:** 2026-06-23
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Class — simple recursive pseudoknot:** "a crosswise interaction of two helices" with three intervening loops u, v, w. The motif (their grammar) is `knot = knt <<< a ~~~ u ~~~ b ~~~ v ~~~ a' ~~~ w ~~~ b'`.
2. **Recursion (the core of this extension), verbatim:** "In simple recursive pseudoknots, we allow the unpaired strands u, v, w in a simple pseudoknot to fold internally in an arbitrary way, **including simple recursive pseudoknots**." So the three loops may themselves contain further pseudoknots — this is what `PredictStructurePseudoknot` (single, loops folded pseudoknot-free) does NOT do and this method does.
3. **Complexity:** O(n⁴) time, O(n²) space to predict the energetically optimal structure possibly containing such pseudoknots; the canonization reduces independent boundaries from 8 to 4 versus Rivas & Eddy's O(n⁶)/O(n⁴).
4. **Canonization rule 1:** "(a) Both strands in a helix must have the same length … (b) Both helices must not have bulges."
5. **Canonization rule 2:** "The helices a, a' and b, b' facing each other must have maximal extent" (compartment v as short as base pairing allows).
6. **Canonization rule 3:** "If two maximal helices would overlap, their boundary is fixed at an arbitrary point between them."
7. **Energy model:** Turner nearest-neighbour stacking for both pseudoknot helices, plus dangling/coaxial terms; pseudoknot-specific penalties: initiation "setting this value to **9 kcal/mole** performs better"; "we penalize each unpaired nucleotide inside a pseudoknot loop with **0.3 kcal/mole**."
8. **Excluded classes, verbatim:** "More complex knotted structures like triple crossing helices or kissing hairpins … are excluded from sr-PK." Chained / complex helix interactions are out of class.

### Reeder & Giegerich 2007 — pknotsRG (Nucleic Acids Research 35:W320)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC1933184/ (open-access full text; DOI 10.1093/nar/gkm258)
**Retrieved via:** WebSearch "pknotsRG ADP grammar pkiss canonical simple recursive pseudoknot knot nonterminal" → WebFetch of the PMC full text.
**Accessed:** 2026-06-23
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Whole-sequence DP / competition:** the algorithm "extends the usual dynamic programming (DP) scheme for RNA folding"; the pseudoknot matrix value "**competes with values of unknotted foldings for the interval (i, j)**." Folding the whole sequence with this per-interval competition is what allows the optimal structure to contain **several** pseudoknots in different regions and pseudoknots nested inside loops.
2. **Recursive class restated, verbatim:** "If the unpaired strands (u,v,w) build secondary structures internally in an arbitrary way, including multiloops and pseudoknots, we call this class simple recursive pseudoknots."
3. **Energy model:** "uses the current thermodynamic energy model by the Turner group, extended by some pseudoknot specific values."

### pknotsRG reference source (jensreeder/pknotsRG, Energy.lhs)

**URL:** https://github.com/jensreeder/pknotsRG
**Retrieved via:** WebSearch (paper search results listed the repo) — reference implementation by the paper's authors.
**Accessed:** 2026-06-23
**Authority rank:** 3 (reference implementation by the paper's authors)

**Key Extracted Points:**

1. Pseudoknot destabilizing penalties (as recorded in the prior RNA-PKPREDICT-001 evidence and re-used unchanged): **"creating a new pseudoknot: 9.0"**, **"not paired base in pk: 0.3"**, **"basepair inside pseudoknot: 0.0"** (kcal/mol). Confirms the two helices are scored with the SAME nearest-neighbour model as nested structures and there is no extra per-base-pair penalty inside the knot.

### Antczak et al. 2018 — crossing condition for dot-bracket

**URL:** https://doi.org/10.1093/bioinformatics/btx783
**Retrieved via:** cited from the prior RNA-PKPREDICT-001 evidence (the crossing test i<k<j<l used by `DetectPseudoknots`).
**Accessed:** 2026-06-23
**Authority rank:** 1

**Key Extracted Points:**

1. Two base pairs (i,j) and (k,l) cross iff i<k<j<l; the two-layer `()` / `[]` annotation places the crossing helix on the second bracket family.

---

## Documented Corner Cases and Failure Modes

### From Reeder & Giegerich 2004 / 2007

1. **Spurious pseudoknots:** the 9 kcal/mol initiation penalty exists precisely so a pseudoknot is taken only when its two crossing helices outweigh both the penalty and the best pseudoknot-free alternative *for that interval*. A pure thermodynamic predictor must NOT report a knot where a competing nested structure is more stable.
2. **Two simultaneous strong G·C knots are thermodynamically degenerate:** two H-type knots need strong (G·C) stems, and the same G·C content almost always offers a more stable *cross-region nested* helix arrangement. Engineering an instance where two side-by-side knots are the genuine MFE therefore requires **isolating** each knot region (e.g. flanking A·U clamps) so the cross-region nested alternative is suppressed. This is a property of the energy model, documented here so the test cases are honest.
3. **Tertiary-stabilised knots (e.g. BWYV / PDB 437D):** held by minor-groove triplexes / ion coordination outside the nearest-neighbour secondary-structure model; not recoverable as the MFE by any NN-only predictor — an energy-model floor, not an algorithm gap.
4. **Excluded classes:** kissing hairpins, triple-crossing / chained ("complex") helix interactions, and bulged/unequal-length pseudoknot helices (rule 1) are NOT in the canonical csr-PK class and are not predicted.

---

## Test Datasets

### Dataset: Knot nested inside an outer (over-arching) helix — fully derivable

**Source:** Constructed from the H-type geometry + an isolating A·U clamp so the outer helix over-arches the knot in its loop. The inner knot is the verified designed H-type `GGGGAACCCCAACCCCAAGGGG`.

| Parameter | Value |
|-----------|-------|
| Sequence | `AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU` (38 nt) |
| Layout | clamp5'=A×8[0–7] · H-type[8–29] · clamp3'=U×8[30–37] |
| Recursive structure | `((((((((((((..[[[[..))))..]]]]))))))))` |
| Outer A·U helix | (0,37)(1,36)(2,35)(3,34)(4,33)(5,32)(6,31)(7,30) — 8 A·U pairs |
| Inner stem 1 (a·a') | (8,23)(9,22)(10,21)(11,20) — 4 G·C pairs |
| Inner stem 2 (b·b', crossing) | (14,29)(15,28)(16,27)(17,26) — 4 C·G pairs |
| Recursive ΔG | **−14.37 kcal/mol**, `HasPseudoknot = true`, 1 crossing |
| Single-knot method (`PredictStructurePseudoknot`) | `......((((((((((((........)))))))))))) `, **−13.05**, `HasPseudoknot = false` (cannot combine outer helix + inner knot) |
| Plain MFE | identical to single: **−13.05**, no knot |
| Verification | recursive ΔG (−14.37) < single/MFE (−13.05): the over-arching knot is recovered ONLY by the recursive method |

### Dataset: Two separate (non-nested) H-type knots — fully derivable

**Source:** Two A·U-clamped copies of the designed H-type, the clamps isolating each region so the cross-region nested alternative is suppressed (see Failure Mode 2).

| Parameter | Value |
|-----------|-------|
| Sequence | `AAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUUAAAAAAAAAAAAGGGGAACCCCAACCCCAAGGGGUUUUUUUU` (80 nt) |
| Layout | clamped-knot[0–37] · A×12 linker[38–49 region overlaps clamp] · clamped-knot[38–79] |
| Recursive structure | `((((((((((((..[[[[..))))..]]]]))))))))((((((((....((((..[[[[..))))..]]]]))))))))` |
| Recovered knots | **two** crossing H-type knots (DetectPseudoknots crossing-count = 32) |
| Recursive ΔG | **−28.74 kcal/mol**, `HasPseudoknot = true` |
| Single-knot method | identical to MFE: **−27.14**, `HasPseudoknot = false` (recovers neither knot) |
| Plain MFE | **−27.14**, no knot |
| Verification | recursive ΔG (−28.74) < single/MFE (−27.14); both knots recovered, single-knot method recovers none |

### Dataset: Plain (non-pseudoknotted) sequences — no spurious knots

**Source:** trivial hairpin and a weak A·U run.

| Parameter | Value |
|-----------|-------|
| Hairpin `GGGGAAAACCCC` | recursive = MFE `((((....))))`, **−5.28**, `HasPseudoknot = false`, crossing = 0 |
| A·U run `AUAUAUAUAUAUAUAU` | recursive = MFE `((((((....))))))`, **−0.26**, no knot, crossing = 0 |
| Random sweep (seed 20260623, 150 seqs, len 12–38) | recursive ΔG ≤ MFE ΔG for all (0 violations); 0 spurious knots reported on random short inputs |

### Dataset: Single canonical H-type (parity with `PredictStructurePseudoknot`)

| Parameter | Value |
|-----------|-------|
| Sequence | `GGGGAACCCCAACCCCAAGGGG` (22 nt) |
| Both methods | `((((..[[[[..))))..]]]]`, **−8.76 kcal/mol**, identical pairs, `HasPseudoknot = true` |

---

## Assumptions

1. **ASSUMPTION: PARTIAL coverage of the full pknotsRG recursion.** This method realises the recursive grammar for the canonical csr-PK class: (i) the three loops u, v, w of a knot fold by the same recursive folder (so a loop may contain a further knot), (ii) the top level chains multiple knots, and (iii) an enclosing helix may over-arch a knot in its loop. To keep within a tractable envelope the H-type helices are enumerated by explicit start/end scan with maximal extension (rules 1–2) rather than the full 4-boundary ADP yield-parser, and a knot component is left-anchored within the interval before chaining; the enclosing-helix production is pursued only when the enclosed region is itself knotted (otherwise the pseudoknot-free Zuker MFE already covers it). The result is the faithful recursive *class* (nested / multiple / over-arching knots) but not a guaranteed bit-identical reproduction of every yield the reference O(n⁴) parser would explore. Documented PARTIAL; no invented parameter.

2. **ASSUMPTION: two-simultaneous-knot test cases are engineered, not random.** Per Failure Mode 2, two strong knots are the genuine MFE only when each knot region is isolated; the test asserts the recovery on the *engineered* isolated-clamp sequence, not a universal "beats single-knot on random input" (which would be false thermodynamics).

---

## Recommendations for Test Coverage

1. **MUST Test:** the over-arching nested-knot sequence → recursive recovers the outer helix + inner crossing knot exactly (pairs as in the dataset), `HasPseudoknot == true`, ΔG = −14.37, strictly below the single-knot method / MFE (−13.05) — Evidence: recursion into loops (Reeder & Giegerich 2004 §recursive class).
2. **MUST Test:** the two-knot sequence → recursive recovers TWO crossing knots (crossing-count = 32), ΔG = −28.74 < single/MFE (−27.14); the single-knot method recovers none — Evidence: per-interval competition (Reeder & Giegerich 2007).
3. **MUST Test (no spurious knots):** plain hairpin / A·U run / random sweep → recursive ΔG ≤ MFE, `HasPseudoknot == false`, no crossing pairs — Evidence: 9 kcal/mol penalty rationale.
4. **MUST Test (invariant):** recursive ΔG ≤ `CalculateMfeStructure` ΔG for any sequence (MFE is the always-available fallback) — Evidence: fallback baseline.
5. **MUST Test (validity):** every index in range, each position paired ≤ once; reported knots have ≥1 genuine crossing — Evidence: crossing condition i<k<j<l (Antczak 2018), rule 1.
6. **SHOULD Test:** null / empty / too-short → empty pseudoknot-free structure (parity with `PredictStructurePseudoknot`).
7. **SHOULD Test:** single canonical H-type → identical result to `PredictStructurePseudoknot` (recursion does not regress the single-knot case).
8. **COULD Test:** DNA spelling (T read as U) folds identically; minLoopSize < 3 clamps to 3.

---

## References

1. Reeder J, Giegerich R. (2004). Design, implementation and evaluation of a practical pseudoknot folding algorithm based on thermodynamics. BMC Bioinformatics 5:104. https://doi.org/10.1186/1471-2105-5-104 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC514697/)
2. Reeder J, Steffen P, Giegerich R. (2007). pknotsRG: RNA pseudoknot folding including near-optimal structures and sliding windows. Nucleic Acids Research 35:W320–W324. https://doi.org/10.1093/nar/gkm258 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC1933184/)
3. Reeder J. pknotsRG source (Energy.lhs — pseudoknot penalties 9.0 / 0.3 / 0.0). https://github.com/jensreeder/pknotsRG
4. Antczak M, et al. (2018). New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation. Bioinformatics 34(8):1304–1312 (crossing condition i<k<j<l). https://doi.org/10.1093/bioinformatics/btx783

---

## Change History

- **2026-06-23**: Initial documentation (RNA-STRUCT-001 limitation fix: recursive pknotsRG grammar — nested / multiple / over-arching pseudoknots).
