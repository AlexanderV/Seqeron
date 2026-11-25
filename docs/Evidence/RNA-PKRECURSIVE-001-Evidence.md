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

## Review 2026-11-25 (Stage A / Stage B)

**Stage A.** The previous `RecursivePkFolder` was not the pknotsRG recursion: (1) an anchored interval memo calling `CalculateMfeStructure` on every interval (6.9 s at 80 nt, three predictors);
(2) a third production "enclosing helix a·a' over a knotted interior" that is not in the grammar of `pknotsRG-mfe.lhs` (a closed stem encloses a knot only through a multiloop with ≥ 2 components);
(3) no knot as a multiloop branch (`pkml` = +600, Algebras.lhs `energyalgext`) and no knot component in knot loops with its own cost; (4) the same energy omissions as RNA-PKPREDICT-001 (flanks, bridging nucleotides, helix-end terms);
(5) a duplicate knot scorer (`EvaluateHTypeRecursive`, `ScoreLoopRecursive`, private `MaxHelixLength`). The old "nested knot under an A·U helix" datasets therefore tested a structure class that is **not** pknotsRG's; they are replaced.

**Stage B.** Same oracle, same corpus and same statistics as RNA-PKPREDICT-001 (Evidence §Stage B): 0 energy differences in the recursive mode over 1 971 comparable sequences (+2 multiloop-knot cases of 62–64 nt), 2 NNDB-bulge explained.
Recursive strictly below single-knot for 98 of 809 designed sequences (95 comparable, all agree). Knot-in-multiloop cases occur in 112 of 187 designed multiloop sequences (long GC stems; ML branch 6.0 + closing 9.3 still pays off).

## Test Datasets (oracle values)

| Case | Sequence | single | recursive | recursive structure |
|------|----------|--------|-----------|---------------------|
| two knots side by side | `GGGGAACCCCAACCCCAAGGGGAAAGGGGAACCCCAACCCCAAGGGG` | −24.46 | −27.52 | `((((..[[[[..))))..]]]]...((((..[[[[..))))..]]]]` |
| two knots, different size | `CCGCUUGGAGUGCGGAAACUCCUUUCGGAAGCCGUCCGUACGGC` | −18.85 | −22.38 | `((((..[[[[[))))..]]]]]...((((.[[[[))))..]]]]` |
| knot inside loop v | `GCCAGCGGCCUUACGGGAUAGGCAUACCCGGGCAACGC` | −11.56 | −12.15 | `(((.[[[((((..[[[[..))))...]]]])))..]]]` |
| knot as multiloop branch | `CGCCCGGGCAAGGCGGGAAACCGCCAAGCGGGAGCCCCCCCCCGCCUGGGGGCAAGCCCGGGCG` | −43.46 | −46.66 | `(((((((((..(((((....)))))..(((((.[[[[[[.)))))..]]]]]]..)))))))))` |
| no knot pays off | `GGGGAAAACCCC` | −5.28 | −5.28 | `((((....))))` |

Plus the 27 rows of `PknotsRgOracleData` and the exhaustive equality "recursive = single-knot optimum for 11–16 nt" (no second/nested knot fits).

## Assumptions

1. As RNA-PKPREDICT-001 (pknotsRG grammar and pseudoknot parameters over the shared Turner 2004 model; strict csr-PK without the relaxed 1-nt bulge).
2. A knot used as a component of a knot loop is charged like a multiloop branch (+6.0), following `block` / `pkml` of `pknotsRG-mfe.lhs`, which uses the same `comps` for multiloops and knot loops.

---

## References

1. Reeder J, Giegerich R. (2004). Design, implementation and evaluation of a practical pseudoknot folding algorithm based on thermodynamics. BMC Bioinformatics 5:104. https://doi.org/10.1186/1471-2105-5-104 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC514697/)
2. Reeder J, Steffen P, Giegerich R. (2007). pknotsRG: RNA pseudoknot folding including near-optimal structures and sliding windows. Nucleic Acids Research 35:W320–W324. https://doi.org/10.1093/nar/gkm258 (full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC1933184/)
3. Reeder J. pknotsRG source (Energy.lhs — pseudoknot penalties 9.0 / 0.3 / 0.0). https://github.com/jensreeder/pknotsRG
4. Antczak M, et al. (2018). New algorithms to represent complex pseudoknotted RNA structures in dot-bracket notation. Bioinformatics 34(8):1304–1312 (crossing condition i<k<j<l). https://doi.org/10.1093/bioinformatics/btx783

---

## Change History

- **2026-06-23**: Initial documentation (RNA-STRUCT-001 limitation fix: recursive pknotsRG grammar — nested / multiple / over-arching pseudoknots).
- **2026-11-25**: Stage A/B review: previous interval-memo recursion (non-grammar production, no multiloop knots, duplicate scorer) replaced by the exact pknotsRG DP shared with RNA-PKPREDICT-001; datasets and assumptions rewritten.
