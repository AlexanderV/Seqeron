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
2. **Tertiary-stabilized knots:** knots like BWYV are held together by triplexes / ion coordination not captured by the secondary-structure energy model; the NN model optimum need not equal the crystal pairing (the model optimum is reported; it was not compared with the crystal pairing).

### From the canonization rules

3. **Equal-length, bulge-free helices (rule 1):** within the canonical class a helix is a contiguous run of pairs of equal strand length; bulged helices are not in the class.
4. **Maximal extent (rule 2):** facing helices are extended to maximal Watson–Crick/GU length, so the predicted helix length is determined by the sequence, not searched independently.

---

## Review 2026-11-25 (Stage A / Stage B) — sources opened and what changed

### Sources opened (raw files, 2026-11-25)

| Source | What was read | Used for |
|--------|---------------|----------|
| pknotsRG `README`, `Foldingspace.lhs`, `Energy.lhs`, `Algebras.lhs`, `pknotsRG-mfe.lhs` (raw.githubusercontent.com/jensreeder/pknotsRG/master) | grammar `struct` / `edangle` / `block` / `comps` / `knot`, `pknot` comprehension (bounds `l <- [i+3 .. j-8]`, `k <- [l+4 .. j-4]`, `front`/`middle`/`back`, `getbeta`), `pk'` = `pkinit + e + 3*npp + fro + mid + bac + dangles`, `pkml` = +600, `npp` = 30, `pkinit` = 900, `pbp` = 0, `ul` = 40, `emptymid` | exact search space and pseudoknot energy terms |
| Reeder & Giegerich 2004 (BMC full text) | rules 1–3, csr-PK, relaxed csr-PK (146/212 = 68 %), complexity O(n⁴)/O(n²) | canonical class |
| ViennaRNA 2.7.2 (`RNA` module: `E_Hairpin`, `E_IntLoop`, `E_MLstem`, `E_ExtLoop`, `fold`, `params_load`) with the NNDB-precision parameter file of RNA-STRUCT-001 §12.2 | nested-loop energies | independent oracle |

No executable pknotsRG exists in the allowed environment (Haskell/ghc not installable; no PyPI/bioconda package; IPknot not on PyPI), so the oracle is
an independent implementation of the documented grammar (below), not pknotsRG output. pknotsRG uses Turner 1999 parameters, so its energies would differ from ours by the
parameter set regardless.

### Stage A findings (docs/claims versus code)

1. The previous doc claimed the unit "exact for the canonical single-knot class … heuristic-free within that class". The code only folded the span [i, r] of the knot: bases 5' of a and 3' of b' stayed unfolded **and unscored**,
   so the knot-only energy competed against the whole-sequence MFE (`GGGCGAAAGCCCAAGGGGAACCCCAACCCCAAGGGG`: knot-only −8.76 vs MFE −15.98 → "no knot"; with the hairpin folded the optimum is −21.70).
2. The knot energy omitted the three bridging nucleotides (3·npp), all helix-end terms and scored the loops as isolated MFEs (hairpin-closed, no exterior terms): designed H-type −8.76 (oracle −13.56).
3. The search space was narrower than pknotsRG's: loops v, w ≥ 1 (pknotsRG: v ≥ 0, w ≥ 2 incl. bridging), candidates with overlapping maximal helices discarded instead of resolved by rule 3.
4. Complexity statement "O(n³)" was wrong (four free boundaries × O(n³) loop MFE per candidate).
5. BWYV statement ("pseudoknot-free hairpin is the NN optimum") held only because of finding 2; with the full model the optimum contains a knot (−10.51 vs −7.90).

### Stage B (oracle) method

`pkoracle.py` (scratchpad, not in the repo): memoised naive grammar; V / WM / W re-implemented with ViennaRNA primitives (`E_Hairpin`, `E_IntLoop` incl. interior-loop mismatch arguments, `E_MLstem`, `E_ExtLoop`),
knots by explicit enumeration over (k, l), helix lengths by direct pair counting, helix stacks by `E_IntLoop(0,0,…)`, no span ordering, no precomputed tables. Parameter file: NNDB-precision (RNA-STRUCT-001 §12.2).
Sequences: random/designed 12–34 nt (3 seeds: 357 + 459 + 452), 809 designed nesting cases ≤ 45 nt, 2 multiloop-knot cases of 62–64 nt; both modes (single / recursive).

| Result | Count |
|--------|-------|
| compared (plain MFE of the sequence equals ViennaRNA) | 1 971 |
| energy differences vs oracle | 0 (|Δ| ≤ 0.005) — except 2, both with a 1-nt bulge inside a knot loop: C# applies the NNDB degeneracy / special-C bulge terms ViennaRNA lacks (verified: `CalculateFreeEnergy` of the loop −5.67 vs `eval_structure` −4.34, Δ = −1.33 = −0.90 special C − 0.43 degeneracy) |
| excluded (plain MFE itself contains an NNDB-only term) | 106 |
| sequences with a knot in the single-knot optimum | 386 of 809 designed; 113 + 140 + 133 of the random sets |
| recursive < single | 98 of 809 designed (95 comparable, all equal to the oracle) |

Exhaustive in-tree check (`PknotsRgBruteForce`, table-free enumeration of all flank structures × all canonical knots × all loop structures, scored with `EvaluateStructureEnergy`): 70 sequences of 11–16 nt, DP = enumeration to 0.006.
Old code on the same inputs: 28 of the first 120 designed sequences had `HasPseudoknot = false` where the optimum contains a knot (e.g. `UCAUGUUAGCUCAUGAAUGGCU` −4.22 → −5.93, empty loop v).

## Test Datasets (oracle values, ViennaRNA primitives + NNDB-precision file; locked in `PknotsRgOracleData`)

| Sequence | plain MFE | single-knot optimum | structure |
|----------|-----------|---------------------|-----------|
| `GGGGAACCCCAACCCCAAGGGG` (designed, a=GGGG u=AA b=CCCC v=AA a'=CCCC w=AA b'=GGGG) | −7.84 | −13.56 | `((((..[[[[..))))..]]]]` = 9.0 + 0.9 + 6×(−3.26) + 0.3 + 0.6 − 4.80 |
| `GGGCGAAAGCCCAAGGGGAACCCCAACCCCAAGGGG` (hairpin + knot) | −15.98 | −21.70 | `((((....))))..((((..[[[[..))))..]]]]` |
| `UCAUGUUAGCUCAUGAAUGGCU` (empty v) | −4.22 | −5.93 | `(((((..[[[[)))))..]]]]` |
| `UCAAGUAUUCCGAAGCUCAACGGGAAAAUGAGCUA` (pknotsRG README) | −6.54 | −12.31 | `.......(((((.[[[[[[.)))))...]]]]]].` |
| `GGCGCGGCACCGUCCGCGGAACAAACGG` (BWYV, PDB 437D sequence) | −7.90 | −10.51 | `..(((((..[[[[)))))......]]]]` (model optimum; crystal pairing not compared) |
| `GGGGAAAACCCC` (hairpin) | −5.28 | −5.28 | `((((....))))` (no knot) |

Plus 27 random/designed rows (hairpin-only optima included) in `PknotsRgOracleData.Rows`.

## Assumptions

1. **ASSUMPTION (adaptation): pknotsRG grammar and pseudoknot parameters over the shared Turner 2004 model.** pknotsRG's loop/stack parameters are Turner 1999; the lead requirement is one energy model. Deviations 1–7 are listed in
   Pseudoknot_Prediction.md §5.3 (loop component and helix-end terms in d2 form, no coaxial stacking, strict helices without the relaxed 1-nt bulge).
2. **ASSUMPTION: relaxed csr-PK (single 1-nt bulge in a helix) is not implemented** — implementable (pknotsRG `stack'`/`getbeta`), recorded as follow-up; strict class csr-PK is the paper's formal class.

---

## Recommendations for Test Coverage

1. MUST: designed H-type, flanked knot, empty-v knot → exact pairs / dot-bracket / ΔG (oracle).
2. MUST: oracle table rows (energy, structure, no-knot ⇒ plain MFE structure, validity, dot-bracket round trip).
3. MUST: DP = exhaustive enumeration (random + H-type-biased 11–16 nt).
4. MUST: invariants ΔG ≤ MFE, recursive ≤ single, HasPseudoknot ⇔ crossing pair.
5. SHOULD: null/empty/short/N/IUPAC/DNA/lower case/minLoopSize/temperature/concurrency/150 nt.

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
- **2026-11-25**: Stage A/B review: engine replaced by the exact pknotsRG DP over the shared Turner 2004 model (flanks, bridging nt, helix ends, v ≥ 0, rule 3); oracle + exhaustive verification; datasets, assumptions, corner cases rewritten; MCP tool `rna_pseudoknot_structure`.
