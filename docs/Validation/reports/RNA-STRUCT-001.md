# Validation Report: RNA-STRUCT-001 — RNA Secondary Structure Prediction (MFE DP + pseudoknot predictors)

- **Validated:** 2026-06-24   **Area:** RnaStructure
- **Canonical method(s):**
  - `RnaSecondaryStructure.CalculateMfeStructure` / `PredictStructureMfe` — MFE-optimal pseudoknot-free
    structure via Zuker–Stiegler (1981) DP traceback (Turner 2004 nearest-neighbour energy model).
  - `RnaSecondaryStructure.CalculateMinimumFreeEnergy` — scalar MFE (same `FillDp`).
  - `RnaSecondaryStructure.PredictStructurePseudoknot` — canonical H-type / csr-PK (Reeder & Giegerich 2004).
  - `RnaSecondaryStructure.PredictStructurePseudoknotRecursive` — recursive pknotsRG grammar (nested/multiple knots).
  - Energy helpers `CalculateStemEnergy` / `CalculateHairpinLoopEnergy`; `PairType` / `CanPair` / `GetBasePairType`;
    `GenerateFullDotBracket` / `GeneratePseudoknotDotBracket` / `DetectPseudoknots`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** CLEAN

---

## Stage A — Description

### Sources opened & what they confirm
- **Zuker & Stiegler (1981), NAR 9(1):133–148** (PMC326673) — the MFE DP: structures decompose into
  hairpin / stacking / bulge–internal / multibranch loops, total ΔG = sum of loop contributions;
  polynomial DP. Matches the V/W/WM decomposition used.
- **Ward, Datta, Wise, Mathews (2017), NAR 45(14):8541** (PMC5737859) — C(i,j)=min(hairpin, interior/bulge
  over inner pair, multiloop); affine multiloop model `a + b·branches + c·unpaired` is O(n³)/O(n²) and "the
  simplest model is best". Confirms the implemented recurrence family and complexity.
- **NNDB Turner 2004** (Mathews et al. 2004, PNAS 101:7287) — stacking, hairpin-initiation, terminal-mismatch,
  AU-end-penalty tables; worked Hairpin Examples 1 (−1.4) and 2 (−1.9) and WC-helix example (−10.13). The
  live `rna.urmc.rochester.edu` server returns HTTP 404 to the fetch tool and web.archive.org is blocked from
  this environment (both noted in `docs/Evidence/RNA-MFE-001-Evidence.md`, which captured the values from a
  Wayback snapshot); the values were re-confirmed here by reproducing those worked examples (below) and by
  hand computation.
- **Reeder & Giegerich (2004), BMC Bioinformatics 5:104** (PMC514697) and **Reeder, Steffen, Giegerich (2007),
  NAR 35:W320** (PMC1933184) — csr-PK: two crossing helices + three intervening loops u/v/w; canonization
  rules (equal-length bulge-free helices, maximal extent); pseudoknot-specific penalties (initiation 9.0
  kcal/mol; 0.3/unpaired loop nt; 0.0/bp inside knot); recursive class allows loops to fold internally
  including further knots; per-interval competition with unknotted foldings; excluded classes = kissing
  hairpins, triple-crossing/complex. All match the implementation and Evidence docs.
- **Antczak et al. (2018)** — crossing condition i<k<j<l (basis of `DetectPseudoknots` / second bracket layer).

### Formula / definitions check
- Allowed pairs A-U/U-A/G-C/C-G (WC) and G-U/U-G (wobble) — `PairType`. Correct.
- Minimum hairpin loop = 3 nt (pair needs j−i ≥ 4); minLoopSize < 3 clamped to 3. Correct (NNDB steric rule).
- Empty/null/too-short → open chain, ΔG = 0; unfoldable → 0 (open chain always available, optimum never > 0).
- Pseudoknot accepted only when ΔG strictly below the pseudoknot-free MFE (the 9.0 initiation penalty
  rationale) — prevents spurious knots. Correct.

### Independent cross-check (numbers)
- **Hairpin GGGAAACCC** by hand: stacks GG/CC = −3.26 (×2) closing G-C, hairpin init(3) = +5.4, no AU end →
  −3.26 −3.26 + 5.4 = **−1.12 kcal/mol**, structure `(((...)))`. (MS-002.)
- **NNDB Hairpin Example 1** `CACAAAAAAAUGUG` → −1.41; **Example 2** `CACAGAAAGUGUG` → −1.91 — reproduced by
  the code in the sibling RNA-MFE-001 validation.
- **Designed H-type** `GGGGAACCCCAACCCCAAGGGG`: two crossing 4-bp G·C helices, `((((..[[[[..))))..]]]]`,
  ΔG = −8.76 < pseudoknot-free MFE −6.94. Knot recovered with lower ΔG (matches Evidence).

### Findings / divergences
None. Description (MFE DP, Turner 2004 NN model, csr-PK predictors) is biologically and mathematically correct.

---

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`
- **Turner stacking table** `StackingEnergies` (l.144) — 16 WC + 20 GU entries; spot-checked against NNDB:
  GC/CG −3.42, GG/CC −3.26, CG/GC −2.36, AU/UA −1.10, UA/AU −1.33, AA/UU −0.93; AU/GU end penalty +0.45
  (l.260). Genuine published values.
- **Hairpin initiation** `HairpinLoopEnergies` (l.175) sizes 3→5.4 … 9→6.4, n>9 extrapolated by
  ΔG(9)+1.75·R·T·ln(n/9). Matches loop.txt.
- **MFE DP fill** `FillDp` (l.1113) — V(i,j) = min(hairpin, stacking, special GGUC 3-stack, internal/bulge,
  multiloop); WM with dangle options; W(j) exterior recurrence; affine multiloop a=9.25, c=−0.63, unpaired=0
  (documented simplification). This is the Zuker–Stiegler recurrence, not a Nussinov approximation.
- **MFE traceback** `TracebackMfe` (l.1536) — a true DP traceback: at each W/V/WM/M sub-problem it
  re-derives the SAME recurrence candidate values and follows the option that attains the stored optimum
  (`|target − candidate| ≤ 1e-6`), recording a pair on each V entry. NOT greedy. Energy returned = `w[n-1]`,
  identical to the scalar MFE — so MS-001 (`FreeEnergy == CalculateMinimumFreeEnergy`) is the locking invariant.
- **H-type predictor** `PredictStructurePseudoknot` (l.2603) — enumerates two crossing helices with maximal
  extension (`MaxHelixLength`, rules 1–2), `EvaluateHType` scores Turner stacking(a)+stacking(b)+9.0 init+
  0.3·(unpaired loop nts)+MFE(u,v,w loop spans); accepts only when ΔG < plain MFE − eps. Matches pknotsRG.
- **Recursive predictor** `PredictStructurePseudoknotRecursive` (l.2851) — recursive folder over loops u/v/w
  (a loop may contain a further knot), top-level chaining of multiple knots, and an enclosing helix that
  over-arches a knot; same energy model/penalties. Matches the recursive csr-PK class (PARTIAL, documented).

### Cross-verification table recomputed vs code (scratch run, then removed)
| Sequence | MFE dot-bracket | MFE ΔG | scalar match | PK | PK ΔG | knot | REC | REC ΔG | knot |
|----------|-----------------|--------|--------------|----|----|------|-----|----|------|
| GGGAAACCC | `(((...)))` | −1.12 | ✓ | same | −1.12 | no | same | −1.12 | no |
| GGGGAAAACCCC | `((((....))))` | −5.28 | ✓ | same | −5.28 | no | same | −5.28 | no |
| GGGGAACCCCAACCCCAAGGGG | `(((....)))..(((....)))` | −6.94 | ✓ | `((((..[[[[..))))..]]]]` | **−8.76** | **yes** | same | −8.76 | yes |
| AAAAAAAAGGGG…GGGGUUUUUUUU (over-arch, 38nt) | `…((((((((((((……))))))))))))` | −13.05 | ✓ | (=MFE) | −13.05 | no | `((((((((((((..[[[[..))))..]]]]))))))))` | **−14.37** | **yes** |
| GGCGCGGCACCGUCCGCGGAACAAACGG (BWYV) | `..(((((......)))))..........` | −8.20 | ✓ | (=MFE) | −8.20 | no | (=MFE) | −8.20 | no |

Key results, all reproduced by the code:
- Reconstructed MFE energy equals the scalar MFE for every input (traceback is optimal).
- Designed H-type knot recovered with ΔG (−8.76) strictly below the pseudoknot-free MFE (−6.94).
- The over-arching nested knot is recovered ONLY by the recursive method (−14.37 < single/MFE −13.05).
- BWYV (tertiary-stabilised, PDB 437D) is NOT recovered as the MFE — the documented energy-model floor.

### Variant/delegate consistency
`CalculateMfeStructure`, `PredictStructureMfe`, and `CalculateMinimumFreeEnergy` share `FillDp` and agree on
energy. `PredictStructurePseudoknotRecursive` reproduces `PredictStructurePseudoknot` on the single canonical
H-type and never returns a structure worse than the plain MFE (fallback baseline). `PairType`/`CanPair`/
`GetBasePairType` agree on the pair set.

### Test quality audit
`RnaSecondaryStructure_MfeStructure_Tests` (MS-001…MS-013), `…_PredictStructurePseudoknot_Tests`, and
`…_PredictStructurePseudoknotRecursive_Tests` ran green (49 in the filtered subset). Assertions check exact
ΔG (−1.12, −5.28, −8.76, −14.37), exact dot-brackets, the energy-consistency invariant, WC/wobble-only pairs,
each position paired ≤ once, ΔG ≤ plain MFE, genuine crossings via `DetectPseudoknots`, no-spurious-knot on
plain/random sequences, and the BWYV non-recovery (documents the limit). Real, exact-value, deterministic.

### Findings / defects
None. The MFE DP traceback is optimal (not greedy), the Turner 2004 NN parameters are the genuine published
values (worked examples + hand calc reproduced exactly), and both pseudoknot predictors recover designed
H-type / nested knots with ΔG strictly below the pseudoknot-free MFE.

---

## Verdict & follow-ups
- **Stage A: PASS** — Zuker–Stiegler MFE DP, Turner 2004 NN model, and the Reeder & Giegerich csr-PK /
  recursive-pknotsRG descriptions are correct against the primary literature.
- **Stage B: PASS** — MFE traceback returns the thermodynamically optimal pseudoknot-free structure (energy ==
  scalar MFE for all tested inputs); H-type and recursive predictors recover designed knots with lower ΔG than
  the MFE; BWYV correctly not recovered.
- **State: CLEAN** — the residual has two distinct, both-documented parts (recorded in
  `docs/Validation/LIMITATIONS.md` §1), sharpened during a focused re-evaluation of part (a):
  - **(a) Non-csr-PK algorithm-class boundary** — kissing hairpins (loop–loop / "H-H" pseudoknots),
    triple-crossing / complex helix interactions, non-canonical bulged or unequal-length helices.
    Reeder & Giegerich (2004), [PMC514697](https://pmc.ncbi.nlm.nih.gov/articles/PMC514697/), state
    *verbatim*: "*More complex knotted structures like triple crossing helices or kissing hairpins …
    are excluded from sr-PK*" — so these classes are genuinely outside the two-crossing-helix csr-PK
    grammar realised by `PredictStructurePseudoknot(Recursive)` (a real extension, not a duplicate of
    the existing predictors). **A faithful kissing-hairpin detector is *not* bundled because no
    sourceable energy model for the loop–loop interaction exists.** The principal algorithm for
    *intramolecular* kissing hairpins — Sperschneider, Datta & Wise (2011), RNA 17:27–38,
    [PMC3004063](https://pmc.ncbi.nlm.nih.gov/articles/PMC3004063/) — states *verbatim*: "*No
    experimentally measured energy parameters for intramolecular kissing hairpins have been established
    to date and thus heuristic energy estimation has to be used*" and "*heuristic methods that allow
    arbitrary pseudoknots strongly depend on the quality of energy parameters, which are not yet
    available for complex pseudoknots*". Its loop–loop model
    `ΔG = Σ stack(S₁,S₂,S₃) + α + β(l₁+l₂+l₄+l₅) + γ·l₃` with **α = 9.0, β = 0.5, γ = 0.0 kcal/mol**
    relies on an *estimated* per-loop-nucleotide penalty (β), which is NOT a Turner-2004
    nearest-neighbour parameter. Adding it would require introducing an unsourced (heuristic,
    non-measured) thermodynamic constant — forbidden under the no-fabrication / cite-verbatim policy —
    so the class is left as an honest documented residual; use the DotKnot / pknotsRG-kissing reference
    tools for it. This boundary is potentially extendable only if/when measured loop–loop parameters
    are published.
  - **(b) Energy-model floor (irreducible)** — tertiary-stabilised knots as the MFE structure (e.g.
    BWYV / PDB 437D) are not recoverable by *any* nearest-neighbour thermodynamic model. This is a
    physics-of-the-model floor, not an algorithm gap, and is not recoverable in this unit.
- No code changed (the re-evaluation confirmed no faithful, citable, deterministically-testable
  kissing-hairpin energy model is available to bundle). No new defects logged. Registry Status unchanged.
