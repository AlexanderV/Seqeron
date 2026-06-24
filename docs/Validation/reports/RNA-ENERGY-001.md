# Validation Report: RNA-ENERGY-001 — RNA Folding Free Energy (Nearest-Neighbor / MFE)

- **Validated:** 2026-06-24   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.CalculateStemEnergy`, `CalculateHairpinLoopEnergy`,
  `CalculateMinimumFreeEnergy`, helper `GetTerminalMismatchEnergy` (+ `CalculateInternalLoopEnergy`,
  `CalculateBulgeLoopEnergy`, `CalculateMultibranchLoopEnergy`, `FillDp` Zuker DP).
  Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** ✅ CLEAN (no defect; genuine Turner-2004 NN model; parameters verified against
  the ViennaRNA `rna_turner2004.par` reference encoding and NNDB; all tests green)

This is an independent re-validation (fresh context). It reaches the same conclusion as the prior
(2026-06-12) report, but the cross-check this time was anchored on the **ViennaRNA Turner-2004
reference parameter file**, which the WebFetch tool *can* retrieve (the NNDB `urmc.rochester.edu`
HTML pages 404/403 to the tool). That file encodes the identical Turner-2004 parameter set, so it is
an authoritative independent witness for the stacking matrix, the hairpin-initiation array, and the
special tri/tetra/hexa-loop tables.

---

## Stage A — Description

### Sources opened & what they confirm
- **ViennaRNA `misc/rna_turner2004.par`** (fetched via `curl`, github raw) — the reference
  encoding of the Mathews/Turner 2004 set. Used as the primary numeric cross-check (values in dcal/mol):
  - `# stack` matrix: CG/CG −240, GC/GC −340, AU −110, UA −130 — consistent with the code's
    NNDB-notation WC values (CG/GC −2.36, GC/CG −3.42, AU/UA −1.10, UA/AU −1.33; the canonical
    Xia-1998 published values, used verbatim by the code).
  - `# hairpin` initiation array, index = loop size: 3→540, 4→560, 5→570, 6→540, 7→600, 8→550,
    9→640 … 30→770. **Exactly** the code's `HairpinLoopEnergies` table (5.4, 5.6, 5.7, **5.4**,
    6.0, 5.5, 6.4 … 7.7), including the non-monotonic dip at size 6.
  - `# Tetraloops`/`# Triloops`/`# Hexaloops`: CCUCGG 250, CUCCGG 270, CUACGG 280, CUGCGG 280,
    CCAAGG 330, CCCAGG 340, CCGAGG 350, CAACGG 550, …; CAACG 680, GUUAC 690; ACAGUGUU 180,
    ACAGUACU 280, ACAGUGCU 290, ACAGUGAU 360. **Exactly** the code's `SpecialHairpinLoops`
    (2.5/2.7/2.8/2.8/3.3/3.4/3.5/5.5/…, 6.8/6.9, 1.8/2.8/2.9/3.6).
  - `# Misc` line `410 360 50 370` → DuplexInit37=4.10, TerminalAU37=**0.50**. ViennaRNA rounds
    the per-AU-end term to 0.50; the genuine NNDB/Xia value is **+0.45** (see below), which is what
    the code uses.
- **Web search (NNDB Turner-2004 WC page + tm-example)** — independently confirms the per-AU-end
  penalty is **0.45 kcal/mol** and the intermolecular helix initiation 4.09 kcal/mol, i.e. the code's
  `TerminalAU_GU_Penalty = 0.45` is the authoritative NNDB value (ViennaRNA's 0.50 is a coarser
  rounding, not the published parameter). The prior report's GU-page quotes (per-GU-end = per-AU-end;
  GGUC/CUGG 3-stack total −4.12) were carried over and are consistent with the GU revision.
- **Xia et al. (1998) Biochemistry 37:14719** (search) — the WC stacking set
  GC/CG −3.42, CG/GC −2.36, GG/CC −3.26, AU/UA −1.10, UA/AU −1.33, AA/UU −0.93 is the foundational
  table assembled into Turner 2004; confirmed as the parameter set the code embeds.

### Formula check (additive NN model)
`ΔG°total = ΔG°init + Σ ΔG°stacking + per-end AU/GU penalty + Σ loop terms` — matches NNDB's
additive helix formula exactly. Stem energy = Σ nearest-neighbour stacks (with the GGUC/CUGG
3-stack special case) + 0.45 per AU/GU helix end. Hairpin = special-loop total OR
init(n) + terminal-mismatch (≥4 nt) + UU/GA & GG bonuses + special-GU closure + all-C penalty,
with Jacobson–Stockmayer extrapolation `ΔG(n)=ΔG(9)+1.75·R·T·ln(n/9)` for n>30. Sign convention
correct (more-negative = more stable; stacks negative, loop init positive).

### Edge-case semantics (sourced)
- Empty / null / too-short (< minLoop+2) sequence → MFE 0. ✅
- Single base pair → stem energy 0 (NN stacking needs ≥2 adjacent pairs). ✅
- Hairpin loop < 3 nt → prohibitive sentinel 100.0, returned *before* any sequence-dependent term
  (avoids the vacuous `"".All(c=='C')` all-C penalty on an empty loop). ✅
- Loop initiation is **non-monotonic** with size (6 = +5.4 < 4 = +5.6) — matches the ViennaRNA
  `# hairpin` array exactly; the spec correctly flags this. ✅

### Independent hand-computation (against ViennaRNA-verified parameters)
- **GGGAAACCC**: 2×GG/CC(−3.26) + init(3)=5.4 = **−1.12**.
- **GGGGAAAACCCC**: 3×GG/CC(−9.78) + init(4)=5.6 + tm(GAAC)=−1.1 = **−5.28**.
- **AU/UA stem stack**: −1.10 + 0.45 + 0.45 (both AU ends) = **−0.20**.

Stage A findings: none. Every checked stacking, hairpin-init, special-loop and penalty parameter
matches the ViennaRNA Turner-2004 reference and/or the canonical Xia-1998 / NNDB values. The model
is the genuine Turner-2004 nearest-neighbour model, not a simplified per-pair score. (A simplified
`PairEnergyByCode` table exists but feeds only the fast heuristic stem-loop scanner, not the energy
methods under test.)

## Stage B — Implementation

### Code path reviewed
`RnaSecondaryStructure.cs`: `StackingEnergies` (144–169, 16 WC + 20 GU), `HairpinLoopEnergies`
(175–183), `SpecialHairpinLoops` (189–216), `TerminalMismatchEnergies` (223–255, 96 entries),
`TerminalAU_GU_Penalty=0.45` (260), bonus/penalty consts (265–270), `SpecialGGUC_CUGG_3Stack=-4.12`
(276); methods `CalculateStemEnergy` (612), `CalculateHairpinLoopEnergy` (706),
`GetTerminalMismatchEnergy` (681), and `CalculateMinimumFreeEnergy` (1062) → `FillDp` (1113),
a genuine Zuker–Stiegler V/WM/W O(n³) DP (hairpin/stack/internal/bulge/multiloop), **not** a
placeholder. MFE and the `MfeStructure` traceback share the same `FillDp`, so the reconstructed
structure's energy equals the scalar MFE by construction.

### Formula realised correctly?
- `CalculateStemEnergy`: sums NN stacks via the `"XY/X'Y'"` key, special-cases the 4-pair
  G-C·G-U·U-G·C-G run to −4.12 (skipping 2 stacks), adds +0.45 at each AU/GU-terminating end. ✅
- `CalculateHairpinLoopEnergy`: returns the experimental total for special tri/tetra/hexaloops,
  else init(n) + (≥4 nt) terminal mismatch + UU/GA(−0.9) & GG(−0.8) bonuses + special-GU closure
  (−2.2) + all-C penalty (3 nt: +1.5; >3 nt: 0.3·n+1.6); JS extrapolation for n>30. ✅
- `CalculateMinimumFreeEnergy`: real Zuker DP over the same verified tables; rounds to 2 dp. ✅

### Cross-verification table recomputed vs code (via tests, exact values)
| Case | Expected (sourced) | Test | Result |
|------|--------------------|------|--------|
| GGGAAACCC MFE | −1.12 | `…_SimpleHairpin_MatchesTurnerManualCalc` | ✅ |
| GGGGAAAACCCC MFE | −5.28 | `…_FourPairGC_MatchesTurner` | ✅ |
| WC stacks (GC/CG,GG/CC,CG/GC,AU/UA,UA/AU,AA/UU) + end penalties | exact | `…_WatsonCrickStacking_MatchesNNDB` (×8) | ✅ |
| GGUC/CUGG 3-stack | −4.12 | `…_GGUC_CUGG_3Stack_MatchesNNDB` | ✅ |
| Terminal AU penalty +0.45 | exact | `…_TerminalAUPenalty_MatchesNNDB` | ✅ |
| Hairpin init 3/4/5/6/9 nt + TM | exact | `…_Initiation_MatchesNNDB` (×5) | ✅ |
| Inner-AU NNDB hairpin example | −3.42 | `…_InnerAUPenalty_NNDBExample_MatchesReference` | ✅ |
| Empty/null/single-pair/no-structure → 0 | 0 | SE-002/004, MFE-002/003 | ✅ |

### Variant/delegate consistency
`PredictStructure`/`PredictMfeStructure` and the pseudoknot energy routines all route stacking
through `CalculateStemEnergy` and loop folding through `CalculateMinimumFreeEnergy`/`CalculateHairpinLoopEnergy`,
so the same verified tables drive every consumer. MFE traceback shares `FillDp`. ✅

### Test quality audit
Assertions check **exact** sourced numbers (`Is.EqualTo(...).Within(0.01)`) with the NNDB derivation
in comments — not tautologies. 128 tests in `RnaSecondaryStructureTests`; 280 across all
`RnaSecondaryStructure*` files. All cover the spec IDs, invariants, and edge cases.

### Findings / defects
None affecting any value under test. (Cosmetic: the `GetInternalLoop2x3MismatchEnergy`
purine/pyrimidine ternaries collapse to equal constants in some RA/YA branches — internal-loop
2×3 detail, outside this unit's scope; noted only.)

## Verdict & follow-ups
- **Stage A: PASS** — genuine Turner-2004 NN model; every checked stacking / hairpin-init /
  special-loop / penalty parameter matches the ViennaRNA `rna_turner2004.par` reference and the
  canonical Xia-1998 / NNDB values; sign convention and non-monotonic loop semantics correct; the
  per-AU/GU-end penalty +0.45 is the published NNDB value (not ViennaRNA's rounded 0.50).
- **Stage B: PASS** — code realises the additive NN model and a real Zuker MFE DP over the verified
  tables; both worked hairpins reproduce exactly (−1.12, −5.28).
- **End state: ✅ CLEAN.** No code changed. Build clean (0 warnings); `RnaSecondaryStructure*`
  tests 280/280 pass (128 in the canonical test class). (Full suite unchanged from prior 4486/0
  baseline — no code touched.)

### Sources
- ViennaRNA `misc/rna_turner2004.par` (raw.githubusercontent.com/ViennaRNA/ViennaRNA/master) —
  reference encoding of the Turner-2004 set: stacking matrix, hairpin-init array, tri/tetra/hexaloop tables.
- NNDB Turner 2004 (rna.urmc.rochester.edu/NNDB/turner04/): wc, gu, hairpin, tm-example — per-AU-end
  +0.45, init 4.09 (HTML pages 404 to the fetch tool; values confirmed via web search of those pages).
- Xia et al. (1998) Biochemistry 37:14719–14735 (WC NN parameters).
- Chen, Turner et al. (2012) Biochemistry doi:10.1021/bi3002709 (revised GU parameters; carried from prior report).
