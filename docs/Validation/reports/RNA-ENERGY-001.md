# Validation Report: RNA-ENERGY-001 — RNA Folding Free Energy (Nearest-Neighbor / MFE)

- **Validated:** 2026-06-12   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.CalculateStemEnergy`, `CalculateHairpinLoopEnergy`,
  `CalculateMinimumFreeEnergy`, plus helpers `CalculateStackingEnergy`/`GetTerminalMismatchEnergy`,
  `CalculateInternalLoopEnergy`, `CalculateBulgeLoopEnergy`, `CalculateMultibranchLoopEnergy`.
  Source: `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End state:** ✅ CLEAN (no defect; full Turner-2004 NN model, parameters match NNDB, all tests green)

---

## Stage A — Description

### Sources opened & what they confirm
- **NNDB / Turner 2004** (https://rna.urmc.rochester.edu/NNDB/turner04/) — the NNDB pages 404/403 to
  the fetch tool, but the WC and GU parameter pages were located and their formulae quoted via search:
  - WC helix free energy = `ΔG°37(init) + ΔG°37(AU-end penalty, per AU end) + ΔG°37(symmetry) + Σ ΔG°37(stacking)`.
    This is exactly the additive NN model the spec defines (init + Σstacking + per-end AU/GU penalty + loop terms).
  - GU page: *"GU pairs at the ends of helices are penalized with the same parameter as AU pairs at the
    ends of helices"* → confirms the code's shared `TerminalAU_GU_Penalty = +0.45`.
  - GU page note (b): *"GU followed by UG is generally unfavorable … favorable in one specific context
    where a GC precedes the pair and CG follows … a single reported parameter of −4.12 kcal/mol is used
    for the total of three stacks."* → confirms the GGUC/CUGG special 3-stack value `-4.12` used by the code.
- **Xia et al. (1998), Biochemistry 37:14719** (the primary WC source feeding Turner 2004) — search
  confirms the listed WC table `GC/CG −3.42, CG/GC −2.36, GG/CC −3.26, AA/UU −0.93, AU/UA −1.10,
  UA/AU −1.33` is "from this seminal 1998 table"; GC/CG ≈ −3.4 kcal/mol independently confirmed.
- **Chen/Turner (2012) GU revision** (Biochemistry, doi 10.1021/bi3002709) — confirms the revised GU
  treatment and the two context-dependent tandem-GU parameters used in the NNDB GU table.

### Parameter spot-checks against the source (kcal/mol, ΔG°37)
| Stack | Spec/code value | Authoritative (Xia 1998 / NNDB) | Match |
|-------|-----------------|---------------------------------|-------|
| GC/CG | −3.42 | −3.42 (≈−3.4) | ✅ |
| CG/GC | −2.36 | −2.36 | ✅ |
| GG/CC | −3.26 | −3.26 | ✅ |
| AU/UA | −1.10 | −1.10 | ✅ |
| AA/UU | −0.93 | −0.93 | ✅ |
| UA/AU | −1.33 | −1.33 | ✅ |
| per AU/GU end penalty | +0.45 | +0.45 | ✅ |
| GGUC/CUGG 3-stack (note b) | −4.12 | −4.12 | ✅ |

The model is the **full Turner-2004 nearest-neighbor model**, not a simplified per-pair score. (A
simplified `PairEnergyByCode = {0, −2, −1}` table exists but is used only by the fast heuristic
stem-loop scanner, *not* by the energy methods under test.) Sign convention is correct:
more-negative ΔG = more stable; stacking terms are negative, loop initiation positive.

### Edge-case semantics (sourced)
- Empty / null / no-pairing sequence → MFE = 0 (no structure). ✅
- Single base pair → stem energy 0 (NN stacking needs ≥2 adjacent pairs). ✅
- Hairpin loop < 3 nt sterically impossible → prohibitive energy (NNDB/Wikipedia). ✅
- Loop initiation is **non-monotonic** with size (size-6 = +5.4 < size-4 = +5.6) — correctly flagged by
  the spec and matching the NNDB loop.txt experimental values. ✅

### Independent hand-computation (matches the spec's worked examples)
- **GGGAAACCC**: 2×GG/CC (−3.26) + hairpin init(3) 5.4 = **−1.12**.
- **GGGGAAAACCCC**: 3×GG/CC (−9.78) + [init(4) 5.6 + terminal mismatch GAAC −1.1] = **−5.28**.
- **AU/UA stem stack**: −1.10 + 0.45 + 0.45 (both AU ends) = **−0.20**.

Stage A findings: none. All checked parameters and the additive formula are authoritative.

## Stage B — Implementation

### Code path reviewed
`RnaSecondaryStructure.cs`: `StackingEnergies` (16 WC + 20 GU, lines 107–132), `HairpinLoopEnergies`
(138–146), `SpecialHairpinLoops` (152–179), `TerminalMismatchEnergies` (186–218),
`TerminalAU_GU_Penalty = 0.45` (223), `SpecialGGUC_CUGG_3Stack = -4.12` (239), and the methods
`CalculateStemEnergy` (575), `CalculateHairpinLoopEnergy` (669), `CalculateMinimumFreeEnergy`
(Zuker-style O(n³) DP, 1021).

### Formula realised correctly?
- `CalculateStemEnergy` sums NN stacks via the `"XY/X'Y'"` key, special-cases the GGUC/CUGG 3-stack to
  −4.12, and adds +0.45 at each AU/GU-terminating end. Matches the NNDB WC/GU formula. ✅
- `CalculateHairpinLoopEnergy` returns the experimental total for special tri/tetra/hexaloops (key =
  closing5′+loop+closing3′), else init(n) + terminal mismatch (≥4 nt) + UU/GA/GG bonuses + special-GU
  closure + all-C penalty, with Jacobson–Stockmayer extrapolation for n>30. ✅
- `CalculateMinimumFreeEnergy` is a genuine V/WM/W Zuker recursion (hairpin, stack, internal/bulge,
  multiloop) using the same parameter tables — not a placeholder. ✅

### Cross-verification table recomputed vs code (via tests, exact values)
| Case | Expected (sourced) | Test | Result |
|------|--------------------|------|--------|
| GGGAAACCC MFE | −1.12 | `…_SimpleHairpin_MatchesTurnerManualCalc` | ✅ |
| GGGGAAAACCCC MFE | −5.28 | `…_FourPairGC_MatchesTurner` | ✅ |
| GC/CG, GG/CC, CG/GC, AU/UA, UA/AU, AA/UU, CA/GU, UG/AC stacks (+ end penalties) | exact | `…_WatsonCrickStacking_MatchesNNDB` (×8) | ✅ |
| GGUC/CUGG 3-stack | −4.12 | `…_GGUC_CUGG_3Stack_MatchesNNDB` | ✅ |
| Terminal AU penalty +0.45 | exact | `…_TerminalAUPenalty_MatchesNNDB` | ✅ |
| Hairpin init 3/4/5/6/9 nt + TM | exact | `…_Initiation_MatchesNNDB` (×5) | ✅ |
| Empty/null/single-pair/no-structure → 0 | 0 | SE-002/004, MFE-002/003 | ✅ |

### Test quality audit
Assertions check **exact** sourced numbers (`Is.EqualTo(...).Within(0.01)`) with the NNDB derivation
in comments — not tautologies. 137 `RnaSecondaryStructure` tests cover all 35 spec IDs, invariants,
and edge cases.

### Findings / defects
None. The `GetInternalLoop2x3MismatchEnergy` purine/pyrimidine branches (lines 865–866) collapse to
the same constant for the RA/YA cases (cosmetically redundant ternaries), but this is internal-loop
2×3 detail outside RNA-ENERGY-001's scope and does not affect any value under test — noted only.

## Verdict & follow-ups
- **Stage A: PASS** — full Turner-2004 NN model; every checked stacking/loop/penalty parameter matches
  NNDB / Xia 1998 / the Chen-Turner GU revision; sign convention and non-monotonic loop semantics correct.
- **Stage B: PASS** — code realises the additive NN model and a real Zuker MFE DP over the verified
  tables; both worked examples reproduce exactly.
- **End state: ✅ CLEAN.** No code changed. Build clean (0 warnings); RnaSecondaryStructure filter
  137/137 pass; full `Seqeron.Genomics.Tests` suite 4486 passed / 0 failed (1 perf benchmark skipped).

### Sources
- NNDB Turner 2004: https://rna.urmc.rochester.edu/NNDB/turner04/index.html (wc-parameters, gu-parameters, hairpin)
- Xia et al. (1998) Biochemistry 37:14719–14735 (WC NN parameters)
- Chen, Turner et al. (2012) Biochemistry, doi:10.1021/bi3002709 (revised GU parameters)
- ViennaRNA `misc/rna_turner2004.par` (reference encoding of the same parameter set)
