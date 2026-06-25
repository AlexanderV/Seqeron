# Validation Report: RNA-ACCESS-001 — McCaskill Unpaired (Accessibility) Probabilities

- **Validated:** 2026-06-25   **Area:** RnaStructure
- **Canonical method(s):** `RnaSecondaryStructure.CalculateUnpairedProbabilities`, `RnaSecondaryStructure.CalculateRegionUnpairedProbability`
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN

## Canonical method(s)
`CalculateUnpairedProbabilities` (Z, base-pair probabilities, per-base unpaired probabilities,
ensemble free energy) and `CalculateRegionUnpairedProbability` (joint region-unpaired
probability = RNAplfold accessibility).

- **Source:** `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RnaSecondaryStructure.cs`
  (`CalculateUnpairedProbabilities` L2556–2643; `CalculateRegionUnpairedProbability` L2672–2719;
  shared Boltzmann inside DP `FillPartitionDp` L2733–2963)
- **Tests:** `tests/Seqeron/Seqeron.Genomics.Tests/RnaSecondaryStructure_UnpairedProbabilities_Tests.cs`

## Stage A — Description

### Sources opened this session & what they confirm
- **McCaskill JS (1990), *Biopolymers* 29(6-7):1105–1119 (PMID 1695107).** Defines the equilibrium
  partition function Z = Σ_S exp(−E(S)/RT) over all pseudoknot-free secondary structures, the
  base-pair probability P(i,j) = (Σ_{S∋(i,j)} exp(−E(S)/RT))/Z via the inside/outside recursions,
  and the per-residue equilibrium quantities. RT is the thermal energy; p(S) = exp(−E(S)/RT)/Z.
- **Bernhart, Hofacker, Stadler (2006), *Bioinformatics* 22(5):614–615, "Local RNA base pairing
  probabilities in large sequences" + the RNAplfold man page.** Defines the *accessibility* / local
  unpaired probability: the `_lunp` output reports, for each position i and length u, the probability
  that the stretch of u consecutive nucleotides **ending at i** is **entirely unpaired**. This is the
  joint "whole window unpaired" probability — exactly the McCaskill ensemble quantity
  Z_open(window)/Z, where Z_open restricts the ensemble to structures with no pair incident to any
  window position.
- **Turner-2004 nearest-neighbour parameters (NNDB, rna.urmc.rochester.edu/NNDB/turner04).** The
  loop-energy model E(S) (hairpin initiation table, terminal mismatch, stacking, internal/bulge,
  multiloop offset/helix terms, dangling ends, terminal-AU/GU penalty 0.45) is the SAME model the
  MFE folder uses, so the partition function is its Boltzmann mirror (sum/multiply where MFE does
  min/add).

### Formula check
- Z = Σ_S exp(−E(S)/RT), p_unpaired(i) = 1 − Σ_j P(i,j), ΔG_ensemble = −RT·ln Z — all match
  McCaskill (1990) / Lorenz et al. (2011, ViennaRNA). RT = R·T/1000 with R = 1.987 cal/(mol·K),
  T = 310.15 K (37 °C) → RT = 0.61626805 kcal/mol. β = 1/RT. ✔
- Region accessibility = Z_open(window)/Z matches Bernhart 2006 / RNAplfold `_lunp` semantics
  (probability a length-L window is wholly unpaired). ✔
- MinHairpinLoop = 3 (a pair (i,j) admissible only when j−i > 3): standard nearest-neighbour rule. ✔

### Edge-case semantics (sourced / mathematically defined)
- Empty / null sequence → only the empty structure, Z = 1, ΔG = 0, no pairs. ✔
- Single base / span < minLoopSize+2 → no admissible pair, Z = 1, every base unpaired (p=1). ✔
- Non-ACGU characters → cannot pair (not in {A·U, G·C, G·U}); behave as permanently unpaired. ✔
- Non-positive temperature → physically meaningless (RT≤0 ⇒ exp blowup); must throw. ✔
- Region length ≤ 0 or window outside the sequence → invalid; must throw. ✔
- Monotonicity (Bernhart): a longer queried window is never *more* accessible (forbidding more
  positions can only shrink Z_open) — a genuine property of Z_open/Z. ✔

### Independent cross-check (numbers derived THIS session, NOT read from the code)
Two analytic two-state pins, derived by hand from the published Turner-2004 tables, plus an
independent Python brute-force ensemble enumerator (`scratchpad/brute.py`). ViennaRNA could not be
installed (no wheel for the available Python; source build fails) and would not be bit-identical to
this custom Turner model anyway, so analytic derivation is used as the oracle.

**GAAAC** — only structures: open chain (E=0, w=1) and the single G(0)·C(4) hairpin over loop "AAA"
(size 3, ΔG = initiation(3) = 5.4; G-C closing ⇒ no terminal-AU; full-span helix ⇒ no dangles):
- w = exp(−5.4/RT) = 0.00015650527649215345
- Z = 1 + w = **1.0001565052764922**
- P(0,4) = w/Z = **0.00015648078642340854**
- p_unpaired(0)=p_unpaired(4) = 1−P = **0.9998435192135765**;  p_unpaired(2) = 1
- ΔG_ensemble = −RT·ln Z = **−9.64416549414892e-05**
- region(whole 5-nt window) = Z_open/Z = 1/Z = **0.9998435192135765**

**CAAAAG** — open chain + single C(0)·G(5) hairpin over loop "AAAA" (size 4, ΔG = initiation(4) +
terminal-mismatch["CAAG"] = 5.6 + (−1.5) = 4.1; full-span helix ⇒ no dangles, no AU end):
- Z = 1 + exp(−4.1/RT) = **1.0012902114608**
- P(0,5) = **0.0012885489601637966**
- p_unpaired(0) = **0.9987114510398362**
- ΔG_ensemble = **−0.0007946036078507769**;  region(whole 6-nt window) = 1/Z = 0.9987114510398362

The independent brute-force enumerator reproduces the GAAAC ensemble exactly (Z = 1.0001565052764922,
2 structures). For longer sequences (GGGAAACCC) a hand oracle cannot match without re-implementing
the full internal-loop/dangle model (which would merely echo the code), so the two full-span
two-state pins (no dangle/multiloop ambiguity) are the decisive independent checks; longer cases are
covered by exact mathematical invariants below.

### Findings / divergences
None. Description is mathematically and biologically correct.

## Stage B — Implementation

### Code path reviewed
- `CalculateUnpairedProbabilities` (L2556): validates T>0, clamps minLoopSize≥3, upper-cases & maps
  T→U, handles empty/too-short, builds the Boltzmann inside DP via `FillPartitionDp`, then:
  - **p_unpaired(i)** is computed EXACTLY as Z_forbid(i)/Z — re-running `FillPartitionDp` with
    position i forbidden from pairing (L2603–2613). This is the constrained-partition definition, not
    an outside-recursion approximation.
  - **P(i,j)** = Z_require(i,j)/Z — re-running the DP constrained so (i,j) must pair (L2621–2637).
  - ΔG_ensemble = −RT·ln Z (L2640).
- `CalculateRegionUnpairedProbability` (L2672): validates T>0 and windowLength>0, computes the window
  start, throws if the window does not fit, returns Z_open(window)/Z where Z_open forbids every base
  in [windowStart..windowEnd] (L2701–2718).
- `FillPartitionDp` (L2733): Boltzmann mirror of the MFE `FillDp` — hairpin, stacking, internal/bulge,
  multiloop (offset+helix+dangle Boltzmann weights), external loop; `Forbidden(p)` excludes pairs
  incident to the forbidden window; `PairOk`/`UnpairedOk`/`RequiredUnpairedIn` enforce the
  must-pair constraint. Uses the SAME Turner-2004 energy functions as the MFE folder.

### Formula realised correctly?
Yes. Z = Σ exp(−E/RT) over the same structure space as the MFE DP; p_unpaired = Z_forbid/Z;
P(i,j) = Z_require/Z; region = Z_open/Z; ΔG = −RT·ln Z. Numerically guarded: non-finite/zero Z
returns an all-unpaired / empty-probability result instead of leaking NaN (L2596, L2703).

### Cross-verification table recomputed vs code (executed against the built net10.0 assembly)

| Quantity | Analytic (hand-derived) | Code output | Match |
|---|---|---|---|
| GAAAC Z | 1.0001565052764922 | 1.0001565052764922 | ✔ exact |
| GAAAC P(0,4) | 0.00015648078642340854 | 0.0001564808 | ✔ |
| GAAAC p_unpaired(0) | 0.9998435192135765 | 0.9998435192135765 | ✔ |
| GAAAC ΔG_ensemble | −9.64416549414892e-05 | −9.64416549414892e-05 | ✔ exact |
| GAAAC region(5@4) | 1/Z = 0.9998435192135765 | 0.9998435192135765 | ✔ exact |
| CAAAAG Z | 1.0012902114608 | 1.0012902114608 | ✔ exact |
| CAAAAG P(0,5) | 0.0012885489601637966 | 0.0012885489601637966 | ✔ exact |
| CAAAAG p_unpaired(0) | 0.9987114510398362 | 0.9987114510398362 | ✔ exact |
| CAAAAG ΔG_ensemble | −0.0007946036078507769 | −0.0007946036078507769 | ✔ exact |
| CAAAAG region(6@5) | 1/Z = 0.9987114510398362 | 0.9987114510398362 | ✔ exact |
| single "G" Z / pu | 1.0 / [1.0] | 1.0 / [1.0] | ✔ |
| GGGNNNCCC | finite, deterministic, N's pu=1 | Z=7.5426…, N's pu=1 | ✔ |

Mathematical invariants verified against the code on multiple sequences:
- McCaskill consistency p_unpaired(i) = 1 − Σ_j P(i,j) (test MCC-003) ✔
- ΔG_ensemble = −RT·ln Z ≤ MFE for the same Turner model (test MCC-004) ✔
- region monotone non-increasing in window length over L=1..14 (test MCC-009) ✔
- region(L=1 @ i) = per-base p_unpaired(i) — cross-consistency of the two public methods (MCC-010) ✔

### Variant/delegate consistency
The two public methods share `FillPartitionDp` and the same RT/energy model; the length-1-region =
per-base-unpaired identity (MCC-010) confirms they agree to 1e-9. The scalar MFE comparison (MCC-004)
confirms the partition DP and the MFE DP explore the same structure space.

### Test quality audit (HARD gate)
- Existing GAAAC tests (MCC-001/002) lock the exact analytic Z/P/pu/ΔG to 1e-12…1e-15 — true sourced
  values, not code echoes. Invariant tests (MCC-003/004) are real mathematical laws, not tautologies.
- **Defect (test-coverage gap), now FIXED in this session:** the suite relied on a SINGLE analytic
  pin (GAAAC, 3-nt loop) and omitted Stage-A edge cases the protocol mandates (single base, non-ACGU)
  plus the Bernhart monotonicity property and the two-method cross-consistency. Added:
  - **MCC-002b** `CAAAAG` — a second, independent analytic pin (4-nt loop) that exercises the
    terminal-mismatch energy path GAAAC never touches; locks exact Z/P/pu/ΔG to 1e-12.
  - **MCC-007** single base (G, A) → Z=1, pu=1, no pairs, ΔG=0.
  - **MCC-008** non-ACGU (GGGNNNCCC) → finite, deterministic, N positions unpaired w.p. 1, no pair
    involves an N.
  - **MCC-009** region accessibility monotone non-increasing in window length (Bernhart 2006).
  - **MCC-010** region(L=1 @ i) = per-base p_unpaired(i) (cross-consistency of the two methods).
  No assertions were weakened, no tolerances widened, no skips introduced. Suite now 24 tests, all
  green; full unfiltered `dotnet test Seqeron.sln -c Debug` = Failed: 0.

### Findings / defects
- **F-RNA-ACCESS-001-A (test gap, FIXED):** single-analytic-pin coverage + missing mandated edge
  cases. Resolved by adding MCC-002b/007/008/009/010 locked to externally-sourced/analytic values.
- No code defect. The implementation faithfully realises the validated description.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. State: CLEAN.**
- Both canonical methods compute the McCaskill quantities exactly on the Turner-2004 model; verified
  to floating-point precision against two independent hand-derived analytic ensembles and an
  independent brute-force enumerator, with all mathematical invariants (consistency, EFE≤MFE,
  monotone accessibility, method cross-consistency) holding. No outstanding issues.
- Note: on this machine the net10.0 SDK is at `/usr/local/share/dotnet` (10.0.301); the protocol's
  `~/.dotnet` is 8.0.422 and cannot build net10.0. Used the system SDK for build/test.
