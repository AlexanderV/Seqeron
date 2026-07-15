---
type: source
title: "Validation report: DISORDER-LC-001 (protein low-complexity region detection — SEG / Wootton & Federhen, DisorderPredictor.PredictLowComplexityRegions)"
tags: [validation, analysis, governance]
doc_path: docs/Validation/reports/DISORDER-LC-001.md
sources:
  - docs/Validation/reports/DISORDER-LC-001.md
source_commit: c9ed6cf3055a7708deeca143f62df61bea0e7263
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: DISORDER-LC-001

The two-stage **validation write-up** for test unit **DISORDER-LC-001** (protein **low-complexity region
detection** — partitioning a protein into compositionally biased low-complexity vs high-complexity
segments via the **SEG algorithm** of Wootton & Federhen), area **ProteinPred**, validated 2026-06-16.
This is the *report* artifact that feeds one row of the [[validation-ledger]]; it records the validator's
independent **verdict** on both the algorithm description (Stage A) and the shipped code (Stage B), inside
the wider [[validation-and-testing]] campaign. The algorithm itself — the Shannon-entropy bits/residue
complexity measure, the two-stage trigger/extend scan, the W=12 / K1=2.2 / K2=2.5 defaults, the oracle
window entropies, and the documented repository extensions — is synthesized in the concept
[[protein-low-complexity-seg]]. [[test-unit-registry]] defines the unit. Despite the `DISORDER-LC` prefix
and the host class name `DisorderPredictor`, this unit is **SEG low-complexity detection**, a *distinct
algorithm* from [[intrinsic-disorder-prediction-top-idp|intrinsic-disorder prediction]] (the TOP-IDP
`PredictDisorder` anchor). Distinct from the pre-implementation evidence artifact
[[disorder-lc-001-evidence]] (sourced from `docs/Evidence/`) — this page is the independent two-stage
re-validation verdict. Source under test:
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/DisorderPredictor.cs`.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: CLEAN / FIXED.** No implementation defect; no
production-code change. The biology/maths is correct and the entropy formula plus the W/K1/K2 defaults are
confirmed against **two independent reference sources**; the code faithfully realises the validated
formula, with every expected test value independently reproduced from the spec. Test-quality gate **PASS**
— three genuine coverage gaps closed with exact sourced values. Full suite **6612 passed / 0 failed**
(1 pre-existing unrelated skipped benchmark); `dotnet build` 0 errors.

## Canonical methods & source under test

In `DisorderPredictor.cs`:

- `PredictLowComplexityRegions(string, int triggerWindow, double triggerThreshold, double extensionThreshold, int minLength)` (`:497`)
  — null-guard → upper-case → short-circuit if `len < W` → Phase 1 trigger (window entropy `≤ K1`) →
  collect spans → Phase 3 greedy extension (whole growing segment entropy `≤ K2`) → merge adjacent
  (`≤ lastEnd+1`) → `minLength` filter + classify.
- `CalculateShannonEntropy` (`:280`, private) — counts A–Z, `p = count/length`, `H = −Σ p·log₂p`; exact
  Shannon entropy normalized by full window length.
- `ClassifyLowComplexityType` (private) — the `"X-rich"` / `"X/Y-rich"` presentation label (>50% dominant
  residue rule).
- Defaults bound to `SegTriggerWindow = 12 / SegTriggerComplexity = 2.2 / SegExtensionComplexity = 2.5`
  (verified verbatim against NCBI). Single public overload; **no `*Fast` variant**.

## Stage A — description (biology faithfulness)

Grounded against sources opened this session, independent of the repo's own Evidence:

- **GCG/Weizmann SEG help** and the **`ncbi-seg` Ubuntu manpage** — confirmed `-WINdow=12`,
  `-LOWcut=2.2` (K1), `-HIGhcut=2.5` (K2); complexity in **bits/residue**; max amino-acid complexity
  `log₂(20) = 4.322`; Stage-1 trigger = "complexity ≤ cutoff", Stage-2 extension merges overlapping
  windows `≤ K2`.
- **NCBI `blast_seg.c`** — constants `kSegWindow = 12`, `kSegLocut = 2.2`, `kSegHicut = 2.5` verbatim;
  carries both an `lnfact[]` factorial-log table (the multinomial WF eq-3/5 complexity used for
  optimization) **and** a Shannon-entropy `s_Entropy` trigger measure.
- **WebSearch (Wootton & Federhen 1993, Comput. Chem. 17(2):149–163)** — confirmed the paper defines
  "local compositional complexity" and the caveat that Shannon entropy "is not a good complexity measure
  for protein sequences" (the reason real SEG adds the P0 optimization).

### Formula check (independent hand computation, L=12)

`H = −Σ pᵢ·log₂pᵢ`, `pᵢ = nᵢ/L`. Reference values hand-derived and matched:

| Composition | Distinct | H (bits) |
|---|---|---|
| `QQQ…` (12×1) | 1 | 0.000000 |
| `AAAAAALLLLLL` (6+6) | 2 | 1.000000 |
| `AAABBBCCCDDD` (3+3+3+3) | 4 | 2.000000 |
| `ACDEFGHIKLMN` (12 distinct) | 12 | 3.584963 = log₂12 |
| 20 equal | 20 | 4.321928 = log₂20 |

Edge-case semantics all defined and source-backed: `len < W` → empty; window of W distinct residues
(`H = log₂W ≈ 3.585 > K2`) → never flagged; homopolymer (`H = 0`) → always triggers.

### Stage A notes (→ PASS-WITH-NOTES, documented, non-blocking)

1. **Complexity measure is the Shannon-entropy form, not WF equation (3) verbatim.** Real SEG's eq (3) is
   the multinomial complexity `K₂ = (1/L)·log_N(N!/Πnᵢ!)`, further reduced by the P0 significance test
   (eq 5). The repo uses the bits/residue Shannon entropy that NCBI `s_Entropy` and the GCG "bits/residue"
   docs describe, with the same 2.2/2.5 calibrated cutoffs, and **explicitly documents** the omission of
   WF eq-3 + P0 (algorithm doc §5.3 "Intentionally simplified / Not implemented"). A legitimate, disclosed
   simplification — boundaries can differ from NCBI SEG on mixed-complexity edges but match on
   homopolymer/biased inputs.
2. **`Type` label (`"X-rich"` / `"X/Y-rich"`)** is a repository presentation extension (>50% dominant
   rule), not part of SEG — registered as Assumption #1; does not affect coordinates.

## Stage B — implementation

The detector was **re-implemented independently in Python** (formula + spec only, not reading the C#) and
reproduced **every** expected value in the test file:

| Case | Input | Expected | Independent recompute |
|---|---|---|---|
| M1 | 26×Q | (0,25) | (0,25) ✓ |
| M2 | all-20 ×2 | empty | empty ✓ |
| M3 | `AAABBBCCCDDD` | (0,11) | (0,11) ✓ |
| M3′ | same, K1=0.5 | empty | empty ✓ |
| M5 | 12A+12L | (0,23) | (0,23) ✓ |
| M6 | 20Q+60 spacer+20A | (0,34),(67,99) | (0,34),(67,99) ✓ |
| S1 | 12Q, minLen=15 | empty | empty ✓ |
| S2 | 6Q | empty | empty ✓ |
| S3 | 12Q | (0,11) | (0,11) ✓ |

**M6 boundaries verified as genuine threshold crossings** (not code echoes): `H(0..34)=2.4939 ≤ 2.5 <
H(0..35)=2.6077`; `H(67..99)=2.4250 ≤ 2.5 < H(66..99)=2.5452`; spacer windows `H=3.585 > K2` (gap
preserved). These trace the test's expected values to the sourced formula, not to the code.

### Test-quality audit (HARD gate) — PASS

Pre-existing fixture had 18 exact-boundary tests (not code echoes). Three coverage gaps closed this
session, each with **exact** sourced/hand-derived values (no `Greater`/`Contains`/range/widened tolerance,
nothing skipped):

- **M8** — dipeptide segment `(0,23)`, A=L=50% → asserts the `"A/L-rich"` label branch (no residue >50%).
- **M9** — custom `triggerWindow` W=4 on `AAABBBCCCDDD` (length-4 windows `H ∈ {0.811,1.0} ≤ K1`) → (0,11).
- **M10** — custom `extensionThreshold` K2=2.0 on the M6 input → (0,28),(71,99); boundaries verified as
  real crossings (`H=1.877 ≤ 2.0 < 2.026`).

Fixture now **21 tests**; all MUST/SHOULD/COULD cases + INV-1/INV-2 properties + null/empty error paths +
both label branches + all three tunable parameters covered. Full unfiltered suite **6612 passed /
0 failed**; changed test file builds warning-free.

## Findings & follow-ups

- **No code defect, no production-code change (State CLEAN / FIXED).** Every expected value independently
  reproduced from the spec; the only "deviations" are the disclosed Shannon-entropy-vs-WF-eq3
  simplification and the cosmetic `Type` label — both documented in the algorithm doc + Assumption
  Register, neither affects the validated coordinate contract.
- **Test coverage gaps fully fixed** (X/Y-rich label, custom `triggerWindow`, custom `extensionThreshold`).
- **No outstanding follow-ups.** A future enhancement — adding the WF eq-3 complexity + P0 optimization for
  NCBI-grade boundaries — is out of scope and already noted as a documented limitation, not a defect.
  Research-grade, not for clinical use.
</content>
