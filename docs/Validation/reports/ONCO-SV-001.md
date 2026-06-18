# Validation Report: ONCO-SV-001 — Somatic Complex Rearrangement Classification (Chromothripsis Inference)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CountCopyNumberStateOscillations(IReadOnlyList<int>)`,
  `OncologyAnalyzer.TestBreakpointClustering(IReadOnlyList<long>)`,
  `OncologyAnalyzer.ClassifyComplexRearrangement(ComplexRearrangementInput)`
  (records/enums `ComplexRearrangementResult` / `ComplexRearrangementType` / `ChromothripsisConfidence` /
  `BreakpointClusteringResult`) — `src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:7192-7432`
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS (after fixing one green-washed assertion)

## Scope note

ONCO-SV-001 is **not** a generic DEL/DUP/INV/TRA/INS breakend-orientation classifier (that is
`StructuralVariantAnalyzer`, out of scope). It is the oncology-specific *complex-rearrangement* layer that
recognises the **chromothripsis** pattern: clustered breakpoints + oscillating copy number between a small
number of states. Confirmed from the implementation, doc, TestSpec, and Evidence.

## Stage A — Description

### Sources opened (independently retrieved this session)

1. **Cortés-Ciriano et al. 2020, Nat Genet 52:331–341 — PMC7058534** (WebFetch). Verbatim extracted:
   - high-confidence calls "display oscillations between two states in **at least seven adjacent segments**";
   - low-confidence calls "involve **between four and six segments**";
   - focal events "comprising **fewer than six SVs**" were excluded;
   - canonical events: "more than **60%** of the CN segments in the affected region oscillated between two states".
2. **Maher & Wilson 2012 "Chromothripsis and beyond" — PMC3861665** (WebFetch), enumerating Korbel & Campbell 2013:
   - hallmark criteria (clustering of rearrangements; interspersed retention/LOH; single chromosome/single
     parental haplotype; random fragment order; random fragment orientations; invariant head/tail alternation
     when walking the chromosome);
   - "The null hypothesis of random breakpoints predicts that the distance between breakpoints should be
     **distributed exponentially**.";
   - first-pass operational thresholds "would require, say, **10, 20, or 50** oscillating copy number changes".
3. **Magrangeas et al. 2011, Blood 118(3):675–678** (WebSearch, ASH/PubMed). Verbatim definition of an oscillating
   CN change and the canonical worked example: chromothripsis is inferred when "**at least ten switches between
   two or three copy number states** are apparent on an individual chromosome, such as a sequence of the states
   '2' and '1' (**'2; 1; 2; 1; 2; 1; 2; 1; 2; 1; 2'**)". This is the exact 11-segment / 10-switch example the
   unit uses (M1/M4).
4. **Exponential distribution CV = 1** (WebSearch, Wikipedia *Coefficient of variation*): for the exponential
   distribution σ = μ, hence CV = 1; it is the boundary between low-variance (CV<1) and over-dispersed (CV>1).

### Formula / threshold check (every value traced to a source above)

| Constant in code | Value | Source | Verdict |
|---|---|---|---|
| `MinOscillatingCopyNumberChanges` | 10 | Magrangeas 2011 "at least ten switches" / Maher&Wilson "10,20, or 50" (lowest) | ✅ |
| `MaxChromothripsisCopyNumberStates` | 3 | Magrangeas 2011 "between **two or three** copy number states" | ✅ |
| `MinChromothripsisSvBurden` | 6 | Cortés-Ciriano 2020 (focal "<six SVs" excluded) | ✅ |
| `HighConfidenceOscillatingSegments` | 7 | Cortés-Ciriano 2020 "≥ seven adjacent segments" | ✅ |
| `LowConfidenceOscillatingSegments` | 4 | Cortés-Ciriano 2020 "between four and six segments" | ✅ |
| Exponential-null CV | 1.0 | Wikipedia CV / exponential (σ=μ) | ✅ |

Oscillation = adjacent per-segment CN-state transition. Magrangeas defines an oscillation as a **switch between
states** with a bounded (two-or-three) alphabet, so the code's "transition count + distinct-state ≤ 3 gate" is a
faithful realisation of the screened quantity. The 11-segment `2,1,…,2` worked example gives exactly 10 switches,
matching the source verbatim.

### Edge-case semantics

- Empty / single segment → 0 oscillations → NotComplex (INV-1 boundary). Defined and sourced.
- Monotone many-state profile → rejected by distinct-state ≤ 3 gate (two-/three-state hallmark). Sourced.
- < 6 SVs → excluded (Cortés-Ciriano focal floor). Sourced.
- < 3 breakpoints → CV undefined (<2 gaps) → IsClustered = false. Reasonable; not a sourced literature value
  but a transparent contract (documented as ASM-02 / Open Question 1).

### Independent cross-check (hand computations)

- M1/M4 profile `2,1,…,2` (11 seg) → 10 switches: matches Magrangeas verbatim example.
- S3 gaps {1,3}: mean 2, population var ((1−2)²+(3−2)²)/2 = 1, sd 1, CV = **0.5** (Python `0.5`).
- M9 gaps {100,100,100,100}: var 0, CV = **0** → not clustered.
- M10 gaps {1,1,1,997}: mean 250, population var 186003, CV = √186003/250 = **1.7251226043386019** (Python).

### Stage A findings / divergences (notes)

- **PASS-WITH-NOTES.** The two-state hallmark is canonically two states; Magrangeas explicitly tolerates "two or
  three", which the code encodes as `MaxChromothripsisCopyNumberStates = 3` — correct. The breakpoint-clustering
  caller is a **transparent CV-vs-exponential-null proxy**, not a formal goodness-of-fit (KS/χ²) test; the
  primary sources fix only the exponential null, not one statistic. Documented as ASM-02. Cortés-Ciriano's
  ">60% of segments oscillating" canonical-fraction criterion is **not** evaluated by the classifier (only the
  ≤3-state gate). Neither is a defect against the stated, sourced scope, but both are genuine simplifications and
  are why Stage A is PASS-WITH-NOTES rather than PASS. Criteria C–F (LOH, haplotype, fragment randomness,
  derivative walk) are explicitly out of scope.

## Stage B — Implementation

### Code path reviewed

`OncologyAnalyzer.cs:7311-7325` (oscillation count), `7338-7368` (clustering CV), `7387-7432` (classifier),
constants `7192-7229`.

- Oscillation count: counts `states[i] != states[i-1]` for i∈[1,n) → [0, n−1]. Correct, O(n).
- Clustering: sorts, computes inter-breakpoint gaps, **population** sd (variance ÷ gapCount), CV = sd/mean,
  `clustered = cv > 1.0`. < 3 breakpoints or mean ≤ 0 → not clustered. Matches the validated description.
- Classifier: `isChromothripsis = distinctStates∈[2,3] ∧ oscillations ≥ 10 ∧ SVcount ≥ 6`; confidence from
  oscillating-segment count (k transitions → k+1 segments; ≥7 High, [4,6] Low, else None). All gates trace to
  sourced thresholds.

### Cross-verification table recomputed vs code (full suite run)

| Case | Input | Expected (sourced/hand) | Code |
|---|---|---|---|
| M1 | `2,1,…,2` (11) | 10 oscillations | 10 ✅ |
| M2 | `2,2,2,2` | 0 | 0 ✅ |
| M3 | `2,3,4,5,6` | 4 transitions | 4 ✅ |
| M4 | profile + 12 SV | Chromothripsis, 2 states | ✅ |
| M5 | `2,1,2,1,2,1` | NotComplex (5<10) | ✅ |
| M6 | `2..12` monotone | NotComplex (11 states) | ✅ |
| M7 | profile + 5 SV | NotComplex (<6 SV) | ✅ |
| M8 | profile | High (11 seg ≥7) | ✅ |
| S1 | `2,1,2,1,2,1` | Low + NotComplex | ✅ |
| S2 | `2,1,2,3,…` (3 states) | Chromothripsis | ✅ |
| M9 | {100..500} | CV 0, not clustered | ✅ |
| M10 | {0,1,2,3,1000} | CV 1.7251226043386019, clustered | ✅ |
| S3 | {0,1,4} | mean 2, CV 0.5, not clustered | ✅ |
| C1/C3/C4 | null | ArgumentNullException | ✅ |
| C2 | []/[2] | 0 | ✅ |

### Variant/delegate consistency

`ClassifyComplexRearrangement` reuses `CountCopyNumberStateOscillations`; oscillating-segment derivation and
distinct-state count are internally consistent. No `*Fast`/delegate variants.

### Test quality audit (HARD gate)

All 17 TestSpec cases (M1–M10, S1–S3, C1–C4) are present across 16 `[Test]` methods. Exact sourced values are
pinned (`Is.EqualTo`), edge/error cases covered (null/empty/single/<3 breakpoints/SV-floor/monotone).

**Defect found & fixed (test-only, 0 production change): green-washed M10.** The MUST clustering test asserted
`CoefficientOfVariation Is.GreaterThan(1.0)` and did **not** check `MeanGap`, although the CV of gaps {1,1,1,997}
is an **exactly computable** value (√186003/250 = 1.7251226043386019). Per the gate ("no Greater/AtLeast where an
exact value is known"), this would pass against a wrong implementation (any CV>1, or a wrong mean). Rewrote M10 to
assert the exact `MeanGap = 250` and `CoefficientOfVariation = √186003/250` `Within(1e-10)` (mirroring S3's exact
style) plus `IsClustered = true`. The value is hand-derived from the population-CV definition, not the code.

No other assertion is weakened: M9 and S3 already lock exact CVs; all classification cases assert exact enum
values; null cases assert the exact exception type. No skips, no widened tolerances, no expected-value-to-output
adjustments.

### Stage B findings / defects

- One green-washed assertion (M10) fixed in-session; no code defect found. Implementation faithfully realises the
  validated, sourced description.

## Verdict & follow-ups

- **Stage A: PASS-WITH-NOTES** — all six numeric thresholds independently confirmed; notes: CV proxy (not a formal
  GoF test) and the unevaluated ">60% canonical fraction" / criteria C–F are documented simplifications, in scope.
- **Stage B: PASS** — code matches the validated description; one green-washed test rewritten to an exact sourced
  value.
- **Test-quality gate: PASS** (after fix).
- **End-state: CLEAN.** Build 0 errors; full unfiltered suite 6694 passed / 0 failed.
- **Follow-up (optional, not a defect):** add the Cortés-Ciriano ">60% of segments oscillate between two states"
  canonical-fraction signal and/or a formal KS/χ² clustering test if a clinical-grade caller is ever required.
