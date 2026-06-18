# Validation Report: ONCO-CNA-002 — Focal Amplification Detection

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.DetectFocalAmplifications(segments, thresholds?)`, `OncologyAnalyzer.IdentifyAmplifiedOncogenes(amplifications)`, `OncologyAnalyzer.IsFocalAmplification(segment, thresholds)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## What the unit computes

A deterministic O(n) filter that implements the GISTIC2.0 length-based focal/broad split combined
with the GISTIC2 amplitude gate, plus an arm-level oncogene-panel mapping:

- **Amplified** ⇔ `Log2Ratio > t_amp` (default 0.1, strict).
- **Focal** ⇔ `Length / ArmLength < broad_len_cutoff` (default 0.98, strict).
- **Focal amplification** ⇔ Amplified ∧ Focal.
- **Oncogene mapping** ⇔ each panel gene whose chromosome arm carries a focal amplification:
  ERBB2 (17q), MYC (8q), EGFR (7p), CCND1 (11q), MDM2 (12q), CDK4 (12q).

## Stage A — Description

### Sources opened this session (independently retrieved)

| Source | Retrieved | Confirms |
|--------|-----------|----------|
| GISTIC2 docs — https://broadinstitute.github.io/gistic2/ (WebFetch) | 2026-06-16 | `t_amp` default **0.1**, "Regions with a copy number gain **above** this positive value are considered amplified"; `broad_len_cutoff` default **0.98**, "distinguish broad from focal events, given in units of fraction of chromosome arm"; `t_del` default 0.1 (out of scope) |
| WebSearch (GISTIC2 parameters) | 2026-06-16 | Independent corroboration: `broad_len_cutoff` 0.98 (fraction of arm), `t_amp` 0.1 (below-threshold filtered as noise) |
| Mermel et al. 2011, GISTIC2.0 — https://pmc.ncbi.nlm.nih.gov/articles/PMC3218867/ (WebFetch) | 2026-06-16 | Length-based classification: events occupying "more than 98% of a chromosome arm" are arm-level; "removing all SCNAs occupying more than 98% of a chromosome arm, leaving only the focal events"; the 98% rule is the deliberate replacement for amplitude thresholding |
| NCBI Gene 2064 (ERBB2) (WebFetch) | 2026-06-16 | Location **17q12** → arm 17q |
| NCBI Gene 4609 (MYC) (WebFetch) | 2026-06-16 | Location **8q24.21** → arm 8q |
| NCBI Gene 1956 (EGFR) (WebFetch) | 2026-06-16 | Location **7p11.2** → arm 7p |
| NCBI Gene 595 (CCND1) (WebFetch) | 2026-06-16 | Location **11q13.3** → arm 11q |
| NCBI Gene 4193 (MDM2) (WebFetch) | 2026-06-16 | Location **12q15** → arm 12q |
| NCBI Gene 1019 (CDK4) (WebFetch) | 2026-06-16 | Location **12q14.1** → arm 12q |

### Formula check
- **Focal/broad split at 0.98 of arm:** matches Mermel 2011 verbatim ("more than 98% of a chromosome
  arm" ⇒ arm-level). Implementation uses strict `<` so exactly 0.98 is arm-level — consistent with
  the paper's "more than 98% ⇒ arm-level" / "less than ⇒ focal" wording and the cutoff semantics.
- **Amplitude gate `t_amp` = 0.1, strict `>`:** matches GISTIC2 docs "**above** this positive value".
  log2 exactly 0.1 is therefore NOT amplified (boundary excluded) — confirmed correct.
- **Single-copy-gain reference value** log2(3/2) = 0.5849625… (hand-computed); 0.585 > 0.1 holds, so
  any genuine gain clears the amplitude gate. Consistent with CNVkit.
- **Oncogene→arm panel:** all six chromosome-arm prefixes match NCBI Gene cytogenetic locations
  retrieved this session.

### Edge-case semantics
- Length exactly at cutoff (0.98) → arm-level (strict `<`). Sourced (paper wording + cutoff).
- log2 exactly at t_amp (0.1) → not amplified (strict `>`). Sourced (GISTIC2 "above").
- Null input → ArgumentNullException; empty → empty; non-positive arm length / End ≤ Start →
  ArgumentException. Reasonable guards, not contradicted by sources.

### Integration assumption (documented, not invented)
GISTIC2.0 itself classifies focal-vs-broad **purely by length**; the amplitude direction/threshold for
calling a gain "amplified" is taken from the GISTIC2 `t_amp` parameter. Combining the paper's length
rule with the docs' `t_amp` into one predicate is the unit's documented integration choice
(Evidence Assumption 1 / Deviation 1). Both halves are individually source-backed. Accepted.

### Findings / divergences
None. Every Stage-A numeric (0.98, 0.1, log2(3/2)=0.585, six cytogenetic arms) traces to a source
retrieved this session. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`:
- L4160–4182 thresholds/constants (`DefaultBroadLengthCutoff` 0.98, `DefaultAmplificationLog2Threshold` 0.1).
- L4194–4206 `CopyNumberArmSegment` (`Length = End − Start`, `ArmFraction = Length / ArmLength`).
- L4217–4226 `IsFocalAmplification`: `amplified = Log2Ratio > t_amp`; `focal = ArmFraction < broad_len_cutoff`; returns conjunction after validation.
- L4239–4256 `DetectFocalAmplifications`: null guard; order-preserving filter via the predicate.
- L4268–4307 `IdentifyAmplifiedOncogenes`: builds case-insensitive arm set, emits panel genes in panel order.
- L4310–4325 `ValidateArmSegment`: ArmLength ≤ 0 and End ≤ Start throw ArgumentException.

### Formula realised correctly?
Yes. Both predicates use strict comparisons exactly as the sources require (`>` for amplitude, `<` for
focal length). The `OncogeneArms` table matches the six NCBI Gene arms verbatim. Order-preserving
subset (INV-03) and arm-set-only mapping (INV-04) hold by construction. ArmFraction is computed in
`double` after a `(double)` cast on the numerator — no integer-division pitfall.

### Cross-verification table recomputed vs code (hand-computed boundaries)

| Segment | ArmFraction | log2 | Amplified (>0.1) | Focal (<0.98) | Focal amp? | Code result |
|---------|-------------|------|------------------|---------------|------------|-------------|
| 17q 0.50, log2 1.0 | 0.50 | 1.0 | yes | yes | **yes** | reported ✓ |
| 8q 0.99, log2 1.5 | 0.99 | 1.5 | yes | no | no | dropped ✓ |
| 11q 0.98, log2 1.0 | 0.98 | 1.0 | yes | no (=cutoff) | no | dropped ✓ |
| 7p 0.30, log2 0.05 | 0.30 | 0.05 | no | yes | no | dropped ✓ |
| 12q 0.10, log2 0.585 | 0.10 | 0.585 | yes | yes | **yes** | reported ✓ |
| 17q 0.10, log2 0.1 (=t_amp) | 0.10 | 0.10 | **no** (strict >) | yes | no | dropped ✓ (new test) |

### Variant/delegate consistency
`DetectFocalAmplifications` delegates to `IsFocalAmplification`; both now have direct test coverage and
agree. Default thresholds and custom-threshold overrides (both amplitude and length) verified.

### Test quality audit (HARD gate)
- **Sourced, not code-echoed:** every expected value (0.98, 0.1, 0.585, six oncogene arms) is asserted
  against an external source, not the implementation output. Boundary tests (length=0.98, log2=0.1)
  would fail a non-strict implementation — they are not green-washes.
- **No weakened assertions:** oncogene tests use exact `Is.EquivalentTo(...)`; counts use exact
  `Has.Count.EqualTo(...)`; no ranges/Greater/Contains where an exact value is known. No tolerances
  widened, nothing skipped/ignored.
- **Coverage gaps found and fixed this session** (added 6 tests, 19 → 25):
  1. `IsFocalAmplification` public predicate had **no direct test** (only indirect). Added 4 cases
     (focal+amplified true; arm-level false; not-amplified false; invalid-segment throws).
  2. **Amplitude boundary at exactly t_amp (log2 = 0.1)** was untested — the analogue of the existing
     length-boundary test. Added `DetectFocalAmplifications_Log2ExactlyAtTamp_NotReported`, sourced to
     GISTIC2 "above this positive value" (strict `>`).
  3. **Custom `broad_len_cutoff` override** was untested (S1 only exercised the amplitude override).
     Added `DetectFocalAmplifications_CustomBroadLengthCutoff_AdmitsLongerSegment`.
- All Stage-A branches/edge cases now exercised: both predicates, both boundaries, custom amplitude
  and custom length thresholds, null/empty/invalid-arm/End≤Start guards, order/subset preservation,
  all six oncogene arms incl. the dual-gene 12q, non-panel arm, and INV-04.

### Findings / defects
No implementation defect. The only deficiency was incomplete test coverage (untested public method +
two untested sourced boundaries/overrides), which is a Stage-B test defect; **fixed in this session**.

## Verdict & follow-ups
- **Stage A: PASS** — description fully matches GISTIC2.0 (Mermel 2011), GISTIC2 docs, NCBI Gene; all
  numerics independently re-derived this session.
- **Stage B: PASS** — code faithfully realises the validated formula with correct strict-comparison
  boundary semantics; tests strengthened to lock all sourced values and exercise every public
  method/branch.
- **Test-quality gate: PASS** (after fix) — 25 sourced tests, exact assertions, full coverage; full
  unfiltered suite `Failed: 0` (6655 passed), build 0 errors.
- **End-state: CLEAN** — no implementation change required; test coverage gap completely fixed.
