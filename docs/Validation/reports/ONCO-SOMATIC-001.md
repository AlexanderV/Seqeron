# Validation Report: ONCO-SOMATIC-001 — Somatic Mutation Calling

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CallSomaticMutations`, `OncologyAnalyzer.Classify`, `OncologyAnalyzer.FilterGermlineVariants`, `OncologyAnalyzer.CalculateSomaticScore`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened & what they confirm (retrieved this session)

1. **Yan et al. 2021, Sci. Rep. 11:11640** — WebFetched `nature.com/articles/s41598-021-91142-1` (via IDP redirect chain). Verbatim: *"WES has a mutation limit of detection (LoD) at variant allele frequencies (VAF) of 5%."* and *"Putative mutations called at ≤ 5% VAF are frequently due to sequencing errors, therefore reporting these subclonal mutations incurs risk of significant false positives."* Of 226 sub-5% calls, 52% were false positives on orthogonal confirmation (82% in cancer genes). → Confirms the τ_t = 0.05 tumor LoD gate and that the gate is applied first.
2. **GATK Mutect2 `mutect.tex`** — WebFetched `raw.githubusercontent.com/broadinstitute/gatk/master/docs/mutect/mutect.tex`. Verbatim: *"If we have no matched normal, ℓ_n = 1."* Germline filter compares three unnormalized probabilities — germline het in normal (∝ 2f(1−f)·ℓ_n·(1−π)), germline hom-alt (∝ f²), and somatic tumor-only (∝ (1−f)²·ℓ_t^S·π); germline error probability = normalized sum of the first two. → Confirms tumor-only mode (f_n = 0 surrogate for ℓ_n = 1) and that a variant clearly present in the matched normal is treated as germline.
3. **Strelka / Saunders 2012; Illumina somatic-caller technotes** — WebSearched. Confirms the somatic concept: somatic = present in tumor and *absent* in matched normal; Strelka's Bayesian model represents continuous allele frequencies for tumor and normal and leverages the expected (ref/ref) normal genotype; VAF = variant-supporting reads / total coverage at the site.

### Formula check
- **VAF** `f = altReads / totalReads`, `totalReads = 0 ⇒ 0` — matches GATK `AD` convention (alt AD / Σ AD) and the Illumina/Strelka definition.
- **Decision rule** `f_t < τ_t → NotDetected`; `f_t ≥ τ_t ∧ f_n ≤ τ_n → Somatic`; else `Germline` — faithfully realises "present in tumor (above LoD), absent in normal." The LoD gate first (INV-2) matches Yan 2021. Boundaries inclusive (`≥ τ_t`, `≤ τ_n`).
- **Score** `max(0, f_t − f_n)` — an explicitly-documented monotone surrogate (ASM-02), bounded [0,1] since both VAFs ∈ [0,1]. Not presented as a caller LOD.

### Edge-case semantics
Tumor-only (normal total 0 ⇒ f_n = 0 ⇒ Somatic, the ℓ_n = 1 analogue), sub-5% tumor (NotDetected), CHIP-like normal above τ_n (Germline), f_n ≥ f_t (score 0), empty input (empty). All defined and sourced.

### Independent cross-check (hand-computed, τ_t=0.05, τ_n=0.01)
| Row | f_t | f_n | Branch | Status | Score |
|-----|-----|-----|--------|--------|-------|
| A 25/100,0/100 | 0.25 | 0.00 | f_t≥τ_t, f_n≤τ_n | Somatic | 0.25 |
| B 48/100,50/100 | 0.48 | 0.50 | f_n>τ_n | Germline | 0 |
| C 2/100,0/100 | 0.02 | 0.00 | f_t<τ_t | NotDetected | 0 |
| D 30/100,3/100 | 0.30 | 0.03 | f_n>τ_n | Germline | 0 |
| E 20/100,0/0 | 0.20 | 0.00 | f_n≤τ_n | Somatic | 0.20 |
| F 5/100,0/100 | 0.05 | 0.00 | f_t=τ_t (incl.) | Somatic | 0.05 |
| G 30/100,1/100 | 0.30 | 0.01 | f_n=τ_n (incl.) | Somatic | 0.29 |

### Findings / divergences
The two ASSUMPTIONs (τ_n = 1% fixed normal ceiling; score = max(0, f_t − f_n)) are honestly labelled as simplifications of the Bayesian callers and are not claimed to be sourced caller outputs. τ_n is exposed as a configurable parameter. No biological or mathematical error. **Stage A: PASS.**

## Stage B — Implementation

### Code path reviewed
- `OncologyAnalyzer.cs:213-229` `CallSomaticMutations` (null + threshold guards, per-variant delegate, order preserved).
- `OncologyAnalyzer.cs:234-267` `Classify` (VAF computation, three-branch decision, score gated on Somatic).
- `OncologyAnalyzer.cs:284-294` `FilterGermlineVariants` (Somatic subset, order preserved).
- `OncologyAnalyzer.cs:309-327` `CalculateSomaticScore` (`max(0, f_t − f_n)`).
- `OncologyAnalyzer.cs:617-643` `CalculateVaf` (alt/total, non-negative + alt≤total guards, 0-coverage⇒0), `ValidateThreshold` ([0,1] + NaN).

### Formula realised correctly?
Yes. The decision order (LoD gate first, then normal-absence) and inclusive boundaries match Stage A exactly. Score gated to 0 for non-Somatic. VAF and validation match the contract.

### Cross-verification table recomputed vs code
All 7 hand-computed rows above reproduce the code's output and the test expectations exactly (verified by running the fixture; boundaries F/G confirm inclusivity). Score M9 = 0.25 − 0.05 = 0.20 and S3 (f_n ≥ f_t) = 0 confirmed.

### Variant/delegate consistency
`CallSomaticMutations` delegates to `Classify`; `FilterGermlineVariants` delegates to `CallSomaticMutations` then filters Status==Somatic; `CalculateSomaticScore(VariantObservation)` and the private `(f_t,f_n)` overload agree. Consistent.

### Test quality audit (HARD gate)
- **Sourced, not code-echoed:** statuses trace to the somatic-state definition (Strelka), the 5% LoD (Yan 2021), and the matched-normal rule (Mutect2); scores to the documented max(0,f_t−f_n). Exact values (`Within(1e-10)`), not ranges. No tautologies.
- **No green-washing:** no assertion weakened, no tolerance widened, no skip/ignore, no expected value bent to actual.
- **Coverage gap found & fixed:** the non-negative-read guard (contract §3.3) was exercised only on the sibling `CalculateVAF`, not on this unit's `Classify`/`CallSomaticMutations` path; and the public single-variant `Classify` overload's three branches were tested only indirectly. Added two tests:
  - `Classify_NegativeReadCounts_Throw` — negative tumor alt, normal alt, and tumor total each throw `ArgumentOutOfRangeException`.
  - `Classify_SingleVariant_CoversAllThreeBranches` — Somatic / Germline / NotDetected directly via `Classify` with sourced thresholds.
  Fixture 18 → 21.
- **Honest green:** full unfiltered `dotnet test` = **6621 passed, 0 failed**; `dotnet build` 0 errors. The 4 build warnings are pre-existing in unrelated `ApproximateMatcher_EditDistance_Tests.cs`; the changed test file is warning-free.

### Findings / defects
No code defect; no description defect. One test-coverage gap (negative-read guard + direct `Classify` branches) — fixed in this session. Logged as **A41** in `FINDINGS_REGISTER.md`. **Stage B: PASS.**

## Verdict & follow-ups
- **Stage A: PASS** — biology/maths correct and sourced; simplifications honestly labelled.
- **Stage B: PASS** — code faithfully realises the validated rule; tests now cover all four public methods, all three classification branches, the boundaries, the score formula, and all documented error cases (null, empty, bad threshold, alt>total, negative reads).
- **End-state: CLEAN.** Test-quality gate: PASS (after adding the two missing-coverage tests).
