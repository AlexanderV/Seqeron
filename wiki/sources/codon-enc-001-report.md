---
type: source
title: "Validation report: CODON-ENC-001 (Effective Number of Codons — ENC/Nc, CodonUsageAnalyzer.CalculateEnc)"
tags: [validation, annotation]
doc_path: docs/Validation/reports/CODON-ENC-001.md
sources:
  - docs/Validation/reports/CODON-ENC-001.md
source_commit: 816a85f76e86d111265c4d6db5e02b68b16f7c07
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CODON-ENC-001

The two-stage **validation write-up** for test unit **CODON-ENC-001** (Effective Number of Codons —
Wright's reference-free codon-usage-bias number in [20, 61]), validated 2026-06-15. This is the
*report* artifact that feeds one row of the [[validation-ledger]]; it records the validator's
independent **verdict** on both the algorithm description (Stage A) and the shipped code (Stage B),
and the wider campaign is [[validation-and-testing]]. The algorithm, its equations, invariants,
oracles and corner-case rules are synthesized in the concept [[effective-number-of-codons]] (the
ENC/Nc anchor in the codon-usage-bias family, built on the [[relative-synonymous-codon-usage]]
synonymous-family frequencies `p_i`); [[test-unit-registry]] defines the unit. Distinct from
[[codon-enc-001-evidence]] — the pre-implementation evidence artifact sourced from `docs/Evidence/` —
this page is the independent two-stage validation verdict.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS-WITH-NOTES · End state: ✅ CLEAN.** The core biology/maths
is correct and sourced, and every defect found was either fixed or fully documented — no half-fix.
Full unfiltered `dotnet test` = **6527 passed, 0 failed** (1 unrelated benchmark skipped); `dotnet
build` 0 errors (4 pre-existing unrelated warnings). The two "NOTES" are documented edge-case
divergences on **degenerate input that real coding sequences never trigger** (see Findings).

## Canonical method & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs:274–360`:

- `CalculateEnc(string)` (core) and `CalculateEnc(DnaSequence)` (null-checks then delegates to the
  string core, `:274–277`); private `CalculateEncCore`. `GetStatistics` reuses `CalculateEncCore`
  (`:393`) — one canonical entry, no divergent re-implementation.
- Eq. (1) codon homozygosity at `:312–318` — `F = (n·Σp² − 1)/(n − 1)`, `p = nᵢ/n`. `n ≤ 1` skipped
  (`:309`); singlets (Met/Trp) skipped (`:306`).
- Eq. (4) within-class averaging via `AverageOrNull` (`:330–355`) — averages only *estimable*
  members. Eq. (5a) isoleucine fallback `F̂₃ = (F̂₂ + F̂₄)/2` at `:337–338`, firing only when F₃
  unestimable *and* F₂,F₄ present. Eq. (3) aggregation with constant `2` at `:343–347`; cap to 61 at
  `:349`. Lower clamp `Math.Max(20, …)` at `:349`.
- `ClassContribution` (`:357–360`) — the divergent path (defect B1, below).

## Stage A — description (algorithm faithfulness)

Confirmed **verbatim** against **Fuglsang (2004) "The 'effective number of codons' revisited", BBRC
317:957–964** (PDF re-fetched, `pdftotext -layout`): Eq. (1) `F̂ = (n·Σpᵢ² − 1)/(n − 1)`; Eq. (3)
`N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆`; Eq. (4) missing-**member** rule (average only the estimable
members of the class); Eq. (5a) isoleucine fallback; the **cap to 61** ("Wright recommends
re-adjusting the result down to 61"); and the 20–61 range. The **codonW reference implementation**
(Peden thesis, §Eq. 2-7/2-8) confirms the formula and, decisively, the **whole-class-empty rule**:
for a sequence where a whole synonymous family is empty (`F̂ₙ = 0`), *"Nc is not calculated"* (gene
assumed too short / extremely skewed), with the isoleucine 3-fold exception. The standard-code
degeneracy partition (NCBI table 1, cross-checked vs the in-code `CodonToAminoAcid` map): 9 doublets,
1 triplet (Ile), 5 quartets, 3 sextets (Leu/Ser/Arg), 2 singlets (Met/Trp) — matches Eq. (3)
constants `2, 9, 1, 5, 3`.

Two documented divergences (both on input real coding sequences never produce):

- **Note A1 (whole-class-absent).** The description originally prescribed "absent class contributes
  its full codon count" for an *entirely empty* 2/4/6-fold class — diverging from codonW's "Nc is
  not calculated" — and even listed that very behaviour as a self-contradictory failure-mode. Spec
  updated this session (§7) to record the divergence honestly.
- **Note A2 (lower clamp at 20).** codonW caps only the top (61); the lower clamp at 20 is a
  defensive bound, not a Wright/codonW instruction. Documented, unchanged.

## Stage B — implementation

The formula is realised correctly (Eq. 1/3/4/5a, cap at 61) and reproduces the independent reference
to **full double precision** on every fully-populated input. Independent Python re-implementation of
the Wright/codonW algorithm (run this session) cross-verification table vs the C# `CalculateEnc`:

| Case | Input | Reference Nc | C# | Match |
|------|-------|--------------|----|-------|
| M1 one-codon-per-aa | each aa ×2, single codon | 20.0 | 20.0 | ✅ |
| M2 near-uniform | all 61 sense codons ×2 | >61 → cap 61.0 | 61.0 | ✅ |
| M3 fully-populated biased | all classes estimable | 41.288461538461526 | 41.288461538461526 | ✅ |
| M5 Ile-absent (Eq. 5a) | all classes but Ile | 39.47394540942927 | 39.47394540942927 | ✅ |
| C1 uniform 2:1 bias | F̂=1/3 every class | 56.0 | 56.0 | ✅ |
| M5b whole-class-absent | only Phe | (undefined per codonW) | 29.0 (library convention) | n/a — divergence pinned |

**Defect B1 (logged, low severity).** `ClassContribution` (`:357–360`): when a whole degeneracy
class has no estimable F̂, it returns the raw codon count (full count) instead of declining to
compute — diverging from codonW's "Nc not calculated". **Reachable only on genes missing an entire
2/4/6-fold class** (never on real coding sequences). A complete fix (return a sentinel/NaN for "Nc
undefined") would change the `double`-returning contract and ripple into `GetStatistics` and the MCP
`EncResult` — a product/contract decision out of scope for one validation session. Documented (spec
§7, FINDINGS_REGISTER) and **pinned by an explicitly LIBRARY-SPECIFIC / NOT-Wright-labelled test
(M5b = 29.0)** rather than masquerading as sourced.

**Test-quality audit (HARD gate, PASS).** Two **code-echo defects fixed**: the original M3
(expected 29.0) and M5 (expected 40.4) asserted values produced by the *unsourced* full-count
fallback (B1) — they would have passed against the very behaviour the spec names as a defect.
**Rewritten** to sourced exact values on fully-populated genes (M3 = 41.288461538461526;
M5 = 39.47394540942927), each traced to the independent reference and each genuinely exercising Eq.
1/3/4 (M5 also Eq. 5a with F₂/F₄/F₆ all estimable). The divergence is pinned by the honest,
explicitly-labelled M5b test. No weakening (no tolerances widened, no assertions softened, no tests
skipped); M1/M2/C1, M4 range-invariant property test, and M6 (null→throws)/M7 (empty→0)/S1
(case)/S2 (invalid codons)/S3 (delegate) all retained. Coverage spans both public overloads, all
formula paths, cap-at-61, the range invariant, null/empty/invalid/case edges, and the divergent
absent-class path.

## Findings & follow-ups

- **End state ✅ CLEAN.** Every defect found was fixed (test code-echoes → sourced exact values;
  added the honest divergence test) or fully documented (B1). Algorithm is fully functional for all
  real coding sequences.
- **Follow-up (optional, not blocking) — FR:** decide a contract for "Nc undefined" (whole class
  empty) to match codonW semantics consistently across `CalculateEnc` / `GetStatistics` / the MCP
  `EncResult`. See [[effective-number-of-codons]] for the standing "lower clamp at 20" nuance
  (defensive bound, not a Wright-prescribed parameter — Note A2).
