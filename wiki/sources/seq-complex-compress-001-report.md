---
type: source
title: "Validation report: SEQ-COMPLEX-COMPRESS-001 (compression-based sequence complexity — Lempel–Ziv 1976)"
tags: [validation, complexity, governance]
doc_path: docs/Validation/reports/SEQ-COMPLEX-COMPRESS-001.md
sources:
  - docs/Validation/reports/SEQ-COMPLEX-COMPRESS-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: SEQ-COMPLEX-COMPRESS-001

The two-stage **validation write-up** for test unit **SEQ-COMPLEX-COMPRESS-001**
(compression-based sequence complexity — the Lempel–Ziv 1976 measure), validated 2026-06-16 in
the Complexity area. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's **verdict** on both the algorithm
description and the shipped code. The algorithm itself is summarized in
[[sequence-complexity-compression-lempel-ziv]]; the two-stage methodology is the
[[validation-protocol]]. Distinct from the pre-implementation
[[seq-complex-compress-001-evidence]] artifact.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: FAIL → FIXED · State: ✅ CLEAN.** Full unfiltered
`Seqeron.Genomics.Tests` = **6603 passed / 0 failed**, `dotnet build` 0 errors, changed files
warning-free (4 pre-existing NUnit2007 warnings in an unrelated file). One defect (the b<2
single-symbol normalization diverging from its own cited reference) was found and completely
fixed in-session; the test was locked to the sourced value rather than the code output.

Canonical methods: `SequenceComplexity.CalculateLempelZivComplexity(string|DnaSequence)`,
`CalculateNormalizedLempelZivComplexity(string|DnaSequence)`,
`EstimateCompressionRatio(string|DnaSequence)` (delegate → normalized LZ). See report at
`docs/Validation/reports/SEQ-COMPLEX-COMPRESS-001.md`. Logged: FINDINGS_REGISTER A29.

## Stage A — description (algorithm faithfulness)

- **Sources retrieved (not trusted from the label):** Naereen/Lempel-Ziv_Complexity reference
  source (docstring "the number of different substrings encountered as the stream is viewed
  from beginning to the end"; doctests `1001111011000010`→8, `1010101010101010`→7,
  `1001111011000010000010`→9, `100111101100001000001010`→10); Wikipedia Lempel–Ziv complexity;
  antropy `lziv_complexity` and entropy `lziv_complexity` (identical normalization block, both
  clamp base to 2); Lempel & Ziv (1976) IEEE TIT 22(1):75–81 (primary, paywalled).
- **Formula check:** raw LZ `c` = exhaustive-history left-to-right set-based parse → number of
  distinct components (matches Naereen exactly). Normalized `LZn = c / (n / log_b n)`, b =
  distinct symbols, matches antropy/entropy exactly for b ≥ 2.
- **DEFECT (Stage A divergence, corrected):** for the b<2 (single-symbol) case, the reference
  **clamps base to 2** and returns `c/(n/log_2 n)`. The TestSpec/Evidence claimed "returns the
  raw count" and mis-attributed it to source #4 — which does the opposite. Description
  corrected in TestSpec §1.3/§4 (M8)/§6 (A2) and Evidence.
- **Edge cases:** empty → 0; homopolymer `"0"×16` → `0/00/000/0000/00000` = 5; single base →
  1; b<2 normalization → clamp base to 2 (corrected from "raw count").
- **Independent cross-check (Python reimplementation of the Naereen set parser):** raw c 8/7/9/10
  for the doctests, `"0"×16`→5, `ACGT`→4, `AAAA`→2, `ACGTACGTACGTACGT`→9 all match. Normalized:
  `1001111011000010`→**2.0**; `ACGTACGTACGTACGT`→**1.125**; `"0"×16` antropy-style clamp→**1.25**
  (NOT 5.0).

## Stage B — implementation (code review + cross-check)

- **Code path:** `SequenceComplexity.cs:460-573`. `CalculateLempelZivComplexityCore` (522) is
  a set-based exhaustive-history parser identical to the Naereen reference; `EstimateCompressionRatio`
  (507/517) is a thin delegate to normalized LZ.
- **DEFECT (found & FIXED):** `CalculateNormalizedLempelZivComplexityCore` (548) computed
  `c/(n/log_b n)` for b ≥ 2 (correct) but **returned the raw count for b<2** — diverging from
  the cited reference. Fix: replaced `if (b < MinAlphabetForNormalization) return c;` with a
  clamp `if (b < MinAlphabetForNormalization) b = MinAlphabetForNormalization;`, keeping the
  `n==1 ⇒ log_b(1)=0` div-by-zero guard. For `"0"×16` the result is now **1.25**
  (= 5/(16/log₂16)), matching the reference.
- **Cross-verification table (recomputed vs fixed code, all PASS):** raw LZ doctests 8/7/9/10;
  `"0"×16`/`ACGT`/`AAAA`/`A`/`""` → 5/4/2/1/0; normalized `1001111011000010`→2.0, `"0"×16`
  (b<2 clamp)→1.25, `ACGTACGTACGTACGT` (b=4)→1.125; `EstimateCompressionRatio` → 2.0 / 1.125.
- **Variant/delegate consistency:** string and DnaSequence overloads agree (raw + normalized);
  `EstimateCompressionRatio` equals `CalculateNormalizedLempelZivComplexity` for both overloads.
- **Test-quality audit (HARD gate) PASS:** expectations are sourced, not code echoes (raw from
  Naereen doctests, normalized from antropy code; the b<2 value corrected from the unsourced 5.0
  to the sourced 1.25 — the **code** was fixed and the test locked to 1.25, no green-washing). A
  deliberately-wrong parser fails M1. Added 5 tests closing untested overloads/branches (fixture
  15→20). **MCP note:** `Seqeron.Mcp.Sequence.Tests/ComplexityCompressionRatioTests.cs` targets
  net9.0 and cannot be built by the installed net8 SDK (pre-existing environmental limit); it is
  an unchanged binding smoke test using b=4 inputs unaffected by the b<2 fix.

## Findings

- **FIXED** — b<2 normalization returned the raw count instead of the reference's clamp-to-2
  normalized value; code corrected to the clamp, tests locked to the sourced 1.25, coverage
  gaps closed.
- **Stage A NOTE** — raw-LZ description and all b ≥ 2 normalization values are correct and
  externally confirmed; the single divergence was the mis-sourced b<2 convention, now corrected.
- **End-state ✅ CLEAN** — defect completely fixed in-session; build + full suite green (6603/0).
