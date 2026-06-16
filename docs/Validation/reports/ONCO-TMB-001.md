# Validation Report: ONCO-TMB-001 — Tumor Mutational Burden (TMB)

- **Validated:** 2026-06-16   **Area:** Oncology
- **Canonical method(s):** `OncologyAnalyzer.CalculateTMB(int mutationCount, double targetRegionMb)`;
  `OncologyAnalyzer.CalculateTMB(IEnumerable<SomaticCall>, double)` (delegate, counts Somatic);
  `OncologyAnalyzer.ClassifyTMB(double tmb)` (Low/High at FDA ≥10 cutoff)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (independently retrieved)

1. **FDA Approval Summary — Pembrolizumab for TMB-High Solid Tumors** (Marcus et al. 2021;
   PMC8416776), via WebSearch. Confirms: FDA approved pembrolizumab **June 16, 2020** for
   **TMB-H = ≥10 mutations/megabase (mut/Mb)** solid tumors, companion diagnostic FoundationOne
   CDx, based on the KEYNOTE-158 trial; "the prespecified definition of tTMB-high status was at
   least 10 mutations per megabase." The cutoff is **inclusive** (≥ 10).
2. **Wikipedia — Tumor mutational burden** (en.wikipedia.org, cited primary refs), via WebFetch.
   Confirms: TMB = "the number of non-inherited [somatic] mutations per million bases (Mb) of
   investigated genomic sequence"; FDA cutoff **≥10 mut/Mb (inclusive)**; mutation-counting
   policy is **platform-dependent** (F1CDx counts synonymous but excludes hotspot drivers;
   MSK-IMPACT includes synonymous + drivers); coefficient of variance rises sharply when panel
   size drops below ~1 Mb.
3. Cross-referenced with the repo Evidence doc's primary sources — Chalmers et al. 2017
   (Genome Medicine 9:34; "TMB was defined as the number of somatic, coding, base substitution,
   and indel mutations per megabase of genome examined"; 315-gene panel = 1.1 Mb) and the FoCR
   Harmonization review (PMC7710563; report in mut/Mb). The WebSearch result for the Frontiers/
   AACR/Lancet-Oncol KEYNOTE-158 set independently restated the same definition and ≥10 cutoff.

### Formula check

`TMB = mutationCount / targetRegionMb`, units mut/Mb. Matches Chalmers 2017 verbatim
("…per megabase of genome examined") and the FoCR harmonized reporting unit. Classification
`High ⇔ TMB ≥ 10` matches the FDA/Marcus 2021 inclusive cutoff exactly.

### Edge-case semantics

- `targetRegionMb = 0` → undefined (division by zero) → throw. Sourced: denominator is Mb.
- Sub-0.5-Mb panel → value still mathematically defined; instability is a documented limitation,
  not an error (Chalmers 2017 variance note; Wikipedia "CV rises below ~1 Mb"). Confirmed correct.
- Boundary `tmb = 10.0` → High (inclusive). Confirmed against FDA "≥10."
- Mutation counting/filtering (synonymous, germline, drivers) is platform-dependent and performed
  upstream (ONCO-SOMATIC-001); the unit correctly counts the caller-supplied somatic count.

### Independent cross-check (hand-computed, traced to sources)

| Input | Expected TMB | Source of expectation |
|-------|-------------|------------------------|
| 11 mut / 1.1 Mb | 10.0 mut/Mb | Chalmers 2017 (315-gene panel = 1.1 Mb); 11/1.1=10.0 |
| 300 mut / 30 Mb | 10.0 mut/Mb | TMB = mut/Mb (FoCR/Chalmers); 300/30=10.0 |
| 150 mut / 10 Mb | 15.0 mut/Mb | TMB = mut/Mb; 150/10=15.0 |
| 2 mut / 0.3 Mb | 6.6667 mut/Mb | TMB = mut/Mb; defined, instability documented |
| tmb 9.9 / 10.0 / 15.0 / 0.0 | Low / High / High / Low | FDA ≥10 inclusive (Marcus 2021) |

### Findings / divergences

- **Threshold conflict (pre-resolved, confirmed correct).** The Registry's old "Low <6 /
  Intermediate 6–20 / High >20" bands have **no retrievable authoritative source**; no source
  found this session defines 6 or 20. The only harmonized, FDA-approved cutoff is **≥10 mut/Mb**.
  The unit correctly implements the source-backed two-tier (Low/High @ ≥10) and does NOT fabricate
  the 6/20 bands. Registry note + TestSpec §7 already document this. No defect.
- Mutation-type policy (synonymous vs nonsynonymous) is platform-dependent per the sources; the
  unit correctly defers it upstream and computes the ratio on the supplied count. No divergence.

Stage A = **PASS**: formula, units, inclusive cutoff, and edge semantics all match independently
retrieved authoritative sources.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs`
- L1475 `TmbHighThreshold = 10.0`
- L1508–1525 `CalculateTMB(int, double)`: validates `mutationCount ≥ 0` and `targetRegionMb`
  finite & > 0 (NaN/Inf/≤0 throw `ArgumentOutOfRangeException`), returns `mutationCount/targetRegionMb`.
- L1540–1546 `CalculateTMB(IEnumerable<SomaticCall>, double)`: null-guards, counts
  `SomaticStatus.Somatic`, delegates to the scalar overload.
- L1557–1565 `ClassifyTMB(double)`: rejects negative/NaN/Inf, returns `High` iff `tmb ≥ 10`.

### Formula realised correctly?

Yes. The division is exact (`count / Mb`), no approximation; the classify uses `>=` against the
constant 10.0 (inclusive). The `SomaticStatus` enum has exactly 3 members (Somatic, Germline,
NotDetected); the overload counts only Somatic — matching "somatic mutations per Mb."

### Cross-verification table recomputed vs code (full suite run)

All five rows above pass against the actual code (verified by the green suite). Hand recompute:
11/1.1=10.0, 300/30=10.0, 150/10=15.0, 2/0.3=6.6667 — all match.

### Variant/delegate consistency

`CalculateTMB(calls, Mb)` reduces to `CalculateTMB(count, Mb)` (same validation, same division);
consistent. `ClassifyTMB` uses the shared `TmbHighThreshold` constant — single source of truth.

### Test quality audit (HARD gate)

- **Sourced, not code-echoes:** M1–M9 assert exact values traceable to Chalmers (mut/Mb) and FDA
  (≥10 inclusive). A deliberately-wrong implementation (e.g. `>` instead of `>=`, or `count*Mb`)
  would fail M7 (10.0→High) and M1 (11/1.1=10.0). Good.
- **No green-washing:** exact `Is.EqualTo(...).Within(1e-10)` for values; exact enum equality for
  classification; the boundary case `10.0` is an exact `High` assertion (not a range). No weakened
  assertions, no widened tolerances, no skips.
- **Coverage:** both `CalculateTMB` overloads + `ClassifyTMB` exercised. All Stage-A branches:
  happy path (M1–M4), zero count (M4), zero region throw (M5), negative count (S2), negative/NaN/Inf
  region (S3), small panel no-throw (S1), somatic-call counting with all 3 enum statuses (M10),
  null calls (S4), classify below/at/above/zero (M6–M9), classify invalid neg/NaN, monotonicity
  invariant sweep (C1), inclusive-boundary flip sweep (C2).
- **Gap found & fixed this session:** the empty-`SomaticCall`-collection path (0 somatic → 0 mut/Mb)
  was untested through the overload. Added `CalculateTMB_FromEmptySomaticCalls_ReturnsZero`
  (exact 0.0 over 1.1 Mb). No other gaps.
- **Honest green:** FULL unfiltered suite = **6629 passed, 0 failed** after the addition (was 6628);
  `dotnet build` 0 errors; the TMB test file builds warning-free (the 4 NUnit2007 warnings are
  pre-existing in unrelated `ApproximateMatcher_EditDistance_Tests.cs`).

### Findings / defects

None. Implementation faithfully realises the validated description; tests are sourced and now
complete. One missing edge case (empty collection) was added and locked.

## Verdict & follow-ups

- **Stage A: PASS** · **Stage B: PASS** · **End-state: CLEAN.**
- **Test-quality gate: PASS** (exact sourced values, all branches + overloads covered, empty-
  collection gap fixed, full suite green).
- No logged defects. The previously-resolved 6/20-band Registry conflict remains correctly handled
  (no fabrication). No follow-ups.
