# Test Specification: CODON-STATS-001

**Test Unit ID:** CODON-STATS-001
**Area:** Codon
**Algorithm:** Codon Usage Statistics (GetStatistics, CalculateCai, EColiOptimalCodons, HumanOptimalCodons)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Sharp & Li (1987), CAI, NAR 15:1281–1295 | 1 | https://doi.org/10.1093/nar/15.3.1281 | 2026-06-13 |
| 2 | Peden (1999) CodonW thesis (GC3s definition §1.8.2.1.3) | 1–3 | https://codonw.sourceforge.net/JohnPedenThesisPressOpt_water.pdf | 2026-06-13 |
| 3 | CodonW codon-usage indices (CAI exclusions) | 3 | https://codonw.sourceforge.net/Indices.html | 2026-06-13 |
| 4 | seqinr `cai` documentation | 3 | https://search.r-project.org/CRAN/refmans/seqinr/html/cai.html | 2026-06-13 |
| 5 | EMBOSS `cusp` (GC1/GC2/GC3) | 3 | https://www.bioinformatics.nl/cgi-bin/emboss/help/cusp | 2026-06-13 |
| 6 | Biopython v1.79 `SharpEcoliIndex` | 3 | https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/SeqUtils/CodonUsageIndices.py | 2026-06-13 |
| 7 | Kazusa H. sapiens [gbpri] codon usage | 5 | https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606 | 2026-06-13 |

### 1.2 Key Evidence Points

1. `w_i = f_i / max(f_j)` over the synonymous family; `CAI = exp[(1/L) Σ ln w_i]` — Sharp & Li 1987 / Wikipedia.
2. CAI excludes non-synonymous (single-codon Met, Trp) and termination codons — seqinr; CodonW Indices.
3. GC3s = "frequency of G or C at the third position of synonymous codons (i.e. excluding Met, Trp and termination codons)" — Peden §1.8.2.1.3.
4. GC1/GC2/GC3 = GC content at codon positions 1/2/3 ("1st/2nd/3rd letter GC") — EMBOSS cusp.
5. E. coli reference w values (CTG=1, GCC=0.122, CGT=1, AGG=0.002 …) — Biopython SharpEcoliIndex.
6. Human RSCU derived from Kazusa per-thousand frequencies (CTG≈2.3713, GCC≈1.5988, Met/Trp=1.0).

### 1.3 Documented Corner Cases

- A sequence of only Met/Trp/stop codons has no scorable codon → CAI 0 and GC3s 0 (empty denominator). Source: seqinr / Peden.
- Empty / null input → zeroed statistics or `ArgumentNullException` (implementation contract; not specified by the literature).

### 1.4 Known Failure Modes / Pitfalls

1. Treating GC3s as GC at *all* third positions (including Met/Trp/stop) — contradicts Peden §1.8.2.1.3.
2. Including Met/Trp/stop codons in CAI — contradicts Sharp & Li 1987 / CodonW.
3. `ln(0)` when a codon's w is 0 — handled by skipping zero-w codons (deviation from the 0.01 floor of seqinr/Bulmer).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `GetStatistics(string)` / `GetStatistics(DnaSequence)` | CodonUsageAnalyzer | Canonical | Aggregates counts, RSCU, ENC, GC1/2/3, GC3s, total codons |
| `CalculateCai(string, refRscu)` / `CalculateCai(DnaSequence, refRscu)` | CodonUsageAnalyzer | Canonical | Sharp & Li 1987 CAI |
| `EColiOptimalCodons` (property) | CodonUsageAnalyzer | Reference | Sharp & Li 1987 w table |
| `HumanOptimalCodons` (property) | CodonUsageAnalyzer | Reference | Kazusa-derived human RSCU |
| `CodonUsageStatistics.OverallGc` | record | Internal | (GC1+GC2+GC3)/3 |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | `0 ≤ CalculateCai ≤ 1` for any sequence and reference | Yes | Sharp & Li 1987 (geometric mean of w∈[0,1]) |
| INV-2 | CAI of an all-optimal sequence (each codon = family's w-max) = 1.0 | Yes | Sharp & Li 1987 |
| INV-3 | GC3s ignores Met (ATG), Trp (TGG) and stop codons; denominator counts only synonymous codons | Yes | Peden §1.8.2.1.3 |
| INV-4 | `0 ≤ GC1, GC2, GC3, GC3s ≤ 100` | Yes | counts/positions ratio ×100 |
| INV-5 | `TotalCodons` = number of valid (ACGT-only) codons in frame | Yes | EMBOSS cusp "Number" |
| INV-6 | `OverallGc` = (GC1+GC2+GC3)/3 | Yes | record definition |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | CAI all-optimal E. coli | `CTGATCGTTGCTCGTAAA` vs `EColiOptimalCodons` | 1.0 (Within 1e-10) | Sharp & Li 1987 |
| M2 | CAI geometric-mean derivation | `GCTGCC` vs E. coli (GCT w=1, GCC w=0.122) | 0.34928498393146 | Sharp & Li 1987 √(1×0.122) |
| M3 | CAI suboptimal derivation | `CTAATAGTC` vs E. coli (0.007,0.003,0.066) | 0.01114947479545 | Sharp & Li 1987 ∛(·) |
| M4 | CAI excludes Met/Trp/stop | `ATGTGGTAA` vs E. coli | 0.0 (no scorable codon) | seqinr / CodonW |
| M5 | GC3s excludes Met/Trp/stop | `ATGGCA`: GC3s vs GC3 | GC3s=0.0, GC3=50.0 | Peden §1.8.2.1.3 |
| M6 | GC3s synonymous fraction | `GCCGCA` (2 Ala; 3rd C,A) | GC3s=50.0 | Peden §1.8.2.1.3 |
| M7 | GC1/GC2/GC3 per position | `CTGGTTAAA` | 66.6667 / 0.0 / 33.3333 | EMBOSS cusp |
| M8 | TotalCodons & counts | `CTGCTGGTTAAA` → 4 codons, CTG=2 | TotalCodons=4, counts["CTG"]=2 | EMBOSS cusp |
| M9 | E. coli table values | `EColiOptimalCodons` | CTG=1.0, GCC=0.122, CGT=1.0, AGG=0.002, TTT=0.296 | Biopython SharpEcoliIndex |
| M10 | Human table values | `HumanOptimalCodons` | CTG=2.3713, GCC=1.5988, ATG=1.0, GTG=1.8517 | Kazusa-derived RSCU |
| M11 | RSCU populated in stats | `GetStatistics("CTGCTG")` Rscu["CTG"] | =6.0 (only Leu codon present, 6-fold) | Sharp et al. 1986 RSCU |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Empty string stats | `GetStatistics("")` | all-zero `CodonUsageStatistics` | input contract |
| S2 | Empty string CAI | `CalculateCai("", EColiOptimalCodons)` | 0.0 | input contract |
| S3 | Null DnaSequence stats | `GetStatistics((DnaSequence)null)` | `ArgumentNullException` | input contract |
| S4 | Null DnaSequence CAI | `CalculateCai((DnaSequence)null, ref)` | `ArgumentNullException` | input contract |
| S5 | Null reference CAI | `CalculateCai(seq, null)` | `ArgumentNullException` | input contract |
| S6 | CAI in [0,1] | arbitrary CDS vs E. coli | 0 ≤ CAI ≤ 1 | INV-1 |
| S7 | Trailing partial codon ignored | `GCCGCAG` (7 nt) | only 2 codons counted | in-frame parsing |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | OverallGc average | `GetStatistics("CTGGTTAAA").OverallGc` | (66.6667+0+33.3333)/3 = 33.3333 | INV-6 |
| C2 | Lowercase normalized | `GetStatistics("ctg")` | same as "CTG" | case-insensitive |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzerTests.cs` — legacy fixture containing weak/permissive tests for `GetStatistics`, `CalculateCai`, `EColiOptimalCodons`, `HumanOptimalCodons` (regions "CAI", "GetStatistics", "Optimal Codon Tables", "Edge Cases"). Also contains tests for `CountCodons`, `CalculateRscu`, `CalculateEnc` belonging to **other** units (CODON-RSCU-001, CODON-ENC-001).
- `tests/Seqeron/Seqeron.Genomics.Tests/Properties/CodonProperties.cs` — property `INV-CAI3` calls `CalculateCai` (bounds in [0,1]); kept (other unit's property file).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| CAI all-optimal (M1) | ⚠ Weak | legacy `CalculateCai_OptimalCodons_ReturnsHigh` used `Within(0.01)` and a non-derived sequence |
| CAI geo-mean (M2) | ❌ Missing | no exact-derivation test existed |
| CAI suboptimal (M3) | ⚠ Weak | legacy used `Is.LessThan(0.5)` (permissive) |
| CAI excludes Met/Trp/stop (M4) | ❌ Missing | not covered (was a defect) |
| GC3s excl. Met/Trp/stop (M5) | ❌ Missing | not covered (was a defect: GC3s==GC3) |
| GC3s synonymous fraction (M6) | ❌ Missing | not covered |
| GC1/GC2/GC3 (M7) | ⚠ Weak | legacy used `Is.GreaterThan(0)` |
| TotalCodons/counts (M8) | ⚠ Weak | legacy checked `Not.Empty`, `EqualTo(5)` on undocumented seq |
| E. coli table (M9) | ⚠ Weak | legacy checked `Count==64` and ordering only |
| Human table (M10) | ⚠ Weak | legacy checked `Count==64` only |
| RSCU in stats (M11) | ❌ Missing | not covered |
| Empty stats (S1) | ⚠ Weak | legacy checked TotalCodons/Enc only |
| Empty CAI (S2) | ✅→re-do | legacy `EqualTo(0)`; re-implemented in canonical file |
| Null stats/CAI/ref (S3–S5) | ✅→re-do | legacy present; re-implemented in canonical file |
| CAI in [0,1] (S6) | ⚠ Weak | legacy range test |
| trailing partial codon (S7) | ❌ Missing | not covered |
| OverallGc (C1) | ⚠ Weak | legacy `Within(0.01)` ok but moved to canonical |
| lowercase (C2) | ❌ Missing | not covered |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/CodonUsageAnalyzer_GetStatistics_Tests.cs` — all CODON-STATS-001 cases (GetStatistics, CalculateCai, EColiOptimalCodons, HumanOptimalCodons).
- **Remove:** the `GetStatistics`/`CAI`/`Optimal Codon Tables`/relevant Edge-Case tests from legacy `CodonUsageAnalyzerTests.cs` (superseded). Keep the legacy `CountCodons`/`CalculateRscu`/`CalculateEnc` tests (other units).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `CodonUsageAnalyzer_GetStatistics_Tests.cs` | canonical CODON-STATS-001 | 21 |
| `CodonUsageAnalyzerTests.cs` | legacy; other-unit tests only | (reduced) |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ⚠ Weak | rewrote with derived seq + Within(1e-10) | ✅ Done |
| 2 | M2 | ❌ Missing | implemented √0.122 derivation | ✅ Done |
| 3 | M3 | ⚠ Weak | rewrote with ∛(·) exact value | ✅ Done |
| 4 | M4 | ❌ Missing | implemented exclusion test | ✅ Done |
| 5 | M5 | ❌ Missing | implemented GC3s≠GC3 test | ✅ Done |
| 6 | M6 | ❌ Missing | implemented GC3s synonymous fraction | ✅ Done |
| 7 | M7 | ⚠ Weak | rewrote with exact per-position values | ✅ Done |
| 8 | M8 | ⚠ Weak | rewrote with derived counts | ✅ Done |
| 9 | M9 | ⚠ Weak | rewrote with exact Sharp&Li w values | ✅ Done |
| 10 | M10 | ⚠ Weak | rewrote with exact Kazusa RSCU values | ✅ Done |
| 11 | M11 | ❌ Missing | implemented RSCU-in-stats test | ✅ Done |
| 12 | S1 | ⚠ Weak | rewrote full zeroed-struct check | ✅ Done |
| 13 | S2 | ✅ | re-implemented in canonical file | ✅ Done |
| 14 | S3 | ✅ | re-implemented in canonical file | ✅ Done |
| 15 | S4 | ✅ | re-implemented in canonical file | ✅ Done |
| 16 | S5 | ✅ | re-implemented in canonical file | ✅ Done |
| 17 | S6 | ⚠ Weak | kept as bounds property test | ✅ Done |
| 18 | S7 | ❌ Missing | implemented partial-codon test | ✅ Done |
| 19 | C1 | ⚠ Weak | rewrote with exact value | ✅ Done |
| 20 | C2 | ❌ Missing | implemented lowercase test | ✅ Done |

**Total items:** 20
**✅ Done:** 20 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact 1.0 Within(1e-10) |
| M2 | ✅ | √0.122 exact |
| M3 | ✅ | ∛(·) exact |
| M4 | ✅ | exclusion → 0.0 |
| M5 | ✅ | GC3s=0, GC3=50 |
| M6 | ✅ | GC3s=50 |
| M7 | ✅ | exact per-position |
| M8 | ✅ | counts/total exact |
| M9 | ✅ | Sharp&Li w exact |
| M10 | ✅ | Kazusa RSCU exact |
| M11 | ✅ | RSCU["CTG"]=6.0 |
| S1 | ✅ | full zeroed struct |
| S2 | ✅ | 0.0 |
| S3 | ✅ | throws |
| S4 | ✅ | throws |
| S5 | ✅ | throws |
| S6 | ✅ | bounds |
| S7 | ✅ | partial codon ignored |
| C1 | ✅ | exact average |
| C2 | ✅ | lowercase normalized |

All in-scope cases ✅. Count of ✅ = 20 = total in-scope cases.

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | GC3s reported as a percentage (×100) for consistency with GC1/GC2/GC3 (CodonW uses a fraction). | M5, M6 expected values |
| 2 | Zero-w codons are skipped (not floored to 0.01 per Bulmer 1988). | CAI of codons absent from reference |

---

## 7. Open Questions / Decisions

1. **Decision:** `EColiOptimalCodons` replaced with Sharp & Li 1987 w-values (Biopython SharpEcoliIndex) and `HumanOptimalCodons` with Kazusa-derived RSCU; prior values were untraceable (defect). CAI rescales by family max, so passing w (max 1.0) reproduces w.
2. **Decision:** Fixed two implementation defects — GC3s now restricted to synonymous third positions (Peden), and CAI now excludes Met/Trp/stop (seqinr/CodonW).
3. None open.
