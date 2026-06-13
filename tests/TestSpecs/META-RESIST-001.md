# Test Specification: META-RESIST-001

**Test Unit ID:** META-RESIST-001
**Area:** Metagenomics
**Algorithm:** Antibiotic Resistance Gene Detection (ResFinder-style acquired-gene detection)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Zankari et al. (2012), JAC 67(11):2640–2644 (ResFinder) | 1 | https://academic.oup.com/jac/article/67/11/2640/707208 | 2026-06-13 |
| 2 | genomicepidemiology/resfinder (reference impl. defaults) | 3 | https://github.com/genomicepidemiology/resfinder | 2026-06-13 |
| 3 | Pipeline validation, Sci Rep (2023) 13 | 1 | https://www.nature.com/articles/s41598-023-42154-6 | 2026-06-13 |
| 4 | Benchmarking AMR-gene ID methods, JAC (2016) 71(9):2484 | 1 | https://academic.oup.com/jac/article/71/9/2484/2238319 | 2026-06-13 |
| 5 | Li H (2018), On the definition of sequence identity | 3 | https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity | 2026-06-13 |
| 6 | CARD RGI (McMaster) | 5 | https://card.mcmaster.ca/analyze/rgi | 2026-06-13 |

### 1.2 Key Evidence Points

1. ResFinder BLASTs each database gene vs the query and reports the **best-matching gene** — Source 1.
2. **%identity (BLAST)** = identical positions / alignment columns; gapless ⇒ denominator = aligned length — Source 5; Source 1 defines %ID as fraction of identical nucleotides between best match and genome.
3. **Coverage** = fraction of the *reference* gene length covered; ResFinder floor ≥ 2/5 (0.60) — Sources 1, 3, 4.
4. **Default thresholds** (web service): 0.90 identity / 0.60 coverage; GitHub README defaults 0.80/0.60; study uses 0.98 — Sources 2, 3, 4.
5. Sub-threshold matches are gene **fragments/noise** and must not be reported — Source 1.
6. 60% coverage floor exists so contig-edge / split genes are still detected — Source 3.
7. Best-hit ranking corroborated by CARD RGI (best bit-score hit); 100% identity full-length = "Perfect" — Source 6.

### 1.3 Documented Corner Cases

- Sub-threshold-identity hits are gene fragments → rejected (Source 1).
- Partial hits below the coverage floor → rejected (Source 1).
- Contig-edge / fragmented genes scored against full reference length, retained if ≥ coverage floor (Source 3).

### 1.4 Known Failure Modes / Pitfalls

1. Treating %identity as motif-length / contig-length is not BLAST identity — Source 5.
2. Computing coverage against the contig instead of the reference gene length — Sources 1, 3.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindAntibioticResistanceGenes(contigs, referenceGenes, identityThreshold, coverageThreshold)` | MetagenomicsAnalyzer | **Canonical** | Best-match detector; deep evidence-based tests |
| `BestUngappedMatch(contig, reference)` | MetagenomicsAnalyzer (private) | **Internal** | Tested indirectly via the canonical method |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | 0 ≤ PercentIdentity ≤ 1 and 0 ≤ Coverage ≤ 1 | Yes | Source 5 (matches ≤ columns); coverage = window/ref ≤ 1 |
| INV-2 | A hit is reported only if PercentIdentity ≥ identityThreshold AND Coverage ≥ coverageThreshold | Yes | Source 1 |
| INV-3 | At most one hit per contig: the best-matching reference gene (max identity, tie → max coverage) | Yes | Source 1 ("best-matching gene"), Source 6 |
| INV-4 | Exact full-length match ⇒ PercentIdentity = 1.0 and Coverage = 1.0 | Yes | Source 6 (Perfect), Source 5 |
| INV-5 | Default thresholds = 0.90 identity, 0.60 coverage | Yes | Sources 2,3,4 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Exact full-length match | Contig contains reference verbatim | %ID=1.0, coverage=1.0, reported | Source 6; Source 5 |
| M2 | Single mismatch, full length | One substituted base over 7-base reference | %ID=6/7≈0.857142857, coverage=1.0 | Source 5 (matches/columns) |
| M3 | Contig-edge partial hit | Only 4 of 7 reference bases present at contig end, default cov=0.6 ⇒ 4/7≈0.571 below floor | not reported | Source 1 (2/5 rule), Source 3 |
| M4 | Partial hit above coverage floor | 5 of 7 bases, coverage 5/7≈0.714 ≥ 0.6, %ID=1.0 | reported, coverage=5/7 | Source 3 |
| M5 | Below identity threshold | match identity below threshold | not reported | Source 1 (fragments/noise) |
| M6 | Best-matching gene only | two reference genes match, one higher identity | only higher-identity gene reported | Source 1; Source 6 |
| M7 | Default thresholds value | constants equal ResFinder values | 0.90 and 0.60 | Sources 2,3,4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Null contigs | contigs=null | ArgumentNullException | contract |
| S2 | Null referenceGenes | referenceGenes=null | ArgumentNullException | contract |
| S3 | Identity threshold out of range | identityThreshold=1.5 | ArgumentOutOfRangeException | contract |
| S4 | Coverage threshold out of range | coverageThreshold=-0.1 | ArgumentOutOfRangeException | contract |
| S5 | Empty contig sequence / no matches | empty/non-matching input | empty result | contract |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Tie-break by coverage | two genes equal identity, different coverage | higher-coverage gene reported | INV-3 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Old `MetagenomicsAnalyzerTests.cs` contains `FindResistanceGenes_*` tests for the legacy motif-containment stub (`FindResistanceGenes`), which is a **different** method (out of scope, MCP-wired). No tests exist for the new canonical `FindAntibioticResistanceGenes`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new unit |
| M2 | ❌ Missing | new unit |
| M3 | ❌ Missing | new unit |
| M4 | ❌ Missing | new unit |
| M5 | ❌ Missing | new unit |
| M6 | ❌ Missing | new unit |
| M7 | ❌ Missing | new unit |
| S1 | ❌ Missing | new unit |
| S2 | ❌ Missing | new unit |
| S3 | ❌ Missing | new unit |
| S4 | ❌ Missing | new unit |
| S5 | ❌ Missing | new unit |
| C1 | ❌ Missing | new unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/MetagenomicsAnalyzer_FindAntibioticResistanceGenes_Tests.cs` — all cases for this unit.
- **Remove:** nothing. Legacy `FindResistanceGenes_*` tests target a separate method and stay.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `MetagenomicsAnalyzer_FindAntibioticResistanceGenes_Tests.cs` | Canonical unit tests | 13 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented | ✅ Done |
| 8 | S1 | ❌ Missing | Implemented | ✅ Done |
| 9 | S2 | ❌ Missing | Implemented | ✅ Done |
| 10 | S3 | ❌ Missing | Implemented | ✅ Done |
| 11 | S4 | ❌ Missing | Implemented | ✅ Done |
| 12 | S5 | ❌ Missing | Implemented | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 13
**✅ Done:** 13 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact match test |
| M2 | ✅ | Single-mismatch identity test |
| M3 | ✅ | Edge partial below floor → not reported |
| M4 | ✅ | Partial above floor reported |
| M5 | ✅ | Below identity threshold → not reported |
| M6 | ✅ | Best-match-only test |
| M7 | ✅ | Default-threshold constants test |
| S1 | ✅ | Null contigs |
| S2 | ✅ | Null referenceGenes |
| S3 | ✅ | Identity out of range |
| S4 | ✅ | Coverage out of range |
| S5 | ✅ | Empty / no match |
| C1 | ✅ | Tie-break by coverage |

**In-scope cases:** 13 | **✅:** 13

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Gapless (ungapped) alignment model (vs full gapped BLAST); identical to BLAST identity formula when no indels | Implementation of `BestUngappedMatch`; M2, M3, M4 |

---

## 7. Open Questions / Decisions

1. Gene database is caller-supplied (no hard-coded CARD/ResFinder gene list), per the unit directive that curated tables must not be fabricated. Decision accepted; Implementation Status = Framework.
