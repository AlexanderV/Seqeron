# Test Specification: ASSEMBLY-MERGE-001

**Test Unit ID:** ASSEMBLY-MERGE-001
**Area:** Assembly
**Algorithm:** Contig Merging (suffix–prefix overlap collapse)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Langmead, *Assembly & shortest common superstring* (JHU notes) | 1 | https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_scs.pdf | 2026-06-13 |
| 2 | Langmead, *Overlap Layout Consensus assembly* (JHU notes) | 1 | https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf | 2026-06-13 |
| 3 | MIT 7.91J Lecture 6, *Genome Assembly* (MIT OCW, 2014) | 1 | https://ocw.mit.edu/courses/7-91j-foundations-of-computational-and-systems-biology-spring-2014/e885f0eb376ea6c2045eb9d8847f106f_MIT7_91JS14_Lecture6.pdf | 2026-06-13 |
| 4 | Compeau, Pevzner & Tesler (2011), Nat Biotechnol 29:987–991 | 1 | https://doi.org/10.1038/nbt.2023 | 2026-06-13 |

### 1.2 Key Evidence Points

1. An overlap is a length-l suffix of X that exactly matches a length-l prefix of Y — Source 1, 3.
2. Merging two overlapping strings keeps a single copy of the overlap: result = X + (Y with its
   length-l prefix removed). Worked: `BAA` + `AAB`, overlap 2 → `BAAB` (len 4 = 3 + 3 − 2) — Source 1.
3. With no overlap the strings are simply concatenated ("just concatenate them"); overlap length 0
   ⇒ result = X + Y — Source 1.
4. A valid overlap length is bounded by `min(|X|, |Y|)` (cannot take a suffix/prefix longer than the
   string); `suffixPrefixMatch` guards `if len(x) < k or len(y) < k: return 0` — Source 1.
5. When several overlaps exist, only the longest suffix/prefix match is reported/used — Source 2.
6. "Fragments are contigs"; OLC layout bundles overlapping reads into contigs — Source 2, 3.

### 1.3 Documented Corner Cases

- Overlap length 0 (no qualifying match) → plain concatenation `X + Y` (Source 1).
- Overlap bounded by the shorter string; oversized request is not a valid overlap (Source 1).

### 1.4 Known Failure Modes / Pitfalls

1. Removing more than `min(|X|,|Y|)` characters would over-trim / corrupt the merge — guarded by the
   shorter-string bound (Source 1).
2. Trusting an unverified overlap length collapses a region that may not actually match; in this
   library the length is produced by `FindOverlap` (Source 2) — documented contract boundary.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `MergeContigs(string contig1, string contig2, int overlapLength)` | SequenceAssembler | Canonical | Suffix/prefix overlap collapse primitive |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | With a valid overlap `0 < l ≤ min(|c1|,|c2|)`, result = `c1 + c2[l..]`; length = `|c1| + |c2| − l` | Yes | Source 1 (merge keeps one copy of overlap) |
| INV-02 | Overlap length 0 ⇒ result = `c1 + c2` (plain concatenation) | Yes | Source 1 ("just concatenate them") |
| INV-03 | Overlap length `l ≤ 0` or `l > min(|c1|,|c2|)` ⇒ result = `c1 + c2` (no usable overlap) | Yes | Source 1 (overlap 0 = concat; bound = min length) |
| INV-04 | Result always contains `c1` as a prefix and the non-overlapped tail of `c2` as a suffix | Yes | Source 1 (collapse preserves both strings' content) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Published two-string merge | `BAA` + `AAB`, overlap 2 | `"BAAB"` (length 4) | Source 1 greedy trace |
| M2 | Chained merge to SCS | merge `AAA`+`AAB` (2) then result+`ABB`(2)… to full chain {AAA,AAB,ABB,BBB,BBA} overlaps 2 | `"AAABBBA"` (length 7) | Source 1 SCS example |
| M3 | No overlap → concatenation | `BAA` + `AAB`, overlap 0 | `"BAAAAB"` (length 6 = `c1+c2`) | Source 1 ("concatenate") |
| M4 | Boundary overlap = min length | `ACGT` + `CGTAA`, overlap 3 (suffix `CGT` = prefix `CGT`) | `"ACGTAA"` (length 6) | Source 1 (overlap ≤ min len; collapse) |
| M5 | Length invariant for valid overlap | `GATTACA` + `ACATGAA`, overlap 3 (suffix `ACA` = prefix `ACA`) | `"GATTACATGAA"` (length 11 = 7 + 7 − 3) | INV-01, Source 1 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Oversized overlap | overlap > `min(|c1|,|c2|)` | result = `c1 + c2` | INV-03 / bound (Source 1) |
| S2 | Negative overlap | overlap = −2 | result = `c1 + c2` | INV-03 (non-positive = no overlap) |
| S3 | Null contig1 | `MergeContigs(null, "AAB", 1)` | `ArgumentNullException` | Sibling null-validation convention |
| S4 | Null contig2 | `MergeContigs("BAA", null, 1)` | `ArgumentNullException` | Sibling null-validation convention |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Empty + non-empty, overlap 0 | `""` + `AAB`, overlap 0 | `"AAB"` | Identity element |
| C2 | Non-empty + empty, overlap 0 | `BAA` + `""`, overlap 0 | `"BAA"` | Identity element |
| C3 | Prefix-containment property | result starts with `c1` and ends with `c2[l..]` | true | INV-04 |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssemblerTests.cs` — legacy fixture for the class;
  searched for `MergeContigs` references.
- `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_AssembleOLC_Tests.cs`,
  `SequenceAssembler_AssembleDeBruijn_Tests.cs` — sibling per-method fixtures (no `MergeContigs`).
- No existing test exercises `MergeContigs`.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | New unit |
| M2 | ❌ Missing | New unit |
| M3 | ❌ Missing | New unit |
| M4 | ❌ Missing | New unit |
| M5 | ❌ Missing | New unit |
| S1 | ❌ Missing | New unit |
| S2 | ❌ Missing | New unit |
| S3 | ❌ Missing | New unit |
| S4 | ❌ Missing | New unit |
| C1 | ❌ Missing | New unit |
| C2 | ❌ Missing | New unit |
| C3 | ❌ Missing | New unit |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_MergeContigs_Tests.cs`
  — all MUST/SHOULD/COULD cases for `MergeContigs`.
- **Remove:** nothing (no pre-existing `MergeContigs` tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceAssembler_MergeContigs_Tests.cs` | Canonical | 12 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented | ✅ Done |
| 6 | S1 | ❌ Missing | Implemented | ✅ Done |
| 7 | S2 | ❌ Missing | Implemented | ✅ Done |
| 8 | S3 | ❌ Missing | Implemented | ✅ Done |
| 9 | S4 | ❌ Missing | Implemented | ✅ Done |
| 10 | C1 | ❌ Missing | Implemented | ✅ Done |
| 11 | C2 | ❌ Missing | Implemented | ✅ Done |
| 12 | C3 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 12
**✅ Done:** 12 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | Exact `"BAAB"` asserted |
| M2 | ✅ Covered | Exact `"AAABBBA"` asserted |
| M3 | ✅ Covered | Exact `"BAAAAB"` asserted |
| M4 | ✅ Covered | Exact `"ACGTAA"` asserted |
| M5 | ✅ Covered | Length-invariant exact value asserted |
| S1 | ✅ Covered | Concatenation asserted |
| S2 | ✅ Covered | Concatenation asserted |
| S3 | ✅ Covered | `ArgumentNullException` asserted |
| S4 | ✅ Covered | `ArgumentNullException` asserted |
| C1 | ✅ Covered | Identity asserted |
| C2 | ✅ Covered | Identity asserted |
| C3 | ✅ Covered | Prefix/suffix containment asserted |

**✅ Covered: 12 / 12 in-scope cases.**

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Caller-supplied overlap length is trusted (verification is `FindOverlap`'s responsibility) — API contract, not a numeric parameter | §1.4, contract |
| 2 | Out-of-range overlap (≤0 or > min length) → plain concatenation, derived from Source 1 facts (overlap 0 = concat; bound = min length) | INV-03, S1, S2 |

---

## 7. Open Questions / Decisions

1. None. The two assumptions are API-contract boundaries derived from the authoritative overlap
   definition, not invented numeric/scoring values; no correctness-affecting assumption remains.
