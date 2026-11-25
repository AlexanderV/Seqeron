# Test Specification: ASSEMBLY-MERGE-001

**Test Unit ID:** ASSEMBLY-MERGE-001
**Area:** Assembly
**Algorithm:** Contig Merging (suffix–prefix overlap collapse)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Langmead, *Assembly & shortest common superstring* (JHU notes) | 1 | https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_scs.pdf | 2026-06-13 |
| 2 | Langmead, *Overlap Layout Consensus assembly* (JHU notes) | 1 | https://www.cs.jhu.edu/~langmea/resources/lecture_notes/assembly_olc.pdf | 2026-06-13 |
| 3 | MIT 7.91J Lecture 6, *Genome Assembly* (MIT OCW, 2014) | 1 | https://ocw.mit.edu/courses/7-91j-foundations-of-computational-and-systems-biology-spring-2014/e885f0eb376ea6c2045eb9d8847f106f_MIT7_91JS14_Lecture6.pdf | 2026-06-13 |
| 4 | Compeau, Pevzner & Tesler (2011), Nat Biotechnol 29:987–991 | 1 | https://doi.org/10.1038/nbt.2023 | 2026-06-13 |
| 5 | Langmead, `CG_SCS.ipynb` (`overlap`, `scs`: `sup += y[olen:]`) | 1 (reference code) | https://raw.githubusercontent.com/BenLangmead/comp-genomics-class/master/notebooks/CG_SCS.ipynb | 2026-11-25 |
| 6 | EMBOSS `merger` source, emboss/merger.c | 1 (reference implementation) | https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/merger.c | 2026-11-25 |

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
| `MergeContigs(string, string, int, bool resolveOverlapConflicts)` | SequenceAssembler | Canonical | Opt-in EMBOSS `merger` overlap conflict rule; MCP `merge_contigs` delegates here |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-01 | With a valid overlap `0 ≤ l ≤ min(|c1|,|c2|)`, result = `c1 + c2[l..]` (default); length = `|c1| + |c2| − l` (both modes) | Yes | Sources 1, 5 |
| INV-02 | Overlap length 0 ⇒ result = `c1 + c2` (plain concatenation) | Yes | Source 1 ("just concatenate them") |
| INV-03 | Overlap length `l < 0` or `l > min(|c1|,|c2|)` ⇒ `ArgumentOutOfRangeException` | Yes | Sources 1, 5 (overlap ∈ [0, min length]); F-MERGE-1 |
| INV-05 | `resolveOverlapConflicts`: overlap columns resolved by the EMBOSS `merger` rule (N yields; quality windows 5/10/20; tie → contig1) | Yes | Source 6 |
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
| S1 | Oversized overlap | overlap > `min(|c1|,|c2|)` (`AC`+`GTAA`,3; `ACGT`+`ACG`,4), both modes | `ArgumentOutOfRangeException` | INV-03 (Sources 1, 5) |
| S2 | Negative overlap | overlap = −2 | `ArgumentOutOfRangeException` | INV-03 |
| S5 | Overlap = both lengths | `ACGT`+`ACGT`,4 | `"ACGT"` | INV-01 |
| S6 | Default mode, mismatched overlap | `AAAACCCCGGGG`+`CCCCGTGGTTTT`,8 | `"AAAACCCCGGGGTTTT"` (contig1 kept) | Source 5 |
| S7 | EMBOSS conflict rule (11 rows: w5/w10/w20 winners each side, full tie, N in contig1/contig2/both, exact overlap) | `MergeContigs(..., resolveOverlapConflicts: true)` | EMBOSS-port values (Evidence table) | Source 6 |
| S8 | EMBOSS mode, case-only difference | `acgtAC`+`GTAcgt`,4 | `"acgtACgt"` (no conflict) | Source 6 (lower-cases inputs) |
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
| S5–S8 | ❌ Missing | Added 2026-11-25 (B10 review) |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/SequenceAssembler_MergeContigs_Tests.cs`
  — all MUST/SHOULD/COULD cases for `MergeContigs`.
- **Remove:** nothing (no pre-existing `MergeContigs` tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `SequenceAssembler_MergeContigs_Tests.cs` | Canonical | 27 |

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
| S1/S2 | ✅ Covered (changed) | Now assert `ArgumentOutOfRangeException` (F-MERGE-1) |
| S5–S8 | ✅ Covered | `MergeContigs_OverlapEqualsBothLengths_ReturnsSingleCopy`, `MergeContigs_DefaultMode_MismatchedOverlapKeepsContig1Bases`, `MergeContigs_Resolve_*` (11), `MergeContigs_ResolveOverlapConflicts_CaseOnlyDifference_IsNotAConflict` |

**✅ Covered: 12 / 12 in-scope cases.**

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Caller-supplied overlap length is trusted (verification is `FindOverlap`'s responsibility) — API contract, not a numeric parameter | §1.4, contract |
| 2 | ~~Out-of-range overlap → plain concatenation~~ withdrawn 2026-11-25: now rejected (F-MERGE-1) | INV-03, S1, S2 |

---

## 7. Open Questions / Decisions

1. None. The two assumptions are API-contract boundaries derived from the authoritative overlap
   definition, not invented numeric/scoring values; no correctness-affecting assumption remains.

## R1-WP1 addendum (2026-11-25): ASCII fold / span k-mers (A15-7)

`SequenceAssembler` private `char.ToUpperInvariant` folds -> `SequenceCase.ToUpperAscii`; `BuildDeBruijnGraph` / `SpellStrandPairsOnce` use `KmerEnumerator` + span dictionary lookups; `ReverseComplement` checks the `TryGetReverseComplement` result. New tests in `Unit/Alignment/SequenceAssembler_AsciiFold_Tests.cs`:

| Test | Locks |
|---|---|
| `ComputeConsensus_LongS_IsNotCountedAsAsciiS` | Biopython `dumb_consensus`: S, U+017F, s -> 2/3 < 0.7 -> N; threshold 0.5 -> S |
| `MergeContigs_ResolveConflicts_LongSDiffersFromS_BestQualityDecides` | EMBOSS merger quality windows: "A\u017F"+"ST" -> "AST" |
| `Scaffold_NegativeGap_LongSOverlapDisagrees_WritesUnknownGap` / `..._LowerCaseOverlapAgrees_Collapsed` | ABySS agreement ignoring ASCII case only; AGP unknown gap 100 |
| `BuildDeBruijnGraph_NonAsciiLetter_KeptVerbatim_SharedPrefixNodeAccumulates` | node/edge multimap unchanged by the span path |
| `AssembleDeBruijn_BothStrands_NonAsciiLetter_StrandPairEmittedOnce` | strand-pair emission unchanged |

