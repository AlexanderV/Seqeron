# Test Specification: PROTMOTIF-COMMON-001

**Test Unit ID:** PROTMOTIF-COMMON-001
**Area:** ProteinMotif
**Algorithm:** Common Motif Finding (`ProteinMotifFinder.FindCommonMotifs`)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | PROSITE PS00001 (ASN_GLYCOSYLATION) | 2 | https://prosite.expasy.org/PS00001 | 2026-06-14 |
| 2 | PROSITE PS00005 (PKC_PHOSPHO_SITE) | 2 | https://prosite.expasy.org/PS00005 | 2026-06-14 |
| 3 | PROSITE PS00006 (CK2_PHOSPHO_SITE) | 2 | https://prosite.expasy.org/PS00006 | 2026-06-14 |
| 4 | PROSITE PS00016 (RGD) | 2 | https://prosite.expasy.org/PS00016 | 2026-06-14 |
| 5 | PROSITE PS00017 (ATP_GTP_A / P-loop) | 2 | https://prosite.expasy.org/PS00017 | 2026-06-14 |
| 6 | ScanProsite documentation (syntax, overlap, coordinates) | 2 | https://prosite.expasy.org/scanprosite/scanprosite_doc.html | 2026-06-14 |
| 7 | Sigrist et al. 2013, New and continuing developments at PROSITE | 1 | https://doi.org/10.1093/nar/gks1067 | 2026-06-14 |

### 1.2 Key Evidence Points

1. `FindCommonMotifs` scans the curated `CommonMotifs` dictionary (PROSITE-style patterns) over the sequence and returns every occurrence of every pattern — PROSITE entries (sources 1–5).
2. PROSITE pattern syntax: `[..]` allowed set, `{..}` excluded set, `x` wildcard, `x(n)` fixed gap, `x(n,m)` variable gap — ScanProsite doc (source 6).
3. PS00001 `N-{P}-[ST]-{P}`: Pro is forbidden at positions 2 and 4 — source 1.
4. ScanProsite default reporting is "greedy, overlaps, no includes": overlapping occurrences are both reported unless one is fully contained in another — source 6.
5. Matched substring content and relative positions are spec-defined; coordinate origin (1-based in PROSITE vs 0-based in repo `MotifMatch`) is a presentation convention — source 6.

### 1.3 Documented Corner Cases

- Proline exclusion in `{P}` positions rejects an otherwise-matching N-glycosylation window (source 1).
- Overlapping occurrences are reported (source 6).

### 1.4 Known Failure Modes / Pitfalls

1. Confusing 1-based PROSITE coordinates with the repository's 0-based `MotifMatch.Start`/`End` — source 6 vs repository convention.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindCommonMotifs(string proteinSequence)` | ProteinMotifFinder | Canonical | Scans all `CommonMotifs` PROSITE-style patterns; aggregates matches |
| `FindMotifByPattern(...)` | ProteinMotifFinder | Internal | Per-pattern engine; canonical under PROTMOTIF-PATTERN-001 — not re-tested here |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every returned `MotifMatch.Sequence` equals `protein.Substring(Start, End-Start+1)` | Yes | ScanProsite coordinate semantics (source 6) |
| INV-2 | `0 ≤ Start ≤ End < protein.Length` for every match | Yes | Source 6 (in-sequence hits) |
| INV-3 | Each match's `MotifName`/`Pattern` equals the `Name`/`Accession` of a `CommonMotifs` entry | Yes | `FindCommonMotifs` aggregation contract |
| INV-4 | `FindCommonMotifs` is deterministic (identical ordered output for identical input) | Yes | Regex scan is deterministic |
| INV-5 | Null or empty input yields an empty result | Yes | Trivial guard |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | NGlyco_FindsSite | `AAAANFTAAAA` contains N-glycosylation window | 1 ASN_GLYCOSYLATION match, Start=4, End=7, Seq="NFTA" | PS00001 |
| M2 | NGlyco_ExcludesProline | `AAAANPSAAAAANPTAAA` has Pro at excluded pos | 0 ASN_GLYCOSYLATION matches | PS00001 {P} |
| M3 | PKC_FindsSite | `AAAAASARKAAA` contains `[ST]-x-[RK]` | 1 PKC_PHOSPHO_SITE match, Start=5, End=7, Seq="SAR" | PS00005 |
| M4 | CK2_FindsTwoSites | `AAAASAAEASDEDAAA` contains two `[ST]-x(2)-[DE]` | 2 CK2_PHOSPHO_SITE matches: (4,7,"SAAE"),(9,12,"SDED") | PS00006 |
| M5 | PLoop_FindsSite | `AAAAAGXXXXGKSAAAA` contains `[AG]-x(4)-G-K-[ST]` | 1 ATP_GTP_A match, Start=5, End=12, Seq="GXXXXGKS" | PS00017 |
| M6 | RGD_FindsSite | `AARGDKK` contains `R-G-D` | 1 RGD match, Start=2, End=4, Seq="RGD" | PS00016 |
| M7 | MultiplePatternTypes | `RGDNFTA` triggers two different patterns (RGD + ASN_GLYCOSYLATION) | Result contains both RGD and ASN_GLYCOSYLATION motif names | whole-dictionary scan (sources 1,4) |
| M8 | MultipleOccurrences_OnePattern | `RGDRGD` contains two RGD sites | 2 RGD matches: (0,2,"RGD"),(3,5,"RGD") | PS00016 + overlap reporting (source 6) |
| M9 | Null_ReturnsEmpty | null input | empty sequence | trivial guard |
| M10 | Empty_ReturnsEmpty | "" input | empty sequence | trivial guard |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | SubstringInvariant | every match Sequence == substring at (Start..End) | INV-1 holds for all matches | content/position consistency |
| S2 | MatchIdentityFromDictionary | RGD match carries Name="RGD", Pattern="PS00016" | exact name + accession | INV-3 |
| S3 | NoMatch_ReturnsEmpty | `AAAAAAAA` (no motif) | empty result | negative control |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Determinism | two calls on same input | identical ordered results | INV-4 |
| C2 | CaseInsensitive | lowercase `aargdkk` finds RGD | 1 RGD match | implementation upcases input |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_MotifSearch_Tests.cs` (PROTMOTIF-FIND-001) already exercises `FindCommonMotifs` per-pattern (N-glyco, PKC, CK2, P-loop) and null input. Those tests belong to FIND-001's canonical file and are not duplicated here.
- This unit's canonical file `ProteinMotifFinder_FindCommonMotifs_Tests.cs` focuses on the aggregation semantics specific to `FindCommonMotifs` (whole-dictionary scan, multi-pattern, multi-occurrence, identity propagation, invariants) plus its own edge cases.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ❌ Missing | new canonical file |
| M2 | ❌ Missing | new canonical file |
| M3 | ❌ Missing | new canonical file |
| M4 | ❌ Missing | new canonical file |
| M5 | ❌ Missing | new canonical file |
| M6 | ❌ Missing | new canonical file |
| M7 | ❌ Missing | new canonical file |
| M8 | ❌ Missing | new canonical file |
| M9 | ❌ Missing | new canonical file |
| M10 | ❌ Missing | new canonical file |
| S1 | ❌ Missing | new canonical file |
| S2 | ❌ Missing | new canonical file |
| S3 | ❌ Missing | new canonical file |
| C1 | ❌ Missing | new canonical file |
| C2 | ❌ Missing | new canonical file |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindCommonMotifs_Tests.cs` — all PROTMOTIF-COMMON-001 cases.
- **Remove:** nothing. Existing FIND-001 tests stay under their own unit; no duplication of the same assertion here.

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `ProteinMotifFinder_FindCommonMotifs_Tests.cs` | PROTMOTIF-COMMON-001 canonical | 15 |

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
| 8 | M8 | ❌ Missing | Implemented | ✅ Done |
| 9 | M9 | ❌ Missing | Implemented | ✅ Done |
| 10 | M10 | ❌ Missing | Implemented | ✅ Done |
| 11 | S1 | ❌ Missing | Implemented | ✅ Done |
| 12 | S2 | ❌ Missing | Implemented | ✅ Done |
| 13 | S3 | ❌ Missing | Implemented | ✅ Done |
| 14 | C1 | ❌ Missing | Implemented | ✅ Done |
| 15 | C2 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 15
**✅ Done:** 15 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | `FindCommonMotifs_NGlycosylationWindow_FindsExactSite` |
| M2 | ✅ Covered | `FindCommonMotifs_ProlineAtExcludedPosition_RejectsSite` |
| M3 | ✅ Covered | `FindCommonMotifs_PkcWindow_FindsExactSite` |
| M4 | ✅ Covered | `FindCommonMotifs_Ck2TwoWindows_FindsBothSites` |
| M5 | ✅ Covered | `FindCommonMotifs_PLoopWindow_FindsExactSite` |
| M6 | ✅ Covered | `FindCommonMotifs_RgdWindow_FindsExactSite` |
| M7 | ✅ Covered | `FindCommonMotifs_TwoDistinctPatterns_ReturnsBothTypes` |
| M8 | ✅ Covered | `FindCommonMotifs_TwoRgdOccurrences_ReturnsBoth` |
| M9 | ✅ Covered | `FindCommonMotifs_NullSequence_ReturnsEmpty` |
| M10 | ✅ Covered | `FindCommonMotifs_EmptySequence_ReturnsEmpty` |
| S1 | ✅ Covered | `FindCommonMotifs_AllMatches_SatisfySubstringInvariant` |
| S2 | ✅ Covered | `FindCommonMotifs_RgdMatch_CarriesDictionaryIdentity` |
| S3 | ✅ Covered | `FindCommonMotifs_SequenceWithNoMotif_ReturnsEmpty` |
| C1 | ✅ Covered | `FindCommonMotifs_SameInput_IsDeterministic` |
| C2 | ✅ Covered | `FindCommonMotifs_LowercaseInput_FindsSite` |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | 0-based inclusive `MotifMatch` coordinates (PROSITE is 1-based) — API-shape convention, not correctness-affecting | M1, M3, M4, M5, M6, M8 |

---

## 7. Open Questions / Decisions

1. `FindAllKnownMotifs` is listed in the Registry method table but does not exist in the codebase and has no authoritative basis distinct from `FindCommonMotifs`. The evidence-derived canonical method for this unit is `FindCommonMotifs`; no aliasing method is invented. Recorded in the algorithm doc and Registry update.
