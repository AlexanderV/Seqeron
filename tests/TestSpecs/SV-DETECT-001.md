# Test Specification: SV-DETECT-001

**Test Unit ID:** SV-DETECT-001
**Area:** StructuralVar
**Algorithm:** Structural Variant Detection from Paired-End Mapping signatures
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Medvedev, Stanciu & Brudno (2009), Nat Methods 6(11s):S13–S20 | 1 | https://doi.org/10.1038/nmeth.1374 | 2026-06-13 |
| 2 | Chen et al. (2009) BreakDancer README (genome/breakdancer) | 3 | https://raw.githubusercontent.com/genome/breakdancer/master/README | 2026-06-13 |
| 3 | Fan, Chen et al. (2014) BreakDancer protocol, PMC3661775 | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC3661775/ | 2026-06-13 |
| 4 | Forward/reverse reads, FR proper pair (SAM FLAG 0x02 / BWA) | 4 | https://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/ | 2026-06-13 |

### 1.2 Key Evidence Points

1. Deletion ⇒ mapped span GREATER than insert size — Medvedev et al. 2009.
2. Insertion ⇒ mapped span SMALLER than insert size — Medvedev et al. 2009.
3. Inversion ⇒ one mate's orientation flipped (intra-chromosomal FF/RR) — Medvedev et al. 2009; cureffi/BWA.
4. Translocation ⇒ mates on different chromosomes (linking/CTX signature) — Medvedev et al. 2009; BreakDancer.
5. Discordant-by-span cutoff: bounds = mean ± c·sd, default c = 3 — BreakDancer README; protocol confirms "3 s.d.".
6. Concordant orientation = FR (one mate +, one mate −); FF/RR/RF is abnormal — cureffi/BWA, SAM proper-pair FLAG 0x02.
7. Minimum supporting read pairs to call an SV: default 2 — BreakDancer README (-r).

### 1.3 Documented Corner Cases

- Insertions larger than the fragment insert size produce no span signature and the inserted sequence is not recovered (Medvedev et al. 2009).
- Clusters below the minimum support are not reported (BreakDancer -r).

### 1.4 Known Failure Modes / Pitfalls

1. Treating a same-orientation inter-chromosomal pair as an inversion — it is a translocation; chromosome difference takes precedence (Medvedev et al. 2009; resolved via ASSUMPTION A1).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `ClassifySV(ReadPairSignature, expectedInsertSize, insertSizeStdDev, cutoff)` | StructuralVariantAnalyzer | Canonical | Maps a discordant signature to an SVType per PEM signature rules. |
| `DetectSVs(readPairs, expectedInsertSize, insertSizeStdDev, cutoff, clusterDistance, minSupport)` | StructuralVariantAnalyzer | Canonical | Orchestrates find-discordant → cluster → classify. |
| `FindDiscordantPairs(readPairs, expectedInsertSize, insertSizeStdDev, cutoff, maxInsertSize)` | StructuralVariantAnalyzer | Internal | Anomaly detection; tested indirectly via DetectSVs and directly for cutoff boundary. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | A pair with chr1 ≠ chr2 classifies as Translocation regardless of orientation/span. | Yes | Medvedev et al. 2009 (linking across chromosomes); ASSUMPTION A1 (precedence) |
| INV-2 | An intra-chromosomal pair with strand1 == strand2 classifies as Inversion. | Yes | Medvedev et al. 2009; cureffi/BWA |
| INV-3 | An intra-chromosomal FR/RF pair with span > mean + c·sd classifies as Deletion. | Yes | Medvedev et al. 2009 |
| INV-4 | An intra-chromosomal FR/RF pair with span < mean − c·sd classifies as Insertion. | Yes | Medvedev et al. 2009 |
| INV-5 | A pair is flagged discordant-by-span iff span < mean − c·sd OR span > mean + c·sd. | Yes | BreakDancer README (bounds = mean ± c·sd) |
| INV-6 | DetectSVs emits an SV for a cluster iff its supporting-pair count ≥ minSupport. | Yes | BreakDancer -r default 2 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | ClassifySV deletion | Same chr, FR, span 5000 (>550) | SVType.Deletion | Medvedev et al. 2009 (span larger) |
| M2 | ClassifySV insertion | Same chr, FR, span 100 (<250) | SVType.Insertion | Medvedev et al. 2009 (span smaller) |
| M3 | ClassifySV inversion | Same chr, FF orientation, span 400 | SVType.Inversion | Medvedev et al. 2009; cureffi/BWA |
| M4 | ClassifySV translocation | chr1="chr1", chr2="chr2" | SVType.Translocation | Medvedev et al. 2009; BreakDancer CTX |
| M5 | ClassifySV inter-chr precedence | chr1≠chr2 AND FF orientation | SVType.Translocation (not Inversion) | Medvedev et al. 2009; ASSUMPTION A1 |
| M6 | FindDiscordantPairs concordant excluded | Same chr, FR, span 400 (in [250,550]) | not returned (concordant) | BreakDancer normal class; FR proper pair |
| M7 | DetectSVs min-support gate (below) | 1 deletion-signature pair, minSupport 2 | no SV emitted | BreakDancer -r default 2 |
| M8 | DetectSVs min-support gate (meets) | 3 clustered deletion pairs, minSupport 2 | 1 SV, Type=Deletion, SupportingReads=3 | BreakDancer -r; Medvedev clustering |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Cutoff lower boundary | span = mean − c·sd (=250) exactly | concordant (not discordant) | bound is inclusive: discordant iff strictly outside |
| S2 | Cutoff just below lower bound | span = 249 | discordant, Insertion | one unit beyond lower bound |
| S3 | Cutoff upper boundary | span = mean + c·sd (=550) exactly | concordant | inclusive upper bound |
| S4 | RF orientation within bounds | same chr, RF, span 400 | concordant (FR/RF both proper) | cureffi: FR/RF pointing-inward proper |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | DetectSVs empty input | no read pairs | empty result | trivial defined behavior |
| C2 | DetectSVs null input | null readPairs | ArgumentNullException | input validation |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing test file for StructuralVariantAnalyzer (`grep` of `tests/Seqeron/Seqeron.Genomics.Tests/` found none). `StructuralVariantAnalyzer.cs` exists with `FindDiscordantPairs`, `ClusterDiscordantPairs`, etc., but no `DetectSVs`/`ClassifySV` (the Registry canonical methods) and no tests.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 | ✅ Covered | ClassifySV_LargeSpanSameChr_ReturnsDeletion |
| M2 | ✅ Covered | ClassifySV_SmallSpanSameChr_ReturnsInsertion |
| M3 | ✅ Covered | ClassifySV_SameOrientationSameChr_ReturnsInversion |
| M4 | ✅ Covered | ClassifySV_DifferentChromosomes_ReturnsTranslocation |
| M5 | ✅ Covered | ClassifySV_DifferentChromosomesSameOrientation_ReturnsTranslocation |
| M6 | ✅ Covered | FindDiscordantPairs_ConcordantFrPair_NotReturned |
| M7 | ✅ Covered | DetectSVs_BelowMinSupport_EmitsNoSv |
| M8 | ✅ Covered | DetectSVs_MeetsMinSupport_EmitsOneDeletion |
| S1 | ✅ Covered | FindDiscordantPairs_SpanAtLowerBound_NotDiscordant |
| S2 | ✅ Covered | FindDiscordantPairs_SpanBelowLowerBound_IsDiscordant |
| S3 | ✅ Covered | FindDiscordantPairs_SpanAtUpperBound_NotDiscordant |
| S4 | ✅ Covered | FindDiscordantPairs_RfOrientationWithinBounds_NotDiscordant |
| C1 | ✅ Covered | DetectSVs_EmptyInput_ReturnsEmpty |
| C2 | ✅ Covered | DetectSVs_NullInput_Throws |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/StructuralVariantAnalyzer_DetectSVs_Tests.cs` — all SV-DETECT-001 cases (ClassifySV, DetectSVs, FindDiscordantPairs boundary).
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `StructuralVariantAnalyzer_DetectSVs_Tests.cs` | Canonical SV-DETECT-001 fixture | 14 |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 2 | M2 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 3 | M3 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 4 | M4 | ❌ Missing | Implemented exact-value test | ✅ Done |
| 5 | M5 | ❌ Missing | Implemented precedence test | ✅ Done |
| 6 | M6 | ❌ Missing | Implemented concordant-exclusion test | ✅ Done |
| 7 | M7 | ❌ Missing | Implemented min-support gate test | ✅ Done |
| 8 | M8 | ❌ Missing | Implemented min-support meets test | ✅ Done |
| 9 | S1 | ❌ Missing | Implemented lower-boundary test | ✅ Done |
| 10 | S2 | ❌ Missing | Implemented below-lower-bound test | ✅ Done |
| 11 | S3 | ❌ Missing | Implemented upper-boundary test | ✅ Done |
| 12 | S4 | ❌ Missing | Implemented RF-within-bounds test | ✅ Done |
| 13 | C1 | ❌ Missing | Implemented empty-input test | ✅ Done |
| 14 | C2 | ❌ Missing | Implemented null-input test | ✅ Done |

**Total items:** 14
**✅ Done:** 14 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | Exact SVType.Deletion asserted |
| M2 | ✅ | Exact SVType.Insertion asserted |
| M3 | ✅ | Exact SVType.Inversion asserted |
| M4 | ✅ | Exact SVType.Translocation asserted |
| M5 | ✅ | Translocation precedence asserted |
| M6 | ✅ | Concordant pair excluded |
| M7 | ✅ | No SV below min support |
| M8 | ✅ | One Deletion, SupportingReads=3 |
| S1 | ✅ | Lower bound inclusive |
| S2 | ✅ | Below bound discordant |
| S3 | ✅ | Upper bound inclusive |
| S4 | ✅ | RF concordant |
| C1 | ✅ | Empty → empty |
| C2 | ✅ | Null → ArgumentNullException |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| A1 | Inter-chromosomal pairs classify as Translocation before orientation/span are considered (chromosome difference has precedence over the inversion rule). | INV-1, M5, ClassifySV |

---

## 7. Open Questions / Decisions

1. None. The s.d. cutoff (default 3), min support (default 2), concordant FR orientation, and the four PEM signatures are all source-traceable; the single precedence ordering is documented as ASSUMPTION A1 with a biological rationale (a cross-chromosome event is a translocation by definition, not an inversion).
