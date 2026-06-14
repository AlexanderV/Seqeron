# Test Specification: ONCO-FUSION-001

**Test Unit ID:** ONCO-FUSION-001
**Area:** Oncology
**Algorithm:** Fusion Gene Detection (candidate fusion calling from breakpoint-supporting reads)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Haas et al. 2019, Genome Biology 20:213 (fusion-detection benchmark) | 1 | https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1842-9 | 2026-06-14 |
| 2 | STAR-Fusion source (defaults MIN_JUNCTION_READS=1, MIN_SUM_FRAGS=2, MIN_SPANNING_FRAGS_ONLY=5) | 3 | https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion | 2026-06-14 |
| 3 | Uhrig et al. 2021, Genome Research 31(3):448 (Arriba) | 1 | https://pmc.ncbi.nlm.nih.gov/articles/PMC7919457/ | 2026-06-14 |
| 4 | Arriba output-file spec (split_reads1/2, discordant_mates, reading_frame) | 3 | https://github.com/suhrig/arriba/wiki/05-Output-files | 2026-06-14 |
| 5 | Genomics England — in-frame/out-of-frame exon-phase rule | 4 | https://www.genomicsengland.co.uk/blog/gene-fusion-reporting | 2026-06-14 |
| 6 | Wikipedia — Reading frame (Badger & Olsen 1999; Lodish 6th ed.) | 4 | https://en.wikipedia.org/wiki/Reading_frame | 2026-06-14 |

### 1.2 Key Evidence Points

1. Fusion support = (# junction/split reads) + (# spanning/discordant fragments) — STAR-Fusion source [2].
2. Default thresholds: ≥1 junction read AND total support ≥2; if 0 junction reads, ≥5 discordant fragments — STAR-Fusion source [2].
3. Total supporting reads = split_reads1 + split_reads2 + discordant_mates — Arriba output spec [4].
4. Split read = read with two segments aligning noncontiguously; discordant mates = paired-end mates aligning nonlinearly — Arriba paper [3].
5. In-frame iff the 3' partner's coding stays in phase across the junction: (5' coding bases − 3' start phase) mod 3 == 0 — Genomics England [5] + reading-frame modulo-3 [6].
6. Detected fusions scored/ordered by abundance of supporting reads — STAR-Fusion / Haas benchmark [1].

### 1.3 Documented Corner Cases

- Spanning-only candidate (0 junction reads): needs ≥5 discordant fragments, else filtered [2].
- Read-through transcripts are common false-positive fusions; distinct-gene and support rules guard against trivial cases [3].
- gene5p == gene3p is not a gene fusion (Registry invariant) [Registry].
- Stop codon before junction (Arriba "stop-codon" value) is out of scope here (requires transcript reconstruction) [4].

### 1.4 Known Failure Modes / Pitfalls

1. Low supporting read count → false positives if threshold too low — STAR-Fusion [2].
2. Counting only split reads or only discordant reads (must sum both classes) — Arriba [4].

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `DetectFusions(IEnumerable<FusionCandidate>, FusionDetectionThresholds?)` | OncologyAnalyzer | Canonical | Rule-based min-support detection; returns ordered fusion calls |
| `IsInFrame(int fivePrimeCodingBases, int threePrimeStartPhase)` | OncologyAnalyzer | Canonical | Codon-phase in-frame test |
| `ComputeTotalSupport(FusionCandidate)` | OncologyAnalyzer | Internal | split1+split2+discordant; exercised via DetectFusions |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported fusion has gene5p ≠ gene3p | Yes | Registry invariant + fusion nomenclature |
| INV-2 | TotalSupport = split_reads1 + split_reads2 + discordant_mates | Yes | Arriba output spec [4] |
| INV-3 | A reported fusion satisfies: junctionReads ≥ minJunctionReads(=1) AND totalSupport ≥ minSumFrags(=2); OR (junctionReads==0 AND discordant ≥ minSpanningFragsOnly(=5)) | Yes | STAR-Fusion source [2] |
| INV-4 | Results are ordered by descending TotalSupport | Yes | STAR-Fusion scoring by abundance [1][2] |
| INV-5 | In-frame ⇔ (fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0 | Yes | exon-phase [5] + modulo-3 [6] |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Junction+sum pass | C1 EML4-ALK split1=3,split2=2,disc=4 | DETECTED, TotalSupport=9 | [2][4] |
| M2 | Spanning-only pass | C4 CD74-ROS1 split=0,disc=5 | DETECTED, TotalSupport=5 | [2] MIN_SPANNING_FRAGS_ONLY=5 |
| M3 | Spanning-only fail | C5 NCOA4-RET split=0,disc=4 | REJECTED (4<5) | [2] |
| M4 | Sum<2 fail | C6 KIF5B-RET split1=1,disc=0 (junc=1,total=1) | REJECTED (1<2) | [2] MIN_SUM_FRAGS=2 |
| M5 | Junc=1,sum=2 pass | C3 TMPRSS2-ERG split1=1,disc=1 | DETECTED, TotalSupport=2 | [2] |
| M6 | Same-gene rejected | C7 ALK-ALK high support | REJECTED (gene5p==gene3p) | INV-1 / Registry |
| M7 | Total support formula | candidate split1=3,split2=2,disc=4 | TotalSupport=9 | [4] |
| M8 | In-frame phase 0 | F1 fivePrimeCodingBases=300, phase=0 | InFrame=true | [5][6] |
| M9 | Out-of-frame phase 1 | F2 fivePrimeCodingBases=301, phase=0 | InFrame=false | [5][6] |
| M10 | In-frame nonzero phase | F4 fivePrimeCodingBases=301, phase=1 | InFrame=true | [5][6] |
| M11 | Ordering by support | mixed candidates | descending TotalSupport order | INV-4 [1][2] |
| M12 | Null candidates | DetectFusions(null) | ArgumentNullException | sibling convention |
| M13 | Negative count | candidate split1=-1 | ArgumentException | sibling convention |
| M14 | Empty input | DetectFusions(empty) | empty result | trivial |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Custom thresholds | min_spanning_frags_only=3, disc=4,split=0 | DETECTED | STAR-Fusion params configurable [2] |
| S2 | Junction boundary | junc=0, disc=5 exactly | DETECTED (≥5) | boundary of MIN_SPANNING_FRAGS_ONLY |
| S3 | Sum boundary | junc=1,total=2 exactly | DETECTED | boundary of MIN_SUM_FRAGS |
| S4 | In-frame phase 2 | fivePrimeCodingBases=302, phase=0 | InFrame=false | [5][6] |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Negative phase guard | IsInFrame(-1,0) | ArgumentException | input validation |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- No existing fusion tests in `tests/Seqeron/Seqeron.Genomics.Tests/`; `FusionDetector`/`DetectFusions` had no prior implementation in `OncologyAnalyzer.cs`. This is a brand-new unit.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M14, S1–S4, C1 | ❌ Missing | New unit; no prior tests |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_DetectFusions_Tests.cs` — all cases.
- **Remove:** none (no prior tests).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| OncologyAnalyzer_DetectFusions_Tests.cs | Canonical fixture | 19 |

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
| 11 | M11 | ❌ Missing | Implemented | ✅ Done |
| 12 | M12 | ❌ Missing | Implemented | ✅ Done |
| 13 | M13 | ❌ Missing | Implemented | ✅ Done |
| 14 | M14 | ❌ Missing | Implemented | ✅ Done |
| 15 | S1 | ❌ Missing | Implemented | ✅ Done |
| 16 | S2 | ❌ Missing | Implemented | ✅ Done |
| 17 | S3 | ❌ Missing | Implemented | ✅ Done |
| 18 | S4 | ❌ Missing | Implemented | ✅ Done |
| 19 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 19
**✅ Done:** 19 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ Covered | DetectFusions_JunctionAndSumPass_Detected |
| M2 | ✅ Covered | DetectFusions_SpanningOnlyFiveFrags_Detected |
| M3 | ✅ Covered | DetectFusions_SpanningOnlyFourFrags_Rejected |
| M4 | ✅ Covered | DetectFusions_SumBelowTwo_Rejected |
| M5 | ✅ Covered | DetectFusions_JunctionOneSumTwo_Detected |
| M6 | ✅ Covered | DetectFusions_SameGene_Rejected |
| M7 | ✅ Covered | ComputeTotalSupport_SumOfThreeClasses (via TotalSupport) |
| M8 | ✅ Covered | IsInFrame_Phase0_True |
| M9 | ✅ Covered | IsInFrame_Phase1_False |
| M10 | ✅ Covered | IsInFrame_NonzeroStartPhase_True |
| M11 | ✅ Covered | DetectFusions_MixedCandidates_OrderedByDescendingSupport |
| M12 | ✅ Covered | DetectFusions_NullInput_Throws |
| M13 | ✅ Covered | DetectFusions_NegativeCount_Throws |
| M14 | ✅ Covered | DetectFusions_EmptyInput_ReturnsEmpty |
| S1 | ✅ Covered | DetectFusions_CustomSpanningThreshold_Detected |
| S2 | ✅ Covered | DetectFusions_SpanningExactlyFive_Detected |
| S3 | ✅ Covered | DetectFusions_SumExactlyTwo_Detected |
| S4 | ✅ Covered | IsInFrame_Phase2_False |
| C1 | ✅ Covered | IsInFrame_NegativeArgument_Throws |

---

## 6. Assumption Register

**Total assumptions:** 2

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Input is candidate-level supporting-read counts (Arriba schema), not raw BAM | Method signature / M1–M7 |
| 2 | In-frame uses codon phase only; premature stop codons (Arriba "stop-codon") out of scope | IsInFrame / M8–M10, S4 |

---

## 7. Open Questions / Decisions

1. Class placement: the session prompt mandates class `OncologyAnalyzer` (the Registry names a future `FusionDetector` class; per the prompt's hard requirement the methods are added to `OncologyAnalyzer`, consistent with the area's existing analyzer). Recorded.
2. `FindChimericReads(bamFile)` and `ValidateFusion(fusion, refGenome)` (Registry "Read extraction"/"Validation" rows) are out of the canonical-threshold scope of this unit and are not implemented here; the count-based decision rule is the formally defined, source-backed core.
