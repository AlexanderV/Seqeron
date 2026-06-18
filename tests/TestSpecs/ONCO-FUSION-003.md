# Test Specification: ONCO-FUSION-003

**Test Unit ID:** ONCO-FUSION-003
**Area:** Oncology
**Algorithm:** Fusion Breakpoint Analysis (junction reading-frame consequence + fusion protein prediction)
**Status:** ☐ In Progress
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-14

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Arriba output-files spec (Uhrig et al. 2021) | 3 | https://github.com/suhrig/arriba/wiki/05-Output-files | 2026-06-14 |
| 2 | AGFusion model.py (Murphy & Elemento 2016) | 3 | https://raw.githubusercontent.com/murphycj/AGFusion/master/agfusion/model.py | 2026-06-14 |
| 3 | Wikipedia "Reading frame" (Badger & Olsen 1999) | 4 | https://en.wikipedia.org/wiki/Reading_frame | 2026-06-14 |

### 1.2 Key Evidence Points

1. Arriba `reading_frame` ∈ {`in-frame`, `out-of-frame`, `stop-codon`, `.`}; the dot means the peptide could not be predicted (breakpoint not in coding context) — Source 1.
2. Arriba `site` categories: `5'UTR`, `3'UTR`, `UTR`, `CDS`, `exon`, `intron`, `intergenic`; a frame call requires both breakpoints in coding context — Source 1.
3. AGFusion frame rule: in-frame iff the chimeric CDS keeps the junction at a codon boundary — the per-segment fractional codon parts must complement; otherwise out-of-frame. The decisive quantity is `(5' coding bases − 3' coding-start phase) mod 3` (reading frames are triplets) — Sources 2, 3.
4. AGFusion protein prediction: chimeric CDS = `cds_5prime_prefix + cds_3prime_suffix`, then `translate()` and truncate at the first stop codon `*` — Source 2.
5. AGFusion out-of-frame: trim the CDS to whole codons (`cds[0:3*(len//3)]`) before translation; the 3' partner is read in a shifted frame — Source 2.

### 1.3 Documented Corner Cases

- Breakpoint outside CDS (UTR/intron/intergenic) → no in/out-of-frame call (Arriba `.`).
- Premature stop in the chimeric ORF → translation truncated at first `*` (AGFusion); Arriba `stop-codon`.
- Out-of-frame fusion → 3' partner read in a frameshifted frame, usually truncated early.

### 1.4 Known Failure Modes / Pitfalls

1. Using gene-level (not CDS-level) lengths for the modulo-3 test gives a wrong frame call — the phase must be the CODING-base count to the breakpoint (Source 2).
2. Forgetting first-stop truncation yields a peptide containing `*` characters (Source 2).

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `AnalyzeBreakpoint(fusion)` | OncologyAnalyzer | Canonical | Classifies the junction: site categories + reading-frame consequence (in/out-of-frame/not-predicted). |
| `PredictFusionProtein(fusion, transcripts)` | OncologyAnalyzer | Canonical | Builds chimeric CDS, translates with the standard genetic code, truncates at first stop. |
| `IsInFrame(fivePrimeCodingBases, threePrimeStartPhase)` | OncologyAnalyzer | Internal | Reused frame primitive (ONCO-FUSION-001); tested indirectly. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | A frame call (`InFrame`/`OutOfFrame`) is made only when both breakpoints lie in CDS; otherwise `NotPredicted`. | Yes | Arriba reading_frame `.` (Source 1) |
| INV-2 | `InFrame` ⟺ `(fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0`. | Yes | AGFusion frame rule + triplet reading (Sources 2, 3) |
| INV-3 | The predicted peptide contains no internal stop `*`; translation ends at the first stop codon. | Yes | AGFusion translate+truncate (Source 2) |
| INV-4 | The chimeric CDS equals 5' CDS prefix (to breakpoint) ++ 3' CDS suffix (from breakpoint). | Yes | AGFusion concatenation (Source 2) |
| INV-5 | An out-of-frame chimeric CDS is trimmed to a whole number of codons before translation. | Yes | AGFusion out-of-frame branch (Source 2) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | InFrame phase 0 | b=9, p=0, both breakpoints CDS | ReadingFrame = InFrame | AGFusion frame rule (Source 2); (9−0)%3=0 |
| M2 | InFrame phase 1 | b=10, p=1 | InFrame | (10−1)%3=0 (Source 2) |
| M3 | InFrame phase 2 | b=11, p=2 | InFrame | (11−2)%3=0 (Source 2) |
| M4 | OutOfFrame | b=10, p=0 | OutOfFrame | (10−0)%3=1≠0 (Source 2) |
| M5 | OutOfFrame phase 1 | b=9, p=1 | OutOfFrame | (9−1)%3=2≠0 (Source 2) |
| M6 | Non-CDS breakpoint → NotPredicted | site5=CDS, site3=5'UTR | ReadingFrame = NotPredicted; no frame call | Arriba `.` (Source 1) |
| M7 | Both UTR → NotPredicted | site5=3'UTR, site3=CDS | NotPredicted | Arriba `.` (Source 1) |
| M8 | Predict protein, in-frame, no stop | cds5=`ATGAAA`, cds3=`GATGGT` | peptide = `MKDG`, Effect=InFrame, PrematureStop=false | AGFusion translate (Source 2) |
| M9 | Predict protein, in-frame, premature stop | cds5=`ATGAAA`, cds3=`GATTAAGGT` | peptide = `MKD`, PrematureStop=true | AGFusion translate+truncate (Source 2) |
| M10 | Predict protein, mid-codon junction, 3' frameshifted | cds5=`ATGA`(4, phase 1), cds3=`AAGGT` (phase 0) | chimeric `ATGAAAGGT` → peptide = `MKG`; Effect=**OutOfFrame** ((4−0)%3=1≠0; the 3' gene is read frameshifted under the Arriba 3'-gene-frame model) | Arriba reading_frame (Source 1); AGFusion concat (Source 2) |
| M10b | Predict protein, mid-codon junction, frames compatible | cds5=`ATGA`(4, phase 1), 3' CDS `TAAGGT` sliced at phase 1 → suffix `AAGGT` | chimeric `ATGAAAGGT` → peptide = `MKG`; Effect=**InFrame** ((4−1)%3=0; 3' read in native frame) | Arriba reading_frame (Source 1); AGFusion concat (Source 2) |
| M11 | Predict protein, out-of-frame trims to whole codons | cds5=`ATGAA`(5), cds3=`GATGGT`(6) | chimeric trimmed `ATGAAGATG` → peptide = `MKM`, Effect=OutOfFrame | AGFusion out-of-frame branch (Source 2) |
| M12 | Chimeric CDS composition | cds5=`ATGAAA`, cds3=`GATGGT` | ChimericCds = `ATGAAAGATGGT` | AGFusion concat (Source 2) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Premature stop → ReadingFrame StopCodon | in-frame fusion whose ORF stops before end | StopCodon reported / PrematureStop flag true | Arriba `stop-codon` value |
| S2 | Invalid 3' phase | phase = 3 | ArgumentOutOfRangeException | mirrors IsInFrame validation |
| S3 | Empty 3' CDS suffix → peptide from 5' only | cds3 contributes nothing after offset | peptide = translation of 5' prefix | concat with empty suffix |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Integration with FusionCall | breakpoint analysis consumes a FusionCall's partners | designation partners preserved | API consistency |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Searched `tests/Seqeron/Seqeron.Genomics.Tests/` for `AnalyzeBreakpoint` / `PredictFusionProtein` — none exist (new unit).
- Sibling fusion fixtures: `OncologyAnalyzer_DetectFusions_Tests.cs` (ONCO-FUSION-001), `OncologyAnalyzer_MatchKnownFusions_Tests.cs` (ONCO-FUSION-002).

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1–M12 | ❌ Missing | New unit; no existing tests. |
| S1–S3 | ❌ Missing | New unit. |
| C1 | ❌ Missing | New unit. |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/OncologyAnalyzer_AnalyzeBreakpoint_Tests.cs` — all ONCO-FUSION-003 tests.
- **Remove:** none (no prior tests for this unit).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `OncologyAnalyzer_AnalyzeBreakpoint_Tests.cs` | Canonical fixture for ONCO-FUSION-003 | 16 |

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
| 13 | S1 | ❌ Missing | Implemented | ✅ Done |
| 14 | S2 | ❌ Missing | Implemented | ✅ Done |
| 15 | S3 | ❌ Missing | Implemented | ✅ Done |
| 16 | C1 | ❌ Missing | Implemented | ✅ Done |

**Total items:** 16
**✅ Done:** 16 | **⛔ Blocked:** 0 | **Remaining:** must be 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | AnalyzeBreakpoint_InFramePhase0_InFrame |
| M2 | ✅ | AnalyzeBreakpoint_InFramePhase1_InFrame |
| M3 | ✅ | AnalyzeBreakpoint_InFramePhase2_InFrame |
| M4 | ✅ | AnalyzeBreakpoint_OutOfFrame_OutOfFrame |
| M5 | ✅ | AnalyzeBreakpoint_OutOfFramePhase1_OutOfFrame |
| M6 | ✅ | AnalyzeBreakpoint_ThreePrimeUtr_NotPredicted |
| M7 | ✅ | AnalyzeBreakpoint_FivePrimeUtr_NotPredicted |
| M8 | ✅ | PredictFusionProtein_InFrameNoStop_TranslatesFullPeptide |
| M9 | ✅ | PredictFusionProtein_PrematureStop_TruncatesAtStop |
| M10 | ✅ | PredictFusionProtein_MidCodonJunctionPhaseMismatch_OutOfFrame |
| M10b | ✅ | PredictFusionProtein_MidCodonJunctionPhaseMatch_InFrame |
| M11 | ✅ | PredictFusionProtein_OutOfFrame_TrimsAndShifts |
| M12 | ✅ | PredictFusionProtein_ChimericCds_ConcatenatesPrefixSuffix |
| S1 | ✅ | PredictFusionProtein_PrematureStop_FlagsStopCodon |
| S2 | ✅ | AnalyzeBreakpoint_InvalidPhase_Throws |
| S3 | ✅ | PredictFusionProtein_EmptyThreePrimeSuffix_PeptideFromFivePrime |
| C1 | ✅ | AnalyzeBreakpoint_PreservesPartners |

---

## 6. Assumption Register

**Total assumptions:** 1

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | Caller supplies partner CDS strings and breakpoint CDS offsets (no genome/GTF DB in repo); API-shape only, no output change. | PredictFusionProtein contract |

---

## 7. Open Questions / Decisions

1. None — frame rule, site categories, protein prediction, and first-stop truncation are all fixed by retrieved sources (Arriba spec, AGFusion source).
</content>
