# Test Specification: GENOMIC-ORF-001

**Test Unit ID:** GENOMIC-ORF-001
**Area:** Analysis
**Algorithm:** Open Reading Frame (ORF) Detection — `GenomicAnalyzer.FindOpenReadingFrames`
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-11-18 (B09 review F21-F26)

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | Rosalind — ORF problem (worked dataset) | 4 | https://rosalind.info/problems/orf/ | 2026-06-14 |
| 2 | Wikipedia — Open reading frame | 4 | https://en.wikipedia.org/wiki/Open_reading_frame | 2026-06-14 |
| 3 | NCBI ORFfinder (tool reference) | 5 | https://www.ncbi.nlm.nih.gov/orffinder/ | 2026-06-14 |
| 4 | NCBI Genetic Codes (Standard, transl_table=1) | 2 | https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi | 2026-06-14 |
| 5 | NCBI C++ Toolkit `COrf::FindOrfs` (`src/algo/sequence/orf.cpp`, `orf.hpp`) — the ORFfinder engine | 1 (tool source) | https://github.com/ncbi/ncbi-cxx-toolkit-public | 2026-11-18 |
| 6 | EMBOSS getorf 6.6.0 (binary run locally, `-find 3`) | 1 (tool) | http://emboss.open-bio.org/rel/rel6/apps/getorf.html | 2026-11-18 |

### 1.2 Key Evidence Points

1. An ORF begins with a start codon and ends with a stop codon, with no internal in-frame stop — source 1, 2.
2. Six reading frames: 3 forward + 3 on the reverse complement — source 1, 2.
3. Standard start codon ATG; stop codons TAA/TAG/TGA — source 1, 4.
4. Every ATG that opens an ORF terminated by an in-frame stop is reported; nested ORFs sharing a stop are all returned (`MGMTPRLGLESLLE` and `MTPRLGLESLLE`) — source 1.
5. A reading begun at an ATG without a downstream in-frame stop yields no protein candidate (translate "until a stop codon") — source 1.
6. `minLength` is a nucleotide threshold; ORFs with length ≥ threshold are kept — source 3. **Correction (B09):** NCBI ORFfinder's own threshold (`COrf`, source 5: "length in base pairs not including the stop") and getorf `-minsize` (source 6) do *not* count the stop codon; this method's default counts it (`Length ≥ minLength`), the NCBI convention is `OrfSearchOptions.MinLengthExcludesStop`.
7. NCBI `COrf` reports every ORF from its most upstream start (`longest_orfs`), or every start (`longest_orfs = false` = Rosalind), reports 3'-open ORFs as partial, takes stop/start codons from the genetic code table (`-g`), and has start modes ATG only / ATG + alternative / any sense codon — source 5.

### 1.3 Documented Corner Cases

- Nested ORFs sharing a stop codon → both reported (source 1).
- ATG with no in-frame stop → not a complete ORF (source 1).
- ORFs on the reverse complement strand (source 1, 2).
- minLength filtering (source 3).

### 1.4 Known Failure Modes / Pitfalls

1. Greedy "first ATG to stop, then skip to next ATG" scanning misses nested ORFs sharing a stop — contradicts source 1 (the pre-existing implementation had this defect).
2. Forgetting the reverse strand misses up to half of ORFs — source 1, 2.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindOpenReadingFrames(DnaSequence, int minLength=100)` | GenomicAnalyzer | Canonical | Six-frame ATG→stop ORF detection; every ATG start reported |
| `FindOpenReadingFrames(string, int minLength[, OrfSearchOptions])` | GenomicAnalyzer | Canonical (A1-21) | IUPAC nucleotide string (DNA/RNA, case-insensitive, `U`→`T`; other symbol → `ArgumentException`): ambiguous codon = stop only if every concrete codon is a stop, never a start; oracle `orf_ambiguous_oracle.tsv` (2844 rows) |
| `FindOpenReadingFrames(DnaSequence, int minLength, OrfSearchOptions)` | GenomicAnalyzer | Canonical (options) | genetic code, start-codon mode, nested vs. longest per stop, open ends, length convention, strands (NCBI COrf / getorf) |
| `OrfScanner.FindRegions` | Core | Shared core | single-strand O(n) scanner used by both and by `Translator.FindOrfs` |
| `OrfInfo.GetForwardStrandRange(int)` | GenomicAnalyzer | Helper | input-strand coordinates of an ORF |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Every reported ORF `Sequence` begins with `ATG` | Yes | Source 1, 4 |
| INV-2 | Every reported ORF `Sequence` ends with TAA/TAG/TGA | Yes | Source 1, 4 |
| INV-3 | Every reported ORF `Length` is divisible by 3 | Yes | Source 2 |
| INV-4 | Every reported ORF `Length` ≥ `minLength` | Yes | Source 3 |
| INV-5 | `Frame` ∈ {1,2,3}; `IsReverseComplement` distinguishes the strand | Yes | Source 1, 2 |
| INV-6 | The stop set is the genetic code's stop set (table 2: AGA/AGG stop, TGA does not) | Yes | Source 4, 5 |
| INV-7 | `GetForwardStrandRange`: reverse complement of `input[Start..End]` = `Sequence` for reverse ORFs, `input[Start..End]` for forward | Yes | definition; Biopython `reverse_complement` |
| INV-8 | Nested = false ⇒ per frame at most one ORF per stop; every such ORF is the longest among the nested set | Yes | Source 5 |
| INV-9 | An open-ended ORF has `IsTerminated = false`, `Length % 3 == 0`, no stop codon inside | Yes | Source 5, 6 |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | Single forward ORF | `ATGAAAAAATAA`, minLength 1 | 1 ORF: Sequence `ATGAAAAAATAA`, Position 0, Frame 1, !RC | Source 4 derivation |
| M2 | Rosalind six-frame dataset | `Rosalind_99` input, minLength 1 | Distinct translated proteins == {`MLLGSFRLIPKETLIQVAGSSPCNLS`,`M`,`MGMTPRLGLESLLE`,`MTPRLGLESLLE`} | Source 1 |
| M3 | Nested ORFs share stop | `ATGGGGATGCCCTAA`, minLength 1 | 2 ORFs at positions 0 and 6, both ending `TAA` | Source 1 |
| M4 | ATG without in-frame stop | `ATGAAAAAAAAA`, minLength 1 | 0 ORFs | Source 1 |
| M5 | Reverse-complement-only ORF | seq whose RC contains `ATG...TAA`, none forward | ORF with IsReverseComplement=true | Source 1, 2 |
| M6 | minLength excludes short ORF | `ATGAAATAA` (9nt), minLength 12 | 0 ORFs | Source 3 |
| M6b | minLength includes exactly-at-threshold | `ATGAAATAA` (9nt), minLength 9 | 1 ORF | Source 3 |
| M7 | INV-1 starts with ATG | Rosalind input | all ORFs start `ATG` | Source 1, 4 |
| M8 | INV-2 ends with stop | Rosalind input | all ORFs end TAA/TAG/TGA | Source 1, 4 |
| M9 | INV-3 length % 3 == 0 | Rosalind input | all lengths divisible by 3 | Source 2 |
| M10 | All three stop codons recognized | `ATGTAA`,`ATGTAG`,`ATGTGA`, minLength 1 | each yields 1 ORF | Source 4 |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | Lowercase input | `atgaaaaaataa` | same single ORF as M1 (case-insensitive) | DnaSequence normalizes case |
| S2 | Multiple ORFs different frames | two ORFs in different frames | both detected with correct frame numbers | Six-frame requirement |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | Too-short sequence | `ATG`, minLength 1 | 0 ORFs (no stop) | Edge guard |

### 4.4 Review additions (B09, F21–F26) — `GenomicAnalyzer_FindOpenReadingFrames_Oracle_Tests.cs`

| ID | Test Case | Expected Outcome | Evidence |
|----|-----------|------------------|----------|
| R1 | 1746-row oracle (all NCBI tables except the dual-coding 27, 28, 31; atg/alt/any; nested/longest; open/closed; min length with/without stop; strands) | exact ordered list `strand,frame,start,end,terminated` | Biopython finder = NCBI COrf port = getorf 6.6.0 (377 rows) |
| R2 | Default overload equals default options on every default row | identical | — |
| R3 | `GetForwardStrandRange`: revcomp of span = ORF (200 random sequences); worked `CC+TTAGGGCAT+GG` → (2,10), `TTAGGGGGGCAT` → (0,11); too-short length throws | exact | INV-7 |
| R4 | String overload (A1-21): 2844 N/IUPAC rows (Biopython definition; getorf 6.6.0 on the 184 A/C/G/T/N rows), worked examples (`ATGNNNTAA`, `TAR`, `NTG`, RNA, reverse strand `TTANNNCAT`), A/C/G/T string == `DnaSequence` overload, invalid symbol / null / dual-coding table | exact | INV-1..7 on IUPAC input |
| R4 | Table 2: `ATGTGAAAATAA` ORF is the whole 12 nt (TGA = Trp), standard code → `ATGTGA`; `ATGCCCAGAGGG` table 2 → `ATGCCCAGA`, table 1 → none | exact | Biopython CodonTable |
| R5 | Tables 27/28/31 throw eagerly; null sequence / options throw | — | dual-coding stops |
| R6 | Starts: `TTG/CTGAAATAA` table 1 alt; `GTG` table 11 only; any sense `CCCAAATAAGGG` (+,1,0,8 ; nested adds 3..8) | exact | NCBI gc.prt, COrf |
| R7 | `NestedOrfs = false` on `ATGGGGATGCCCTAA` → 1 ORF at 0 | exact | COrf longest_orfs |
| R8 | `IncludeOpenEnded` on `ATGAAACCCGG` → `ATGAAACCC`, `IsTerminated = false` | exact | COrf partial, getorf |
| R9 | `MinLengthExcludesStop`: `ATGAAATAA` min 6 kept, 7 dropped; default min 9 kept, 10 dropped | exact | COrf, getorf |
| R10 | Differential vs `Translator.FindOrfs` (tables 1, 2, 4, 11; 240 sequences) = outermost starts of the code, open ends, aa-length min | identical spans/frames | B02 getorf-confirmed engine |
| R11 | Differential vs `GenomeAnnotator.FindOrfs` (150 sequences) ATG-initiated ORFs = ours (forward coords, half-open) | identical | B11 engine |
| R12 | Complexity: `ATG`×300 000 + `TAA`, min 899 000 → 335 ORFs in ≪ 30 s; 5 Mbp random six frames < 60 s | counts | Rosalind definition |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzerTests.cs` contained 3 pre-existing smoke tests under `#region Open Reading Frames` (`FindOpenReadingFrames_SimpleOrf_FindsIt`, `_MultipleFrames_FindsAll`, `_NoOrf_ReturnsEmpty`). These used permissive assertions (`Has.Count.GreaterThanOrEqualTo`, `o => o.Sequence.StartsWith`) and did not pin exact ORFs/positions.
- ANNOT-ORF-001 (`GenomeAnnotator.FindOrfs`) is the related but separate canonical ORF unit; it classified `GenomicAnalyzer.FindOpenReadingFrames` as an "Alternate (smoke)". GENOMIC-ORF-001 now makes `FindOpenReadingFrames` the canonical method for its own unit with deep evidence-based tests in a dedicated fixture.

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Notes |
|---------------------|--------|-------|
| M1 single ORF | ❌ Missing | new fixture |
| M2 Rosalind dataset | ❌ Missing | new fixture |
| M3 nested ORFs | ❌ Missing | pre-existing greedy impl could not pass this |
| M4 no in-frame stop | ❌ Missing | |
| M5 reverse strand | ❌ Missing | old `_MultipleFrames` was weak |
| M6 / M6b minLength | ❌ Missing | |
| M7–M9 invariants | ❌ Missing | |
| M10 stop codons | ❌ Missing | |
| S1 lowercase | ❌ Missing | |
| S2 multiple frames | ⚠ Weak | old `_MultipleFrames_FindsAll` only asserts StartsWith ATG |
| C1 too short | ⚠ Weak | old `_NoOrf_ReturnsEmpty` close but not exact-case |

### 5.3 Consolidation Plan

- **Canonical file:** `tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_FindOpenReadingFrames_Tests.cs` — all deep evidence-based tests for this unit.
- **Remove:** the 3 weak smoke tests and the `#region Open Reading Frames` block from `GenomicAnalyzerTests.cs` (superseded by the canonical fixture; avoids duplicate/weak coverage).

### 5.4 Final State After Consolidation

| File | Role | Test Count |
|------|------|------------|
| `GenomicAnalyzer_FindOpenReadingFrames_Tests.cs` | Canonical | 15 (+ `…_Oracle_Tests.cs`, 18, B09) |
| `GenomicAnalyzerTests.cs` | Other GenomicAnalyzer methods | ORF region removed |

### 5.5 Phase 7 Work Queue

| # | Test Case ID | §5.2 Status | Action Taken | Final Status |
|---|-------------|-------------|--------------|--------------|
| 1 | M1 | ❌ | implemented | ✅ Done |
| 2 | M2 | ❌ | implemented | ✅ Done |
| 3 | M3 | ❌ | implemented | ✅ Done |
| 4 | M4 | ❌ | implemented | ✅ Done |
| 5 | M5 | ❌ | implemented | ✅ Done |
| 6 | M6/M6b | ❌ | implemented | ✅ Done |
| 7 | M7–M9 | ❌ | implemented | ✅ Done |
| 8 | M10 | ❌ | implemented | ✅ Done |
| 9 | S1 | ❌ | implemented | ✅ Done |
| 10 | S2 | ⚠ | rewritten exact | ✅ Done |
| 11 | C1 | ⚠ | rewritten exact | ✅ Done |

**Total items:** 11
**✅ Done:** 11 | **⛔ Blocked:** 0 | **Remaining:** 0

### 5.6 Post-Implementation Coverage

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1 | ✅ | exact Sequence/Position/Frame |
| M2 | ✅ | exact 4-protein set |
| M3 | ✅ | exact 2 positions, shared stop |
| M4 | ✅ | empty |
| M5 | ✅ | IsReverseComplement true |
| M6 | ✅ | empty |
| M6b | ✅ | 1 ORF |
| M7 | ✅ | invariant |
| M8 | ✅ | invariant |
| M9 | ✅ | invariant |
| M10 | ✅ | TestCase per stop codon |
| S1 | ✅ | lowercase equals uppercase result |
| S2 | ✅ | two frames |
| C1 | ✅ | empty |

All in-scope cases ✅ (14 tests).

---

## 6. Assumption Register

**Total assumptions:** 3

| # | Assumption | Used In |
|---|-----------|---------|
| 1 | ORF `Sequence` includes the terminating stop codon (protein candidate excludes it) | M1, M3, M6b, INV-2/3 |
| 2 | `minLength` is a nucleotide threshold, inclusive lower bound, counted **with** the stop codon by default (NCBI/getorf count without: option) | M6, M6b, INV-4 |
| 3 | Standard start codon ATG only (alternative initiation codons out of scope) | all |

These affect reported span/filtering only, not which start sites are detected; each is backed by an authoritative convention (NCBI nucleotide length filter; Wikipedia "length divisible by three, bounded by stop codons"; NCBI default "ATG only").

---

## 7. Open Questions / Decisions

1. Relationship to ANNOT-ORF-001: kept as a separate canonical unit (`GenomicAnalyzer.FindOpenReadingFrames`). The pre-existing greedy implementation was a defect (missed nested ORFs) and was corrected to the canonical every-ATG-to-first-in-frame-stop definition. Not consolidated away because the Registry lists it as its own unit with a distinct struct/API (`OrfInfo`).
2. Suffix tree: not applicable — this is a fixed-stride six-frame codon scan (O(n)), not query-occurrence enumeration. Decision recorded in algorithm doc §5.2.
