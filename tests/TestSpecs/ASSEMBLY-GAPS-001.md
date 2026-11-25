# Test Specification: ASSEMBLY-GAPS-001

**Test Unit ID:** ASSEMBLY-GAPS-001
**Area:** Assembly
**Algorithm:** Assembly gaps, scaffold structure and contig extraction (N-run detection, gap summary, scaffold → contigs, length filter / sort)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect (B10 sweep U1)
**Last Updated:** 2026-11-25

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | NCBI AGP Specification v2.1 | 2 (format standard) | https://www.ncbi.nlm.nih.gov/genbank/genome_agp_specification/ | 2026-11-25 |
| 2 | NCBI GenBank "Gapped Format for Genome Submissions" | 2 | https://www.ncbi.nlm.nih.gov/genbank/wgs_gapped/ | 2026-11-25 |
| 3 | NCBI GenBank "Genome Submission Guide" (gap information) | 2 | https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/ | 2026-11-25 |
| 4 | QUAST `quast_libs/qutils.py` (`split_by_ns`, `broke_scaffolds`, `correct_fasta`, `correct_seq`), `qconfig.py` | 3 (reference implementation) | https://raw.githubusercontent.com/ablab/quast/master/quast_libs/qutils.py | 2026-11-25 |
| 5 | seqtk `seqtk.c` (`stk_gap`, `stk_cutN`, `seq_nt6_table`) | 3 | https://raw.githubusercontent.com/lh3/seqtk/master/seqtk.c | 2026-11-25 |
| 6 | numpy 2.5.3 `numpy.median` / `numpy.mean` | 3 | pip-installed | 2026-11-25 |

### 1.2 Key Evidence Points

1. A contig contains no gaps "although there may be short runs of Ns due to ambiguous base calls"; AGP coordinates are 1-based inclusive; gap lines must have positive length; column 7b gap types (scaffold, contig, centromere, …) encode linkage evidence, not length (source 1).
2. "Every run of 10 or more Ns is recognized as a gap when the assembly statistics are calculated"; shorter runs are ambiguous bases (`table2asn -gaps-min 10`); 100 N is the suggested unknown-length gap (sources 2, 3).
3. QUAST breaks scaffolds into contigs at maximal N-runs of length ≥ `Ns_break_threshold = 10` (`split_by_ns`); shorter N-runs stay inside the piece; broken pieces are then filtered by `len(seq) >= min_contig` (source 4).
4. QUAST upper-cases (`n` = `N`) and maps IUPAC codes to N (`correct_seq`); seqtk `gap` treats every non-ACGT character as gap and reports BED `[i−l, i)` for runs ≥ `-l` (default 50) (sources 4, 5).
5. seqtk `cutN` (default `-n 1000`, `-p 10`) is a scored N-tract cutter, not exact runs — a different model, not used (source 5).
6. Median = sample median (numpy) (source 6; same convention as ASSEMBLY-STATS-001 F11).

### 1.3 Documented Corner Cases

- Leading / trailing gap: no empty contig is emitted (QUAST emits `seq[0:0]` and drops it via `min_contig`).
- All-N scaffold of length ≥ minGapLength: no contig; shorter than minGapLength: one contig made of Ns (QUAST).
- Lowercase `n` is a gap symbol; IUPAC codes (R, Y, …) are bases (deviation from QUAST/seqtk, consistent with ASSEMBLY-STATS-001).
- Even gap count → median is the mean of the two central lengths.

### 1.4 Known Failure Modes / Pitfalls

1. Splitting a scaffold at N-runs that are below the gap threshold while not recording them as gaps breaks `ContigLength + GapLength = TotalLength` (old AnalyzeScaffolds).
2. Upper median `lengths[count/2]` for even counts (old AnalyzeGapDistribution).
3. Mixing coordinate conventions: Seqeron 0-based inclusive end vs AGP 1-based inclusive vs BED half-open.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindGaps(sequences, minGapLength = 1)` | GenomeAssemblyAnalyzer | Canonical | Maximal N/n runs; also listed in ASSEMBLY-STATS-001 (M10–M13) |
| `AnalyzeGapDistribution(gaps)` | GenomeAssemblyAnalyzer | Canonical | Median via `StatisticsHelper.Median` |
| `AnalyzeScaffolds(scaffolds, minGapLength = 10)` | GenomeAssemblyAnalyzer | Canonical | QUAST `split_by_ns` |
| `ExtractContigs(scaffolds, minContigLength = 200, minGapLength = 10)` | GenomeAssemblyAnalyzer | Canonical | QUAST broken scaffolds + `--min-contig` |
| `FilterByLength(sequences, min, max)` | GenomeAssemblyAnalyzer | Utility | Inclusive bounds |
| `SortByLength(sequences, descending)` | GenomeAssemblyAnalyzer | Utility | Stable |
| MCP `find_gaps`, `gap_distribution`, `analyze_scaffolds`, `extract_contigs` | ChromosomeTools | Wrapper | Delegate |

All four N-run scanners (and `CalculateStatistics` gap fields) share one private `EnumerateNRuns`; scaffold splitting is one private `SplitAtGaps`.

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Gaps = maximal N/n runs of length ≥ minGapLength, 0-based, End inclusive, Length = End − Start + 1 | Yes | Source 5 (`stk_gap`) |
| INV-2 | Scaffold contigs = QUAST `split_by_ns` pieces (non-empty) | Yes | Source 4 |
| INV-3 | ContigLength + GapLength = TotalLength | Yes | Sources 1, 4 |
| INV-4 | ExtractContigs = AnalyzeScaffolds contigs (same minGapLength) with length ≥ minContigLength | Yes | Source 4 |
| INV-5 | FindGaps(default) count / total length = CalculateStatistics TotalGaps / TotalGapLength | Yes | ASSEMBLY-STATS-001 |
| INV-6 | MedianLength = numpy median of gap lengths | Yes | Source 6 |

---

## 4. Test Cases

### 4.1 MUST Tests

| ID | Test Case | Expected Outcome | Test |
|----|-----------|------------------|------|
| M1 | FindGaps vs seqtk `gap` (10 cases: lowercase, IUPAC, leading/trailing, threshold) | seqtk table (Evidence) | `GenomeAssemblyAnalyzer_Gaps_Tests.FindGaps_MatchesSeqtkGap` |
| M2 | IUPAC codes are bases, `n` is a gap | "ARYnNKMA" → (3,4,2) | `FindGaps_IupacCodesAreBases_LowercaseNIsAGap` |
| M3 | Default FindGaps = CalculateStatistics gaps | equal counts / lengths | `FindGaps_DefaultThreshold_CountsTheGapsOfCalculateStatistics` |
| M4 | Gap median / mean = numpy | {2,4} → 3; {1,5,100} → 5; {10,100,1000,5} → 55 | `AnalyzeGapDistribution_MedianAndMean_MatchNumpy` (4), `GenomeAssemblyAnalyzerMutationTests.AnalyzeGapDistribution_ExactSummary` |
| M5 | AnalyzeScaffolds contigs = QUAST `split_by_ns` (9 cases) | pieces table (Evidence) | `AnalyzeScaffolds_ContigsMatchQuastSplitByNs` |
| M6 | Short N-run stays inside contig | "AAAAANNNNNTTTTT",10 → [0,14], 15/15/0 | `AnalyzeScaffolds_NRunBelowMinGap_StaysInsideTheContig`; MCP `AnalyzeScaffolds_ShortGapBelowThreshold_NotSplit` |
| M7 | ExtractContigs default threshold 10 | 300A+9N+300C → one 609 bp; 10 N → 300 + 300; minGap 1 → 300 + 300 | `ExtractContigs_DefaultBreaksOnlyAtTenOrMoreNs`; MCP `ExtractContigs_MinGapLength_DelegatesQuastBreakThreshold` |
| M8a | A10-2: contigs numbered from the full scaffold contig list (dropped short ones keep their number), same names as AnalyzeScaffolds | 5/300/2/250-bp pieces → `s_contig2`, `s_contig4` | `ExtractContigs_NumbersFromTheFullContigList_AgreeingWithAnalyzeScaffolds`, `ExtractContigs_NamesCountDroppedContigs_LikeAnalyzeScaffolds` (MCP) |
| M8 | INV-4 | equal lists | `ExtractContigs_EqualsAnalyzeScaffoldsContigsFilteredByLength` |
| M9 | Validation | null → ArgumentNullException; minGapLength < 1, minContigLength < 0 → ArgumentOutOfRangeException; null element → ArgumentException | `FindGaps_InvalidArguments_Throw`, `AnalyzeScaffolds_InvalidArguments_Throw`, `ExtractContigs_InvalidArguments_Throw`, `AnalyzeGapDistribution_EmptyAndNull`, `FilterByLength_MinAboveMax_IsEmpty_AndNullsThrow` |
| M10 | Random cross-check | regex oracle of `split_by_ns` / `stk_gap` 2000/2000 | Fuzz `AssemblyGapsFuzzTests.GapsScaffoldsContigs_RandomSequences_MatchRegexOracle` |

### 4.2 SHOULD Tests

| ID | Test Case | Expected Outcome | Test |
|----|-----------|------------------|------|
| S1 | ExtractContigs with minGapLength 1 (every N cuts), inclusive length bound, renumbering | as QUAST `split_by_ns(threshold=1)` | `ExtractContigs_SplitsOnNRunsAndFiltersByLength`, `ExtractContigs_LengthBoundaryIsInclusive` |
| S2 | Case preserved, minContigLength 0 keeps all | "acgt"+10n+"g" → acgt, g | `ExtractContigs_LowercaseCasePreserved_MinContigZeroKeepsAll` |
| S3 | SortByLength stable on ties | b,e,a,c,d / d,a,c,b,e | `SortByLength_IsStableOnTies` |
| S4 | FilterByLength inclusive bounds; min > max → empty | | `FilterByLength_InclusiveBounds`, `FilterByLength_MinAboveMax_IsEmpty_AndNullsThrow` |

### 4.3 COULD Tests

| ID | Test Case | Expected Outcome | Test |
|----|-----------|------------------|------|
| C1 | Empty gap list | all zeros, empty map | `AnalyzeGapDistribution_EmptyAndNull` |
| C2 | Leading gap opens no empty contig | | `GenomeAssemblyAnalyzer_MutationKillers7_Tests.AnalyzeScaffolds_LeadingGap_NoEmptyContigAtStart` |

---

## 5. Audit of Existing Tests

| Test | Status | Action |
|------|--------|--------|
| `GenomeAssemblyAnalyzerMutationTests.AnalyzeGapDistribution_ExactSummary` | Locked the upper median (4.0) | Corrected to numpy 3.0 |
| `ExtractContigs_SplitsOnNRunsAndFiltersByLength`, `ExtractContigs_LengthBoundaryIsInclusive` | Relied on cut-at-every-N | Pass `minGapLength: 1` explicitly (QUAST threshold 1); default-10 assertion added |
| MCP `GapDistribution_Binding_InvokesSuccessfully` | Locked median 100 for {5,100} | Corrected to 52.5 |
| Other GenomeAssemblyAnalyzer gap/scaffold tests | Consistent with sources | Kept |

---

## 6. Deviations and Assumptions

| # | Item | Source behaviour | Seqeron behaviour |
|---|------|------------------|-------------------|
| 1 | IUPAC codes | QUAST → N; seqtk → gap | Bases (AGP; consistent with ASSEMBLY-STATS-001) |
| 2 | GapType | AGP 7b = evidence type | Seqeron length class (Short < 10, Medium < 100, Long < 1000, Scaffold); documented, not an AGP type |
| 3 | Contig naming | QUAST numbers dropped pieces too | resolved (A10-2): ExtractContigs numbers every piece, as AnalyzeScaffolds |
| 4 | FindGaps default | NCBI 10, seqtk 50 | 1 (every N-run; = CalculateStatistics); 10 via `DefaultScaffoldGapMinLength` |
| 5 | Coordinates | AGP 1-based inclusive; BED half-open | 0-based, inclusive End (documented conversions) |

## R1-WP4 addendum (2026-11-25): `nonAcgtAsGap` (A3-2)

Deviation 1 (IUPAC codes are bases) is resolved to the sourced reference: QUAST `qutils.correct_seq` turns IUPAC codes into N (seqtk `gap` treats every non-ACGT as a gap), so `FindGaps`, `AnalyzeScaffolds`, `ExtractContigs` (and `CalculateStatistics`) take `nonAcgtAsGap` (default true); false restores the AGP/NCBI reading. QUAST skips a contig file that still has non-ACGTN symbols after `correct_seq`; Seqeron reads them as N (documented lenient deviation). Oracle: `wp4_oracle.stats_expected` (QUAST `correct_seq` + `split_by_ns` verbatim) equals the C# gaps / contigs / extracted pieces on 1500 random cases x 7 seeds. Tests: `FindGaps_MatchesSeqtkGap` rows `ACRYNNKMN`, `ACRYACKMN`, `ACGU-*T` (seqtk transcription), `FindGaps_Default_IupacCodesAndOtherNonAcgtSymbolsAreGapSymbols`, `FindGaps_NonAcgtAsGapFalse_IupacCodesAreBases_LowercaseNIsAGap`, `ScaffoldAndContigApis_NonAcgtAsGap_AppliesConsistently`, fuzz `GapsScaffoldsContigs_RandomSequences_MatchRegexOracle` (both modes), MCP `FindGaps_NonAcgtAsGap_DefaultTrueFollowsQuastCorrectSeq`, `AnalyzeScaffolds_NonAcgtAsGap_SplitsAtIupacRuns`, `ExtractContigs_NonAcgtAsGap_CutsAtIupacRuns`.
