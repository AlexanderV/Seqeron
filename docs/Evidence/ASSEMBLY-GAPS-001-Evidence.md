# Evidence Artifact: ASSEMBLY-GAPS-001

**Test Unit ID:** ASSEMBLY-GAPS-001
**Algorithm:** Assembly gaps, scaffold structure and contig extraction (`GenomeAssemblyAnalyzer.FindGaps`, `AnalyzeGapDistribution`, `AnalyzeScaffolds`, `ExtractContigs`, `FilterByLength`, `SortByLength`; MCP `find_gaps`, `gap_distribution`, `analyze_scaffolds`, `extract_contigs`)
**Date Collected:** 2026-11-25 (B10 sweep U1)

---

## Online Sources

### NCBI — AGP Specification v2.1

**URL:** https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/ (redirects to https://www.ncbi.nlm.nih.gov/genbank/genome_agp_specification/)
**Accessed:** 2026-11-25 (curl, neutral User-Agent; HTML stripped to text and read)
**Authority rank:** 2 (format standard)

**Key Extracted Points:**

1. Contig definition (verbatim): "There should be no gaps in a sequence contig (although there may be short runs of Ns due to ambiguous base calls)." → short N-runs belong to the contig.
2. "All coordinates given in the file are 1-based inclusive (not 0-based)." "Gap lengths must be positive. Negative gaps and gap lines with zero length" are invalid.
3. Column 7b `gap_type` = scaffold / contig / centromere / short_arm / heterochromatin / telomere / repeat …, combined with column 8b linkage evidence — the type is evidence-based, not a length class.
4. "A length of 100 must be used for all U type gaps."

### NCBI GenBank — "Gapped Format for Genome Submissions"

**URL:** https://www.ncbi.nlm.nih.gov/genbank/wgs_gapped/
**Accessed:** 2026-11-25 (curl, text read)
**Authority rank:** 2

1. (verbatim) "Note that every run of 10 or more Ns is recognized as a gap when the assembly statistics are calculated in NCBI's Assembly resource."
2. `table2asn -gaps-min`: "the minimum number of Ns in a row that represents a gap"; example "if runs of 10 or more N's are estimated gaps, and shorter runs of N's are just ambiguous bases … -gaps-min 10 -gaps-unknown 100".
3. "Do not add assembly_gaps for Ns that represent ambiguous base calls".

### NCBI GenBank — "Genome Submission Guide"

**URL:** https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/
**Accessed:** 2026-11-25

1. (verbatim) "Be aware that the assembly statistics are always calculated using 10 or more Ns as a gap, regardless of the presence/absence of gaps in the final genome sequence."

### QUAST — `quast_libs/qutils.py`, `quast_libs/qconfig.py` (ablab/quast master)

**URL:** https://raw.githubusercontent.com/ablab/quast/master/quast_libs/qutils.py, …/qconfig.py
**Accessed:** 2026-11-25 (raw source downloaded, functions read in full)
**Authority rank:** 3 (reference implementation)

1. `qconfig.Ns_break_threshold = 10`; `broke_scaffolds` calls `split_by_ns(seq, name, broken_scaffolds_fasta, qconfig.Ns_break_threshold, qconfig.min_contig)`.
2. `split_by_ns` (transcribed verbatim into the oracle): finds each maximal `N` run `[start, end)`; `if end - start >= Ns_break_threshold:` the piece `seq[cur_contig_start:start]` is emitted and `cur_contig_start = end`; shorter runs stay in the piece; the last piece is emitted if `len(seq) - cur_contig_start >= min_contig`.
3. The broken file is re-read by `correct_fasta` with `min_threshold = qconfig.min_contig` → `if (len(seq) >= min_contig)` (inclusive length filter; drops the empty piece produced by a leading gap).
4. `correct_seq`: `seq.upper()` and IUPAC M/K/R/Y/W/S/V/B/H/D → N.

### seqtk — `seqtk.c` (lh3/seqtk master)

**URL:** https://raw.githubusercontent.com/lh3/seqtk/master/seqtk.c
**Accessed:** 2026-11-25

1. `stk_gap` (`seqtk gap [-l 50]`): runs of characters with `seq_nt6_table` ∉ {1..4} (i.e. not A/C/G/T, either case) of length `l ≥ min_size` are printed as `name  i-l  i` (0-based half-open).
2. `stk_cutN` (`-n 1000` min N tract, `-p 10` non-N penalty) extends N tracts by a score that tolerates isolated non-N bases — a different, scored model; not adopted.

### numpy 2.5.3 — `numpy.median`, `numpy.mean` (pip-installed)

Sample median: mean of the two central values for an even count.

---

## Reference numbers (oracle `tools/review/oracles/B10/gaps_oracle.py`, `python3 -I`)

### FindGaps vs seqtk `stk_gap` (End inclusive; IUPAC mapped to bases before seqtk)

| Sequence | minGap | Gaps (Start, End, Length) |
|---|---|---|
| NNAANNNN | 1 | (0,1,2) (4,7,4) |
| AAANNNNNGGG | 1 / 10 | (3,7,5) / — |
| nnACnNnGT | 1 | (0,1,2) (4,6,3) |
| ACRYNNKMN | 1 / 2 | (4,5,2) (8,8,1) / (4,5,2) |
| NNNNNNNNNN | 1 | (0,9,10) |
| "" ; ACGT | 1 | — |
| A+10N+C+9N+G | 10 | (1,10,10) |

### AnalyzeScaffolds / ExtractContigs vs QUAST `split_by_ns` (pieces, 0-based half-open)

| Sequence | minGap | QUAST pieces | Old Seqeron contigs |
|---|---|---|---|
| AAAAANNNNNTTTTT | 10 | [0,15) | [0,5) [10,15), no gap, ContigLength 10 ≠ 15 |
| AAANNNNNNNNNNAAA | 10 | [0,3) [13,16) | same |
| NNNNNNNNNNAAA | 10 | [10,13) | same |
| AAAA+10N | 10 | [0,4) | same |
| ACnnGT+12N+acNNNNNNNNNt | 10 | [0,6) [18,30) | split also at `nn` and the 9-N run |
| NNN | 10 | [0,3) | none |
| 10N | 10 | — | — |
| AAANNAAAA | 1 / 3 | [0,3) [5,9) / [0,9) | [0,3) [5,9) for both |
| 300A+9N+300C (ExtractContigs, 200) | 10 | one 609-bp contig | two 300-bp contigs |

### AnalyzeGapDistribution vs numpy

| Lengths | numpy median | numpy mean | Old median |
|---|---|---|---|
| {2,4} | 3.0 | 3.0 | 4 |
| {5,100} | 52.5 | 52.5 | 100 |
| {1,5,100} | 5.0 | 35.333333333333336 | 5 |
| {10,100,1000,5} | 55.0 | 278.75 | 100 |

### Random cross-check

3000 random sequences (seed 2026; length 0–60; N/n runs of 1–25; ACGTacgt; minGap ∈ {1,2,3,5,10,11,20}; minContig ∈ {0,1,2,3,5,8,20}; 1393 with an N-run shorter than minGap) — oracle JSON compared against the C# methods (FindGaps, AnalyzeScaffolds contigs + gaps, ExtractContigs, median, mean, INV-3): **3000/3000**. The oracle asserts its coordinate transcription equals the verbatim `split_by_ns` slices on every case. A C# regex oracle (alphabet with IUPAC R/Y) is locked as fuzz test `AssemblyGapsFuzzTests` (2000 cases).

---

## Findings (B10 F36–F40)

- **F36** AnalyzeScaffolds cut contigs at every N-run but recorded only runs ≥ minGapLength as gaps → sub-threshold N-runs vanished (ContigLength + GapLength < TotalLength). Fixed to QUAST `split_by_ns`.
- **F37** ExtractContigs cut at every single N (QUAST threshold 1) while NCBI/QUAST break scaffolds at ≥ 10 N; new `minGapLength` parameter, default `DefaultScaffoldGapMinLength = 10`; consistent with AnalyzeScaffolds.
- **F38** AnalyzeGapDistribution median was the upper median → `StatisticsHelper.Median`.
- **F39** Validation: null inputs threw NRE lazily, minGapLength ≤ 0 accepted silently → eager ArgumentNullException / ArgumentOutOfRangeException, null element → ArgumentException.
- **F40** MCP docs/descriptions (upper median, "gap-free runs", coordinates, GapType meaning) corrected; `extract_contigs` gains `minGapLength`.

## Deviations kept (documented)

IUPAC codes are bases (QUAST → N, seqtk → gap); GapType is a Seqeron length class, not an AGP 7b type (the sequence carries no linkage evidence); FindGaps default 1 (every N-run, matches CalculateStatistics); ExtractContigs numbers only returned contigs.

**R1-WP2 (2026-11-25):** QUAST `qutils.py` `split_by_ns` re-opened: `total_contigs` is incremented for every piece appended and `--min-contig` filters afterwards, so `ExtractContigs` now names `{id}_contig{n}` with n = position among all contigs of the scaffold (same as `AnalyzeScaffolds`); locked by `ExtractContigs_NumbersFromTheFullContigList_AgreeingWithAnalyzeScaffolds` (pieces 5/300/2/250 bp → `s_contig2`, `s_contig4`). Oracle: `tools/review/oracles/B10/stats_oracle.py` header.
