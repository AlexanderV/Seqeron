# bio-assembly tool map (assembly-relevant subset)

Grouped by sub-task. Each row: **[MCP] tool** · `Method ID` · one-line purpose. Open the linked
per-tool doc for the full I/O schema — do not guess parameters.

**Server split.** This skill's owned slice is **Analysis** (k-mer + repeat/complexity subset) +
**Core** (suffix-tree repeats). The assembly **engine** tools are on the **Alignment** server
(`SequenceAssembler.*`) — listed here because they drive the workflow, but their reference lives with
[`bio-alignment`](../../bio-alignment/SKILL.md). The Analysis server is **shared with `bio-annotation`**: everything
motif/variant/ORF/RNA/protein on Analysis is `bio-annotation`'s, not listed here.

## Assembly engine (Alignment server — `SequenceAssembler.*`)

| Tool | Method ID | Purpose |
|---|---|---|
| [`assemble_de_bruijn`](../../../../docs/mcp/tools/alignment/assemble_de_bruijn.md) | `SequenceAssembler.AssembleDeBruijn` | Shred reads into k-mers, build (k−1)-mer graph, Eulerian walk → contigs + N50. |
| [`assemble_olc`](../../../../docs/mcp/tools/alignment/assemble_olc.md) | `SequenceAssembler.AssembleOLC` | Overlap-layout-consensus assembly (long clear overlaps) → contigs + N50. |
| [`find_overlap`](../../../../docs/mcp/tools/alignment/find_overlap.md) | `SequenceAssembler.FindOverlap` | Best suffix-prefix overlap between two reads. |
| [`find_all_overlaps`](../../../../docs/mcp/tools/alignment/find_all_overlaps.md) | `SequenceAssembler.FindAllOverlaps` | All pairwise overlaps in a read set (inspect the graph pre-assembly). |
| [`merge_contigs`](../../../../docs/mcp/tools/alignment/merge_contigs.md) | `SequenceAssembler.MergeContigs` | Merge overlapping contigs. |
| [`scaffold_contigs`](../../../../docs/mcp/tools/alignment/scaffold_contigs.md) | `SequenceAssembler.Scaffold` | Order/orient contigs into scaffolds. |
| [`error_correct_reads`](../../../../docs/mcp/tools/alignment/error_correct_reads.md) | `SequenceAssembler.ErrorCorrectReads` | k-mer-spectrum read error correction (pre-assembly). |
| [`quality_trim_reads`](../../../../docs/mcp/tools/alignment/quality_trim_reads.md) | `SequenceAssembler.QualityTrimReads` | Quality-based read trimming (pre-assembly). |
| [`compute_consensus`](../../../../docs/mcp/tools/alignment/compute_consensus.md) | `SequenceAssembler.ComputeConsensus` | Majority-vote consensus from pre-aligned equal-length reads. |

## Coverage & assembly stats (Alignment server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`assembly_stats`](../../../../docs/mcp/tools/alignment/assembly_stats.md) | `SequenceAssembler.CalculateStats` | N50 / longest / total length / read accounting for a contig set. |
| [`calculate_coverage`](../../../../docs/mcp/tools/alignment/calculate_coverage.md) | `SequenceAssembler.CalculateCoverage` | Per-base coverage depth array of reads over a reference. |

> Chromosome-scale assembly QC (auN, Nx, L50, gaps, BUSCO-like completeness, `extract_contigs`) lives
> on the **Chromosome** server — see [`bio-chromosome`](../../bio-chromosome/SKILL.md).

## k-mer (Analysis server — `KmerAnalyzer.*`)

| Tool | Method ID | Purpose |
|---|---|---|
| [`kmer_spectrum`](../../../../docs/mcp/tools/analysis/kmer_spectrum.md) | `KmerAnalyzer.GetKmerSpectrum` | Frequency-of-frequencies; separates error k-mers from genomic, hints coverage. |
| [`count_kmers`](../../../../docs/mcp/tools/analysis/count_kmers.md) | `KmerAnalyzer.CountKmers` | Exact per-k-mer counts (single strand). |
| [`count_kmers_both_strands`](../../../../docs/mcp/tools/analysis/count_kmers_both_strands.md) | `KmerAnalyzer.CountKmersBothStrands` | Canonical both-strand k-mer counts (correct for dsDNA). |
| [`most_frequent_kmers`](../../../../docs/mcp/tools/analysis/most_frequent_kmers.md) | `KmerAnalyzer.FindMostFrequentKmers` | Top-count k-mers (repeats / adapters / contamination). |
| [`kmers_with_min_count`](../../../../docs/mcp/tools/analysis/kmers_with_min_count.md) | `KmerAnalyzer.FindKmersWithMinCount` | k-mers occurring ≥ minCount times. |
| [`unique_kmers`](../../../../docs/mcp/tools/analysis/unique_kmers.md) | `KmerAnalyzer.FindUniqueKmers` | k-mers occurring exactly once (unique anchors/markers). |
| [`find_clumps`](../../../../docs/mcp/tools/analysis/find_clumps.md) | `KmerAnalyzer.FindClumps` | k-mers clumping ≥ minOccurrences within a window (ori/satellite arrays). |
| [`kmer_frequencies`](../../../../docs/mcp/tools/analysis/kmer_frequencies.md) | `KmerAnalyzer.GetKmerFrequencies` | Full k-mer → frequency table. |
| [`kmer_positions`](../../../../docs/mcp/tools/analysis/kmer_positions.md) | `KmerAnalyzer.FindKmerPositions` | 0-based positions of a given k-mer. |
| [`generate_all_kmers`](../../../../docs/mcp/tools/analysis/generate_all_kmers.md) | `KmerAnalyzer.GenerateAllKmers` | Enumerate all k-mers over the alphabet. |
| [`analyze_kmers`](../../../../docs/mcp/tools/analysis/analyze_kmers.md) | `KmerAnalyzer.AnalyzeKmers` | Composite k-mer summary for a sequence. |
| [`kmer_distance`](../../../../docs/mcp/tools/analysis/kmer_distance.md) | `KmerAnalyzer.KmerDistance` | k-mer-profile distance between two sequences. |

## Repeats & low-complexity (Analysis server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`find_repeats`](../../../../docs/mcp/tools/analysis/find_repeats.md) | `GenomicAnalyzer.FindRepeats` | General repeated-region detection. |
| [`find_tandem_repeats`](../../../../docs/mcp/tools/analysis/find_tandem_repeats.md) | `GenomicAnalyzer.FindTandemRepeats` | Tandem repeat arrays. |
| [`find_direct_repeats`](../../../../docs/mcp/tools/analysis/find_direct_repeats.md) | `RepeatFinder.FindDirectRepeats` | Dispersed direct repeats. |
| [`find_inverted_repeats`](../../../../docs/mcp/tools/analysis/find_inverted_repeats.md) | `RepeatFinder.FindInvertedRepeats` | Inverted repeats (branch-inducing). |
| [`find_palindromes`](../../../../docs/mcp/tools/analysis/find_palindromes.md) | `RepeatFinder.FindPalindromes` | DNA palindromes. |
| [`find_microsatellites`](../../../../docs/mcp/tools/analysis/find_microsatellites.md) | `RepeatFinder.FindMicrosatellites` | Microsatellites / SSRs. |
| [`tandem_repeat_summary`](../../../../docs/mcp/tools/analysis/tandem_repeat_summary.md) | `RepeatFinder.GetTandemRepeatSummary` | Compact per-sequence tandem-repeat report. |
| [`find_low_complexity_regions`](../../../../docs/mcp/tools/analysis/find_low_complexity_regions.md) | `SequenceComplexity.FindLowComplexityRegions` | Low-complexity region intervals. |
| [`dust_score`](../../../../docs/mcp/tools/analysis/dust_score.md) | `SequenceComplexity.CalculateDustScore` | DUST low-complexity score. |
| [`mask_low_complexity`](../../../../docs/mcp/tools/analysis/mask_low_complexity.md) | `SequenceComplexity.MaskLowComplexity` | Mask low-complexity before assembly. |

## Repeats — independent suffix-tree path (Core server)

| Tool | Method ID | Purpose |
|---|---|---|
| [`find_longest_repeat`](../../../../docs/mcp/tools/core/find_longest_repeat.md) | `GenomicAnalyzer.FindLongestRepeat` | Longest repeated region in a DNA sequence (with position). |
| [`suffix_tree_lrs`](../../../../docs/mcp/tools/core/suffix_tree_lrs.md) | `SuffixTree.LongestRepeatedSubstring` | Longest repeated substring of a raw text (suffix-tree path — cross-check). |

> **Deferred to [`bio-annotation`](../../bio-annotation/SKILL.md)** (Analysis-server overlap): all
> motif discovery/scanning (`discover_motifs`, `create_pwm`, `scan_with_pwm`, `find_*_motif*`),
> ORFs/genes/promoters, variant calling/annotation/effect, RNA secondary structure, and protein
> features. Those are annotation-type analyses, not assembly — use `bio-annotation`, not this skill.
