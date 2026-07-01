# bio-chromosome — tool map (all 32)

Server: **chromosome**. Two backing classes: `ChromosomeAnalyzer.*` (cytogenetics / SV / telomere /
ploidy) and `GenomeAssemblyAnalyzer.*` (assembly QC / contiguity / gaps / completeness).
Verify schemas in `docs/mcp/tools/chromosome/<tool>.md`. GC-skew & replication-origin are on the
**analysis** server (see `bio-annotation`) — listed at the bottom for convenience.

> Coordinates: gap / band / region outputs are typically **1-based inclusive `[start,end]`** as
> documented per tool. Depth-bin positions are `position / binSize`. Always confirm in the tool doc.

## Karyotype / ploidy / aneuploidy

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `analyze_karyotype` | `ChromosomeAnalyzer.AnalyzeKaryotype` | Group autosomes by base name, flag ISCN aneuploidy vs `expectedPloidyLevel` from descriptors. | `analyze_karyotype.md` |
| `detect_aneuploidy` | `ChromosomeAnalyzer.DetectAneuploidy` | Per-bin `logRatio`/`copyNumber`/confidence from depth samples + genome median. | `detect_aneuploidy.md` |
| `identify_whole_chromosome_aneuploidy` | `ChromosomeAnalyzer.IdentifyWholeChromosomeAneuploidy` | Dominant per-chromosome CN ≠ 2 covering ≥ `minFraction` → ISCN aneuploidy. | `identify_whole_chromosome_aneuploidy.md` |
| `detect_ploidy` | `ChromosomeAnalyzer.DetectPloidy` | Genome-wide `ploidy = round(median/expectedDiploidDepth × 2)`, clamped [1,8]. | `detect_ploidy.md` |

## Centromere / arm / bands / heterochromatin

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `analyze_centromere` | `ChromosomeAnalyzer.AnalyzeCentromere` | Alpha-satellite window scan → centromere `start`/`end` + Levan type. | `analyze_centromere.md` |
| `arm_ratio` | `ChromosomeAnalyzer.CalculateArmRatio` | `p=centromerePos`, `q=len-pos`, `armRatio=q/p`. | `arm_ratio.md` |
| `classify_chromosome_by_arm_ratio` | `ChromosomeAnalyzer.ClassifyChromosomeByArmRatio` | long/short ratio → Metacentric…Acrocentric/Telocentric. | `classify_chromosome_by_arm_ratio.md` |
| `predict_g_bands` | `ChromosomeAnalyzer.PredictGBands` | Windowed GC → gpos100 (dark) / gpos50 / gneg cytobands. | `predict_g_bands.md` |
| `find_heterochromatin_regions` | `ChromosomeAnalyzer.FindHeterochromatinRegions` | Merge repeat-rich windows → Telomeric/Centromeric/Constitutive. | `find_heterochromatin_regions.md` |

## Telomere

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `analyze_telomeres` | `ChromosomeAnalyzer.AnalyzeTelomeres` | Scan both ends for TTAGGG/CCCTAA tracts → present? length? critically short? | `analyze_telomeres.md` |
| `estimate_telomere_length_from_ts_ratio` | `ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio` | qPCR T/S ratio → length via reference scaling. | `estimate_telomere_length_from_ts_ratio.md` |
| `estimate_cell_divisions_from_telomere_length` | `ChromosomeAnalyzer.EstimateCellDivisionsFromTelomereLength` | Replicative divisions from telomere shortening. | `estimate_cell_divisions_from_telomere_length.md` |

## Synteny / structural rearrangements

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `find_synteny_blocks` | `ChromosomeAnalyzer.FindSyntenyBlocks` | Ortholog pairs → collinear blocks (≥ `minGenes`, `maxGap` Mb) + strand. | `find_synteny_blocks.md` |
| `detect_rearrangements` | `ChromosomeAnalyzer.DetectRearrangements` | Blocks → inversion / translocation / deletion / duplication + breakpoints. | `detect_rearrangements.md` |
| `find_syntenic_blocks_assemblies` | `GenomeAssemblyAnalyzer.FindSyntenicBlocks` | K-mer-anchored collinear blocks between two assemblies (≥ `minBlockSize`), orientation flags. | `find_syntenic_blocks_assemblies.md` |

## Assembly QC — contiguity / gaps / completeness

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `assembly_statistics` | `GenomeAssemblyAnalyzer.CalculateStatistics` | N50/L50, N90/L90, GC, largest/mean/median, gap counts/lengths. | `assembly_statistics.md` |
| `nx_statistics` | `GenomeAssemblyAnalyzer.CalculateNx` | Nx/Lx at one `threshold%` from **descending-sorted lengths + totalLength**. | `nx_statistics.md` |
| `nx_curve` | `GenomeAssemblyAnalyzer.CalculateNxCurve` | Nx/Lx across thresholds (default deciles 10…90). | `nx_curve.md` |
| `au_n` | `GenomeAssemblyAnalyzer.CalculateAuN` | Threshold-free contiguity `auN = Σlᵢ²/Σlᵢ`. | `au_n.md` |
| `length_distribution` | `GenomeAssemblyAnalyzer.CalculateLengthDistribution` | Bucket lengths into size bins. | `length_distribution.md` |
| `find_gaps` | `GenomeAssemblyAnalyzer.FindGaps` | Maximal N-runs ≥ `minGapLength` → `[start,end]` + length class. | `find_gaps.md` |
| `gap_distribution` | `GenomeAssemblyAnalyzer.AnalyzeGapDistribution` | Gap count, mean/median/max, per-type counts. | `gap_distribution.md` |
| `extract_contigs` | `GenomeAssemblyAnalyzer.ExtractContigs` | Split scaffolds at N-runs → `{id}_contig{n}` ≥ `minContigLength`. | `extract_contigs.md` |
| `analyze_scaffolds` | `GenomeAssemblyAnalyzer.AnalyzeScaffolds` | Scaffold→contig breakdown with gap runs ≥ `minGapLength`. | `analyze_scaffolds.md` |
| `assess_completeness` | `GenomeAssemblyAnalyzer.AssessCompleteness` | Marker-gene (BUSCO-like) complete/duplicated/fragmented/missing. | `assess_completeness.md` |
| `estimate_completeness_from_kmers` | `GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers` | K-mer spectrum → completeness / coverage peak (auto if `expectedCoverage=0`). | `estimate_completeness_from_kmers.md` |
| `local_quality` | `GenomeAssemblyAnalyzer.CalculateLocalQuality` | Per-window GC / N-count / 4-mer linguistic complexity. | `local_quality.md` |
| `find_suspicious_regions` | `GenomeAssemblyAnalyzer.FindSuspiciousRegions` | Windows with GC deviating > `gcDeviation` from global or low complexity. | `find_suspicious_regions.md` |
| `find_tandem_repeats` | `GenomeAssemblyAnalyzer.FindTandemRepeats` | Tandem arrays (unit `minUnitLength`–`maxUnitLength`, ≥ `minCopies`) + purity. | `find_tandem_repeats.md` |
| `find_repetitive_regions` | `GenomeAssemblyAnalyzer.FindRepetitiveRegions` | K-mers recurring ≥ `minCopies`, merged within `windowSize`. | `find_repetitive_regions.md` |
| `repeat_content` | `GenomeAssemblyAnalyzer.CalculateRepeatContent` | Total repeat bp × 100 / genomeLength, grouped by `repeatClass`. | `repeat_content.md` |
| `compare_assemblies` | `GenomeAssemblyAnalyzer.CompareAssemblies` | Shared-k-mer fractions + identity proxy (structural counts reported 0). | `compare_assemblies.md` |

## Off-server (analysis) — GC-skew & replication origin (use via `bio-annotation`)

| Tool | Method ID | Purpose | Doc |
|---|---|---|---|
| `gc_skew` | `GcSkewCalculator.CalculateGcSkew` | (G−C)/(G+C) skew. | `docs/mcp/tools/analysis/gc_skew.md` |
| `cumulative_gc_skew` | `GcSkewCalculator.CalculateCumulativeGcSkew` | Running cumulative skew (min→origin, max→terminus). | `docs/mcp/tools/analysis/cumulative_gc_skew.md` |
| `windowed_gc_skew` | `GcSkewCalculator.CalculateWindowedGcSkew` | Sliding-window skew profile. | `docs/mcp/tools/analysis/windowed_gc_skew.md` |
| `predict_replication_origin` | `GcSkewCalculator.PredictReplicationOrigin` | Origin/terminus from cumulative-skew extrema. | `docs/mcp/tools/analysis/predict_replication_origin.md` |
