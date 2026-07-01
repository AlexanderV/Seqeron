using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using ModelContextProtocol.Server;
using Seqeron.Genomics.Core;
using Seqeron.Mcp.Chromosome.Models;
using SourceCa = Seqeron.Genomics.Chromosome.ChromosomeAnalyzer;
using SourceGa = Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Mcp.Chromosome.Tools;

/// <summary>
/// MCP tools for chromosome operations.
/// </summary>
[McpServerToolType]
public class ChromosomeTools
{
    #region ChromosomeAnalyzer tools

    [McpServerTool(Name = "analyze_karyotype", Title = "Karyotype — Analyze", ReadOnly = true)]
    [Description("Analyze karyotype from chromosome data; reports total/autosome counts, sex chromosomes, ploidy and aneuploidy abnormalities (ISCN nomenclature).")]
    public static SourceCa.Karyotype AnalyzeKaryotype(
        [Description("Chromosomes (name, length, isSexChromosome).")] IReadOnlyList<ChromosomeInput> chromosomes,
        [Description("Expected ploidy level (default 2).")] int expectedPloidyLevel = 2)
    {
        if (chromosomes is null)
            throw new ArgumentException("Chromosomes cannot be null", nameof(chromosomes));
        if (expectedPloidyLevel <= 0)
            throw new ArgumentException("Expected ploidy level must be positive", nameof(expectedPloidyLevel));

        var tuples = chromosomes.Select(c => (c.Name, c.Length, c.IsSexChromosome));
        return SourceCa.AnalyzeKaryotype(tuples, expectedPloidyLevel);
    }

    [McpServerTool(Name = "detect_ploidy", Title = "Ploidy — Detect from depth", ReadOnly = true)]
    [Description("Detect ploidy level from normalized read depth values; returns ploidy and confidence.")]
    public static PloidyResult DetectPloidy(
        [Description("Normalized read depths.")] IReadOnlyList<double> normalizedDepths,
        [Description("Expected diploid depth (default 1.0).")] double expectedDiploidDepth = 1.0)
    {
        var (ploidy, confidence) = SourceCa.DetectPloidy(normalizedDepths, expectedDiploidDepth);
        return new PloidyResult(ploidy, confidence);
    }

    [McpServerTool(Name = "analyze_telomeres", Title = "Telomeres — Analyze ends", ReadOnly = true)]
    [Description("Analyze 5'/3' telomere repeat tracts on a chromosome sequence; reports lengths, repeat purity, and whether critically short.")]
    public static SourceCa.TelomereResult AnalyzeTelomeres(
        [Description("Chromosome name.")] string chromosomeName,
        [Description("Chromosome sequence.")] string sequence,
        [Description("Telomere repeat unit (default TTAGGG, vertebrate).")] string telomereRepeat = "TTAGGG",
        [Description("End-region search length in bp (default 10000).")] int searchLength = 10000,
        [Description("Minimum length to call a telomere present (default 500).")] int minTelomereLength = 500,
        [Description("Critical telomere shortening threshold (default 3000).")] int criticalLength = 3000)
    {
        if (string.IsNullOrEmpty(chromosomeName))
            throw new ArgumentException("Chromosome name cannot be null or empty", nameof(chromosomeName));
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (string.IsNullOrEmpty(telomereRepeat))
            throw new ArgumentException("Telomere repeat cannot be null or empty", nameof(telomereRepeat));
        if (searchLength <= 0)
            throw new ArgumentException("Search length must be positive", nameof(searchLength));

        return SourceCa.AnalyzeTelomeres(chromosomeName, sequence, telomereRepeat, searchLength, minTelomereLength, criticalLength);
    }

    [McpServerTool(Name = "estimate_telomere_length_from_ts_ratio", Title = "Telomeres — Length from qPCR T/S ratio", ReadOnly = true)]
    [Description("Estimate telomere length in bp from a qPCR Telomere/Single-copy gene (T/S) ratio against a reference.")]
    public static TelomereLengthEstimate EstimateTelomereLengthFromTsRatio(
        [Description("Sample T/S ratio.")] double tsRatio,
        [Description("Reference T/S ratio (default 1.0).")] double referenceRatio = 1.0,
        [Description("Telomere length corresponding to reference ratio in bp (default 7000).")] double referenceLength = 7000)
        => new(SourceCa.EstimateTelomereLengthFromTSRatio(tsRatio, referenceRatio, referenceLength));

    [McpServerTool(Name = "analyze_centromere", Title = "Centromere — Analyze", ReadOnly = true)]
    [Description("Locate the centromere region by alpha-satellite-like repeat content and classify its type per Levan et al. (1964).")]
    public static SourceCa.CentromereResult AnalyzeCentromere(
        [Description("Chromosome name.")] string chromosomeName,
        [Description("Chromosome sequence.")] string sequence,
        [Description("Scan window size in bp (default 100000).")] int windowSize = 100000,
        [Description("Minimum alpha-satellite-like content (default 0.3).")] double minAlphaSatelliteContent = 0.3)
    {
        if (string.IsNullOrEmpty(chromosomeName))
            throw new ArgumentException("Chromosome name cannot be null or empty", nameof(chromosomeName));
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (windowSize <= 0)
            throw new ArgumentException("Window size must be positive", nameof(windowSize));
        if (minAlphaSatelliteContent is < 0 or > 1)
            throw new ArgumentException("Minimum alpha-satellite content must be in [0, 1]", nameof(minAlphaSatelliteContent));

        return SourceCa.AnalyzeCentromere(chromosomeName, sequence, windowSize, minAlphaSatelliteContent);
    }

    [McpServerTool(Name = "predict_g_bands", Title = "Cytogenetic — Predict G-bands", ReadOnly = true)]
    [Description("Predict cytogenetic G-band pattern from sequence GC content (gpos100/gpos50/gneg, simplified).")]
    public static CytogeneticBandsResult PredictGBands(
        [Description("Chromosome name.")] string chromosomeName,
        [Description("Chromosome sequence.")] string sequence,
        [Description("Band size in bp (default 5,000,000).")] int bandSize = 5000000,
        [Description("GC threshold for dark band (default 0.37).")] double darkBandGcThreshold = 0.37,
        [Description("GC threshold for light band (default 0.45).")] double lightBandGcThreshold = 0.45)
        => new(SourceCa.PredictGBands(chromosomeName, sequence, bandSize, darkBandGcThreshold, lightBandGcThreshold).ToList());

    [McpServerTool(Name = "find_heterochromatin_regions", Title = "Heterochromatin — Find regions", ReadOnly = true)]
    [Description("Identify heterochromatin regions by k-mer repeat content; classifies as Telomeric, Centromeric, or Constitutive.")]
    public static HeterochromatinRegionsResult FindHeterochromatinRegions(
        [Description("Chromosome sequence.")] string sequence,
        [Description("Window size in bp (default 100000).")] int windowSize = 100000,
        [Description("Minimum repeat content fraction (default 0.5).")] double minRepeatContent = 0.5)
        => new(SourceCa.FindHeterochromatinRegions(sequence, windowSize, minRepeatContent)
            .Select(r => new HeterochromatinRegion(r.Start, r.End, r.Type))
            .ToList());

    [McpServerTool(Name = "find_synteny_blocks", Title = "Synteny — From orthologs", ReadOnly = true)]
    [Description("Identify collinear synteny blocks between two genomes from a list of ortholog gene pairs.")]
    public static SyntenyBlocksResult FindSyntenyBlocks(
        [Description("Ortholog gene pairs across the two genomes.")] IReadOnlyList<OrthologPair> orthologPairs,
        [Description("Minimum genes per block (default 3).")] int minGenes = 3,
        [Description("Maximum gap between consecutive genes in Mb (default 10).")] int maxGap = 10)
    {
        var tuples = orthologPairs.Select(p =>
            (p.Chr1, p.Start1, p.End1, p.Gene1, p.Chr2, p.Start2, p.End2, p.Gene2));
        return new(SourceCa.FindSyntenyBlocks(tuples, minGenes, maxGap).ToList());
    }

    [McpServerTool(Name = "detect_rearrangements", Title = "Rearrangements — Detect from synteny", ReadOnly = true)]
    [Description("Detect chromosomal rearrangements (inversions, translocations, deletions, duplications) from synteny blocks.")]
    public static RearrangementsResult DetectRearrangements(
        [Description("Synteny blocks (typically from find_synteny_blocks).")] IReadOnlyList<SourceCa.SyntenyBlock> syntenyBlocks)
        => new(SourceCa.DetectRearrangements(syntenyBlocks).ToList());

    [McpServerTool(Name = "detect_aneuploidy", Title = "Aneuploidy — Detect from depth", ReadOnly = true)]
    [Description("Detect copy-number states from binned read-depth data; reports per-bin copy number, log2 ratio and confidence.")]
    public static CopyNumberStatesResult DetectAneuploidy(
        [Description("Read-depth samples per chromosome position.")] IReadOnlyList<DepthSample> depthData,
        [Description("Genome-wide median depth.")] double medianDepth,
        [Description("Bin size in bp (default 1,000,000).")] int binSize = 1000000)
    {
        var tuples = depthData.Select(d => (d.Chromosome, d.Position, d.Depth));
        return new(SourceCa.DetectAneuploidy(tuples, medianDepth, binSize).ToList());
    }

    [McpServerTool(Name = "identify_whole_chromosome_aneuploidy", Title = "Aneuploidy — Whole-chromosome", ReadOnly = true)]
    [Description("Identify whole-chromosome aneuploidies (Monosomy/Trisomy/Tetrasomy/etc.) from per-bin copy-number states.")]
    public static WholeChromosomeAneuploidiesResult IdentifyWholeChromosomeAneuploidy(
        [Description("Per-bin copy-number states.")] IReadOnlyList<SourceCa.CopyNumberState> copyNumberStates,
        [Description("Minimum dominant fraction (default 0.8).")] double minFraction = 0.8)
        => new(SourceCa.IdentifyWholeChromosomeAneuploidy(copyNumberStates, minFraction)
            .Select(a => new WholeChromosomeAneuploidy(a.Chromosome, a.CopyNumber, a.Type))
            .ToList());

    [McpServerTool(Name = "arm_ratio", Title = "Chromosome — p/q arm ratio", ReadOnly = true)]
    [Description("Compute chromosome arm ratio (p/q) from centromere position and total length.")]
    public static ArmRatioResult ArmRatio(
        [Description("Centromere position in bp.")] int centromerePosition,
        [Description("Total chromosome length in bp.")] int chromosomeLength)
        => new(SourceCa.CalculateArmRatio(centromerePosition, chromosomeLength));

    [McpServerTool(Name = "classify_chromosome_by_arm_ratio", Title = "Chromosome — Classify by arm ratio", ReadOnly = true)]
    [Description("Classify a chromosome (Metacentric/Submetacentric/Acrocentric/Telocentric) from its p/q arm ratio per Levan et al. (1964).")]
    public static ChromosomeClassification ClassifyChromosomeByArmRatio(
        [Description("Arm ratio (p/q).")] double armRatio)
        => new(SourceCa.ClassifyChromosomeByArmRatio(armRatio));

    [McpServerTool(Name = "estimate_cell_divisions_from_telomere_length", Title = "Telomeres — Estimate cell divisions", ReadOnly = true)]
    [Description("Estimate the number of cell divisions from current telomere length given birth length and per-division loss rate.")]
    public static CellDivisionsEstimate EstimateCellDivisionsFromTelomereLength(
        [Description("Current telomere length in bp.")] int currentLength,
        [Description("Birth telomere length in bp (default 15000).")] int birthLength = 15000,
        [Description("Telomere loss per division in bp (default 50).")] int lossPerDivision = 50)
        => new(SourceCa.EstimateCellDivisionsFromTelomereLength(currentLength, birthLength, lossPerDivision));

    #endregion

    #region GenomeAssemblyAnalyzer tools

    [McpServerTool(Name = "assembly_statistics", Title = "Assembly — Statistics (N50, L50, GC, gaps)", ReadOnly = true)]
    [Description("Compute comprehensive assembly statistics: total length, N50/L50/N90/L90, largest/smallest, mean/median, GC, and gap metrics.")]
    public static SourceGa.AssemblyStatistics AssemblyStatistics(
        [Description("Assembled sequences (id, sequence).")] IReadOnlyList<NamedSequence> sequences)
        => SourceGa.CalculateStatistics(sequences.Select(s => (s.Id, s.Sequence)));

    [McpServerTool(Name = "nx_statistics", Title = "Assembly — Nx/Lx (single threshold)", ReadOnly = true)]
    [Description("Compute Nx and Lx for a single threshold given pre-sorted (descending) sequence lengths and total length.")]
    public static SourceGa.NxStatistics NxStatistics(
        [Description("Lengths sorted descending.")] IReadOnlyList<int> sortedLengths,
        [Description("Sum of all lengths.")] long totalLength,
        [Description("Threshold percent (e.g. 50 for N50).")] int threshold)
        => SourceGa.CalculateNx(sortedLengths, totalLength, threshold);

    [McpServerTool(Name = "nx_curve", Title = "Assembly — Nx curve", ReadOnly = true)]
    [Description("Compute Nx/Lx for multiple thresholds (default 10..90 step 10).")]
    public static NxStatisticsListResult NxCurve(
        [Description("Sequence lengths.")] IReadOnlyList<int> lengths,
        [Description("Thresholds (empty → 10..90 step 10).")] int[]? thresholds = null)
        => new(SourceGa.CalculateNxCurve(lengths, thresholds ?? System.Array.Empty<int>()).ToList());

    [McpServerTool(Name = "au_n", Title = "Assembly — auN (area under Nx curve)", ReadOnly = true)]
    [Description("Compute auN (area under Nx curve) — a length-weighted contiguity metric robust to outliers.")]
    public static AuNResult AuN(
        [Description("Sequence lengths.")] IReadOnlyList<int> lengths)
        => new(SourceGa.CalculateAuN(lengths));

    [McpServerTool(Name = "find_gaps", Title = "Assembly — Find gaps (N runs)", ReadOnly = true)]
    [Description("Find gaps (N-runs) in assembled sequences; returns position, length and length-class type.")]
    public static GapsResult FindGaps(
        [Description("Assembled sequences (id, sequence).")] IReadOnlyList<NamedSequence> sequences,
        [Description("Minimum gap length (default 1).")] int minGapLength = 1)
        => new(SourceGa.FindGaps(sequences.Select(s => (s.Id, s.Sequence)), minGapLength).ToList());

    [McpServerTool(Name = "gap_distribution", Title = "Assembly — Gap distribution", ReadOnly = true)]
    [Description("Summarize a list of gaps: count, mean/median/max length, and counts by gap-length type.")]
    public static GapDistributionResult GapDistribution(
        [Description("Gaps (typically from find_gaps).")] IReadOnlyList<SourceGa.GapInfo> gaps)
    {
        var (count, mean, median, max, types) = SourceGa.AnalyzeGapDistribution(gaps);
        return new GapDistributionResult(count, mean, median, max, types.ToDictionary(kv => kv.Key, kv => kv.Value));
    }

    [McpServerTool(Name = "analyze_scaffolds", Title = "Assembly — Scaffold structure", ReadOnly = true)]
    [Description("Decompose scaffolds into contigs and gaps; returns per-scaffold contig list, gap list, and length totals.")]
    public static ScaffoldStructuresResult AnalyzeScaffolds(
        [Description("Scaffolds (id, sequence).")] IReadOnlyList<NamedSequence> scaffolds,
        [Description("Minimum gap length to split contigs (default 10).")] int minGapLength = 10)
    {
        if (scaffolds is null)
            throw new ArgumentException("Scaffolds cannot be null", nameof(scaffolds));
        if (minGapLength <= 0)
            throw new ArgumentException("Minimum gap length must be positive", nameof(minGapLength));

        return new(SourceGa.AnalyzeScaffolds(scaffolds.Select(s => (s.Id, s.Sequence)), minGapLength)
            .Select(ss => new ScaffoldStructureItem(
                ss.ScaffoldId,
                ss.Contigs.Select(c => new ScaffoldContig(c.ContigId, c.Start, c.End)).ToList(),
                ss.Gaps,
                ss.TotalLength,
                ss.ContigLength,
                ss.GapLength))
            .ToList());
    }

    [McpServerTool(Name = "extract_contigs", Title = "Assembly — Extract contigs", ReadOnly = true)]
    [Description("Extract contigs (gap-free runs) from scaffolds with a minimum length filter.")]
    public static NamedSequencesResult ExtractContigs(
        [Description("Scaffolds (id, sequence).")] IReadOnlyList<NamedSequence> scaffolds,
        [Description("Minimum contig length (default 200).")] int minContigLength = 200)
        => new(SourceGa.ExtractContigs(scaffolds.Select(s => (s.Id, s.Sequence)), minContigLength)
            .Select(c => new NamedSequence(c.Id, c.Sequence))
            .ToList());

    [McpServerTool(Name = "assess_completeness", Title = "Assembly — BUSCO-like completeness", ReadOnly = true)]
    [Description("Assess assembly completeness by aligning marker genes (BUSCO-like, k-mer based); reports complete/single/duplicated/fragmented/missing.")]
    public static SourceGa.CompletenessResult AssessCompleteness(
        [Description("Assembled sequences (id, sequence).")] IReadOnlyList<NamedSequence> assembly,
        [Description("Marker genes (geneId, sequence).")] IReadOnlyList<NamedSequence> markerGenes,
        [Description("Identity threshold (default 0.9).")] double identityThreshold = 0.9,
        [Description("Coverage threshold (default 0.9).")] double coverageThreshold = 0.9)
        => SourceGa.AssessCompleteness(
            assembly.Select(s => (s.Id, s.Sequence)),
            markerGenes.Select(s => (s.Id, s.Sequence)),
            identityThreshold,
            coverageThreshold);

    [McpServerTool(Name = "estimate_completeness_from_kmers", Title = "Assembly — Completeness from k-mer spectrum", ReadOnly = true)]
    [Description("Estimate genome completeness, error rate and genome size from a k-mer count spectrum.")]
    public static KmerCompletenessResult EstimateCompletenessFromKmers(
        [Description("K-mer spectrum (kmer, count).")] IReadOnlyList<KmerCount> kmerSpectrum,
        [Description("Expected coverage (0 = auto-detect peak).")] int expectedCoverage = 0)
    {
        var (completeness, errorRate, size) =
            SourceGa.EstimateCompletenessFromKmers(kmerSpectrum.Select(k => (k.Kmer, k.Count)), expectedCoverage);
        return new KmerCompletenessResult(completeness, errorRate, size);
    }

    [McpServerTool(Name = "find_repetitive_regions", Title = "Assembly — Find repetitive regions", ReadOnly = true)]
    [Description("Identify repetitive regions across assembled sequences using k-mer copy-number frequency.")]
    public static RepetitiveRegionsResult FindRepetitiveRegions(
        [Description("Sequences (id, sequence).")] IReadOnlyList<NamedSequence> sequences,
        [Description("K-mer size (default 15).")] int kmerSize = 15,
        [Description("Minimum k-mer copies (default 3).")] int minCopies = 3,
        [Description("Window size (default 100).")] int windowSize = 100)
        => new(SourceGa.FindRepetitiveRegions(
                sequences.Select(s => (s.Id, s.Sequence)), kmerSize, minCopies, windowSize)
            .Select(r => new RepetitiveRegion(r.SequenceId, r.Start, r.End, r.Copies))
            .ToList());

    [McpServerTool(Name = "find_tandem_repeats", Title = "Assembly — Find tandem repeats", ReadOnly = true)]
    [Description("Identify tandem repeats; returns repeat unit, copy number and purity per occurrence.")]
    public static TandemRepeatsResult FindTandemRepeats(
        [Description("Sequences (id, sequence).")] IReadOnlyList<NamedSequence> sequences,
        [Description("Minimum unit length (default 2).")] int minUnitLength = 2,
        [Description("Maximum unit length (default 50).")] int maxUnitLength = 50,
        [Description("Minimum copies (default 3).")] int minCopies = 3)
        => new(SourceGa.FindTandemRepeats(
                sequences.Select(s => (s.Id, s.Sequence)), minUnitLength, maxUnitLength, minCopies)
            .Select(r => new TandemRepeat(r.SequenceId, r.Start, r.End, r.Unit, r.Copies, r.Purity))
            .ToList());

    [McpServerTool(Name = "repeat_content", Title = "Assembly — Repeat content from annotations", ReadOnly = true)]
    [Description("Compute total repeat length, repeat percentage and per-class lengths from repeat annotations.")]
    public static RepeatContentResult RepeatContent(
        [Description("Repeat annotations.")] IReadOnlyList<SourceGa.RepeatAnnotation> repeats,
        [Description("Genome length in bp.")] long genomeLength)
    {
        var (total, percent, classLens) = SourceGa.CalculateRepeatContent(repeats, genomeLength);
        return new RepeatContentResult(total, percent, classLens.ToDictionary(kv => kv.Key, kv => kv.Value));
    }

    [McpServerTool(Name = "compare_assemblies", Title = "Assembly — Compare two assemblies", ReadOnly = true)]
    [Description("Compare two assemblies by shared k-mer content; reports aligned fractions and an identity proxy.")]
    public static SourceGa.AssemblyComparison CompareAssemblies(
        [Description("Assembly 1 sequences.")] IReadOnlyList<NamedSequence> assembly1,
        [Description("Assembly 2 sequences.")] IReadOnlyList<NamedSequence> assembly2,
        [Description("Name for assembly 1 (default Assembly1).")] string name1 = "Assembly1",
        [Description("Name for assembly 2 (default Assembly2).")] string name2 = "Assembly2",
        [Description("K-mer size (default 21).")] int kmerSize = 21)
        => SourceGa.CompareAssemblies(
            assembly1.Select(s => (s.Id, s.Sequence)),
            assembly2.Select(s => (s.Id, s.Sequence)),
            name1, name2, kmerSize);

    [McpServerTool(Name = "find_syntenic_blocks_assemblies", Title = "Assembly — Syntenic blocks (k-mer based)", ReadOnly = true)]
    [Description("Find syntenic blocks between two assemblies via k-mer anchor clustering; flags inverted blocks.")]
    public static SyntenicBlocksResult FindSyntenicBlocksAssemblies(
        [Description("Assembly 1 sequences.")] IReadOnlyList<NamedSequence> assembly1,
        [Description("Assembly 2 sequences.")] IReadOnlyList<NamedSequence> assembly2,
        [Description("Minimum block size (default 1000).")] int minBlockSize = 1000,
        [Description("K-mer size (default 21).")] int kmerSize = 21)
        => new(SourceGa.FindSyntenicBlocks(
                assembly1.Select(s => (s.Id, s.Sequence)),
                assembly2.Select(s => (s.Id, s.Sequence)),
                minBlockSize, kmerSize)
            .Select(b => new SyntenicBlockItem(b.Seq1, b.Start1, b.End1, b.Seq2, b.Start2, b.End2, b.IsInverted))
            .ToList());

    [McpServerTool(Name = "local_quality", Title = "Assembly — Local quality windows", ReadOnly = true)]
    [Description("Compute per-window local quality metrics (GC content, N count, linguistic complexity).")]
    public static LocalQualityResult LocalQuality(
        [Description("Sequences (id, sequence).")] IReadOnlyList<NamedSequence> sequences,
        [Description("Window size (default 1000).")] int windowSize = 1000)
        => new(SourceGa.CalculateLocalQuality(sequences.Select(s => (s.Id, s.Sequence)), windowSize)
            .Select(q => new LocalQualityWindow(q.SequenceId, q.Position, q.WindowSize, q.GcContent, q.NCount, q.Complexity))
            .ToList());

    [McpServerTool(Name = "find_suspicious_regions", Title = "Assembly — Suspicious regions", ReadOnly = true)]
    [Description("Flag potentially misassembled regions by GC deviation, low complexity, and high N content.")]
    public static SuspiciousRegionsResult FindSuspiciousRegions(
        [Description("Sequences (id, sequence).")] IReadOnlyList<NamedSequence> sequences,
        [Description("Allowed GC deviation from global GC (default 0.15).")] double gcDeviation = 0.15,
        [Description("Minimum linguistic complexity (default 0.3).")] double minComplexity = 0.3)
        => new(SourceGa.FindSuspiciousRegions(sequences.Select(s => (s.Id, s.Sequence)), gcDeviation, minComplexity)
            .Select(r => new SuspiciousRegion(r.SequenceId, r.Start, r.End, r.Reason, r.Score))
            .ToList());

    [McpServerTool(Name = "length_distribution", Title = "Assembly — Length distribution", ReadOnly = true)]
    [Description("Bucket sequence lengths into bins (default: 100..1,000,000 powers).")]
    public static LengthDistributionResult LengthDistribution(
        [Description("Sequence lengths.")] IReadOnlyList<int> lengths,
        [Description("Bin upper bounds (empty → defaults).")] int[]? bins = null)
    {
        var dist = SourceGa.CalculateLengthDistribution(lengths, bins ?? System.Array.Empty<int>());
        return new LengthDistributionResult(dist.ToDictionary(kv => kv.Key, kv => kv.Value));
    }

    #endregion
}
