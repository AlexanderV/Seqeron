using SourceCa = Seqeron.Genomics.Chromosome.ChromosomeAnalyzer;
using SourceGa = Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Mcp.Chromosome.Models;

#region Input DTOs

/// <summary>
/// Named DNA/RNA sequence (id + sequence). Replaces value-tuple inputs.
/// </summary>
public record NamedSequence(string Id, string Sequence);

/// <summary>
/// Chromosome descriptor used for karyotype analysis.
/// </summary>
public record ChromosomeInput(string Name, long Length, bool IsSexChromosome);

/// <summary>
/// Ortholog gene pair between two species used for synteny block detection.
/// </summary>
public record OrthologPair(
    string Chr1, int Start1, int End1, string Gene1,
    string Chr2, int Start2, int End2, string Gene2);

/// <summary>
/// Read-depth sample at a single chromosome position.
/// </summary>
public record DepthSample(string Chromosome, int Position, double Depth);

/// <summary>
/// K-mer with its observed count from a k-mer spectrum.
/// </summary>
public record KmerCount(string Kmer, int Count);

#endregion

#region Output DTOs (wrappers for tuples / scalars / nested-tuple records)

public record PloidyResult(int PloidyLevel, double Confidence);

public record TelomereLengthEstimate(double TelomereLength);

public record HeterochromatinRegion(int Start, int End, string Type);

public record HeterochromatinRegionsResult(IReadOnlyList<HeterochromatinRegion> Items);

public record SyntenyBlocksResult(IReadOnlyList<SourceCa.SyntenyBlock> Items);

public record CytogeneticBandsResult(IReadOnlyList<SourceCa.CytogeneticBand> Items);

public record RearrangementsResult(IReadOnlyList<SourceCa.ChromosomalRearrangement> Items);

public record CopyNumberStatesResult(IReadOnlyList<SourceCa.CopyNumberState> Items);

public record WholeChromosomeAneuploidy(string Chromosome, int CopyNumber, string Type);

public record WholeChromosomeAneuploidiesResult(IReadOnlyList<WholeChromosomeAneuploidy> Items);

public record ArmRatioResult(double ArmRatio);

public record ChromosomeClassification(string Classification);

public record CellDivisionsEstimate(double CellDivisions);

public record NxStatisticsListResult(IReadOnlyList<SourceGa.NxStatistics> Items);

public record AuNResult(double AuN);

public record GapsResult(IReadOnlyList<SourceGa.GapInfo> Items);

public record GapDistributionResult(
    int Count,
    double MeanLength,
    double MedianLength,
    int MaxLength,
    IReadOnlyDictionary<string, int> TypeCounts);

/// <summary>
/// Contig entry inside a scaffold (replaces value-tuple in source ScaffoldStructure).
/// </summary>
public record ScaffoldContig(string ContigId, int Start, int End);

/// <summary>
/// Scaffold structure with named contigs (MCP-serializable variant of source ScaffoldStructure).
/// </summary>
public record ScaffoldStructureItem(
    string ScaffoldId,
    IReadOnlyList<ScaffoldContig> Contigs,
    IReadOnlyList<SourceGa.GapInfo> Gaps,
    int TotalLength,
    int ContigLength,
    int GapLength);

public record ScaffoldStructuresResult(IReadOnlyList<ScaffoldStructureItem> Items);

public record NamedSequencesResult(IReadOnlyList<NamedSequence> Items);

public record KmerCompletenessResult(
    double Completeness,
    double ErrorRate,
    long EstimatedGenomeSize);

public record RepetitiveRegion(string SequenceId, int Start, int End, int Copies);

public record RepetitiveRegionsResult(IReadOnlyList<RepetitiveRegion> Items);

public record TandemRepeat(
    string SequenceId,
    int Start,
    int End,
    string Unit,
    int Copies,
    double Purity);

public record TandemRepeatsResult(IReadOnlyList<TandemRepeat> Items);

public record RepeatContentResult(
    long TotalRepeatLength,
    double RepeatPercentage,
    IReadOnlyDictionary<string, long> RepeatClassLengths);

public record SyntenicBlockItem(
    string Seq1,
    int Start1,
    int End1,
    string Seq2,
    int Start2,
    int End2,
    bool IsInverted);

public record SyntenicBlocksResult(IReadOnlyList<SyntenicBlockItem> Items);

public record LocalQualityWindow(
    string SequenceId,
    int Position,
    int WindowSize,
    double GcContent,
    int NCount,
    double Complexity);

public record LocalQualityResult(IReadOnlyList<LocalQualityWindow> Items);

public record SuspiciousRegion(
    string SequenceId,
    int Start,
    int End,
    string Reason,
    double Score);

public record SuspiciousRegionsResult(IReadOnlyList<SuspiciousRegion> Items);

public record LengthDistributionResult(IReadOnlyDictionary<string, int> Distribution);

#endregion
