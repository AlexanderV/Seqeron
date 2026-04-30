using System.Collections.Generic;
using System.Text.Json.Serialization;
using Seqeron.Genomics.Chromosome;
using Seqeron.Mcp.Chromosome.Models;

namespace Seqeron.Mcp.Chromosome;

[JsonSourceGenerationOptions(WriteIndented = false)]
[JsonSerializable(typeof(string))]
[JsonSerializable(typeof(int))]
[JsonSerializable(typeof(long))]
[JsonSerializable(typeof(double))]
[JsonSerializable(typeof(bool))]
[JsonSerializable(typeof(int[]))]
[JsonSerializable(typeof(double[]))]

// ===== Input DTOs =====
[JsonSerializable(typeof(NamedSequence))]
[JsonSerializable(typeof(IReadOnlyList<NamedSequence>))]
[JsonSerializable(typeof(ChromosomeInput))]
[JsonSerializable(typeof(IReadOnlyList<ChromosomeInput>))]
[JsonSerializable(typeof(OrthologPair))]
[JsonSerializable(typeof(IReadOnlyList<OrthologPair>))]
[JsonSerializable(typeof(DepthSample))]
[JsonSerializable(typeof(IReadOnlyList<DepthSample>))]
[JsonSerializable(typeof(KmerCount))]
[JsonSerializable(typeof(IReadOnlyList<KmerCount>))]
[JsonSerializable(typeof(IReadOnlyList<int>))]
[JsonSerializable(typeof(IReadOnlyList<double>))]

// ===== ChromosomeAnalyzer source records =====
[JsonSerializable(typeof(ChromosomeAnalyzer.Karyotype))]
[JsonSerializable(typeof(ChromosomeAnalyzer.TelomereResult))]
[JsonSerializable(typeof(ChromosomeAnalyzer.CentromereResult))]
[JsonSerializable(typeof(ChromosomeAnalyzer.CytogeneticBand))]
[JsonSerializable(typeof(ChromosomeAnalyzer.SyntenyBlock))]
[JsonSerializable(typeof(IReadOnlyList<ChromosomeAnalyzer.SyntenyBlock>))]
[JsonSerializable(typeof(ChromosomeAnalyzer.ChromosomalRearrangement))]
[JsonSerializable(typeof(ChromosomeAnalyzer.CopyNumberState))]
[JsonSerializable(typeof(IReadOnlyList<ChromosomeAnalyzer.CopyNumberState>))]

// ===== GenomeAssemblyAnalyzer source records =====
[JsonSerializable(typeof(GenomeAssemblyAnalyzer.AssemblyStatistics))]
[JsonSerializable(typeof(GenomeAssemblyAnalyzer.NxStatistics))]
[JsonSerializable(typeof(GenomeAssemblyAnalyzer.GapInfo))]
[JsonSerializable(typeof(IReadOnlyList<GenomeAssemblyAnalyzer.GapInfo>))]
[JsonSerializable(typeof(GenomeAssemblyAnalyzer.RepeatAnnotation))]
[JsonSerializable(typeof(IReadOnlyList<GenomeAssemblyAnalyzer.RepeatAnnotation>))]
[JsonSerializable(typeof(GenomeAssemblyAnalyzer.CompletenessResult))]
[JsonSerializable(typeof(GenomeAssemblyAnalyzer.AssemblyComparison))]

// ===== Output wrapper DTOs =====
[JsonSerializable(typeof(PloidyResult))]
[JsonSerializable(typeof(TelomereLengthEstimate))]
[JsonSerializable(typeof(HeterochromatinRegion))]
[JsonSerializable(typeof(HeterochromatinRegionsResult))]
[JsonSerializable(typeof(SyntenyBlocksResult))]
[JsonSerializable(typeof(CytogeneticBandsResult))]
[JsonSerializable(typeof(RearrangementsResult))]
[JsonSerializable(typeof(CopyNumberStatesResult))]
[JsonSerializable(typeof(WholeChromosomeAneuploidy))]
[JsonSerializable(typeof(WholeChromosomeAneuploidiesResult))]
[JsonSerializable(typeof(ArmRatioResult))]
[JsonSerializable(typeof(ChromosomeClassification))]
[JsonSerializable(typeof(CellDivisionsEstimate))]
[JsonSerializable(typeof(NxStatisticsListResult))]
[JsonSerializable(typeof(AuNResult))]
[JsonSerializable(typeof(GapsResult))]
[JsonSerializable(typeof(GapDistributionResult))]
[JsonSerializable(typeof(ScaffoldContig))]
[JsonSerializable(typeof(ScaffoldStructureItem))]
[JsonSerializable(typeof(ScaffoldStructuresResult))]
[JsonSerializable(typeof(NamedSequencesResult))]
[JsonSerializable(typeof(KmerCompletenessResult))]
[JsonSerializable(typeof(RepetitiveRegion))]
[JsonSerializable(typeof(RepetitiveRegionsResult))]
[JsonSerializable(typeof(TandemRepeat))]
[JsonSerializable(typeof(TandemRepeatsResult))]
[JsonSerializable(typeof(RepeatContentResult))]
[JsonSerializable(typeof(SyntenicBlockItem))]
[JsonSerializable(typeof(SyntenicBlocksResult))]
[JsonSerializable(typeof(LocalQualityWindow))]
[JsonSerializable(typeof(LocalQualityResult))]
[JsonSerializable(typeof(SuspiciousRegion))]
[JsonSerializable(typeof(SuspiciousRegionsResult))]
[JsonSerializable(typeof(LengthDistributionResult))]
public partial class AppJsonContext : JsonSerializerContext
{
}
