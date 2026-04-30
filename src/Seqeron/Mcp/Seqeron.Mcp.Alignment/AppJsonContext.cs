using System.Text.Json.Serialization;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment;

[JsonSourceGenerationOptions(WriteIndented = false)]
// Primitives & arrays used at the tool boundary
[JsonSerializable(typeof(string))]
[JsonSerializable(typeof(int))]
[JsonSerializable(typeof(double))]
[JsonSerializable(typeof(bool))]
[JsonSerializable(typeof(string[]))]
[JsonSerializable(typeof(int[]))]
// Inputs
[JsonSerializable(typeof(ContigLinkInput))]
[JsonSerializable(typeof(ContigLinkInput[]))]
[JsonSerializable(typeof(QualityReadInput))]
[JsonSerializable(typeof(QualityReadInput[]))]
// SequenceAligner outputs
[JsonSerializable(typeof(AlignmentResultDto))]
[JsonSerializable(typeof(AlignmentStatisticsDto))]
[JsonSerializable(typeof(FormatAlignmentResult))]
[JsonSerializable(typeof(MultipleAlignmentResultDto))]
// ApproximateMatcher outputs
[JsonSerializable(typeof(ApproximateMatchDto))]
[JsonSerializable(typeof(ApproximateMatchDto[]))]
[JsonSerializable(typeof(ApproximateMatchListResult))]
[JsonSerializable(typeof(FindBestMatchResult))]
[JsonSerializable(typeof(FrequentKmerItem))]
[JsonSerializable(typeof(FrequentKmerItem[]))]
[JsonSerializable(typeof(FrequentKmersResult))]
// SequenceAssembler outputs
[JsonSerializable(typeof(AssemblyResultDto))]
[JsonSerializable(typeof(OverlapDto))]
[JsonSerializable(typeof(OverlapDto[]))]
[JsonSerializable(typeof(OverlapListResult))]
[JsonSerializable(typeof(FindOverlapInfo))]
[JsonSerializable(typeof(FindOverlapResult))]
[JsonSerializable(typeof(IdentityResult))]
[JsonSerializable(typeof(MergedContigResult))]
[JsonSerializable(typeof(ScaffoldsResult))]
[JsonSerializable(typeof(CoverageResult))]
[JsonSerializable(typeof(ConsensusResult))]
[JsonSerializable(typeof(TrimmedReadsResult))]
[JsonSerializable(typeof(CorrectedReadsResult))]
public partial class AppJsonContext : JsonSerializerContext
{
}
