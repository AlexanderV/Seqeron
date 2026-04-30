using System.Collections.Generic;
using System.Text.Json.Serialization;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population;

[JsonSourceGenerationOptions(WriteIndented = false)]
[JsonSerializable(typeof(string))]
// --- Allele Frequency
[JsonSerializable(typeof(AlleleFrequenciesResult))]
[JsonSerializable(typeof(MinorAlleleFrequencyResult))]
[JsonSerializable(typeof(int[]))]
[JsonSerializable(typeof(VariantItem))]
[JsonSerializable(typeof(VariantItem[]))]
[JsonSerializable(typeof(FilterVariantsByMafResult))]
// --- Diversity Statistics
[JsonSerializable(typeof(string[]))]
[JsonSerializable(typeof(NucleotideDiversityResult))]
[JsonSerializable(typeof(WattersonsThetaResult))]
[JsonSerializable(typeof(TajimasDResult))]
[JsonSerializable(typeof(DiversityStatisticsResult))]
// --- Hardy-Weinberg
[JsonSerializable(typeof(HardyWeinbergTestResult))]
// --- F-statistics
[JsonSerializable(typeof(AlleleItem))]
[JsonSerializable(typeof(AlleleItem[]))]
[JsonSerializable(typeof(PopulationItem))]
[JsonSerializable(typeof(PopulationItem[]))]
[JsonSerializable(typeof(FstResult))]
[JsonSerializable(typeof(PairwiseFstResult))]
[JsonSerializable(typeof(VariantDataItem))]
[JsonSerializable(typeof(VariantDataItem[]))]
[JsonSerializable(typeof(FStatisticsResult))]
// --- Linkage Disequilibrium
[JsonSerializable(typeof(GenotypePairItem))]
[JsonSerializable(typeof(GenotypePairItem[]))]
[JsonSerializable(typeof(LinkageDisequilibriumResult))]
[JsonSerializable(typeof(VariantGenotypesItem))]
[JsonSerializable(typeof(VariantGenotypesItem[]))]
[JsonSerializable(typeof(HaplotypeFrequencyItem))]
[JsonSerializable(typeof(HaplotypeBlockItem))]
[JsonSerializable(typeof(HaplotypeBlocksResult))]
// --- Selection
[JsonSerializable(typeof(double[]))]
[JsonSerializable(typeof(IhsResult))]
[JsonSerializable(typeof(RegionItem))]
[JsonSerializable(typeof(RegionItem[]))]
[JsonSerializable(typeof(SelectionSignalItem))]
[JsonSerializable(typeof(ScanSelectionSignalsResult))]
[JsonSerializable(typeof(IndividualItem))]
[JsonSerializable(typeof(IndividualItem[]))]
[JsonSerializable(typeof(RefPopItem))]
[JsonSerializable(typeof(RefPopItem[]))]
[JsonSerializable(typeof(AncestryProportionItem))]
[JsonSerializable(typeof(EstimateAncestryResult))]
[JsonSerializable(typeof(Dictionary<string, double>))]
[JsonSerializable(typeof(IReadOnlyDictionary<string, double>))]
// --- Inbreeding / ROH
[JsonSerializable(typeof(RohSegmentItem))]
[JsonSerializable(typeof(RohSegmentItem[]))]
[JsonSerializable(typeof(InbreedingFromRohResult))]
[JsonSerializable(typeof(GenotypePositionItem))]
[JsonSerializable(typeof(GenotypePositionItem[]))]
[JsonSerializable(typeof(RohRegionItem))]
[JsonSerializable(typeof(RunsOfHomozygosityResult))]
public partial class AppJsonContext : JsonSerializerContext
{
}
