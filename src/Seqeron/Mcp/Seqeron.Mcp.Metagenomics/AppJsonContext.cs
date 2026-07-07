using System.Text.Json.Serialization;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics;

[JsonSourceGenerationOptions(WriteIndented = false)]
[JsonSerializable(typeof(string))]
[JsonSerializable(typeof(int))]
[JsonSerializable(typeof(double))]
[JsonSerializable(typeof(bool))]

// ----- Inputs (root) -----
[JsonSerializable(typeof(IReadOnlyList<ReadInput>))]
[JsonSerializable(typeof(IReadOnlyList<KmerDatabaseEntry>))]
[JsonSerializable(typeof(IReadOnlyList<ReferenceGenomeInput>))]
[JsonSerializable(typeof(IReadOnlyList<TaxonNodeInput>))]
[JsonSerializable(typeof(IReadOnlyList<AbundanceItem>))]
[JsonSerializable(typeof(IReadOnlyList<AbundanceSample>))]
[JsonSerializable(typeof(IReadOnlyList<ContigInput>))]
[JsonSerializable(typeof(IReadOnlyList<ProteinInput>))]
[JsonSerializable(typeof(IReadOnlyList<FunctionDatabaseEntry>))]
[JsonSerializable(typeof(IReadOnlyList<GeneInput>))]
[JsonSerializable(typeof(IReadOnlyList<ResistanceDatabaseEntry>))]
[JsonSerializable(typeof(IReadOnlyList<GenomeInput>))]
[JsonSerializable(typeof(IReadOnlyList<MetagenomicsAnalyzer.TaxonomicClassification>))]
[JsonSerializable(typeof(IReadOnlyList<MetagenomicsAnalyzer.FunctionalAnnotation>))]
[JsonSerializable(typeof(IReadOnlyList<PanGenomeAnalyzer.GeneCluster>))]

// ----- Element types (registered explicitly for AOT closure) -----
[JsonSerializable(typeof(ReadInput))]
[JsonSerializable(typeof(KmerDatabaseEntry))]
[JsonSerializable(typeof(ReferenceGenomeInput))]
[JsonSerializable(typeof(TaxonNodeInput))]
[JsonSerializable(typeof(AbundanceItem))]
[JsonSerializable(typeof(AbundanceSample))]
[JsonSerializable(typeof(ContigInput))]
[JsonSerializable(typeof(ProteinInput))]
[JsonSerializable(typeof(FunctionDatabaseEntry))]
[JsonSerializable(typeof(GeneInput))]
[JsonSerializable(typeof(ResistanceDatabaseEntry))]
[JsonSerializable(typeof(GenomeInput))]
[JsonSerializable(typeof(MetagenomicsAnalyzer.TaxonomicClassification))]
[JsonSerializable(typeof(MetagenomicsAnalyzer.FunctionalAnnotation))]
[JsonSerializable(typeof(MetagenomicsAnalyzer.GenomeBin))]
[JsonSerializable(typeof(MetagenomicsAnalyzer.AlphaDiversity))]
[JsonSerializable(typeof(MetagenomicsAnalyzer.BetaDiversity))]
[JsonSerializable(typeof(PanGenomeAnalyzer.GeneCluster))]

// ----- Output result DTOs (root) -----
[JsonSerializable(typeof(ClassifyReadsResult))]
[JsonSerializable(typeof(BuildKmerDatabaseResult))]
[JsonSerializable(typeof(TaxonomicProfileResult))]
[JsonSerializable(typeof(BinContigsResult))]
[JsonSerializable(typeof(PredictFunctionsResult))]
[JsonSerializable(typeof(FunctionalDiversityResult))]
[JsonSerializable(typeof(FindResistanceGenesResult))]
[JsonSerializable(typeof(DifferentialAbundanceResult))]
[JsonSerializable(typeof(PanGenomeResultDto))]
[JsonSerializable(typeof(ClusterGenesResult))]
[JsonSerializable(typeof(PresenceAbsenceMatrixResult))]
[JsonSerializable(typeof(HeapsLawFitDto))]
[JsonSerializable(typeof(CoreGeneClustersResult))]
[JsonSerializable(typeof(CoreGenomeAlignmentResult))]
[JsonSerializable(typeof(AccessoryGenesResult))]
[JsonSerializable(typeof(GenomeSpecificGenesResult))]
[JsonSerializable(typeof(SelectPhylogeneticMarkersResult))]

// ----- Nested item DTOs referenced from result DTOs -----
[JsonSerializable(typeof(PathwayCountItem))]
[JsonSerializable(typeof(ResistanceGeneItem))]
[JsonSerializable(typeof(DifferentialAbundanceItem))]
[JsonSerializable(typeof(GenomeToGenesItem))]
[JsonSerializable(typeof(PanGenomeStatisticsDto))]
[JsonSerializable(typeof(GenePresenceItem))]
[JsonSerializable(typeof(GenePresenceRowDto))]
[JsonSerializable(typeof(AccessoryGeneItem))]
[JsonSerializable(typeof(GenomeSpecificGeneItem))]
internal partial class AppJsonContext : JsonSerializerContext
{
}
