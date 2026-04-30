using System.Collections.Generic;

namespace Seqeron.Mcp.Metagenomics.Tools;

// =====================================================================
// Shared input DTOs
// =====================================================================

/// <summary>Input read for taxonomic classification.</summary>
public record ReadInput(string Id, string Sequence);

/// <summary>K-mer → taxon mapping entry (flattened from a dictionary).</summary>
public record KmerDatabaseEntry(string Kmer, string TaxonId);

/// <summary>Reference genome (taxon id + nucleotide sequence) for k-mer DB construction.</summary>
public record ReferenceGenomeInput(string TaxonId, string Sequence);

/// <summary>Single (species/taxon → fractional abundance) entry.</summary>
public record AbundanceItem(string Name, double Fraction);

/// <summary>One sample's abundance vector for differential-abundance analysis.</summary>
public record AbundanceSample(IReadOnlyList<AbundanceItem> Items);

/// <summary>Contig with coverage for genome binning.</summary>
public record ContigInput(string ContigId, string Sequence, double Coverage);

/// <summary>Predicted protein for functional annotation.</summary>
public record ProteinInput(string GeneId, string ProteinSequence);

/// <summary>Function-database entry (motif → function/pathway/KO).</summary>
public record FunctionDatabaseEntry(string Motif, string Function, string Pathway, string Ko);

/// <summary>Gene (id + nucleotide sequence) used as input by several tools.</summary>
public record GeneInput(string GeneId, string Sequence);

/// <summary>Resistance-gene database entry (motif → name/antibiotic class).</summary>
public record ResistanceDatabaseEntry(string Motif, string Name, string AntibioticClass);

/// <summary>Genome (id + ordered list of genes) for pan-genome tools.</summary>
public record GenomeInput(string GenomeId, IReadOnlyList<GeneInput> Genes);

// =====================================================================
// MetagenomicsAnalyzer result DTOs
// =====================================================================

/// <summary>Result of classify_reads.</summary>
public record ClassifyReadsResult(
    IReadOnlyList<global::Seqeron.Genomics.Metagenomics.MetagenomicsAnalyzer.TaxonomicClassification> Items);

/// <summary>Result of build_kmer_database.</summary>
public record BuildKmerDatabaseResult(IReadOnlyList<KmerDatabaseEntry> Entries, int Count);

/// <summary>Result of taxonomic_profile (dictionaries flattened to lists).</summary>
public record TaxonomicProfileResult(
    IReadOnlyList<AbundanceItem> KingdomAbundance,
    IReadOnlyList<AbundanceItem> PhylumAbundance,
    IReadOnlyList<AbundanceItem> GenusAbundance,
    IReadOnlyList<AbundanceItem> SpeciesAbundance,
    double ShannonDiversity,
    double SimpsonDiversity,
    int TotalReads,
    int ClassifiedReads);

/// <summary>Result of bin_contigs.</summary>
public record BinContigsResult(
    IReadOnlyList<global::Seqeron.Genomics.Metagenomics.MetagenomicsAnalyzer.GenomeBin> Items);

/// <summary>Result of predict_functions.</summary>
public record PredictFunctionsResult(
    IReadOnlyList<global::Seqeron.Genomics.Metagenomics.MetagenomicsAnalyzer.FunctionalAnnotation> Items);

/// <summary>Pathway-occurrence count item used by functional_diversity.</summary>
public record PathwayCountItem(string Pathway, int Count);

/// <summary>Result of functional_diversity.</summary>
public record FunctionalDiversityResult(
    double FunctionalRichness,
    double FunctionalDiversity,
    IReadOnlyList<PathwayCountItem> PathwayCounts);

/// <summary>One predicted resistance-gene hit.</summary>
public record ResistanceGeneItem(
    string GeneId,
    string ResistanceGene,
    string AntibioticClass,
    double Identity);

/// <summary>Result of find_resistance_genes.</summary>
public record FindResistanceGenesResult(IReadOnlyList<ResistanceGeneItem> Items);

/// <summary>One taxon's differential-abundance result (Welch t-test).</summary>
public record DifferentialAbundanceItem(
    string Taxon,
    double FoldChange,
    double PValue,
    bool Significant);

/// <summary>Result of differential_abundance.</summary>
public record DifferentialAbundanceResult(IReadOnlyList<DifferentialAbundanceItem> Items);

// =====================================================================
// PanGenomeAnalyzer result DTOs
// =====================================================================

/// <summary>(genome → its gene ids) mapping entry.</summary>
public record GenomeToGenesItem(string GenomeId, IReadOnlyList<string> GeneIds);

/// <summary>Pan-genome statistics (enum serialized as string).</summary>
public record PanGenomeStatisticsDto(
    int TotalGenomes,
    int TotalGenes,
    int CoreGeneCount,
    int AccessoryGeneCount,
    int UniqueGeneCount,
    double CoreFraction,
    double GenomeFluidity,
    string Type);

/// <summary>Result of construct_pangenome (with dictionaries flattened).</summary>
public record PanGenomeResultDto(
    IReadOnlyList<string> CoreGenes,
    IReadOnlyList<string> AccessoryGenes,
    IReadOnlyList<string> UniqueGenes,
    IReadOnlyList<GenomeToGenesItem> GenomeToGenes,
    PanGenomeStatisticsDto Statistics);

/// <summary>Result of cluster_genes.</summary>
public record ClusterGenesResult(
    IReadOnlyList<global::Seqeron.Genomics.Metagenomics.PanGenomeAnalyzer.GeneCluster> Items);

/// <summary>One (cluster → present?) entry within a presence/absence row.</summary>
public record GenePresenceItem(string ClusterId, bool Present);

/// <summary>One row of the presence/absence matrix (per genome).</summary>
public record GenePresenceRowDto(
    string GenomeId,
    IReadOnlyList<GenePresenceItem> GenePresence,
    int TotalGenes,
    int PresentGenes);

/// <summary>Result of gene_presence_absence_matrix.</summary>
public record PresenceAbsenceMatrixResult(IReadOnlyList<GenePresenceRowDto> Items);

/// <summary>Result of fit_heaps_law (Func predictor field dropped).</summary>
public record HeapsLawFitDto(double K, double Gamma, double RSquared);

/// <summary>Result of core_gene_clusters.</summary>
public record CoreGeneClustersResult(
    IReadOnlyList<global::Seqeron.Genomics.Metagenomics.PanGenomeAnalyzer.GeneCluster> Items);

/// <summary>Result of core_genome_alignment (concatenated alignment).</summary>
public record CoreGenomeAlignmentResult(string Result);

/// <summary>One accessory-cluster summary.</summary>
public record AccessoryGeneItem(
    string ClusterId,
    IReadOnlyList<string> GenomesWithGene,
    double Frequency);

/// <summary>Result of accessory_genes.</summary>
public record AccessoryGenesResult(IReadOnlyList<AccessoryGeneItem> Items);

/// <summary>Genes unique to a single genome.</summary>
public record GenomeSpecificGeneItem(string GenomeId, IReadOnlyList<string> UniqueGeneIds);

/// <summary>Result of find_genome_specific_genes.</summary>
public record GenomeSpecificGenesResult(IReadOnlyList<GenomeSpecificGeneItem> Items);

/// <summary>Result of select_phylogenetic_markers.</summary>
public record SelectPhylogeneticMarkersResult(
    IReadOnlyList<global::Seqeron.Genomics.Metagenomics.PanGenomeAnalyzer.GeneCluster> Items);
