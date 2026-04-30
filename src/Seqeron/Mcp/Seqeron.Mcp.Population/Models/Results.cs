namespace Seqeron.Mcp.Population.Tools;

// ================================
// Allele Frequency
// ================================

/// <summary>Result of allele_frequencies operation.</summary>
public record AlleleFrequenciesResult(double MajorFreq, double MinorFreq);

/// <summary>Result of minor_allele_frequency operation.</summary>
public record MinorAlleleFrequencyResult(double Maf);

/// <summary>
/// Genetic variant with allele frequency data. Mirrors
/// <c>PopulationGeneticsAnalyzer.Variant</c> as an MCP-friendly DTO.
/// </summary>
public record VariantItem(
    string Id,
    string Chromosome,
    int Position,
    string ReferenceAllele,
    string AlternateAllele,
    double AlleleFrequency,
    int SampleCount);

/// <summary>Result of filter_variants_by_maf operation.</summary>
public record FilterVariantsByMafResult(IReadOnlyList<VariantItem> Items);

// ================================
// Diversity Statistics
// ================================

/// <summary>Result of nucleotide_diversity operation.</summary>
public record NucleotideDiversityResult(double Pi);

/// <summary>Result of wattersons_theta operation.</summary>
public record WattersonsThetaResult(double Theta);

/// <summary>Result of tajimas_d operation.</summary>
public record TajimasDResult(double TajimasD);

/// <summary>Result of diversity_statistics operation.</summary>
public record DiversityStatisticsResult(
    double NucleotideDiversity,
    double WattersonTheta,
    double TajimasD,
    int SegregatingSites,
    int SampleSize,
    double HeterozygosityObserved,
    double HeterozygosityExpected);

// ================================
// Hardy-Weinberg
// ================================

/// <summary>Result of hardy_weinberg_test operation.</summary>
public record HardyWeinbergTestResult(
    string VariantId,
    int ObservedAA,
    int ObservedAa,
    int Observedaa,
    double ExpectedAA,
    double ExpectedAa,
    double Expectedaa,
    double ChiSquare,
    double PValue,
    bool InEquilibrium);

// ================================
// F-statistics
// ================================

/// <summary>Per-variant allele frequency entry for Fst computation.</summary>
public record AlleleItem(double AlleleFreq, int SampleSize);

/// <summary>Per-population allele frequency vector for pairwise Fst.</summary>
public record PopulationItem(string PopulationId, IReadOnlyList<AlleleItem> Variants);

/// <summary>Result of fst operation.</summary>
public record FstResult(double Fst);

/// <summary>Result of pairwise_fst operation.</summary>
public record PairwiseFstResult(
    IReadOnlyList<string> PopulationIds,
    IReadOnlyList<IReadOnlyList<double>> Matrix);

/// <summary>Per-variant input row for f_statistics operation.</summary>
public record VariantDataItem(
    int HetObs1,
    int N1,
    int HetObs2,
    int N2,
    double AlleleFreq1,
    double AlleleFreq2);

/// <summary>Result of f_statistics operation.</summary>
public record FStatisticsResult(
    double Fst,
    double Fis,
    double Fit,
    string Population1,
    string Population2);

// ================================
// Linkage Disequilibrium
// ================================

/// <summary>Pair of genotypes (0/1/2 encoding) at two variants.</summary>
public record GenotypePairItem(int Geno1, int Geno2);

/// <summary>Result of linkage_disequilibrium operation.</summary>
public record LinkageDisequilibriumResult(
    string Variant1,
    string Variant2,
    double DPrime,
    double RSquared,
    double Distance);

/// <summary>Variant with positional and per-individual genotype data.</summary>
public record VariantGenotypesItem(
    string VariantId,
    int Position,
    IReadOnlyList<int> Genotypes);

/// <summary>Haplotype frequency entry inside a haplotype block.</summary>
public record HaplotypeFrequencyItem(string Haplotype, double Frequency);

/// <summary>A haplotype block discovered by haplotype_blocks.</summary>
public record HaplotypeBlockItem(
    int Start,
    int End,
    IReadOnlyList<string> Variants,
    IReadOnlyList<HaplotypeFrequencyItem> Haplotypes);

/// <summary>Result of haplotype_blocks operation.</summary>
public record HaplotypeBlocksResult(IReadOnlyList<HaplotypeBlockItem> Items);

// ================================
// Selection
// ================================

/// <summary>Result of integrated_haplotype_score operation.</summary>
public record IhsResult(double Ihs);

/// <summary>Per-region input for scan_selection_signals.</summary>
public record RegionItem(
    string Region,
    int Start,
    int End,
    double TajimaD,
    double Fst,
    double IHS);

/// <summary>A selection signal emitted by scan_selection_signals.</summary>
public record SelectionSignalItem(
    string Region,
    int Start,
    int End,
    double Score,
    string TestType,
    double PValue,
    string Interpretation);

/// <summary>Result of scan_selection_signals operation.</summary>
public record ScanSelectionSignalsResult(IReadOnlyList<SelectionSignalItem> Items);

/// <summary>Individual with per-SNP genotypes (0/1/2).</summary>
public record IndividualItem(string IndividualId, IReadOnlyList<int> Genotypes);

/// <summary>Reference population with per-SNP allele frequencies.</summary>
public record RefPopItem(string PopulationId, IReadOnlyList<double> AlleleFrequencies);

/// <summary>Ancestry proportion for one individual across reference populations.</summary>
public record AncestryProportionItem(
    string IndividualId,
    IReadOnlyDictionary<string, double> Proportions);

/// <summary>Result of estimate_ancestry operation.</summary>
public record EstimateAncestryResult(IReadOnlyList<AncestryProportionItem> Items);

// ================================
// Inbreeding / ROH
// ================================

/// <summary>An ROH segment with start/end coordinates.</summary>
public record RohSegmentItem(int Start, int End);

/// <summary>Result of inbreeding_from_roh operation.</summary>
public record InbreedingFromRohResult(double InbreedingCoefficient);

/// <summary>Per-position genotype call for runs_of_homozygosity input.</summary>
public record GenotypePositionItem(int Position, int Genotype);

/// <summary>An ROH region produced by runs_of_homozygosity.</summary>
public record RohRegionItem(int Start, int End, int SnpCount);

/// <summary>Result of runs_of_homozygosity operation.</summary>
public record RunsOfHomozygosityResult(IReadOnlyList<RohRegionItem> Items);
