namespace Seqeron.Mcp.Annotation.Tools;

// ================================
// Shared DTOs
// ================================

/// <summary>Open reading frame record.</summary>
public record OrfDto(
    int Start,
    int End,
    int Frame,
    bool IsReverseComplement,
    string Sequence,
    string ProteinSequence);

/// <summary>Longest ORF per reading frame entry.</summary>
public record LongestOrfFrame(int Frame, OrfDto? Orf);

/// <summary>Predicted gene annotation (0-based half-open Start/End).</summary>
/// <param name="Source">GFF3 column 2; null or empty → ".".</param>
/// <param name="Score">GFF3 column 6; null → ".".</param>
/// <param name="AttributeValues">Multi-valued GFF3 attributes (e.g. Parent → [m1, m2]), written comma-separated.</param>
public record GeneAnnotationDto(
    string GeneId,
    int Start,
    int End,
    char Strand,
    string Type,
    string Product,
    Dictionary<string, string> Attributes,
    string? Source = null,
    double? Score = null,
    Dictionary<string, List<string>>? AttributeValues = null);

/// <summary>Generic GFF3 genomic feature (1-based inclusive Start/End).</summary>
/// <param name="Seqid">GFF3 column 1, percent-decoded.</param>
/// <param name="Source">GFF3 column 2, percent-decoded.</param>
/// <param name="AttributeValues">Column 9 as tag → value list (split on the unescaped ',', each value decoded).</param>
public record GenomicFeatureDto(
    string FeatureId,
    string Type,
    int Start,
    int End,
    char Strand,
    double? Score,
    int? Phase,
    Dictionary<string, string> Attributes,
    string? Seqid = null,
    string? Source = null,
    Dictionary<string, List<string>>? AttributeValues = null);

/// <summary>Ribosome binding site (Shine–Dalgarno) hit.</summary>
public record RibosomeBindingSite(int Position, string Sequence, double Score);

/// <summary>Promoter motif hit (-10/-35 box).</summary>
public record PromoterMotif(int Position, string Type, string Sequence, double Score);

/// <summary>Repetitive element (tandem or inverted repeat).</summary>
public record RepetitiveElement(int Start, int End, string Type, string Sequence);

// ================================
// Tool result wrappers
// ================================

/// <summary>Result of find_orfs.</summary>
public record FindOrfsResult(IReadOnlyList<OrfDto> Orfs);

/// <summary>Result of longest_orfs_per_frame.</summary>
public record LongestOrfsPerFrameResult(IReadOnlyList<LongestOrfFrame> Frames);

/// <summary>Result of find_ribosome_binding_sites.</summary>
public record FindRibosomeBindingSitesResult(IReadOnlyList<RibosomeBindingSite> Sites);

/// <summary>Result of predict_genes.</summary>
public record PredictGenesResult(IReadOnlyList<GeneAnnotationDto> Genes);

/// <summary>Result of parse_gff3.</summary>
public record ParseGff3Result(IReadOnlyList<GenomicFeatureDto> Features);

/// <summary>Result of to_gff3.</summary>
public record ToGff3Result(IReadOnlyList<string> Lines);

/// <summary>Result of find_promoter_motifs.</summary>
public record FindPromoterMotifsResult(IReadOnlyList<PromoterMotif> Motifs);

/// <summary>Result of coding_potential.</summary>
public record CodingPotentialResult(double Score);

/// <summary>Result of fickett_score.</summary>
public record FickettScoreResult(double Score, string Call);

/// <summary>One fickett_scan window (0-based half-open).</summary>
public record FickettWindowDto(int Start, int End, double Score, string Call);

/// <summary>Result of fickett_scan.</summary>
public record FickettScanResult(IReadOnlyList<FickettWindowDto> Windows);

/// <summary>One scored CPAT ORF candidate (CPAT 1-based coordinates; '-' ORFs have Start &gt; End).</summary>
public record CpatOrfDto(string Strand, int Frame, int Start, int End, int OrfLength, string OrfSequence,
    int MrnaLength, double FickettScore, double HexamerScore, double CodingProbability, bool IsCoding);

/// <summary>Result of cpat_coding_probability: candidates in CPAT order and the index of the best one (-1 if none).</summary>
public record CpatCodingProbabilityResult(IReadOnlyList<CpatOrfDto> Candidates, int BestIndex);

/// <summary>Result of find_repetitive_elements.</summary>
public record FindRepetitiveElementsResult(IReadOnlyList<RepetitiveElement> Repeats);

/// <summary>Result of codon_usage.</summary>
public record CodonUsageResult(Dictionary<string, int> Usage);

/// <summary>Result of pooled_rscu: sense codon (upper-case DNA) to RSCU.</summary>
public record PooledRscuResult(Dictionary<string, double> Rscu);

// ================================
// VariantCaller — shared DTOs
// ================================

/// <summary>A detected genetic variant.</summary>
public record VariantDto(
    int Position,
    string ReferenceAllele,
    string AlternateAllele,
    string Type,
    int QueryPosition);

/// <summary>A variant with functional annotation.</summary>
public record AnnotatedVariantDto(
    VariantDto Variant,
    string Effect,
    string MutationType);

/// <summary>Summary variant statistics between two sequences.</summary>
public record VariantStatisticsDto(
    int TotalVariants,
    int Snps,
    int Insertions,
    int Deletions,
    double TiTvRatio,
    double VariantDensity,
    int ReferenceLength,
    int QueryLength);

// ================================
// VariantCaller — tool result wrappers
// ================================

/// <summary>Result of call_variants / call_variants_from_alignment / find_* / find_indels.</summary>
public record CallVariantsResult(IReadOnlyList<VariantDto> Variants);

/// <summary>Result of classify_mutation.</summary>
public record ClassifyMutationResult(string MutationType);

/// <summary>Result of titv_ratio.</summary>
public record TiTvRatioResult(double Ratio);

/// <summary>Result of predict_variant_effect.</summary>
public record PredictVariantEffectResult(string Effect);

/// <summary>Result of annotate_variants.</summary>
public record AnnotateVariantsResult(IReadOnlyList<AnnotatedVariantDto> Annotated);

/// <summary>Result of variants_to_vcf.</summary>
public record VariantsToVcfResult(IReadOnlyList<string> Lines);

// ================================
// VariantAnnotator — shared DTOs
// ================================

/// <summary>A genetic variant in genomic coordinates (annotator flavour).</summary>
public record AnnotatorVariantDto(
    string Chromosome,
    int Position,
    string Reference,
    string Alternate,
    string Type,
    double? Quality,
    string? Id);

/// <summary>Half-open genomic interval (inclusive Start/End in source).</summary>
public record GenomicIntervalDto(int Start, int End);

/// <summary>Transcript model used for variant consequence prediction.</summary>
public record TranscriptDto(
    string TranscriptId,
    string GeneId,
    string GeneName,
    string Chromosome,
    int Start,
    int End,
    char Strand,
    IReadOnlyList<GenomicIntervalDto> Exons,
    IReadOnlyList<GenomicIntervalDto> CodingExons,
    int? CdsStart,
    int? CdsEnd);

/// <summary>VEP-like variant annotation against a single transcript.</summary>
public record VariantAnnotationDto(
    AnnotatorVariantDto Variant,
    string TranscriptId,
    string GeneId,
    string GeneName,
    string Consequence,
    string Impact,
    string? CodonChange,
    string? AminoAcidChange,
    int? ProteinPosition,
    int? CdsPosition,
    double? SiftScore,
    double? PolyphenScore,
    double? CaddScore,
    string? ExistingVariation,
    Dictionary<string, double>? PopulationFrequencies);

/// <summary>ACMG-like pathogenicity prediction.</summary>
public record PathogenicityPredictionDto(
    string Classification,
    double ConfidenceScore,
    IReadOnlyList<string> EvidenceCriteria,
    double? ClinicalSignificance,
    bool IsActionable);

/// <summary>Per-position multi-species conservation score input.</summary>
public record ConservationPositionInputDto(
    string Chromosome,
    int Position,
    string SpeciesAlleles);

/// <summary>Per-position conservation scores (PhyloP, PhastCons, GERP).</summary>
public record ConservationScoreDto(
    string Chromosome,
    int Position,
    double PhyloP,
    double PhastCons,
    double Gerp,
    int ConservedSpeciesCount);

/// <summary>A conserved genomic element (run of high-scoring positions).</summary>
public record ConservedElementDto(
    string Chromosome,
    int Start,
    int End,
    double Score);

/// <summary>Regulatory region input for overlap annotation.</summary>
public record RegulatoryRegionInputDto(
    string Chromosome,
    int Start,
    int End,
    string Type,
    string? CellType,
    double? Score,
    IReadOnlyList<string> TranscriptionFactors);

/// <summary>Regulatory region annotated as overlapping a variant.</summary>
public record RegulatoryAnnotationDto(
    string Chromosome,
    int Start,
    int End,
    string FeatureType,
    string? CellType,
    double? Score,
    IReadOnlyList<string> TranscriptionFactors);

/// <summary>Transcription-factor binding motif input.</summary>
public record TfMotifInputDto(string TfName, string Motif, double Threshold);

/// <summary>Transcription-factor binding score change induced by a SNV.</summary>
public record TfBindingChangeDto(
    string TfName,
    double RefScore,
    double AltScore,
    double ScoreDifference);

// ================================
// VariantAnnotator — tool result wrappers
// ================================

/// <summary>Result of classify_variant.</summary>
public record ClassifyVariantResult(string VariantType);

/// <summary>Result of normalize_variant / parse_vcf_variant.</summary>
public record AnnotatorVariantResult(AnnotatorVariantDto Variant);

/// <summary>Result of annotate_variant_on_transcripts.</summary>
public record AnnotateVariantOnTranscriptsResult(IReadOnlyList<VariantAnnotationDto> Annotations);

/// <summary>Result of impact_level.</summary>
public record ImpactLevelResult(string Impact);

/// <summary>Result of predict_pathogenicity.</summary>
public record PredictPathogenicityResult(PathogenicityPredictionDto Prediction);

/// <summary>Result of calculate_conservation.</summary>
public record CalculateConservationResult(IReadOnlyList<ConservationScoreDto> Scores);

/// <summary>Result of find_conserved_elements.</summary>
public record FindConservedElementsResult(IReadOnlyList<ConservedElementDto> Elements);

/// <summary>Result of annotate_regulatory_elements.</summary>
public record AnnotateRegulatoryElementsResult(IReadOnlyList<RegulatoryAnnotationDto> Annotations);

/// <summary>Result of predict_tf_binding_change.</summary>
public record PredictTfBindingChangeResult(IReadOnlyList<TfBindingChangeDto> Changes);

/// <summary>Result of format_vcf_info.</summary>
public record FormatVcfInfoResult(string Info);

// ================================
// EpigeneticsAnalyzer — shared DTOs
// ================================

/// <summary>Methylation site (used for both input and output).</summary>
public record MethylationSiteDto(
    int Position,
    string Type,
    string Context,
    double MethylationLevel,
    int Coverage);

/// <summary>Position-level methylation entry within a profile.</summary>
public record MethylationByPositionDto(int Position, double Level);

/// <summary>Aggregated methylation profile.</summary>
public record MethylationProfileDto(
    double GlobalMethylation,
    double CpGMethylation,
    double CHGMethylation,
    double CHHMethylation,
    int TotalCpGSites,
    int MethylatedCpGSites,
    IReadOnlyList<MethylationByPositionDto> MethylationByPosition);

/// <summary>CpG island record.</summary>
public record CpGIslandDto(int Start, int End, double GcContent, double CpGRatio);

/// <summary>Bisulfite read input.</summary>
public record BisulfiteReadDto(string ReadSequence, int StartPosition);

/// <summary>Differentially methylated region.</summary>
public record DmrDto(
    int Start,
    int End,
    double MeanDifference,
    double PValue,
    int CpGCount,
    string Annotation);

/// <summary>Histone modification interval (input to annotation).</summary>
public record HistoneModificationInputDto(int Start, int End, string Mark, double Signal);

/// <summary>Histone modification with predicted chromatin state.</summary>
public record HistoneModificationAnnotationDto(
    int Start,
    int End,
    string Mark,
    double Signal,
    string PredictedState);

/// <summary>Per-position chromatin accessibility signal (input).</summary>
public record AccessibilitySignalDto(int Position, double Signal);

/// <summary>Accessible chromatin region (output).</summary>
public record AccessibilityRegionDto(
    int Start,
    int End,
    double AccessibilityScore,
    string PeakType,
    IReadOnlyList<string> NearbyGenes);

/// <summary>Allele-specific methylation gene input for imprinting prediction.</summary>
public record ImprintedGeneInputDto(
    string GeneId,
    int Start,
    int End,
    double MaternalMethylation,
    double PaternalMethylation);

/// <summary>Predicted imprinted gene output.</summary>
public record ImprintedGeneDto(
    string GeneId,
    int Start,
    int End,
    double ImprintingScore,
    string ParentalOrigin,
    bool HasDMR);

// ================================
// EpigeneticsAnalyzer — tool result wrappers
// ================================

/// <summary>Result of find_cpg_sites.</summary>
public record FindCpGSitesResult(IReadOnlyList<int> Positions);

/// <summary>Result of find_methylation_sites / methylation_from_bisulfite.</summary>
public record FindMethylationSitesResult(IReadOnlyList<MethylationSiteDto> Sites);

/// <summary>Result of cpg_observed_expected.</summary>
public record CpGObservedExpectedResult(double Ratio);

/// <summary>Result of find_cpg_islands.</summary>
public record FindCpGIslandsResult(IReadOnlyList<CpGIslandDto> Islands);

/// <summary>Result of simulate_bisulfite_conversion.</summary>
public record SimulateBisulfiteConversionResult(string Converted);

/// <summary>Result of find_dmrs.</summary>
public record FindDmrsResult(IReadOnlyList<DmrDto> Regions);

/// <summary>Result of predict_chromatin_state.</summary>
public record PredictChromatinStateResult(string State);

/// <summary>Result of annotate_histone_modifications.</summary>
public record AnnotateHistoneModificationsResult(IReadOnlyList<HistoneModificationAnnotationDto> Annotations);

/// <summary>Result of find_accessible_regions.</summary>
public record FindAccessibleRegionsResult(IReadOnlyList<AccessibilityRegionDto> Regions);

/// <summary>Result of predict_imprinted_genes.</summary>
public record PredictImprintedGenesResult(IReadOnlyList<ImprintedGeneDto> Imprinted);

/// <summary>Result of epigenetic_age.</summary>
public record EpigeneticAgeResult(double Age);

// ================================
// MiRnaAnalyzer — shared DTOs
// ================================

/// <summary>microRNA record (T→U normalized) with seed metadata.</summary>
public record MiRnaDto(
    string Name,
    string Sequence,
    string SeedSequence,
    int SeedStart,
    int SeedEnd);

/// <summary>Predicted miRNA target site.</summary>
public record TargetSiteDto(
    int Start,
    int End,
    string TargetSequence,
    string MiRnaName,
    string Type,
    int SeedMatchLength,
    double Score,
    double FreeEnergy,
    string Alignment);

/// <summary>Pre-miRNA hairpin candidate.</summary>
public record PreMiRnaHairpinDto(
    int Start,
    int End,
    string Sequence,
    string MatureSequence,
    string StarSequence,
    string Structure,
    double FreeEnergy);

/// <summary>Group of miRNAs sharing an identical seed sequence.</summary>
public record SeedFamilyDto(string SeedFamily, IReadOnlyList<MiRnaDto> Members);

// ================================
// MiRnaAnalyzer — tool result wrappers
// ================================

/// <summary>Result of mirna_seed_sequence.</summary>
public record MiRnaSeedResult(string Seed);

/// <summary>Result of create_mirna.</summary>
public record CreateMiRnaResult(MiRnaDto MiRna);

/// <summary>Result of compare_seed_regions.</summary>
public record CompareSeedRegionsResult(int Matches, int Mismatches, bool IsSameFamily);

/// <summary>Result of find_mirna_target_sites.</summary>
public record FindMiRnaTargetSitesResult(IReadOnlyList<TargetSiteDto> Sites);

/// <summary>Result of rna_reverse_complement.</summary>
public record RnaReverseComplementResult(string ReverseComplement);

/// <summary>Result of can_pair.</summary>
public record CanPairResult(bool CanPair);

/// <summary>Result of is_wobble_pair.</summary>
public record IsWobblePairResult(bool IsWobble);

/// <summary>Result of align_mirna_to_target.</summary>
public record AlignMiRnaToTargetResult(
    string MiRnaSequence,
    string TargetSequence,
    string AlignmentString,
    int Matches,
    int Mismatches,
    int GuWobbles,
    int Gaps,
    double FreeEnergy);

/// <summary>Result of find_pre_mirna_hairpins.</summary>
public record FindPreMiRnaHairpinsResult(IReadOnlyList<PreMiRnaHairpinDto> Hairpins);

/// <summary>Result of calculate_randfold_pvalue. <c>Found</c> is false (P = 1) for an empty sequence.</summary>
public record CalculateRandfoldPValueResult(bool Found, double FreeEnergy, double PValue, int Iterations, int LowerOrEqualCount, int Seed);

/// <summary>Features fed to the in-repo pre-miRNA logistic model.</summary>
public record PreMiRnaFeaturesDto(
    double FreeEnergy, double Amfe, double Mfei, double GcContent, double PairedFraction,
    int StemBasePairs, int LoopSize, int Length, double BasePairingPropensity);

/// <summary>Result of classify_pre_mirna. <c>Found</c> is false (<c>Features</c> null) for an empty sequence.</summary>
public record ClassifyPreMiRnaResult(bool Found, PreMiRnaFeaturesDto? Features, double NaturalProbability, bool IsNatural);

/// <summary>One stem-loop scored by the Triplet-SVM model (1-based start/end in the query).</summary>
public record TripletSvmCandidateDto(
    int Start, int End, string Sequence, string DotBracket, IReadOnlyList<double> Features, double DecisionValue, bool IsPreMiRna);

/// <summary>Result of classify_pre_mirna_triplet_svm. <c>Found</c> is false when no stem-loop was retrieved.</summary>
public record ClassifyPreMiRnaTripletSvmResult(bool Found, IReadOnlyList<TripletSvmCandidateDto> Candidates);

/// <summary>
/// Result of classify_pre_mirna_mipred. <c>Found</c> is false when nothing remains of the sequence. When the package
/// pre-judgement fails (<c>IsHairpinLike</c> false) <c>PValue</c>, <c>Features</c> and <c>RealProbability</c> are null.
/// </summary>
public record ClassifyPreMiRnaMiPredResult(
    bool Found, string Sequence, string DotBracket, double FreeEnergy,
    bool IsHairpinLike, int BasePairs, bool HasMultiloop, int Length,
    double? PValue, IReadOnlyList<double>? Features, int RealVotes, int TreeCount,
    double? RealProbability, bool IsPreMiRna, int Seed, int ShuffleCount);

/// <summary>Result of compute_mipred_features: the 34 features V1..V34.</summary>
public record ComputeMiPredFeaturesResult(IReadOnlyList<double> Features);

/// <summary>Hairpin assessed from its real MFE structure (Ambros 2003 / Bartel 2004 / Zhang 2006); 0-based inclusive span in the scanned sequence.</summary>
public record PreMiRnaMfeDto(
    int Start, int End, string Sequence, string DotBracket, int StemBasePairs,
    int TerminalLoopStart, int TerminalLoopSize, double FreeEnergy, double Amfe, double Mfei);

/// <summary>Result of assess_hairpin_by_mfe. <c>Found</c> is false (<c>Hairpin</c> null) when the candidate is rejected.</summary>
public record AssessHairpinByMfeResult(bool Found, PreMiRnaMfeDto? Hairpin);

/// <summary>Result of find_pre_mirna_hairpins_by_mfe.</summary>
public record FindPreMiRnaHairpinsByMfeResult(IReadOnlyList<PreMiRnaMfeDto> Hairpins);

/// <summary>Result of calculate_mfe_index.</summary>
public record CalculateMfeIndexResult(double Mfei);

/// <summary>Result of dinucleotide_shuffle: <c>Count</c> shuffles drawn in order from one System.Random(<c>Seed</c>).</summary>
public record DinucleotideShuffleResult(IReadOnlyList<string> Shuffles, int Seed, int Count);

/// <summary>Result of score_pre_mirna_features.</summary>
public record ScorePreMiRnaFeaturesResult(double NaturalProbability);

/// <summary>Result of compute_triplet_features. <c>Found</c> is false (<c>Features</c> null) when the structure is not a single arm pair or the sequence has non-ACGU characters.</summary>
public record ComputeTripletFeaturesResult(bool Found, IReadOnlyList<double>? Features);

/// <summary>Result of triplet_svm_decision_value.</summary>
public record TripletSvmDecisionValueResult(double DecisionValue, bool IsPreMiRna);

/// <summary>Base pair of a miRNA·target duplex (0-based indices, miRNA 5'→3' and target 5'→3').</summary>
public record DuplexBasePairDto(int MiRnaIndex, int TargetIndex);

/// <summary>Result of fold_mirna_duplex. <c>Found</c> is false (and <c>Duplex</c> null) when no base pair is possible.</summary>
public record FoldMiRnaDuplexResult(bool Found, MiRnaDuplexStructureDto? Duplex);

/// <summary>Minimum-free-energy miRNA·target duplex (RNAduplex model, Turner 2004).</summary>
public record MiRnaDuplexStructureDto(
    string MiRnaSequence,
    string TargetSequence,
    int MiRnaStart,
    int MiRnaEnd,
    int TargetStart,
    int TargetEnd,
    string Structure,
    IReadOnlyList<DuplexBasePairDto> BasePairs,
    double FreeEnergy);

/// <summary>Result of predict_drosha_dicer_cleavage. <c>Predicted</c> is false (and <c>Cleavage</c> null) when the stem or junction does not allow a prediction.</summary>
public record PredictDroshaDicerCleavageResult(bool Predicted, DroshaDicerCleavageDto? Cleavage);

/// <summary>Predicted Drosha/Dicer cut coordinates (0-based, inclusive spans) and 5p/3p products.</summary>
public record DroshaDicerCleavageDto(
    int BasalJunction,
    int ApicalJunction,
    int DroshaCut5Prime,
    int DroshaCut3Prime,
    int MatureStart,
    int MatureEnd,
    int StarStart,
    int StarEnd,
    string MatureSequence,
    string StarSequence,
    string Sequence,
    string Structure,
    int ThreePrimeOverhang,
    bool HasCnncMotif,
    bool HasBasalUgMotif,
    string DicerRule,
    string DroshaRule);

/// <summary>Result of find_centered_sites.</summary>
public record FindCenteredSitesResult(IReadOnlyList<TargetSiteDto> Sites);

/// <summary>TargetScan 5.0 / Grimson 2007 context score with its contributions (log2 units).</summary>
public record GrimsonContextScoreDto(
    string Type,
    double SiteTypeContribution,
    double ThreePrimePairingContribution,
    double LocalAuContribution,
    double PositionContribution,
    double ThreePrimePairingScore,
    double LocalAuFraction,
    double? ContextScore);

/// <summary>Result of score_target_site_grimson. <c>Scored</c> is false (and <c>Score</c> null) for site types TargetScan 5.0 has no coefficients for.</summary>
public record ScoreTargetSiteGrimsonResult(bool Scored, GrimsonContextScoreDto? Score);

/// <summary>Optional multi-species conservation input of score_target_site_context_plus_plus (tree + species with the site; PCT sigmoid row optional).</summary>
public record ContextConservationInputDto(
    string? Newick,
    int? TargetScanBin,
    IReadOnlyList<string> SpeciesWithSite,
    PctSigmoidParametersDto? SigmoidParameters);

/// <summary>PCT sigmoid parameters PCT(Bls) = b0 + b1 / (1 + exp(-b2*Bls + b3)).</summary>
public record PctSigmoidParametersDto(double B0, double B1, double B2, double B3);

/// <summary>TargetScan 7 context++ (Agarwal 2015) per-feature contributions of one site.</summary>
public record ContextPlusPlusScoreDto(
    string SiteType,
    double Intercept,
    double LocalAuContribution,
    double SRna1Contribution,
    double SRna8Contribution,
    double Site8Contribution,
    double SaContribution,
    double ThreePrimePairingContribution,
    double MinDistContribution,
    double Len3UtrContribution,
    double Off6mContribution,
    double SpsContribution,
    double TaContribution,
    double LenOrfContribution,
    double Orf8mContribution,
    double PctContribution,
    double BranchLengthScore,
    double Pct,
    bool PctIsNa,
    double ContextScorePartial,
    double? ContextScore,
    IReadOnlyList<string> OmittedFeatures);

/// <summary>Result of score_target_site_context_plus_plus.</summary>
public record ScoreTargetSiteContextPlusPlusResult(ContextPlusPlusScoreDto Score);

/// <summary>Result of compute_ta_3utr: the site count N and TA_3UTR = log10(N) (0 when N = 0).</summary>
public record ComputeTa3UtrResult(long SiteCount, double Ta, int UtrCount);

/// <summary>One seed-match site of a TargetScan 3'UTR alignment (1-based coordinates, as targetscan_70.pl).</summary>
public record TargetScanAlignedSiteDto(
    string Species,
    int MsaStart,
    int MsaEnd,
    int UtrStart,
    int UtrEnd,
    string SiteType,
    bool OrfOverlap);

/// <summary>Species of a site group that have a given site type (or a type containing it).</summary>
public record TargetScanSiteTypeSpeciesDto(string SiteType, IReadOnlyList<string> Species);

/// <summary>A TargetScan site group: overlapping seed-match sites across species.</summary>
public record TargetScanSiteGroupDto(
    int GroupNumber,
    string GroupType,
    IReadOnlyList<string> SiteTypes,
    IReadOnlyList<string> Species,
    IReadOnlyList<TargetScanAlignedSiteDto> Sites,
    IReadOnlyList<TargetScanSiteTypeSpeciesDto> SpeciesBySiteType);

/// <summary>Result of find_targetscan_site_groups.</summary>
public record FindTargetScanSiteGroupsResult(IReadOnlyList<TargetScanSiteGroupDto> Groups);

/// <summary>One site with TargetScan branch length, PCT (null = NA) and conservation call.</summary>
public record TargetScanSiteConservationDto(
    TargetScanAlignedSiteDto Site,
    int GroupNumber,
    IReadOnlyList<string> SpeciesWithSite,
    double BranchLength,
    double? Pct,
    bool Conserved);

/// <summary>Result of score_targetscan_conservation. <c>MedianBranchLength</c> is null when the bin was supplied by the caller.</summary>
public record ScoreTargetScanConservationResult(
    int UtrBin,
    double? MedianBranchLength,
    IReadOnlyList<TargetScanSiteGroupDto> Groups,
    IReadOnlyList<TargetScanSiteConservationDto> Sites);

/// <summary>Result of analyze_target_context.</summary>
public record AnalyzeTargetContextResult(
    double AuContent,
    bool NearStart,
    bool NearEnd,
    double ContextScore);

/// <summary>Result of site_accessibility.</summary>
public record SiteAccessibilityResult(double Accessibility);

/// <summary>Result of group_by_seed_family.</summary>
public record GroupBySeedFamilyResult(IReadOnlyList<SeedFamilyDto> Families);

/// <summary>Result of find_similar_mirnas.</summary>
public record FindSimilarMiRnasResult(IReadOnlyList<MiRnaDto> Matches);

/// <summary>Result of generate_seed_variants.</summary>
public record GenerateSeedVariantsResult(IReadOnlyList<string> Variants);

// ================================
// SpliceSitePredictor — shared DTOs
// ================================

/// <summary>Predicted splice site (donor, acceptor, branch point, or U12 variant).</summary>
public record SpliceSiteDto(
    int Position,
    string Type,
    string Motif,
    double Score,
    double Confidence);

/// <summary>Predicted intron with its flanking splice sites and optional branch point.</summary>
public record IntronDto(
    int Start,
    int End,
    int Length,
    SpliceSiteDto DonorSite,
    SpliceSiteDto AcceptorSite,
    SpliceSiteDto? BranchPoint,
    string Sequence,
    string Type,
    double Score);

/// <summary>Predicted exon segment.</summary>
public record ExonDto(
    int Start,
    int End,
    int Length,
    string Type,
    int? Phase,
    string Sequence);

/// <summary>Candidate alternative-splicing event.</summary>
public record AlternativeSplicingEventDto(
    string Type,
    int Position,
    string Description);

// ================================
// SpliceSitePredictor — tool result wrappers
// ================================

/// <summary>Result of find_donor_sites.</summary>
public record FindDonorSitesResult(IReadOnlyList<SpliceSiteDto> Sites);

/// <summary>Result of find_acceptor_sites.</summary>
public record FindAcceptorSitesResult(IReadOnlyList<SpliceSiteDto> Sites);

/// <summary>Result of find_branch_points.</summary>
public record FindBranchPointsResult(IReadOnlyList<SpliceSiteDto> Sites);

/// <summary>Result of predict_introns.</summary>
public record PredictIntronsResult(IReadOnlyList<IntronDto> Introns);

/// <summary>Result of predict_gene_structure.</summary>
public record PredictGeneStructureResult(
    IReadOnlyList<ExonDto> Exons,
    IReadOnlyList<IntronDto> Introns,
    string SplicedSequence,
    double OverallScore,
    int? CdsStart = null,
    string? CodingSequence = null,
    double? IntronScore = null);

/// <summary>Result of detect_alternative_splicing.</summary>
public record DetectAlternativeSplicingResult(IReadOnlyList<AlternativeSplicingEventDto> Events);

/// <summary>Result of find_retained_intron_candidates.</summary>
public record FindRetainedIntronCandidatesResult(IReadOnlyList<IntronDto> Introns);

/// <summary>Result of maxent_score.</summary>
public record MaxEntScoreResult(double Score);

/// <summary>Result of is_within_coding_region.</summary>
public record IsWithinCodingRegionResult(bool IsCoding);

// ================================
// StructuralVariantAnalyzer — shared DTOs
// ================================

/// <summary>A structural variant call.</summary>
public record StructuralVariantDto(
    string Id,
    string Chromosome,
    int Start,
    int End,
    string Type,
    int Length,
    double Quality,
    int SupportingReads,
    string? InsertedSequence);

/// <summary>A breakpoint joining two genomic loci.</summary>
public record BreakpointDto(
    string Chromosome1,
    int Position1,
    char Strand1,
    string Chromosome2,
    int Position2,
    char Strand2,
    int SupportingReads,
    double Quality);

/// <summary>A copy-number segment.</summary>
public record CopyNumberSegmentDto(
    string Chromosome,
    int Start,
    int End,
    double LogRatio,
    int CopyNumber,
    double BAlleleFrequency,
    int ProbeCount);

/// <summary>A read-pair signature (potentially discordant).</summary>
public record ReadPairSignatureDto(
    string ReadId,
    string Chromosome1,
    int Position1,
    char Strand1,
    string Chromosome2,
    int Position2,
    char Strand2,
    int InsertSize,
    bool IsDiscordant);

/// <summary>A split-read alignment record.</summary>
public record SplitReadDto(
    string ReadId,
    string Chromosome,
    int PrimaryPosition,
    int SupplementaryPosition,
    int ClipLength,
    string ClippedSequence);

/// <summary>SV annotation against gene/exon models.</summary>
public record SvAnnotationDto(
    string SvId,
    IReadOnlyList<string> AffectedGenes,
    IReadOnlyList<string> AffectedExons,
    string FunctionalImpact,
    double PopulationFrequency,
    bool IsPathogenic);

// ================================
// StructuralVariantAnalyzer — input DTOs
// ================================

/// <summary>A read-pair input for discordant-pair detection.</summary>
public record ReadPairInputDto(
    string ReadId,
    string Chr1,
    int Pos1,
    char Strand1,
    string Chr2,
    int Pos2,
    char Strand2,
    int InsertSize);

/// <summary>An alignment record (CIGAR + sequence) for split-read detection.</summary>
public record AlignmentInputDto(
    string ReadId,
    string Chromosome,
    int Position,
    string Cigar,
    string Sequence);

/// <summary>A copy-number probe observation.</summary>
public record CopyNumberProbeDto(
    string Chromosome,
    int Position,
    double LogRatio,
    double Baf);

/// <summary>An exon interval (used by SV gene-overlap input).</summary>
public record SvExonIntervalDto(int Start, int End);

/// <summary>A gene model with exons used for SV functional annotation.</summary>
public record SvGeneInputDto(
    string GeneId,
    string Chromosome,
    int Start,
    int End,
    IReadOnlyList<SvExonIntervalDto> Exons);

// ================================
// StructuralVariantAnalyzer — tool result wrappers
// ================================

/// <summary>Result of find_discordant_pairs.</summary>
public record FindDiscordantPairsResult(IReadOnlyList<ReadPairSignatureDto> Pairs);

/// <summary>Result of cluster_discordant_pairs.</summary>
public record ClusterDiscordantPairsResult(IReadOnlyList<StructuralVariantDto> Variants);

/// <summary>Result of find_split_reads.</summary>
public record FindSplitReadsResult(IReadOnlyList<SplitReadDto> Reads);

/// <summary>Result of cluster_split_reads.</summary>
public record ClusterSplitReadsResult(IReadOnlyList<BreakpointDto> Breakpoints);

/// <summary>Result of segment_copy_number.</summary>
public record SegmentCopyNumberResult(IReadOnlyList<CopyNumberSegmentDto> Segments);

/// <summary>Result of identify_cnvs.</summary>
public record IdentifyCnvsResult(IReadOnlyList<StructuralVariantDto> Variants);

/// <summary>Result of merge_overlapping_svs.</summary>
public record MergeOverlappingSvsResult(IReadOnlyList<StructuralVariantDto> Merged);

/// <summary>Result of filter_svs.</summary>
public record FilterSvsResult(IReadOnlyList<StructuralVariantDto> Variants);

/// <summary>Result of annotate_svs.</summary>
public record AnnotateSvsResult(IReadOnlyList<SvAnnotationDto> Annotations);

/// <summary>Result of genotype_sv.</summary>
public record GenotypeSvResult(string Genotype, double Quality);

/// <summary>Result of assemble_breakpoint_sequence.</summary>
public record AssembleBreakpointSequenceResult(string? Sequence);

/// <summary>Result of find_microhomology.</summary>
public record FindMicrohomologyResult(int MicrohomologyLength, string Sequence);

// ================================
// TranscriptomeAnalyzer — shared / input DTOs
// ================================

/// <summary>Raw count input for TPM/FPKM quantification.</summary>
public record GeneCountInputDto(string GeneId, double RawCount, int Length);

/// <summary>Per-gene expression record (TPM, FPKM, raw count, length).</summary>
public record GeneExpressionDto(
    string GeneId,
    double RawCount,
    double Tpm,
    double Fpkm,
    int Length);

/// <summary>One row of a differential-expression input matrix.</summary>
public record DiffExprInputDto(
    string GeneId,
    IReadOnlyList<double> Group1,
    IReadOnlyList<double> Group2);

/// <summary>Differential-expression result for a single gene.</summary>
public record DifferentialExpressionDto(
    string GeneId,
    double Log2FoldChange,
    double PValue,
    double AdjustedPValue,
    bool IsSignificant,
    string Regulation);

/// <summary>Pathway / gene-set definition for ORA.</summary>
public record PathwayInputDto(
    string PathwayId,
    string PathwayName,
    IReadOnlyList<string> Genes);

/// <summary>Gene-set / pathway enrichment result.</summary>
public record EnrichmentResultDto(
    string PathwayId,
    string PathwayName,
    int GenesInPathway,
    int OverlappingGenes,
    double EnrichmentScore,
    double PValue,
    IReadOnlyList<string> Genes);

/// <summary>Inclusion/skipping read counts for a candidate skipped exon.</summary>
public record SkippedExonInputDto(
    string GeneId,
    int ExonStart,
    int ExonEnd,
    double InclusionReads,
    double SkippingReads);

/// <summary>Two-condition PSI input for differential splicing.</summary>
public record SplicingComparisonInputDto(
    string GeneId,
    int Start,
    int End,
    double PsiCondition1,
    double PsiCondition2);

/// <summary>Splicing event with PSI / delta-PSI.</summary>
public record SplicingEventDto(
    string GeneId,
    string EventType,
    int Start,
    int End,
    double InclusionLevel,
    double DeltaPsi);

/// <summary>Exon coordinate range used by transcript isoforms.</summary>
public record ExonRangeDto(int Start, int End);

/// <summary>Transcript isoform descriptor.</summary>
public record TranscriptIsoformDto(
    string TranscriptId,
    string GeneId,
    int Length,
    int ExonCount,
    double Expression,
    bool IsProteinCoding,
    IReadOnlyList<ExonRangeDto> Exons);

/// <summary>Dominant isoform per gene with its dominance ratio.</summary>
public record DominantIsoformDto(
    string GeneId,
    TranscriptIsoformDto DominantIsoform,
    double DominanceRatio);

/// <summary>Per-isoform expression in two conditions for switching analysis.</summary>
public record IsoformExpressionInputDto(
    TranscriptIsoformDto Isoform,
    double Expression1,
    double Expression2);

/// <summary>Detected isoform-usage switch between two conditions.</summary>
public record IsoformSwitchDto(
    string GeneId,
    string TranscriptId1,
    string TranscriptId2,
    double SwitchScore);

/// <summary>Gene expression profile across samples / conditions.</summary>
public record GeneProfileInputDto(
    string GeneId,
    IReadOnlyList<double> Expression);

/// <summary>Edge of a co-expression network.</summary>
public record CoExpressionEdgeDto(string Gene1, string Gene2, double Correlation);

/// <summary>Co-expression cluster of genes.</summary>
public record CoExpressionClusterDto(
    int ClusterId,
    IReadOnlyList<string> Genes,
    double MeanCorrelation,
    string RepresentativeGene,
    IReadOnlyList<string> EnrichedFunctions);

/// <summary>Per-sample expression vector for PCA.</summary>
public record SampleExpressionInputDto(
    string SampleId,
    IReadOnlyList<double> Expression);

/// <summary>Sample projection on the first two principal components.</summary>
public record PcaPointDto(string SampleId, double Pc1, double Pc2);

// ================================
// TranscriptomeAnalyzer — tool result wrappers
// ================================

/// <summary>Result of calculate_tpm.</summary>
public record CalculateTpmResult(IReadOnlyList<GeneExpressionDto> Expressions);

/// <summary>Result of quantile_normalize.</summary>
public record QuantileNormalizeResult(IReadOnlyList<IReadOnlyList<double>> Normalized);

/// <summary>Result of log2_transform.</summary>
public record Log2TransformResult(IReadOnlyList<double> Transformed);

/// <summary>Result of differential_expression.</summary>
public record DifferentialExpressionResult(IReadOnlyList<DifferentialExpressionDto> Results);

/// <summary>Result of over_representation_analysis.</summary>
public record OverRepresentationAnalysisResult(IReadOnlyList<EnrichmentResultDto> Results);

/// <summary>Result of enrichment_score.</summary>
public record EnrichmentScoreResult(double Score);

/// <summary>Result of find_skipped_exon_events.</summary>
public record FindSkippedExonEventsResult(IReadOnlyList<SplicingEventDto> Events);

/// <summary>Result of detect_differential_splicing.</summary>
public record DetectDifferentialSplicingResult(IReadOnlyList<SplicingEventDto> Events);

/// <summary>Result of find_dominant_isoforms.</summary>
public record FindDominantIsoformsResult(IReadOnlyList<DominantIsoformDto> Dominants);

/// <summary>Result of detect_isoform_switching.</summary>
public record DetectIsoformSwitchingResult(IReadOnlyList<IsoformSwitchDto> Switches);

/// <summary>Result of pearson_correlation.</summary>
public record PearsonCorrelationResult(double Correlation);

/// <summary>Result of build_coexpression_network.</summary>
public record BuildCoexpressionNetworkResult(IReadOnlyList<CoExpressionEdgeDto> Edges);

/// <summary>Result of cluster_genes_by_expression.</summary>
public record ClusterGenesByExpressionResult(IReadOnlyList<CoExpressionClusterDto> Clusters);

/// <summary>Result of rnaseq_quality_metrics.</summary>
public record RnaseqQualityMetricsResult(
    double MappingRate,
    double ExonicRate,
    double RRnaRate,
    int DetectedGenes);

/// <summary>Result of perform_pca.</summary>
public record PerformPcaResult(IReadOnlyList<PcaPointDto> Points);

/// <summary>Result of assign_targetscan_utr_bin: bin 1-10, the median per-position branch length (null when the reference row has no nucleotide) and the 10 bin lower bounds (bundled generic tree).</summary>
public record AssignTargetScanUtrBinResult(int Bin, double? MedianBranchLength, IReadOnlyList<double> BinThresholds);

/// <summary>Result of compute_branch_length_score: Friedman 2009 branch-length score (4 decimals) and the tree source used.</summary>
public record ComputeBranchLengthScoreResult(double BranchLengthScore, string TreeSource);

/// <summary>Result of pct_from_branch_length. <c>Pct</c> is null (TargetScan "NA") when the seed has no bundled PCT family row; <c>Parameters</c> are the sigmoid row used (null when none).</summary>
public record PctFromBranchLengthResult(double? Pct, bool IsNa, PctSigmoidParametersDto? Parameters);

/// <summary>Result of is_conserved_by_branch_length: TargetScan "Conserved" call and the threshold applied (null for site types that are never conserved).</summary>
public record IsConservedByBranchLengthResult(bool Conserved, double? Threshold);

/// <summary>Result of get_targetscan_parameters: PCT family membership, PCT sigmoid row for the site type, and the bundled TA / SPS row.</summary>
public record GetTargetScanParametersResult(
    string SeedRegion,
    bool IsPctSeedFamily,
    PctSigmoidParametersDto? PctParameters,
    bool HasTaSps,
    double? Ta,
    double? Sps);

/// <summary>Result of count_targetscan_orf_8mer_sites: ORF8m of targetscan_count_8mers.pl.</summary>
public record CountTargetScanOrf8merSitesResult(int Count);
