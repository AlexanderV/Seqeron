namespace Seqeron.Mcp.Analysis.Tools;

// ================================
// KmerAnalyzer Results
// ================================

/// <summary>k-mer count map (k-mer → occurrence count).</summary>
public record KmerCountsResult(Dictionary<string, int> Counts);

/// <summary>k-mer frequency-of-frequencies spectrum.</summary>
public record KmerSpectrumResult(Dictionary<int, int> Spectrum);

/// <summary>List of k-mers.</summary>
public record KmerListResult(string[] Kmers);

/// <summary>Normalized k-mer frequencies.</summary>
public record KmerFrequenciesResult(Dictionary<string, double> Frequencies);

/// <summary>Euclidean distance between two k-mer frequency vectors.</summary>
public record KmerDistanceResult(double Distance);

/// <summary>k-mers paired with their occurrence count.</summary>
public record KmerCountItem(string Kmer, int Count);

/// <summary>Result of <c>kmers_with_min_count</c>.</summary>
public record KmersWithMinCountResult(KmerCountItem[] Items);

/// <summary>Positions of all (overlapping) k-mer occurrences.</summary>
public record KmerPositionsResult(int[] Positions);

/// <summary>Aggregate k-mer statistics.</summary>
public record AnalyzeKmersResult(
    int TotalKmers,
    int UniqueKmers,
    int MaxCount,
    int MinCount,
    double AverageCount,
    double Entropy);

// ================================
// SequenceStatistics Results
// ================================

/// <summary>Sliding-window double-valued profile.</summary>
public record DoubleProfileResult(double[] Values);

/// <summary>Dinucleotide frequency map.</summary>
public record DinucleotideFrequenciesResult(Dictionary<string, double> Frequencies);

/// <summary>Observed/expected dinucleotide ratios.</summary>
public record DinucleotideRatiosResult(Dictionary<string, double> Ratios);

/// <summary>Codon frequency map.</summary>
public record CodonFrequenciesResult(Dictionary<string, double> Frequencies);

/// <summary>Per-window Chou-Fasman propensities.</summary>
public record ChouFasmanItem(double Helix, double Sheet, double Turn);

/// <summary>Result of <c>predict_chou_fasman</c>.</summary>
public record PredictChouFasmanResult(ChouFasmanItem[] Items);

// ================================
// GenomicAnalyzer Results
// ================================

/// <summary>A repeated region with all its occurrence positions.</summary>
public record RepeatItem(string Sequence, int[] Positions, int Length, int Count);

/// <summary>Result of <c>find_repeats</c>.</summary>
public record FindRepeatsResult(RepeatItem[] Items);

/// <summary>A consecutive tandem repeat.</summary>
public record TandemRepeatItem(string Unit, int Position, int Repetitions, int TotalLength);

/// <summary>Result of <c>find_tandem_repeats</c>.</summary>
public record FindTandemRepeatsResult(TandemRepeatItem[] Items);

/// <summary>Positions of motif occurrences.</summary>
public record MotifPositionsResult(int[] Positions);

/// <summary>Map of motif → matching positions for hits with at least one match.</summary>
public record FindKnownMotifsResult(Dictionary<string, int[]> Matches);

/// <summary>A common region between two sequences.</summary>
public record CommonRegionItem(string Sequence, int PositionInFirst, int PositionInSecond, int Length);

/// <summary>Result of <c>find_common_regions</c>.</summary>
public record FindCommonRegionsResult(CommonRegionItem[] Items);

/// <summary>An open reading frame.</summary>
public record OrfItem(
    string Sequence,
    int Position,
    int Frame,
    bool IsReverseComplement,
    int Length,
    int CodonCount);

/// <summary>Result of <c>find_open_reading_frames</c>.</summary>
public record FindOpenReadingFramesResult(OrfItem[] Items);

// ================================
// RepeatFinder Results
// ================================

/// <summary>A microsatellite (STR) hit.</summary>
public record MicrosatelliteItem(
    int Position,
    string RepeatUnit,
    int RepeatCount,
    int TotalLength,
    string RepeatType);

/// <summary>Result of <c>find_microsatellites</c>.</summary>
public record FindMicrosatellitesResult(MicrosatelliteItem[] Items);

/// <summary>An inverted repeat candidate (potential hairpin).</summary>
public record InvertedRepeatItem(
    int LeftArmStart,
    int RightArmStart,
    int ArmLength,
    int LoopLength,
    string LeftArm,
    string RightArm,
    string Loop,
    bool CanFormHairpin,
    int TotalLength);

/// <summary>Result of <c>find_inverted_repeats</c>.</summary>
public record FindInvertedRepeatsResult(InvertedRepeatItem[] Items);

/// <summary>A direct repeat (two identical occurrences with a spacer).</summary>
public record DirectRepeatItem(
    int FirstPosition,
    int SecondPosition,
    string RepeatSequence,
    int Length,
    int Spacing);

/// <summary>Result of <c>find_direct_repeats</c>.</summary>
public record FindDirectRepeatsResult(DirectRepeatItem[] Items);

/// <summary>A palindrome (sequence equal to its reverse complement).</summary>
public record PalindromeItem(int Position, string Sequence, int Length);

/// <summary>Result of <c>find_palindromes</c>.</summary>
public record FindPalindromesResult(PalindromeItem[] Items);

/// <summary>Aggregate statistics across all microsatellites.</summary>
public record TandemRepeatSummaryResult(
    int TotalRepeats,
    int TotalRepeatBases,
    double PercentageOfSequence,
    int MononucleotideRepeats,
    int DinucleotideRepeats,
    int TrinucleotideRepeats,
    int TetranucleotideRepeats,
    MicrosatelliteItem? LongestRepeat,
    string? MostFrequentUnit);

// ================================
// MotifFinder Results
// ================================

/// <summary>A motif match (position, matched substring, pattern, score).</summary>
public record MotifMatchItem(int Position, string MatchedSequence, string Pattern, double Score);

/// <summary>Result of <c>find_degenerate_motif</c>.</summary>
public record FindDegenerateMotifResult(MotifMatchItem[] Items);

/// <summary>PWM input: jagged 4×L matrix (rows A,C,G,T) and length L.</summary>
public record PwmInput(double[][] Matrix, int Length);

/// <summary>Result of <c>create_pwm</c>: log-odds 4×L matrix plus consensus and score bounds.</summary>
public record PwmResult(double[][] Matrix, int Length, string Consensus, double MaxScore, double MinScore);

/// <summary>Result of <c>scan_with_pwm</c>.</summary>
public record ScanWithPwmResult(MotifMatchItem[] Items);

/// <summary>Result of <c>generate_consensus</c>.</summary>
public record ConsensusResult(string Consensus);

/// <summary>A discovered (overrepresented) k-mer motif.</summary>
public record DiscoveredMotifItem(string Sequence, int Count, int[] Positions, double Enrichment);

/// <summary>Result of <c>discover_motifs</c>.</summary>
public record DiscoverMotifsResult(DiscoveredMotifItem[] Items);

/// <summary>A motif shared between sequences.</summary>
public record SharedMotifItem(string Sequence, int[] SequenceIndices, double Prevalence);

/// <summary>Result of <c>find_shared_motifs</c>.</summary>
public record FindSharedMotifsResult(SharedMotifItem[] Items);

/// <summary>A regulatory element occurrence.</summary>
public record RegulatoryElementItem(string Name, int Position, string Sequence, string Pattern, string Description);

/// <summary>Result of <c>find_regulatory_elements</c>.</summary>
public record FindRegulatoryElementsResult(RegulatoryElementItem[] Items);

// ================================
// ProteinMotifFinder Results
// ================================

/// <summary>A protein motif match.</summary>
public record ProteinMotifMatchItem(
    int Start,
    int End,
    string Sequence,
    string MotifName,
    string Pattern,
    double Score,
    double EValue);

/// <summary>Result of <c>find_protein_motifs</c>, <c>find_motif_by_pattern</c>, and <c>find_motif_by_prosite</c>.</summary>
public record FindProteinMotifsResult(ProteinMotifMatchItem[] Items);

/// <summary>Result of <c>prosite_to_regex</c>.</summary>
public record PrositeRegexResult(string Regex);

/// <summary>Signal-peptide cleavage-site prediction (von Heijne 1986 weight-matrix method).</summary>
public record SignalPeptideResult(
    bool Found,
    int CleavagePosition,
    double Score,
    string SignalSequence,
    string WindowSequence,
    bool IsLikelySignalPeptide);

/// <summary>A scored region (start, end, score) — used for TM helices and coiled-coils.</summary>
public record RegionScoreItem(int Start, int End, double Score);

/// <summary>Result of <c>predict_transmembrane_helices</c>.</summary>
public record PredictTransmembraneHelicesResult(RegionScoreItem[] Items);

/// <summary>Result of <c>predict_coiled_coils</c>.</summary>
public record PredictCoiledCoilsResult(RegionScoreItem[] Items);

/// <summary>A protein low-complexity region (SEG): 0-based inclusive span and minimum window complexity in bits/residue.</summary>
public record ProteinLowComplexityItem(int Start, int End, double Complexity);

/// <summary>Result of <c>find_protein_low_complexity_regions</c>.</summary>
public record FindProteinLowComplexityRegionsResult(ProteinLowComplexityItem[] Items);

/// <summary>A protein domain hit.</summary>
public record ProteinDomainItem(
    string Name,
    string Accession,
    int Start,
    int End,
    double Score,
    string Description);

/// <summary>Result of <c>find_protein_domains</c>.</summary>
public record FindProteinDomainsResult(ProteinDomainItem[] Items);

// ================================
// SequenceComplexity Results
// ================================

/// <summary>Sliding-window complexity point.</summary>
public record ComplexityPointItem(
    int Position,
    double ShannonEntropy,
    double LinguisticComplexity,
    int WindowStart,
    int WindowEnd);

/// <summary>Result of <c>windowed_complexity</c>.</summary>
public record WindowedComplexityResult(ComplexityPointItem[] Items);

/// <summary>A low-complexity DNA region (entropy-based).</summary>
public record DnaLowComplexityItem(int Start, int End, int Length, double MinEntropy, string Sequence);

/// <summary>Result of <c>find_low_complexity_regions</c> (DNA).</summary>
public record FindLowComplexityRegionsResult(DnaLowComplexityItem[] Items);

/// <summary>Result of <c>dust_score</c>.</summary>
public record DustScoreResult(double Score);

/// <summary>Result of <c>mask_low_complexity</c>.</summary>
public record MaskLowComplexityResult(string Masked);

/// <summary>Result of <c>compression_ratio</c>.</summary>
public record CompressionRatioResult(double Ratio);

// ================================
// ComparativeGenomics Results
// ================================

/// <summary>Gene input for comparative analyses.</summary>
public record GeneInput(string Id, string GenomeId, int Start, int End, char Strand, string? Sequence = null);

/// <summary>A syntenic block between two genomes.</summary>
public record SyntenicBlockItem(
    string Genome1Id,
    int Start1,
    int End1,
    string Genome2Id,
    int Start2,
    int End2,
    bool IsInverted,
    int GeneCount,
    double Identity);

/// <summary>Result of <c>find_syntenic_blocks</c>.</summary>
public record FindSyntenicBlocksResult(SyntenicBlockItem[] Items);

/// <summary>An ortholog gene pair.</summary>
public record OrthologPairItem(
    string Gene1Id,
    string Gene2Id,
    double Identity,
    double Coverage,
    int AlignmentLength);

/// <summary>Result of <c>find_orthologs</c>.</summary>
public record FindOrthologsResult(OrthologPairItem[] Items);

/// <summary>Result of <c>find_reciprocal_best_hits</c>.</summary>
public record FindReciprocalBestHitsResult(OrthologPairItem[] Items);

/// <summary>A genome rearrangement event (type as string: Inversion, Translocation, Deletion, Insertion, Duplication, Transposition).</summary>
public record RearrangementEventItem(
    string Type,
    string GenomeId,
    int Position,
    int Length,
    string? TargetPosition);

/// <summary>Result of <c>detect_rearrangements</c>.</summary>
public record DetectRearrangementsResult(RearrangementEventItem[] Items);

/// <summary>Result of <c>compare_genomes</c>.</summary>
public record CompareGenomesResult(
    SyntenicBlockItem[] SyntenicBlocks,
    OrthologPairItem[] Orthologs,
    RearrangementEventItem[] Rearrangements,
    double OverallSynteny,
    int ConservedGenes,
    int GenomeSpecificGenes1,
    int GenomeSpecificGenes2);

/// <summary>Result of <c>reversal_distance</c>.</summary>
public record ReversalDistanceResult(int Distance);

/// <summary>Result of <c>find_conserved_clusters</c>: each cluster is a list of ortholog-group IDs.</summary>
public record FindConservedClustersResult(string[][] Clusters);

/// <summary>Result of <c>calculate_ani</c>.</summary>
public record AniResult(double Ani);

/// <summary>A dot-plot match coordinate.</summary>
public record DotPlotPoint(int X, int Y);

/// <summary>Result of <c>generate_dot_plot</c>.</summary>
public record GenerateDotPlotResult(DotPlotPoint[] Points);

// ================================
// DisorderPredictor Results
// ================================

/// <summary>Per-residue disorder prediction.</summary>
public record DisorderResiduePredictionItem(int Position, char Residue, double DisorderScore, bool IsDisordered);

/// <summary>Contiguous disordered region.</summary>
public record DisorderRegionItem(int Start, int End, double MeanScore, double Confidence, string RegionType);

/// <summary>Result of <c>predict_disorder</c>.</summary>
public record PredictDisorderResult(
    string Sequence,
    DisorderResiduePredictionItem[] ResiduePredictions,
    DisorderRegionItem[] DisorderedRegions,
    double OverallDisorderContent,
    double MeanDisorderScore);

/// <summary>SEG low-complexity region item.</summary>
public record SegRegionItem(int Start, int End, string Type);

/// <summary>Result of <c>predict_low_complexity_seg</c>.</summary>
public record PredictLowComplexitySegResult(SegRegionItem[] Items);

/// <summary>Predicted MoRF (Molecular Recognition Feature) item.</summary>
public record MorfItem(int Start, int End, double Score);

/// <summary>Result of <c>predict_morfs</c>.</summary>
public record PredictMorfsResult(MorfItem[] Items);

/// <summary>Result of <c>disorder_propensity</c>.</summary>
public record DisorderPropensityResult(double Propensity);

/// <summary>Result of <c>is_disorder_promoting</c>.</summary>
public record IsDisorderPromotingResult(bool Result);

// ================================
// GcSkewCalculator Results
// ================================

/// <summary>Result of <c>gc_skew</c>.</summary>
public record GcSkewResult(double GcSkew);

/// <summary>Result of <c>at_skew</c>.</summary>
public record AtSkewResult(double AtSkew);

/// <summary>A windowed GC skew point.</summary>
public record GcSkewPointItem(int Position, double GcSkew, int WindowStart, int WindowEnd);

/// <summary>Result of <c>windowed_gc_skew</c>.</summary>
public record WindowedGcSkewResult(GcSkewPointItem[] Items);

/// <summary>A cumulative GC skew point.</summary>
public record CumulativeGcSkewPointItem(int Position, double GcSkew, double CumulativeGcSkew);

/// <summary>Result of <c>cumulative_gc_skew</c>.</summary>
public record CumulativeGcSkewResult(CumulativeGcSkewPointItem[] Items);

/// <summary>A windowed GC content point.</summary>
public record GcContentPointItem(int Position, double GcContent, int WindowStart, int WindowEnd);

/// <summary>Result of <c>predict_replication_origin</c>.</summary>
public record PredictReplicationOriginResult(
    int PredictedOrigin,
    int PredictedTerminus,
    double OriginSkew,
    double TerminusSkew,
    bool IsSignificant);

/// <summary>Result of <c>analyze_gc_content</c>.</summary>
public record AnalyzeGcContentResult(
    double OverallGcContent,
    double OverallGcSkew,
    double OverallAtSkew,
    double GcContentVariance,
    double GcSkewVariance,
    GcSkewPointItem[] WindowedGcSkew,
    GcContentPointItem[] WindowedGcContent,
    int SequenceLength);

// ================================
// RnaSecondaryStructure Results
// ================================

/// <summary>RNA base pair (renamed to avoid colliding with source <c>BasePair</c> nested record).</summary>
public record BasePairItem(int Position1, int Position2, char Base1, char Base2, string Type);

/// <summary>Bare base-pair coordinate (used by <c>parse_dot_bracket</c>).</summary>
public record BasePairCoord(int Position1, int Position2);

/// <summary>RNA stem (helix region).</summary>
public record StemItem(
    int Start5Prime,
    int End5Prime,
    int Start3Prime,
    int End3Prime,
    int Length,
    BasePairItem[] BasePairs,
    double FreeEnergy);

/// <summary>RNA loop region.</summary>
public record LoopItem(string Type, int Start, int End, int Size, string Sequence);

/// <summary>RNA stem-loop (hairpin) item.</summary>
public record StemLoopItem(
    int Start,
    int End,
    StemItem Stem,
    LoopItem Loop,
    double TotalFreeEnergy,
    string DotBracketNotation);

/// <summary>RNA pseudoknot item.</summary>
public record PseudoknotItem(int Start1, int End1, int Start2, int End2, BasePairItem[] CrossingPairs);

/// <summary>Result of <c>can_pair</c>.</summary>
public record CanPairResult(bool Result);

/// <summary>Result of <c>base_pair_type</c>; <c>Type</c> is null if bases cannot pair.</summary>
public record BasePairTypeResult(string? Type);

/// <summary>Result of <c>rna_complement_base</c>.</summary>
public record RnaComplementBaseResult(char Complement);

/// <summary>Result of <c>find_stem_loops</c>.</summary>
public record FindStemLoopsResult(StemLoopItem[] Items);

/// <summary>Result of <c>stem_energy</c>.</summary>
public record StemEnergyResult(double Energy);

/// <summary>Result of <c>terminal_mismatch_energy</c>.</summary>
public record TerminalMismatchEnergyResult(double Energy);

/// <summary>Result of <c>dangling_end_energy</c>.</summary>
public record DanglingEndEnergyResult(double Energy);

/// <summary>Result of <c>hairpin_loop_energy</c>.</summary>
public record HairpinLoopEnergyResult(double Energy);

/// <summary>Result of <c>internal_loop_energy</c>.</summary>
public record InternalLoopEnergyResult(double Energy);

/// <summary>Result of <c>bulge_loop_energy</c>.</summary>
public record BulgeLoopEnergyResult(double Energy);

/// <summary>Result of <c>multibranch_loop_energy</c>.</summary>
public record MultibranchLoopEnergyResult(double Energy);

/// <summary>Result of <c>flush_coaxial_stacking</c>.</summary>
public record FlushCoaxialStackingResult(double Energy);

/// <summary>Result of <c>mismatch_coaxial_stacking</c>.</summary>
public record MismatchCoaxialStackingResult(double Energy);

/// <summary>Result of <c>minimum_free_energy</c>.</summary>
public record MinimumFreeEnergyResult(double Mfe);

/// <summary>Result of <c>predict_rna_structure</c>.</summary>
public record PredictRnaStructureResult(
    string Sequence,
    string DotBracket,
    BasePairItem[] BasePairs,
    StemLoopItem[] StemLoops,
    PseudoknotItem[] Pseudoknots,
    double MinimumFreeEnergy);

/// <summary>Result of <c>detect_pseudoknots</c>.</summary>
public record DetectPseudoknotsResult(PseudoknotItem[] Items);

/// <summary>Result of <c>parse_dot_bracket</c>.</summary>
public record ParseDotBracketResult(BasePairCoord[] Pairs);

/// <summary>Result of <c>validate_dot_bracket</c>.</summary>
public record ValidateDotBracketResult(bool Result);

/// <summary>Inverted-repeat item for RNA structure (antiparallel complementary regions).</summary>
public record FindRnaInvertedRepeatsItem(int Start1, int End1, int Start2, int End2, int Length);

/// <summary>Result of <c>find_rna_inverted_repeats</c>.</summary>
public record FindRnaInvertedRepeatsResult(FindRnaInvertedRepeatsItem[] Items);
