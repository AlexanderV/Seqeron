namespace Seqeron.Mcp.Analysis.Tools;

// ================================
// KmerAnalyzer Results
// ================================

/// <summary>k-mer count map (k-mer → occurrence count).</summary>
public record KmerCountsResult(Dictionary<string, int> Counts);

/// <summary>
/// k-mer frequency-of-frequencies spectrum; <paramref name="Histogram"/> holds the jellyfish-histo rows when any
/// histo option (low/high/increment/full) was given, otherwise null.
/// </summary>
public record KmerSpectrumResult(Dictionary<int, int> Spectrum, KmerHistogramRow[]? Histogram = null);

/// <summary>One jellyfish <c>histo</c> row: bucket label (low end point) and number of distinct k-mers in it.</summary>
public record KmerHistogramRow(long Bin, long Frequency);

/// <summary>List of k-mers.</summary>
public record KmerListResult(string[] Kmers);

/// <summary>Normalized k-mer frequencies.</summary>
public record KmerFrequenciesResult(Dictionary<string, double> Frequencies);

/// <summary>Euclidean distance between two k-mer frequency vectors.</summary>
public record KmerDistanceResult(double Distance);

/// <summary>
/// Result of <c>kmer_d2_statistics</c>: raw D2* / D2S, the dissimilarities d2* / d2S, the Markov orders used and each
/// sequence's BIC for orders 0..min(k−1, 10) (index = order).
/// </summary>
public record KmerD2StatisticsResult(
    double D2Star,
    double D2Shepherd,
    double D2StarDistance,
    double D2ShepherdDistance,
    int MarkovOrder1,
    int MarkovOrder2,
    double[] Bic1,
    double[] Bic2);

/// <summary>
/// k-mer Jaccard index (fraction in [0,1]) and the Mash distance derived from it — exact, or estimated from Mash
/// MinHash sketches when <c>sketchSize</c> &gt; 0 (then <c>SharedHashes</c>/<c>SketchDenominator</c> = Mash "x/s" and
/// <c>PValue</c> = Mash p-value; null in exact mode) — plus the exact containment indices |A∩B|/|A| and |A∩B|/|B|.
/// With <c>scaled</c> &gt; 0 every value comes from sourmash FracMinHash sketches: <c>SharedHashes</c> = |A∩B|,
/// <c>SketchDenominator</c> = |A∪B|, the containments are sourmash <c>contained_by</c> (bias-corrected) and
/// <c>MaxContainment</c> = <c>max_containment</c> (null otherwise); with <c>trackAbundance</c>, <c>AngularSimilarity</c> =
/// sourmash <c>angular_similarity</c> (null otherwise).
/// </summary>
public record KmerJaccardResult(
    double Jaccard,
    double MashDistance,
    double ContainmentSeq1InSeq2 = 0,
    double ContainmentSeq2InSeq1 = 0,
    int? SharedHashes = null,
    int? SketchDenominator = null,
    double? PValue = null,
    double? MaxContainment = null,
    double? AngularSimilarity = null);

/// <summary>A maximal run of consecutive qualifying window starts (0-based, both ends inclusive).</summary>
public record ClumpWindowRunItem(int FirstWindowStart, int LastWindowStart);

/// <summary>An (L, t)-clump k-mer, its leftmost qualifying window start and its maximal window runs.</summary>
public record KmerClumpItem(string Kmer, int FirstWindowStart, ClumpWindowRunItem[] WindowRuns);

/// <summary>Result of <c>find_clump_windows</c>.</summary>
public record KmerClumpWindowsResult(KmerClumpItem[] Clumps);

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
    double Entropy)
{
    /// <summary>Number of distinct k-mers (Jellyfish "Distinct"); same value as <see cref="UniqueKmers"/>.</summary>
    public int DistinctKmers { get; init; }

    /// <summary>Number of k-mers occurring exactly once (Jellyfish "Unique").</summary>
    public int SingletonKmers { get; init; }
}

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

/// <summary>A MISA compound microsatellite (types c / c*).</summary>
public record CompoundMicrosatelliteItem(
    int Start,
    int End,
    int Length,
    string Type,
    string Notation,
    MicrosatelliteItem[] Components);

/// <summary>Result of <c>find_microsatellites</c>.</summary>
public record FindMicrosatellitesResult(MicrosatelliteItem[] Items)
{
    /// <summary>MISA compound microsatellites (only when <c>maxCompoundInterruption</c> ≥ 0; otherwise null).</summary>
    public CompoundMicrosatelliteItem[]? Compounds { get; init; }
}

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
    int TotalLength)
{
    /// <summary>Mismatched pairs inside the stem (0 unless <c>maxMismatches</c> &gt; 0).</summary>
    public int Mismatches { get; init; }
}

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
    string? MostFrequentUnit)
{
    /// <summary>Count of pentanucleotide (5 bp unit) STRs.</summary>
    public int PentanucleotideRepeats { get; init; }

    /// <summary>Count of hexanucleotide (6 bp unit) STRs.</summary>
    public int HexanucleotideRepeats { get; init; }

    /// <summary>STR count per searched unit size in bp (1-6, or every size of misaDefinition, incl. sizes above 6); sums to TotalRepeats (misa.pl .statistics "Distribution to different repeat type classes").</summary>
    public Dictionary<int, int>? CountsByUnitLength { get; init; }

    /// <summary>STR counts per MISA repeat-type class (rotations + reverse complement, e.g. "AC/GT").</summary>
    public Dictionary<string, int>? CanonicalMotifCounts { get; init; }

    /// <summary>STR counts per Krait standard motif (only when <c>standardMotifLevel</c> is 0-4; otherwise null).</summary>
    public Dictionary<string, int>? StandardMotifCounts { get; init; }
}

/// <summary>An EMBOSS einverted scored inverted repeat (0-based inclusive arm coordinates).</summary>
public record ScoredInvertedRepeatItem(
    int LeftArmStart,
    int LeftArmEnd,
    int RightArmStart,
    int RightArmEnd,
    int Score,
    int Matches,
    int Mismatches,
    int Gaps,
    string LeftArmAlignment,
    string MatchLine,
    string RightArmAlignment,
    int LeftArmLength,
    int RightArmLength,
    int LoopLength,
    double PercentMatches);

/// <summary>Result of <c>find_inverted_repeats_scored</c>.</summary>
public record FindInvertedRepeatsScoredResult(ScoredInvertedRepeatItem[] Items);

/// <summary>A maximal reverse-complement repeat pair.</summary>
public record ReverseComplementRepeatItem(
    int FirstPosition,
    int SecondPosition,
    string RepeatSequence,
    string SecondSequence,
    int Length,
    int Spacing);

/// <summary>Result of <c>find_reverse_complement_repeats</c>.</summary>
public record FindReverseComplementRepeatsResult(ReverseComplementRepeatItem[] Items);

/// <summary>A maximal k-mismatch direct repeat.</summary>
public record ApproximateDirectRepeatItem(
    int FirstPosition,
    int SecondPosition,
    int Length,
    int Mismatches,
    int Spacing,
    string FirstCopy,
    string SecondCopy);

/// <summary>Result of <c>find_approximate_direct_repeats</c>.</summary>
public record FindApproximateDirectRepeatsResult(ApproximateDirectRepeatItem[] Items);

/// <summary>A maximal degenerate (k-differences / k-mismatches) repeat.</summary>
public record DegenerateRepeatItem(
    int FirstPosition,
    int FirstLength,
    int SecondPosition,
    int SecondLength,
    int Distance,
    int Spacing,
    string FirstCopy,
    string SecondCopy,
    bool IsReverseComplement);

/// <summary>Result of <c>find_degenerate_repeats</c>.</summary>
public record FindDegenerateRepeatsResult(DegenerateRepeatItem[] Items);

/// <summary>A supermaximal repeat with all its occurrences.</summary>
public record SupermaximalRepeatItem(string Sequence, int Length, int[] Positions);

/// <summary>Result of <c>find_supermaximal_repeats</c>.</summary>
public record FindSupermaximalRepeatsResult(SupermaximalRepeatItem[] Items);

/// <summary>A Tandem Repeats Finder row (0-based start; TRF .dat columns).</summary>
public record ApproximateTandemRepeatItem(
    int Start,
    int SpanLength,
    int Period,
    int ConsensusSize,
    string Consensus,
    double CopyNumber,
    double PercentMatches,
    double PercentIndels,
    int AlignmentScore)
{
    /// <summary>Percentage of A in the repeat.</summary>
    public double PercentA { get; init; }

    /// <summary>Percentage of C in the repeat.</summary>
    public double PercentC { get; init; }

    /// <summary>Percentage of G in the repeat.</summary>
    public double PercentG { get; init; }

    /// <summary>Percentage of T in the repeat.</summary>
    public double PercentT { get; init; }

    /// <summary>Shannon entropy (bits) over A/C/G/T only.</summary>
    public double Entropy { get; init; }

    /// <summary>TRF's entropy column (denominator = all non-gap symbols).</summary>
    public double EntropyTrf { get; init; }

    /// <summary>Final alignment row of the repeat.</summary>
    public string? AlignedSequence { get; init; }

    /// <summary>Final alignment row of the consensus copies.</summary>
    public string? AlignedConsensus { get; init; }

    /// <summary>Left flank (only when flankLength &gt; 0).</summary>
    public string? LeftFlank { get; init; }

    /// <summary>Right flank (only when flankLength &gt; 0).</summary>
    public string? RightFlank { get; init; }

    /// <summary>Matching column pairs between adjacent copies (numerator of PercentMatches).</summary>
    public int CopyMatches { get; init; }

    /// <summary>Mismatching column pairs between adjacent copies.</summary>
    public int CopyMismatches { get; init; }

    /// <summary>Indel columns between adjacent copies (numerator of PercentIndels).</summary>
    public int CopyIndels { get; init; }

    /// <summary>1-based report rank before MaxPeriod filtering / redundancy elimination (TRF OUTPUTcount; TRF row order).</summary>
    public int OutputIndex { get; init; }

    /// <summary>Alignments reported for the whole sequence before MaxPeriod filtering / redundancy elimination (TRF final OUTPUTcount).</summary>
    public int OutputCount { get; init; }

    /// <summary>0-based position at which the candidate was detected (TRF alignment page "Found at i:" = this + 1).</summary>
    public int DetectionPosition { get; init; }

    /// <summary>Candidate distance at detection (TRF "original size"; 0 = not set).</summary>
    public int DetectionDistance { get; init; }
}

/// <summary>One TRF HTML page: TRF's file name and the page content.</summary>
public record TrfHtmlPageItem(string FileName, string Html);

/// <summary>Result of <c>find_approximate_tandem_repeats</c>.</summary>
public record FindApproximateTandemRepeatsResult(ApproximateTandemRepeatItem[] Items)
{
    /// <summary>TRF .dat (format 'dat') or -ngs (format 'ngs') text; null for 'json' / 'html'.</summary>
    public string? Formatted { get; init; }

    /// <summary>TRF repeat-table HTML pages (format 'html'); null otherwise.</summary>
    public TrfHtmlPageItem[]? HtmlPages { get; init; }

    /// <summary>TRF alignment pages (<c>.N.txt.html</c>, the targets of the table links; format 'html'); null otherwise.</summary>
    public TrfHtmlPageItem[]? AlignmentPages { get; init; }
}

/// <summary>Result of <c>mask_approximate_tandem_repeats</c>.</summary>
public record MaskApproximateTandemRepeatsResult(string Masked);

/// <summary>Result of <c>tandem_repeat_bernoulli_statistics</c>.</summary>
public record TandemRepeatBernoulliStatisticsResult(
    int Period,
    int AdjacentCopyPairs,
    int BernoulliTrials,
    int Matches,
    int Mismatches,
    int Indels,
    double MatchProbability,
    double IndelProbability,
    double PercentMatches,
    double PercentIndels,
    double ExpectedMatches,
    bool MeetsExpectedMatchProbability);

/// <summary>Result of <c>standardize_repeat_motif</c>.</summary>
public record StandardizeRepeatMotifResult(string CanonicalClass, string StandardMotif, int Level);

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

/// <summary>A both-strand PWM hit (Biopython <c>search(both=True)</c> coordinates in <c>BiopythonPosition</c>).</summary>
public record PwmStrandMatchItem(int Position, int BiopythonPosition, string Strand, string MatchedSequence, string Pattern, double Score);

/// <summary>Result of <c>scan_with_pwm_both_strands</c>.</summary>
public record ScanWithPwmBothStrandsResult(PwmStrandMatchItem[] Items);

/// <summary>Result of <c>pwm_score_thresholds</c>: score-distribution grid and derived thresholds.</summary>
public record PwmScoreThresholdsResult(
    double MinScore,
    double Step,
    int PointCount,
    double MeanScore,
    double ThresholdFpr,
    double ThresholdFnr,
    double ThresholdBalanced,
    double BalancedFalsePositiveRate,
    double ThresholdPatser);

/// <summary>Result of <c>pwm_score_pvalue</c>: score threshold and its exact p-value P(S &gt;= Score).</summary>
/// <param name="Score">The requested score, or the computed threshold; null when no word score qualifies (ScoreAboveMaximum).</param>
/// <param name="ScoreAboveMaximum">True when even the best word has P(S &gt;= max) &gt; pValue (threshold above the maximum, p-value 0).</param>
/// <param name="PValue">P(S &gt;= Score); the conservative upper bound when IsExact is false.</param>
/// <param name="PValueLowerBound">Certified lower bound.</param>
/// <param name="PValueUpperBound">Certified upper bound.</param>
/// <param name="IsExact">True when the p-value is exact.</param>
/// <param name="Granularity">Scale g of the integer-rounded matrix that resolved it (0 = decided without rounding).</param>
public record PwmScorePValueResult(
    double? Score,
    bool ScoreAboveMaximum,
    double PValue,
    double PValueLowerBound,
    double PValueUpperBound,
    bool IsExact,
    double Granularity);

/// <summary>A promoter-element weight-matrix hit.</summary>
public record PromoterMatrixHitItem(string Name, string MatrixId, int Position, string Strand, string Sequence, double Score, double Threshold);

/// <summary>Result of <c>find_promoter_elements_by_matrix</c>.</summary>
public record FindPromoterElementsByMatrixResult(PromoterMatrixHitItem[] Items);

/// <summary>Result of <c>create_alphabet_pwm</c>: K×L log-odds (rows in alphabet order; null = −∞) and summary statistics.</summary>
public record AlphabetPwmResult(
    string Alphabet,
    double?[][] Matrix,
    int Length,
    string Consensus,
    string Anticonsensus,
    double? MaxScore,
    double? MinScore,
    double Mean,
    double Std);

/// <summary>Result of <c>scan_with_alphabet_pwm</c>: every window score (null when not finite), NaN windows, and hits.</summary>
public record ScanWithAlphabetPwmResult(
    string Consensus,
    int Length,
    double?[] Scores,
    int[] InvalidWindows,
    MotifMatchItem[] Hits);

/// <summary>A paired σ70 −35/−10 consensus candidate.</summary>
public record Sigma70CandidateItem(
    string Strand,
    int Minus35Start,
    string Minus35,
    int Minus10Start,
    string Minus10,
    int Spacer,
    int Mismatches35,
    int Mismatches10,
    int TotalMismatches,
    int SpacerDeviation);

/// <summary>Result of <c>find_sigma70_promoters</c>.</summary>
public record FindSigma70PromotersResult(Sigma70CandidateItem[] Items);

/// <summary>Promoter Calculator v1.0 minimum-ΔG configuration of one TSS.</summary>
public record Sigma70PredictionItem(
    string Strand,
    int Tss,
    string PromoterSequence,
    string Up,
    string Minus35,
    string Spacer,
    string Minus10,
    string Discriminator,
    string Itr,
    int UpStart,
    int Minus35Start,
    int SpacerStart,
    int Minus10Start,
    int DiscriminatorStart,
    double DeltaGTotal,
    double DeltaG10,
    double DeltaG35,
    double DeltaGDiscriminator,
    double DeltaGItr,
    double DeltaGExtended10,
    double DeltaGSpacer,
    double DeltaGUp,
    double DeltaGBind,
    double TranscriptionRate);

/// <summary>Result of <c>predict_sigma70_promoters</c>: per-TSS predictions and the overall lowest-ΔG one.</summary>
public record PredictSigma70PromotersResult(Sigma70PredictionItem[] Items, Sigma70PredictionItem? Best);

/// <summary>A strand-annotated regulatory element occurrence.</summary>
public record StrandedRegulatoryElementItem(string Name, int Position, string Sequence, string Pattern, string Description, string Strand);

/// <summary>Result of <c>find_regulatory_elements_both_strands</c>.</summary>
public record FindRegulatoryElementsBothStrandsResult(StrandedRegulatoryElementItem[] Items);

/// <summary>A k-mer (or reverse-complement pair) scored by RSAT oligo-analysis occurrence statistics.</summary>
public record OligoMotifItem(
    string Sequence,
    string? ReverseComplement,
    int Count,
    int[] Positions,
    double ExpectedFrequency,
    double ExpectedOccurrences,
    double? Ratio,
    double OccurrenceProbability,
    double OccurrenceEValue,
    double OccurrenceSignificance)
{
    /// <summary>Sequence index of each position (0 = <c>sequence</c>); RSAT-option runs only.</summary>
    public int[]? SequenceIndices { get; init; }

    /// <summary>Discarded overlapping windows (<c>-noov</c>, RSAT <c>ovl_occ</c>); RSAT-option runs only.</summary>
    public int? Overlaps { get; init; }

    /// <summary>occ / total occurrences (RSAT <c>obs_freq</c>); RSAT-option runs only.</summary>
    public double? ObservedFrequency { get; init; }

    /// <summary>RSAT <c>exp_var</c> (with zscore or a calibration).</summary>
    public double? ExpectedVariance { get; init; }

    /// <summary>RSAT <c>ovlp</c> overlap coefficient (with zscore).</summary>
    public double? OverlapCoefficient { get; init; }

    /// <summary>RSAT z-score (occ − exp_occ)/√exp_var (with zscore; null when the variance is ≤ 0).</summary>
    public double? ZScore { get; init; }

    /// <summary>Distribution behind occ_P: Binomial, Poisson or NegativeBinomial; RSAT-option runs only.</summary>
    public string? FittedDistribution { get; init; }

    /// <summary>Best lexicon segmentation "prefix|suffix" (background 'lexicon').</summary>
    public string? LexiconSegmentation { get; init; }
}

/// <summary>Result of <c>oligo_analysis</c>.</summary>
public record OligoAnalysisResultDto(
    OligoMotifItem[] Motifs,
    int OligoLength,
    string Strands,
    bool CountOverlapping,
    long TotalOccurrences,
    int TestedPatterns,
    double? PossibleOligos)
{
    /// <summary>Number of analysed sequences; RSAT-option runs only.</summary>
    public int? SequenceCount { get; init; }

    /// <summary>Degenerate-word mode ('none', 'oneN', 'onedeg'); RSAT-option runs only.</summary>
    public string? Degenerate { get; init; }

    /// <summary>Sequence type ('protein' or 'other', RSAT -seqtype); null for DNA.</summary>
    public string? SequenceType { get; init; }

    /// <summary>RSAT alphabet_size (20 for protein, number of distinct residues for other); null for DNA.</summary>
    public int? AlphabetSize { get; init; }
}

/// <summary>A spaced dyad (or reverse-complement pair) scored by RSAT dyad-analysis.</summary>
public record DyadItem(
    string Pattern,
    string FirstMonad,
    int Spacing,
    string SecondMonad,
    string? ReverseComplement,
    int Occurrences,
    int Overlaps,
    int[] SequenceIndices,
    int[] Positions,
    double ObservedFrequency,
    double? ExpectedFrequency,
    double? ExpectedOccurrences,
    double? ExpectedVariance,
    double OverlapCoefficient,
    double? ZScore,
    double Ratio,
    double? OccurrenceProbability,
    double? OccurrenceEValue,
    double? OccurrenceSignificance,
    bool IsDirectRepeat,
    bool IsReversePalindrome,
    bool ExpectedFromMonads);

/// <summary>Per-spacing positions and counts of <c>dyad_analysis</c>.</summary>
public record DyadSpacingItem(int Spacing, long PossiblePositions, long Occurrences, long Overlaps);

/// <summary>Result of <c>dyad_analysis</c>.</summary>
public record DyadAnalysisResultDto(
    DyadItem[] Dyads,
    int MonadLength,
    int MinSpacing,
    int MaxSpacing,
    string DyadType,
    string Strands,
    bool CountOverlapping,
    int SequenceCount,
    long MonadOccurrences,
    DyadSpacingItem[] Spacings,
    int TestedPatterns,
    double? PossibleDyads);

/// <summary>A k-mer (or pair) scored by RSAT oligo-analysis matching-sequence statistics.</summary>
public record SharedMotifSignificanceItem(
    string Sequence,
    string? ReverseComplement,
    int[] SequenceIndices,
    double Prevalence,
    double ExpectedFrequency,
    double ExpectedMatchingSequences,
    double MatchingSequenceProbability,
    double? MatchingSequenceEValue,
    double? MatchingSequenceSignificance);

/// <summary>Result of <c>shared_motifs_significance</c>.</summary>
public record SharedMotifSignificanceResult(
    SharedMotifSignificanceItem[] Motifs,
    int OligoLength,
    string Strands,
    int SequenceCount,
    long PossiblePositions,
    double? PossibleOligos)
{
    /// <summary>Degenerate-word mode ('oneN' or 'onedeg', RSAT -oneN / -onedeg); null for plain words.</summary>
    public string? Degenerate { get; init; }
}

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

/// <summary>A 0-based half-open [Start, End) low-complexity interval.</summary>
public record ComplexityIntervalItem(int Start, int End, int Length);

/// <summary>Result of <c>find_low_complexity_intervals</c>.</summary>
public record FindLowComplexityIntervalsResult(ComplexityIntervalItem[] Items);

/// <summary>Result of <c>longdust_score</c>: S_L(x) and the number of k-mer positions l(x).</summary>
public record LongdustScoreResult(double Score, int KmerPositions);

/// <summary>Result of <c>find_longdust_regions</c>.</summary>
public record FindLongdustRegionsResult(ComplexityIntervalItem[] Items);

/// <summary>Result of <c>lempel_ziv_complexity</c>.</summary>
public record LempelZivComplexityResult(int Complexity, double Normalized);

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
