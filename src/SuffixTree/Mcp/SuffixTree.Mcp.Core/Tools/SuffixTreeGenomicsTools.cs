using System.ComponentModel;
using ModelContextProtocol.Server;

namespace SuffixTree.Mcp.Core.Tools;

/// <summary>
/// MCP tools for genomic analysis operations backed by suffix tree algorithms.
/// </summary>
[McpServerToolType]
public class SuffixTreeGenomicsTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private SuffixTreeGenomicsTools() { }

    /// <summary>
    /// Find the longest repeated region in a DNA sequence.
    /// </summary>
    [McpServerTool(Name = "find_longest_repeat", Title = "Genomics — Longest Repeated Region", ReadOnly = true)]
    [Description("Find the longest repeated region in a nucleotide sequence (overlapping occurrences allowed). Returns the repeat sequence and all positions; on a length tie the substring whose first occurrence is leftmost. Case-insensitive; only A/C/G/T match, N / IUPAC codes / other symbols never repeat and end a repeat. bothStrands=true: the longest repeat over both strands (a string and its reverse complement are one class; reversePositions lists the windows that spell the reverse complement).")]
    public static FindLongestRepeatResult FindLongestRepeat(
        [Description("The nucleotide sequence to analyze (A/C/G/T plus IUPAC codes such as N, case-insensitive)")] string sequence,
        [Description("Search both strands (default false = forward strand only)")] bool bothStrands = false,
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        cancellationToken.ThrowIfCancellationRequested();
        if (bothStrands)
        {
            var all = GenomicAnalyzer.FindAllLongestRepeatsBothStrands(sequence, cancellationToken);
            var best = all.Count == 0 ? StrandedRepeatInfo.None : all[0];
            return new FindLongestRepeatResult(best.Sequence, best.Positions.ToArray(), best.Length, best.ReversePositions.ToArray());
        }

        var result = GenomicAnalyzer.FindLongestRepeat(sequence);
        cancellationToken.ThrowIfCancellationRequested();
        return new FindLongestRepeatResult(result.Sequence, result.Positions.ToArray(), result.Length);
    }

    /// <summary>
    /// Find the longest common region between two DNA sequences.
    /// </summary>
    [McpServerTool(Name = "find_longest_common_region", Title = "Genomics — Longest Common Region", ReadOnly = true)]
    [Description("Find the longest common region (longest common substring) between two nucleotide sequences. Case-insensitive; only A/C/G/T match, N / IUPAC codes / other symbols match nothing and end a region. strands: Forward (default) | ReverseComplement | Both (the region may be shared with the reverse complement of the second sequence; position2 is then the forward-strand start of the window that spells it and isReverseComplement is true). maxMismatches > 0: the longest pair of equal-length windows that differ in at most that many positions (Hamming, no gaps; Flouri et al. 2015), returned with region2 and the mismatch count. Cancellation is checked before and after the search (the search itself is not interruptible).")]
    public static FindLongestCommonRegionResult FindLongestCommonRegion(
        [Description("The first nucleotide sequence (A/C/G/T plus IUPAC codes such as N, case-insensitive)")] string sequence1,
        [Description("The second nucleotide sequence")] string sequence2,
        [Description("Which strands of the second sequence (case-insensitive): Forward (default) | ReverseComplement | Both")] string strands = "Forward",
        [Description("Largest number of mismatching positions of the two windows (default 0 = exact common substring)")] int maxMismatches = 0,
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));
        if (maxMismatches < 0)
            throw new ArgumentOutOfRangeException(nameof(maxMismatches), "maxMismatches must be >= 0");
        if (!Enum.TryParse<StrandSelection>(strands, ignoreCase: true, out var selection) || !Enum.IsDefined(selection) || selection == 0 || int.TryParse(strands, out _))
            throw new ArgumentException($"strands must be Forward, ReverseComplement or Both, got '{strands}'", nameof(strands));

        cancellationToken.ThrowIfCancellationRequested();
        FindLongestCommonRegionResult result;
        if (maxMismatches > 0)
        {
            var a = GenomicAnalyzer.FindLongestCommonRegionWithMismatches(sequence1, sequence2, maxMismatches, selection);
            result = new FindLongestCommonRegionResult(a.Sequence1, a.PositionInFirst, a.PositionInSecond, a.Length, a.IsReverseComplement, a.Mismatches, a.Sequence2);
        }
        else
        {
            var c = GenomicAnalyzer.FindLongestCommonRegion(sequence1, sequence2, selection);
            result = new FindLongestCommonRegionResult(c.Sequence, c.PositionInFirst, c.PositionInSecond, c.Length, c.IsReverseComplement);
        }

        cancellationToken.ThrowIfCancellationRequested();
        return result;
    }

    /// <summary>
    /// Calculate similarity between two DNA sequences using k-mer Jaccard index.
    /// </summary>
    [McpServerTool(Name = "calculate_similarity", Title = "Genomics — Sequence Similarity", ReadOnly = true)]
    [Description("Calculate similarity between two DNA sequences using the k-mer Jaccard index |A∩B|/|A∪B| of the distinct k-mer sets, as a percentage (0-100 scale; 0 when both sets are empty). Default: the exact index; compares literal (strand-specific) k-mers like a Python set / scipy jaccard; canonical=true uses min(k-mer, reverse complement) as Mash and sourmash do, so a sequence and its reverse complement score 100. sketchSize = s estimates it from a Mash bottom-s MinHash sketch (MurmurHash3-x64-128 seed 42, as 'mash dist'; k <= 32; canonical k-mers are Mash's default); scaled = S from a sourmash FracMinHash sketch (also gives the containment indices and containment ANI). The result always adds jaccardDistance (1 - J), mashDistance (-ln(2J/(1+J))/k) and jaccardAni = (2J/(1+J))^(1/k); the exact and FracMinHash estimators add containmentAInB / containmentBInA / maxContainment and containmentAniAInB / containmentAniBInA (C^(1/k), sourmash containment ANI point estimate); sketches add sharedHashes / denominatorHashes (and mashPValue for a bottom-s sketch). With scaled the result also carries aniEstimates: the sourmash ANI results containment_ani (both ways), max_containment_ani and jaccard_ani (sourmash.distance_utils; n_unique_kmers = sketch hashes x scaled) with distance, the probability that two such sketches share nothing by chance (probabilityNothingInCommon / probabilityExceedsThreshold, 1e-3), sizeIsInaccurate (a sketch size estimate off by more than 20 % with probability above 5 % withholds the ANI), the Jaccard-ANI error bound (jaccardError, withheld above 1e-4) and, with estimateAniInterval = true, the confidence interval of the containment ANI at the given confidence (distanceLow / distanceHigh, aniLow / aniHigh; Hera et al. 2023, Brent root search; null when the sketch is too small). Plain A/C/G/T input only (sequences with N: the Analysis server accepts IUPAC input). Verified against sourmash 4.9.4 (1 120 ANI / interval rows) and mash 2.3 (400 distance / p-value rows).")]
    public static CalculateSimilarityResult CalculateSimilarity(
        [Description("The first nucleotide sequence (N runs, IUPAC codes and lower case accepted: a k-mer window with a symbol other than A/C/G/T is skipped)")] string sequence1,
        [Description("The second nucleotide sequence")] string sequence2,
        [Description("K-mer size (default: 5)")] int kmerSize = 5,
        [Description("Use canonical (strand-agnostic) k-mers min(w, reverse complement) as Mash/sourmash (default: false)")] bool canonical = false,
        [Description("Estimate the index from a Mash bottom-s MinHash sketch of this size (>= 1; k <= 32). Cannot be combined with scaled.")] int? sketchSize = null,
        [Description("Estimate from a sourmash FracMinHash sketch with this scaled factor (1..4294967295). Cannot be combined with sketchSize.")] long? scaled = null,
        [Description("With scaled: compute the confidence interval of the containment ANI (sourmash estimate_ci; default false).")] bool estimateAniInterval = false,
        [Description("Confidence level of the ANI interval in (0, 1) (default 0.95).")] double confidence = 0.95,
        CancellationToken cancellationToken = default)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));
        cancellationToken.ThrowIfCancellationRequested();
        // strings: N runs, IUPAC codes and lower case are accepted; a k-mer window with a symbol other than A/C/G/T is skipped (Mash -n / Jellyfish rule)
        var e = GenomicAnalyzer.CompareSimilarity(sequence1, sequence2, kmerSize, canonical, sketchSize, scaled, new AniOptions(estimateAniInterval, confidence));
        cancellationToken.ThrowIfCancellationRequested();
        static AniEstimateDto Dto(AniEstimate a)
        {
            bool? jaccardExceeds = null;
            if (a.JaccardError is not null)
                jaccardExceeds = a.JaccardErrorExceedsThreshold;
            return new AniEstimateDto(a.Distance, a.DistanceLow, a.DistanceHigh, a.Ani, a.AniLow, a.AniHigh, a.ProbabilityNothingInCommon, a.ProbabilityExceedsThreshold, a.SizeIsInaccurate, a.JaccardError, jaccardExceeds);
        }

        var ani = e.AniEstimates is { } all
            ? new AniEstimatesDto(Dto(all.ContainmentAInB), Dto(all.ContainmentBInA), Dto(all.MaxContainment), Dto(all.Jaccard))
            : null;
        return new CalculateSimilarityResult(
            e.Jaccard * 100.0, e.Estimator.ToString(), e.JaccardDistance, e.MashDistance, e.JaccardAni, e.ContainmentAInB, e.ContainmentBInA, e.MaxContainment,
            e.ContainmentAniAInB, e.ContainmentAniBInA, e.SharedHashes, e.DenominatorHashes, e.MashPValue, ani);
    }

    /// <summary>
    /// Calculate Hamming distance between two sequences of equal length.
    /// </summary>
    [McpServerTool(Name = "hamming_distance", Title = "Genomics — Hamming Distance", ReadOnly = true)]
    [Description("Calculate Hamming distance between two sequences of equal length. Returns number of positions with different characters.")]
    public static HammingDistanceResult HammingDistance(
        [Description("The first sequence")] string sequence1,
        [Description("The second sequence")] string sequence2)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));
        if (sequence1.Length != sequence2.Length)
            throw new ArgumentException("Sequences must have equal length for Hamming distance");

        var distance = ApproximateMatcher.HammingDistance(sequence1, sequence2);
        return new HammingDistanceResult(distance);
    }

    /// <summary>
    /// Calculate edit distance (Levenshtein distance) between two sequences.
    /// </summary>
    [McpServerTool(Name = "edit_distance", Title = "Genomics — Edit Distance", ReadOnly = true)]
    [Description("Calculate edit distance (Levenshtein distance) between two sequences. Returns minimum number of edits needed. Optional insertionCost/deletionCost/substitutionCost (default 1) give the weighted Levenshtein distance with rapidfuzz Levenshtein.distance(s1, s2, weights=(insertion, deletion, substitution)) semantics: transforming sequence1 into sequence2, an insertion adds a character of sequence2, a deletion removes one of sequence1.")]
    public static EditDistanceResult EditDistance(
        [Description("The first sequence")] string sequence1,
        [Description("The second sequence")] string sequence2,
        [Description("Cost of inserting a character of sequence2 (>= 0; default 1)")] int insertionCost = 1,
        [Description("Cost of deleting a character of sequence1 (>= 0; default 1)")] int deletionCost = 1,
        [Description("Cost of substituting a character (>= 0; default 1)")] int substitutionCost = 1)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));
        if (insertionCost < 0)
            throw new ArgumentOutOfRangeException(nameof(insertionCost), "insertionCost must be >= 0");
        if (deletionCost < 0)
            throw new ArgumentOutOfRangeException(nameof(deletionCost), "deletionCost must be >= 0");
        if (substitutionCost < 0)
            throw new ArgumentOutOfRangeException(nameof(substitutionCost), "substitutionCost must be >= 0");

        var distance = insertionCost == 1 && deletionCost == 1 && substitutionCost == 1
            ? ApproximateMatcher.EditDistance(sequence1, sequence2)
            : ApproximateMatcher.EditDistance(sequence1, sequence2, insertionCost, deletionCost, substitutionCost);
        return new EditDistanceResult(distance);
    }

    /// <summary>
    /// Count approximate occurrences of a pattern in a sequence.
    /// </summary>
    [McpServerTool(Name = "count_approximate_occurrences", Title = "Genomics — Approximate Pattern Count", ReadOnly = true)]
    [Description("Count approximate occurrences of a pattern in a sequence, allowing up to maxMismatches substitutions.")]
    public static CountApproximateOccurrencesResult CountApproximateOccurrences(
        [Description("The sequence to search in")] string sequence,
        [Description("The pattern to find")] string pattern,
        [Description("Maximum number of allowed mismatches")] int maxMismatches)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException("Pattern cannot be null or empty", nameof(pattern));
        if (maxMismatches < 0)
            throw new ArgumentException("MaxMismatches cannot be negative", nameof(maxMismatches));

        var count = ApproximateMatcher.CountApproximateOccurrences(sequence, pattern, maxMismatches);
        return new CountApproximateOccurrencesResult(count);
    }
}
