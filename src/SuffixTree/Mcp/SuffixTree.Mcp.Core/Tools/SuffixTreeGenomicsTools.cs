using System.ComponentModel;
using ModelContextProtocol.Server;

namespace SuffixTree.Mcp.Core.Tools;

/// <summary>
/// MCP tools for genomic analysis operations backed by suffix tree algorithms.
/// </summary>
[McpServerToolType]
public class SuffixTreeGenomicsTools
{
    /// <summary>
    /// Find the longest repeated region in a DNA sequence.
    /// </summary>
    [McpServerTool(Name = "find_longest_repeat", Title = "Genomics — Longest Repeated Region", ReadOnly = true)]
    [Description("Find the longest repeated region in a DNA sequence. Returns the repeat sequence and all positions.")]
    public static FindLongestRepeatResult FindLongestRepeat(
        [Description("The DNA sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));

        var dna = new global::Seqeron.Genomics.Core.DnaSequence(sequence);
        var result = GenomicAnalyzer.FindLongestRepeat(dna);

        return new FindLongestRepeatResult(result.Sequence, result.Positions.ToArray(), result.Length);
    }

    /// <summary>
    /// Find the longest common region between two DNA sequences.
    /// </summary>
    [McpServerTool(Name = "find_longest_common_region", Title = "Genomics — Longest Common Region", ReadOnly = true)]
    [Description("Find the longest common region between two DNA sequences.")]
    public static FindLongestCommonRegionResult FindLongestCommonRegion(
        [Description("The first DNA sequence")] string sequence1,
        [Description("The second DNA sequence")] string sequence2)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));

        var dna1 = new global::Seqeron.Genomics.Core.DnaSequence(sequence1);
        var dna2 = new global::Seqeron.Genomics.Core.DnaSequence(sequence2);
        var result = GenomicAnalyzer.FindLongestCommonRegion(dna1, dna2);

        return new FindLongestCommonRegionResult(result.Sequence, result.PositionInFirst, result.PositionInSecond, result.Length);
    }

    /// <summary>
    /// Calculate similarity between two DNA sequences using k-mer Jaccard index.
    /// </summary>
    [McpServerTool(Name = "calculate_similarity", Title = "Genomics — Sequence Similarity", ReadOnly = true)]
    [Description("Calculate similarity between two DNA sequences using k-mer Jaccard index (0-100 percentage scale).")]
    public static CalculateSimilarityResult CalculateSimilarity(
        [Description("The first DNA sequence")] string sequence1,
        [Description("The second DNA sequence")] string sequence2,
        [Description("K-mer size (default: 5)")] int kmerSize = 5)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));

        var dna1 = new global::Seqeron.Genomics.Core.DnaSequence(sequence1);
        var dna2 = new global::Seqeron.Genomics.Core.DnaSequence(sequence2);
        var similarity = GenomicAnalyzer.CalculateSimilarity(dna1, dna2, kmerSize);

        return new CalculateSimilarityResult(similarity);
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
    [Description("Calculate edit distance (Levenshtein distance) between two sequences. Returns minimum number of edits needed.")]
    public static EditDistanceResult EditDistance(
        [Description("The first sequence")] string sequence1,
        [Description("The second sequence")] string sequence2)
    {
        if (string.IsNullOrEmpty(sequence1))
            throw new ArgumentException("Sequence1 cannot be null or empty", nameof(sequence1));
        if (string.IsNullOrEmpty(sequence2))
            throw new ArgumentException("Sequence2 cannot be null or empty", nameof(sequence2));

        var distance = ApproximateMatcher.EditDistance(sequence1, sequence2);
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
