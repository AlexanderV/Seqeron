using System.ComponentModel;
using ModelContextProtocol.Server;

namespace SuffixTree.Mcp.Core.Tools;

/// <summary>
/// MCP tools for core suffix tree operations.
/// </summary>
[McpServerToolType]
public class SuffixTreeTools
{
    // Compatibility wrappers preserved for direct static calls from tests/consumers.
    // MCP tool registration for suffix-tree core operations is hosted in SuffixTreeCoreTools.

    /// <summary>
    /// Check if a pattern exists in text using suffix tree.
    /// </summary>
    /// <param name="text">The text to search in.</param>
    /// <param name="pattern">The pattern to search for.</param>
    /// <returns>Result indicating whether pattern was found.</returns>
    public static SuffixTreeContainsResult SuffixTreeContains(
        string text,
        string pattern) => SuffixTreeCoreTools.SuffixTreeContains(text, pattern);

    /// <summary>
    /// Count occurrences of a pattern in text using suffix tree.
    /// </summary>
    public static SuffixTreeCountResult SuffixTreeCount(
        string text,
        string pattern) => SuffixTreeCoreTools.SuffixTreeCount(text, pattern);

    /// <summary>
    /// Find all positions where a pattern occurs in text.
    /// </summary>
    public static SuffixTreeFindAllResult SuffixTreeFindAll(
        string text,
        string pattern) => SuffixTreeCoreTools.SuffixTreeFindAll(text, pattern);

    /// <summary>
    /// Find the longest repeated substring in text.
    /// </summary>
    public static SuffixTreeLrsResult SuffixTreeLrs(
        string text) => SuffixTreeCoreTools.SuffixTreeLrs(text);

    /// <summary>
    /// Find the longest common substring between two texts.
    /// </summary>
    public static SuffixTreeLcsResult SuffixTreeLcs(
        string text1,
        string text2) => SuffixTreeCoreTools.SuffixTreeLcs(text1, text2);

    /// <summary>
    /// Get statistics about a suffix tree built from text.
    /// </summary>
    public static SuffixTreeStatsResult SuffixTreeStats(
        string text) => SuffixTreeCoreTools.SuffixTreeStats(text);

    /// <summary>
    /// Find the longest repeated region in a DNA sequence.
    /// </summary>
    [McpServerTool(Name = "find_longest_repeat")]
    [Description("Find the longest repeated region in a DNA sequence. Returns the repeat sequence and all positions.")]
    public static FindLongestRepeatResult FindLongestRepeat(
        [Description("The DNA sequence to analyze")] string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            throw new ArgumentException("Sequence cannot be null or empty", nameof(sequence));
        }

        var dna = new global::Seqeron.Genomics.Core.DnaSequence(sequence);
        var result = GenomicAnalyzer.FindLongestRepeat(dna);

        return new FindLongestRepeatResult(result.Sequence, result.Positions.ToArray(), result.Length);
    }

    /// <summary>
    /// Find the longest common region between two DNA sequences.
    /// </summary>
    [McpServerTool(Name = "find_longest_common_region")]
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
    [McpServerTool(Name = "calculate_similarity")]
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
    [McpServerTool(Name = "hamming_distance")]
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
    [McpServerTool(Name = "edit_distance")]
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
    [McpServerTool(Name = "count_approximate_occurrences")]
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

/// <summary>
/// Result of suffix_tree_contains operation.
/// </summary>
/// <param name="Found">Whether the pattern was found in the text.</param>
public record SuffixTreeContainsResult(bool Found);

/// <summary>
/// Result of suffix_tree_count operation.
/// </summary>
/// <param name="Count">Number of pattern occurrences in the text.</param>
public record SuffixTreeCountResult(int Count);

/// <summary>
/// Result of suffix_tree_find_all operation.
/// </summary>
/// <param name="Positions">Array of positions where pattern was found.</param>
public record SuffixTreeFindAllResult(int[] Positions);

/// <summary>
/// Result of suffix_tree_lrs operation.
/// </summary>
public record SuffixTreeLrsResult(string Substring, int Length);

/// <summary>
/// Result of suffix_tree_lcs operation.
/// </summary>
public record SuffixTreeLcsResult(string Substring, int Length);

/// <summary>
/// Result of suffix_tree_stats operation.
/// </summary>
public record SuffixTreeStatsResult(int NodeCount, int LeafCount, int MaxDepth, int TextLength);

/// <summary>
/// Result of find_longest_repeat operation.
/// </summary>
public record FindLongestRepeatResult(string Repeat, int[] Positions, int Length);

/// <summary>
/// Result of find_longest_common_region operation.
/// </summary>
public record FindLongestCommonRegionResult(string Region, int Position1, int Position2, int Length);

/// <summary>
/// Result of calculate_similarity operation.
/// </summary>
public record CalculateSimilarityResult(double Similarity);

/// <summary>
/// Result of hamming_distance operation.
/// </summary>
public record HammingDistanceResult(int Distance);

/// <summary>
/// Result of edit_distance operation.
/// </summary>
public record EditDistanceResult(int Distance);

/// <summary>
/// Result of count_approximate_occurrences operation.
/// </summary>
public record CountApproximateOccurrencesResult(int Count);
