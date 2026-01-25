using System.ComponentModel;
using ModelContextProtocol.Server;

namespace SuffixTree.Mcp.Core.Tools;

/// <summary>
/// MCP tools for core suffix tree operations.
/// </summary>
[McpServerToolType]
public static class SuffixTreeTools
{
    /// <summary>
    /// Check if a pattern exists in text using suffix tree.
    /// </summary>
    /// <param name="text">The text to search in.</param>
    /// <param name="pattern">The pattern to search for.</param>
    /// <returns>Result indicating whether pattern was found.</returns>
    [McpServerTool(Name = "suffix_tree_contains")]
    [Description("Check if a pattern exists in text using suffix tree. Returns true if pattern is found, false otherwise.")]
    public static SuffixTreeContainsResult SuffixTreeContains(
        [Description("The text to search in")] string text,
        [Description("The pattern to search for")] string pattern)
    {
        if (string.IsNullOrEmpty(text))
        {
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        }

        if (pattern == null)
        {
            throw new ArgumentException("Pattern cannot be null", nameof(pattern));
        }

        var tree = global::SuffixTree.SuffixTree.Build(text);
        var found = tree.Contains(pattern);

        return new SuffixTreeContainsResult(found);
    }

    /// <summary>
    /// Count occurrences of a pattern in text using suffix tree.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_count")]
    [Description("Count the number of occurrences of a pattern in text using suffix tree.")]
    public static SuffixTreeCountResult SuffixTreeCount(
        [Description("The text to search in")] string text,
        [Description("The pattern to count")] string pattern)
    {
        if (string.IsNullOrEmpty(text))
        {
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        }

        if (pattern == null)
        {
            throw new ArgumentException("Pattern cannot be null", nameof(pattern));
        }

        var tree = global::SuffixTree.SuffixTree.Build(text);
        var count = tree.CountOccurrences(pattern);

        return new SuffixTreeCountResult(count);
    }

    /// <summary>
    /// Find all positions where a pattern occurs in text.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_find_all")]
    [Description("Find all positions where a pattern occurs in text using suffix tree.")]
    public static SuffixTreeFindAllResult SuffixTreeFindAll(
        [Description("The text to search in")] string text,
        [Description("The pattern to find")] string pattern)
    {
        if (string.IsNullOrEmpty(text))
        {
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        }

        if (pattern == null)
        {
            throw new ArgumentException("Pattern cannot be null", nameof(pattern));
        }

        var tree = global::SuffixTree.SuffixTree.Build(text);
        var positions = tree.FindAllOccurrences(pattern);

        return new SuffixTreeFindAllResult(positions.ToArray());
    }

    /// <summary>
    /// Find the longest repeated substring in text.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_lrs")]
    [Description("Find the longest repeated substring in text using suffix tree.")]
    public static SuffixTreeLrsResult SuffixTreeLrs(
        [Description("The text to analyze")] string text)
    {
        if (string.IsNullOrEmpty(text))
        {
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        }

        var tree = global::SuffixTree.SuffixTree.Build(text);
        var lrs = tree.LongestRepeatedSubstring();

        return new SuffixTreeLrsResult(lrs, lrs.Length);
    }

    /// <summary>
    /// Find the longest common substring between two texts.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_lcs")]
    [Description("Find the longest common substring between two texts using suffix tree.")]
    public static SuffixTreeLcsResult SuffixTreeLcs(
        [Description("The first text")] string text1,
        [Description("The second text")] string text2)
    {
        if (string.IsNullOrEmpty(text1))
        {
            throw new ArgumentException("Text1 cannot be null or empty", nameof(text1));
        }

        if (string.IsNullOrEmpty(text2))
        {
            throw new ArgumentException("Text2 cannot be null or empty", nameof(text2));
        }

        var tree = global::SuffixTree.SuffixTree.Build(text1);
        var lcs = tree.LongestCommonSubstring(text2);

        return new SuffixTreeLcsResult(lcs, lcs.Length);
    }

    /// <summary>
    /// Get statistics about a suffix tree built from text.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_stats")]
    [Description("Get statistics about a suffix tree: node count, leaf count, max depth, and text length.")]
    public static SuffixTreeStatsResult SuffixTreeStats(
        [Description("The text to analyze")] string text)
    {
        if (string.IsNullOrEmpty(text))
        {
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        }

        var tree = global::SuffixTree.SuffixTree.Build(text);

        return new SuffixTreeStatsResult(tree.NodeCount, tree.LeafCount, tree.MaxDepth, text.Length);
    }

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

        var dna = new global::SuffixTree.Genomics.DnaSequence(sequence);
        var result = global::SuffixTree.Genomics.GenomicAnalyzer.FindLongestRepeat(dna);

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

        var dna1 = new global::SuffixTree.Genomics.DnaSequence(sequence1);
        var dna2 = new global::SuffixTree.Genomics.DnaSequence(sequence2);
        var result = global::SuffixTree.Genomics.GenomicAnalyzer.FindLongestCommonRegion(dna1, dna2);

        return new FindLongestCommonRegionResult(result.Sequence, result.PositionInFirst, result.PositionInSecond, result.Length);
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
