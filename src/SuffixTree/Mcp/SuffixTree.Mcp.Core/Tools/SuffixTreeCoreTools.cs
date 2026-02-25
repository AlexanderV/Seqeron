using System.ComponentModel;
using ModelContextProtocol.Server;

namespace SuffixTree.Mcp.Core.Tools;

/// <summary>
/// MCP tools for pure suffix-tree operations.
/// </summary>
[McpServerToolType]
public class SuffixTreeCoreTools
{
    /// <summary>
    /// Check if a pattern exists in text using suffix tree.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_contains")]
    [Description("Check if a pattern exists in text using suffix tree. Returns true if pattern is found, false otherwise.")]
    public static SuffixTreeContainsResult SuffixTreeContains(
        [Description("The text to search in")] string text,
        [Description("The pattern to search for")] string pattern)
    {
        if (string.IsNullOrEmpty(text))
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        if (pattern == null)
            throw new ArgumentException("Pattern cannot be null", nameof(pattern));

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
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        if (pattern == null)
            throw new ArgumentException("Pattern cannot be null", nameof(pattern));

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
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        if (pattern == null)
            throw new ArgumentException("Pattern cannot be null", nameof(pattern));

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
            throw new ArgumentException("Text cannot be null or empty", nameof(text));

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
            throw new ArgumentException("Text1 cannot be null or empty", nameof(text1));
        if (string.IsNullOrEmpty(text2))
            throw new ArgumentException("Text2 cannot be null or empty", nameof(text2));

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
            throw new ArgumentException("Text cannot be null or empty", nameof(text));

        var tree = global::SuffixTree.SuffixTree.Build(text);
        return new SuffixTreeStatsResult(tree.NodeCount, tree.LeafCount, tree.MaxDepth, text.Length);
    }
}
