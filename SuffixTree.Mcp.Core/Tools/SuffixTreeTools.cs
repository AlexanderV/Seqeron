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
