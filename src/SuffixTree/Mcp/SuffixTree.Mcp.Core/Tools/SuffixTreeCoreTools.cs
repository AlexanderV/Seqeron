using System.ComponentModel;
using ModelContextProtocol.Server;

namespace SuffixTree.Mcp.Core.Tools;

/// <summary>
/// MCP tools for pure suffix-tree operations.
/// </summary>
[McpServerToolType]
public class SuffixTreeCoreTools
{
    // Utility holder for static MCP tools; never instantiated (S1118).
    private SuffixTreeCoreTools() { }

    /// <summary>
    /// Check if a pattern exists in text using suffix tree.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_contains", Title = "Suffix Tree — Contains", ReadOnly = true)]
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
    [McpServerTool(Name = "suffix_tree_count", Title = "Suffix Tree — Count Occurrences", ReadOnly = true)]
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
    [McpServerTool(Name = "suffix_tree_find_all", Title = "Suffix Tree — Find All Positions", ReadOnly = true)]
    [Description("Find all 0-based start positions (including overlapping occurrences) where a pattern occurs in text using suffix tree, in ascending order.")]
    public static SuffixTreeFindAllResult SuffixTreeFindAll(
        [Description("The text to search in")] string text,
        [Description("The pattern to find")] string pattern)
    {
        if (string.IsNullOrEmpty(text))
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        if (pattern == null)
            throw new ArgumentException("Pattern cannot be null", nameof(pattern));

        var tree = global::SuffixTree.SuffixTree.Build(text);
        // The tree reports leaves in traversal order; sort so the tool output is deterministic
        // and ascending (Rosalind SUBS convention).
        var positions = tree.FindAllOccurrences(pattern).ToArray();
        Array.Sort(positions);
        return new SuffixTreeFindAllResult(positions);
    }

    /// <summary>
    /// Find the longest repeated substring in text.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_lrs", Title = "Suffix Tree — Longest Repeated Substring", ReadOnly = true)]
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
    [McpServerTool(Name = "suffix_tree_lcs", Title = "Suffix Tree — Longest Common Substring", ReadOnly = true)]
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
    /// Find every longest repeated substring (all length ties) with all positions.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_all_lrs", Title = "Suffix Tree — All Longest Repeated Substrings", ReadOnly = true)]
    [Description("Find every distinct longest repeated substring of text (all length ties, unlike suffix_tree_lrs which returns one representative), each with all 0-based start positions in ascending order (occurrences may overlap); substrings ordered by first occurrence. Empty when no character repeats.")]
    public static SuffixTreeAllLrsResult SuffixTreeAllLrs(
        [Description("The text to analyze")] string text)
    {
        if (string.IsNullOrEmpty(text))
            throw new ArgumentException("Text cannot be null or empty", nameof(text));

        var tree = global::SuffixTree.SuffixTree.Build(text);
        var all = tree.FindAllLongestRepeatedSubstrings();
        var items = all.Select(r => new RepeatedSubstringItem(r.Substring, r.Positions.ToArray())).ToArray();
        return new SuffixTreeAllLrsResult(items, items.Length == 0 ? 0 : items[0].Substring.Length);
    }

    /// <summary>
    /// Find all maximal exact matches (MEMs) between a reference text and a query (MUMmer -maxmatch).
    /// </summary>
    [McpServerTool(Name = "suffix_tree_find_mems", Title = "Suffix Tree — Maximal Exact Matches (MEMs)", ReadOnly = true)]
    [Description("Find all maximal exact matches (MEMs) of length >= minLength between a reference text and a query: every left- and right-maximal match with every reference occurrence (MUMmer 3 'mummer -maxmatch -l minLength', forward strand). Positions are 0-based; sorted by query position, then text position.")]
    public static SuffixTreeMaximalMatchesResult SuffixTreeFindMems(
        [Description("The reference text (the suffix tree is built on it)")] string text,
        [Description("The query string matched against the reference")] string query,
        [Description("Minimum match length (>= 1; MUMmer default 20)")] int minLength = 20)
    {
        if (string.IsNullOrEmpty(text))
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));
        if (minLength < 1)
            throw new ArgumentOutOfRangeException(nameof(minLength), "minLength must be >= 1");

        var tree = global::SuffixTree.SuffixTree.Build(text);
        return ToMaximalMatchesResult(tree.FindMaximalExactMatches(query, minLength));
    }

    /// <summary>
    /// Find maximal unique matches (MUMs) between a reference text and a query (MUMmer -mum / -mumreference).
    /// </summary>
    [McpServerTool(Name = "suffix_tree_find_mums", Title = "Suffix Tree — Maximal Unique Matches (MUMs)", ReadOnly = true)]
    [Description("Find maximal unique matches (MUMs) of length >= minLength: MEMs whose string occurs exactly once in the reference and, with uniqueness 'both', exactly once in the query (MUMmer 3 'mummer -mum'); 'reference' requires uniqueness in the reference only ('mummer -mumreference'). Forward strand, 0-based; sorted by query position, then text position.")]
    public static SuffixTreeMaximalMatchesResult SuffixTreeFindMums(
        [Description("The reference text (the suffix tree is built on it)")] string text,
        [Description("The query string matched against the reference")] string query,
        [Description("Minimum match length (>= 1; MUMmer default 20)")] int minLength = 20,
        [Description("Uniqueness requirement: 'both' (unique in reference and query, -mum; default) or 'reference' (unique in reference only, -mumreference)")] string uniqueness = "both")
    {
        if (string.IsNullOrEmpty(text))
            throw new ArgumentException("Text cannot be null or empty", nameof(text));
        if (string.IsNullOrEmpty(query))
            throw new ArgumentException("Query cannot be null or empty", nameof(query));
        if (minLength < 1)
            throw new ArgumentOutOfRangeException(nameof(minLength), "minLength must be >= 1");

        var mode = (uniqueness ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "both" => MumUniqueness.Both,
            "reference" => MumUniqueness.Reference,
            _ => throw new ArgumentException("Uniqueness must be 'both' or 'reference'", nameof(uniqueness)),
        };

        var tree = global::SuffixTree.SuffixTree.Build(text);
        return ToMaximalMatchesResult(tree.FindMaximalUniqueMatches(query, minLength, mode));
    }

    private static SuffixTreeMaximalMatchesResult ToMaximalMatchesResult(
        IReadOnlyList<(int PositionInText, int PositionInQuery, int Length)> matches)
        => new(matches.Select(m => new MaximalMatchItem(m.PositionInText, m.PositionInQuery, m.Length)).ToArray());

    /// <summary>
    /// Get statistics about a suffix tree built from text.
    /// </summary>
    [McpServerTool(Name = "suffix_tree_stats", Title = "Suffix Tree — Statistics", ReadOnly = true)]
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
