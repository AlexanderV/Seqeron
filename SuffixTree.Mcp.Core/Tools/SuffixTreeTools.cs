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
}

/// <summary>
/// Result of suffix_tree_contains operation.
/// </summary>
/// <param name="Found">Whether the pattern was found in the text.</param>
public record SuffixTreeContainsResult(bool Found);
