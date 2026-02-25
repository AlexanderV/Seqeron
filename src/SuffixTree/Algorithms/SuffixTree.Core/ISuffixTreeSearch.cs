using System;
using System.Collections.Generic;

namespace SuffixTree;

/// <summary>
/// Search and basic metadata contract for suffix trees.
/// </summary>
public interface ISuffixTreeSearch
{
    /// <summary>
    /// Gets the original text source that this suffix tree was built from.
    /// </summary>
    ITextSource Text { get; }

    /// <summary>
    /// Gets the total number of nodes in the tree (including root).
    /// </summary>
    int NodeCount { get; }

    /// <summary>
    /// Gets the number of leaf nodes in the tree.
    /// </summary>
    int LeafCount { get; }

    /// <summary>
    /// Gets the maximum depth of the tree (longest path from root to leaf in characters).
    /// </summary>
    int MaxDepth { get; }

    /// <summary>
    /// Gets a value indicating whether this tree is empty (built from an empty string).
    /// </summary>
    bool IsEmpty { get; }

    /// <summary>
    /// Checks if the specified string is a substring of the tree content.
    /// </summary>
    bool Contains(string value);

    /// <summary>
    /// Checks if the specified character span is a substring of the tree content.
    /// </summary>
    bool Contains(ReadOnlySpan<char> value);

    /// <summary>
    /// Finds all starting positions where the pattern occurs in the original string.
    /// </summary>
    IReadOnlyList<int> FindAllOccurrences(string pattern);

    /// <summary>
    /// Finds all starting positions where the pattern occurs in the original string.
    /// </summary>
    IReadOnlyList<int> FindAllOccurrences(ReadOnlySpan<char> pattern);

    /// <summary>
    /// Counts the number of occurrences of a pattern in the text.
    /// </summary>
    int CountOccurrences(string pattern);

    /// <summary>
    /// Counts the number of occurrences of a pattern in the text.
    /// </summary>
    int CountOccurrences(ReadOnlySpan<char> pattern);
}
