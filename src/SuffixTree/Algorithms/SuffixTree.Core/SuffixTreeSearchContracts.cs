using System;
using System.Collections.Generic;

namespace SuffixTree;

/// <summary>
/// Shared contract helpers for suffix-tree search overloads.
/// Keeps null/empty-pattern behavior aligned across implementations.
/// </summary>
public static class SuffixTreeSearchContracts
{
    /// <summary>
    /// Ensures that a string argument is not null.
    /// </summary>
    public static void EnsureNotNull(string? value, string paramName)
        => ArgumentNullException.ThrowIfNull(value, paramName);

    /// <summary>
    /// Builds the canonical index list for empty-pattern matches: all valid start positions.
    /// </summary>
    public static IReadOnlyList<int> BuildAllStartPositions(int textLength)
    {
        var all = new List<int>(textLength);
        for (int i = 0; i < textLength; i++)
            all.Add(i);
        return all;
    }
}
