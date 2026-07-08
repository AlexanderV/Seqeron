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
    // S3236: paramName is deliberately forwarded from this contract helper's own caller, so the
    // thrown exception names the original argument — not overriding CallerArgumentExpression by accident.
#pragma warning disable S3236
    public static void EnsureNotNull(string? value, string paramName)
        => ArgumentNullException.ThrowIfNull(value, paramName);
#pragma warning restore S3236

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
