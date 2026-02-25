using System;

namespace SuffixTree;

/// <summary>
/// Optional contract for text sources that can compare a pattern at a given offset
/// without allocating an intermediate slice.
/// </summary>
public interface ITextPatternMatcher
{
    /// <summary>
    /// Compares <paramref name="pattern"/> with source content starting at <paramref name="start"/>.
    /// Returns false when the requested range is out of bounds.
    /// </summary>
    bool SequenceEqualAt(int start, ReadOnlySpan<char> pattern);
}
