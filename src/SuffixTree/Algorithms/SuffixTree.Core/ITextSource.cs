using System;
using System.Collections.Generic;

namespace SuffixTree;

/// <summary>
/// Defines an abstraction for a source of text characters.
/// Allows the suffix tree to operate on strings, memory-mapped files, or other storage.
/// </summary>
public interface ITextSource : IEnumerable<char>
{
    /// <summary>
    /// Gets the number of characters in the source.
    /// </summary>
    int Length { get; }

    /// <summary>
    /// Gets the character at the specified index.
    /// </summary>
    /// <param name="index">The zero-based index of the character to retrieve.</param>
    /// <returns>The character at the specified index.</returns>
    char this[int index] { get; }

    /// <summary>
    /// Retrieves a substring from this instance.
    /// </summary>
    /// <param name="start">The zero-based starting character position of a substring in this instance.</param>
    /// <param name="length">The number of characters in the substring.</param>
    /// <returns>A string that is equivalent to the substring of length <paramref name="length"/> that begins at <paramref name="start"/>.</returns>
    string Substring(int start, int length);

    /// <summary>
    /// Retrieves a character span from this instance.
    /// </summary>
    /// <param name="start">The zero-based starting character position.</param>
    /// <param name="length">The number of characters.</param>
    /// <returns>A read-only span of characters.</returns>
    ReadOnlySpan<char> Slice(int start, int length);
}
