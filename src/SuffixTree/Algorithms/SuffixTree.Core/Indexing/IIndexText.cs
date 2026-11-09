using System;

namespace SuffixTree.Indexing;

/// <summary>64-bit counterpart of <see cref="ITextSource"/>: read access to the indexed text of an engine.</summary>
public interface IIndexText
{
    /// <summary>Number of symbols n.</summary>
    long Length { get; }

    /// <summary>
    /// The UTF-16 code unit at <paramref name="position"/> (0..65535), or -1 (the terminator) for every <paramref name="position"/> at or beyond
    /// <see cref="Length"/>.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="position"/> is negative.</exception>
    int SymbolAt(long position);

    /// <summary>
    /// Copies the symbols <c>[start, start + destination.Length)</c> clipped to the text into <paramref name="destination"/>.
    /// </summary>
    /// <returns>The number of code units written: <c>min(destination.Length, Length - start)</c> (0 when <paramref name="start"/> equals <see cref="Length"/>).</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="start"/> is negative or larger than <see cref="Length"/>.</exception>
    int CopyTo(long start, Span<char> destination);
}
