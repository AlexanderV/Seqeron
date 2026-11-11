using SuffixTree.Indexing;

namespace SuffixTree.Esa;

/// <summary>
/// <see cref="IIndexText"/> over the string of an <see cref="EsaIndex"/> (work package 1b.5): the 64-bit view of the text that the index contract hands out. The string is shared, not copied.
/// </summary>
internal sealed class EsaIndexText(string text) : IIndexText
{
    /// <inheritdoc />
    public long Length => text.Length;

    /// <inheritdoc />
    public int SymbolAt(long position)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(position);
        return position >= text.Length ? -1 : text[(int)position];
    }

    /// <inheritdoc />
    public int CopyTo(long start, Span<char> destination)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(start);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(start, text.Length);
        int count = (int)Math.Min(destination.Length, text.Length - start);
        text.AsSpan((int)start, count).CopyTo(destination);
        return count;
    }
}
