using SuffixTree.Indexing;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>
/// Test-only oracle engine: the in-memory Ukkonen tree seen through the index contract. It uses nothing but the public <see cref="ISuffixTree"/> API,
/// so every answer is the answer of the existing, independently validated tree. A locus is a text occurrence: <see cref="IndexLocus.Lo"/> is the
/// leftmost start position of the locus string, <see cref="IndexLocus.Hi"/> the end (exclusive) of that occurrence, <see cref="IndexLocus.Length"/> the
/// length. Only the four core operations are supported (Count, ExtendRight, LazyLocate, TextAccess); every other operation is
/// <see cref="IndexSupport.Unsupported"/> and throws. Cost is O(n) per call, which is irrelevant for an oracle.
/// </summary>
public sealed class UkkonenIndexAdapter : ISuffixIndex
{
    private readonly string _text;
    private readonly UkkonenTree _tree;
    private readonly StringIndexText _indexText;

    /// <summary>Builds the Ukkonen tree of <paramref name="text"/>.</summary>
    public UkkonenIndexAdapter(string text)
    {
        ArgumentNullException.ThrowIfNull(text);
        _text = text;
        _tree = UkkonenTree.Build(text);
        _indexText = new StringIndexText(text);
    }

    /// <inheritdoc />
    public string Engine => "ukkonen-adapter";

    /// <inheritdoc />
    public long TextLength => _text.Length;

    /// <inheritdoc />
    public IIndexText Text => _indexText;

    /// <inheritdoc />
    public IndexLocus Root => new(0, 0, 0, 0);

    /// <inheritdoc />
    public IndexSupport Support(IndexOperation operation)
    {
        ThrowIfUndefined(operation);
        return operation is IndexOperation.Count or IndexOperation.ExtendRight or IndexOperation.LazyLocate or IndexOperation.TextAccess
            ? IndexSupport.Native
            : IndexSupport.Unsupported;
    }

    /// <inheritdoc />
    public string? SupportNote(IndexOperation operation)
    {
        ThrowIfUndefined(operation);
        return Support(operation) == IndexSupport.Unsupported ? $"the Ukkonen adapter does not provide {operation}" : null;
    }

    /// <inheritdoc />
    public long Count(in IndexLocus locus) => _tree.CountOccurrences(PatternOf(locus));

    /// <inheritdoc />
    public bool TryExtendRight(in IndexLocus locus, int symbol, out IndexLocus result)
    {
        if ((uint)symbol > char.MaxValue)
        {
            result = default;
            return false;
        }
        return TryExtend(PatternOf(locus), [(char)symbol], out result);
    }

    /// <inheritdoc />
    public bool TryExtendRight(in IndexLocus locus, ReadOnlySpan<char> pattern, out IndexLocus result)
        => TryExtend(PatternOf(locus), pattern, out result);

    /// <inheritdoc />
    public bool TryExtendLeft(in IndexLocus locus, int symbol, out IndexLocus result)
        => throw new NotSupportedException($"engine '{Engine}' does not support {IndexOperation.ExtendLeft}");

    /// <inheritdoc />
    public int Locate(in IndexLocus locus, long skip, Span<long> positions)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(skip);
        IReadOnlyList<int> all = _tree.FindAllOccurrences(PatternOf(locus));
        int written = 0;
        for (long i = skip; i < all.Count && written < positions.Length; i++)
            positions[written++] = all[(int)i];
        return written;
    }

    private string PatternOf(in IndexLocus locus) => _text.Substring((int)locus.Lo, (int)locus.Length);

    private bool TryExtend(string head, ReadOnlySpan<char> tail, out IndexLocus result)
    {
        string pattern = string.Concat(head, tail);
        IReadOnlyList<int> positions = _tree.FindAllOccurrences(pattern);
        if (positions.Count == 0 && pattern.Length > 0)
        {
            result = default;
            return false;
        }
        int first = pattern.Length == 0 ? 0 : positions.Min();
        result = new IndexLocus(first, first + pattern.Length, pattern.Length, 0);
        return true;
    }

    private static void ThrowIfUndefined(IndexOperation operation)
    {
        if (!Enum.IsDefined(operation)) throw new ArgumentOutOfRangeException(nameof(operation), operation, "not an IndexOperation");
    }
}

/// <summary>An <see cref="IIndexText"/> over a string (test support).</summary>
public sealed class StringIndexText(string text) : IIndexText
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
