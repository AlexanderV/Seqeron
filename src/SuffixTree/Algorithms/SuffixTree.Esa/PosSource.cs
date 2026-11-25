namespace SuffixTree.Esa;

/// <summary>
/// An array of positions, ranks or LCP values at 64-bit indexes, whichever way it is stored (work package S2.5): a view of a suffix array or of a section of a container (<see cref="SaView"/>: an <c>int[]</c>, 4-byte or 5-byte entries) or
/// any <see cref="IPosReader"/> (a closed form of a test). The long navigation engine reads the suffix array, PSI and the persisted RMQ levels through it: not the hot path of any gate row (those are the 32-bit navigator and the descents), so
/// one virtual read costs nothing that matters and one engine serves every width. The checks are the reader's (Checked mode bounds every read).
/// </summary>
internal abstract class PosSource
{
    /// <summary>Number of elements.</summary>
    internal abstract long Count { get; }

    /// <summary>The element at <paramref name="index"/>.</summary>
    /// <exception cref="InvalidDataException">The index is outside the section (a checked reader).</exception>
    internal abstract long this[long index] { get; }

    /// <summary>The elements of a suffix-array-shaped view.</summary>
    internal static PosSource Of(SaView view) => new ViewSource(view);

    /// <summary>The elements of a reader (its access mode is the reader's).</summary>
    internal static PosSource Of<TPos>(TPos reader)
        where TPos : struct, IPosReader => new ReaderSource<TPos>(reader);

    private sealed class ViewSource(SaView view) : PosSource
    {
        internal override long Count => view.Length;

        internal override long this[long index] => view[index];
    }

    private sealed class ReaderSource<TPos>(TPos reader) : PosSource
        where TPos : struct, IPosReader
    {
        internal override long Count => reader.Count;

        internal override long this[long index] => (long)reader.Get(index);
    }
}
