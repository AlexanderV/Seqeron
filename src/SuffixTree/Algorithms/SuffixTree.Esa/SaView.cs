using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// Read access to the suffix array of an <see cref="EsaIndex"/>, whichever way it is stored (work package 1b.7f): a managed <c>int[]</c> (an index built in memory) or <c>n + 1</c> little-endian 32-bit entries inside a mapped
/// SESA container (an index opened from a file). A struct, so that it costs the caller one reference, one pointer and a count. The container's entries are untrusted: in <see cref="EsaAccessMode.Checked"/> mode (the default) an index outside
/// <c>0..n</c> is an <see cref="InvalidDataException"/> and never a read outside the section; a value that a read returns can still be garbage (it is a position or a rank that the caller validates against its own bounds).
/// </summary>
internal readonly unsafe struct SaView
{
    private readonly int[]? _array;
    private readonly int* _pointer;
    private readonly long _count;
    private readonly bool _checked;

    /// <summary>A view of an array.</summary>
    internal SaView(int[] array)
    {
        _array = array;
        _pointer = null;
        _count = array.Length;
        _checked = true;
    }

    /// <summary>A view of <paramref name="count"/> entries at <paramref name="pointer"/> (the caller keeps the memory alive and unchanged).</summary>
    internal SaView(int* pointer, long count, bool isChecked)
    {
        _array = null;
        _pointer = pointer;
        _count = count;
        _checked = isChecked;
    }

    /// <summary>Number of entries (<c>n + 1</c>).</summary>
    internal long Length => _count;

    /// <summary>The entry at <paramref name="rank"/>.</summary>
    /// <exception cref="InvalidDataException">The rank is outside the section (a damaged table sent the read there); checked mode only for a file.</exception>
    internal int this[int rank]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            if (_array is not null) return _array[rank];
            if (_checked && (uint)rank >= (uint)_count) ThrowOutside(rank);
            return _pointer[rank];
        }
    }

    /// <summary>The same read for a 64-bit rank (the contract's).</summary>
    internal int At(long rank)
    {
        if ((ulong)rank >= (ulong)_count) ThrowOutside(rank);
        return this[(int)rank];
    }

    /// <summary>Copies <c>destination.Length</c> entries that start at <paramref name="start"/>; the range is always checked (it is a bulk read).</summary>
    internal void CopyTo(long start, Span<int> destination)
    {
        if (start < 0 || destination.Length > _count - start) ThrowOutside(start);
        if (_array is not null) _array.AsSpan((int)start, destination.Length).CopyTo(destination);
        else new ReadOnlySpan<int>(_pointer + start, destination.Length).CopyTo(destination);
    }

    /// <summary>The entries <c>[start, start + count)</c> as a new array.</summary>
    internal int[] ToArray(long start, int count)
    {
        var result = new int[count];
        CopyTo(start, result);
        return result;
    }

#pragma warning disable S2325 // stubs of the tests-first step
    // ---- work package S2.3, the long members (STUBS of the tests-first step)

    /// <summary>A view of <paramref name="count"/> entries of <paramref name="posWidth"/> (4 or 5) bytes at <paramref name="section"/>.</summary>
    internal static SaView FromMapped(byte* section, long count, int posWidth, bool isChecked) => throw new NotImplementedException();

    /// <summary>The entry at a 64-bit rank (the stored width of the view; the checks of the access mode).</summary>
    internal long this[long rank] => throw new NotImplementedException();

    /// <summary>The stored width in bytes: 4 for an array and a 32-bit section, 5 for a 40-bit section.</summary>
    internal int StoredWidth => throw new NotImplementedException();

    /// <summary>Copies <c>destination.Length</c> entries that start at the 64-bit rank <paramref name="start"/>.</summary>
    internal void CopyTo(long start, Span<long> destination) => throw new NotImplementedException();

    /// <summary>The entries <c>[start, start + count)</c> as a new array of 64-bit positions.</summary>
    internal long[] ToLongArray(long start, int count) => throw new NotImplementedException();

    /// <summary>The pointer to the first entry (null for an array); the navigation pack reads through it.</summary>
    internal int* Pointer => _pointer;

    /// <summary>The wrapped array, or null for a mapped section.</summary>
    internal int[]? Array => _array;

    [MethodImpl(MethodImplOptions.NoInlining)]
    private void ThrowOutside(long index) => throw new InvalidDataException($"the suffix array is read at rank {index}, but the section has {_count} entries: the file is damaged");
}

/// <summary>
/// The suffix array over a position reader (ADR-005 section 3.3, work package S2.3): the generic form of <see cref="SaView"/> whose every number is 64-bit. The stored width is the type argument (<see cref="Pos32{TAccess}"/>,
/// <see cref="Pos40{TAccess}"/>, <see cref="PosInt"/>, a closed form in the tests), so that a loop over the view is specialised for the width.
/// </summary>
internal readonly struct SaView<TPos>
    where TPos : struct, IPosReader
{
    private readonly TPos _reader;

    internal SaView(TPos reader) => _reader = reader;

    /// <summary>Number of entries (<c>n + 1</c>).</summary>
    internal long Length => _reader.Count;

    /// <summary>The reader.</summary>
    internal TPos Reader => _reader;

    /// <summary>The entry at <paramref name="rank"/>.</summary>
    internal long this[long rank] => throw new NotImplementedException();

    /// <summary>The same read, always bounded (the contract's).</summary>
    internal long At(long rank) => throw new NotImplementedException();

    /// <summary>Copies <c>destination.Length</c> entries that start at <paramref name="start"/>.</summary>
    internal void CopyTo(long start, Span<long> destination) => throw new NotImplementedException();

    /// <summary>Copies <c>destination.Length</c> entries that start at <paramref name="start"/> as 32-bit positions.</summary>
    internal void CopyTo(long start, Span<int> destination) => throw new NotImplementedException();

    /// <summary>The entries <c>[start, start + count)</c> as a new array.</summary>
    internal long[] ToArray(long start, int count) => throw new NotImplementedException();
}
