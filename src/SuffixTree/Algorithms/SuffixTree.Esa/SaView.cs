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
    private readonly int _count;
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
    internal SaView(int* pointer, int count, bool isChecked)
    {
        _array = null;
        _pointer = pointer;
        _count = count;
        _checked = isChecked;
    }

    /// <summary>Number of entries (<c>n + 1</c>).</summary>
    internal int Length => _count;

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
    internal void CopyTo(int start, Span<int> destination)
    {
        if (start < 0 || destination.Length > _count - start) ThrowOutside(start);
        if (_array is not null) _array.AsSpan(start, destination.Length).CopyTo(destination);
        else new ReadOnlySpan<int>(_pointer + start, destination.Length).CopyTo(destination);
    }

    /// <summary>The entries <c>[start, start + count)</c> as a new array.</summary>
    internal int[] ToArray(int start, int count)
    {
        var result = new int[count];
        CopyTo(start, result);
        return result;
    }

    /// <summary>The pointer to the first entry (null for an array); the navigation pack reads through it.</summary>
    internal int* Pointer => _pointer;

    /// <summary>The wrapped array, or null for a mapped section.</summary>
    internal int[]? Array => _array;

    [MethodImpl(MethodImplOptions.NoInlining)]
    private void ThrowOutside(long index) => throw new InvalidDataException($"the suffix array is read at rank {index}, but the section has {_count} entries: the file is damaged");
}
