using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// Read access to the suffix array of an <see cref="EsaIndex"/>, whichever way it is stored (work package 1b.7f, widened by S2.3): a managed <c>int[]</c> (an index built in memory) or <c>n + 1</c> little-endian entries of 4 bytes
/// (<c>POS_WIDTH</c> 4) or of 5 bytes (<c>POS_WIDTH</c> 5, ADR-005 section 6.1) inside a mapped SESA container (an index opened from a file). A struct, so that it costs the caller one reference, one pointer and a count. The container's entries
/// are untrusted: in <see cref="EsaAccessMode.Checked"/> mode (the default) an index outside <c>0..n</c> is an <see cref="InvalidDataException"/> and never a read outside the section; a value that a read returns can still be garbage (it is a position or a
/// rank that the caller validates against its own bounds). Every number of the long members is 64-bit (<see cref="Length"/>, <c>this[long]</c>, <see cref="At"/>, the copies); the <see cref="int"/> members are what the 32-bit engine of S1 reads until
/// S2.4 and S2.5 and are unchanged on an array and on a 4-byte section of at most <see cref="int.MaxValue"/> entries (the hot path: the same test of the flags tells the bridged forms apart); on a 5-byte section and on a 4-byte section of more entries they are the checked bridge (a position above <see cref="int.MaxValue"/> is an <see cref="InvalidDataException"/>, never a wrapped number). The generic form over a
/// position reader is <see cref="SaView{TPos}"/>; this one is the non-generic holder that a container reader hands out (a switch on the form on every long member, none on the 32-bit read of an array or a 4-byte section).
/// </summary>
internal readonly unsafe struct SaView
{
    private const byte Checked = 1;
    private const byte FiveBytes = 2;

    /// <summary>A 4-byte section of more than <see cref="int.MaxValue"/> entries: an entry may be a position above <see cref="int.MaxValue"/> (a text of 2^31 symbols and more), which the 32-bit read must refuse, not wrap.</summary>
    private const byte FourBytesWide = 4;

    /// <summary>The sections whose 32-bit read is the checked bridge (the same test of the flags as the one that told the 5-byte section from the 4-byte one).</summary>
    private const byte Bridged = FiveBytes | FourBytesWide;

    private readonly int[]? _array;
    private readonly int* _pointer;
    private readonly ulong _countFlags; // bits 0-55: the number of entries; bits 56-63: the flags (24 bytes in all, as in S1: the struct is copied by value on every call that reads the suffix array)

    private const int FlagShift = 56;
    private const ulong CountMask = (1UL << FlagShift) - 1;

    private long Count => (long)(_countFlags & CountMask);

    private byte Flags => (byte)(_countFlags >> FlagShift);

    private static ulong Pack(long count, byte flags) => (ulong)count | ((ulong)flags << FlagShift);

    /// <summary>A view of an array.</summary>
    internal SaView(int[] array)
    {
        ArgumentNullException.ThrowIfNull(array);
        _array = array;
        _pointer = null;
        _countFlags = Pack(array.Length, Checked);
    }

    /// <summary>A view of <paramref name="count"/> 4-byte entries at <paramref name="pointer"/> (the caller keeps the memory alive and unchanged).</summary>
    internal SaView(int* pointer, long count, bool isChecked)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        _array = null;
        _pointer = pointer;
        _countFlags = Pack(count, FlagsOf(isChecked, 4, count));
    }

    private static byte FlagsOf(bool isChecked, int posWidth, long count)
    {
        int flags = isChecked ? Checked : 0;
        if (posWidth == 5) flags |= FiveBytes;
        else if (count > int.MaxValue) flags |= FourBytesWide;
        return (byte)flags;
    }

    private SaView(int* pointer, long count, byte flags)
    {
        _array = null;
        _pointer = pointer;
        _countFlags = Pack(count, flags);
    }

    /// <summary>A view of <paramref name="count"/> entries of <paramref name="posWidth"/> (4 or 5) bytes at <paramref name="section"/> (the caller keeps the memory alive and unchanged; a 5-byte section has its 3 bytes of padding).</summary>
    /// <exception cref="ArgumentOutOfRangeException">The width is not 4 or 5, or the count is negative.</exception>
    internal static SaView FromMapped(byte* section, long count, int posWidth, bool isChecked)
    {
        if (posWidth is not (4 or 5)) throw new ArgumentOutOfRangeException(nameof(posWidth), posWidth, "the stored width of a position is 4 or 5 bytes");
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        return new SaView((int*)section, count, FlagsOf(isChecked, posWidth, count));
    }

    /// <summary>Number of entries (<c>n + 1</c>).</summary>
    internal long Length => Count;

    /// <summary>The stored width in bytes: 4 for an array and a 4-byte section, 5 for a 5-byte section.</summary>
    internal int StoredWidth => (Flags & FiveBytes) != 0 ? 5 : 4;

    /// <summary>The entry at <paramref name="rank"/> (the 32-bit engine's read).</summary>
    /// <exception cref="InvalidDataException">The rank is outside the section (a damaged table sent the read there; checked mode only for a file), or the entry of a 5-byte section is not a 32-bit number.</exception>
    internal int this[int rank]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            if (_array is not null) return _array[rank];
            if ((Flags & Checked) != 0 && (ulong)rank >= (ulong)Count) ThrowOutside(rank);
            if ((Flags & Bridged) != 0) return BridgedRead(rank);
            return _pointer[rank];
        }
    }

    /// <summary>The entry at a 64-bit rank: the checks of the access mode, the stored width of the view.</summary>
    /// <exception cref="InvalidDataException">The rank is outside the section (checked mode only for a file).</exception>
    internal long this[long rank]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get
        {
            if (_array is not null) return (uint)_array[rank];
            if ((Flags & Checked) != 0 && (ulong)rank >= (ulong)Count) ThrowOutside(rank);
            return (Flags & FiveBytes) != 0 ? (long)Load40(rank) : (uint)_pointer[rank];
        }
    }

    /// <summary>The same read, always bounded in both modes (the contract's).</summary>
    /// <exception cref="InvalidDataException">The rank is outside <c>0..n</c>.</exception>
    internal long At(long rank)
    {
        if ((ulong)rank >= (ulong)Count) ThrowOutside(rank);
        return this[rank];
    }

    /// <summary>Copies <c>destination.Length</c> entries that start at <paramref name="start"/>; the range is always checked (it is a bulk read).</summary>
    /// <exception cref="InvalidDataException">The range is outside the section, or an entry of a 5-byte section is not a 32-bit number.</exception>
    internal void CopyTo(long start, Span<int> destination)
    {
        if (start < 0 || destination.Length > Count - start) ThrowOutside(start);
        if (_array is not null) _array.AsSpan((int)start, destination.Length).CopyTo(destination);
        else if ((Flags & Bridged) == 0) new ReadOnlySpan<int>(_pointer + start, destination.Length).CopyTo(destination);
        else CopyBridged(start, destination);
    }

    /// <summary>The bulk copy of a bridged section (out of line, so that the copy of an array and of a 4-byte section stays the small code of S1: the locate of a short pattern ends in it).</summary>
    [MethodImpl(MethodImplOptions.NoInlining)]
    private void CopyBridged(long start, Span<int> destination)
    {
        for (int i = 0; i < destination.Length; i++) destination[i] = IntBridge.Value(ReadUnchecked(start + i), "a suffix array entry");
    }

    /// <summary>Copies <c>destination.Length</c> entries that start at the 64-bit rank <paramref name="start"/>; the range is always checked.</summary>
    /// <exception cref="InvalidDataException">The range is outside the section.</exception>
    internal void CopyTo(long start, Span<long> destination)
    {
        if (start < 0 || destination.Length > Count - start) ThrowOutside(start);
        for (int i = 0; i < destination.Length; i++) destination[i] = ReadUnchecked(start + i);
    }

    /// <summary>The entries <c>[start, start + count)</c> as a new array.</summary>
    internal int[] ToArray(long start, int count)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        int[] result = new int[count];
        CopyTo(start, result);
        return result;
    }

    /// <summary>
    /// <see cref="ToArray"/> for a long list, without zeroing the array first: every entry is written by the copy (a failed copy drops the array), which saves the clearing of the 6 KB that the locate of 1,500 hits returns. Out of line:
    /// the locate of a short pattern (a few entries, below <see cref="UninitializedFrom"/>) calls the small <see cref="ToArray"/> of S1, and a runtime that allocates a zeroed array below 2,048 bytes anyway would only add its dispatch.
    /// </summary>
    [MethodImpl(MethodImplOptions.NoInlining)]
    internal int[] ToArrayUninitialized(long start, int count)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        int[] result = GC.AllocateUninitializedArray<int>(count);
        CopyTo(start, result);
        return result;
    }

    /// <summary>The number of entries from which <see cref="ToArrayUninitialized"/> pays (2,048 bytes: the size below which the runtime zeroes the array of <see cref="GC.AllocateUninitializedArray{T}(int, bool)"/> anyway).</summary>
    internal const int UninitializedFrom = 512;

    /// <summary>The entries <c>[start, start + count)</c> as a new array of 64-bit positions.</summary>
    internal long[] ToLongArray(long start, int count)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        var result = new long[count];
        CopyTo(start, result);
        return result;
    }

    /// <summary>
    /// The pointer to the first entry (null for an array); the navigation pack and <see cref="MappedSearch"/> read through it.
    /// </summary>
    /// <exception cref="NotSupportedException">The section has 5-byte entries: no pointer to 32-bit entries exists (the 32-bit engine does not read it).</exception>
    internal int* Pointer => (Flags & FiveBytes) == 0 ? _pointer : throw new NotSupportedException("the suffix array has 5-byte entries (POS_WIDTH 5): there is no pointer to 32-bit entries; use the 64-bit members");

    /// <summary>The wrapped array, or null for a mapped section.</summary>
    internal int[]? Array => _array;

    private long ReadUnchecked(long rank)
    {
        if (_array is not null) return (uint)_array[rank];
        return (Flags & FiveBytes) != 0 ? (long)Load40(rank) : (uint)_pointer[rank];
    }

    /// <summary>The 32-bit read of a bridged section (out of line, so that the indexer that the navigator inlines stays the code of S1).</summary>
    [MethodImpl(MethodImplOptions.NoInlining)]
    private int BridgedRead(int rank) => IntBridge.Value(ReadUnchecked(rank), "a suffix array entry");

    private ulong Load40(long rank) => Unsafe.ReadUnaligned<ulong>((byte*)_pointer + (5 * rank)) & PosReader.Mask40;

    [MethodImpl(MethodImplOptions.NoInlining)]
    private void ThrowOutside(long index) => SaViewErrors.Outside(index, Count);
}

/// <summary>The failure that every suffix-array view shares.</summary>
internal static class SaViewErrors
{
    [MethodImpl(MethodImplOptions.NoInlining)]
    internal static void Outside(long index, long count) =>
        throw new InvalidDataException($"the suffix array is read at rank {index}, but the section has {count} entries: the file is damaged");
}

/// <summary>
/// The suffix array over a position reader (ADR-005 section 3.3, work package S2.3): the generic form of <see cref="SaView"/> whose every number is 64-bit. The stored width is the type argument (<see cref="Pos32{TAccess}"/>,
/// <see cref="Pos40{TAccess}"/>, <see cref="PosInt"/>, a closed form in the tests), so that a loop over the view is specialised for the width and the access mode and the read costs the load, the decode and, in Checked mode, one comparison.
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

    /// <summary>The entry at <paramref name="rank"/> (the checks are the reader's: Checked mode bounds it).</summary>
    internal long this[long rank]
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        get => (long)_reader.Get(rank);
    }

    /// <summary>The same read, always bounded (the contract's).</summary>
    /// <exception cref="InvalidDataException">The rank is outside <c>0..n</c>.</exception>
    internal long At(long rank)
    {
        if ((ulong)rank >= (ulong)_reader.Count) SaViewErrors.Outside(rank, _reader.Count);
        return (long)_reader.Get(rank);
    }

    /// <summary>Copies <c>destination.Length</c> entries that start at <paramref name="start"/>; the range is always checked.</summary>
    /// <exception cref="InvalidDataException">The range is outside the section.</exception>
    internal void CopyTo(long start, Span<long> destination)
    {
        CheckRange(start, destination.Length);
        for (int i = 0; i < destination.Length; i++) destination[i] = (long)_reader.Get(start + i);
    }

    /// <summary>Copies <c>destination.Length</c> entries that start at <paramref name="start"/> as 32-bit positions (the bridge: a position that is not a 32-bit number is damage, not a wrapped number).</summary>
    /// <exception cref="InvalidDataException">The range is outside the section, or an entry does not fit.</exception>
    internal void CopyTo(long start, Span<int> destination)
    {
        CheckRange(start, destination.Length);
        for (int i = 0; i < destination.Length; i++) destination[i] = IntBridge.Value((long)_reader.Get(start + i), "a suffix array entry");
    }

    /// <summary>The entries <c>[start, start + count)</c> as a new array.</summary>
    internal long[] ToArray(long start, int count)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(count);
        var result = new long[count];
        CopyTo(start, result);
        return result;
    }

    private void CheckRange(long start, int length)
    {
        if (start < 0 || length > _reader.Count - start) SaViewErrors.Outside(start, _reader.Count);
    }
}
