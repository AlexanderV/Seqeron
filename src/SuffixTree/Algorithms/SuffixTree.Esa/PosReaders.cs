using System.Runtime.CompilerServices;

namespace SuffixTree.Esa;

/// <summary>
/// The stored width of an array of positions, ranks or LCP values as a struct strategy (ADR-005 section 3.3, work package S2.3): one element at a <c>long</c> index, returned as an unsigned 64-bit number in a register. A reader is a
/// <c>readonly struct</c> that holds a pointer (or an array) and the element count of its section; a generic view over it is specialised by the JIT for each reader (as <see cref="MappedSearch"/> is for the access mode), so a read costs the
/// load and the decode and nothing else. The readers of a container are <see cref="Pos32{TAccess}"/> (<c>POS_WIDTH</c> 4), <see cref="Pos40{TAccess}"/> (<c>POS_WIDTH</c> 5) and, for the LCP array, <see cref="Pos8{TAccess}"/>,
/// <see cref="Pos16{TAccess}"/> and <see cref="Pos64{TAccess}"/>; <see cref="PosInt"/> reads an <c>int[]</c> of the in-memory index; the tests add a reader that computes its elements (<c>PosGen</c>).
/// </summary>
internal interface IPosReader
{
    /// <summary>Bytes of one element of the stored form (4 for <see cref="Pos32{TAccess}"/> and <see cref="PosInt"/>, 5 for <see cref="Pos40{TAccess}"/>); a closed form that computes its elements says the width it stands for.</summary>
    static abstract int ElementBytes { get; }

    /// <summary>The number of elements of the section (a reader never reads the element <c>Count</c> or beyond in <see cref="EsaAccessMode.Checked"/> mode).</summary>
    long Count { get; }

    /// <summary>
    /// The element at <paramref name="index"/>. <see cref="EsaAccessMode.Checked"/> (<see cref="CheckedAccess"/>): an index outside <c>0..Count-1</c> is an <see cref="InvalidDataException"/> (a damaged table sent the read there) and never a
    /// read outside the section; <see cref="EsaAccessMode.Trusted"/>: no check.
    /// </summary>
    ulong Get(long index);
}

/// <summary>The bounds check and the failure that every pointer reader shares.</summary>
internal static class PosReader
{
    /// <summary>The mask of the 40 low bits that a <see cref="Pos40{TAccess}"/> keeps of its unaligned 8-byte load.</summary>
    internal const ulong Mask40 = 0xFF_FFFF_FFFFUL;

    /// <summary>The padding that a 5-byte section (and an <c>LCP_OVERFLOW</c> section of pairs of 5-byte fields) carries after its last element, inside its payload: the load at the last element reads 8 bytes.</summary>
    internal const int Padding40 = 3;

    [MethodImpl(MethodImplOptions.NoInlining)]
    internal static void ThrowOutside(long index, long count) =>
        throw new InvalidDataException($"a section of {count} elements is read at element {index}: the file is damaged");

    internal static void CheckCount(long count) => ArgumentOutOfRangeException.ThrowIfNegative(count);
}

/// <summary>One byte per element (the LCP array of kind U8).</summary>
internal readonly unsafe struct Pos8<TAccess> : IPosReader
    where TAccess : struct, IAccessChecks
{
    private readonly byte* _p;
    private readonly long _count;

    internal Pos8(byte* section, long count)
    {
        PosReader.CheckCount(count);
        _p = section;
        _count = count;
    }

    public static int ElementBytes => 1;

    public long Count => _count;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long index)
    {
        if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);
        return _p[index];
    }
}

/// <summary>Two bytes per element, little endian (the LCP array of kind U16).</summary>
internal readonly unsafe struct Pos16<TAccess> : IPosReader
    where TAccess : struct, IAccessChecks
{
    private readonly byte* _p;
    private readonly long _count;

    internal Pos16(byte* section, long count)
    {
        PosReader.CheckCount(count);
        _p = section;
        _count = count;
    }

    public static int ElementBytes => 2;

    public long Count => _count;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long index)
    {
        if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);
        return Unsafe.ReadUnaligned<ushort>(_p + (2 * index));
    }
}

/// <summary>Four bytes per element, little endian, unsigned (<c>POS_WIDTH</c> 4; the LCP array of kind U32).</summary>
internal readonly unsafe struct Pos32<TAccess> : IPosReader
    where TAccess : struct, IAccessChecks
{
    private readonly byte* _p;
    private readonly long _count;

    internal Pos32(byte* section, long count)
    {
        PosReader.CheckCount(count);
        _p = section;
        _count = count;
    }

    public static int ElementBytes => 4;

    public long Count => _count;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long index)
    {
        if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);
        return Unsafe.ReadUnaligned<uint>(_p + (4 * index));
    }
}

/// <summary>
/// Five bytes per element, little endian (<c>POS_WIDTH</c> 5, ADR-005 section 6.1): one unaligned 8-byte load at <c>5 * index</c> and a mask of the 40 low bits. The section ends with 3 zero bytes that belong to its payload, so that the load
/// at the last element ends inside it; the caller of the constructor guarantees <c>5 * count + 3</c> readable bytes.
/// </summary>
internal readonly unsafe struct Pos40<TAccess> : IPosReader
    where TAccess : struct, IAccessChecks
{
    private readonly byte* _p;
    private readonly long _count;

    internal Pos40(byte* section, long count)
    {
        PosReader.CheckCount(count);
        _p = section;
        _count = count;
    }

    public static int ElementBytes => 5;

    public long Count => _count;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long index)
    {
        if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);
        return Unsafe.ReadUnaligned<ulong>(_p + (5 * index)) & PosReader.Mask40;
    }
}

/// <summary>Eight bytes per element, little endian (the LCP array of kind U64: a text above 2^32 with a repeat longer than 2^32 - 2).</summary>
internal readonly unsafe struct Pos64<TAccess> : IPosReader
    where TAccess : struct, IAccessChecks
{
    private readonly byte* _p;
    private readonly long _count;

    internal Pos64(byte* section, long count)
    {
        PosReader.CheckCount(count);
        _p = section;
        _count = count;
    }

    public static int ElementBytes => 8;

    public long Count => _count;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long index)
    {
        if (TAccess.IsChecked && (ulong)index >= (ulong)_count) PosReader.ThrowOutside(index, _count);
        return Unsafe.ReadUnaligned<ulong>(_p + (8 * index));
    }
}

/// <summary>The <c>int[]</c> of an index built in memory, read as unsigned 32-bit numbers (a position or a rank is never negative); the bounds are the array's own (an <see cref="IndexOutOfRangeException"/>).</summary>
internal readonly struct PosInt : IPosReader
{
    private readonly int[] _values;

    internal PosInt(int[] values)
    {
        ArgumentNullException.ThrowIfNull(values);
        _values = values;
    }

    public static int ElementBytes => 4;

    public long Count => _values.Length;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public ulong Get(long index) => (uint)_values[index];
}
