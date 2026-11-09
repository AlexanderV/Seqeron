using System.Buffers.Binary;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

namespace SuffixTree.EmProto.Text;

/// <summary>
/// A text of dense symbols stored with <see cref="Bits"/> bits per symbol (ADR-004 section 6a, Q3: the code is parameterised over b; 2-bit DNA, 3- and 4-bit with N and IUPAC,
/// 8-bit bytes and, for the reduced strings of the difference-cover recursion, up to 32-bit names).
/// <para>
/// Layout: one big-endian bit stream, symbol <c>i</c> occupies stream bits <c>[i*b, (i+1)*b)</c>, the first bit of a byte is its most significant bit; <see cref="PadBytes"/> zero bytes follow
/// the data, so a 64-bit window may start at any symbol up to the end and beyond (the padding reads as the symbol 0, which is why every user clamps comparisons to the length of the
/// shorter suffix). The numeric order of <see cref="Key"/> is then the lexicographic order of the first <see cref="KeySymbols"/> symbols (the key of a suffix shorter than that is padded
/// with zero symbols, hence <c>suffix(i) &lt; suffix(j)</c> implies <c>Key(i) &lt;= Key(j)</c>: a monotone, not strict, key; ties go to the comparator).
/// </para>
/// </summary>
internal sealed unsafe class PackedText : IDisposable
{
    /// <summary>Zero bytes after the last data byte (a 256-bit window plus the carry byte reads at most 41 bytes past the start byte).</summary>
    internal const int PadBytes = 64;

    private static readonly Vector256<byte> ReverseInLane = Vector256.Create(
        (byte)7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8,
        7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);

    private readonly byte* _p;
    private readonly IMemoryBlock? _owner;

    /// <summary>Wraps an existing packed block (owned by <paramref name="owner"/>, which is disposed with this object when <paramref name="ownsOwner"/>).</summary>
    internal PackedText(byte* data, long length, int bits, IMemoryBlock? owner, bool ownsOwner)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(bits, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(bits, 32);
        _p = data;
        Length = length;
        Bits = bits;
        _owner = ownsOwner ? owner : null;
        KeySymbols = 64 / bits;
        KeyShift = 64 - (KeySymbols * bits);
        VectorSymbols = 256 / bits;
    }

    /// <summary>Number of symbols n.</summary>
    internal long Length { get; }

    /// <summary>Bits per symbol b (1..32).</summary>
    internal int Bits { get; }

    /// <summary>Symbols in a 64-bit key: <c>64 / b</c>.</summary>
    internal int KeySymbols { get; }

    internal int KeyShift { get; }

    /// <summary>Symbols compared per step of the vector LCE: <c>256 / b</c>.</summary>
    internal int VectorSymbols { get; }

    /// <summary>When false (tests), the 64-bit scalar LCE is used even where AVX2 exists.</summary>
    internal bool UseVectors { get; set; } = Avx2.IsSupported;

    internal byte* Pointer => _p;

    /// <summary>Bytes of the data part (without padding) for <paramref name="n"/> symbols of <paramref name="bits"/> bits.</summary>
    internal static long DataBytes(long n, int bits) => ((n * bits) + 7) / 8;

    /// <summary>Total bytes of a packed block including the padding.</summary>
    internal static long TotalBytes(long n, int bits) => DataBytes(n, bits) + PadBytes;

    /// <summary>Smallest b with <c>2^b &gt;= sigma</c> (at least 1).</summary>
    internal static int BitsFor(long sigma)
    {
        int b = 1;
        while ((1L << b) < sigma) b++;
        return b;
    }

    /// <summary>Packs dense symbols (each below <c>2^bits</c>) into freshly allocated native memory.</summary>
    internal static PackedText Pack(ReadOnlySpan<byte> symbols, int bits)
    {
        var block = new NativeBlock(TotalBytes(symbols.Length, bits));
        PackInto(symbols, bits, block.Pointer);
        return new PackedText(block.Pointer, symbols.Length, bits, block, ownsOwner: true);
    }

    /// <summary>Packs 32-bit symbols (names of the difference-cover recursion).</summary>
    internal static PackedText Pack(ReadOnlySpan<int> symbols, int bits)
    {
        var block = new NativeBlock(TotalBytes(symbols.Length, bits));
        byte* dst = block.Pointer;
        ulong acc = 0;
        int accBits = 0;
        long o = 0;
        foreach (int s in symbols)
        {
            if (((ulong)(uint)s >> bits) != 0) throw new ArgumentOutOfRangeException(nameof(symbols), $"symbol {s} does not fit in {bits} bits");
            acc = (acc << bits) | (uint)s;
            accBits += bits;
            while (accBits >= 8)
            {
                dst[o++] = (byte)(acc >> (accBits - 8));
                accBits -= 8;
            }
        }
        if (accBits > 0) dst[o] = (byte)(acc << (8 - accBits));
        return new PackedText(block.Pointer, symbols.Length, bits, block, ownsOwner: true);
    }

    /// <summary>Writes the packed form of <paramref name="symbols"/> to <paramref name="dst"/>, which must hold <see cref="TotalBytes"/> zeroed bytes.</summary>
    internal static void PackInto(ReadOnlySpan<byte> symbols, int bits, byte* dst)
    {
        ulong acc = 0;
        int accBits = 0;
        long o = 0;
        foreach (byte s in symbols)
        {
            acc = (acc << bits) | s;
            accBits += bits;
            while (accBits >= 8)
            {
                dst[o++] = (byte)(acc >> (accBits - 8));
                accBits -= 8;
            }
        }
        if (accBits > 0) dst[o] = (byte)(acc << (8 - accBits));
    }

    /// <summary>The 64 stream bits that start at symbol <paramref name="i"/>; the first symbol is in the top bits.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal ulong Window(long i)
    {
        long bit = i * Bits;
        byte* q = _p + (bit >> 3);
        int s = (int)(bit & 7);
        ulong w = BinaryPrimitives.ReverseEndianness(Unsafe.ReadUnaligned<ulong>(q));
        return s == 0 ? w : (w << s) | ((ulong)q[8] >> (8 - s));
    }

    /// <summary>The first <see cref="KeySymbols"/> symbols of the suffix at <paramref name="i"/> as one number (monotone in the suffix order; see the class remarks).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal ulong Key(long i) => Window(i) >> KeyShift;

    /// <summary>The symbol at <paramref name="i"/> (0 for <c>i &gt;= n</c> inside the padding).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint Symbol(long i) => (uint)(Window(i) >> (64 - Bits));

    /// <summary>
    /// Length of the common prefix of the suffixes at <paramref name="i"/> and <paramref name="j"/>, at most <paramref name="max"/>. The caller guarantees
    /// <c>max &lt;= n - max(i, j)</c> (the padding would otherwise read as symbol 0). Vector path: 256 bits (<c>256 / b</c> symbols) per step; scalar path: <c>64 / b</c> symbols.
    /// </summary>
    internal long Lce(long i, long j, long max, ref long steps)
    {
        long k = 0;
        if (UseVectors && max >= VectorSymbols)
        {
            int vs = VectorSymbols;
            int b = Bits;
            while (k + vs <= max)
            {
                steps++;
                Vector256<ulong> a = LoadShifted((i + k) * b);
                Vector256<ulong> c = LoadShifted((j + k) * b);
                uint mask = Vector256.Equals(a.AsByte(), c.AsByte()).ExtractMostSignificantBits();
                if (mask == uint.MaxValue)
                {
                    k += vs;
                    continue;
                }
                int lane = BitOperations.TrailingZeroCount(~mask) >> 3;
                ulong x = a.GetElement(lane) ^ c.GetElement(lane);
                int matched = ((lane * 64) + BitOperations.LeadingZeroCount(x)) / b;
                if (matched < vs) return Math.Min(k + matched, max);
                k += vs; // the difference lies in the partial symbol at the end of the window: continue after the whole symbols
            }
        }
        int ks = KeySymbols;
        int bits = Bits;
        while (k < max)
        {
            steps++;
            ulong x = Window(i + k) ^ Window(j + k);
            if (x != 0)
            {
                int matched = BitOperations.LeadingZeroCount(x) / bits;
                if (matched < ks) return Math.Min(k + matched, max);
            }
            k += ks;
        }
        return max;
    }

    /// <summary>Lexicographic comparison of the suffixes at <paramref name="i"/> and <paramref name="j"/> (the shorter of two equal-prefix suffixes is smaller; the empty suffix at n is the smallest).</summary>
    internal int CompareSuffixes(long i, long j, long skip, ref long steps, out long lce)
    {
        long maxLen = Length - Math.Max(i, j);
        skip = Math.Min(skip, maxLen);
        long m = skip + Lce(i + skip, j + skip, maxLen - skip, ref steps);
        lce = m;
        if (m == maxLen) return i > j ? -1 : 1; // one suffix ended: the one that starts later is shorter, hence smaller
        return Symbol(i + m) < Symbol(j + m) ? -1 : 1;
    }

    public void Dispose() => _owner?.Dispose();

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private Vector256<ulong> LoadShifted(long bit)
    {
        byte* q = _p + (bit >> 3);
        ulong s = (ulong)(bit & 7);
        Vector256<ulong> lo = Avx2.Shuffle(Vector256.Load(q), ReverseInLane).AsUInt64();
        Vector256<ulong> hi = Avx2.Shuffle(Vector256.Load(q + 8), ReverseInLane).AsUInt64();
        return Avx2.ShiftLeftLogical(lo, Vector128.CreateScalar(s)) | Avx2.ShiftRightLogical(hi, Vector128.CreateScalar(64 - s));
    }
}
