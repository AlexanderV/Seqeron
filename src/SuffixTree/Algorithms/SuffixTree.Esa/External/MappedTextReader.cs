using System.Buffers.Binary;
using System.IO.MemoryMappedFiles;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;
using System.Security.Cryptography;

namespace SuffixTree.Esa.External;

/// <summary>
/// Long-indexed random reads of a packed text (TEXT codec 4, <see cref="PackedTextFormat"/>) through a memory-mapped file: the symbol at a position, the 64-bit window and the key (the first <c>64 / b</c> symbols as one
/// number, in the order of the symbols, so the numeric order of keys is the lexicographic order of the first symbols), a bulk read of a range, and the common prefix of two suffixes with an AVX2 path (256 bits per step,
/// guarded by <see cref="Avx2.IsSupported"/>) and a scalar 64-bit path that gives identical results and is the one that runs on arm64. The behaviour is that of the prototype <c>PackedText</c> of E3-E8 (the padding reads
/// as the symbol 0, so every user clamps comparisons to the shorter suffix); what is added is the checking: a position below 0 is an <see cref="ArgumentOutOfRangeException"/>, a position at or beyond the end reads
/// the zero padding (the hot path is one unsigned comparison more than the prototype), a common-prefix limit beyond the text is clamped, and nothing here can read outside the mapping.
/// <para>
/// Opening checks the layout in O(1) (the region has exactly <see cref="PackedTextFormat.TotalBytes"/>, the unused bits of the last data byte and the padding are zero) and reports damage as <see cref="InvalidDataException"/>
/// (values that come from a file are never an <see cref="ArgumentException"/>); <see cref="Verify"/> is the O(n) check of the symbols against the recorded alphabet, counts and SHA-256. Reads are thread-safe and
/// need no lock; a read that races <see cref="Dispose"/> is undefined (as for any mapping), a read after it throws <see cref="ObjectDisposedException"/>.
/// </para>
/// </summary>
internal sealed unsafe class MappedTextReader : IDisposable
{
    private static readonly Vector256<byte> ReverseInLane = Vector256.Create(
        (byte)7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8,
        7, 6, 5, 4, 3, 2, 1, 0, 15, 14, 13, 12, 11, 10, 9, 8);

    /// <summary>Bits per symbol of a reduced string of the difference-cover recursion (see <see cref="OpenReduced"/>).</summary>
    internal const int ReducedBits = 32;

    private readonly byte* _p;
    private readonly long _n;
    private readonly int _bits;
    private readonly int _keySymbols;
    private readonly int _keyShift;
    private readonly int _vectorSymbols;
    private readonly MemoryMappedViewAccessor? _view;
    private readonly MemoryMappedFile? _map;
    private long _limit;
    private bool _useVectors = Avx2.IsSupported;
    private bool _disposed;
#if ESA_BOUNDS
    private readonly long _bytes; // the checked build: the bytes behind _p (data and padding); every raw read below lies inside them
#endif

    private MappedTextReader(byte* data, long n, int bits, MemoryMappedFile? map, MemoryMappedViewAccessor? view)
    {
        _p = data;
        _n = n;
        _limit = n;
        _bits = bits;
        _keySymbols = 64 / bits;
        _keyShift = 64 - (_keySymbols * bits);
        _vectorSymbols = 256 / bits;
        _map = map;
        _view = view;
#if ESA_BOUNDS
        _bytes = PackedTextFormat.TotalBytes(n, bits);
#endif
    }

    /// <summary>Number of symbols n.</summary>
    internal long Length => _n;

    /// <summary>Bits per symbol b (1..8).</summary>
    internal int Bits => _bits;

    /// <summary>The first byte of the packed payload (the pattern descent of <see cref="MappedSearch"/> reads 64-bit windows through it; the 64 bytes of padding make a window at any symbol of the text safe).</summary>
    internal byte* Pointer => _p;

    /// <summary>Symbols in a 64-bit key: <c>64 / b</c>.</summary>
    internal int KeySymbols => _keySymbols;

    /// <summary>Low bits of the 64-bit window that a key drops: <c>64 - (64 / b) * b</c>.</summary>
    internal int KeyShift => _keyShift;

    /// <summary>Symbols compared per step of the vector common prefix: <c>256 / b</c>.</summary>
    internal int VectorSymbols => _vectorSymbols;

    /// <summary>Whether the AVX2 common prefix is used; it cannot be switched on where the hardware has no AVX2 (arm64, old x86), so the scalar path is what runs there. Tests switch it off to compare the two.</summary>
    internal bool UseVectors
    {
        get => _useVectors;
        set => _useVectors = value && Avx2.IsSupported;
    }

    /// <summary>Opens a file that is exactly one packed text (the length of the file is <see cref="PackedTextFormat.TotalBytes"/>).</summary>
    internal static MappedTextReader Open(string path, long n, int bits) => OpenCore(path, 0, null, n, bits);

    /// <summary>Opens the packed text that is a section of a larger file: <paramref name="length"/> (the LENGTH of the section entry) must be <see cref="PackedTextFormat.TotalBytes"/> and the section must lie inside the file.</summary>
    internal static MappedTextReader OpenRegion(string path, long offset, long length, long n, int bits) => OpenCore(path, offset, length, n, bits);

    /// <summary>
    /// The reduced string of the difference-cover recursion (ADR-004 PH1, work package 1b.7d): <paramref name="n"/> symbols of <see cref="ReducedBits"/> = 32 bits each (names and ranks, big endian, the layout of the prototype's
    /// <c>PackedText</c> with 32 bits per symbol) followed by the 64 zero bytes of padding, in a file that is exactly that. Every reading operation is that of the byte-symbol reader; only <see cref="ReadSymbols"/> and the codec
    /// checks of <see cref="Verify"/> are for the 1 to 8 bits of the file codec and refuse such a reader.
    /// </summary>
    internal static MappedTextReader OpenReduced(string path, long n) => OpenCore(path, 0, null, n, ReducedBits, reduced: true);

    /// <summary>A reduced string in memory that the caller owns (tests): the same checks, nothing is mapped or freed.</summary>
    internal static MappedTextReader FromMemoryReduced(byte* data, long length, long n)
    {
        CheckShape(n, ReducedBits, reduced: true);
        CheckLayout(data, length, n, ReducedBits);
        return new MappedTextReader(data, n, ReducedBits, null, null);
    }

    /// <summary>A reader over memory that the caller owns (tests, in-RAM texts): the same checks, nothing is mapped or freed.</summary>
    internal static MappedTextReader FromMemory(byte* data, long length, long n, int bits)
    {
        CheckShape(n, bits);
        CheckLayout(data, length, n, bits);
        return new MappedTextReader(data, n, bits, null, null);
    }

    /// <summary>The 64 stream bits that start at symbol <paramref name="i"/>, the first symbol in the top bits; zero for <c>i &gt;= n</c>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal ulong Window(long i) => (ulong)i < (ulong)_limit ? WindowAt(i) : WindowOutside(i);

    /// <summary>The first <see cref="KeySymbols"/> symbols of the suffix at <paramref name="i"/> as one number (monotone in the suffix order, not strict: the suffixes that end within the key are padded with zero symbols, ties go to the comparator).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal ulong Key(long i) => Window(i) >> _keyShift;

    /// <summary>Asks the CPU to load the cache line of the window at <paramref name="i"/> (a hint for a read that comes a few iterations later; no effect for a position outside the text or where the hardware has no prefetch, never a fault).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal void Prefetch(long i)
    {
        if (Sse.IsSupported && (ulong)i < (ulong)_limit) Sse.Prefetch0(_p + ((i * _bits) >> 3));
    }

    /// <summary>The symbol at <paramref name="i"/>; 0 for <c>i &gt;= n</c>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint Symbol(long i) => (uint)(Window(i) >> (64 - _bits));

    /// <summary>Unpacks the symbols <c>[start, start + destination.Length)</c> into <paramref name="destination"/>, one byte each (the dense symbol).</summary>
    internal void ReadSymbols(long start, Span<byte> destination)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        if (_bits > PackedTextFormat.MaxBits) throw new NotSupportedException($"a reduced string of {_bits} bits per symbol has no one-byte symbols to unpack");
        if (start < 0 || start > _n || destination.Length > _n - start) throw new ArgumentOutOfRangeException(nameof(start), $"range [{start}, {start + destination.Length}) outside [0, {_n})");
#if ESA_BOUNDS
        CheckRead(start * _bits >> 3, (((long)destination.Length * _bits) + 7 >> 3) + 2);
#endif
        if (_bits == 8)
        {
            new ReadOnlySpan<byte>(_p + start, destination.Length).CopyTo(destination);
            return;
        }
        int bits = _bits;
        int mask = (1 << bits) - 1;
        long bit = start * bits;
        for (int k = 0; k < destination.Length; k++, bit += bits)
        {
            byte* q = _p + (bit >> 3);
            int shift = 16 - (int)(bit & 7) - bits;
            destination[k] = (byte)((((q[0] << 8) | q[1]) >> shift) & mask);
        }
    }

    /// <summary>
    /// Length of the common prefix of the suffixes at <paramref name="i"/> and <paramref name="j"/>, at most <paramref name="max"/> (clamped to the shorter suffix, so the padding is never compared); a position at or beyond n
    /// is the empty suffix. Vector path: 256 bits (<c>256 / b</c> symbols) per step; scalar path: <c>64 / b</c> symbols. <paramref name="steps"/> counts the steps (the cost model of the comparator).
    /// </summary>
    internal long Lce(long i, long j, long max, ref long steps)
    {
        if ((i | j | max) < 0) ThrowNegative(i, j, max);
        long room = _limit - Math.Max(i, j); // _limit is 0 after disposal
        if (room <= 0)
        {
            ObjectDisposedException.ThrowIf(_disposed, this);
            return 0;
        }
        if (max > room) max = room;
        long k = 0;
        if (_useVectors && max >= _vectorSymbols)
        {
            int vs = _vectorSymbols;
            int b = _bits;
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
        int ks = _keySymbols;
        int bits = _bits;
        while (k < max)
        {
            steps++;
            ulong x = WindowAt(i + k) ^ WindowAt(j + k);
            if (x != 0)
            {
                int matched = BitOperations.LeadingZeroCount(x) / bits;
                if (matched < ks) return Math.Min(k + matched, max);
            }
            k += ks;
        }
        return max;
    }

    /// <summary>
    /// Lexicographic comparison of the suffixes at <paramref name="i"/> and <paramref name="j"/> (the shorter of two suffixes with equal prefix is the smaller; the empty suffix at n, and any position beyond it, is the smallest),
    /// the first <paramref name="skip"/> symbols being known equal (clamped to the shorter suffix). Returns -1, 0 (the same position) or 1; <paramref name="lce"/> is the common prefix length including the skipped symbols.
    /// </summary>
    internal int CompareSuffixes(long i, long j, long skip, ref long steps, out long lce)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        ArgumentOutOfRangeException.ThrowIfNegative(i);
        ArgumentOutOfRangeException.ThrowIfNegative(j);
        ArgumentOutOfRangeException.ThrowIfNegative(skip);
        i = Math.Min(i, _n);
        j = Math.Min(j, _n);
        long maxLen = _n - Math.Max(i, j);
        if (i == j)
        {
            lce = _n - i;
            return 0;
        }
        skip = Math.Min(skip, maxLen);
        long m = skip + Lce(i + skip, j + skip, maxLen - skip, ref steps);
        lce = m;
        if (m == maxLen) return i > j ? -1 : 1; // one suffix ended: the one that starts later is shorter, hence smaller
        return Symbol(i + m) < Symbol(j + m) ? -1 : 1;
    }

    /// <summary>
    /// Reads every symbol once and checks it against what was recorded about the text: the length and the width, the alphabet size, every symbol below sigma, the symbol counts and (with <paramref name="checkHash"/>) the SHA-256 of
    /// the text that the alphabet maps the symbols back to. O(n) with sequential reads. Any difference is an <see cref="InvalidDataException"/> that names it.
    /// </summary>
    internal void Verify(PackedTextInfo info, bool checkHash = true)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        ArgumentNullException.ThrowIfNull(info);
        if (info.N != _n) throw new InvalidDataException($"the text has {_n} symbols but {info.N} were recorded: wrong length");
        if (info.Bits != _bits) throw new InvalidDataException($"the text is read with {_bits} bits per symbol but {info.Bits} bits were recorded");
        int sigma = info.Letters.Length;
        if (info.Counts.Length != sigma) throw new InvalidDataException($"the alphabet has {sigma} letters but {info.Counts.Length} symbol counts were recorded");
        if (sigma > (1 << _bits)) throw new InvalidDataException($"the alphabet has {sigma} letters, more than the 2^{_bits} symbols of {_bits} bits");
        var seen = new long[sigma];
        var block = new byte[1 << 16];
        using IncrementalHash? hash = checkHash ? IncrementalHash.CreateHash(HashAlgorithmName.SHA256) : null;
        for (long at = 0; at < _n; at += block.Length)
        {
            int take = (int)Math.Min(block.Length, _n - at);
            ReadSymbols(at, block.AsSpan(0, take));
            for (int k = 0; k < take; k++)
            {
                byte s = block[k];
                if (s >= sigma) throw new InvalidDataException($"symbol {s} at position {at + k} is not below sigma = {sigma}");
                seen[s]++;
                block[k] = info.Letters[s];
            }
            hash?.AppendData(block.AsSpan(0, take));
        }
        for (int s = 0; s < sigma; s++)
        {
            if (seen[s] != info.Counts[s]) throw new InvalidDataException($"symbol count of letter {s} is {seen[s]} but {info.Counts[s]} was recorded");
        }
        if (hash is not null && !string.Equals(Convert.ToHexStringLower(hash.GetHashAndReset()), info.Sha256, StringComparison.OrdinalIgnoreCase)) throw new InvalidDataException("the SHA-256 of the text differs from the recorded one");
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        _limit = 0;
        if (_view is not null)
        {
            _view.SafeMemoryMappedViewHandle.ReleasePointer();
            _view.Dispose();
        }
        _map?.Dispose();
    }

    private static MappedTextReader OpenCore(string path, long offset, long? sectionLength, long n, int bits, bool reduced = false)
    {
        ArgumentNullException.ThrowIfNull(path);
        CheckShape(n, bits, reduced);
        long total = PackedTextFormat.TotalBytes(n, bits);
        if (offset < 0) throw new InvalidDataException($"the section starts at offset {offset}: negative");
        if (sectionLength is long declared && declared != total) throw new InvalidDataException($"the section has length {declared} but {n} symbols of {bits} bits need {total} bytes ({PackedTextFormat.DataBytes(n, bits)} data bytes and {PackedTextFormat.PadBytes} of padding)");
        using var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite | FileShare.Delete);
        long fileLength = stream.Length;
        if (offset > fileLength - total) throw new InvalidDataException($"truncated: the file has {fileLength} bytes but the text needs {total} bytes at offset {offset}");
        if (sectionLength is null && fileLength != total) throw new InvalidDataException($"the file has length {fileLength} but {n} symbols of {bits} bits need {total} bytes");
        MemoryMappedFile map = MemoryMappedFile.CreateFromFile(stream, null, 0, MemoryMappedFileAccess.Read, HandleInheritability.None, leaveOpen: true);
        MemoryMappedViewAccessor? view = null;
        try
        {
            view = map.CreateViewAccessor(offset, total, MemoryMappedFileAccess.Read);
            byte* p = null;
            view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
            byte* data = p + view.PointerOffset;
            try
            {
                CheckLayout(data, total, n, bits);
            }
            catch
            {
                view.SafeMemoryMappedViewHandle.ReleasePointer();
                throw;
            }
            return new MappedTextReader(data, n, bits, map, view);
        }
        catch
        {
            view?.Dispose();
            map.Dispose();
            throw;
        }
    }

    private static void CheckShape(long n, int bits, bool reduced = false)
    {
        if (reduced ? bits != ReducedBits : bits < PackedTextFormat.MinBits || bits > PackedTextFormat.MaxBits) throw new InvalidDataException($"{bits} bits per symbol: the packed text has {PackedTextFormat.MinBits} to {PackedTextFormat.MaxBits}");
        if (n < 0 || n > long.MaxValue / PackedTextFormat.MaxBits - PackedTextFormat.PadBytes) throw new InvalidDataException($"{n} symbols: not a length of a text");
    }

    private static void CheckLayout(byte* data, long length, long n, int bits)
    {
        long total = PackedTextFormat.TotalBytes(n, bits);
        long dataBytes = PackedTextFormat.DataBytes(n, bits);
        if (length < total) throw new InvalidDataException($"truncated: the packed text has {length} bytes but {n} symbols of {bits} bits need {total}");
        if (length > total) throw new InvalidDataException($"the packed text has length {length} but {n} symbols of {bits} bits need {total}");
        int unused = (int)((dataBytes * 8) - (n * bits));
        if (unused > 0 && (data[dataBytes - 1] & ((1 << unused) - 1)) != 0) throw new InvalidDataException($"the {unused} unused bits of the last data byte are not zero");
        int bad = new ReadOnlySpan<byte>(data + dataBytes, PackedTextFormat.PadBytes).IndexOfAnyExcept((byte)0);
        if (bad >= 0) throw new InvalidDataException($"the padding byte {bad} after the data is not zero");
    }

    /// <summary>The window of a position inside the text (the caller checked <c>0 &lt;= i &lt; n</c>): the eight bytes from the byte of the symbol, shifted by the bit offset, and the carry from the ninth.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private ulong WindowAt(long i)
    {
        long bit = i * _bits;
#if ESA_BOUNDS
        CheckRead(bit >> 3, 9);
#endif
        byte* q = _p + (bit >> 3);
        int s = (int)(bit & 7);
        ulong w = BinaryPrimitives.ReverseEndianness(Unsafe.ReadUnaligned<ulong>(q));
        return s == 0 ? w : (w << s) | ((ulong)q[8] >> (8 - s));
    }

#if ESA_BOUNDS
    /// <summary>The checked build: a raw read of <paramref name="count"/> bytes from byte <paramref name="offset"/> of the mapped text must lie inside the region (data and the 64 bytes of padding).</summary>
    private void CheckRead(long offset, long count)
    {
        if (offset < 0 || count > _bytes - offset) ThrowOutOfBounds(offset, count);
    }

    private void ThrowOutOfBounds(long offset, long count)
        => throw new EsaMemoryViolationException(
            MemoryViolationKind.OutOfBounds,
            $"out-of-bounds read of the mapped text: {count} bytes from byte {offset} of a region of {_bytes} bytes ({_n} symbols of {_bits} bits)");
#endif

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static void ThrowNegative(long i, long j, long max)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(i);
        ArgumentOutOfRangeException.ThrowIfNegative(j);
        ArgumentOutOfRangeException.ThrowIfNegative(max);
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private ulong WindowOutside(long i)
    {
        ObjectDisposedException.ThrowIf(_disposed, this);
        ArgumentOutOfRangeException.ThrowIfNegative(i);
        return 0;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private Vector256<ulong> LoadShifted(long bit)
    {
#if ESA_BOUNDS
        CheckRead(bit >> 3, 40);
#endif
        byte* q = _p + (bit >> 3);
        ulong s = (ulong)(bit & 7);
        Vector256<ulong> lo = Avx2.Shuffle(Vector256.Load(q), ReverseInLane).AsUInt64();
        Vector256<ulong> hi = Avx2.Shuffle(Vector256.Load(q + 8), ReverseInLane).AsUInt64();
        return Avx2.ShiftLeftLogical(lo, Vector128.CreateScalar(s)) | Avx2.ShiftRightLogical(hi, Vector128.CreateScalar(64 - s));
    }
}
