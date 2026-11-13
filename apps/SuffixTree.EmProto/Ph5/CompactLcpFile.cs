using System.Runtime.InteropServices;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Ph5;

/// <summary>Receives blocks of resolved LCP values in the order a file reader produces them.</summary>
internal delegate void LcpBlockHandler(long firstRank, ReadOnlySpan<uint> values);

/// <summary>
/// The compact LCP of PH4 (ADR-004 section 4.2, B3; the files of <c>e4</c>: a base array of <c>n + 1</c> entries of 1, 2 or 4 bytes, little endian, the escape value 255 / 65,535 in the narrow kinds, and an overflow file of
/// pairs <c>(rank, value)</c>, 4 bytes each, ascending ranks) read in rank order or backwards in blocks, the escapes resolved through a moving cursor over the overflow pairs (no binary search, no random access).
/// This is the "given LCP file" form of the PH5 input and the source of layout F's backward scan. Every read goes through <see cref="CountedFile"/>, so the volume is counted. Not thread-safe.
/// </summary>
internal sealed class CompactLcpFile : IDisposable
{
    private readonly CountedFile _base;
    private readonly CountedFile _overflow;
    private readonly long _n;
    private readonly int _blockEntries;
    private readonly uint _escape;
    private readonly byte[] _bytes;
    private readonly uint[] _values;
    private readonly uint[] _pairs;
    private long _pairsStart = -1;
    private int _pairsCount;

    /// <param name="basePath">The base array.</param>
    /// <param name="overflowPath">The overflow pairs.</param>
    /// <param name="n">Text length; the base array has <c>n + 1</c> entries.</param>
    /// <param name="bufferBytes">Bytes of the base array read per block (and of the overflow pairs per read).</param>
    /// <exception cref="ArgumentException">The file lengths do not fit <c>n</c>.</exception>
    internal CompactLcpFile(string basePath, string overflowPath, long n, int bufferBytes = 1 << 20)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        _n = n;
        long baseLength = new FileInfo(basePath).Length;
        long overflowLength = new FileInfo(overflowPath).Length;
        Width = (int)(baseLength / (n + 1));
        if (Width is not (1 or 2 or 4) || baseLength != (n + 1) * Width) throw new ArgumentException($"the base file has {baseLength} bytes, not (n + 1) x 1, 2 or 4 for n = {n}");
        if (overflowLength % 8 != 0) throw new ArgumentException($"the overflow file has {overflowLength} bytes, not a whole number of pairs");
        OverflowCount = overflowLength / 8;
        if (Width == 4 && OverflowCount != 0) throw new ArgumentException("a base array of 4-byte entries has no overflow list");
        _escape = Width switch { 1 => byte.MaxValue, 2 => ushort.MaxValue, _ => uint.MaxValue };
        _blockEntries = Math.Max(1, bufferBytes / Width);
        _bytes = new byte[_blockEntries * Width];
        _values = new uint[_blockEntries];
        _pairs = new uint[2 * Math.Max(1, bufferBytes / 8)];
        _base = new CountedFile(basePath, FileMode.Open, FileAccess.Read);
        _overflow = new CountedFile(overflowPath, FileMode.Open, FileAccess.Read);
    }

    /// <summary>Bytes per entry of the base array.</summary>
    internal int Width { get; }

    /// <summary>Pairs in the overflow file.</summary>
    internal long OverflowCount { get; }

    /// <summary>Calls <paramref name="handler"/> with consecutive blocks from rank 0 to rank n (resolved values; the span is valid for the call).</summary>
    /// <exception cref="InvalidDataException">The overflow pairs do not match the escapes of the base array.</exception>
    internal void ForEachForward(LcpBlockHandler handler)
    {
        long pair = 0;
        _pairsStart = -1;
        for (long first = 0; first <= _n; first += _blockEntries)
        {
            int count = (int)Math.Min(_blockEntries, _n + 1 - first);
            Load(first, count);
            for (int i = 0; i < count; i++)
            {
                if (_values[i] != _escape || Width == 4) continue;
                if (pair >= OverflowCount) throw new InvalidDataException($"rank {first + i} is an escape but the overflow list has only {OverflowCount} pairs");
                (uint rank, uint value) = Pair(pair++, backward: false);
                if (rank != first + i || value < _escape) throw new InvalidDataException($"overflow pair {pair - 1} is ({rank}, {value}), expected rank {first + i} and a value of at least {_escape}");
                _values[i] = value;
            }
            handler(first, _values.AsSpan(0, count));
        }
        if (pair != OverflowCount) throw new InvalidDataException($"{OverflowCount} overflow pairs, {pair} escapes in the base array");
    }

    /// <summary>Calls <paramref name="handler"/> with consecutive blocks from the end of the array to rank 0 (each block in ascending rank order inside, blocks in descending order).</summary>
    /// <exception cref="InvalidDataException">The overflow pairs do not match the escapes of the base array.</exception>
    internal void ForEachBackward(LcpBlockHandler handler)
    {
        long pair = OverflowCount - 1;
        _pairsStart = -1;
        long total = _n + 1;
        long blocks = (total + _blockEntries - 1) / _blockEntries;
        for (long b = blocks - 1; b >= 0; b--)
        {
            long first = b * _blockEntries;
            int count = (int)Math.Min(_blockEntries, total - first);
            Load(first, count);
            for (int i = count - 1; i >= 0; i--)
            {
                if (_values[i] != _escape || Width == 4) continue;
                if (pair < 0) throw new InvalidDataException($"rank {first + i} is an escape but the overflow list is exhausted");
                (uint rank, uint value) = Pair(pair--, backward: true);
                if (rank != first + i || value < _escape) throw new InvalidDataException($"overflow pair {pair + 1} is ({rank}, {value}), expected rank {first + i} and a value of at least {_escape}");
                _values[i] = value;
            }
            handler(first, _values.AsSpan(0, count));
        }
        if (pair != -1) throw new InvalidDataException($"{OverflowCount} overflow pairs, {OverflowCount - 1 - pair} escapes in the base array");
    }

    public void Dispose()
    {
        _base.Dispose();
        _overflow.Dispose();
    }

    private void Load(long firstRank, int count)
    {
        Span<byte> bytes = _bytes.AsSpan(0, count * Width);
        _base.ReadExactAt(bytes, firstRank * Width);
        switch (Width)
        {
            case 1:
                for (int i = 0; i < count; i++) _values[i] = bytes[i];
                break;
            case 2:
                ReadOnlySpan<ushort> u16 = MemoryMarshal.Cast<byte, ushort>(bytes);
                for (int i = 0; i < count; i++) _values[i] = u16[i];
                break;
            default:
                MemoryMarshal.Cast<byte, uint>(bytes).CopyTo(_values);
                break;
        }
    }

    /// <summary>Pair number <paramref name="index"/>, from a window of the overflow file that is read when the index leaves it (sequential: the window starts at the index going forwards and ends at it going backwards).</summary>
    private (uint Rank, uint Value) Pair(long index, bool backward)
    {
        int capacity = _pairs.Length / 2;
        if (_pairsStart < 0 || index < _pairsStart || index >= _pairsStart + _pairsCount)
        {
            long start = backward ? Math.Max(0, index - capacity + 1) : index;
            int take = (int)Math.Min(capacity, OverflowCount - start);
            _overflow.ReadExactAt(MemoryMarshal.AsBytes(_pairs.AsSpan(0, take * 2)), start * 8);
            _pairsStart = start;
            _pairsCount = take;
        }
        int at = (int)(index - _pairsStart) * 2;
        return (_pairs[at], _pairs[at + 1]);
    }
}
