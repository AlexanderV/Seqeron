using System.Buffers.Binary;

namespace SuffixTree.Esa.External;

/// <summary>
/// The writer of the RMQ section (kind 9: the persisted LCP pyramid of the navigation pack; ADR-005 sections 3.3 and 6.1, work package S2.5): the block minima of the LCP section as it lies in the container, level by level. Level 0 holds
/// the minimum of every <see cref="SesaFormat.RmqFanout"/> consecutive LCP entries (the last block may be shorter), level k + 1 the minima of the blocks of level k, until a level has one block; an LCP array of a single entry (the empty
/// text) has no level and gets no section. The levels are written one after the other, level 0 first, as 4-byte little-endian numbers (<c>POS_WIDTH</c> 4; the width-5 form is the one of ADR-005 section 6.1 and belongs to the 40-bit
/// builder of S2.7). The LCP section is read once, sequentially, in blocks of 1 MiB with its overflow list merged by rank (both are in rank order), so the memory is the section itself, 0.13 bytes per character, and nothing is built twice.
/// </summary>
internal static class RmqWriter
{
    private const int ReadBlockEntries = 1 << 20;

    /// <summary>The number of entries of all levels of the pyramid of an LCP array of <paramref name="rows"/> entries (0 for a single entry).</summary>
    internal static long CountFor(long rows, int fanout = SesaFormat.RmqFanout)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(rows, 1);
        CheckFanout(fanout);
        long total = 0;
        for (long size = rows; size > 1;)
        {
            size = (size + fanout - 1) / fanout;
            total += size;
        }
        return total;
    }

    /// <summary>The payload bytes of the section (4 per entry at <c>POS_WIDTH</c> 4).</summary>
    internal static long LengthFor(long rows, int fanout = SesaFormat.RmqFanout) => 4 * CountFor(rows, fanout);

    /// <summary>
    /// The levels of the pyramid of the LCP section that lies at <paramref name="lcpOffset"/> of the file <paramref name="path"/>: <paramref name="rows"/> entries of the width of <paramref name="kind"/> and, for an 8- or 16-bit array
    /// with escapes, the <paramref name="overflowCount"/> pairs (rank, value) of 4-byte numbers at <paramref name="overflowOffset"/>.
    /// </summary>
    /// <exception cref="InvalidDataException">The overflow list is not in ascending rank order, or does not match the escapes of the base array: the container that was just written is inconsistent.</exception>
    internal static uint[] Build(string path, LcpKind kind, long lcpOffset, long rows, long overflowOffset, long overflowCount, IoCounters? counters, int fanout = SesaFormat.RmqFanout)
    {
        ArgumentException.ThrowIfNullOrEmpty(path);
        CheckFanout(fanout);
        if (kind == LcpKind.U64) throw new NotSupportedException("the RMQ section of a text with a 64-bit LCP array is written by the 40-bit builder of S2.7");
        long count = CountFor(rows, fanout);
        var levels = new uint[checked((int)count)];
        if (count == 0) return levels;

        int width = LcpWidth.Bytes(kind);
        uint escape = kind switch { LcpKind.U8 => byte.MaxValue, LcpKind.U16 => ushort.MaxValue, _ => uint.MaxValue };
        long blocks = (rows + fanout - 1) / fanout;
        using var file = new CountedFile(path, FileMode.Open, FileAccess.Read, counters);
        int entries = (int)Math.Min(rows, (long)ReadBlockEntries / fanout * fanout);
        var buffer = new byte[checked(entries * width)];
        var pairs = new OverflowCursor(file, overflowOffset, kind == LcpKind.U32 ? 0 : overflowCount);

        long block = 0;
        uint min = uint.MaxValue;
        for (long first = 0; first < rows; first += entries)
        {
            int take = (int)Math.Min(entries, rows - first);
            file.ReadExactAt(buffer.AsSpan(0, take * width), lcpOffset + (first * width));
            for (int i = 0; i < take; i++)
            {
                uint v = width switch
                {
                    1 => buffer[i],
                    2 => BinaryPrimitives.ReadUInt16LittleEndian(buffer.AsSpan(2 * i)),
                    _ => BinaryPrimitives.ReadUInt32LittleEndian(buffer.AsSpan(4 * i)),
                };
                if (kind != LcpKind.U32 && v == escape) v = pairs.ValueAt(first + i);
                if (v < min) min = v;
                if ((((first + i + 1) % fanout) == 0) || first + i + 1 == rows)
                {
                    levels[block++] = min;
                    min = uint.MaxValue;
                }
            }
        }
        if (block != blocks || pairs.Remaining != 0) throw new InvalidDataException("the overflow list of the LCP section does not match the escapes of its base array");

        // the upper levels from the level below: the same minimum over groups of the fanout
        long below = 0, belowCount = blocks, here = blocks;
        while (belowCount > 1)
        {
            long hereCount = (belowCount + fanout - 1) / fanout;
            for (long b = 0; b < hereCount; b++)
            {
                uint m = uint.MaxValue;
                for (long j = b * fanout, end = Math.Min(belowCount, (b + 1) * fanout); j < end; j++) m = Math.Min(m, levels[below + j]);
                levels[here + b] = m;
            }
            below = here;
            here += hereCount;
            belowCount = hereCount;
        }
        return levels;
    }

    private static void CheckFanout(int fanout)
    {
        if (fanout < 2 || (fanout & (fanout - 1)) != 0) throw new ArgumentOutOfRangeException(nameof(fanout), fanout, "the fanout must be a power of two of at least 2");
    }

    /// <summary>The overflow pairs <c>(rank, value)</c> of the LCP section in rank order, read in blocks; <see cref="ValueAt"/> must be asked for ascending ranks.</summary>
    private sealed class OverflowCursor
    {
        private const int PairsPerBlock = 1 << 16;
        private readonly CountedFile _file;
        private readonly long _offset;
        private readonly long _count;
        private readonly byte[] _block = new byte[8 * PairsPerBlock];
        private long _loaded;     // pairs read so far
        private int _have;        // pairs in the block
        private int _at;          // next pair in the block

        internal OverflowCursor(CountedFile file, long offset, long count)
        {
            _file = file;
            _offset = offset;
            _count = count;
        }

        /// <summary>Pairs not yet consumed.</summary>
        internal long Remaining => _count - (_loaded - _have + _at);

        internal uint ValueAt(long rank)
        {
            if (_at == _have) Refill();
            if (_at == _have) throw new InvalidDataException($"the overflow list ends before the escape at rank {rank}");
            uint r = BinaryPrimitives.ReadUInt32LittleEndian(_block.AsSpan(8 * _at));
            if (r != rank) throw new InvalidDataException($"the escape at rank {rank} has no pair in the overflow list (the next pair is for rank {r})");
            uint v = BinaryPrimitives.ReadUInt32LittleEndian(_block.AsSpan((8 * _at) + 4));
            _at++;
            return v;
        }

        private void Refill()
        {
            int take = (int)Math.Min(PairsPerBlock, _count - _loaded);
            if (take <= 0) return;
            _file.ReadExactAt(_block.AsSpan(0, 8 * take), _offset + (8 * _loaded));
            _loaded += take;
            _have = take;
            _at = 0;
        }
    }
}
