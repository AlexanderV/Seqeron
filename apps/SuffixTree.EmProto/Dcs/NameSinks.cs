using System.Buffers.Binary;
using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Dcs;

/// <summary>Where the names go once they are in slot order: consecutive windows of slots, each call continuing where the previous one ended.</summary>
internal interface INameSink : IDisposable
{
    /// <summary>Takes the names of the slots <c>[firstSlot, firstSlot + names.Length)</c>.</summary>
    void Accept(long firstSlot, ReadOnlySpan<uint> names);
}

/// <summary>
/// All names distinct: the names are the ranks, and the names in slot order are the rank file (little-endian <c>uint</c> per slot; the unused slots of the last block stay zero). Written through
/// bounded mapped windows. This machine is little endian, as the in-RAM table and the E3 files assume.
/// </summary>
internal sealed class RankFileSink : INameSink
{
    private readonly MappedWindowWriter _writer;

    internal RankFileSink(string path, long slots, int windowBytes, bool scratch)
    {
        _writer = new MappedWindowWriter(path, slots * 4, windowBytes);
        if (scratch) ScratchTracker.Add(slots * 4);
    }

    public void Accept(long firstSlot, ReadOnlySpan<uint> names) => _writer.WriteAt(firstSlot * 4, MemoryMarshal.AsBytes(names));

    public void Dispose() => _writer.Dispose();
}

/// <summary>
/// The names go to the reduced string as a file: class-major (class c occupies the indices <c>RedStart[c] .. RedStart[c] + ClassCount[c] - 1</c>, a zero separator follows it), every symbol a 32-bit
/// <em>big-endian</em> integer so that the file is a <see cref="PackedText"/> of 32 bits per symbol, followed by the padding of 64 zero bytes. The names arrive in slot order, i.e. block by block;
/// every class gets one write cursor with its own buffer (the sigma-way distribution of ADR-004 section 3.3, here with |D| streams), so every class is written sequentially.
/// </summary>
internal sealed class ReducedFileSink : INameSink
{
    private readonly SampleGeometry _g;
    private readonly CountedFile _file;
    private readonly int _bufEntries;
    private readonly Buf<byte> _buffers;
    private readonly int[] _filled;
    private readonly long[] _flushed;
    internal ReducedFileSink(string path, SampleGeometry geometry, int ioBufferBytes, BudgetedArena arena)
    {
        _g = geometry;
        _bufEntries = Math.Max(1, ioBufferBytes / 4);
        _buffers = arena.Allocate<byte>(_bufEntries * 4 * geometry.Classes, zero: false, "reduced.cursors");
        _filled = new int[geometry.Classes];
        _flushed = new long[geometry.Classes];
        _file = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite);
        long bytes = (geometry.RedLength * 4) + PackedText.PadBytes;
        _file.SetLength(bytes);
        ScratchTracker.Add(bytes);
    }

    public void Accept(long firstSlot, ReadOnlySpan<uint> names)
    {
        int m = _g.Classes;
        int c = (int)(firstSlot % m);
        Span<byte> all = _buffers.AsSpan(_bufEntries * 4 * m);
        for (int i = 0; i < names.Length; i++)
        {
            Span<byte> buf = all.Slice(c * _bufEntries * 4, _bufEntries * 4);
            BinaryPrimitives.WriteUInt32BigEndian(buf[(_filled[c] * 4)..], names[i]);
            if (++_filled[c] == _bufEntries) Flush(c, buf);
            if (++c == m) c = 0;
        }
    }

    /// <summary>Writes the cursors that still hold names.</summary>
    internal void Complete()
    {
        int m = _g.Classes;
        Span<byte> all = _buffers.AsSpan(_bufEntries * 4 * m);
        for (int c = 0; c < m; c++)
        {
            if (_filled[c] > 0) Flush(c, all.Slice(c * _bufEntries * 4, _bufEntries * 4));
        }
    }

    public void Dispose() => _file.Dispose();

    private void Flush(int c, Span<byte> buf)
    {
        _file.WriteAt(buf[..(_filled[c] * 4)], (_g.RedStart[c] + _flushed[c]) * 4);
        _flushed[c] += _filled[c];
        _filled[c] = 0;
    }
}

/// <summary>The names go to the reduced string in RAM (the base case): a zeroed <c>int</c> array, class-major, separators 0.</summary>
internal sealed class ReducedRamSink : INameSink
{
    private readonly SampleGeometry _g;

    internal ReducedRamSink(SampleGeometry geometry, BudgetedArena arena)
    {
        _g = geometry;
        Red = arena.Allocate<int>(checked((int)geometry.RedLength), zero: true, "reduced.string");
    }

    /// <summary>The reduced string.</summary>
    internal Buf<int> Red { get; }

    public void Accept(long firstSlot, ReadOnlySpan<uint> names)
    {
        int m = _g.Classes;
        long block = firstSlot / m;
        int c = (int)(firstSlot % m);
        Span<int> red = Red.AsSpan(checked((int)_g.RedLength));
        for (int i = 0; i < names.Length; i++)
        {
            red[(int)(_g.RedStart[c] + block)] = (int)names[i];
            if (++c == m)
            {
                c = 0;
                block++;
            }
        }
    }

    public void Dispose()
    {
    }
}
