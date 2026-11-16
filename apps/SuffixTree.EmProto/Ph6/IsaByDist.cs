using System.Runtime.InteropServices;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.Esa;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Ph6;

internal sealed record DistOptions
{
    /// <summary>Directory of the chunk files (the tuples <c>(SA[i], i)</c>); created if missing, empty again afterwards.</summary>
    internal string ScratchDir { get; init; } = "";

    /// <summary>R: positions per chunk (the ADR: 8 M; P = (n + 1) / R chunks).</summary>
    internal long ChunkPositions { get; init; } = 8_000_000;

    /// <summary>Write buffer of one chunk (the ADR: 64-256 KB); a chunk that never outgrows it stays in memory.</summary>
    internal int IoBufferBytes { get; init; } = 64 * 1024;

    /// <summary>Bytes of the suffix array read at a time.</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    internal int WindowBytes { get; init; } = 1 << 20;

    internal ArenaBudget? Budget { get; init; }
}

internal sealed record DistStats(long Chunks, long TupleBytesWritten, long TupleBytesRead);

/// <summary>
/// The inverse permutation of the suffix array by one-level radix distribution (the DIST primitive of ADR-004 section 3.3): <c>ISA[SA[i]] = i</c>, 4 (n + 1) bytes. Pass 1 reads the SA in rank order and distributes the
/// tuples <c>(SA[i], i)</c> (8 bytes) into P = ceil((n + 1) / R) chunk files by <c>SA[i] / R</c> through P write buffers (a chunk that never outgrows its buffer stays in memory); pass 2 takes the chunks in position
/// order, scatters each into an array of R entries and appends it to the ISA file through the bounded mapped window writer. I/O: SA read 4 + tuples written 8 + tuples read 8 + ISA written 4 = 24 B per entry.
/// </summary>
internal static class IsaByDist
{
    internal static DistStats Build(string saPath, long n, string isaPath, DistOptions o)
    {
        ArgumentNullException.ThrowIfNull(o);
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.ChunkPositions, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.IoBufferBytes, 8);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.ReadBytes, 4);
        ArgumentOutOfRangeException.ThrowIfLessThan(o.WindowBytes, 1);
        if (n >= uint.MaxValue - 2) throw new NotSupportedException("positions are 32-bit in the prototype");
        long rows = n + 1;
        long chunkPositions = o.ChunkPositions;
        long chunks = (rows + chunkPositions - 1) / chunkPositions;
        int bufEntries = o.IoBufferBytes / 8;
        int readEntries = (int)Math.Min(o.ReadBytes / 4, rows);
        using var sa = new CountedFile(saPath, FileMode.Open, FileAccess.Read);
        if (sa.Length != 4 * rows) throw new InvalidDataException($"the suffix array file has {sa.Length} bytes, expected {4 * rows} for n = {n}");
        Directory.CreateDirectory(o.ScratchDir);

        var filled = new int[chunks];
        var flushed = new long[chunks];
        long bytesWritten = 0, bytesRead = 0;
        var arenaBudget = o.Budget ?? new ArenaBudget(long.MaxValue);
        try
        {
            using var arena = new BudgetedArena(arenaBudget);
            Buf<ulong> tuples = arena.Allocate<ulong>(checked((int)(chunks * bufEntries)), zero: false, "dist.tuples");
            Buf<uint> readBuf = arena.Allocate<uint>(readEntries, zero: false, "dist.sa");
            Span<ulong> all = tuples.AsSpan(checked((int)(chunks * bufEntries)));

            // pass 1: the suffix array in rank order
            for (long first = 0; first < rows; first += readEntries)
            {
                int take = (int)Math.Min(readEntries, rows - first);
                Span<uint> entries = readBuf.AsSpan(take);
                sa.ReadExactAt(MemoryMarshal.AsBytes(entries), first * 4);
                for (int i = 0; i < take; i++)
                {
                    uint p = entries[i];
                    if (p > n) throw new InvalidDataException($"SA[{first + i}] = {p} lies beyond n = {n}");
                    long chunk = p / chunkPositions;
                    Span<ulong> buf = all.Slice(checked((int)(chunk * bufEntries)), bufEntries);
                    buf[filled[chunk]++] = ((ulong)p << 32) | (uint)(first + i);
                    if (filled[chunk] == bufEntries) Flush(o, chunk, buf, filled, flushed, ref bytesWritten);
                }
            }

            // pass 2: chunk by chunk in position order
            using var writer = new MappedWindowWriter(isaPath, 4 * rows, o.WindowBytes);
            Buf<uint> chunkArray = arena.Allocate<uint>((int)Math.Min(chunkPositions, rows), zero: false, "dist.chunk");
            Buf<ulong> chunkRead = arena.Allocate<ulong>(bufEntries, zero: false, "dist.read");
            for (long chunk = 0; chunk < chunks; chunk++)
            {
                long startPos = chunk * chunkPositions;
                int length = (int)Math.Min(chunkPositions, rows - startPos);
                Span<uint> array = chunkArray.AsSpan(length);
                array.Fill(uint.MaxValue);
                string path = ChunkPath(o, chunk);
                if (flushed[chunk] > 0)
                {
                    using (var f = new CountedFile(path, FileMode.Open, FileAccess.Read))
                    {
                        for (long at = 0; at < flushed[chunk]; at += bufEntries)
                        {
                            int take = (int)Math.Min(bufEntries, flushed[chunk] - at);
                            Span<ulong> part = chunkRead.AsSpan(take);
                            f.ReadExactAt(MemoryMarshal.AsBytes(part), at * 8);
                            bytesRead += take * 8L;
                            Scatter(part, array, startPos);
                        }
                    }
                    ScratchTracker.Remove(flushed[chunk] * 8);
                    File.Delete(path);
                    flushed[chunk] = -flushed[chunk]; // done: the file is gone
                }
                Scatter(all.Slice(checked((int)(chunk * bufEntries)), filled[chunk]), array, startPos);
                for (int k = 0; k < length; k++)
                {
                    if (array[k] == uint.MaxValue) throw new InvalidDataException($"the suffix array is not a permutation: position {startPos + k} does not occur");
                }
                writer.Append(MemoryMarshal.AsBytes(array));
            }
            return new DistStats(chunks, bytesWritten, bytesRead);
        }
        finally
        {
            for (long chunk = 0; chunk < chunks; chunk++)
            {
                if (flushed[chunk] > 0)
                {
                    ScratchTracker.Remove(flushed[chunk] * 8);
                    File.Delete(ChunkPath(o, chunk));
                }
            }
        }
    }

    private static void Scatter(ReadOnlySpan<ulong> tuples, Span<uint> array, long startPos)
    {
        foreach (ulong v in tuples)
        {
            ref uint slot = ref array[(int)((long)(v >> 32) - startPos)];
            if (slot != uint.MaxValue) throw new InvalidDataException($"the suffix array is not a permutation: position {(long)(v >> 32)} occurs twice");
            slot = (uint)v;
        }
    }

    private static string ChunkPath(DistOptions o, long chunk) => Path.Combine(o.ScratchDir, $"dist-tup-{chunk:D6}.bin");

    private static void Flush(DistOptions o, long chunk, ReadOnlySpan<ulong> buf, int[] filled, long[] flushed, ref long bytesWritten)
    {
        using (var f = new CountedFile(ChunkPath(o, chunk), FileMode.OpenOrCreate, FileAccess.Write))
        {
            f.WriteAt(MemoryMarshal.AsBytes(buf[..filled[chunk]]), flushed[chunk] * 8);
        }
        ScratchTracker.Add(filled[chunk] * 8L);
        bytesWritten += filled[chunk] * 8L;
        flushed[chunk] += filled[chunk];
        filled[chunk] = 0;
    }
}

/// <summary>The other route to Psi: <c>PSI[i] = ISA[SA[i] + 1]</c> from the suffix array and the ISA file (PSI[0] = 0): the SA is read sequentially, the ISA at random through its mapping.</summary>
internal static class PsiFromIsa
{
    internal static void Build(string saPath, string isaPath, long n, string psiPath, int readBytes = 1 << 20, int windowBytes = 1 << 20, ArenaBudget? budget = null)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        ArgumentOutOfRangeException.ThrowIfLessThan(readBytes, 4);
        ArgumentOutOfRangeException.ThrowIfLessThan(windowBytes, 1);
        if (n >= uint.MaxValue - 2) throw new NotSupportedException("positions are 32-bit in the prototype");
        long rows = n + 1;
        using var sa = new CountedFile(saPath, FileMode.Open, FileAccess.Read);
        if (sa.Length != 4 * rows) throw new InvalidDataException($"the suffix array file has {sa.Length} bytes, expected {4 * rows} for n = {n}");
        if (new FileInfo(isaPath).Length != 4 * rows) throw new InvalidDataException($"the ISA file has {new FileInfo(isaPath).Length} bytes, expected {4 * rows} for n = {n}");
        using var isa = new MappedFile(isaPath, writable: false);
        int entries = (int)Math.Min(readBytes / 4, rows);
        using var arena = new BudgetedArena(budget ?? new ArenaBudget(long.MaxValue));
        Buf<uint> saBuf = arena.Allocate<uint>(entries, zero: false, "psiisa.sa");
        Buf<uint> outBuf = arena.Allocate<uint>(entries, zero: false, "psiisa.out");
        using var writer = new MappedWindowWriter(psiPath, 4 * rows, windowBytes);
        unsafe
        {
            uint* table = (uint*)isa.Pointer;
            for (long first = 0; first < rows; first += entries)
            {
                int take = (int)Math.Min(entries, rows - first);
                Span<uint> read = saBuf.AsSpan(take);
                sa.ReadExactAt(MemoryMarshal.AsBytes(read), first * 4);
                Span<uint> output = outBuf.AsSpan(take);
                for (int i = 0; i < take; i++)
                {
                    if (first + i == 0)
                    {
                        output[i] = 0;
                        continue;
                    }
                    uint p = read[i];
                    if (p >= n) throw new InvalidDataException($"SA[{first + i}] = {p} is not a position below n = {n}");
                    output[i] = table[p + 1];
                }
                writer.Append(MemoryMarshal.AsBytes(output));
            }
        }
    }
}
