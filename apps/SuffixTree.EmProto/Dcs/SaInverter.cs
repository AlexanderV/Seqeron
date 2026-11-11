using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Dcs;

/// <summary>
/// The inverse permutation of the suffix array of the reduced string, straight into the rank file of the level above (the DIST primitive of ADR-004 section 3.3 with a rank file as the target): pass 1 reads
/// the suffix array in rank order and distributes the tuples <c>(slot, rank)</c> of the real (non-separator) reduced positions into P = ceil(Slots / R) chunk files by <c>slot / R</c> (P write buffers; a
/// chunk that never outgrows its buffer stays in memory); pass 2 takes the chunks in slot order, scatters each into an array of R entries and appends it to the rank file through the bounded mapped
/// window writer. The rank of a sample suffix is its index in the suffix array of the reduced string, separators included (they occupy the first |D| ranks), the convention of <see cref="DcsBuilder"/>.
/// </summary>
internal static class SaInverter
{
    internal static void Invert(string saPath, long reducedLength, SampleGeometry g, string rankPath, bool scratchOutput, ExtDcsOptions o, string tag, ArenaBudget budget)
    {
        long slots = g.Slots;
        long chunkSlots = Math.Max(1, o.ChunkSlots);
        long chunks = (slots + chunkSlots - 1) / chunkSlots;
        int bufEntries = Math.Max(1, o.IoBufferBytes / 8);
        using var arena = new BudgetedArena(budget, null, o.TouchBuffers);
        Buf<ulong> tuples = arena.Allocate<ulong>(checked((int)(chunks * bufEntries)), zero: false, "invert.tuples");
        var filled = new int[chunks];
        var flushed = new long[chunks];

        // pass 1: the suffix array in rank order
        int readEntries = Math.Max(1, o.IoBufferBytes / 4);
        Buf<uint> readBuf = arena.Allocate<uint>(readEntries, zero: false, "invert.sa");
        long total = reducedLength + 1;
        using (var sa = new CountedFile(saPath, FileMode.Open, FileAccess.Read))
        {
            for (long first = 0; first < total; first += readEntries)
            {
                int take = (int)Math.Min(readEntries, total - first);
                Span<uint> entries = readBuf.AsSpan(take);
                sa.ReadExactAt(MemoryMarshal.AsBytes(entries), first * 4);
                for (int i = 0; i < take; i++)
                {
                    long rank = first + i;
                    if (rank == 0) continue; // the empty suffix
                    if (!g.TrySlotOfRed(entries[i], out long slot)) continue; // a separator
                    long chunk = slot / chunkSlots;
                    Span<ulong> buf = tuples.AsSpan(checked((int)(chunks * bufEntries))).Slice(checked((int)(chunk * bufEntries)), bufEntries);
                    buf[filled[chunk]++] = ((ulong)slot << 32) | (uint)rank;
                    if (filled[chunk] == bufEntries) Flush(chunk, buf, filled, flushed, o, tag);
                }
            }
        }

        // pass 2: chunk by chunk in slot order
        using var writer = new MappedWindowWriter(rankPath, slots * 4, o.WindowBytes);
        long outputBytes = scratchOutput ? slots * 4 : 0;
        if (outputBytes > 0) ScratchTracker.Add(outputBytes);
        Buf<uint> chunkArray = arena.Allocate<uint>((int)Math.Min(chunkSlots, slots), zero: false, "invert.chunk");
        Buf<ulong> chunkRead = arena.Allocate<ulong>(bufEntries, zero: false, "invert.read");
        for (long chunk = 0; chunk < chunks; chunk++)
        {
            long start = chunk * chunkSlots;
            Span<uint> array = chunkArray.AsSpan((int)Math.Min(chunkSlots, slots - start));
            array.Clear();
            string path = ChunkPath(o, tag, chunk);
            if (flushed[chunk] > 0)
            {
                using (var f = new CountedFile(path, FileMode.Open, FileAccess.Read))
                {
                    for (long at = 0; at < flushed[chunk]; at += bufEntries)
                    {
                        int take = (int)Math.Min(bufEntries, flushed[chunk] - at);
                        Span<ulong> part = chunkRead.AsSpan(take);
                        f.ReadExactAt(MemoryMarshal.AsBytes(part), at * 8);
                        foreach (ulong v in part) array[(int)((long)(v >> 32) - start)] = (uint)v;
                    }
                }
                ScratchTracker.Remove(flushed[chunk] * 8);
                File.Delete(path);
            }
            ReadOnlySpan<ulong> rest = tuples.AsSpan(checked((int)(chunks * bufEntries))).Slice(checked((int)(chunk * bufEntries)), filled[chunk]);
            foreach (ulong v in rest) array[(int)((long)(v >> 32) - start)] = (uint)v;
            writer.Append(MemoryMarshal.AsBytes(array));
        }
    }

    private static string ChunkPath(ExtDcsOptions o, string tag, long chunk) => Path.Combine(o.ScratchDir, $"{tag}-tup-{chunk:D6}.bin");

    private static void Flush(long chunk, ReadOnlySpan<ulong> buf, int[] filled, long[] flushed, ExtDcsOptions o, string tag)
    {
        using (var f = new CountedFile(ChunkPath(o, tag, chunk), FileMode.OpenOrCreate, FileAccess.Write))
        {
            f.WriteAt(MemoryMarshal.AsBytes(buf[..filled[chunk]]), flushed[chunk] * 8);
        }
        ScratchTracker.Add(filled[chunk] * 8L);
        flushed[chunk] += filled[chunk];
        filled[chunk] = 0;
    }
}
