using System.Diagnostics;
using System.Runtime.InteropServices;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Ph5;

/// <summary>What the NEXTSYM pass did.</summary>
internal sealed class NextSymResult
{
    internal double Seconds { get; set; }

    /// <summary>Random reads of the text (one per rank: <c>n + 1</c>).</summary>
    internal long TextReads { get; set; }

    internal WindowStats Windows { get; set; }
}

/// <summary>
/// NEXTSYM of ADR-004 section 4.2 (row PH5): <c>NEXTSYM[r] = text[SA[r] + LCP[r]]</c>, one byte per rank, the symbol that follows the longest common prefix with the previous suffix (the dense symbol code; <see cref="NoSymbol"/> when
/// <c>SA[r] + LCP[r] = n</c>, i.e. at the end of the text, which is the case of rank 0). It is a separate pass of the prototype, report-only: it reads the suffix array file and the compact LCP file in rank order and the text at
/// random (<c>n + 1</c> reads through the read-only mapping of the page cache), and writes the output sequentially through bounded mapped windows. Symbol codes must be below 255.
/// </summary>
internal static class NextSymScan
{
    /// <summary>The stored byte when the position is the end of the text.</summary>
    internal const byte NoSymbol = 0xFF;

    internal static NextSymResult Run(PackedText text, string saPath, CompactLcpFile lcp, string outPath, ChildBuildOptions options)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(lcp);
        ArgumentNullException.ThrowIfNull(options);
        long n = text.Length;
        if (text.Bits > 8) throw new ArgumentException("NEXTSYM holds one byte per rank: the text must have at most 8 bits per symbol");
        using var sa = new CountedFile(saPath, FileMode.Open, FileAccess.Read);
        if (sa.Length != (n + 1) * 4) throw new ArgumentException($"the suffix array file has {sa.Length} bytes, expected {(n + 1) * 4} for n = {n}");
        using var writer = new AlignedWindowWriter(outPath, n + 1, options.WindowBytes, 2, options.FlushWindows);
        var positions = new uint[Math.Max(1, options.ReadBufferBytes / 4)];
        var bytes = new byte[positions.Length];
        long t0 = Stopwatch.GetTimestamp();
        lcp.ForEachForward((first, values) =>
        {
            for (int at = 0; at < values.Length; at += positions.Length)
            {
                int count = Math.Min(positions.Length, values.Length - at);
                sa.ReadExactAt(MemoryMarshal.AsBytes(positions.AsSpan(0, count)), (first + at) * 4);
                Fill(text, n, positions, values.Slice(at, count), bytes);
                writer.Write(first + at, bytes.AsSpan(0, count));
            }
        });
        writer.Dispose();
        return new NextSymResult { Seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency, TextReads = n + 1, Windows = writer.Stats };
    }

    /// <summary>The bytes of a block: <c>text[SA + LCP]</c> or <see cref="NoSymbol"/> at the end of the text.</summary>
    internal static void Fill(PackedText text, long n, ReadOnlySpan<uint> sa, ReadOnlySpan<uint> lcp, Span<byte> destination)
    {
        for (int i = 0; i < lcp.Length; i++)
        {
            long p = (long)sa[i] + lcp[i];
            destination[i] = p < n ? (byte)text.Symbol(p) : NoSymbol;
        }
    }
}
