using System.Diagnostics;
using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Ph5;

/// <summary>What the backward scan did.</summary>
internal sealed class NsvResult
{
    internal long StackPeak { get; set; }

    internal SpillStats Spill { get; set; }

    internal WindowStats Windows { get; set; }

    internal double Seconds { get; set; }
}

/// <summary>An entry of the NSV stack: a value and the rank where it stands (8 bytes).</summary>
[StructLayout(LayoutKind.Sequential)]
internal struct ValueEntry
{
    internal uint Value;
    internal uint Index;
}

/// <summary>
/// The backward half of layout F (ADR-004 section 3.2: "NSV by a backward scan of the LCP stream with the same stack"): next smaller value (strict) of every rank, <c>NSV[i] = min { j &gt; i : LCP[j] &lt; LCP[i] }</c> or
/// <c>n + 1</c>, written in descending rank order, so every write is sequential. The stack holds the distinct values seen so far (strictly increasing from the bottom, at most <c>1 + max LCP</c> entries) in a
/// <see cref="SpillStack{T}"/>; the LCP is read backwards from the compact file.
/// </summary>
internal static class NsvScan
{
    /// <summary>Scans <paramref name="lcp"/> backwards and writes NSV[0..n] (32-bit little endian) to <paramref name="nsvPath"/>.</summary>
    internal static NsvResult Run(long n, CompactLcpFile lcp, string nsvPath, ChildBuildOptions options)
    {
        ArgumentNullException.ThrowIfNull(lcp);
        ArgumentNullException.ThrowIfNull(options);
        Directory.CreateDirectory(options.ScratchDir);
        uint none = (uint)(n + 1);
        using var arena = new BudgetedArena(options.Budget ?? new ArenaBudget(long.MaxValue), options.ArenaGuards, options.TouchBuffers);
        using var stack = new SpillStack<ValueEntry>(arena, Path.Combine(options.ScratchDir, "ph5-nsv-" + Guid.NewGuid().ToString("N") + ".spill"), options.BlockEntries, options.SpillFlushBytes);
        using var writer = new AlignedWindowWriter(nsvPath, 4 * (n + 1), options.WindowBytes, 2, options.FlushWindows);
        long t0 = Stopwatch.GetTimestamp();
        lcp.ForEachBackward((first, values) =>
        {
            for (int j = values.Length - 1; j >= 0; j--)
            {
                uint x = values[j];
                while (stack.Height > 0 && stack.Top.Value >= x) stack.Pop();
                uint rank = (uint)(first + j);
                writer.Write32(rank, stack.Height > 0 ? stack.Top.Index : none);
                stack.Push(new ValueEntry { Value = x, Index = rank });
            }
        });
        writer.Dispose();
        return new NsvResult { StackPeak = stack.PeakHeight, Spill = stack.Stats, Windows = writer.Stats, Seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency };
    }
}
