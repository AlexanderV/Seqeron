using SuffixTree.Esa;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Ph4;

/// <summary>How PH4 gets from the Phi tuples to the LCP in rank order (ADR-004 section 4.2, row PH4 and the bounded-scratch paragraph).</summary>
internal enum LcpMode
{
    /// <summary>Base mode of the ADR: three distributions, tuples <c>(SA[i], SA[i-1], i)</c> by <c>SA[i] / R</c>, then <c>(i, PLCP[j])</c> by <c>i / R</c>; scratch peaks at 12 B/char.</summary>
    Distributed,

    /// <summary>Bounded-scratch mode: text-range rounds, Phi tuples <c>(SA[i], SA[i-1])</c> of one round only, PLCP of every round into a text-order file, then a gather <c>PLCP[SA[i]]</c> in rank order.</summary>
    Bounded,
}

/// <summary>Receives the dense LCP values chunk by chunk in strictly ascending rank order, once the width of the LCP is known (the stream point of the child-table scan PH5, which is not part of this prototype).</summary>
internal interface ILcpChunkConsumer
{
    /// <summary>Called with the LCP values of ranks <c>firstRank .. firstRank + lcp.Length - 1</c> (resolved values, never escapes); the span is valid for the call only.</summary>
    void Consume(long firstRank, ReadOnlySpan<uint> lcp);
}

/// <summary>Parameters of PH4 (ADR-004 section 4.2): everything that sets a chunk size or a buffer size. The output never depends on any of them (owner answer Q5).</summary>
internal sealed record ExtLcpOptions
{
    /// <summary>R: text positions per Phi chunk and ranks per LCP chunk (ADR: 8 M, P = n / R = 375 chunks at 3x10^9 symbols).</summary>
    internal long ChunkPositions { get; init; } = 8_000_000;

    /// <summary>Upper bound of the parallel B2 tasks.</summary>
    internal int Threads { get; init; } = 4;

    /// <summary>B2 tasks run in parallel; 0 = the number the arena budget allows (each holds R Phi entries and R rank entries).</summary>
    internal int Workers { get; init; }

    /// <summary>Size of a sequential read buffer (suffix array blocks, tuple file blocks).</summary>
    internal int ReadBufferBytes { get; init; } = 1 << 20;

    /// <summary>Size of the write buffer kept per chunk by each distribution.</summary>
    internal int TupleBufferBytes { get; init; } = 32 * 1024;

    /// <summary>Size of one sliding window of the mapped output writers.</summary>
    internal int WindowBytes { get; init; } = 1 << 20;

    /// <summary>The mode of the construction.</summary>
    internal LcpMode Mode { get; init; } = LcpMode.Distributed;

    /// <summary>Bounded mode: text chunks per round (0 = all chunks in one round).</summary>
    internal long RoundChunks { get; init; }

    /// <summary>Makes every native buffer resident when it is allocated (worst case of the working set, for the RssAnon figures of the scale-down emulation).</summary>
    internal bool TouchBuffers { get; init; }

    /// <summary>Shared cap of the native buffers (null = unlimited).</summary>
    internal ArenaBudget? Budget { get; init; }

    /// <summary>Directory of the scratch files (created when missing; every file of a build is removed when it ends).</summary>
    internal required string ScratchDir { get; init; }
}

/// <summary>What one PH4 build did.</summary>
internal sealed class ExtLcpStats
{
    /// <summary>The width chosen before any LCP entry was written.</summary>
    internal LcpKind Kind { get; set; }

    internal long N { get; set; }

    /// <summary>Values at or above 255 (counted on PLCP, equal to the count on LCP).</summary>
    internal long Escapes8 { get; set; }

    /// <summary>Values at or above 65,535.</summary>
    internal long Escapes16 { get; set; }

    /// <summary>Entries of the overflow list.</summary>
    internal long OverflowCount { get; set; }

    /// <summary>Phi chunks (<c>ceil(n / R)</c>).</summary>
    internal long TextChunks { get; set; }

    /// <summary>LCP chunks (<c>ceil((n + 1) / R)</c>).</summary>
    internal long RankChunks { get; set; }

    internal int Rounds { get; set; }

    /// <summary>B2 tasks that ran in parallel.</summary>
    internal int Workers { get; set; }

    internal double B1Seconds { get; set; }

    internal double B2Seconds { get; set; }

    internal double B3Seconds { get; set; }

    /// <summary>Vector (or word) compare steps of the PLCP loops of all chunks.</summary>
    internal long CompareSteps { get; set; }

    /// <summary>Sum over the chunks of the PLCP value at the first position of the chunk: the symbols the restart without a carry compares.</summary>
    internal long RestartSymbols { get; set; }

    /// <summary>Largest PLCP value.</summary>
    internal long MaxPlcp { get; set; }
}

/// <summary>PH4 of ADR-004: the LCP array of a text from its suffix array file with bounded RAM (prototype E4).</summary>
internal static class ExternalLcp
{
    /// <summary>
    /// Builds the compact LCP of <paramref name="text"/> from the suffix array file <paramref name="saPath"/> (<c>n + 1</c> little-endian <c>uint</c> entries, <c>SA[0] = n</c>) into
    /// <paramref name="basePath"/> (<c>n + 1</c> entries of 1, 2 or 4 bytes, little endian, escape value 255 / 65,535 in the narrow kinds) and <paramref name="overflowPath"/> (pairs
    /// <c>(rank, value)</c> of 4 bytes each, ascending ranks, ADR-003 section 4 kind 5). The output files are created after the width is known and only then; the scratch files of the build are
    /// removed when it ends, also when it fails.
    /// </summary>
    /// <exception cref="ArgumentException">The suffix array file has not <c>n + 1</c> entries, <c>SA[0] != n</c>, or an entry is outside <c>[0, n)</c>.</exception>
    /// <exception cref="InvalidOperationException">The suffix array is not a permutation (a chunk holds the wrong number of tuples).</exception>
    internal static ExtLcpStats Build(string saPath, PackedText text, string basePath, string overflowPath, ExtLcpOptions options, ILcpChunkConsumer? consumer = null)
        => new LcpBuild(saPath, text, basePath, overflowPath, options, consumer).Run();

    /// <summary>
    /// Bytes of native buffers one B2 task holds for a text of <paramref name="n"/> symbols under <paramref name="options"/> (the arithmetic that limits the number of parallel tasks under a budget):
    /// the Phi chunk, in the base mode the rank chunk and the write buffers of the second distribution, the tuple read buffer.
    /// </summary>
    internal static long WorkerBytes(long n, ExtLcpOptions options) => LcpBuild.WorkerBytes(n, options);
}
