using SuffixTree.EmProto.Harness;

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

/// <summary>The inverse permutation of the suffix array by one-level radix distribution (the DIST primitive of ADR-004 section 3.3): <c>ISA[SA[i]] = i</c>, 4 (n + 1) bytes.</summary>
internal static class IsaByDist
{
    internal static DistStats Build(string saPath, long n, string isaPath, DistOptions o) => throw new NotImplementedException();
}

/// <summary>The other route to Psi: <c>PSI[i] = ISA[SA[i] + 1]</c> from the suffix array and the ISA file (PSI[0] = 0).</summary>
internal static class PsiFromIsa
{
    internal static void Build(string saPath, string isaPath, long n, string psiPath, int readBytes = 1 << 20, int windowBytes = 1 << 20, ArenaBudget? budget = null) => throw new NotImplementedException();
}
