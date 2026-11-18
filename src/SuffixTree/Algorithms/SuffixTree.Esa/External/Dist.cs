#pragma warning disable S2325, CA1822, S1144, IDE0051, IDE0060, S1172 // red stub of work package 1b.7a: every member throws NotImplementedException; the implementation replaces this file

using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>Parameters of <see cref="Dist"/>.</summary>
internal sealed record DistOptions
{
    /// <summary>Directory of the chunk files; created if missing, empty again afterwards.</summary>
    internal string ScratchDir { get; init; } = "";

    /// <summary>R: slots per chunk (ADR-004: 8 M).</summary>
    internal long ChunkEntries { get; init; } = 8_000_000;

    /// <summary>Write buffer of one chunk (ADR-004: 64-256 KB).</summary>
    internal int IoBufferBytes { get; init; } = 64 * 1024;

    /// <summary>Bytes of the source read at a time.</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    /// <summary>Window of the mapped target writer.</summary>
    internal int WindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    internal ArenaBudget? Budget { get; init; }

    internal IoCounters? Counters { get; init; }

    /// <summary>Guard bands for the tests.</summary>
    internal ArenaOptions? Arena { get; init; }
}

/// <summary>What one run wrote and read of its chunk files.</summary>
internal readonly record struct DistStats(long Chunks, long TupleBytesWritten, long TupleBytesRead);

/// <summary>Maps entry <c>index</c> of the source (a file of 32-bit values) to a tuple, or skips it.</summary>
internal interface IDistMapper<TValue>
    where TValue : unmanaged
{
    bool TryMap(long index, uint sourceValue, out uint key, out TValue value);
}

/// <summary>The target array: <paramref name="Slots"/> entries of <typeparamref name="TValue"/>.</summary>
internal readonly record struct DistTarget<TValue>(long Slots, TValue Fill, bool ExactlyOnce)
    where TValue : unmanaged;

/// <summary>A (key, value) tuple of the chunk files; 8 bytes for a 32-bit value.</summary>
[StructLayout(LayoutKind.Sequential, Pack = 1)]
internal readonly struct DistEntry<TValue>(uint key, TValue value)
    where TValue : unmanaged
{
    internal readonly uint Key = key;
    internal readonly TValue Value = value;
}

/// <summary>ISA: key <c>SA[i]</c>, value <c>i</c>.</summary>
internal readonly struct IsaMapper(long n) : IDistMapper<uint>
{
    internal long N { get; } = n;

    public bool TryMap(long index, uint sourceValue, out uint key, out uint value) => throw new NotImplementedException();
}

/// <summary>Phi: key <c>SA[i]</c>, value <c>SA[i-1]</c> for i &gt;= 1 (the entry 0, the empty suffix, only starts the chain).</summary>
internal struct PhiMapper(long n) : IDistMapper<uint>
{
    internal long N { get; } = n;

    public bool TryMap(long index, uint sourceValue, out uint key, out uint value) => throw new NotImplementedException();
}

/// <summary>The one-level radix distribute-and-scatter primitive.</summary>
internal static class Dist
{
    /// <summary>Number of chunks of <paramref name="chunkEntries"/> slots needed for <paramref name="slots"/> slots.</summary>
    internal static long ChunksFor(long slots, long chunkEntries) => throw new NotImplementedException();

    internal static DistStats Scatter<TValue, TMapper>(string sourcePath, long sourceEntries, ref TMapper mapper, string targetPath, DistTarget<TValue> target, DistOptions options)
        where TValue : unmanaged, IEquatable<TValue>
        where TMapper : IDistMapper<TValue>
        => throw new NotImplementedException();

    /// <summary>ISA by DIST: <c>ISA[SA[i]] = i</c> for the suffix array file of a text of n symbols (n + 1 entries), 4 (n + 1) bytes.</summary>
    internal static DistStats BuildIsa(string saPath, long n, string isaPath, DistOptions options) => throw new NotImplementedException();

    /// <summary>Phi by DIST: <c>Phi[SA[i]] = SA[i-1]</c> for i = 1..n, 4 n bytes (the suffix array file has n + 1 entries, SA[0] = n).</summary>
    internal static DistStats BuildPhi(string saPath, long n, string phiPath, DistOptions options) => throw new NotImplementedException();
}
