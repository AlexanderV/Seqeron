using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>Parameters of <see cref="Dist"/>.</summary>
internal sealed record DistOptions
{
    /// <summary>Directory of the chunk files; created if missing, free of this run's files again afterwards.</summary>
    internal string ScratchDir { get; init; } = "";

    /// <summary>R: slots per chunk, P = ceil(slots / R) chunk files (ADR-004: R = 8 M, P = 375 at 3x10^9).</summary>
    internal long ChunkEntries { get; init; } = 8_000_000;

    /// <summary>Write buffer of one chunk (ADR-004: 64-256 KB); a chunk that never outgrows it stays in memory.</summary>
    internal int IoBufferBytes { get; init; } = 64 * 1024;

    /// <summary>Bytes of the source read at a time.</summary>
    internal int ReadBytes { get; init; } = 1 << 20;

    /// <summary>Window of the mapped target writer (at least 1 MiB: see <see cref="MappedWindowWriter"/>).</summary>
    internal int WindowBytes { get; init; } = MappedWindowWriter.DefaultWindowBytes;

    /// <summary>The shared RAM cap; none = unlimited.</summary>
    internal ArenaBudget? Budget { get; init; }

    /// <summary>Receives the explicit and mapped I/O and the scratch bytes of the run.</summary>
    internal IoCounters? Counters { get; init; }

    /// <summary>Arena options: guard bands in the tests, none in production.</summary>
    internal ArenaOptions? Arena { get; init; }

    /// <summary>How the target file is created (reservation of the disk space through <see cref="MappedOutputFile.Create"/>: a full disk is an <see cref="IOException"/> at creation); null = created sparse, as in 1b.7a (work package 1b.7d added the option for the rank file of PH1).</summary>
    internal MappedFileOptions? Files { get; init; }
}

/// <summary>What one run wrote to and read from its chunk files (a chunk that never outgrew its buffer is in neither).</summary>
internal readonly record struct DistStats(long Chunks, long TupleBytesWritten, long TupleBytesRead, ReservationOutcome Reservation = ReservationOutcome.Sparse);

/// <summary>Maps the entry <c>index</c> of the source (a file of 32-bit values, in file order) to a tuple <c>(key, value)</c>, or skips it. A mapper may keep state (it is passed by reference).</summary>
internal interface IDistMapper<TValue>
    where TValue : unmanaged
{
    bool TryMap(long index, uint sourceValue, out uint key, out TValue value);
}

/// <summary>
/// The target array of a scatter: <paramref name="Slots"/> entries of <typeparamref name="TValue"/>. A slot nobody writes keeps <paramref name="Fill"/>. With <paramref name="ExactlyOnce"/> every slot must be
/// written exactly once (a permutation: ISA, Phi); a repeat or a miss throws <see cref="InvalidDataException"/>, which needs a <paramref name="Fill"/> that is not a legal value.
/// </summary>
internal readonly record struct DistTarget<TValue>(long Slots, TValue Fill, bool ExactlyOnce)
    where TValue : unmanaged;

/// <summary>A (key, value) tuple of the chunk files; 4 + sizeof(value) bytes (8 for a 32-bit value).</summary>
[StructLayout(LayoutKind.Sequential, Pack = 1)]
internal readonly struct DistEntry<TValue>(uint key, TValue value)
    where TValue : unmanaged
{
    internal readonly uint Key = key;
    internal readonly TValue Value = value;
}

/// <summary>
/// <c>key / R</c> for a 32-bit key and a fixed R, without a hardware division (it runs once per source entry, and a 64-bit division costs 25-90 cycles on the machines of this class): for R from 2 to 2^32 - 1 two
/// multiplications of Lemire, Kaser and Kurz ("Faster remainder by direct computation", 2019: with M = ceil(2^64 / R) the high word of M x key is exactly key / R for every 32-bit key), R = 1 is the identity and an R
/// above 32 bits puts every key in chunk 0.
/// </summary>
internal readonly struct ChunkDivider
{
    private readonly ulong _m;
    private readonly uint _mode; // 0 = single chunk, 1 = identity, 2 = multiply

    internal ChunkDivider(long chunkEntries)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(chunkEntries, 1);
        if (chunkEntries > uint.MaxValue)
        {
            _mode = 0;
            _m = 0;
        }
        else if (chunkEntries == 1)
        {
            _mode = 1;
            _m = 0;
        }
        else
        {
            _mode = 2;
            _m = (ulong.MaxValue / (uint)chunkEntries) + 1;
        }
    }

    /// <summary>The chunk of <paramref name="key"/>: <c>key / R</c>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal uint Of(uint key)
    {
        if (_mode == 2) return (uint)Math.BigMul(_m, key, out _);
        return _mode == 1 ? key : 0;
    }
}

/// <summary>ISA: key <c>SA[i]</c>, value <c>i</c>; refuses an entry beyond n.</summary>
internal readonly struct IsaMapper(long n) : IDistMapper<uint>
{
    internal long N { get; } = n;

    public bool TryMap(long index, uint sourceValue, out uint key, out uint value)
    {
        if (sourceValue > N) ThrowBeyondN(index, sourceValue);
        key = sourceValue;
        value = (uint)index;
        return true;
    }

    [DoesNotReturn]
    private void ThrowBeyondN(long index, uint sourceValue) => throw new InvalidDataException($"SA[{index}] = {sourceValue} lies beyond n = {N}");
}

/// <summary>Phi: key <c>SA[i]</c>, value <c>SA[i-1]</c> for i &gt;= 1; the entry 0 (the empty suffix, which must be n) only starts the chain.</summary>
internal struct PhiMapper(long n) : IDistMapper<uint>
{
    private uint _previous;

    internal long N { get; } = n;

    public bool TryMap(long index, uint sourceValue, out uint key, out uint value)
    {
        key = 0;
        value = 0;
        if (index == 0)
        {
            if (sourceValue != N) throw new InvalidDataException($"SA[0] must be n = {N}, the empty suffix, it is {sourceValue}");
            _previous = sourceValue;
            return false;
        }
        if (sourceValue >= N) throw new InvalidDataException($"SA[{index}] = {sourceValue} is not a position below n = {N}");
        key = sourceValue;
        value = _previous;
        _previous = sourceValue;
        return true;
    }
}

/// <summary>
/// DIST, the one-level radix distribute-and-scatter primitive of ADR-004 section 3.3, generalised in the tuple and in the target as section 4.4 asks. Pass 1 reads a file of 32-bit values in file order, maps every entry
/// to a tuple <c>(key, value)</c> and distributes it into P = ceil(slots / R) chunk files by <c>key / R</c> through P write buffers (a chunk that never outgrows its buffer stays in memory); pass 2 takes the chunks in
/// key order, scatters each into an array of R entries of the target type and appends it to the target file through a bounded mapped window. Used by PH6 (ISA, <see cref="BuildIsa"/>), PH4 (Phi, <see cref="BuildPhi"/>;
/// the per-chunk PLCP work of PH4 uses <see cref="DistScratch{T}"/> and <see cref="DistWriter{T}"/> directly) and PH1 (the rank file of the sample, a skipping mapper and a fill of zero). I/O of ISA by DIST: source read 4 +
/// tuples written 8 + tuples read 8 + target written 4 = 24 B per entry. Memory: the P buffers, a read buffer in pass 1, the chunk array and a read buffer in pass 2, all from one arena under the budget. Keys are 32-bit.
/// </summary>
internal static class Dist
{
    /// <summary>Number of chunks of <paramref name="chunkEntries"/> slots needed for <paramref name="slots"/> slots.</summary>
    internal static long ChunksFor(long slots, long chunkEntries)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(slots);
        ArgumentOutOfRangeException.ThrowIfLessThan(chunkEntries, 1);
        return (slots / chunkEntries) + (slots % chunkEntries == 0 ? 0 : 1);
    }

    /// <summary>
    /// Runs both passes: the source file of <paramref name="sourceEntries"/> 32-bit values is mapped to tuples, the target file of <c>target.Slots x sizeof(TValue)</c> bytes is written. The call returns
    /// when the target is complete; on an exception the target file may be partial and must not be used, the chunk files are removed and the arena budget is given back in every case.
    /// </summary>
    internal static unsafe DistStats Scatter<TValue, TMapper>(string sourcePath, long sourceEntries, ref TMapper mapper, string targetPath, DistTarget<TValue> target, DistOptions options)
        where TValue : unmanaged, IEquatable<TValue>
        where TMapper : IDistMapper<TValue>
    {
        ArgumentNullException.ThrowIfNull(options);
        ArgumentOutOfRangeException.ThrowIfNegative(sourceEntries);
        ArgumentOutOfRangeException.ThrowIfNegative(target.Slots);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ChunkEntries, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.IoBufferBytes, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.ReadBytes, 4);
        ArgumentOutOfRangeException.ThrowIfLessThan(options.WindowBytes, 1);
        ArgumentException.ThrowIfNullOrEmpty(options.ScratchDir);
        EnsureKeysFit(target.Slots);

        long slots = target.Slots;
        long chunkEntries = options.ChunkEntries;
        long chunks = ChunksFor(slots, chunkEntries);
        int tupleBytes = sizeof(DistEntry<TValue>);
        int bufEntries = Math.Max(1, options.IoBufferBytes / tupleBytes);
        int readEntries = (int)Math.Min(options.ReadBytes / 4, sourceEntries);
        IoCounters? counters = options.Counters;
        var divide = new ChunkDivider(chunkEntries);

        using var source = new CountedFile(sourcePath, FileMode.Open, FileAccess.Read, counters);
        if (source.Length != 4 * sourceEntries) throw new InvalidDataException($"the source file has {source.Length} bytes, expected {4 * sourceEntries} for {sourceEntries} entries");
        Directory.CreateDirectory(options.ScratchDir);

        using var scratch = new DistScratch<DistEntry<TValue>>(options.ScratchDir, "dist", chunks, counters);
        using var arena = new BudgetedArena(options.Budget ?? new ArenaBudget(long.MaxValue), options.Arena);
        var writer = new DistWriter<DistEntry<TValue>>(scratch, arena, 0, chunks, options.IoBufferBytes);

        // pass 1: the source in file order; its read buffer is released before pass 2 allocates the chunk array
        int mark = arena.Mark();
        Buf<uint> readBuf = arena.Allocate<uint>(readEntries, zero: false, "dist.src");
        Distribute(source, sourceEntries, ref mapper, readBuf, readEntries, writer, divide, slots);
        arena.ReleaseTo(mark);

        // pass 2: chunk by chunk in key order
        using MappedOutputFile? created = options.Files is null ? null : MappedOutputFile.Create(targetPath, slots * sizeof(TValue), options.Files);
        using MappedWindowWriter window = created is null ? new MappedWindowWriter(targetPath, slots * sizeof(TValue), options.WindowBytes, counters) : created.CreateWriter(options.WindowBytes, counters);
        Buf<TValue> chunkArray = arena.Allocate<TValue>((int)Math.Min(chunkEntries, slots), zero: false, "dist.chunk");
        Buf<DistEntry<TValue>> chunkRead = arena.Allocate<DistEntry<TValue>>(bufEntries, zero: false, "dist.read");
        long bytesWritten = 0;
        long bytesRead = 0;
        for (long chunk = 0; chunk < chunks; chunk++)
        {
            long startSlot = chunk * chunkEntries;
            int length = (int)Math.Min(chunkEntries, slots - startSlot);
            Span<TValue> array = chunkArray.AsSpan(length);
            array.Fill(target.Fill);
            long flushed = scratch.Entries(chunk);
            if (flushed > 0)
            {
                using (DistReader<DistEntry<TValue>> reader = scratch.OpenReader(chunk))
                {
                    int got;
                    while ((got = reader.Read(chunkRead.AsSpan(bufEntries))) > 0)
                    {
                        bytesRead += (long)got * tupleBytes;
                        ScatterInto(chunkRead.AsSpan(got), array, startSlot, target);
                    }
                }
                bytesWritten += flushed * tupleBytes;
                scratch.Delete(chunk);
            }
            ScatterInto(writer.PendingSpan(chunk), array, startSlot, target);
            if (target.ExactlyOnce)
            {
                for (int k = 0; k < length; k++)
                {
                    if (array[k].Equals(target.Fill)) throw new InvalidDataException($"slot {startSlot + k} is never written: the source is not a permutation");
                }
            }
            window.Append(MemoryMarshal.AsBytes(array));
        }
        arena.VerifyCanaries();
        return new DistStats(chunks, bytesWritten, bytesRead);
    }

    /// <summary>Keys are 32-bit: a target of more than 2^32 slots is refused.</summary>
    internal static void EnsureKeysFit(long slots)
    {
        if (slots > 1L << 32) throw new NotSupportedException("keys are 32-bit: at most 2^32 slots");
    }

    /// <summary>ISA by DIST: <c>ISA[SA[i]] = i</c> for the suffix array file of a text of n symbols (n + 1 entries, SA[0] = n), 4 (n + 1) bytes; the file must be a permutation of 0..n.</summary>
    internal static DistStats BuildIsa(string saPath, long n, string isaPath, DistOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        if (n >= uint.MaxValue - 2) throw new NotSupportedException("positions are 32-bit: n must be below 2^32 - 2");
        var mapper = new IsaMapper(n);
        return Scatter(saPath, n + 1, ref mapper, isaPath, new DistTarget<uint>(n + 1, uint.MaxValue, true), options);
    }

    /// <summary>Phi by DIST: <c>Phi[SA[i]] = SA[i-1]</c> for i = 1..n, 4 n bytes (the suffix array file has n + 1 entries and SA[0] = n).</summary>
    internal static DistStats BuildPhi(string saPath, long n, string phiPath, DistOptions options)
    {
        ArgumentNullException.ThrowIfNull(options);
        ArgumentOutOfRangeException.ThrowIfNegative(n);
        if (n >= uint.MaxValue - 2) throw new NotSupportedException("positions are 32-bit: n must be below 2^32 - 2");
        var mapper = new PhiMapper(n);
        return Scatter(saPath, n + 1, ref mapper, phiPath, new DistTarget<uint>(n, uint.MaxValue, true), options);
    }

    /// <summary>
    /// Pass 1: the hot loop of the primitive, one iteration per source entry (10^8 at the size of the E9 runs). It is a method of its own and compiled fully optimised at its first call: inside the large
    /// <see cref="Scatter{TValue,TMapper}"/> the runtime would run it as an on-stack replacement of an unoptimised frame, which keeps the loop variables in memory (measured: about 10 % of the time of ISA by DIST).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static void Distribute<TValue, TMapper>(CountedFile source, long sourceEntries, ref TMapper mapper, Buf<uint> readBuf, int readEntries, DistWriter<DistEntry<TValue>> writer, ChunkDivider divide, long slots)
        where TValue : unmanaged
        where TMapper : IDistMapper<TValue>
    {
        for (long first = 0; first < sourceEntries; first += readEntries)
        {
            int take = (int)Math.Min(readEntries, sourceEntries - first);
            Span<uint> entries = readBuf.AsSpan(take);
            source.ReadExactAt(MemoryMarshal.AsBytes(entries), first * 4);
            for (int i = 0; i < take; i++)
            {
                if (!mapper.TryMap(first + i, entries[i], out uint key, out TValue value)) continue;
                if (key >= slots) ThrowKeyOutside(first + i, key, slots);
                writer.Add(divide.Of(key), new DistEntry<TValue>(key, value));
            }
        }
    }

    [DoesNotReturn]
    private static void ThrowKeyOutside(long index, uint key, long slots)
        => throw new InvalidDataException($"entry {index} maps to key {key}, outside the {slots} slots of the target");

    [DoesNotReturn]
    private static void ThrowChunkDamaged(uint key, long startSlot)
        => throw new InvalidDataException($"key {key} does not belong to the chunk that starts at slot {startSlot}: the chunk file is damaged");

    [DoesNotReturn]
    private static void ThrowWrittenTwice(uint key)
        => throw new InvalidDataException($"slot {key} is written twice: the source is not a permutation");

    /// <summary>Scatters a run of tuples into the array of the chunk that starts at slot <paramref name="startSlot"/>; a key outside the chunk (a damaged chunk file) or, with <c>ExactlyOnce</c>, a slot written twice throws <see cref="InvalidDataException"/>.</summary>
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    internal static void ScatterInto<TValue>(ReadOnlySpan<DistEntry<TValue>> tuples, Span<TValue> array, long startSlot, DistTarget<TValue> target)
        where TValue : unmanaged, IEquatable<TValue>
    {
        bool exactlyOnce = target.ExactlyOnce;
        TValue fill = target.Fill;
        foreach (DistEntry<TValue> e in tuples)
        {
            long at = e.Key - startSlot;
            if ((ulong)at >= (ulong)array.Length) ThrowChunkDamaged(e.Key, startSlot);
            ref TValue slot = ref array[(int)at];
            if (exactlyOnce && !slot.Equals(fill)) ThrowWrittenTwice(e.Key);
            slot = e.Value;
        }
    }
}
