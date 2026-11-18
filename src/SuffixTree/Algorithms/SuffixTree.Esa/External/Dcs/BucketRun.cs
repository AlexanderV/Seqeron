using System.Runtime.InteropServices;

namespace SuffixTree.Esa.External;

/// <summary>
/// The scratch bytes that one build of PH1 or PH2 / PH3 holds on disk, mirrored into the shared <see cref="IoCounters"/> (live and peak: ADR-004 Q1 reports the peak scratch in B/char next to the I/O volume). The ledger knows
/// what this build took, so that a failed build gives back exactly that and not what another stage of the same counters holds. Thread-safe.
/// </summary>
internal sealed class ScratchLedger(IoCounters? counters)
{
    private long _live;

    /// <summary>Bytes this build holds now.</summary>
    internal long Live => Interlocked.Read(ref _live);

    internal void Add(long bytes)
    {
        Interlocked.Add(ref _live, bytes);
        counters?.AddScratch(bytes);
    }

    internal void Remove(long bytes)
    {
        Interlocked.Add(ref _live, -bytes);
        counters?.RemoveScratch(bytes);
    }

    /// <summary>Gives back everything that is still held (the files are deleted by the failure path).</summary>
    internal void ReleaseAll() => Remove(Interlocked.Read(ref _live));
}

/// <summary>Where the runs of the sample buckets go: the directory, the tag of the build, the size up to which a run stays in memory, the numbering of the files (per build, so that two builds never share a name) and the ledger and counters.</summary>
internal sealed class RunStore(string directory, string tag, int spillBytes, ScratchLedger ledger, IoCounters? counters)
{
    private int _next;

    internal ScratchLedger Ledger { get; } = ledger;

    /// <summary>Creates a run from records that are already sorted by slot.</summary>
    internal BucketRun Create(ReadOnlySpan<ulong> sorted)
    {
        if ((long)sorted.Length * 8 <= spillBytes) return new BucketRun(sorted.ToArray(), null, sorted.Length, this);
        string path = Path.Combine(directory, $"{tag}-run-{Interlocked.Increment(ref _next):D8}.bin");
        long bytes = (long)sorted.Length * 8;
        Ledger.Add(bytes); // taken with the reservation, so that a failed write gives back exactly what was taken
        try
        {
            using var f = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite, counters);
            f.WriteAt(MemoryMarshal.AsBytes(sorted), 0);
        }
        catch
        {
            Ledger.Remove(bytes);
            DeleteQuietly(path);
            throw;
        }
        return new BucketRun(null, path, sorted.Length, this);
    }

    internal CountedFile Open(string path) => new(path, FileMode.Open, FileAccess.Read, counters);

    internal static void DeleteQuietly(string path)
    {
        try
        {
            File.Delete(path);
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException)
        {
            // best effort: the owner of the scratch directory removes what is left
        }
    }
}

/// <summary>
/// The output of one sample bucket: its records <c>(slot, local name)</c> packed as <c>slot &lt;&lt; 32 | name</c>, sorted by slot (hence by text position). A run of at most the spill size stays in a managed array, a
/// larger one is one scratch file of 8-byte little-endian records (counted by <see cref="IoCounters"/>, accounted in the <see cref="ScratchLedger"/>, deleted on disposal).
/// </summary>
internal sealed class BucketRun : IDisposable
{
    private readonly string? _path;
    private readonly RunStore _store;
    private ulong[]? _memory;
    private int _disposed;

    internal BucketRun(ulong[]? memory, string? path, long count, RunStore store)
    {
        _memory = memory;
        _path = path;
        _store = store;
        Count = count;
    }

    /// <summary>Records in the run.</summary>
    internal long Count { get; }

    internal bool InMemory => _memory is not null;

    /// <summary>The records of an in-memory run.</summary>
    internal ReadOnlySpan<ulong> Memory => _memory;

    /// <summary>Reads up to <paramref name="dest"/>.Length records from record index <paramref name="start"/> of a file run (the file is opened for the call and closed again: a merge of thousands of runs must not hold thousands of descriptors).</summary>
    internal int ReadAt(long start, Span<ulong> dest)
    {
        long left = Count - start;
        int take = (int)Math.Min(dest.Length, left);
        if (take <= 0) return 0;
        using CountedFile f = _store.Open(_path!);
        f.ReadExactAt(MemoryMarshal.AsBytes(dest[..take]), start * 8);
        return take;
    }

    public void Dispose()
    {
        if (Interlocked.Exchange(ref _disposed, 1) != 0) return;
        _memory = null;
        if (_path is null) return;
        _store.Ledger.Remove(Count * 8);
        RunStore.DeleteQuietly(_path);
    }
}

/// <summary>Reads a <see cref="BucketRun"/> in order through a bounded buffer and scatters its records into the name window of the merge.</summary>
internal sealed class RunReader(BucketRun run, Buf<ulong> buffer, int bufferRecords)
{
    private long _consumed;
    private int _pos;
    private int _len;

    /// <summary>Records scattered so far.</summary>
    internal long Scattered { get; private set; }

    /// <summary>
    /// Writes <c>offset + local name</c> into <c>window[slot - windowStart]</c> for every following record whose slot is below <paramref name="limit"/>; stops at the first record at or above it.
    /// Returns the number of records scattered.
    /// </summary>
    internal long Drain(uint limit, uint offset, Span<uint> window, long windowStart)
    {
        long drained = 0;
        while (true)
        {
            if (_pos == _len && !Refill()) return drained;
            ReadOnlySpan<ulong> cur = run.InMemory ? run.Memory : buffer.AsSpan(_len);
            int i = _pos;
            while (i < _len)
            {
                ulong v = cur[i];
                uint slot = (uint)(v >> 32);
                if (slot >= limit) break;
                window[(int)(slot - windowStart)] = offset + (uint)v;
                i++;
            }
            drained += i - _pos;
            Scattered += i - _pos;
            _pos = i;
            if (_pos < _len) return drained;
        }
    }

    private bool Refill()
    {
        if (_consumed >= run.Count) return false;
        if (run.InMemory)
        {
            _len = (int)run.Count;
            _consumed = run.Count;
        }
        else
        {
            _len = run.ReadAt(_consumed, buffer.AsSpan(bufferRecords));
            _consumed += _len;
        }
        _pos = 0;
        return _len > 0;
    }
}
