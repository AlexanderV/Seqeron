using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Dcs;

/// <summary>
/// The output of one sample bucket: its records <c>(slot, local name)</c> packed as <c>slot &lt;&lt; 32 | name</c>, sorted by slot (hence by text position). A run of at most the spill size stays in
/// a managed array, a larger one is one scratch file of 8-byte little-endian records (counted by <see cref="IoCounters"/>, tracked by <see cref="ScratchTracker"/>, deleted on disposal).
/// </summary>
internal sealed class BucketRun : IDisposable
{
    private static long s_counter;
    private ulong[]? _memory;
    private readonly string? _path;

    private BucketRun(ulong[]? memory, string? path, long count)
    {
        _memory = memory;
        _path = path;
        Count = count;
    }

    /// <summary>Records in the run.</summary>
    internal long Count { get; }

    internal bool InMemory => _memory is not null;

    /// <summary>Creates a run from records that are already sorted by slot.</summary>
    internal static BucketRun Create(ReadOnlySpan<ulong> sorted, string dir, string tag, int spillBytes)
    {
        if ((long)sorted.Length * 8 <= spillBytes) return new BucketRun(sorted.ToArray(), null, sorted.Length);
        string path = Path.Combine(dir, $"{tag}-run-{Interlocked.Increment(ref s_counter):D8}.bin");
        using (var f = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite))
        {
            f.WriteAt(MemoryMarshal.AsBytes(sorted), 0);
        }
        ScratchTracker.Add((long)sorted.Length * 8);
        return new BucketRun(null, path, sorted.Length);
    }

    /// <summary>The records of an in-memory run.</summary>
    internal ReadOnlySpan<ulong> Memory => _memory;

    /// <summary>Reads up to <paramref name="dest"/>.Length records from record index <paramref name="start"/> of a file run (the file is opened for the call and closed again: a merge of thousands of runs must not hold thousands of descriptors).</summary>
    internal int ReadAt(long start, Span<ulong> dest)
    {
        long left = Count - start;
        int take = (int)Math.Min(dest.Length, left);
        if (take <= 0) return 0;
        using var f = new CountedFile(_path!, FileMode.Open, FileAccess.Read);
        f.ReadExactAt(MemoryMarshal.AsBytes(dest[..take]), start * 8);
        return take;
    }

    public void Dispose()
    {
        _memory = null;
        if (_path is null) return;
        ScratchTracker.Remove(Count * 8);
        try
        {
            File.Delete(_path);
        }
        catch (IOException)
        {
            // best effort: the owner of the scratch directory removes what is left
        }
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
