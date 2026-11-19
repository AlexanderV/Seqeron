using System.Runtime.InteropServices;
using System.Text;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>A scratch directory of one test (deleted with it).</summary>
internal sealed class TempDir : IDisposable
{
    internal TempDir()
    {
        Path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "esa-ext-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(Path);
    }

    internal string Path { get; }

    /// <summary>A file path inside the directory (nothing is created).</summary>
    internal string File(string name) => System.IO.Path.Combine(Path, name);

    /// <summary>A subdirectory path (nothing is created).</summary>
    internal string Sub(string name) => System.IO.Path.Combine(Path, name);

    public void Dispose()
    {
        if (Directory.Exists(Path)) Directory.Delete(Path, true);
    }
}

/// <summary>An <see cref="IProcSource"/> that serves fixed texts, or fails the way a platform without <c>/proc</c> does.</summary>
internal sealed class FakeProcSource : IProcSource
{
    private readonly Dictionary<string, Func<byte[]>> _files = [];

    /// <summary>Number of reads, answered or not.</summary>
    internal int Calls { get; private set; }

    /// <summary>When set, every read throws this exception (a hostile source).</summary>
    internal Exception? Throws { get; set; }

    internal FakeProcSource With(string path, string text)
    {
        byte[] bytes = Encoding.ASCII.GetBytes(text);
        _files[path] = () => bytes;
        return this;
    }

    /// <summary>The content is produced anew at every read.</summary>
    internal FakeProcSource WithSeries(string path, Func<int, string> text)
    {
        int call = 0;
        _files[path] = () => Encoding.ASCII.GetBytes(text(call++));
        return this;
    }

    public int TryRead(string path, Span<byte> buffer)
    {
        Calls++;
        if (Throws is not null) throw Throws;
        if (!_files.TryGetValue(path, out Func<byte[]>? make)) return -1;
        byte[] bytes = make();
        int take = Math.Min(bytes.Length, buffer.Length);
        bytes.AsSpan(0, take).CopyTo(buffer);
        return take;
    }
}

internal static class FileBits
{
    internal static void WriteUInts(string path, ReadOnlySpan<uint> values) => File.WriteAllBytes(path, MemoryMarshal.AsBytes(values).ToArray());

    internal static uint[] ReadUInts(string path) => MemoryMarshal.Cast<byte, uint>(File.ReadAllBytes(path)).ToArray();

    /// <summary>Every permutation of 0..count-1 (Heap's algorithm order is not needed: lexicographic by recursion).</summary>
    internal static IEnumerable<uint[]> Permutations(int count)
    {
        var current = Enumerable.Range(0, count).Select(i => (uint)i).ToArray();
        return Permute(current, 0);
    }

    private static IEnumerable<uint[]> Permute(uint[] a, int k)
    {
        if (k == a.Length)
        {
            yield return (uint[])a.Clone();
            yield break;
        }
        for (int i = k; i < a.Length; i++)
        {
            (a[k], a[i]) = (a[i], a[k]);
            foreach (uint[] p in Permute(a, k + 1)) yield return p;
            (a[k], a[i]) = (a[i], a[k]);
        }
    }
}

/// <summary>An <see cref="IReservationSource"/> that behaves like a platform: with or without fallocate, with a known or unknown amount of free space, or like a disk that is full.</summary>
internal sealed class FakeReservationSource : IReservationSource
{
    public bool CanPreallocate { get; set; }

    /// <summary>Free bytes reported (-1: unknown).</summary>
    public long Free { get; set; } = -1;

    /// <summary>When set, the preallocating creation throws this exception.</summary>
    public Exception? PreallocateThrows { get; set; }

    public int PreallocateCalls { get; private set; }

    public int FreeCalls { get; private set; }

    public string? LastDirectory { get; private set; }

    /// <summary>Every directory that the source was asked about, in order.</summary>
    public List<string> Asked { get; } = [];

    public FileStream CreatePreallocated(string path, long length)
    {
        PreallocateCalls++;
        if (PreallocateThrows is not null) throw PreallocateThrows;
        return new FileStream(path, FileMode.Create, FileAccess.ReadWrite, FileShare.ReadWrite);
    }

    public long AvailableFreeBytes(string directory)
    {
        FreeCalls++;
        Asked.Add(directory);
        LastDirectory = directory;
        return Free;
    }
}

/// <summary>What the kernel says about the mappings of this process (Linux: <c>/proc/self/maps</c>): lets a test see that a disposed window or reader really unmapped its view.</summary>
internal static class ProcMaps
{
    /// <summary>True when a mapping of <paramref name="path"/> exists, false when none does, null when this platform has no <c>/proc/self/maps</c>.</summary>
    internal static bool? Contains(string path)
    {
        const string maps = "/proc/self/maps";
        if (!OperatingSystem.IsLinux() || !File.Exists(maps)) return null;
        string full = System.IO.Path.GetFullPath(path);
        return File.ReadLines(maps).Any(l => l.EndsWith(full, StringComparison.Ordinal));
    }
}
