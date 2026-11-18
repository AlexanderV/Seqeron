using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>A scratch directory of one test.</summary>
internal sealed class DiffDir : IDisposable
{
    internal DiffDir()
    {
        Path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "emproto-1b7a-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(Path);
    }

    internal string Path { get; }

    internal string File(string name) => System.IO.Path.Combine(Path, name);

    public void Dispose()
    {
        if (Directory.Exists(Path)) Directory.Delete(Path, true);
    }
}

/// <summary>The inputs of the differential tests: generated texts, the data files of the benchmark directory (skipped when absent) and the suffix-array files made from them.</summary>
internal static class DiffInputs
{
    /// <summary>The benchmark directory: <c>SEQERON_BENCH_DIR</c> (default <c>/tmp/stbench</c>) + <c>/data</c>.</summary>
    internal static string DataDir => System.IO.Path.Combine(Environment.GetEnvironmentVariable("SEQERON_BENCH_DIR") ?? "/tmp/stbench", "data");

    /// <summary>The text of a data file (<c>dna_10m</c>, <c>prot_1m</c>, <c>rep_10m</c>), or ignores the test when the file is not there.</summary>
    internal static string DataText(string name)
    {
        string path = System.IO.Path.Combine(DataDir, name + ".txt");
        if (!File.Exists(path)) Assert.Ignore($"{path} is not present (generate it with analysis/bench/gen/gen_inputs.py)");
        return File.ReadAllText(path).TrimEnd('\n', '\r');
    }

    internal static IEnumerable<(string Name, string Text)> Families(int n) =>
    [
        ("a^n", StringFamilies.AllSame(n)),
        ("fibonacci", StringFamilies.Fibonacci(n)),
        ("thue-morse", StringFamilies.ThueMorse(n)),
    ];

    /// <summary>The suffix array of the text with the empty suffix first (n + 1 entries), written as 32-bit little-endian values.</summary>
    internal static string WriteSa(DiffDir dir, string text, string name = "sa.bin")
    {
        int[] sa = EsaIndex.Build(text, ChildLayoutKind.Tagged).SuffixArray;
        string path = dir.File(name);
        File.WriteAllBytes(path, MemoryMarshal.AsBytes(sa.AsSpan()).ToArray());
        return path;
    }

    internal static bool SameBytes(string a, string b)
    {
        if (new FileInfo(a).Length != new FileInfo(b).Length) return false;
        using var fa = new FileStream(a, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 16);
        using var fb = new FileStream(b, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 16);
        var ba = new byte[1 << 16];
        var bb = new byte[1 << 16];
        int got;
        while ((got = fa.Read(ba, 0, ba.Length)) > 0)
        {
            fb.ReadExactly(bb, 0, got);
            if (!ba.AsSpan(0, got).SequenceEqual(bb.AsSpan(0, got))) return false;
        }
        return true;
    }
}
