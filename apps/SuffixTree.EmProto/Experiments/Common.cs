using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Gen;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Experiments;

/// <summary>Aborts the process when a run exceeds its time limit (E1: a run above 10 minutes counts as a failure); prints one row first.</summary>
internal sealed class Watchdog : IDisposable
{
    private readonly Timer? _timer;

    internal Watchdog(double seconds, Report report, string experiment, string caseName)
    {
        if (seconds <= 0) return;
        _timer = new Timer(
            _ =>
            {
                report.Row(experiment, caseName, "status", "TIMEOUT", "", $"aborted after {seconds:F0} s (the limit of the gate)");
                Environment.Exit(3);
            },
            null,
            TimeSpan.FromSeconds(seconds),
            Timeout.InfiniteTimeSpan);
    }

    public void Dispose() => _timer?.Dispose();
}

/// <summary>Loads the benchmark inputs and generates the families by a short specification (<c>dna_10m</c>, <c>a:1000000</c>, <c>fib:1000000</c>, <c>tm:1000000</c>, <c>tandem:171x10000</c>, <c>nrun:10000000:1000000</c>, <c>genome:10000000</c>, <c>dna:100000000</c>).</summary>
internal static class FamilyLoader
{
    internal static SymbolText Load(string spec, string dataDir)
    {
        string[] p = spec.Split(':');
        switch (p[0])
        {
            case "a": return Families.AllSame(long.Parse(p[1], CultureInfo.InvariantCulture));
            case "fib": return Families.Fibonacci(long.Parse(p[1], CultureInfo.InvariantCulture));
            case "tm": return Families.ThueMorse(long.Parse(p[1], CultureInfo.InvariantCulture));
            case "dna": return Families.RandomDna(long.Parse(p[1], CultureInfo.InvariantCulture), 20240001);
            case "genome": return Families.GenomeLike(long.Parse(p[1], CultureInfo.InvariantCulture), 20240002);
            case "tandem":
                string[] u = p[1].Split('x');
                return Families.Tandem(int.Parse(u[0], CultureInfo.InvariantCulture), int.Parse(u[1], CultureInfo.InvariantCulture), 20240003);
            case "nrun":
                long n = long.Parse(p[1], CultureInfo.InvariantCulture);
                long len = long.Parse(p[2], CultureInfo.InvariantCulture);
                return Families.NRun(n, (n - len) / 2, len, 20240004);
            default:
                return Families.FromFile(Path.Combine(dataDir, spec + ".txt"), spec);
        }
    }

    /// <summary>The corpus name of <c>analysis/sa_vectors.tsv</c> that belongs to a specification, or null when libsais has no hash for it.</summary>
    internal static string? VectorName(string spec)
    {
        string[] p = spec.Split(':');
        return p[0] switch
        {
            "a" => $"all-same-n{p[1]}",
            "fib" => $"fibonacci-n{p[1]}",
            "tm" => $"thue-morse-n{p[1]}",
            "dna_1m" or "dna_10m" or "rep_10m" or "prot_1m" => p[0],
            _ => null,
        };
    }
}

/// <summary>The committed libsais hashes (<c>analysis/sa_vectors.tsv</c>): corpus name to SHA-256 of SA[1..n] as little-endian int32.</summary>
internal static class VectorTable
{
    internal static Dictionary<string, string> Load(string path)
    {
        var d = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (string line in File.ReadLines(path))
        {
            if (line.Length == 0 || line[0] == '#') continue;
            string[] c = line.Split('\t');
            if (c.Length >= 6 && c[0] != "corpus") d[c[0]] = c[4];
        }
        return d;
    }
}

/// <summary>Correctness checks of a finished suffix array against the oracles allowed by the rules: the committed libsais hash, or <see cref="SaisKernel"/> run on the same text.</summary>
internal static class Oracle
{
    /// <summary>The in-repository kernel on the same symbols; also reports its time as a speed reference.</summary>
    internal static int[] SaisKernelSa(SymbolText t, out double seconds)
    {
        var prepared = new PreparedText(t.Letters.ToCharArray(), t.Symbols, null);
        long t0 = Stopwatch.GetTimestamp();
        int[] sa = SaisKernel.Build(prepared);
        seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        return sa;
    }

    internal static unsafe bool Equal(MemorySaSink sink, int[] reference)
    {
        if (reference.Length != sink.Entries) return false;
        return new ReadOnlySpan<uint>(sink.Data, reference.Length).SequenceEqual(System.Runtime.InteropServices.MemoryMarshal.Cast<int, uint>(reference));
    }
}
