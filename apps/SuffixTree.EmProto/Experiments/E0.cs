using System.Diagnostics;
using System.Globalization;
using SuffixTree.Esa;
using SuffixTree.Esa.Bench;
using SuffixTree.EmProto.Gen;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Experiments;

/// <summary>
/// E0 of ADR-004 section 5: the measurement harness itself. Gates: the RssAnon sampler (10 ms) costs less than 2 % on a CPU-bound workload (paired, interleaved runs with an A/A calibration
/// in the same session), and the explicit I/O counters agree with <c>/proc/self/io</c> within 1 %. The <see cref="ArenaBudget"/> is exercised for enforcement and for release.
/// </summary>
internal static class E0
{
    internal static int Run(Options o, Report r)
    {
        string dir = o.Require("dir");
        Directory.CreateDirectory(dir);
        int pairs = o.Int("pairs", 15);
        bool ok = true;
        ok &= SamplerBasics(r);
        ok &= SamplerOverhead(r, pairs, o.Int("n", 10_000_000));
        ok &= IoCounters(r, dir, o.Int("io-mb", 256));
        ok &= Budget(r);
        r.Row("E0", "all", "gate", ok ? "PASS" : "FAIL", "", "sampler overhead < 2 %, counters within 1 % of /proc/self/io, budget enforced");
        return ok ? 0 : 1;
    }

    private static unsafe bool SamplerBasics(Report r)
    {
        if (!RssSampler.IsSupported)
        {
            r.Row("E0", "sampler", "supported", "no", "", "no /proc/self/status RssAnon on this platform");
            return false;
        }
        using var sampler = new RssSampler(10, keepSeries: true);
        sampler.Start();
        long beforeKb = Environment.ProcessId > 0 ? ReadAnonKb() : 0;
        const long mb = 200;
        byte* p = (byte*)System.Runtime.InteropServices.NativeMemory.Alloc((nuint)(mb << 20));
        for (long i = 0; i < mb << 20; i += 4096) p[i] = 1; // touch every page
        Thread.Sleep(120);
        System.Runtime.InteropServices.NativeMemory.Free(p);
        Thread.Sleep(60);
        RssStats s = sampler.Stop();
        double deltaMb = (s.MaxRssAnonKb - beforeKb) * 1024.0 / 1e6;
        bool ok = deltaMb is > 190 and < 230;
        r.Row("E0", "sampler", "samples", s.Samples, "n", "10 ms interval over about 0.2 s");
        r.Row("E0", "sampler", "touched_alloc_mb", mb * 1.048576, "MB", "200 MiB of touched native memory");
        r.Row("E0", "sampler", "rssanon_delta", deltaMb, "MB", ok ? "within 190..230 MB: PASS" : "FAIL");
        return ok;
    }

    private static long ReadAnonKb()
    {
        foreach (string line in File.ReadLines("/proc/self/status"))
        {
            if (!line.StartsWith("RssAnon:", StringComparison.Ordinal)) continue;
            string[] parts = line.Split([' ', '\t'], StringSplitOptions.RemoveEmptyEntries);
            return long.Parse(parts[1], CultureInfo.InvariantCulture);
        }
        return 0;
    }

    private static bool SamplerOverhead(Report r, int pairs, int n)
    {
        SymbolText text = Families.RandomDna(n, 7);
        var prepared = new PreparedText(text.Letters.ToCharArray(), text.Symbols, null);
        _ = SaisKernel.Build(prepared); // JIT warm-up, excluded
        _ = SaisKernel.Build(prepared);

        double Time(bool withSampler)
        {
            using var sampler = new RssSampler(10);
            if (withSampler) sampler.Start();
            long t0 = Stopwatch.GetTimestamp();
            _ = SaisKernel.Build(prepared);
            double s = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
            if (withSampler) sampler.Stop();
            return s;
        }

        var aaA = new List<double>();
        var aaB = new List<double>();
        for (int i = 0; i < pairs; i++)
        {
            (double a, double b) = i % 2 == 0 ? (Time(false), Time(false)) : Swap(Time(false), Time(false));
            aaA.Add(a);
            aaB.Add(b);
        }
        var abA = new List<double>();
        var abB = new List<double>();
        for (int i = 0; i < pairs; i++)
        {
            if (i % 2 == 0)
            {
                double a = Time(false);
                double b = Time(true);
                abA.Add(a);
                abB.Add(b);
            }
            else
            {
                double b = Time(true);
                double a = Time(false);
                abA.Add(a);
                abB.Add(b);
            }
        }
        List<double> aa = PairedStats.Ratios(aaA, aaB);
        List<double> ab = PairedStats.Ratios(abA, abB);
        (double lo, double hi) = PairedStats.NoiseBand(aa);
        (double ilo, double ihi) = PairedStats.MedianInterval(ab);
        double median = PairedStats.Median(ab);
        double overhead = median - 1;
        bool ok = overhead < 0.02;
        string case1 = $"sais-{n / 1_000_000}m";
        r.Row("E0", case1, "workload_median", PairedStats.Median(abA), "s", $"SaisKernel on {n} random DNA symbols, unsampled side of {pairs} pairs");
        r.Row("E0", case1, "aa_median_ratio", PairedStats.Median(aa), "ratio", $"A/A, {pairs} pairs; noise band [{lo:F4}, {hi:F4}]");
        r.Row("E0", case1, "ab_median_ratio", median, "ratio", $"sampled/unsampled, {pairs} interleaved pairs; 95% interval [{ilo:F4}, {ihi:F4}]");
        r.Row("E0", case1, "sampler_overhead", overhead * 100, "%", ok ? "< 2 %: PASS" : "FAIL");
        return ok;
    }

    private static (double, double) Swap(double a, double b) => (b, a);

    private static bool IoCounters(Report r, string dir, int mb)
    {
        string path = Path.Combine(dir, "io-test.bin");
        using var sampler = new RssSampler(10);
        sampler.Start();
        Harness.IoCounters.Reset();
        ProcIoSnapshot before = ProcIoSnapshot.Read();
        long total = (long)mb << 20;
        using (var f = new CountedFile(path, FileMode.Create, FileAccess.ReadWrite))
        {
            var chunk = new byte[1 << 20];
            new Random(1).NextBytes(chunk);
            for (long off = 0; off < total; off += chunk.Length) f.WriteAt(chunk, off);
            var small = new byte[4096];
            var rng = new SplitMix64(3);
            for (int i = 0; i < 8192; i++) f.WriteAt(small, (long)rng.NextInt((int)(total / 4096)) * 4096);
            var big = new byte[256 * 1024];
            for (long off = 0; off + big.Length <= total; off += big.Length) f.ReadExactAt(big, off);
            for (int i = 0; i < 16384; i++) f.ReadExactAt(small, (long)rng.NextInt((int)(total / 4096)) * 4096);
        }
        ProcIoSnapshot after = ProcIoSnapshot.Read();
        sampler.Stop();
        ProcIoSnapshot d = after.Since(before);
        long myRead = Harness.IoCounters.ReadBytes, myWrite = Harness.IoCounters.WriteBytes;
        double relRead = Math.Abs(myRead - d.RChar) / (double)d.RChar;
        double relWrite = Math.Abs(myWrite - d.WChar) / (double)d.WChar;
        bool ok = relRead < 0.01 && relWrite < 0.01;
        r.Row("E0", "io", "counted_write", myWrite, "bytes", $"{Harness.IoCounters.WriteCalls} calls");
        r.Row("E0", "io", "proc_wchar", d.WChar, "bytes", "delta of /proc/self/io");
        r.Row("E0", "io", "write_difference", relWrite * 100, "%", relWrite < 0.01 ? "< 1 %: PASS" : "FAIL");
        r.Row("E0", "io", "counted_read", myRead, "bytes", $"{Harness.IoCounters.ReadCalls} calls");
        r.Row("E0", "io", "proc_rchar", d.RChar, "bytes", "delta of /proc/self/io (includes the sampler reads of /proc/self/status)");
        r.Row("E0", "io", "read_difference", relRead * 100, "%", relRead < 0.01 ? "< 1 %: PASS" : "FAIL");
        r.Row("E0", "io", "proc_read_bytes", d.ReadBytes, "bytes", "block layer (page-cache misses), informational");
        r.Row("E0", "io", "proc_write_bytes", d.WriteBytes, "bytes", "block layer (writeback done so far), informational");
        File.Delete(path);
        return ok;
    }

    private static bool Budget(Report r)
    {
        var budget = new ArenaBudget(10 << 20);
        bool refused = false;
        using (var arena = new BudgetedArena(budget))
        {
            _ = arena.Allocate<int>(1 << 20, zero: true, "a"); // 4 MiB
            int mark = arena.Mark();
            _ = arena.Allocate<int>(1 << 20, zero: true, "b"); // 8 MiB live
            try
            {
                _ = arena.Allocate<int>(1 << 20, zero: true, "c"); // would be 12 MiB
            }
            catch (ArenaBudgetExceededException)
            {
                refused = true;
            }
            arena.ReleaseTo(mark);
            _ = arena.Allocate<int>(1 << 20, zero: true, "d"); // fits again
        }
        bool ok = refused && budget.LiveBytes == 0 && budget.PeakBytes == 8 << 20;
        r.Row("E0", "budget", "enforced", ok ? "yes" : "no", "", $"cap 10 MiB, refused={refused}, live after dispose={budget.LiveBytes}, peak={budget.PeakBytes}");
        return ok;
    }
}
