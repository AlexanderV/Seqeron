using System.Diagnostics;
using System.Globalization;
using System.Reflection;
using System.Runtime.Loader;
using System.Security.Cryptography;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Bench;

/// <summary>
/// Work package 1b.9e: the resumable build (<c>SesaBuildOptions.Resume</c>) as commands. <c>sesa resume build</c>: one build in this process, in a work directory that is NOT cleaned (so that a build that was killed leaves its
/// temporary file and state for the next call), with <c>M</c> lines on stdout at every boundary and phase that a driver can kill the process at (<c>--stop-at</c> blocks the process at a marker until it is killed), and every
/// figure of the result as rows. <c>sesa resume ab</c>: the same build with Resume off and on, alternating in one process (the cost of recording the state, where 2 % can be told from noise). <c>sesa resume abhead</c>: the
/// build of the commit before the package (loaded from its own directory into a second load context) against the build of this commit with Resume off, alternating in one process: the default must not slow down.
/// The drivers are <c>analysis/em/run_1b9e_kill.py</c> and <c>run_1b9e_cost.sh</c>.
/// </summary>
internal static class ResumeCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private sealed record Config(string Input, string Work, string Label, int Threads, bool Psi, bool Fsync, bool Resume, bool Markers, string? StopAt, bool Hash, bool Cold, string Device, int Pairs, string? Old, bool Keep);

    internal static void Run(string[] args)
    {
        if (args.Length == 0 || args[0] is not ("build" or "ab" or "abhead")) throw new ArgumentException("usage: sesa resume build|ab|abhead --input <text file> --work <directory> [--threads 4] [--psi true] [--fsync true] [--resume true] [--markers false] [--stop-at <marker>] [--hash false] [--cold false] [--pairs 5] [--old <directory of the old build>] [--label <text>]");
        string mode = args[0];
        string? input = null, work = null, label = null, stopAt = null, old = null;
        string device = "vda";
        int threads = 4, pairs = 5;
        bool psi = true, fsync = true, resume = true, markers = false, hash = false, cold = false, keep = false;
        for (int i = 1; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            switch (args[i])
            {
                case "--input": input = value; break;
                case "--work": work = value; break;
                case "--label": label = value; break;
                case "--stop-at": stopAt = value; break;
                case "--old": old = value; break;
                case "--device": device = value; break;
                case "--threads": threads = int.Parse(value, Inv); break;
                case "--pairs": pairs = int.Parse(value, Inv); break;
                case "--psi": psi = bool.Parse(value); break;
                case "--fsync": fsync = bool.Parse(value); break;
                case "--resume": resume = bool.Parse(value); break;
                case "--markers": markers = bool.Parse(value); break;
                case "--hash": hash = bool.Parse(value); break;
                case "--cold": cold = bool.Parse(value); break;
                case "--keep": keep = bool.Parse(value); break;
                default: throw new ArgumentException($"unknown option '{args[i]}'");
            }
        }
        if (input is null || work is null) throw new ArgumentException("--input and --work are required");
        Directory.CreateDirectory(work);
        var c = new Config(Path.GetFullPath(input), Path.GetFullPath(work), label ?? Path.GetFileNameWithoutExtension(input), threads, psi, fsync, resume, markers, stopAt, hash, cold, device, pairs, old is null ? null : Path.GetFullPath(old), keep);
        switch (mode)
        {
            case "build": Build(c); break;
            case "ab": Ab(c); break;
            default: AbHead(c); break;
        }
    }

    private static double Seconds(long start) => (Stopwatch.GetTimestamp() - start) / (double)Stopwatch.Frequency;

    private static string Sha256Of(string path)
    {
        using var fs = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite, 1 << 20);
        return Convert.ToHexStringLower(SHA256.HashData(fs));
    }

    // ------------------------------------------------------------------ one build, as the unit of a kill campaign

    private static void Build(Config c)
    {
        string target = Path.Combine(c.Work, "out.sesa");
        long t0 = Stopwatch.GetTimestamp();
        var gate = new object();
        void Mark(string name)
        {
            if (!c.Markers && c.StopAt is null) return;
            lock (gate)
            {
                Console.WriteLine(string.Join('\t', "M", name, Seconds(t0).ToString("F4", Inv)));
                Console.Out.Flush();
            }
            if (name == c.StopAt) Thread.Sleep(Timeout.Infinite); // the driver kills the process here
        }
        var counters = new IoCounters();
        var options = new SesaBuildOptions
        {
            Threads = c.Threads,
            IncludePsi = c.Psi,
            Fsync = c.Fsync,
            Resume = c.Resume,
            Counters = counters,
            Hook = Mark,
            SaTrace = p =>
            {
                if (p is "ph1" or "ph23" or "rank-written") Mark(p);
            },
            LcpTrace = Mark,
        };
        if (c.Cold)
        {
            Sync();
            DropCaches();
        }
        (long r0, long w0) = Sectors(c.Device);
        long start = Stopwatch.GetTimestamp();
        SesaBuildResult r;
        try
        {
            r = SesaBuilder.Build(c.Input, target, options);
        }
        catch (IOException e)
        {
            Emit.Row("resume", "exception", $"{e.GetType().Name}: {e.Message.Replace('\t', ' ').Replace('\n', ' ')}", "text", c.Label);
            Environment.ExitCode = 3;
            return;
        }
        double wall = Seconds(start);
        (long r1, long w1) = Sectors(c.Device);
        Sync();
        (_, long w2) = Sectors(c.Device);
        void Row(string metric, double v, string unit, string notes = "") => Emit.Row("resume", metric, v, unit, string.IsNullOrEmpty(notes) ? c.Label : $"{c.Label} {notes}");
        void RowL(string metric, long v, string unit, string notes = "") => Emit.Row("resume", metric, v, unit, string.IsNullOrEmpty(notes) ? c.Label : $"{c.Label} {notes}");
        Row("total_s", r.Seconds.Total, "s", $"resume {c.Resume} fsync {c.Fsync} threads {c.Threads}");
        Row("wall_s", wall, "s");
        Row("scan_s", r.Seconds.Scan, "s");
        Row("pack_s", r.Seconds.Pack, "s");
        Row("ph1_s", r.Seconds.Ph1, "s");
        Row("ph23_s", r.Seconds.Ph23, "s");
        Row("ph45_s", r.Seconds.Ph45, "s");
        Row("bwt_s", r.Seconds.Bwt, "s");
        Row("psi_s", r.Seconds.Psi, "s");
        Row("crc_s", r.Seconds.Crc, "s");
        Row("seal_s", r.Seconds.Seal, "s");
        RowL("n", r.N, "symbols", $"sigma {r.Sigma}");
        Row("sort_partition_thread_s", r.Sa.Sort.PartitionSeconds, "s", "PH2: seconds of the partition steps, summed over the worker threads");
        Row("sort_leaf_thread_s", r.Sa.Sort.LeafSeconds, "s", "PH3: seconds of the leaf sorts (and their slice writes), summed over the worker threads");
        RowL("sort_buckets", r.Sa.Sort.Buckets, "buckets", $"{r.Sa.Sort.Leaves} leaves, {r.Sa.Sort.Partitions} partitions, depth {r.Sa.Sort.MaxDepth}");
        RowL("dcs_levels", r.Sa.Dcs.Levels.Count, "levels", $"names unique at level 0: {(r.Sa.Dcs.Levels.Count > 0 && r.Sa.Dcs.Levels[0].NamesUnique)}");
        RowL("file_size", r.FileSize, "bytes", $"LCP {r.LcpKind}");
        RowL("counted_read", r.Io.ReadBytes, "bytes");
        RowL("counted_write", r.Io.WriteBytes, "bytes");
        RowL("counted_mapped_write", r.Io.MappedWriteBytes, "bytes");
        RowL("counted_flushes", counters.Flushes, "calls", "fsync calls counted (header, magic, and every checkpoint with Resume)");
        RowL("dev_sectors_written", w1 - w0, "sectors", "over the call (512 B)");
        RowL("dev_sectors_written_synced", w2 - w0, "sectors", "until the end of a sync after the call");
        RowL("dev_sectors_read", r1 - r0, "sectors");
        if (r.RssAnonPeakMb is { } rss) Row("rssanon_peak_mb", rss, "MB", $"verdict {r.RssVerdict}");
        if (r.Resume is { } rp)
        {
            Emit.Row("resume", "skipped", string.Join(',', rp.Skipped), "stages", c.Label);
            Emit.Row("resume", "state_found", rp.StateFound ? "yes" : "no", "bool", $"{c.Label} {rp.Reason}");
            Row("verify_s", rp.VerifySeconds, "s", "state parsed and the outputs of the skipped stages verified (CRC-32C)");
            RowL("clear_scanned_bytes", rp.ClearScannedBytes, "bytes", "read and compared with zero: the container outside the ranges of the kept stages");
            RowL("clear_bytes", rp.ClearedBytes, "bytes", "of those, written as zeros (units of 4096 bytes that held a non-zero byte)");
            Row("clear_s", rp.ClearSeconds, "s", "the normalisation of the reopened container");
            RowL("state_writes", rp.StateWrites, "files");
            RowL("state_bytes", rp.StateBytes, "bytes");
            foreach (SesaCheckpointRow k in rp.Checkpoints)
            {
                Row($"ckpt_{k.Stage}_fsync_s", k.FsyncSeconds, "s");
                Row($"ckpt_{k.Stage}_crc_s", k.CrcSeconds, "s", $"{k.CrcBytes} bytes");
                Row($"ckpt_{k.Stage}_state_s", k.StateSeconds, "s");
            }
            Row("ckpt_total_s", rp.Checkpoints.Sum(k => k.FsyncSeconds + k.CrcSeconds + k.StateSeconds), "s", "recording all stages");
        }
        if (c.Hash) Emit.Row("resume", "sha256", Sha256Of(target), "hex", c.Label);
        if (c.Keep) return;
        File.Delete(target);
    }

    // ------------------------------------------------------------------ Resume off against on, in one process

    private static void Ab(Config c)
    {
        string target = Path.Combine(c.Work, "out.sesa");
        var ratios = new List<double>();
        var synced = new List<double>();
        for (int pair = 0; pair <= c.Pairs; pair++)
        {
            bool offFirst = pair % 2 == 1 || pair == 0;
            double[] seconds = new double[2];
            long[] sectors = new long[2];
            long[] flushes = new long[2];
            foreach (bool on in offFirst ? new[] { false, true } : new[] { true, false })
            {
                Cleanup(target);
                if (c.Cold)
                {
                    Sync();
                    DropCaches();
                }
                var counters = new IoCounters();
                (_, long w0) = Sectors(c.Device);
                SesaBuildResult r = SesaBuilder.Build(c.Input, target, new SesaBuildOptions { Threads = c.Threads, IncludePsi = c.Psi, Fsync = c.Fsync, Resume = on, Counters = counters, MeasureRss = false });
                Sync();
                (_, long w1) = Sectors(c.Device);
                int k = on ? 1 : 0;
                seconds[k] = r.Seconds.Total;
                sectors[k] = w1 - w0;
                flushes[k] = counters.Flushes;
                string label = $"{c.Label}.{pair}.{(on ? "on" : "off")}";
                Emit.Row("resume-ab", "total_s", r.Seconds.Total, "s", label);
                foreach ((string name, double v) in new[] { ("scan", r.Seconds.Scan), ("create", r.Seconds.Create), ("pack", r.Seconds.Pack), ("ph1", r.Seconds.Ph1), ("ph23", r.Seconds.Ph23), ("ph45", r.Seconds.Ph45), ("bwt", r.Seconds.Bwt), ("psi", r.Seconds.Psi), ("crc", r.Seconds.Crc), ("seal", r.Seconds.Seal) })
                {
                    Emit.Row("resume-ab", "phase_" + name + "_s", v, "s", label);
                }
                Emit.Row("resume-ab", "dev_sectors_written_synced", w1 - w0, "sectors", label);
                Emit.Row("resume-ab", "flushes", counters.Flushes, "calls", label);
                if (r.Resume is { } rp) Emit.Row("resume-ab", "ckpt_total_s", rp.Checkpoints.Sum(x => x.FsyncSeconds + x.CrcSeconds + x.StateSeconds), "s", label);
            }
            Cleanup(target);
            if (pair == 0) continue;
            ratios.Add(seconds[1] / seconds[0]);
            synced.Add(sectors[1] / (double)Math.Max(1, sectors[0]));
            Emit.Row("resume-ab", "pair_ratio_total", ratios[^1], "on / off", $"{c.Label} pair {pair}: {seconds[1]:F3} / {seconds[0]:F3} s, first {(offFirst ? "off" : "on")}");
            Emit.Row("resume-ab", "pair_ratio_device_sectors", synced[^1], "on / off", $"{c.Label} pair {pair}: {sectors[1]} / {sectors[0]} sectors written until the end of a sync");
        }
        ratios.Sort();
        synced.Sort();
        Emit.Row("resume-ab", "ratio_total_median", ratios[ratios.Count / 2], "on / off", $"{c.Label}: min {ratios[0]:F4} max {ratios[^1]:F4} over {ratios.Count} pairs, fsync {c.Fsync}");
        Emit.Row("resume-ab", "ratio_device_sectors_median", synced[synced.Count / 2], "on / off", $"{c.Label}: min {synced[0]:F4} max {synced[^1]:F4}");
    }

    // ------------------------------------------------------------------ the build before the package against this one, in one process

#pragma warning disable S3011 // the builder and its options of the commit before the package are internal: reaching them by reflection is what this instrument is for
    private static void AbHead(Config c)
    {
        if (c.Old is null) throw new ArgumentException("--old <directory with the SuffixTree.Esa.dll of the commit before the package> is required");
        string target = Path.Combine(c.Work, "out.sesa");
        var alc = new AssemblyLoadContext("before-1b9e");
        alc.Resolving += (ctx, name) =>
        {
            string path = Path.Combine(c.Old, name.Name + ".dll");
            return File.Exists(path) ? ctx.LoadFromAssemblyPath(path) : null;
        };
        Assembly oldEsa = alc.LoadFromAssemblyPath(Path.Combine(c.Old, "SuffixTree.Esa.dll"));
        Type builder = oldEsa.GetType("SuffixTree.Esa.External.SesaBuilder", true)!;
        Type optionsType = oldEsa.GetType("SuffixTree.Esa.External.SesaBuildOptions", true)!;
        MethodInfo build = builder.GetMethod("Build", BindingFlags.NonPublic | BindingFlags.Public | BindingFlags.Static)!;
        const BindingFlags any = BindingFlags.NonPublic | BindingFlags.Public | BindingFlags.Instance;
        double RunOld()
        {
            object o = Activator.CreateInstance(optionsType, true)!;
            optionsType.GetProperty("Threads", any)!.SetValue(o, c.Threads);
            optionsType.GetProperty("IncludePsi", any)!.SetValue(o, c.Psi);
            optionsType.GetProperty("Fsync", any)!.SetValue(o, c.Fsync);
            optionsType.GetProperty("MeasureRss", any)!.SetValue(o, false);
            object r = build.Invoke(null, [c.Input, target, o])!;
            object seconds = r.GetType().GetProperty("Seconds", any)!.GetValue(r)!;
            return (double)seconds.GetType().GetProperty("Total", any)!.GetValue(seconds)!;
        }
        double RunNew()
        {
            SesaBuildResult r = SesaBuilder.Build(c.Input, target, new SesaBuildOptions { Threads = c.Threads, IncludePsi = c.Psi, Fsync = c.Fsync, MeasureRss = false });
            return r.Seconds.Total;
        }
        var ratios = new List<double>();
        for (int pair = 0; pair <= c.Pairs; pair++)
        {
            bool oldFirst = pair % 2 == 1 || pair == 0;
            double a = 0, b = 0;
            foreach (bool useOld in oldFirst ? new[] { true, false } : new[] { false, true })
            {
                Cleanup(target);
                if (c.Cold)
                {
                    Sync();
                    DropCaches();
                }
                double s = useOld ? RunOld() : RunNew();
                if (useOld) a = s;
                else b = s;
                Emit.Row("resume-abhead", "total_s", s, "s", $"{c.Label}.{pair}.{(useOld ? "before" : "after")}");
            }
            Cleanup(target);
            if (pair == 0) continue;
            ratios.Add(b / a);
            Emit.Row("resume-abhead", "pair_ratio_total", ratios[^1], "after / before", $"{c.Label} pair {pair}: {b:F3} / {a:F3} s, first {(oldFirst ? "before" : "after")}");
        }
        ratios.Sort();
        Emit.Row("resume-abhead", "ratio_total_median", ratios[ratios.Count / 2], "after / before", $"{c.Label}: min {ratios[0]:F4} max {ratios[^1]:F4} over {ratios.Count} pairs, Resume off");
    }

#pragma warning restore S3011

    private static (long Read, long Written) Sectors(string device) => SesaCommands.Sectors(device);

    private static void Sync() => SesaCommands.Sync();

    private static void DropCaches() => SesaCommands.DropCaches();

    private static void Cleanup(string target)
    {
        File.Delete(target);
        File.Delete(target + ".tmp");
        File.Delete(target + ".lock");
        if (Directory.Exists(target + ".tmp.d")) Directory.Delete(target + ".tmp.d", true);
    }
}
