using System.Diagnostics;
using System.Globalization;
using Microsoft.Win32.SafeHandles;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Bench;

/// <summary>
/// Work package 1b.9d: the production container build (<c>SesaBuilder.Build</c>, default parameters unless an option says otherwise) as a command, for the two measurements of the package. <c>sesa build</c>: one build in
/// this process, all the figures of the result and of the operating system as rows (the page-cache-limited runs of <c>analysis/em/run_1b9d_cache.sh</c> start it inside a cgroup). <c>sesa ab</c>: the same build several
/// times in this one process, alternating the two ways of treating the growth of the container behind its reserved length, the sparse growth that the builder had before this package (<c>A</c>: the range that PH4 writes the
/// LCP into has no blocks) and the reserved growth (<c>B</c>: zeros are written into it first), so that the difference is measured inside one process, where two percent can be told from noise.
/// </summary>
internal static class SesaCommands
{
    private static readonly CultureInfo Inv = CultureInfo.InvariantCulture;

    private static double Seconds(long start) => (Stopwatch.GetTimestamp() - start) / (double)Stopwatch.Frequency;

    /// <summary>The system's reservation, but the growth is sparse: the behaviour before work package 1b.9d (variant A of the A/B).</summary>
    private sealed class SparseTailSource : IReservationSource
    {
        public bool CanPreallocate => SystemReservationSource.Instance.CanPreallocate;

        public FileStream CreatePreallocated(string path, long length) => SystemReservationSource.Instance.CreatePreallocated(path, length);

        public long AvailableFreeBytes(string directory) => SystemReservationSource.Instance.AvailableFreeBytes(directory);

        public void ReserveRange(SafeFileHandle file, long offset, long length, IoCounters? counters)
        {
            // nothing: the range is a hole until PH4 writes into it
        }
    }

    /// <summary>Experiment of the tail measurement: one zero byte at the start of every 4 KiB page of the range instead of all of it (the pages are dirty either way).</summary>
    private sealed class StrideSource : IReservationSource
    {
        public bool CanPreallocate => SystemReservationSource.Instance.CanPreallocate;

        public FileStream CreatePreallocated(string path, long length) => SystemReservationSource.Instance.CreatePreallocated(path, length);

        public long AvailableFreeBytes(string directory) => SystemReservationSource.Instance.AvailableFreeBytes(directory);

        public void ReserveRange(SafeFileHandle file, long offset, long length, IoCounters? counters)
        {
            ReadOnlySpan<byte> one = [0];
            for (long at = offset; at < offset + length; at += 4096)
            {
                RandomAccess.Write(file, one, at);
                counters?.AddReserved(1);
            }
        }
    }

    /// <summary>Experiment of the tail measurement: the zeros, then a flush of the file to the device before the mapped writes (which adds the device writes of the zeros).</summary>
    private sealed class SyncedZeroSource : IReservationSource
    {
        public bool CanPreallocate => SystemReservationSource.Instance.CanPreallocate;

        public FileStream CreatePreallocated(string path, long length) => SystemReservationSource.Instance.CreatePreallocated(path, length);

        public long AvailableFreeBytes(string directory) => SystemReservationSource.Instance.AvailableFreeBytes(directory);

        public void ReserveRange(SafeFileHandle file, long offset, long length, IoCounters? counters)
        {
            ZeroFill.Write(file, offset, length, counters);
            RandomAccess.FlushToDisk(file);
        }
    }

    /// <summary>The system's reservation with the seconds of the reservation of the growth (variant B of the A/B).</summary>
    private sealed class TimedSource : IReservationSource
    {
        internal double ReserveSeconds { get; private set; }

        public bool CanPreallocate => SystemReservationSource.Instance.CanPreallocate;

        public FileStream CreatePreallocated(string path, long length) => SystemReservationSource.Instance.CreatePreallocated(path, length);

        public long AvailableFreeBytes(string directory) => SystemReservationSource.Instance.AvailableFreeBytes(directory);

        public void ReserveRange(SafeFileHandle file, long offset, long length, IoCounters? counters)
        {
            long t = Stopwatch.GetTimestamp();
            ZeroFill.Write(file, offset, length, counters);
            ReserveSeconds += Seconds(t);
        }
    }

    private sealed record Config(string Input, string Work, string Label, int Threads, long ScaleN, bool Psi, bool Fsync, bool Sparse, bool Cold, bool Keep, string Device, string? RowsFile, int Pairs, bool MeasureRss, bool Lean, string? FillAt);

    internal static void Run(string[] args)
    {
        if (args.Length > 0 && args[0] == "tail")
        {
            Tail(args[1..]);
            return;
        }
        if (args.Length > 0 && args[0] == "resume")
        {
            ResumeCommands.Run(args[1..]);
            return;
        }
        if (args.Length > 0 && args[0] == "handoff")
        {
            HandoffCommands.Run(args[1..]);
            return;
        }
        if (args.Length == 0 || args[0] is not ("build" or "ab")) throw new ArgumentException("usage: sesa build|ab --input <text file> --work <directory> [options] | sesa tail --work <directory> --mb <n> [options] | sesa resume ... | sesa handoff env|open|queries ...");
        string mode = args[0];
        string? input = null, work = null, label = null, rows = null, fillAt = null;
        string device = "vda";
        int threads = 4, pairs = 3;
        long scaleN = 0;
        bool psi = true, fsync = true, sparse = false, cold = false, keep = false, rss = true, lean = false;
        for (int i = 1; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            switch (args[i])
            {
                case "--input": input = value; break;
                case "--work": work = value; break;
                case "--label": label = value; break;
                case "--threads": threads = int.Parse(value, Inv); break;
                case "--scale-n": scaleN = long.Parse(value, Inv); break;
                case "--psi": psi = bool.Parse(value); break;
                case "--fsync": fsync = bool.Parse(value); break;
                case "--sparse": sparse = bool.Parse(value); break;
                case "--cold": cold = bool.Parse(value); break;
                case "--keep": keep = bool.Parse(value); break;
                case "--device": device = value; break;
                case "--rows": rows = value; break;
                case "--pairs": pairs = int.Parse(value, Inv); break;
                case "--rss": rss = bool.Parse(value); break;
                case "--lean": lean = bool.Parse(value); break;
                case "--fill-at": fillAt = value; break;
                default: throw new ArgumentException($"unknown option '{args[i]}'");
            }
        }
        if (input is null || work is null) throw new ArgumentException("--input and --work are required");
        Directory.CreateDirectory(work);
        var c = new Config(Path.GetFullPath(input), Path.GetFullPath(work), label ?? Path.GetFileNameWithoutExtension(input), threads, scaleN, psi, fsync, sparse, cold, keep, device, rows, pairs, rss, lean, fillAt);
        if (mode == "build") Build(c);
        else Ab(c);
    }

    /// <summary>
    /// <c>sesa tail</c>: the growth of a container alone, the cost of its three ways, with and without the page cache. A file of 1 MiB is created (reserved) and grown by <c>--mb</c> MiB, then the whole growth is written
    /// through 1 MiB mapped windows (what PH4 does with the LCP), then the file is completed and synced. Variants: <c>sparse</c> (the growth is a hole: the behaviour before 1b.9d), <c>zero</c> (the growth is reserved
    /// by writing zeros: 1b.9d), <c>worst</c> (the file is created with the growth reserved by <c>fallocate</c> and the unused part would be cut off afterwards: the alternative of reserving the worst case up front; the
    /// length is not cut here, only the creation is timed, and the file is deleted at once, so the disk it holds is the transient cost). Rows: seconds of the growth, of the write, of the sync, the process's I/O and the device's.
    /// </summary>
    private static void Tail(string[] args)
    {
        string? work = null;
        string variant = "zero", label = "tail";
        int mb = 100;
        bool cold = false;
        string device = "vda";
        for (int i = 0; i < args.Length; i += 2)
        {
            string value = i + 1 < args.Length ? args[i + 1] : throw new ArgumentException($"option {args[i]} needs a value");
            switch (args[i])
            {
                case "--work": work = value; break;
                case "--mb": mb = int.Parse(value, Inv); break;
                case "--variant": variant = value; break;
                case "--cold": cold = bool.Parse(value); break;
                case "--device": device = value; break;
                case "--label": label = value; break;
                default: throw new ArgumentException($"unknown option '{args[i]}'");
            }
        }
        if (work is null) throw new ArgumentException("--work is required");
        Directory.CreateDirectory(work);
        string path = Path.Combine(Path.GetFullPath(work), "tail.bin");
        File.Delete(path);
        const long head = 1 << 20;
        long growth = (long)mb << 20;
        var counters = new IoCounters();
        if (cold)
        {
            Sync();
            DropCaches();
        }
        ProcIoSnapshot? p0 = ProcIoSnapshot.TryRead();
        (long r0, long w0) = Sectors(device);
        long t0 = Stopwatch.GetTimestamp();
        double growS, writeS, syncS;
        long t1;
        if (variant == "worst")
        {
            using MappedOutputFile big = MappedOutputFile.Create(path, head + growth);
            growS = Seconds(t0);
            t1 = Stopwatch.GetTimestamp();
            using (MappedWindowWriter w = big.CreateWriter(0, head + growth, 1 << 20, counters)) w.Append(new byte[head + growth]);
            big.Complete(false);
            writeS = Seconds(t1);
        }
        else
        {
            IReservationSource source = variant switch
            {
                "sparse" => new SparseTailSource(),
                "stride" => new StrideSource(),
                "zerosync" => new SyncedZeroSource(),
                _ => SystemReservationSource.Instance,
            };
            using MappedOutputFile file = MappedOutputFile.Create(path, head);
            file.Grow(head + growth, new MappedFileOptions { Source = source }, counters);
            growS = Seconds(t0);
            t1 = Stopwatch.GetTimestamp();
            using (MappedWindowWriter w = file.CreateWriter(head, growth, 1 << 20, counters)) w.Append(new byte[growth]);
            file.Complete(false);
            writeS = Seconds(t1);
        }
        long t2 = Stopwatch.GetTimestamp();
        ProcIoSnapshot? p1 = ProcIoSnapshot.TryRead();
        (long r1, long w1) = Sectors(device);
        Sync();
        syncS = Seconds(t2);
        (_, long w2) = Sectors(device);
        string tag = $"{label} {variant} {mb} MiB {(cold ? "cold" : "warm")}";
        void Row(string metric, double v, string unit, string notes = "") => Emit.Row("sesa-tail", metric, v, unit, $"{tag} {notes}".TrimEnd());
        Row("grow_s", growS, "s", "create or grow, with the reservation");
        Row("write_s", writeS, "s", "the mapped windows over the whole growth");
        Row("sync_s", syncS, "s", "sync after the write");
        Row("total_s", growS + writeS + syncS, "s");
        Emit.Row("sesa-tail", "reserved_bytes", counters.ReservedWriteBytes, "bytes", tag);
        if (p0 is { } a && p1 is { } b)
        {
            ProcIoSnapshot d = b.Since(a);
            Emit.Row("sesa-tail", "proc_wchar", d.WChar, "bytes", tag);
            Emit.Row("sesa-tail", "proc_write_bytes", d.WriteBytes, "bytes", $"{tag} pages dirtied");
        }
        Emit.Row("sesa-tail", "dev_read_bytes", (r1 - r0) * 512, "bytes", tag);
        Emit.Row("sesa-tail", "dev_written_bytes_before_sync", (w1 - w0) * 512, "bytes", tag);
        Emit.Row("sesa-tail", "dev_written_bytes", (w2 - w0) * 512, "bytes", $"{tag} until the end of the sync");
        File.Delete(path);
    }

    private static void Build(Config c)
    {
        Measured m;
        try
        {
            m = BuildOnce(c, c.Sparse ? new SparseTailSource() : null);
        }
        catch (IOException e)
        {
            // a build that cannot get the disk space ends here with its exception (a signal would end the process before this line): exit code 3
            Emit.Row("sesa", "io_exception", e.Message.Replace('\t', ' ').Replace('\n', ' '), "text", c.Label);
            Emit.Row("sesa", "left_behind", Directory.GetFileSystemEntries(c.Work).Count(x => Path.GetFileName(x) is not "ballast"), "entries", $"{c.Label} entries of the work directory besides the ballast");
            Environment.ExitCode = 3;
            return;
        }
        Report(c, c.Label, m);
    }

    private static void Ab(Config c)
    {
        var ratios = new List<double>();
        var ratios45 = new List<double>();
        var reserve = new List<double>();
        // pair 0 is the warm-up of the process (JIT, the first touch of every buffer): it is run and reported like the others but left out of the ratios
        for (int pair = 0; pair <= c.Pairs; pair++)
        {
            bool aFirst = pair % 2 == 1 || pair == 0;
            Measured? a = null, b = null;
            TimedSource timed = new();
            foreach (char which in aFirst ? "AB" : "BA")
            {
                string label = $"{c.Label}.{pair}.{which}";
                if (which == 'A')
                {
                    a = BuildOnce(c, new SparseTailSource());
                    Report(c, label, a);
                }
                else
                {
                    b = BuildOnce(c, timed);
                    Report(c, label, b);
                    Emit.Row("sesa", "reserve_s", timed.ReserveSeconds, "s", $"{label} the seconds inside ReserveRange (the zeros of the growth)");
                }
            }
            if (pair == 0) continue;
            ratios.Add(b!.Result.Seconds.Total / a!.Result.Seconds.Total);
            ratios45.Add(b.Result.Seconds.Ph45 / a.Result.Seconds.Ph45);
            reserve.Add(timed.ReserveSeconds);
            Emit.Row("sesa", "pair_ratio_total", ratios[^1], "B / A", $"{c.Label} pair {pair}: total {b.Result.Seconds.Total:F3} / {a.Result.Seconds.Total:F3} s, first {(aFirst ? "A" : "B")}");
            Emit.Row("sesa", "pair_ratio_ph45", ratios45[^1], "B / A", $"{c.Label} pair {pair}: PH4 + PH5 {b.Result.Seconds.Ph45:F3} / {a.Result.Seconds.Ph45:F3} s");
        }
        ratios.Sort();
        ratios45.Sort();
        reserve.Sort();
        Emit.Row("sesa", "ratio_total_median", ratios[ratios.Count / 2], "B / A", $"{c.Label}: min {ratios[0]:F4} max {ratios[^1]:F4} over {ratios.Count} pairs");
        Emit.Row("sesa", "ratio_ph45_median", ratios45[ratios45.Count / 2], "B / A", $"{c.Label}: min {ratios45[0]:F4} max {ratios45[^1]:F4}");
        Emit.Row("sesa", "reserve_s_median", reserve[reserve.Count / 2], "s", $"{c.Label}: min {reserve[0]:F4} max {reserve[^1]:F4}");
    }

    private sealed record Measured(SesaBuildResult Result, ProcIoSnapshot? Proc, long SectorsRead, long SectorsWritten, long SectorsWrittenAfterSync, double Wall, long Reserved);

    private static Measured BuildOnce(Config c, IReservationSource? source)
    {
        string target = Path.Combine(c.Work, "out.sesa");
        Cleanup(target);
        if (c.Cold)
        {
            Sync();
            DropCaches();
        }
        var counters = new IoCounters();
        var options = new SesaBuildOptions
        {
            Threads = c.Threads,
            ScaleN = c.ScaleN,
            IncludePsi = c.Psi,
            Fsync = c.Fsync,
            MeasureRss = c.MeasureRss,
            // the lean settings of the cache tests of E7 (section 8 of RESULTS-E7.md: --budget-mb 64 --ph5-budget-mb 24 --spill-kb 16 --tuple-kb 8 --b2-workers 2), which cut RssAnon by more than half
            ArenaCapBytes = c.Lean ? 64L << 20 : new SesaBuildOptions().ArenaCapBytes,
            ChildArenaCapBytes = c.Lean ? 24L << 20 : new SesaBuildOptions().ChildArenaCapBytes,
            SaTuning = c.Lean ? new SaBuildTuning { SpillBytes = 16 << 10 } : new SaBuildTuning(),
            LcpTuning = c.Lean ? new LcpChildTuning { LcpTupleBytes = 8 << 10, LcpWorkers = 2 } : new LcpChildTuning(),
            Counters = counters,
            Files = source is null ? null : new MappedFileOptions { Source = source },
            Hook = c.FillAt is null ? null : point =>
            {
                if (point == c.FillAt) FillVolume(c.Work);
            },
        };
        ProcIoSnapshot? p0 = ProcIoSnapshot.TryRead();
        (long r0, long w0) = Sectors(c.Device);
        long t0 = Stopwatch.GetTimestamp();
        SesaBuildResult r = SesaBuilder.Build(c.Input, target, options);
        double wall = Seconds(t0);
        ProcIoSnapshot? p1 = ProcIoSnapshot.TryRead();
        (long r1, long w1) = Sectors(c.Device);
        Sync();
        (_, long w2) = Sectors(c.Device);
        if (!c.Keep) Cleanup(target);
        return new Measured(r, p0 is { } a && p1 is { } b ? b.Since(a) : null, r1 - r0, w1 - w0, w2 - w0, wall, counters.ReservedWriteBytes);
    }

    private static void Report(Config c, string label, Measured m)
    {
        SesaBuildResult r = m.Result;
        void Row(string metric, double v, string unit, string notes = "") => Emit.Row("sesa", metric, v, unit, string.IsNullOrEmpty(notes) ? label : $"{label} {notes}");
        void RowL(string metric, long v, string unit, string notes = "") => Emit.Row("sesa", metric, v, unit, string.IsNullOrEmpty(notes) ? label : $"{label} {notes}");
        Row("total_s", r.Seconds.Total, "s", "SesaBuilder.Build, inside the call");
        Row("wall_s", m.Wall, "s", "around the call");
        Row("scan_s", r.Seconds.Scan, "s");
        Row("create_s", r.Seconds.Create, "s");
        Row("pack_s", r.Seconds.Pack, "s");
        Row("ph1_s", r.Seconds.Ph1, "s");
        Row("ph23_s", r.Seconds.Ph23, "s");
        Row("ph45_s", r.Seconds.Ph45, "s");
        Row("bwt_s", r.Seconds.Bwt, "s");
        Row("psi_s", r.Seconds.Psi, "s");
        Row("crc_s", r.Seconds.Crc, "s");
        Row("seal_s", r.Seconds.Seal, "s");
        RowL("n", r.N, "symbols", $"sigma {r.Sigma}, {r.Bits} bits");
        RowL("file_size", r.FileSize, "bytes", $"LCP {r.LcpKind}, {r.OverflowCount} overflow pairs");
        RowL("grown_bytes", r.GrownBytes, "bytes", $"growth {(r.GrowthReservation?.ToString() ?? "none")}, creation {r.Reservation}");
        RowL("counted_read", r.Io.ReadBytes, "bytes", "explicit reads of the build (IoCounters)");
        RowL("counted_write", r.Io.WriteBytes, "bytes", "explicit writes of the build (IoCounters)");
        RowL("counted_mapped_write", r.Io.MappedWriteBytes, "bytes", "writes through mapped windows (IoCounters)");
        RowL("counted_reserved_write", m.Reserved, "bytes", "zeros written to reserve the growth");
        RowL("scratch_peak", r.ScratchPeakBytes, "bytes");
        RowL("arena_peak", r.ArenaPeakBytes, "bytes");
        if (m.Proc is { } p)
        {
            RowL("proc_rchar", p.RChar, "bytes", "/proc/self/io over the call");
            RowL("proc_wchar", p.WChar, "bytes");
            RowL("proc_read_bytes", p.ReadBytes, "bytes", "block-layer reads charged to the process");
            RowL("proc_write_bytes", p.WriteBytes, "bytes", "pages dirtied (a page that is written again counts again)");
        }
        RowL("dev_sectors_read", m.SectorsRead, "sectors", $"/sys/block/{c.Device}/stat over the call (512 B)");
        RowL("dev_sectors_written", m.SectorsWritten, "sectors", "over the call");
        RowL("dev_sectors_written_synced", m.SectorsWrittenAfterSync, "sectors", "until the end of a sync after the call");
        if (r.RssAnonPeakMb is { } rss) Row("rssanon_peak_mb", rss, "MB", $"verdict {r.RssVerdict}");
        foreach ((string name, double? v) in new[] { ("ph0", r.Rss.Ph0), ("ph1", r.Rss.Ph1), ("ph23", r.Rss.Ph23), ("ph45", r.Rss.Ph45), ("ph6", r.Rss.Ph6), ("seal", r.Rss.Seal) })
        {
            if (v is { } x) Row("rssanon_" + name + "_mb", x, "MB");
        }
        if (c.RowsFile is not null) File.AppendAllText(c.RowsFile, $"#{label}\ttotal {r.Seconds.Total:F3}\tph45 {r.Seconds.Ph45:F3}\n");
    }

    internal static (long Read, long Written) Sectors(string device)
    {
        try
        {
            string[] f = File.ReadAllText($"/sys/block/{device}/stat").Split(' ', StringSplitOptions.RemoveEmptyEntries);
            return (long.Parse(f[2], Inv), long.Parse(f[6], Inv));
        }
        catch (Exception e) when (e is IOException or UnauthorizedAccessException or IndexOutOfRangeException or FormatException)
        {
            return (0, 0);
        }
    }

    /// <summary>
    /// The disk of a full-volume run is taken at the point <c>--fill-at</c> names (a hook of <see cref="SesaHooks"/>, e.g. <c>lcp-placed</c>: after the late pre-flight and the growth of the container): a file of the name
    /// <c>ballast</c> in the directory is written until the volume refuses the next byte, and a thread goes on trying for the rest of the process, so that every block that the build frees afterwards (it deletes its scratch
    /// files as it goes) is taken at once: that is what another process that writes without end does to the volume in the window of work package 1b.9d.
    /// </summary>
    private static void FillVolume(string directory)
    {
        SafeFileHandle handle = File.OpenHandle(Path.Combine(directory, "ballast"), FileMode.Create, FileAccess.Write, FileShare.None);
        long offset = 0;
        offset = TakeWhatIsFree(handle, offset);
        Console.Error.WriteLine($"filled {directory}: ballast of {offset} bytes");
        var thief = new Thread(() =>
        {
            while (true)
            {
                offset = TakeWhatIsFree(handle, offset);
                Thread.Sleep(1);
            }
        })
        { IsBackground = true, Name = "ballast" };
        thief.Start();
    }

    private static long TakeWhatIsFree(SafeFileHandle handle, long offset)
    {
        var chunk = new byte[1 << 20];
        foreach (int size in new[] { 1 << 20, 4096, 1 })
        {
            try
            {
                while (true)
                {
                    RandomAccess.Write(handle, chunk.AsSpan(0, size), offset);
                    offset += size;
                }
            }
            catch (IOException)
            {
                offset = RandomAccess.GetLength(handle);
            }
        }
        return offset;
    }

    internal static void Sync()
    {
        using System.Diagnostics.Process? p = System.Diagnostics.Process.Start(new ProcessStartInfo("/usr/bin/sync"));
        p?.WaitForExit();
    }

    internal static void DropCaches() => File.WriteAllText("/proc/sys/vm/drop_caches", "3");

    private static void Cleanup(string target)
    {
        File.Delete(target);
        File.Delete(target + ".tmp");
        if (Directory.Exists(target + ".tmp.d")) Directory.Delete(target + ".tmp.d", true);
    }
}
