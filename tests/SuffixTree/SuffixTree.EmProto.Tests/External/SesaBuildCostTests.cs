using System.Diagnostics;
using System.Globalization;
using SuffixTree.EmProto.Dcs;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph6;
using SuffixTree.EmProto.Text;
using Prod = SuffixTree.Esa.External;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>
/// Gates G2 (identity at 10^8), G3 (cost) and G5 (reported costs) of work package 1b.7e: one run per process, driven by <c>docs/refactoring/suffix-tree-upgrade/analysis/em/run_1b7e_g2.sh</c> and
/// <c>run_1b7e_g3.sh</c> (cold or warm cache, processes one after the other, environment variables name the files). Not part of any run (<c>[Explicit]</c>). Variants: <c>container</c> (the whole build
/// <c>SesaBuilder.Build</c>, with or without PSI), <c>separate</c> (the separate production phases that the earlier packages measured, each into a file of its own: <c>PackedTextWriter.WriteFile</c>,
/// <c>SaBuilder.Build</c>, <c>LcpChildBuilder.Build</c> and, with PSI, <c>BwtPass</c> + <c>PsiStream</c>, then the checksum pass over the files: the sum of the phases), <c>ph6proto</c> (the prototype of E9,
/// code unchanged, on the suffix array and the packed text of E6) and <c>ph6prod</c> (the production PH6 on the same files, into a file of its own).
/// </summary>
[TestFixture]
[Explicit("gates G2 / G3 / G5 of 1b.7e: run by analysis/em/run_1b7e_g2.sh and run_1b7e_g3.sh, which set the G3_* variables")]
[Category("G3Cost")]
public unsafe class SesaBuildCostTests
{
    private static string Env(string name) => Environment.GetEnvironmentVariable(name) ?? throw new InvalidOperationException($"{name} is not set");

    private static string EnvOr(string name, string fallback) => Environment.GetEnvironmentVariable(name) ?? fallback;

    private static string F(double v) => v.ToString("F4", CultureInfo.InvariantCulture);

    private static string Mb(double? v) => v is { } x ? x.ToString("F1", CultureInfo.InvariantCulture) : "na";

    [Test]
    public void RunOnce()
    {
        string variant = Env("G3_VARIANT");
        string rows = Env("G3_ROWS");
        string label = Env("G3_LABEL");
        long scale = long.Parse(EnvOr("G3_SCALE", "3000000000"), CultureInfo.InvariantCulture);
        int threads = int.Parse(EnvOr("G3_THREADS", "4"), CultureInfo.InvariantCulture);
        bool psi = EnvOr("G3_PSI", "1") == "1";
        var rss = new Prod.RssSampler(new Prod.RssReader(), 10);
        long before = new Prod.RssReader().Read().RssAnonKb;
        rss.Start();
        string row = variant switch
        {
            "container" => RunContainer(scale, threads, psi),
            "separate" => RunSeparate(scale, threads, psi),
            "ph6proto" => RunPh6Proto(),
            "ph6prod" => RunPh6Prod(),
            "crc" => RunCrc(),
            _ => throw new InvalidOperationException("G3_VARIANT must be container, separate, ph6proto, ph6prod or crc"),
        };
        Prod.RssStats s = rss.Stop();
        File.AppendAllText(rows, string.Join('\t', label, variant, row, before, s.MaxRssAnonKb, s.Samples) + "\n");
        TestContext.Out.WriteLine(row);
        Assert.That(row, Is.Not.Empty);
    }

    private static string RunContainer(long scale, int threads, bool psi)
    {
        string input = Env("G3_INPUT"), target = Env("G3_TARGET");
        var options = new Prod.SesaBuildOptions
        {
            Threads = threads,
            ScaleN = scale,
            IncludePsi = psi,
            Fsync = EnvOr("G3_FSYNC", "1") == "1",
            MeasureRss = EnvOr("G3_SAMPLER", "1") == "1",
            Files = EnvOr("G3_RESERVATION", "auto") == "sparse" ? new Prod.MappedFileOptions { Reservation = Prod.ReservationMode.Sparse } : null,
        };
        Prod.SesaBuildResult r = Prod.SesaBuilder.Build(input, target, options);
        File.WriteAllText(Env("G3_ROWS") + ".crc." + Env("G3_LABEL"), string.Join('\n', r.CrcRows.Select(c => string.Join('\t', c.Kind, c.Bytes, F(c.Seconds), c.Hardware))) + "\n");
        Prod.SesaPhaseSeconds s = r.Seconds;
        Prod.SesaPhaseRss p = r.Rss;
        return string.Join('\t', F(s.Scan), F(s.Create), F(s.Pack), F(s.Ph1), F(s.Ph23), F(s.Ph45), F(s.Bwt), F(s.Psi), F(s.Crc), F(s.Seal), F(s.Total), r.Io.ReadBytes, r.Io.WriteBytes, r.Io.MappedWriteBytes, r.ScratchPeakBytes, r.ArenaPeakBytes, Mb(p.Ph0), Mb(p.Ph1), Mb(p.Ph23), Mb(p.Ph45), Mb(p.Ph6), Mb(p.Seal), Mb(r.RssAnonPeakMb), r.RssVerdict, r.FileSize, r.GrownBytes, r.Reservation, r.OverflowCount, r.LcpKind);
    }

    private static string RunSeparate(long scale, int threads, bool psi)
    {
        string input = Env("G3_INPUT"), work = Env("G3_WORK");
        Directory.CreateDirectory(work);
        string pk = Path.Combine(work, "sep.pk"), sa = Path.Combine(work, "sep.sa"), scratch = Path.Combine(work, "scratch");
        var counters = new Prod.IoCounters();
        long t0 = Stopwatch.GetTimestamp();
        Prod.PackedTextInfo info = Prod.PackedTextWriter.WriteFile(input, pk, new Prod.PackedTextWriterOptions { Counters = counters, Fsync = false });
        double ph0 = Seconds(t0);
        using Prod.MappedTextReader text = Prod.MappedTextReader.Open(pk, info.N, info.Bits);
        Prod.SaBuildResult saRes = Prod.SaBuilder.Build(text, sa, new Prod.SaBuildOptions { ScratchDir = scratch, Threads = threads, ScaleN = scale, Counters = counters });
        Prod.LcpChildResult lcp = Prod.LcpChildBuilder.Build(sa, text, new Prod.LcpChildOutputs(Path.Combine(work, "sep.lcp"), Path.Combine(work, "sep.ovf"), Path.Combine(work, "sep.cld")), new Prod.LcpChildOptions { ScratchDir = scratch, Threads = threads, ScaleN = scale, Counters = counters });
        double bwtS = 0, psiS = 0;
        if (psi)
        {
            string bwt = Path.Combine(work, "sep.bwt"), psiPath = Path.Combine(work, "sep.psi");
            t0 = Stopwatch.GetTimestamp();
            Prod.BwtResult b = Prod.BwtPass.Build(sa, 0, text, info.Sigma, bwt, new Prod.BwtOptions { Counters = counters });
            bwtS = Seconds(t0);
            t0 = Stopwatch.GetTimestamp();
            using (Prod.MappedOutputFile file = Prod.MappedOutputFile.Create(psiPath, 4 * (info.N + 1)))
            {
                Prod.PsiStream.Build(bwt, info.N, b.PrimaryRow, info.Counts, file, 0, new Prod.PsiOptions { Counters = counters });
                file.Complete(false);
            }
            psiS = Seconds(t0);
        }
        string[] files = [pk, sa, Path.Combine(work, "sep.lcp"), Path.Combine(work, "sep.ovf"), Path.Combine(work, "sep.cld"), .. psi ? new[] { Path.Combine(work, "sep.psi") } : []];
        t0 = Stopwatch.GetTimestamp();
        (uint Crc, Prod.SesaCrcRow Row)[] crc = Prod.SesaChecksums.Compute(Path.Combine(work, "sep.sa"), [(3, 0, new FileInfo(sa).Length)]); // one file for the structure; the others below
        double crcS = Seconds(t0);
        foreach (string f in files.Where(f => f != sa))
        {
            long t1 = Stopwatch.GetTimestamp();
            _ = Prod.SesaChecksums.Compute(f, [(0, 0, new FileInfo(f).Length)]);
            crcS += Seconds(t1);
        }
        _ = crc;
        double total = Seconds(t0) + ph0 + saRes.Ph1Seconds + saRes.Ph23Seconds + lcp.Seconds + bwtS + psiS;
        double copyS = 0;
        long copied = 0;
        if (EnvOr("G3_COPY", "0") == "1")
        {
            // what building into the sections avoids: the finished outputs are copied one after the other into one container-sized file (1 MiB buffers, fsync at the end)
            string container = Path.Combine(work, "sep.container");
            t0 = Stopwatch.GetTimestamp();
            using (var outFile = new FileStream(container, FileMode.Create, FileAccess.Write, FileShare.None, 1 << 20))
            {
                var buffer = new byte[1 << 20];
                foreach (string f in files)
                {
                    using var inFile = new FileStream(f, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 20);
                    int got;
                    while ((got = inFile.Read(buffer, 0, buffer.Length)) > 0)
                    {
                        outFile.Write(buffer, 0, got);
                        copied += got;
                    }
                }
                outFile.Flush(true);
            }
            copyS = Seconds(t0);
        }
        return string.Join('\t', F(ph0), F(saRes.Ph1Seconds), F(saRes.Ph23Seconds), F(lcp.Seconds), F(bwtS), F(psiS), F(crcS), F(total), counters.ReadBytes, counters.WriteBytes, counters.MappedWriteBytes, counters.ScratchPeakBytes, files.Sum(f => new FileInfo(f).Length), F(copyS), copied);
    }

    private static string RunPh6Proto()
    {
        string pk = Env("G3_PK"), sa = Env("G3_SA"), work = Env("G3_WORK");
        long n = long.Parse(Env("G3_N"), CultureInfo.InvariantCulture);
        int bits = int.Parse(Env("G3_BITS"), CultureInfo.InvariantCulture), sigma = int.Parse(Env("G3_SIGMA"), CultureInfo.InvariantCulture);
        int window = int.Parse(EnvOr("G3_WINDOW", "1048576"), CultureInfo.InvariantCulture);
        Directory.CreateDirectory(work);
        string bwt = Path.Combine(work, "p.bwt"), psi = Path.Combine(work, "p.psi");
        using var mapped = new MappedFile(pk, writable: false);
        using var text = new PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        IoCounters.Reset();
        MappedIo.Reset();
        long t0 = Stopwatch.GetTimestamp();
        BwtMeta meta = BwtPass.Build(sa, text, sigma, bwt, new BwtOptions { WindowBytes = window });
        double bwtS = Seconds(t0);
        long r1 = IoCounters.ReadBytes, w1 = IoCounters.WriteBytes, m1 = MappedIo.WriteBytes;
        t0 = Stopwatch.GetTimestamp();
        PsiStats ps = PsiStream.Build(bwt, meta, psi, new PsiOptions { WindowBytes = window });
        double psiS = Seconds(t0);
        return string.Join('\t', F(bwtS), F(psiS), r1, w1, m1, IoCounters.ReadBytes - r1, IoCounters.WriteBytes - w1, MappedIo.WriteBytes - m1, ps.Windows, meta.PrimaryRow);
    }

    private static string RunPh6Prod()
    {
        string pk = Env("G3_PK"), sa = Env("G3_SA"), work = Env("G3_WORK");
        long n = long.Parse(Env("G3_N"), CultureInfo.InvariantCulture);
        int bits = int.Parse(Env("G3_BITS"), CultureInfo.InvariantCulture), sigma = int.Parse(Env("G3_SIGMA"), CultureInfo.InvariantCulture);
        int window = int.Parse(EnvOr("G3_WINDOW", "1048576"), CultureInfo.InvariantCulture);
        long[] counts = [.. Env("G3_COUNTS").Split(',').Select(c => long.Parse(c, CultureInfo.InvariantCulture))];
        Directory.CreateDirectory(work);
        string bwt = Path.Combine(work, "q.bwt"), psi = Path.Combine(work, "q.psi");
        using Prod.MappedTextReader text = Prod.MappedTextReader.Open(pk, n, bits);
        var bwtCounters = new Prod.IoCounters();
        var psiCounters = new Prod.IoCounters();
        long t0 = Stopwatch.GetTimestamp();
        Prod.BwtResult r = Prod.BwtPass.Build(sa, 0, text, sigma, bwt, new Prod.BwtOptions { WindowBytes = window, Counters = bwtCounters });
        double bwtS = Seconds(t0);
        t0 = Stopwatch.GetTimestamp();
        Prod.PsiStats ps;
        using (Prod.MappedOutputFile file = Prod.MappedOutputFile.Create(psi, 4 * (n + 1)))
        {
            ps = Prod.PsiStream.Build(bwt, n, r.PrimaryRow, counts, file, 0, new Prod.PsiOptions { WindowBytes = window, Counters = psiCounters });
            file.Complete(false);
        }
        double psiS = Seconds(t0);
        return string.Join('\t', F(bwtS), F(psiS), bwtCounters.ReadBytes, bwtCounters.WriteBytes, bwtCounters.MappedWriteBytes, psiCounters.ReadBytes, psiCounters.WriteBytes, psiCounters.MappedWriteBytes, ps.Windows, r.PrimaryRow);
    }

    /// <summary>G5: the checksum pass over the sections of a finished container (warm cache), the hardware path and the managed table path, per section.</summary>
    private static string RunCrc()
    {
        string path = Env("G3_SESA");
        using var stream = new FileStream(path, FileMode.Open, FileAccess.Read, FileShare.ReadWrite);
        var page = new byte[4096];
        stream.ReadExactly(page);
        int count = BitConverter.ToInt32(page, 48);
        var rows = new List<string>();
        using var map = System.IO.MemoryMappedFiles.MemoryMappedFile.CreateFromFile(stream, null, 0, System.IO.MemoryMappedFiles.MemoryMappedFileAccess.Read, HandleInheritability.None, leaveOpen: true);
        for (int i = 0; i < count; i++)
        {
            int kind = BitConverter.ToInt32(page, 128 + (64 * i));
            long offset = BitConverter.ToInt64(page, 128 + (64 * i) + 8), length = BitConverter.ToInt64(page, 128 + (64 * i) + 16);
            if (length < (1 << 20)) continue; // the small sections are noise
            foreach (string pass in new[] { "hardware", "table", "hardware" })
            {
                uint crc = 0;
                long t0 = Stopwatch.GetTimestamp();
                for (long done = 0; done < length; done += 16 << 20)
                {
                    long take = Math.Min(16 << 20, length - done);
                    using System.IO.MemoryMappedFiles.MemoryMappedViewAccessor view = map.CreateViewAccessor(offset + done, take, System.IO.MemoryMappedFiles.MemoryMappedFileAccess.Read);
                    byte* p = null;
                    view.SafeMemoryMappedViewHandle.AcquirePointer(ref p);
                    crc = pass == "table" ? Prod.Crc32C.AppendTable(crc, new ReadOnlySpan<byte>(p + view.PointerOffset, (int)take)) : Prod.Crc32C.Append(crc, p + view.PointerOffset, take);
                    view.SafeMemoryMappedViewHandle.ReleasePointer();
                }
                double s = Seconds(t0);
                rows.Add(string.Join('\t', "crc", kind, length, pass, F(s), F(length / 1e9 / s), crc.ToString("x8", CultureInfo.InvariantCulture)));
            }
        }
        return string.Join('|', rows);
    }

    private static double Seconds(long t0) => (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
}
