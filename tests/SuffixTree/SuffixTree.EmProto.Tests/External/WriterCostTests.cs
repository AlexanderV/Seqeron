using System.Diagnostics;
using System.Globalization;
using Prod = SuffixTree.Esa.External;
using ProtoDcs = SuffixTree.EmProto.Dcs;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>
/// Gate G3 of work package 1b.7b, part 2 (the mapped writer against the prototype's, and the reservation of the disk space): one process writes a 400 MB output (4 bytes for each of 10^8 symbols) in the pattern of E9 (PH6, Psi):
/// <c>seq</c> is one cursor writing sequentially (the SA, LCP, CHILD writers), <c>psi</c> is one cursor per symbol of a BWT file, each writing its own contiguous range through its own window (the Psi writer); variants are
/// the prototype <c>MappedWindowWriter</c> (one file per cursor, it cannot share a file), the production <see cref="Prod.MappedOutputFile"/> sparse and with <c>fallocate</c>. The time of the creation and the total are measured here;
/// the device writes (block layer, after <c>sync</c>) and the cold cache are the script's (<c>analysis/em/run_1b7b_g3.sh</c>). Not part of any run (<c>[Explicit]</c>).
/// </summary>
[TestFixture]
[Explicit("gate G3 of 1b.7b: run by analysis/em/run_1b7b_g3.sh, which sets the G3W_* variables")]
[Category("G3Cost")]
public class WriterCostTests
{
    private static string Env(string name) => Environment.GetEnvironmentVariable(name) ?? throw new InvalidOperationException($"{name} is not set");

    [Test]
    public void RunOnce()
    {
        string variant = Env("G3W_VARIANT"); // proto | sparse | fallocate
        string pattern = Env("G3W_PATTERN"); // seq | psi
        string bwtPath = Env("G3W_BWT");
        string outDir = Env("G3W_OUT");
        int window = int.Parse(Environment.GetEnvironmentVariable("G3W_WINDOW") ?? "1048576", CultureInfo.InvariantCulture);
        string rows = Env("G3W_ROWS");
        string label = Env("G3W_LABEL");
        Directory.CreateDirectory(outDir);

        // the symbols of the BWT file (outside the timed region: the counts are what the cursors' ranges need)
        var cnt = new long[4];
        long length = new FileInfo(bwtPath).Length;
        var block = new byte[1 << 20];
        using (var f = new FileStream(bwtPath, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 16))
        {
            int got;
            while ((got = f.Read(block, 0, block.Length)) > 0)
            {
                for (int i = 0; i < got; i++) cnt[block[i] & 3]++; // the file stays in the page cache: only the output is cold
            }
        }
        var start = new long[4];
        for (int c = 1; c < 4; c++) start[c] = start[c - 1] + (4 * cnt[c - 1]);
        long total = 4 * length;
        var rss = new Prod.RssSampler(new Prod.RssReader(), 10);
        long before = new Prod.RssReader().Read().RssAnonKb;
        rss.Start();
        var counters = new Prod.IoCounters();
        var chunk = new byte[1 << 16];
        for (int i = 0; i < chunk.Length; i++) chunk[i] = (byte)((i * 131) + 7);
        long t0 = Stopwatch.GetTimestamp();
        double created;
        long writtenBytes = 0;
        if (variant == "proto")
        {
            if (pattern == "seq")
            {
                using var w = new ProtoDcs.MappedWindowWriter(Path.Combine(outDir, "out.bin"), total, window);
                created = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
                for (long at = 0; at < total; at += chunk.Length) w.Append(chunk.AsSpan(0, (int)Math.Min(chunk.Length, total - at)));
                writtenBytes = total;
            }
            else
            {
                var w = new ProtoDcs.MappedWindowWriter[4];
                for (int c = 0; c < 4; c++) w[c] = new ProtoDcs.MappedWindowWriter(Path.Combine(outDir, $"out{c}.bin"), 4 * cnt[c], window);
                created = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
                Span<byte> four = stackalloc byte[4];
                using var f = new FileStream(bwtPath, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 16);
                long at = 0;
                int got;
                while ((got = f.Read(block, 0, block.Length)) > 0)
                {
                    for (int i = 0; i < got; i++)
                    {
                        BitConverter.TryWriteBytes(four, (uint)(at + i));
                        w[block[i] & 3].Append(four);
                    }
                    at += got;
                }
                writtenBytes = total;
                foreach (ProtoDcs.MappedWindowWriter x in w) x.Dispose();
            }
        }
        else
        {
            var options = new Prod.MappedFileOptions { Reservation = variant == "fallocate" ? Prod.ReservationMode.Preallocate : Prod.ReservationMode.Sparse };
            using Prod.MappedOutputFile file = Prod.MappedOutputFile.Create(Path.Combine(outDir, "out.bin"), total, options);
            created = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
            if (pattern == "seq")
            {
                using Prod.MappedWindowWriter w = file.CreateWriter(window, counters);
                for (long at = 0; at < total; at += chunk.Length) w.Append(chunk.AsSpan(0, (int)Math.Min(chunk.Length, total - at)));
                writtenBytes = total;
            }
            else
            {
                var w = new Prod.MappedWindowWriter[4];
                for (int c = 0; c < 4; c++) w[c] = file.CreateWriter(start[c], 4 * cnt[c], window, counters);
                Span<byte> four = stackalloc byte[4];
                using var f = new FileStream(bwtPath, FileMode.Open, FileAccess.Read, FileShare.Read, 1 << 16);
                long at = 0;
                int got;
                while ((got = f.Read(block, 0, block.Length)) > 0)
                {
                    for (int i = 0; i < got; i++)
                    {
                        BitConverter.TryWriteBytes(four, (uint)(at + i));
                        w[block[i] & 3].Append(four);
                    }
                    at += got;
                }
                writtenBytes = total;
                foreach (Prod.MappedWindowWriter x in w) x.Dispose();
            }
            file.Complete(fsync: false);
        }
        double seconds = (Stopwatch.GetTimestamp() - t0) / (double)Stopwatch.Frequency;
        Prod.RssStats stats = rss.Stop();
        string row = string.Join('\t', label, variant, pattern, window.ToString(CultureInfo.InvariantCulture), created.ToString("F4", CultureInfo.InvariantCulture), seconds.ToString("F4", CultureInfo.InvariantCulture),
            before.ToString(CultureInfo.InvariantCulture), stats.MaxRssAnonKb.ToString(CultureInfo.InvariantCulture), writtenBytes.ToString(CultureInfo.InvariantCulture), counters.MappedWindows.ToString(CultureInfo.InvariantCulture));
        File.AppendAllText(rows, row + "\n");
        TestContext.Out.WriteLine(row);
        Assert.That(writtenBytes, Is.EqualTo(total));
    }
}
