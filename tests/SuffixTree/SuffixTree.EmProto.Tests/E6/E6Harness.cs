using System.Globalization;
using System.Security.Cryptography;
using System.Text;
using SuffixTree.EmProto.Pipeline;
using SuffixTree.EmProto.Ph5;
using SuffixTree.EmProto.Tests.Ph5;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.WholePipeline;

/// <summary>
/// The oracles and the harness of gate G1 of the prototype E6 (ADR-004 section 5): the whole pipeline (PH0 to PH5) on a text, every output compared with the in-RAM structures of the product
/// (<see cref="ChildReference"/>: suffix array of SA-IS, compact LCP of the kernel, packed child table of the product, statistics), the packed text with an independent bit-by-bit packer,
/// the metadata with direct counts, and the clean-up of the scratch directory.
/// </summary>
internal static class PipelineChecks
{
    /// <summary>Tiny parameters that rotate with <paramref name="v"/>: K = 3..50, F = 2..4, X = 7..64 (<see cref="Period"/>), blocks and windows of a few bytes, so that every multi-bucket, multi-level, multi-chunk and spill path runs on small texts.</summary>
    internal static EndToEndTuning Tiny(int v) => new()
    {
        BucketCapacity = 3 + ((v * 7) % 48),
        FanOut = 2 + (v % 3),
        Oversample = 4 + (v % 5),
        SpillBytes = 64 + (16 * (v % 4)),
        DcsBucketCapacity = 3 + ((v * 5) % 48),
        DcsBaseCase = 4 + ((v * 29) % 197),
        DcsChunkSlots = 16 + ((v * 3) % 60),
        DcsSortCapacity = 3 + ((v * 11) % 46),
        DcsRunSpillBytes = 16 + ((v * 13) % 200),
        DcsFanOut = 2 + (v % 5),
        DcsReservoir = 8 + ((v * 17) % 120),
        DcsIoBytes = 16 + ((v * 19) % 300),
        DcsMergeBudget = 8 + ((v * 23) % 2000),
        LcpChunk = 3 + ((v * 11) % 48),
        LcpWorkers = 1 + (v % 3),
        LcpReadBytes = 16 + (8 * (v % 8)),
        LcpTupleBytes = 48 + (4 * (v % 5)),
        BlockEntries = 2 + (v % 5),
        ChildWindowBytes = 16 << (v % 3),
        ChildSlots = 1 + (v % 3),
        SpillFlushBytes = 64,
        WindowBytes = 64 + (32 * (v % 6)),
        PrepareBlockSymbols = 8 << (v % 4),
    };

    internal static int Period(int v) => 7 + ((v * 13) % 58);

    internal static byte[] Latin1(string text) => Encoding.Latin1.GetBytes(text);

    /// <summary>Options of a tiny run in <paramref name="dir"/>: the input file is written here, outputs are <c>dir/out.*</c>, the scratch directory is <c>dir/scratch</c> (it does not exist yet and must not exist afterwards).</summary>
    internal static EndToEndOptions Options(string dir, string text, int variant, int threads = 1)
    {
        string input = Path.Combine(dir, "in.txt");
        File.WriteAllBytes(input, Latin1(text));
        return new EndToEndOptions
        {
            InputPath = input,
            OutPrefix = Path.Combine(dir, "out"),
            ScratchDir = Path.Combine(dir, "scratch"),
            Threads = threads,
            Period = Period(variant),
            ExactPeriod = true,
            Tuning = Tiny(variant),
            ArenaGuards = ArenaOptions.Guarded(0xCD),
            SyncOutputs = false,
        };
    }

    /// <summary>Dense rank of every character of <paramref name="text"/> among its distinct characters (increasing code unit order).</summary>
    internal static (byte[] Symbols, char[] Letters) Dense(string text)
    {
        char[] letters = [.. text.Distinct().OrderBy(c => c)];
        var symbols = new byte[text.Length];
        for (int i = 0; i < text.Length; i++) symbols[i] = (byte)Array.BinarySearch(letters, text[i]);
        return (symbols, letters);
    }

    /// <summary>The packed text written bit by bit (first bit of a byte is its most significant bit), 64 zero bytes of padding: independent of <c>PackedText.PackInto</c>.</summary>
    internal static byte[] ReferencePacked(byte[] symbols, int sigma)
    {
        int bits = 1;
        while ((1 << bits) < sigma) bits++;
        var bytes = new byte[(((long)symbols.Length * bits) + 7) / 8 + 64];
        for (int i = 0; i < symbols.Length; i++)
        {
            for (int b = 0; b < bits; b++)
            {
                long pos = ((long)i * bits) + b;
                if (((symbols[i] >> (bits - 1 - b)) & 1) != 0) bytes[pos / 8] |= (byte)(0x80 >> (int)(pos % 8));
            }
        }
        return bytes;
    }

    internal static Dictionary<string, string> ReadKeyValues(string path) => File.ReadAllLines(path).Where(l => l.Contains('=')).Select(l => l.Split('=', 2)).ToDictionary(p => p[0], p => p[1]);

    /// <summary>Runs the pipeline on <paramref name="text"/> with <paramref name="options"/> and compares every output and every side effect with the in-RAM structures.</summary>
    internal static EndToEndResult RunAndCompare(ChildReference r, EndToEndOptions options, string label, bool childTraversal = false)
    {
        string dir = Path.GetDirectoryName(options.OutPrefix)!;
        EndToEndResult res = EndToEnd.Run(options);
        EndToEndFiles f = res.Files;
        int n = r.N;
        (byte[] symbols, char[] letters) = Dense(r.Text);

        // the packed text
        Assert.That(File.ReadAllBytes(f.Packed), Is.EqualTo(ReferencePacked(symbols, letters.Length)), $"{label}: packed text");

        // the suffix array
        uint[] sa = ChildFiles.ReadUInts(f.Sa);
        int[] expectedSa = r.Index.SuffixArray;
        Assert.That(sa, Has.Length.EqualTo(n + 1), $"{label}: SA length");
        for (int i = 0; i <= n; i++)
        {
            if (sa[i] != (uint)expectedSa[i]) Assert.Fail($"{label}: SA[{i}] = {sa[i]}, expected {expectedSa[i]}");
        }

        // the compact LCP, byte for byte, in the width the kernel chose
        string refDir = Path.Combine(dir, "ref");
        Directory.CreateDirectory(refDir);
        (string refBase, string refOverflow) = r.WriteLcpFiles(refDir, "ref");
        Assert.That(File.ReadAllBytes(f.Lcp), Is.EqualTo(File.ReadAllBytes(refBase)), $"{label}: base LCP ({r.Lcp.Kind})");
        Assert.That(File.ReadAllBytes(f.Overflow), Is.EqualTo(File.ReadAllBytes(refOverflow)), $"{label}: overflow list");

        // the child table, entry by entry, in the packed codec of the product
        ChildChecks.CompareChild(r, f.Child, options.Codec, label);

        // the statistics
        Dictionary<string, string> st = ReadKeyValues(f.Stats);
        Assert.That(long.Parse(st["nodecount"], CultureInfo.InvariantCulture), Is.EqualTo(r.Statistics.NodeCount), $"{label}: NODE_COUNT");
        Assert.That(long.Parse(st["internal"], CultureInfo.InvariantCulture), Is.EqualTo(r.Statistics.InternalCount), $"{label}: internal nodes");
        Assert.That(long.Parse(st["lrsdepth"], CultureInfo.InvariantCulture), Is.EqualTo(r.Statistics.LrsDepth), $"{label}: LRS_DEPTH");
        Assert.That(long.Parse(st["lrsrank"], CultureInfo.InvariantCulture), Is.EqualTo(r.Statistics.LrsRank), $"{label}: LRS_RANK");

        // the metadata of PH0
        Dictionary<string, string> meta = ReadKeyValues(f.PackedMeta);
        Assert.That(long.Parse(meta["n"], CultureInfo.InvariantCulture), Is.EqualTo(n), $"{label}: meta n");
        Assert.That(int.Parse(meta["bits"], CultureInfo.InvariantCulture), Is.EqualTo(Math.Max(1, (int)Math.Ceiling(Math.Log2(Math.Max(1, letters.Length))))), $"{label}: meta bits");
        Assert.That(int.Parse(meta["sigma"], CultureInfo.InvariantCulture), Is.EqualTo(letters.Length), $"{label}: meta sigma");
        Assert.That(meta["letters_hex"], Is.EqualTo(Convert.ToHexString(letters.Select(c => (byte)c).ToArray())), $"{label}: meta letters");
        long[] counts = letters.Select(c => (long)r.Text.Count(x => x == c)).ToArray();
        Assert.That(meta["counts"], Is.EqualTo(string.Join(',', counts)), $"{label}: meta counts");
        Assert.That(meta["sha256"], Is.EqualTo(Convert.ToHexStringLower(SHA256.HashData(Latin1(r.Text)))), $"{label}: meta SHA-256");

        // final lengths, written out again
        int width = r.Lcp.Kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 };
        Assert.That(new FileInfo(f.Sa).Length, Is.EqualTo(4L * (n + 1)), $"{label}: SA bytes");
        Assert.That(new FileInfo(f.Lcp).Length, Is.EqualTo((long)width * (n + 1)), $"{label}: LCP bytes");
        Assert.That(new FileInfo(f.Overflow).Length, Is.EqualTo(8L * r.Lcp.OverflowCount), $"{label}: overflow bytes");
        Assert.That(new FileInfo(f.Child).Length, Is.EqualTo(4L * (n + 1)), $"{label}: CHILD bytes");

        // nothing is left behind: no scratch entry, no scratch directory, no tracked scratch byte, no live arena byte; and no stray output file
        Assert.That(Directory.Exists(options.ScratchDir), Is.False, $"{label}: scratch directory still exists");
        Assert.That(res.Clean, Is.True, $"{label}: scratch entries {res.ScratchEntriesLeft}, tracked scratch {res.ScratchLiveBytes} bytes, live arena {res.ArenaLiveBytes} bytes");
        Assert.That(Directory.GetFiles(dir, "out*").Select(Path.GetFileName).OrderBy(x => x), Is.EqualTo(new[] { "out.cld", "out.lcp", "out.ovf", "out.pk", "out.pk.meta", "out.sa", "out.stats" }), $"{label}: the output files");
        if (childTraversal) CompareTraversal(r, res, label);
        Directory.Delete(refDir, true);
        return res;
    }

    /// <summary>The child enumeration through the packed CHILD file equals the in-RAM index's, and the logical hash equals the Ukkonen tree's (tiny texts).</summary>
    internal static void CompareTraversal(ChildReference r, EndToEndResult res, string label)
    {
        var expected = new EventRecorder();
        r.Index.Traverse(expected);
        var view = new PackedView(ChildFiles.ReadUInts(res.Files.Child), r.Lcp);
        var got = new EventRecorder();
        ChildTraversal.Traverse(r.Text, r.Index.SuffixArray, r.Lcp, view, got);
        Assert.That(got.Events, Is.EqualTo(expected.Events), $"{label}: enumeration through CHILD");
        byte[] ukkonen = Persistent.SuffixTreeSerializer.CalculateLogicalHash(global::SuffixTree.SuffixTree.Build(r.Text));
        byte[] through = Persistent.SuffixTreeSerializer.CalculateLogicalHash(new WalkTree(r.Text, v => ChildTraversal.Traverse(r.Text, r.Index.SuffixArray, r.Lcp, view, v)));
        Assert.That(Convert.ToHexString(through), Is.EqualTo(Convert.ToHexString(ukkonen)), $"{label}: logical hash through CHILD against the Ukkonen tree");
    }

    /// <summary>SHA-256 of every output file of a run, by extension: the identity of two runs for the thread and K comparisons.</summary>
    internal static Dictionary<string, string> OutputHashes(string prefix)
        => new[] { "pk", "pk.meta", "sa", "lcp", "ovf", "cld", "stats" }.ToDictionary(e => e, e => Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(prefix + "." + e))));
}

/// <summary>Per-test scratch directory under the system temp directory.</summary>
public abstract class TempDirFixture
{
    protected string Dir { get; private set; } = string.Empty;

    /// <summary>
    /// Where the runs of a test live: a memory file system when there is one (<c>/dev/shm</c>; the tiny-parameter runs msync and fsync tens of small files each, which on a disk costs 50 to 100 ms per run), else the system temp directory;
    /// <c>EMPROTO_TEST_TMP</c> overrides. The tests that need a real file system (dirty pages of a mapping) use <see cref="Path.GetTempPath"/> themselves.
    /// </summary>
    internal static string TempRoot()
    {
        string? env = Environment.GetEnvironmentVariable("EMPROTO_TEST_TMP");
        if (!string.IsNullOrEmpty(env)) return env;
        return Directory.Exists("/dev/shm") ? "/dev/shm" : Path.GetTempPath();
    }

    [SetUp]
    public void SetUpDir()
    {
        Dir = Path.Combine(TempRoot(), "emproto-e6-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(Dir);
    }

    [TearDown]
    public void TearDownDir() => Directory.Delete(Dir, true);

    /// <summary>A fresh sub-directory (one run each).</summary>
    protected string Run(string name)
    {
        string d = Path.Combine(Dir, name);
        Directory.CreateDirectory(d);
        return d;
    }
}
