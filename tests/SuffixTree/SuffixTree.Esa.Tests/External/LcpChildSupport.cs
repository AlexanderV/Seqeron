using System.Globalization;
using System.Reflection;
using System.Runtime.InteropServices;
using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The in-RAM oracles of PH4 + PH5 for one text (gate G1 of work package 1b.7c): the <see cref="EsaIndex"/> of the product (SA-IS, the LCP kernel, the child-table scan), its compact LCP, the raw field arrays of
/// its two single-array child layouts (read through reflection: the fields are private) and the header statistics. Everything the production stage writes is compared with these.
/// </summary>
internal sealed class LcpChildCase
{
    private LcpChildCase(string text, byte[] symbols, int sigma, EsaIndex packed, EsaIndex tagged)
    {
        Text = text;
        Symbols = symbols;
        Sigma = sigma;
        Bits = TextOracle.BitsOf(Math.Max(1, sigma));
        N = text.Length;
        Index = packed;
        Packed = Raw(packed.Child.Packed);
        Tagged = Raw(tagged.Child.Tagged);
        int[] dense = packed.Lcp.ToArray();
        MaxLcp = dense.Length == 0 ? 0 : dense.Max();
    }

    internal string Text { get; }

    /// <summary>The dense symbols: the rank of each character among the distinct characters (increasing code unit order).</summary>
    internal byte[] Symbols { get; }

    internal int Sigma { get; }

    internal int Bits { get; }

    internal int N { get; }

    internal EsaIndex Index { get; }

    internal int[] Sa => Index.SuffixArray;

    internal CompactLcp Lcp => Index.Lcp;

    internal EsaStatistics Statistics => Index.Statistics;

    /// <summary>The raw array of the in-RAM packed layout (the codec <see cref="ChildCodec.Packed"/>).</summary>
    internal uint[] Packed { get; }

    /// <summary>The raw array of the in-RAM tagged layout (the codec <see cref="ChildCodec.Tagged"/>).</summary>
    internal uint[] Tagged { get; }

    internal int MaxLcp { get; }

    internal uint[] Child(ChildCodec codec) => codec == ChildCodec.Tagged ? Tagged : Packed;

    internal static LcpChildCase Of(string text)
    {
        char[] letters = [.. text.Distinct().OrderBy(c => c)];
        var symbols = new byte[text.Length];
        for (int i = 0; i < text.Length; i++) symbols[i] = (byte)Array.BinarySearch(letters, text[i]);
        EsaIndex packed = EsaIndex.Build(text, ChildLayoutKind.Packed);
        return new LcpChildCase(text, symbols, letters.Length, packed, packed.WithLayout(ChildLayoutKind.Tagged));
    }

    /// <summary>The expected base array of the compact LCP, as the bytes of the file.</summary>
    internal byte[] LcpBaseBytes() => Lcp.Kind switch
    {
        LcpKind.U8 => Lcp.BaseU8!,
        LcpKind.U16 => MemoryMarshal.AsBytes(Lcp.BaseU16.AsSpan()).ToArray(),
        _ => MemoryMarshal.AsBytes(Lcp.BaseU32.AsSpan()).ToArray(),
    };

    /// <summary>The expected overflow list: pairs <c>(rank, value)</c> of two uints.</summary>
    internal byte[] OverflowBytes()
    {
        var pairs = new uint[Lcp.OverflowCount * 2];
        for (int i = 0; i < Lcp.OverflowCount; i++)
        {
            pairs[2 * i] = (uint)Lcp.OverflowRanks[i];
            pairs[(2 * i) + 1] = (uint)Lcp.OverflowValues[i];
        }
        return MemoryMarshal.AsBytes(pairs.AsSpan()).ToArray();
    }

    /// <summary>The dense LCP values as uints (rank 0 included).</summary>
    internal uint[] DenseLcp() => Lcp.ToArray().Select(v => (uint)v).ToArray();

    /// <summary>The text of the statistics file, written out from the definition in the test (the format of the prototype's <c>.stats</c>).</summary>
    internal string ExpectedStats(int width) => string.Create(CultureInfo.InvariantCulture,
        $"nodecount={Statistics.NodeCount}\ninternal={Statistics.InternalCount}\nlrsdepth={Statistics.LrsDepth}\nlrsrank={Statistics.LrsRank}\nn={N}\nsigma={Sigma}\nlcpwidth={width}\noverflow={Lcp.OverflowCount}\nmaxlcp={MaxLcp}\n");

    private static uint[] Raw<TLayout>(TLayout layout)
        where TLayout : struct
    {
        object field = typeof(TLayout).GetField("_cld", BindingFlags.NonPublic | BindingFlags.Instance)!.GetValue(layout)!;
        return field switch
        {
            uint[] u => u,
            int[] i => MemoryMarshal.Cast<int, uint>(i).ToArray(),
            _ => throw new InvalidOperationException("unexpected field type"),
        };
    }
}

/// <summary>The open file descriptors of this process (Linux: <c>/proc/self/fd</c>): lets a test see that a build closed every file it opened.</summary>
internal static class ProcFds
{
    /// <summary>The paths the descriptors point to (a deleted file without the kernel's " (deleted)" suffix), or null where the platform has no <c>/proc/self/fd</c>.</summary>
    internal static List<string>? Open()
    {
        const string dir = "/proc/self/fd";
        if (!OperatingSystem.IsLinux() || !Directory.Exists(dir)) return null;
        var result = new List<string>();
        foreach (string entry in Directory.GetFileSystemEntries(dir))
        {
            try
            {
                string? target = new FileInfo(entry).LinkTarget;
                if (target is null) continue;
                result.Add(target.EndsWith(" (deleted)", StringComparison.Ordinal) ? target[..^" (deleted)".Length] : target);
            }
            catch (IOException)
            {
                // the descriptor of the directory listing itself, or one that closed in the meantime
            }
        }
        return result;
    }

    /// <summary>Fails when a descriptor points to one of <paramref name="paths"/> or below one of <paramref name="directories"/>.</summary>
    internal static void AssertClosed(string label, IEnumerable<string?> paths, IEnumerable<string?> directories)
    {
        List<string>? open = Open();
        if (open is null) return;
        string[] files = [.. paths.Where(p => p is not null).Select(p => Path.GetFullPath(p!))];
        string[] dirs = [.. directories.Where(d => d is not null).Select(d => Path.GetFullPath(d!) + "/")];
        foreach (string fd in open)
        {
            if (files.Contains(fd) || dirs.Any(d => fd.StartsWith(d, StringComparison.Ordinal))) Assert.Fail($"{label}: a file descriptor is still open on {fd}");
        }
    }
}

/// <summary>One run of the production stage on a case, and the comparison of everything it wrote with the oracles.</summary>
internal static class LcpChildRun
{
    /// <summary>The files and the result of one run.</summary>
    internal sealed record Outcome(LcpChildResult Result, LcpChildOutputs Files, string OutDir, string ScratchDir, string? SaPath = null);

    internal static int Width(LcpKind kind) => kind switch { LcpKind.U8 => 1, LcpKind.U16 => 2, _ => 4 };

    /// <summary>Tiny parameters that rotate with <paramref name="v"/>: R = 3..50, tiny buffers, blocks of 2..6 entries, windows of 16..64 bytes, so that every chunk, spill and window path runs on small texts.</summary>
    internal static LcpChildTuning TinyTuning(int v) => new()
    {
        LcpChunk = 3 + ((v * 11) % 48),
        LcpWorkers = 1 + (v % 3),
        LcpReadBytes = 16 + (8 * (v % 8)),
        LcpTupleBytes = 48 + (4 * (v % 5)),
        WindowBytes = 64 + (32 * (v % 6)),
        BlockEntries = 2 + (v % 5),
        ChildWindowBytes = 16 << (v % 3),
        ChildSlots = 1 + (v % 3),
        SpillFlushBytes = 64,
    };

    internal static LcpChildOptions Options(string scratch, LcpChildTuning tuning, int threads = 1, ChildCodec codec = ChildCodec.Packed, int? sigma = null) => new()
    {
        ScratchDir = scratch,
        Threads = threads,
        Codec = codec,
        Tuning = tuning,
        Sigma = sigma,
        ArenaGuards = ArenaOptions.Guarded(0xCD), // every buffer that is read before it is written shows as 0xCDCDCDCD
    };

    /// <summary>Writes the suffix array (as uints) and the packed text of <paramref name="c"/> into <paramref name="dir"/> and returns the paths.</summary>
    internal static (string Sa, string Text) WriteInputs(TempDir dir, LcpChildCase c, string tag = "")
    {
        string sa = dir.File($"sa{tag}.bin");
        FileBits.WriteUInts(sa, c.Sa.Select(v => (uint)v).ToArray());
        string text = TextOracle.WriteFile(dir, $"text{tag}.pk", c.Symbols, c.Bits);
        return (sa, text);
    }

    /// <summary>Runs the stage with the outputs in <c>dir/out</c> (the scratch directory is <c>dir/scratch</c>, it does not exist before and must not exist afterwards).</summary>
    internal static Outcome Run(TempDir dir, LcpChildCase c, LcpChildOptions options, bool stats = true, string tag = "", Action<MappedTextReader>? configure = null)
    {
        (string sa, string textPath) = WriteInputs(dir, c, tag);
        return RunOn(dir, sa, textPath, c.N, c.Bits, options, stats, tag, configure);
    }

    internal static Outcome RunOn(TempDir dir, string saPath, string textPath, long n, int bits, LcpChildOptions options, bool stats = true, string tag = "", Action<MappedTextReader>? configure = null)
    {
        string outDir = dir.Sub("out" + tag);
        Directory.CreateDirectory(outDir);
        var files = new LcpChildOutputs(Path.Combine(outDir, "o.lcp"), Path.Combine(outDir, "o.ovf"), Path.Combine(outDir, "o.cld"), stats ? Path.Combine(outDir, "o.stats") : null);
        using MappedTextReader reader = MappedTextReader.Open(textPath, n, bits);
        configure?.Invoke(reader);
        LcpChildResult result = LcpChildBuilder.Build(saPath, reader, files, options);
        return new Outcome(result, files, outDir, options.ScratchDir, saPath);
    }

    /// <summary>Fails with a description of the first difference between everything the run wrote and the in-RAM structures; also checks that nothing is left behind.</summary>
    internal static void Compare(LcpChildCase c, Outcome o, ChildCodec codec, string label)
    {
        LcpChildResult r = o.Result;
        Assert.That(r.Lcp.Kind, Is.EqualTo(c.Lcp.Kind), $"{label}: LCP width");
        int width = Width(c.Lcp.Kind);
        byte[] baseBytes = File.ReadAllBytes(o.Files.Lcp);
        Assert.That(baseBytes, Has.Length.EqualTo((long)width * (c.N + 1)), $"{label}: base LCP bytes");
        Assert.That(baseBytes.AsSpan().SequenceEqual(c.LcpBaseBytes()), Is.True, $"{label}: base LCP ({c.Lcp.Kind})");
        Assert.That(File.ReadAllBytes(o.Files.Overflow).AsSpan().SequenceEqual(c.OverflowBytes()), Is.True, $"{label}: overflow list");
        uint[] child = FileBits.ReadUInts(o.Files.Child);
        uint[] expected = c.Child(codec);
        Assert.That(child, Has.Length.EqualTo(c.N + 1), $"{label}: CHILD entries");
        for (int k = 0; k < child.Length; k++)
        {
            if (child[k] != expected[k]) Assert.Fail($"{label}: CHILD[{k}] = 0x{child[k]:X8}, expected 0x{expected[k]:X8} ({codec})");
        }
        Assert.That(r.Scan.NodeCount, Is.EqualTo(c.Statistics.NodeCount), $"{label}: NODE_COUNT");
        Assert.That(r.Scan.InternalCount, Is.EqualTo(c.Statistics.InternalCount), $"{label}: internal nodes");
        Assert.That((int)r.Scan.LrsDepth, Is.EqualTo(c.Statistics.LrsDepth), $"{label}: LRS_DEPTH");
        Assert.That((int)r.Scan.LrsRank, Is.EqualTo(c.Statistics.LrsRank), $"{label}: LRS_RANK");
        Assert.That(r.Scan.MaxLcp, Is.EqualTo(c.MaxLcp), $"{label}: max LCP");
        Assert.That(r.Lcp.OverflowCount, Is.EqualTo(c.Lcp.OverflowCount), $"{label}: overflow count");
        if (o.Files.Stats is not null) Assert.That(File.ReadAllText(o.Files.Stats), Is.EqualTo(c.ExpectedStats(width)), $"{label}: statistics file");
        AssertNothingLeft(o, label);
    }

    internal static void AssertNothingLeft(Outcome o, string label)
    {
        Assert.That(Directory.Exists(o.ScratchDir), Is.False, $"{label}: the scratch directory is left behind");
        string[] expected = new[] { o.Files.Lcp, o.Files.Overflow, o.Files.Child, o.Files.Stats }.Where(p => p is not null).Select(p => Path.GetFileName(p!)).Order().ToArray();
        Assert.That(Directory.GetFiles(o.OutDir).Select(Path.GetFileName).Order(), Is.EqualTo(expected), $"{label}: the output directory holds exactly the outputs");
        foreach (string? path in new[] { o.Files.Lcp, o.Files.Overflow, o.Files.Child })
            Assert.That(ProcMaps.Contains(path!) ?? false, Is.False, $"{label}: {Path.GetFileName(path)} is still mapped");
        ProcFds.AssertClosed(label, [o.Files.Lcp, o.Files.Overflow, o.Files.Child, o.SaPath], [o.ScratchDir]);
    }

    /// <summary>A run that must end with an exception: nothing is left behind (no output, no scratch).</summary>
    internal static void AssertFailsCleanly(TempDir dir, string saPath, string textPath, long n, int bits, LcpChildOptions options, Type exception, string label, bool stats = true)
    {
        string outDir = dir.Sub("out-fail");
        Directory.CreateDirectory(outDir);
        var files = new LcpChildOutputs(Path.Combine(outDir, "o.lcp"), Path.Combine(outDir, "o.ovf"), Path.Combine(outDir, "o.cld"), stats ? Path.Combine(outDir, "o.stats") : null);
        using MappedTextReader reader = MappedTextReader.Open(textPath, n, bits);
        Exception? thrown = null;
        try
        {
            LcpChildBuilder.Build(saPath, reader, files, options);
        }
        catch (Exception e)
        {
            thrown = e;
        }
        Assert.That(thrown, Is.Not.Null, $"{label}: the build must fail");
        Assert.That(thrown, Is.InstanceOf(exception), $"{label}: {thrown}");
        Assert.That(Directory.GetFileSystemEntries(outDir), Is.Empty, $"{label}: partial outputs are left behind");
        Assert.That(Directory.Exists(options.ScratchDir), Is.False, $"{label}: scratch is left behind");
        ProcFds.AssertClosed(label, [saPath], [options.ScratchDir, outDir]);
        Directory.Delete(outDir);
    }

    /// <summary>SHA-256 hex of the files of an outcome, in the order lcp, ovf, cld, stats.</summary>
    internal static string[] Hashes(Outcome o)
        => new[] { o.Files.Lcp, o.Files.Overflow, o.Files.Child, o.Files.Stats }.Where(p => p is not null)
            .Select(p => Convert.ToHexString(System.Security.Cryptography.SHA256.HashData(File.ReadAllBytes(p!)))).ToArray();

    /// <summary>The families of the contract tests, each as the text of n characters.</summary>
    internal static IEnumerable<(string Name, string Text)> Families(int n)
    {
        yield return ("a^n", TestSupport.StringFamilies.AllSame(n));
        yield return ("fibonacci", TestSupport.StringFamilies.Fibonacci(n));
        yield return ("thue-morse", TestSupport.StringFamilies.ThueMorse(n));
        yield return ("random-dna", TestSupport.StringFamilies.RandomDna(n, 17));
        yield return ("random-protein", TestSupport.StringFamilies.RandomProtein(n, 18));
        yield return ("repeat-rich", RepeatRich(n, 19));
        yield return ("n-runs", NRuns(n, 20));
    }

    /// <summary>Random DNA in which blocks of the text are copied with a few substitutions (like <c>rep_10m</c>).</summary>
    internal static string RepeatRich(int n, ulong seed)
    {
        char[] t = TestSupport.StringFamilies.RandomDna(n, seed).ToCharArray();
        var rng = new TestSupport.SplitMix64(seed + 1000);
        for (int copy = 0; copy < 6 && n > 40; copy++)
        {
            int len = Math.Max(8, n / 12);
            int from = rng.NextInt(n - len), to = rng.NextInt(n - len);
            for (int k = 0; k < len; k++) t[to + k] = t[from + k];
            for (int s = 0; s < len / 50; s++) t[to + rng.NextInt(len)] = "ACGT"[rng.NextInt(4)];
        }
        return new string(t);
    }

    /// <summary>Random DNA with runs of N (a 4-letter text with a long run of one more letter).</summary>
    internal static string NRuns(int n, ulong seed)
    {
        char[] t = TestSupport.StringFamilies.RandomDna(n, seed).ToCharArray();
        var rng = new TestSupport.SplitMix64(seed + 2000);
        for (int run = 0; run < 3 && n > 30; run++)
        {
            int len = Math.Max(4, n / 15);
            int at = rng.NextInt(n - len);
            for (int k = 0; k < len; k++) t[at + k] = 'N';
        }
        return new string(t);
    }
}
