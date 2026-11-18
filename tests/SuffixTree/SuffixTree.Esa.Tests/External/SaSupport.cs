using System.Runtime.InteropServices;
using System.Security.Cryptography;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The in-RAM oracles of PH1 and PH2 / PH3 for one text (gate G1 of work package 1b.7d): the suffix array from the definition (a sort of the suffixes with an ordinal span comparison) for the short texts and from the
/// kernel of the product (<see cref="SaisKernel"/>) above that, and the dense symbols that the packed text holds.
/// </summary>
internal sealed class SaCase
{
    private const int NaiveLimit = 3000;
    private uint[]? _sa;

    private SaCase(string text, byte[] symbols, int sigma)
    {
        Text = text;
        Symbols = symbols;
        Sigma = sigma;
        Bits = TextOracle.BitsOf(Math.Max(1, sigma));
        N = text.Length;
    }

    internal string Text { get; }

    /// <summary>The dense symbols: the rank of each character among the distinct characters (increasing code unit order).</summary>
    internal byte[] Symbols { get; }

    internal int Sigma { get; }

    internal int Bits { get; }

    internal int N { get; }

    /// <summary>The suffix array, <c>n + 1</c> entries, <c>SA[0] = n</c>.</summary>
    internal uint[] Sa => _sa ??= BuildSa();

    internal static SaCase Of(string text)
    {
        char[] letters = [.. text.Distinct().OrderBy(c => c)];
        var symbols = new byte[text.Length];
        for (int i = 0; i < text.Length; i++) symbols[i] = (byte)Array.BinarySearch(letters, text[i]);
        return new SaCase(text, symbols, letters.Length);
    }

    private uint[] BuildSa()
    {
        int[] sa = N <= NaiveLimit ? NaiveSuffixArray.Build(Text) : SaisKernel.Build(TextPreparation.Prepare(Text));
        return [.. sa.Select(v => (uint)v)];
    }
}

/// <summary>The rank file that PH1 must write, from the definition: the sample positions in slot order, the rank of every sample suffix among the sample suffixes by a sort with the ordinal comparison.</summary>
internal static class DcsOracle
{
    /// <summary>
    /// The expected contents of the rank file: one <c>uint</c> per slot (block-major, the unused slots of the last block zero). A sample suffix's value is its 1-based rank among the sample suffixes (the empty suffix at n, when it is
    /// a sample position, is the smallest) plus |D| when the X-prefixes of the sample are not all distinct (then the reduced string of the recursion holds the |D| separators, which take the first ranks).
    /// </summary>
    internal static uint[] ExpectedRanks(byte[] s, DifferenceCover cover, out bool namesUnique)
    {
        int n = s.Length, x = cover.Period, m = cover.Residues.Length;
        long slots = (((long)n / x) + 1) * m;
        var sample = new List<(int Pos, long Slot)>();
        for (long block = 0; block <= n / x; block++)
        {
            for (int c = 0; c < m; c++)
            {
                long p = (block * x) + cover.Residues[c];
                if (p <= n) sample.Add(((int)p, (block * m) + c));
            }
        }
        int[] order = [.. sample.Select(t => t.Pos)];
        Array.Sort(order, (a, b) => Compare(s, a, b));
        var rank = new Dictionary<int, uint>();
        for (int i = 0; i < order.Length; i++) rank[order[i]] = (uint)(i + 1);
        var prefixes = new HashSet<string>();
        namesUnique = true;
        foreach ((int pos, _) in sample)
        {
            if (n - pos < x) continue;
            if (!prefixes.Add(Convert.ToBase64String(s.AsSpan(pos, x)))) namesUnique = false;
        }
        var file = new uint[slots];
        foreach ((int pos, long slot) in sample) file[slot] = rank[pos] + (namesUnique ? 0u : (uint)m);
        return file;
    }

    /// <summary>The definition of the order of two suffixes (a position at n is the empty suffix).</summary>
    internal static int Compare(byte[] s, int a, int b) => Math.Sign(s.AsSpan(a).SequenceCompareTo(s.AsSpan(b)));
}

/// <summary>Tiny parameters that rotate with a variant number, so that every multi-bucket, multi-level, spill, merge and window path runs on texts of 10^1-10^4 symbols.</summary>
internal static class SaTiny
{
    /// <summary>The periods used with an exact (pruned or Singer) cover: Singer 7, 13, 31, 57, pruned 2, 3, 5, 9, 16, 21, 25 (the tiny periods make the reduced string as long as the text).</summary>
    internal static readonly int[] Periods = [7, 3, 9, 13, 2, 16, 21, 5, 25, 31, 57];

    internal static SaBuildTuning Tuning(int v) => new()
    {
        DcsBucketCapacity = 3 + (v * 7 % 20),
        DcsBaseCase = 4 + (v * 5 % 40),
        DcsChunkSlots = 16 + (v * 13 % 60),
        DcsSortCapacity = 3 + (v * 3 % 30),
        DcsRunSpillBytes = v % 4 * 40,
        DcsFanOut = 2 + (v % 4),
        DcsReservoir = 3 + (v * 11 % 40),
        DcsIoBytes = 16 + (8 * (v % 9)),
        DcsMergeBudget = 8 + (16 * (v % 7)),
        BucketCapacity = 3 + (v * 17 % 48),
        FanOut = 2 + (v % 3),
        Oversample = 1 + (v % 5),
        SpillBytes = 16 + (16 * (v % 5)),
        WindowBytes = 1 + (v * 37 % 200),
    };

    internal static SaBuildOptions Build(string scratch, int v, int threads = 1) => new()
    {
        ScratchDir = scratch,
        Threads = threads,
        Period = Periods[v % Periods.Length],
        ExactPeriod = true,
        Tuning = Tuning(v),
        Seed = (ulong)(v + 1),
        ArenaGuards = ArenaOptions.Guarded(0xCD), // every buffer that is read before it is written shows as 0xCDCDCDCD
    };

    /// <summary>The options of PH1 alone for the variant (the same sizes as <see cref="Build"/>).</summary>
    internal static ExtDcsOptions Dcs(string scratch, int v, int threads = 1)
    {
        SaBuildTuning t = Tuning(v);
        return new ExtDcsOptions
        {
            ScratchDir = scratch,
            Threads = threads,
            Seed = (ulong)(v + 1),
            BucketCapacity = t.DcsBucketCapacity!.Value,
            BaseCase = t.DcsBaseCase!.Value,
            ChunkSlots = t.DcsChunkSlots!.Value,
            SortLeafCapacity = t.DcsSortCapacity!.Value,
            RunSpillBytes = t.DcsRunSpillBytes!.Value,
            FanOut = t.DcsFanOut!.Value,
            ReservoirSize = t.DcsReservoir!.Value,
            IoBufferBytes = t.DcsIoBytes!.Value,
            MergeBufferBudget = t.DcsMergeBudget!.Value,
            WindowBytes = t.WindowBytes!.Value,
            ArenaGuards = ArenaOptions.Guarded(0xCD),
        };
    }

    /// <summary>The options of PH2 / PH3 alone for the variant.</summary>
    internal static SortOptions Sort(string scratch, int v, int threads = 1)
    {
        SaBuildTuning t = Tuning(v);
        return new SortOptions
        {
            ScratchDir = scratch,
            Threads = threads,
            Seed = (ulong)(v + 1),
            LeafCapacity = t.BucketCapacity!.Value,
            FanOut = t.FanOut!.Value,
            Oversample = t.Oversample!.Value,
            SpillBytes = t.SpillBytes!.Value,
            ArenaGuards = ArenaOptions.Guarded(0xCD),
        };
    }
}

/// <summary>One run of the stage on a case, and the comparison of everything it wrote with the oracles.</summary>
internal static class SaRun
{
    /// <summary>The result and the places of one run.</summary>
    internal sealed record Outcome(SaBuildResult Result, string SaPath, string OutDir, string ScratchDir, string TextPath, string RankPath);

    internal static string WriteText(TempDir dir, SaCase c, string tag = "") => TextOracle.WriteFile(dir, $"text{tag}.pk", c.Symbols, c.Bits);

    /// <summary>Runs the stage with the output in <c>dir/out</c> (the scratch directory is that of the options, it must not exist before and must not exist afterwards).</summary>
    internal static Outcome Run(TempDir dir, SaCase c, SaBuildOptions options, string tag = "", Action<MappedTextReader>? configure = null)
    {
        string textPath = WriteText(dir, c, tag);
        string outDir = dir.Sub("out" + tag);
        Directory.CreateDirectory(outDir);
        string saPath = Path.Combine(outDir, "o.sa");
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        configure?.Invoke(reader);
        SaBuildResult result = SaBuilder.Build(reader, saPath, options);
        return new Outcome(result, saPath, outDir, options.ScratchDir, textPath, options.RankPath ?? Path.Combine(options.ScratchDir, "ranks.bin"));
    }

    /// <summary>Fails with a description of the first difference between the suffix array file and the oracle; also checks that nothing is left behind.</summary>
    internal static void Compare(SaCase c, Outcome o, string label)
    {
        Assert.That(new FileInfo(o.SaPath).Length, Is.EqualTo(4L * (c.N + 1)), $"{label}: length of the suffix array file");
        uint[] sa = FileBits.ReadUInts(o.SaPath);
        uint[] expected = c.Sa;
        for (int k = 0; k < sa.Length; k++)
        {
            if (sa[k] != expected[k]) Assert.Fail($"{label}: SA[{k}] = {sa[k]}, expected {expected[k]} (text '{(c.Text.Length > 60 ? c.Text[..60] + "..." : c.Text)}')");
        }
        AssertNothingLeft(o, label);
    }

    internal static void AssertNothingLeft(Outcome o, string label)
    {
        Assert.That(Directory.Exists(o.ScratchDir), Is.False, $"{label}: the scratch directory is left behind");
        Assert.That(File.Exists(o.RankPath), Is.False, $"{label}: the rank file is left behind");
        Assert.That(Directory.GetFiles(o.OutDir).Select(Path.GetFileName), Is.EqualTo(new[] { Path.GetFileName(o.SaPath) }), $"{label}: the output directory holds exactly the suffix array");
        Assert.That(ProcMaps.Contains(o.SaPath) ?? false, Is.False, $"{label}: the suffix array is still mapped");
        ProcFds.AssertClosed(label, [o.SaPath, o.RankPath], [o.ScratchDir]);
    }

    /// <summary>A run that must end with an exception: nothing is left behind (no output, no rank file, no scratch).</summary>
    internal static void AssertFailsCleanly(TempDir dir, SaCase c, SaBuildOptions options, Type exception, string label, Action<MappedTextReader>? configure = null, string? textPath = null, long? n = null, int? bits = null)
    {
        textPath ??= WriteText(dir, c, "-fail");
        string outDir = dir.Sub("out-fail");
        Directory.CreateDirectory(outDir);
        string saPath = Path.Combine(outDir, "o.sa");
        using MappedTextReader reader = MappedTextReader.Open(textPath, n ?? c.N, bits ?? c.Bits);
        configure?.Invoke(reader);
        Exception? thrown = null;
        try
        {
            SaBuilder.Build(reader, saPath, options);
        }
        catch (Exception e)
        {
            thrown = e;
        }
        Assert.That(thrown, Is.Not.Null, $"{label}: the build must fail");
        Assert.That(thrown, Is.InstanceOf(exception), $"{label}: {thrown}");
        Assert.That(Directory.GetFileSystemEntries(outDir), Is.Empty, $"{label}: partial outputs are left behind");
        Assert.That(Directory.Exists(options.ScratchDir), Is.False, $"{label}: scratch is left behind");
        Assert.That(File.Exists(options.RankPath ?? Path.Combine(options.ScratchDir, "ranks.bin")), Is.False, $"{label}: the rank file is left behind");
        ProcFds.AssertClosed(label, [saPath, options.RankPath], [options.ScratchDir, outDir]);
        Assert.That(ProcMaps.Contains(saPath) ?? false, Is.False, $"{label}: the suffix array is still mapped");
        Directory.Delete(outDir);
    }

    internal static string Sha(string path) => Convert.ToHexString(SHA256.HashData(File.ReadAllBytes(path)));

    /// <summary>The families of the contract tests, each as the text of n characters.</summary>
    internal static IEnumerable<(string Name, string Text)> Families(int n)
    {
        yield return ("a^n", StringFamilies.AllSame(n));
        yield return ("fibonacci", StringFamilies.Fibonacci(n));
        yield return ("thue-morse", StringFamilies.ThueMorse(n));
        yield return ("tandem", StringFamilies.Periodic(n, 37, 5));
        yield return ("random-dna", StringFamilies.RandomDna(n, 17));
        yield return ("random-protein", StringFamilies.RandomProtein(n, 18));
        yield return ("repeat-rich", LcpChildRun.RepeatRich(n, 19));
        yield return ("n-runs", LcpChildRun.NRuns(n, 20));
        yield return ("genome-like", GenomeLike(n, 21));
    }

    /// <summary>DNA with satellites (a short unit repeated), a few long duplications and runs of N: tie-rich in the way of a real genome.</summary>
    internal static string GenomeLike(int n, ulong seed)
    {
        char[] t = StringFamilies.RandomDna(n, seed).ToCharArray();
        var rng = new SplitMix64(seed + 3000);
        for (int sat = 0; sat < 5 && n > 200; sat++)
        {
            int unit = 3 + rng.NextInt(12), len = Math.Max(30, n / 20);
            int at = rng.NextInt(n - len);
            for (int k = 0; k < len; k++) t[at + k] = t[at + (k % unit)];
        }
        for (int copy = 0; copy < 4 && n > 400; copy++)
        {
            int len = Math.Max(20, n / 10);
            int from = rng.NextInt(n - len), to = rng.NextInt(n - len);
            for (int k = 0; k < len; k++) t[to + k] = t[from + k];
        }
        for (int run = 0; run < 2 && n > 100; run++)
        {
            int len = Math.Max(5, n / 30);
            int at = rng.NextInt(n - len);
            for (int k = 0; k < len; k++) t[at + k] = 'N';
        }
        return new string(t);
    }
}

/// <summary>A rank file written from an oracle, and the packed text, for the tests of the parts that take them as input.</summary>
internal static class SaInputs
{
    /// <summary>Writes the oracle rank file of <paramref name="c"/> for the cover and returns its path.</summary>
    internal static string WriteRankFile(TempDir dir, SaCase c, DifferenceCover cover, string name = "ranks.bin")
    {
        uint[] ranks = DcsOracle.ExpectedRanks(c.Symbols, cover, out _);
        string path = dir.File(name);
        File.WriteAllBytes(path, MemoryMarshal.AsBytes(ranks.AsSpan()).ToArray());
        return path;
    }
}
