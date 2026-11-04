using System.Diagnostics;
using System.Globalization;
using System.Text;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Review;

/// <summary>
/// Independent review of the optimized kernels (docs/refactoring/suffix-tree-upgrade/adr/REVIEW-S1a-opt.md), item 3: worst-case complexity of the LCP in
/// suffix-array order (iteration 1o.3b, "bounded fallback to the three-pass kernel").
/// <para>
/// <see cref="WorstCaseTimings"/> is an explicit measurement (not a test of correctness): for every adversarial family it times the SA-IS pointer kernel and the three LCP
/// implementations (safe reference, three-pass pointer kernel, direct scan with fallback) and the bare direct attempt, and writes one TSV row per family to the file named by
/// <c>ESA_REVIEW_OUT</c>. Parameters: <c>ESA_REVIEW_N</c> (default 1,000,000), <c>ESA_REVIEW_FAMILIES</c> (comma separated names, default all), <c>ESA_REVIEW_REPEATS</c> (default 3),
/// <c>ESA_REVIEW_DUMP</c> (directory: byte-alphabet families are also written there as Latin-1 files for the libsais driver and the Bench app).
/// </para>
/// <para>
/// The default tests pin the complexity contract itself: whatever the bound constants are mutated to, a text whose resolution is quadratic without the step budget must still
/// finish in linear time, and the result must stay the reference result.
/// </para>
/// </summary>
[TestFixture]
public class ReviewS1aOptWorstCaseTests
{
    internal static readonly string[] AllFamilies =
    [
        "random-dna", "random-binary", "random-protein", "all-same", "ab-repeat", "fibonacci", "thue-morse",
        "period2-noise", "period1000-noise", "rep-100x-1pct", "debruijn-2", "debruijn-4",
        "dup400-p10", "dup400-p15", "dup400-p20", "dup400-p25", "embedded-015pct", "embedded-040pct", "embedded-1pct", "embedded-5pct", "embedded-20pct", "two-copies", "three-copies", "tandem-1000",
        "u16-random-1000", "u16-bmp", "u16-all-same", "u16-period1000-noise", "akb-k100", "runs-L256",
    ];

    internal static string Family(string name, int n)
    {
        const ulong Seed = 20261104;
        switch (name)
        {
            case "random-dna": return StringFamilies.RandomDna(n, Seed);
            case "random-binary": return StringFamilies.Random(n, 2, Seed);
            case "random-protein": return StringFamilies.RandomProtein(n, Seed);
            case "all-same": return StringFamilies.AllSame(n);
            case "ab-repeat": return StringFamilies.AbRepeat(n);
            case "fibonacci": return StringFamilies.Fibonacci(n);
            case "thue-morse": return StringFamilies.ThueMorse(n);
            case "period2-noise": return Noisy(Repeat("ab", n), 0.001, Seed);
            case "period1000-noise": return Noisy(StringFamilies.Periodic(n, 1000, Seed), 0.001, Seed + 1);
            case "rep-100x-1pct": return RepeatedBlock(n, 100, 0.01, Seed);
            case "debruijn-2": return StringFamilies.DeBruijn(2, OrderFor(2, n), n);
            case "debruijn-4": return StringFamilies.DeBruijn(4, OrderFor(4, n), n);
            case "dup400-p10": return DuplicatedBlocks(n, 400, 0.10, Seed);
            case "dup400-p15": return DuplicatedBlocks(n, 400, 0.15, Seed);
            case "dup400-p20": return DuplicatedBlocks(n, 400, 0.20, Seed);
            case "dup400-p25": return DuplicatedBlocks(n, 400, 0.255, Seed);
            case "embedded-015pct": return EmbeddedRepeat(n, 0.0015, Seed);
            case "embedded-040pct": return EmbeddedRepeat(n, 0.004, Seed);
            case "embedded-1pct": return EmbeddedRepeat(n, 0.01, Seed);
            case "embedded-5pct": return EmbeddedRepeat(n, 0.05, Seed);
            case "embedded-20pct": return EmbeddedRepeat(n, 0.20, Seed);
            case "two-copies": return Copies(n, 2, Seed);
            case "three-copies": return Copies(n, 3, Seed);
            case "tandem-1000": return StringFamilies.Periodic(n, 1000, Seed);
            case "u16-random-1000": return StringFamilies.Random(n, 1000, Seed);
            case "u16-bmp": return StringFamilies.Bmp(n, Seed);
            case "u16-all-same": return Wide(300) + new string('a', n - 300);
            case "u16-period1000-noise": return Wide(Noisy(StringFamilies.Periodic(n, 1000, Seed), 0.001, Seed + 2), 1000, Seed + 3);
            case "akb-k100": return StringFamilies.AkB(n, 100);
            case "runs-L256": return StringFamilies.Runs(n, 256, Seed);
            default: throw new ArgumentException("unknown family " + name);
        }
    }

    [Test]
    [Explicit("Measurement for the review: ESA_REVIEW_N, ESA_REVIEW_OUT, ESA_REVIEW_FAMILIES, ESA_REVIEW_REPEATS, ESA_REVIEW_DUMP; see the class remarks")]
    [Category("ReviewTiming")]
    public void WorstCaseTimings()
    {
        int n = int.Parse(Environment.GetEnvironmentVariable("ESA_REVIEW_N") ?? "1000000", CultureInfo.InvariantCulture);
        int repeats = int.Parse(Environment.GetEnvironmentVariable("ESA_REVIEW_REPEATS") ?? "3", CultureInfo.InvariantCulture);
        string? outPath = Environment.GetEnvironmentVariable("ESA_REVIEW_OUT");
        string? dump = Environment.GetEnvironmentVariable("ESA_REVIEW_DUMP");
        string[] families = (Environment.GetEnvironmentVariable("ESA_REVIEW_FAMILIES") ?? string.Join(',', AllFamilies)).Split(',', StringSplitOptions.RemoveEmptyEntries);
        var rows = new StringBuilder();
        rows.AppendLine("family\tn\tsigma\tsa_s\tlcp_ref_s\tlcp_3pass_s\tlcp_direct_s\ttrydirect_s\tdirect_ok\tkind\tescapes\tsa_plus_direct_s\tsa_plus_3pass_s");
        foreach (string family in families)
        {
            string text = Family(family, n);
            if (dump is not null && text.All(c => c <= 0xFF))
            {
                Directory.CreateDirectory(dump);
                File.WriteAllBytes(Path.Combine(dump, $"{family}_{n}.txt"), [.. text.Select(c => (byte)c)]);
            }

            PreparedText prepared = TextPreparation.Prepare(text);
            double saS = double.MaxValue;
            int[] sa = [];
            for (int r = 0; r < Math.Max(1, repeats - 1); r++)
            {
                GC.Collect();
                var sw = Stopwatch.StartNew();
                sa = SaisKernel.Build(SaisKernelKind.Pointer, prepared);
                saS = Math.Min(saS, sw.Elapsed.TotalSeconds);
            }

            double refS = Best(repeats, () => LcpKernel.Build(LcpKernelKind.Reference, prepared, sa));
            double passS = Best(repeats, () => LcpKernel.Build(LcpKernelKind.Pointer, prepared, sa));
            CompactLcp? last = null;
            double directS = Best(repeats, () => last = LcpKernel.Build(LcpKernelKind.Direct, prepared, sa));
            double tryS = Best(repeats, () => TryDirectOnce(prepared, sa, out _));
            TryDirectOnce(prepared, sa, out bool ok);
            Assert.That(last, Is.Not.Null);
            rows.AppendLine(string.Create(CultureInfo.InvariantCulture,
                $"{family}\t{prepared.Length}\t{prepared.AlphabetSize}\t{saS:F4}\t{refS:F4}\t{passS:F4}\t{directS:F4}\t{tryS:F4}\t{(ok ? 1 : 0)}\t{last!.Kind}\t{last.OverflowCount}\t{saS + directS:F4}\t{saS + passS:F4}"));
            Console.WriteLine(rows.ToString().Split('\n')[^2]);
        }

        if (outPath is not null) File.WriteAllText(outPath, rows.ToString());
    }

    private static double Best(int repeats, Action run)
    {
        double best = double.MaxValue;
        for (int r = 0; r < repeats; r++)
        {
            GC.Collect();
            var sw = Stopwatch.StartNew();
            run();
            best = Math.Min(best, sw.Elapsed.TotalSeconds);
        }

        return best;
    }

    private static void TryDirectOnce(PreparedText prepared, int[] sa, out bool ok)
    {
        using var arena = new NativeArena();
        int n = prepared.Length;
        Buf<int> saBuf = arena.Allocate<int>(n + 1, zero: false, "sa");
        sa.CopyTo(saBuf.AsSpan(n + 1));
        CompactLcp? result;
        if (prepared.Width == SymbolWidth.Byte)
        {
            Buf<byte> t = arena.Allocate<byte>(n, zero: false, "text");
            prepared.Bytes.CopyTo(t.AsSpan(n));
            result = LcpKernel.TryDirect(t, n, saBuf, arena);
        }
        else
        {
            Buf<ushort> t = arena.Allocate<ushort>(n, zero: false, "text");
            prepared.UShorts.CopyTo(t.AsSpan(n));
            result = LcpKernel.TryDirect(t, n, saBuf, arena);
        }

        ok = result is not null;
    }

    // ---- the complexity contract (default tests) ---------------------------------------------------------------------------------

    [Test]
    public void Direct_OnAnEmbeddedLongRepeat_IsLinear_AndEqualsTheReference()
    {
        // random text with one block of 20 % of the length occurring twice: 20 % of the ranks are long (below the quarter limit, so the scan completes and the early exit
        // does not fire: the long ranks are spread over the whole suffix array) and their LCP values average 400,000, i.e. 25,000 vector steps each. Without the step budget
        // (2n steps, then the three-pass fallback) the resolution would take 2e10 steps (tens of seconds); with it, milliseconds. A mutation of ResolveStepsPerChar or
        // of the budget test shows up here as a timeout, not as a wrong result (the fallback result equals the direct result).
        const int N = 4_000_000;
        string text = EmbeddedRepeat(N, 0.20, 7);
        PreparedText prepared = TextPreparation.Prepare(text);
        int[] sa = SaisKernel.Build(SaisKernelKind.Pointer, prepared);
        string warm = new(text.AsSpan(0, 2000));
        _ = LcpKernel.Build(LcpKernelKind.Direct, TextPreparation.Prepare(warm), SaisReference.BuildSuffixArray(warm)); // JIT
        var sw = Stopwatch.StartNew();
        CompactLcp direct = LcpKernel.Build(LcpKernelKind.Direct, prepared, sa);
        double seconds = sw.Elapsed.TotalSeconds;
        CompactLcp reference = LcpKernel.Build(LcpKernelKind.Reference, prepared, sa);
        Assert.That(seconds, Is.LessThan(3.0), "direct LCP on a text with an embedded long repeat must be linear (normally ~0.1 s)");
        Assert.That(direct.Kind, Is.EqualTo(reference.Kind));
        Assert.That(direct.OverflowCount, Is.EqualTo(reference.OverflowCount));
        Assert.That(direct.ToArray(), Is.EqualTo(reference.ToArray()));
    }

    [Test]
    public void Direct_OnDuplicatedBlocks_JustBelowTheQuarterLimit_IsLinear_AndEqualsTheReference()
    {
        // about 24 % of the ranks long, average LCP 200 (12 vector steps each): the scan runs to the end, the resolution exhausts its budget, the fallback runs
        const int N = 400_000;
        string text = DuplicatedBlocks(N, 400, 0.20, 11);
        PreparedText prepared = TextPreparation.Prepare(text);
        int[] sa = SaisKernel.Build(SaisKernelKind.Pointer, prepared);
        CompactLcp direct = LcpKernel.Build(LcpKernelKind.Direct, prepared, sa);
        CompactLcp reference = LcpKernel.Build(LcpKernelKind.Reference, prepared, sa);
        Assert.That(direct.Kind, Is.EqualTo(reference.Kind));
        Assert.That(direct.OverflowCount, Is.EqualTo(reference.OverflowCount));
        Assert.That(direct.ToArray(), Is.EqualTo(reference.ToArray()));
    }

    [Test]
    [Explicit("KNOWN LIMITATION, fails today (REVIEW-S1a-opt MAJ-1): one repeat of 1 % of the text length already exhausts the 2n-step resolution budget, the direct attempt gives up and the LCP phase runs at three-pass speed (plus the wasted attempt)")]
    public void Direct_Succeeds_OnRandomDnaWithOneEmbeddedRepeatOfOnePercent()
    {
        // n = 10^6 random DNA; the block of 10,000 letters occurs twice. 10,000 ranks are long (1 %, far below the quarter limit) and their LCP values average 5,000,
        // i.e. 312 vector steps each: 3.1e6 steps in total against a budget of 2e6, so ResolveLong returns -1. In general a single repeat of length L has about L^2 / 32 steps
        // and the budget is 2n, so every repeat longer than about 8 * sqrt(n) symbols (8,000 at 10^6, 25,000 at 10^7, 80,000 at 10^8) defeats the direct scan. Real genomes
        // contain such repeats (segmental duplications, satellite arrays).
        const int N = 1_000_000;
        string text = EmbeddedRepeat(N, 0.01, 3);
        PreparedText prepared = TextPreparation.Prepare(text);
        int[] sa = SaisKernel.Build(SaisKernelKind.Pointer, prepared);
        TryDirectOnce(prepared, sa, out bool ok);
        Assert.That(ok, Is.True, "the direct scan should resolve a single 1 % repeat");
    }

    [Test]
    public void InvalidSuffixArrays_InRange_AreProcessedInLinearTime_ByBothLcpKernels()
    {
        // The kernels trust a suffix array whose entries are in range. For a wrong one the result is garbage, but the time must stay linear and nothing may be read out of
        // bounds: the Kasai-style loop carries l, which never exceeds n and drops by at most 1 per step, so its total work is <= 2n whatever the entries are; the direct scan
        // gives up through its caps. A random permutation and an all-zero array over a^n and over a periodic text (the worst shapes for an unverified carried length).
        const int N = 400_000;
        foreach (string text in new[] { StringFamilies.AllSame(N), StringFamilies.Periodic(N, 3, 9) })
        {
            PreparedText prepared = TextPreparation.Prepare(text);
            var rng = new SplitMix64(5);
            int[] permutation = new int[N + 1];
            permutation[0] = N;
            for (int i = 1; i <= N; i++) permutation[i] = i - 1;
            for (int i = N; i > 1; i--)
            {
                int j = 1 + rng.NextInt(i);
                (permutation[i], permutation[j]) = (permutation[j], permutation[i]);
            }

            int[] zeros = new int[N + 1];
            zeros[0] = N;
            foreach (int[] sa in new[] { permutation, zeros })
                foreach (LcpKernelKind kind in new[] { LcpKernelKind.Direct, LcpKernelKind.Pointer })
                {
                    var sw = Stopwatch.StartNew();
                    CompactLcp lcp = LcpKernel.Build(kind, prepared, sa);
                    Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(5.0), $"{kind} on an invalid suffix array of {N} entries");
                    Assert.That(lcp.Length, Is.EqualTo(N + 1));
                }
        }
    }

    // ---- generators --------------------------------------------------------------------------------------------------------------

    private static int OrderFor(int k, int n)
    {
        int order = 1;
        long size = k;
        while (size < n)
        {
            size *= k;
            order++;
        }

        return order;
    }

    private static string Repeat(string unit, int n)
    {
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = unit[i % unit.Length];
        return new string(c);
    }

    /// <summary>Replaces a fraction <paramref name="rate"/> of the positions (chosen at random, with repetition) by a random letter of 'a'..'d'.</summary>
    internal static string Noisy(string text, double rate, ulong seed)
    {
        var rng = new SplitMix64(seed);
        char[] c = text.ToCharArray();
        long edits = (long)(c.Length * rate);
        for (long k = 0; k < edits; k++) c[rng.NextInt(c.Length)] = (char)('a' + rng.NextInt(4));
        return new string(c);
    }

    /// <summary>A random block of n / copies letters of 'a'..'d' repeated <paramref name="copies"/> times, each copy with a fraction <paramref name="rate"/> of substitutions.</summary>
    internal static string RepeatedBlock(int n, int copies, double rate, ulong seed)
    {
        var rng = new SplitMix64(seed);
        int block = n / copies;
        var b = new char[block];
        for (int i = 0; i < block; i++) b[i] = (char)('a' + rng.NextInt(4));
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = b[i % block];
        long edits = (long)(n * rate);
        for (long k = 0; k < edits; k++) c[rng.NextInt(n)] = (char)('a' + rng.NextInt(4));
        return new string(c);
    }

    /// <summary>
    /// Random DNA of length n in which one block of <paramref name="fraction"/> * n letters occurs twice (random, block, random, the same block, random): about
    /// <paramref name="fraction"/> of the ranks have an LCP between 16 and the block length, spread over the whole suffix array.
    /// </summary>
    internal static string EmbeddedRepeat(int n, double fraction, ulong seed)
    {
        var rng = new SplitMix64(seed);
        int block = (int)(n * fraction);
        int gap = (n - (2 * block)) / 3;
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = (char)('a' + rng.NextInt(4));
        int first = gap, second = gap + block + gap;
        for (int i = 0; i < block; i++) c[second + i] = c[first + i];
        return new string(c);
    }

    /// <summary>Exactly <paramref name="copies"/> copies of one random block (the last one cut to length n).</summary>
    internal static string Copies(int n, int copies, ulong seed)
    {
        var rng = new SplitMix64(seed);
        int block = (n + copies - 1) / copies;
        var b = new char[block];
        for (int i = 0; i < block; i++) b[i] = (char)('a' + rng.NextInt(4));
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = b[i % block];
        return new string(c);
    }

    /// <summary>
    /// Blocks of <paramref name="blockLength"/> random DNA letters; each block after the first is, with probability <paramref name="copyProbability"/>, a copy of a random
    /// earlier block (its end differs from the original's neighbours, so the common prefix of a copy with its original is the block length minus the offset): about
    /// <c>copyProbability</c> of the ranks then have an LCP between 16 and <paramref name="blockLength"/> and the rest stay at about log4(n).
    /// </summary>
    internal static string DuplicatedBlocks(int n, int blockLength, double copyProbability, ulong seed)
    {
        var rng = new SplitMix64(seed);
        var c = new char[n];
        int blocks = (n + blockLength - 1) / blockLength;
        for (int b = 0; b < blocks; b++)
        {
            int start = b * blockLength;
            int len = Math.Min(blockLength, n - start);
            bool copy = b > 0 && (double)(rng.Next() >> 11) / (1UL << 53) < copyProbability;
            if (copy)
            {
                int from = rng.NextInt(b) * blockLength;
                for (int i = 0; i < len; i++) c[start + i] = c[from + i];
            }
            else
            {
                for (int i = 0; i < len; i++) c[start + i] = (char)('a' + rng.NextInt(4));
            }
        }

        return new string(c);
    }

    /// <summary>The first <paramref name="count"/> symbols are distinct code units above 255, which forces the 16-bit representation.</summary>
    private static string Wide(int count)
    {
        var sb = new StringBuilder(count);
        for (int i = 0; i < count; i++) sb.Append((char)(0x100 + i));
        return sb.ToString();
    }

    /// <summary>Maps the letters of <paramref name="text"/> onto a random-looking 16-bit alphabet of <paramref name="sigma"/> symbols (letters keep their identity, noise symbols are drawn from the whole alphabet).</summary>
    private static string Wide(string text, int sigma, ulong seed)
    {
        var rng = new SplitMix64(seed);
        var map = new char[4];
        for (int i = 0; i < 4; i++) map[i] = (char)(0x100 + rng.NextInt(sigma));
        var c = new char[text.Length];
        for (int i = 0; i < c.Length; i++) c[i] = map[text[i] - 'a'];
        // widen the alphabet: every 997th symbol becomes a distinct wide symbol (keeps the repeat structure, adds sigma)
        for (int i = 0; i < c.Length; i += 997) c[i] = (char)(0x100 + rng.NextInt(sigma));
        return new string(c);
    }
}
