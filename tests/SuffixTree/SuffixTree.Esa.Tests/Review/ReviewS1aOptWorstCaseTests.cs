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
        "dup400-p10", "dup400-p15", "dup400-p20", "dup400-p25", "two-copies", "three-copies", "tandem-1000",
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
    public void Direct_OnTwoCopiesOfARandomBlock_IsLinear_AndEqualsTheReference()
    {
        // B B with |B| = n / 2: the LCP of suffix i of the second copy with suffix i of the first is n/2 - i, so the sum of the LCP values is n^2 / 8 and an unbounded
        // continuation would take n^2 / 128 = 7.8e9 vector steps at n = 10^6 (seconds to minutes); the step budget (2n steps, then the three-pass fallback) makes it ~ms.
        // A mutation of ResolveStepsPerChar, of the early exit or of the quarter limit that removes the bound shows up here as a timeout, not as a wrong result.
        const int N = 1_000_000;
        string text = Copies(N, 2, 7);
        PreparedText prepared = TextPreparation.Prepare(text);
        int[] sa = SaisKernel.Build(SaisKernelKind.Pointer, prepared);
        string warm = new(text.AsSpan(0, 2000));
        _ = LcpKernel.Build(LcpKernelKind.Direct, TextPreparation.Prepare(warm), SaisReference.BuildSuffixArray(warm)); // JIT
        var sw = Stopwatch.StartNew();
        CompactLcp direct = LcpKernel.Build(LcpKernelKind.Direct, prepared, sa);
        double seconds = sw.Elapsed.TotalSeconds;
        CompactLcp reference = LcpKernel.Build(LcpKernelKind.Reference, prepared, sa);
        Assert.That(seconds, Is.LessThan(5.0), "direct LCP on B B must be linear (normally ~20 ms)");
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
