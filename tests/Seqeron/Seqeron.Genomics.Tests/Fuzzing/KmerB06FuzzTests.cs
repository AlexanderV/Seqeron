namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// B06 heavy tier — fuzz tests for the public <see cref="KmerAnalyzer"/> APIs added or changed by batch B06
/// (docs/Validation/review-2026-09/B06.md, F1–F38): counting options (F10), parallel counting (F14), Jellyfish
/// <c>histo</c> (F12), clump windows (F12), <c>dump -L/-U</c> (F9), option-aware frequency / most-frequent /
/// both-strand / statistics overloads (F7, F25, F26), the twelve distance metrics incl. D2*/D2S (F15, F19, F31),
/// Jensen–Shannon and <c>spaced -d EV</c> (F22, F28), spaced words (F22, F27, F29), exact Jaccard / containment / Mash
/// distance (F11, F24), Mash MinHash sketches and p-value (F23, F31, F35, F38), FracMinHash sketches, comparison and
/// downsampling (F30, F32–F34), and the BIC Markov-order helpers (F20).
///
/// Contract checked for every call on random strings (upper/lower case, N, IUPAC, U, gaps, '*', digits, blanks,
/// control and non-ASCII characters, null, empty, very short):
///   • a call whose arguments satisfy the documented preconditions never throws;
///   • a call that violates a documented precondition throws an <see cref="ArgumentException"/> (or a subclass:
///     <see cref="ArgumentNullException"/>, <see cref="ArgumentOutOfRangeException"/>) — never NullReference,
///     IndexOutOfRange, InvalidOperation, DivideByZero, Overflow, KeyNotFound …;
///   • results lie in their documented ranges.
/// Inputs are generated from fixed seeds, so every failure is reproducible.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class KmerB06FuzzTests
{
    #region Helpers

    private const string Pool = "ACGTACGTACGTACGTacgtNnRYKMSWBDHVUu-*. 09\t\0é漢\u2013";

    private static readonly KmerCountingOptions[] AllOptions =
    [
        KmerCountingOptions.Default, new(AcgtOnly: true), new(Canonical: true), new(Canonical: true, AcgtOnly: true),
    ];

    private static string? RandomInput(Random rng)
    {
        int kind = rng.Next(10);
        if (kind == 0) return null;
        if (kind == 1) return string.Empty;
        int len = kind <= 3 ? rng.Next(1, 4) : rng.Next(4, 120);
        string alphabet = kind >= 7 ? "ACGT" : Pool;
        var chars = new char[len];
        for (int i = 0; i < len; i++) chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    private static int RandomK(Random rng, string? s) => rng.Next(8) switch
    {
        0 => -rng.Next(0, 3),
        1 => (s?.Length ?? 0) + rng.Next(1, 3),
        _ => rng.Next(1, 7),
    };

    /// <summary>
    /// Runs <paramref name="call"/>: valid = true → must not throw; valid = false → must throw an ArgumentException
    /// family exception; valid = null → either. Any other exception type fails the test.
    /// </summary>
    private static bool Probe<T>(Func<T> call, bool? valid, string ctx, out T result)
    {
        try
        {
            result = call();
        }
        catch (ArgumentException ex)
        {
            if (valid == true)
                Assert.Fail($"{ctx}: documented-valid input threw {ex.GetType().Name}: {ex.Message}");
            result = default!;
            return false;
        }
        catch (Exception ex) when (ex is not NUnit.Framework.AssertionException)
        {
            Assert.Fail($"{ctx}: undocumented exception {ex.GetType().Name}: {ex.Message}");
            throw;
        }
        if (valid == false)
            Assert.Fail($"{ctx}: documented-invalid input did not throw");
        return true;
    }

    private static string Show(string? s) => s is null ? "null" : $"\"{s.Replace("\0", "\\0")}\"";

    private static bool KValid(string? s, int k) => string.IsNullOrEmpty(s) || k > 0;

    private static bool IsAcgtWord(string w) => w.All(c => c is 'A' or 'C' or 'G' or 'T');

    private static int Windows(string? s, int k) => s is null || k <= 0 ? 0 : Math.Max(0, s.Length - k + 1);

    #endregion

    /// <summary>
    /// Counting family: CountKmers / CountKmersParallel / DistinctKmers / GetKmerSpectrum / GetKmerFrequencies /
    /// FindMostFrequentKmers / FindUniqueKmers / FindKmersWithMinCount / CountKmersBothStrands / AnalyzeKmers under every
    /// option: k ≤ 0 with non-empty input → ArgumentOutOfRangeException; negative maxCount → ArgumentOutOfRangeException;
    /// degree 0 or &lt; −1 → ArgumentOutOfRangeException; canonical both-strand → ArgumentException; otherwise results obey
    /// Σ counts ≤ L − k + 1, keys of length k, frequencies in [0, 1] summing to 1, ranges/orders as documented.
    /// </summary>
    [Test]
    public void CountingFamily_RandomInputs_OnlyDocumentedExceptions_AndResultsInRange()
    {
        var rng = new Random(60601);
        for (int iter = 0; iter < 1500; iter++)
        {
            string? s = RandomInput(rng);
            int k = RandomK(rng, s);
            var options = AllOptions[rng.Next(AllOptions.Length)];
            bool kOk = KValid(s, k);
            string ctx = $"s={Show(s)} k={k} {options}";

            if (Probe(() => KmerAnalyzer.CountKmers(s!, k, options), kOk, "CountKmers " + ctx, out var counts))
            {
                counts.Values.Sum().Should().BeLessThanOrEqualTo(Windows(s, k), ctx);
                counts.Values.Should().OnlyContain(c => c > 0, ctx);
                counts.Keys.Should().OnlyContain(w => w.Length == k, ctx);
                if (options.SkipsNonAcgt)
                    counts.Keys.Should().OnlyContain(w => IsAcgtWord(w), ctx);
            }

            int degree = new[] { -2, -1, 0, 1, 3 }[rng.Next(5)];
            bool degreeOk = degree == -1 || degree > 0;
            if (Probe(() => KmerAnalyzer.CountKmersParallel(s!, k, options, degree), degreeOk && kOk,
                    $"CountKmersParallel degree={degree} " + ctx, out var parallel))
                parallel.Should().BeEquivalentTo(counts, ctx);

            if (Probe(() => KmerAnalyzer.DistinctKmers(s!, k, options), kOk, "DistinctKmers " + ctx, out var distinct))
                distinct.Should().BeEquivalentTo(counts!.Keys, ctx);

            if (Probe(() => KmerAnalyzer.GetKmerSpectrum(s!, k, options), kOk, "GetKmerSpectrum " + ctx, out var spectrum))
                spectrum.Sum(kv => (long)kv.Key * kv.Value).Should().Be(counts!.Values.Sum(), ctx);

            if (Probe(() => KmerAnalyzer.GetKmerFrequencies(s!, k, options), kOk, "GetKmerFrequencies " + ctx, out var freqs))
            {
                freqs.Values.Should().OnlyContain(f => f > 0 && f <= 1, ctx);
                if (freqs.Count > 0)
                    freqs.Values.Sum().Should().BeApproximately(1.0, 1e-9, ctx);
            }

            if (Probe(() => KmerAnalyzer.FindMostFrequentKmers(s!, k, options).ToList(), kOk, "FindMostFrequentKmers " + ctx, out var top)
                && counts!.Count > 0)
                top.Should().OnlyContain(w => counts[w] == counts.Values.Max(), ctx);

            if (Probe(() => KmerAnalyzer.FindUniqueKmers(s!, k, options).ToList(), kOk, "FindUniqueKmers " + ctx, out var unique))
                unique.Should().OnlyContain(w => counts![w] == 1, ctx).And.BeInAscendingOrder(StringComparer.Ordinal);

            int min = rng.Next(-3, 6), max = rng.Next(-2, 8);
            if (Probe(() => KmerAnalyzer.FindKmersWithMinCount(s!, k, min, max, options).ToList(), kOk && max >= 0,
                    $"FindKmersWithMinCount [{min},{max}] " + ctx, out var ranged))
                ranged.Should().OnlyContain(p => p.Count >= min && p.Count <= max && counts![p.Kmer] == p.Count, ctx);

            if (Probe(() => KmerAnalyzer.CountKmersBothStrands(s!, k, options), !options.Canonical && kOk,
                    "CountKmersBothStrands " + ctx, out var both))
                both.Values.Sum().Should().Be(2 * counts!.Values.Sum(), "each counted window contributes once per strand; " + ctx);

            int lower = rng.Next(-2, 4), upper = rng.Next(-1, 6);
            if (Probe(() => KmerAnalyzer.AnalyzeKmers(s!, k, options, lower, upper), null,
                    $"AnalyzeKmers [{lower},{upper}] " + ctx, out var stats) && kOk)
            {
                stats.TotalKmers.Should().BeGreaterThanOrEqualTo(0, ctx);
                stats.SingletonKmers.Should().BeLessThanOrEqualTo(stats.DistinctKmers, ctx);
                double.IsFinite(stats.AverageCount).Should().BeTrue(ctx);
                stats.Entropy.Should().BeGreaterThanOrEqualTo(0, ctx);
                if (stats.DistinctKmers > 0)
                    stats.AverageCount.Should().BeApproximately((double)stats.TotalKmers / stats.DistinctKmers, 1e-12, "F6 exact mean; " + ctx);
            }
        }
    }

    /// <summary>
    /// Jellyfish <c>histo</c> parameters (F12): low &lt; 0, high &lt; low or increment &lt; 1 → ArgumentOutOfRangeException;
    /// otherwise bins ascending, frequencies positive (or ≥ 0 with --full) and summing to the number of distinct k-mers.
    /// FindClumpWindows never throws (documented empty results) and reports the FindClumps set.
    /// </summary>
    [Test]
    public void HistogramAndClumpWindows_RandomParameters_OnlyDocumentedExceptions()
    {
        var rng = new Random(60602);
        for (int iter = 0; iter < 1500; iter++)
        {
            string? s = RandomInput(rng);
            int k = RandomK(rng, s);
            var options = AllOptions[rng.Next(AllOptions.Length)];
            long low = rng.Next(-1, 6), high = new long[] { -1, 0, 1, 3, 100, 10000 }[rng.Next(6)], inc = rng.Next(0, 8);
            bool full = rng.Next(2) == 0;
            bool valid = KValid(s, k) && low >= 0 && high >= low && inc >= 1;
            string ctx = $"s={Show(s)} k={k} {options} l={low} h={high} i={inc} f={full}";
            if (Probe(() => KmerAnalyzer.GetKmerHistogram(s!, k, options, low, high, inc, full), valid, "GetKmerHistogram " + ctx, out var hist))
            {
                hist.Select(b => b.Bin).Should().BeInAscendingOrder(ctx);
                hist.Should().OnlyContain(b => full ? b.Frequency >= 0 : b.Frequency > 0, ctx);
                hist.Sum(b => b.Frequency).Should().Be(KmerAnalyzer.DistinctKmers(s!, k, options).Count, ctx);
            }

            int window = rng.Next(-1, 40), t = rng.Next(-1, 5);
            string cctx = $"s={Show(s)} k={k} L={window} t={t}";
            if (Probe(() => KmerAnalyzer.FindClumpWindows(s!, k, window, t), true, "FindClumpWindows " + cctx, out var clumps))
            {
                Probe(() => KmerAnalyzer.FindClumps(s!, k, window, t).ToList(), true, "FindClumps " + cctx, out var set);
                clumps.Select(c => c.Kmer).Should().BeEquivalentTo(set, cctx);
                int lastStart = (s?.Length ?? 0) - window;
                clumps.Should().OnlyContain(c => c.WindowRuns.Count > 0
                    && c.WindowRuns.All(r => r.FirstWindowStart >= 0 && r.FirstWindowStart <= r.LastWindowStart && r.LastWindowStart <= lastStart), cctx);
            }
        }
    }

    /// <summary>
    /// Distance metrics (F11, F15, F19, F22, F28, F31): every <see cref="KmerDistanceMetric"/> with random k, Markov
    /// order and strand mode. Documented invalid arguments (k ≤ 0; k &gt; 12 or an order outside [−1, k) for D2*/D2S; an
    /// order ≠ 0 or bothStrands for a metric without background; D2*/D2S on a sequence with no ACGT k-mer; EV on a
    /// sequence with fewer than k letters) throw ArgumentException-family exceptions; valid calls return values in the
    /// metric's documented range (Euclidean ≤ √2, Manhattan ≤ 2, Chebyshev/JS/d2*/d2S ∈ [0, 1], cosine ∈ [0, 2],
    /// count metrics and D2 ≥ 0; d2*/d2S may be NaN when a normaliser is 0; EV may be NaN/∞ as in <c>spaced</c>).
    /// </summary>
    [Test]
    public void DistanceMetrics_RandomInputs_OnlyDocumentedExceptions_AndValuesInRange()
    {
        var rng = new Random(60603);
        var metrics = Enum.GetValues<KmerDistanceMetric>();
        for (int iter = 0; iter < 1500; iter++)
        {
            string? a = RandomInput(rng), b = rng.Next(3) == 0 ? a : RandomInput(rng);
            int k = rng.Next(10) == 0 ? rng.Next(-1, 1) : rng.Next(1, 7);
            if (rng.Next(40) == 0) k = KmerAnalyzer.MaxBackgroundAdjustedK + 1;
            var metric = metrics[rng.Next(metrics.Length)];
            int order = rng.Next(4) == 0 ? rng.Next(-2, k + 1) : 0;
            bool both = rng.Next(3) == 0;
            bool background = metric is KmerDistanceMetric.D2Star or KmerDistanceMetric.D2Shepherd;
            bool? valid = k <= 0 ? false
                : !background && order != 0 ? false
                : !background && metric != KmerDistanceMetric.SpacedEvolutionary && both ? false
                : background && (k > KmerAnalyzer.MaxBackgroundAdjustedK || order < -1 || order >= k) ? false
                : background || metric == KmerDistanceMetric.SpacedEvolutionary ? null // data-dependent (no ACGT k-mer / too short)
                : true;
            string ctx = $"a={Show(a)} b={Show(b)} k={k} {metric} r={order} both={both}";
            if (!Probe(() => KmerAnalyzer.KmerDistance(a!, b!, k, metric, order, both), valid, "KmerDistance " + ctx, out double d))
                continue;
            switch (metric)
            {
                case KmerDistanceMetric.Euclidean: d.Should().BeInRange(0, Math.Sqrt(2) + 1e-12, ctx); break;
                case KmerDistanceMetric.Manhattan: d.Should().BeInRange(0, 2 + 1e-12, ctx); break;
                case KmerDistanceMetric.Chebyshev:
                case KmerDistanceMetric.JensenShannon: d.Should().BeInRange(0, 1 + 1e-12, ctx); break;
                case KmerDistanceMetric.Cosine: d.Should().BeInRange(0, 2, ctx); break;
                case KmerDistanceMetric.D2Star:
                case KmerDistanceMetric.D2Shepherd:
                    if (!double.IsNaN(d)) d.Should().BeInRange(0, 1, "F31 clamp; " + ctx);
                    break;
                case KmerDistanceMetric.SpacedEvolutionary: break; // spaced's own NaN/∞/1.2 conventions
                default:
                    double.IsFinite(d).Should().BeTrue(ctx);
                    d.Should().BeGreaterThanOrEqualTo(0, ctx);
                    break;
            }
        }
    }

    /// <summary>
    /// D2*/D2S and BIC helpers (F15, F20): <c>BackgroundAdjustedD2</c> with a null sequence → ArgumentNullException
    /// (documented), the metric overload treats null as empty → ArgumentException (F21); MarkovOrderBic(s) / Bics /
    /// SelectMarkovOrder: null → ArgumentNullException, negative order → ArgumentOutOfRangeException, no ACGT (r+1)-mer →
    /// ArgumentException; otherwise finite BICs and a selected order in [0, maxOrder]; ParseDistanceMetric accepts only
    /// the documented names.
    /// </summary>
    [Test]
    public void BackgroundModelHelpers_RandomInputs_OnlyDocumentedExceptions()
    {
        var rng = new Random(60604);
        Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.BackgroundAdjustedD2(null!, "ACGTACGT", 2));
        Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.BackgroundAdjustedD2("ACGTACGT", null!, 2));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.KmerDistance(null!, "ACGTACGT", 2, KmerDistanceMetric.D2Star));
        for (int iter = 0; iter < 600; iter++)
        {
            string? s = RandomInput(rng);
            int order = rng.Next(-1, 4);
            string ctx = $"s={Show(s)} r={order}";
            bool? valid = s is null || order < 0 ? false : null;
            if (Probe(() => KmerAnalyzer.MarkovOrderBic(s!, order), valid, "MarkovOrderBic " + ctx, out double bic))
                double.IsNaN(bic).Should().BeFalse(ctx);
            if (Probe(() => KmerAnalyzer.MarkovOrderBics(s!, order), valid, "MarkovOrderBics " + ctx, out var bics))
                bics.Should().HaveCount(order + 1, ctx);
            if (Probe(() => KmerAnalyzer.SelectMarkovOrder(s!, order), valid, "SelectMarkovOrder " + ctx, out int chosen))
                chosen.Should().BeInRange(0, order, ctx);

            string name = new[] { "euclidean", " D2S ", "js", "EV", "cosine", "d2star", "", "bogus", "d3", "manhattan2" }[rng.Next(10)];
            bool known = name.Trim().ToLowerInvariant() is "" or "euclidean" or "d2s" or "js" or "ev" or "cosine" or "d2star";
            Probe(() => KmerAnalyzer.ParseDistanceMetric(name), known, $"ParseDistanceMetric '{name}'", out _);
        }
    }

    /// <summary>
    /// Spaced words (F22, F27, F29): random patterns (valid, empty, not starting/ending with '1', foreign symbols,
    /// unequal weights, null) through CountSpacedWords and SpacedWordDistance in both strand modes — malformed patterns,
    /// canonical options and D2*/D2S metrics → ArgumentException family; literal counts sum to L − ℓ + 1; ACGT-only words
    /// are upper-case ACGT; distances obey the metric ranges.
    /// </summary>
    [Test]
    public void SpacedWords_RandomPatternsAndInputs_OnlyDocumentedExceptions()
    {
        var rng = new Random(60605);
        string?[] patternPool = ["1", "11", "101", "1101", "1011", "11011", "10101", "", "0", "10", "01", "1x1", "1 1", null];
        KmerDistanceMetric[] metrics =
        [
            KmerDistanceMetric.Euclidean, KmerDistanceMetric.JensenShannon, KmerDistanceMetric.EuclideanCounts,
            KmerDistanceMetric.Cosine, KmerDistanceMetric.SpacedEvolutionary, KmerDistanceMetric.D2Star, KmerDistanceMetric.Manhattan,
        ];
        for (int iter = 0; iter < 1200; iter++)
        {
            string? a = RandomInput(rng), b = RandomInput(rng);
            string? pattern = patternPool[rng.Next(patternPool.Length)];
            bool patternOk = pattern is { Length: > 0 } && pattern[0] == '1' && pattern[^1] == '1' && pattern.All(c => c is '0' or '1');
            var options = AllOptions[rng.Next(AllOptions.Length)];
            string ctx = $"a={Show(a)} p={pattern ?? "null"} {options}";

            if (Probe(() => KmerAnalyzer.CountSpacedWords(a!, pattern!, options), patternOk && !options.Canonical,
                    "CountSpacedWords " + ctx, out var words))
            {
                int windows = a is null || pattern!.Length > a.Length ? 0 : a.Length - pattern.Length + 1;
                if (options.AcgtOnly)
                {
                    words.Values.Sum().Should().BeLessThanOrEqualTo(windows, ctx);
                    words.Keys.Should().OnlyContain(w => IsAcgtWord(w), ctx);
                }
                else
                {
                    words.Values.Sum().Should().Be(windows, ctx);
                }
            }

            var patterns = Enumerable.Range(0, rng.Next(1, 4)).Select(_ => patternPool[rng.Next(patternPool.Length)]).ToList();
            var metric = metrics[rng.Next(metrics.Length)];
            bool both = rng.Next(2) == 0;
            bool allOk = patterns.All(p => p is { Length: > 0 } && p[0] == '1' && p[^1] == '1' && p.All(c => c is '0' or '1'));
            bool sameWeight = allOk && patterns.Select(p => p!.Count(c => c == '1')).Distinct().Count() == 1;
            bool? valid = !allOk || !sameWeight || options.Canonical || metric == KmerDistanceMetric.D2Star ? false
                : metric == KmerDistanceMetric.SpacedEvolutionary ? null // EV: equal lengths and long enough reads required
                : true;
            string dctx = $"a={Show(a)} b={Show(b)} P=[{string.Join(',', patterns.Select(p => p ?? "null"))}] {metric} {options} both={both}";
            if (Probe(() => KmerAnalyzer.SpacedWordDistance(a!, b!, patterns!, metric, options, both), valid,
                    "SpacedWordDistance " + dctx, out double d) && metric != KmerDistanceMetric.SpacedEvolutionary)
            {
                double.IsFinite(d).Should().BeTrue(dctx);
                d.Should().BeGreaterThanOrEqualTo(0, dctx);
                if (metric == KmerDistanceMetric.JensenShannon)
                    d.Should().BeLessThanOrEqualTo(1 + 1e-12, dctx);
            }
        }
        Assert.Throws<ArgumentNullException>(() => KmerAnalyzer.SpacedWordDistance("ACGT", "ACGT", null!));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.SpacedWordDistance("ACGT", "ACGT", []));
    }

    /// <summary>
    /// Exact Jaccard / containment / Mash distance (F11, F24): k ≤ 0 → ArgumentOutOfRangeException; otherwise values in
    /// [0, 1] for any input (null = empty); MashDistanceFromJaccard rejects J outside [0, 1] and NaN.
    /// </summary>
    [Test]
    public void ExactSetSimilarities_RandomInputs_InUnitInterval()
    {
        var rng = new Random(60606);
        for (int iter = 0; iter < 1500; iter++)
        {
            string? a = RandomInput(rng), b = RandomInput(rng);
            int k = rng.Next(8) == 0 ? rng.Next(-2, 1) : rng.Next(1, 10);
            var options = AllOptions[rng.Next(AllOptions.Length)];
            string ctx = $"a={Show(a)} b={Show(b)} k={k} {options}";
            foreach (var (name, f) in new (string, Func<double>)[]
                     {
                         ("Jaccard", () => KmerAnalyzer.JaccardSimilarity(a!, b!, k, options)),
                         ("Containment", () => KmerAnalyzer.ContainmentIndex(a!, b!, k, options)),
                         ("Mash", () => KmerAnalyzer.MashDistance(a!, b!, k, options)),
                     })
            {
                if (Probe(f, k > 0, name + " " + ctx, out double v))
                    v.Should().BeInRange(0.0, 1.0, name + " " + ctx);
            }

            double j = new[] { -0.1, 0.0, 0.3, 1.0, 1.5, double.NaN, double.Epsilon }[rng.Next(7)];
            Probe(() => KmerAnalyzer.MashDistanceFromJaccard(j, k), k > 0 && j >= 0 && j <= 1, $"MashDistanceFromJaccard J={j} k={k}", out _);
        }
    }

    /// <summary>
    /// Mash sketches (F23, F31, F35, F38): CreateMinHashSketch (k outside 1..32 or s &lt; 1 → ArgumentOutOfRangeException),
    /// CompareMinHashSketches (different k / seed / canonical, or a malformed sketch → ArgumentException), MashPValue
    /// (negative inputs or k outside 1..32 → ArgumentOutOfRangeException, x &gt; s or x ≥ 1 with an empty set →
    /// ArgumentException) and MinHashSketch.FromHashes: sketches sorted and distinct, Jaccard/distance/p-value in [0, 1].
    /// </summary>
    [Test]
    public void MinHashSketches_RandomInputs_OnlyDocumentedExceptions()
    {
        var rng = new Random(60607);
        for (int iter = 0; iter < 800; iter++)
        {
            string? a = RandomInput(rng), b = RandomInput(rng);
            int k = rng.Next(10) == 0 ? new[] { 0, -1, 33 }[rng.Next(3)] : rng.Next(1, 33);
            int s = rng.Next(10) == 0 ? 0 : new[] { 1, 3, 50, 1000 }[rng.Next(4)];
            bool canonical = rng.Next(2) == 0;
            uint seed = rng.Next(4) == 0 ? (uint)rng.Next() : KmerAnalyzer.DefaultMashSeed;
            bool valid = k >= 1 && k <= 32 && s >= 1;
            string ctx = $"a={Show(a)} k={k} s={s} canonical={canonical} seed={seed}";
            if (!Probe(() => KmerAnalyzer.CreateMinHashSketch(a!, k, s, canonical, seed), valid, "CreateMinHashSketch " + ctx, out var sa))
                continue;
            sa.Hashes.Count.Should().BeLessThanOrEqualTo(s, ctx);
            sa.Hashes.Should().BeInAscendingOrder(ctx).And.OnlyHaveUniqueItems(ctx);
            if (k <= 16)
                sa.Hashes.Should().OnlyContain(h => h <= uint.MaxValue, "32-bit Mash hashes for k ≤ 16; " + ctx);

            int kb = rng.Next(5) == 0 ? Math.Max(1, k - 1) : k;
            bool cb = rng.Next(5) == 0 ? !canonical : canonical;
            var sb = KmerAnalyzer.CreateMinHashSketch(b!, kb, s, cb, seed);
            bool same = kb == k && cb == canonical;
            if (Probe(() => KmerAnalyzer.CompareMinHashSketches(sa, sb), same, "CompareMinHashSketches " + ctx, out var cmp))
            {
                cmp.SharedHashes.Should().BeInRange(0, Math.Max(0, cmp.Denominator), ctx);
                cmp.Jaccard.Should().BeInRange(0.0, 1.0, ctx);
                cmp.Distance.Should().BeInRange(0.0, 1.0, ctx);
                cmp.PValue.Should().BeInRange(0.0, 1.0, ctx);
            }

            // Malformed sketches (F38): Use64 inconsistent with K, unsorted hashes.
            if (sa.Hashes.Count >= 2)
            {
                var unsorted = sa with { Hashes = sa.Hashes.Reverse().ToArray() };
                Probe(() => KmerAnalyzer.CompareMinHashSketches(unsorted, sa), false, "unsorted sketch " + ctx, out _);
            }
            Probe(() => KmerAnalyzer.CompareMinHashSketches(sa with { Use64 = !sa.Use64 }, sa with { Use64 = !sa.Use64 }), false,
                "Use64 ≠ (K > 16) " + ctx, out _);

            long x = rng.Next(-1, 6), l1 = rng.Next(-1, 3) * rng.Next(0, 5000), l2 = rng.Next(0, 3) * rng.Next(0, 5000), ss = rng.Next(-1, 6);
            int pk = rng.Next(8) == 0 ? new[] { 0, 33, 600 }[rng.Next(3)] : rng.Next(1, 33);
            bool pValid = x >= 0 && l1 >= 0 && l2 >= 0 && ss >= 0 && pk >= 1 && pk <= 32 && x <= ss && (x == 0 || (l1 > 0 && l2 > 0));
            if (Probe(() => KmerAnalyzer.MashPValue(x, l1, l2, pk, ss), pValid, $"MashPValue x={x} l1={l1} l2={l2} k={pk} s={ss}", out double p))
                p.Should().BeInRange(0.0, 1.0);

            var raw = Enumerable.Range(0, rng.Next(0, 8)).Select(_ => rng.Next(3) == 0 ? (ulong)rng.NextInt64() << 1 : (ulong)rng.Next()).ToList();
            int fk = rng.Next(1, 33);
            bool rawOk = MinHashSketchFitsWidth(raw, fk, s);
            Probe(() => MinHashSketch.FromHashes(fk, s, canonical, seed, 100, raw), valid && rawOk ? true : s < 1 ? false : null,
                $"FromHashes k={fk} s={s}", out _);
        }
    }

    private static bool MinHashSketchFitsWidth(List<ulong> raw, int k, int s) =>
        k > 16 || s < 1 || raw.Distinct().Order().Take(s).All(h => h <= uint.MaxValue);

    /// <summary>
    /// FracMinHash (F30, F32–F34): CreateFracMinHashSketch (k &lt; 1 or scaled outside 1..2^32 − 1 →
    /// ArgumentOutOfRangeException), DownsampleFracMinHash (smaller scaled → ArgumentException), CompareFracMinHashSketches
    /// (different scaled without downsample, different k/seed/canonical, malformed sketch → ArgumentException): every kept
    /// hash ≤ max_hash, abundances positive, comparison values in [0, 1].
    /// </summary>
    [Test]
    public void FracMinHashSketches_RandomInputs_OnlyDocumentedExceptions()
    {
        var rng = new Random(60608);
        long[] scaledPool = [0, 1, 2, 7, 1000, uint.MaxValue, (long)uint.MaxValue + 1, -5];
        for (int iter = 0; iter < 800; iter++)
        {
            string? a = RandomInput(rng), b = RandomInput(rng);
            int k = rng.Next(10) == 0 ? rng.Next(-1, 1) : rng.Next(1, 40);
            long scaled = scaledPool[rng.Next(scaledPool.Length)];
            bool track = rng.Next(2) == 0, canonical = rng.Next(3) != 0;
            bool valid = k >= 1 && scaled >= 1 && scaled <= uint.MaxValue;
            string ctx = $"a={Show(a)} k={k} S={scaled} track={track} canonical={canonical}";
            if (!Probe(() => KmerAnalyzer.CreateFracMinHashSketch(a!, k, scaled, canonical, trackAbundance: track), valid,
                    "CreateFracMinHashSketch " + ctx, out var fa))
                continue;
            fa.Hashes.Should().BeInAscendingOrder(ctx).And.OnlyHaveUniqueItems(ctx).And.OnlyContain(h => h <= fa.MaxHash, ctx);
            fa.MaxHash.Should().Be(KmerAnalyzer.FracMinHashMaxHash(scaled), ctx);
            if (track)
                fa.Abundances.Should().HaveCount(fa.Hashes.Count, ctx).And.OnlyContain(v => v >= 1, ctx);
            else
                fa.Abundances.Should().BeNull(ctx);

            long newScaled = new[] { scaled - 1, scaled, scaled * 3, (long)uint.MaxValue, (long)uint.MaxValue + 1 }[rng.Next(5)];
            bool dValid = newScaled >= scaled && newScaled >= 1 && newScaled <= uint.MaxValue;
            if (Probe(() => KmerAnalyzer.DownsampleFracMinHash(fa, newScaled), dValid, $"Downsample {scaled}->{newScaled} " + ctx, out var down))
                down.Hashes.Should().BeSubsetOf(fa.Hashes, ctx);

            long scaledB = rng.Next(3) == 0 ? new[] { 1L, 5L, 1000L }[rng.Next(3)] : scaled;
            int kb = rng.Next(6) == 0 ? k + 1 : k;
            var fb = KmerAnalyzer.CreateFracMinHashSketch(b!, kb, scaledB, canonical, trackAbundance: rng.Next(2) == 0);
            bool downsample = rng.Next(2) == 0;
            bool cValid = kb == k && (scaledB == scaled || downsample);
            if (Probe(() => KmerAnalyzer.CompareFracMinHashSketches(fa, fb, downsample), cValid,
                    $"CompareFracMinHash S={scaled}/{scaledB} k={k}/{kb} downsample={downsample} " + ctx, out var cmp))
            {
                cmp.SharedHashes.Should().BeLessThanOrEqualTo(cmp.UnionHashes, ctx);
                foreach (double v in new[] { cmp.Jaccard, cmp.ContainmentAInB, cmp.ContainmentBInA, cmp.MaxContainment })
                    v.Should().BeInRange(0.0, 1.0, ctx);
                foreach (double? v in new[] { cmp.AngularSimilarity, cmp.WeightedContainmentAInB, cmp.WeightedContainmentBInA })
                    if (v is { } x) x.Should().BeInRange(0.0, 1.0, ctx);
            }

            if (fa.Hashes.Count > 0)
            {
                var badMax = fa with { MaxHash = fa.MaxHash - 1 };
                Probe(() => KmerAnalyzer.CompareFracMinHashSketches(badMax, badMax), false, "MaxHash ≠ max_hash(S) " + ctx, out _);
            }
        }
    }
}
