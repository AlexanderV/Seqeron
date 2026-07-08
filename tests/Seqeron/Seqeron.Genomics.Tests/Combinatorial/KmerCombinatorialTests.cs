namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the K-mer area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Kmer")]
public class KmerCombinatorialTests
{
    /// <summary>Deterministic well-mixed ACGT sequence (LCG).</summary>
    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    /// <summary>Independent brute-force k-mer counter (ground truth, no production code).</summary>
    private static Dictionary<string, int> BruteCount(string s, int k)
    {
        var d = new Dictionary<string, int>();
        for (int i = 0; i + k <= s.Length; i++)
        {
            string km = s.Substring(i, k);
            d[km] = d.GetValueOrDefault(km) + 1;
        }
        return d;
    }

    // A sequence with a 16-nt block repeated 4×, so k-mers (k ≤ 11) inside the block recur
    // ≥ 4 times — guaranteeing non-vacuous results up to minFreq 3 across all four k values.
    private static readonly string KmerSeq =
        DiverseDna(120, 0xBEEFu) + string.Join("A", Enumerable.Repeat("ACGTTGCAACGATCGT", 4)) + DiverseDna(120, 0xFADEu);

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: KMER-FIND-001 — K-mer search: most-frequent / min-count (K-mer)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 34.
    // Spec: tests/TestSpecs/KMER-FIND-001.md (canonical FindMostFrequentKmers,
    //       FindKmersWithMinCount, FindUniqueKmers).
    // Dimensions: k(4) × topN(3) × minFreq(3). Grid 4×3×3 = 36.
    //
    // Model (Rosalind BA1B; Wikipedia K-mer): the recurrent-k-mer query returns the k-mers
    // whose occurrence count is at least minFreq, ordered by count descending; taking the
    // first topN yields the "top-N most frequent k-mers occurring ≥ minFreq times".
    //
    // The combinatorial point: k, the count floor (minFreq) and the result cap (topN)
    // interact. Every returned k-mer must (a) actually occur ≥ minFreq times (verified by
    // an independent counter), (b) appear in non-increasing count order, (c) number at most
    // topN, and (d) genuinely be the highest-count eligible k-mers (top-N correctness).
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void KmerFind_TopNRecurrentKmers_AreCorrectlySelectedAndCounted(
        [Values(3, 5, 7, 11)] int k,
        [Values(1, 3, 5)] int topN,
        [Values(1, 2, 3)] int minFreq)
    {
        var results = KmerAnalyzer.FindKmersWithMinCount(KmerSeq, k, minFreq).Take(topN).ToList();
        var brute = BruteCount(KmerSeq, k);
        var eligibleCountsDesc = brute.Values.Where(v => v >= minFreq).OrderByDescending(v => v).ToList();

        results.Should().HaveCount(Math.Min(topN, eligibleCountsDesc.Count), "topN caps the eligible set");

        foreach (var (kmer, count) in results)
        {
            kmer.Length.Should().Be(k);
            brute[kmer].Should().Be(count, "the reported count matches an independent recount");
            count.Should().BeGreaterThanOrEqualTo(minFreq, "every result clears the minFreq floor");
        }

        results.Select(r => r.Count).Should().BeInDescendingOrder("results are ordered by count");
        results.Select(r => r.Count).Should().Equal(eligibleCountsDesc.Take(topN),
            "the returned counts are exactly the top-N highest eligible counts");
    }

    /// <summary>
    /// Interaction witness: minFreq is a monotone floor — raising it can only drop k-mers,
    /// so the recurrent set at a higher minFreq is a subset of the set at a lower one.
    /// </summary>
    [Test]
    public void KmerFind_MinFreq_IsMonotoneFilter()
    {
        const int k = 5;
        var at1 = KmerAnalyzer.FindKmersWithMinCount(KmerSeq, k, 1).Select(r => r.Kmer).ToHashSet();
        var at2 = KmerAnalyzer.FindKmersWithMinCount(KmerSeq, k, 2).Select(r => r.Kmer).ToHashSet();
        var at3 = KmerAnalyzer.FindKmersWithMinCount(KmerSeq, k, 3).Select(r => r.Kmer).ToHashSet();

        at3.Should().BeSubsetOf(at2);
        at2.Should().BeSubsetOf(at1);
        at1.Should().NotBeEmpty();
    }

    /// <summary>
    /// Interaction witness: FindMostFrequentKmers returns exactly the maximum-count group,
    /// which equals FindKmersWithMinCount at minFreq = the maximum count.
    /// </summary>
    [Test]
    public void KmerFind_MostFrequent_EqualsMaxCountGroup()
    {
        const int k = 5;
        var counts = BruteCount(KmerSeq, k);
        int maxCount = counts.Values.Max();

        var mostFrequent = KmerAnalyzer.FindMostFrequentKmers(KmerSeq, k).ToHashSet();
        var expectedMax = counts.Where(kvp => kvp.Value == maxCount).Select(kvp => kvp.Key).ToHashSet();

        mostFrequent.Should().BeEquivalentTo(expectedMax);
        KmerAnalyzer.FindKmersWithMinCount(KmerSeq, k, maxCount).Select(r => r.Kmer)
            .Should().BeEquivalentTo(expectedMax, "minFreq = maxCount isolates the most-frequent group");
    }

    /// <summary>
    /// Worked example (Rosalind BA1B): the most frequent 4-mers of the sample text are CATG
    /// and GCAT, each occurring three times.
    /// </summary>
    [Test]
    public void KmerFind_RosalindBa1b_WorkedExample()
    {
        const string text = "ACGTTGCATGTCGCATGATGCATGAGAGCT";

        KmerAnalyzer.FindMostFrequentKmers(text, 4).Should().BeEquivalentTo("CATG", "GCAT");
        var counts = KmerAnalyzer.CountKmers(text, 4);
        counts["CATG"].Should().Be(3);
        counts["GCAT"].Should().Be(3);
    }

    /// <summary>A sequence carrying both recurrent (4× repeated block) and unique (diverse pad) k-mers.</summary>
    private static string RecurrentAndUniqueSeq(int padLen) =>
        DiverseDna(padLen, 0x1234u) + string.Join("A", Enumerable.Repeat("ACGTTGCAACGATCGT", 4)) + DiverseDna(padLen, 0x5678u);

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: KMER-ASYNC-001 — Asynchronous / cancellable k-mer counting (K-mer)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 156.
    // Spec: canonical CountKmersAsync / CountKmers(string,k,CancellationToken). ADVANCED §10.
    // Dimensions: k(4) × seqLen(3) × parallelism(2). Grid 4×3×2 = 24 (full, exhaustive ⊇ pairwise).
    //
    // Model (Wikipedia K-mer; Compeau & Pevzner Count(Text,Pattern)): k-mer counting is the
    // overlapping-window occurrence histogram — a pure, deterministic function of (sequence, k).
    // The async overload merely runs the same computation on a Task; it MUST yield bit-identical
    // counts (parallelism is an execution detail, not a semantic one).
    //
    // Axis mapping (documented): parallelism → {Sync = CountKmers, Async = CountKmersAsync}. The
    // combinatorial point: across every k and length, both execution paths reproduce the independent
    // brute-force histogram exactly — the async path adds no nondeterminism.
    // ═══════════════════════════════════════════════════════════════════════

    public enum CountMode { Sync, Async }

    private static Dictionary<string, int> CountVia(CountMode mode, string seq, int k) => mode == CountMode.Sync
        ? KmerAnalyzer.CountKmers(seq, k)
        : KmerAnalyzer.CountKmersAsync(seq, k).GetAwaiter().GetResult();

    [Test, Combinatorial]
    public void KmerAsync_CountsMatchGroundTruth_AcrossKLengthParallelism(
        [Values(3, 5, 7, 11)] int k,
        [Values(50, 120, 240)] int seqLen,
        [Values(CountMode.Sync, CountMode.Async)] CountMode parallelism)
    {
        string seq = DiverseDna(seqLen, 0xC0FFEEu);

        var counts = CountVia(parallelism, seq, k);
        var brute = BruteCount(seq, k);

        counts.Should().BeEquivalentTo(brute, "k-mer counting is a deterministic histogram, independent of execution mode");
        counts.Values.Sum().Should().Be(seqLen - k + 1, "every overlapping window is counted exactly once");
    }

    /// <summary>
    /// Interaction witness — the async path is semantically identical to the sync path, and an
    /// already-cancelled token aborts the count (cooperative cancellation).
    /// </summary>
    [Test]
    public void KmerAsync_MatchesSync_AndHonoursCancellation()
    {
        string seq = RecurrentAndUniqueSeq(80);
        KmerAnalyzer.CountKmersAsync(seq, 6).GetAwaiter().GetResult()
            .Should().BeEquivalentTo(KmerAnalyzer.CountKmers(seq, 6), "async == sync");

        using var cts = new System.Threading.CancellationTokenSource();
        cts.Cancel();
        Action cancelled = () => KmerAnalyzer.CountKmersAsync(seq, 6, cts.Token).GetAwaiter().GetResult();
        cancelled.Should().Throw<OperationCanceledException>("a cancelled token aborts the count");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: KMER-DIST-001 — Alignment-free k-mer distance (K-mer)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 158.
    // Spec: canonical KmerDistance (Euclidean over normalized k-mer frequency vectors). ADVANCED §10.
    // Dimensions: k(4) × metric(2) × seqLen(3). Grid 4×2×3 = 24 (full, exhaustive ⊇ pairwise).
    //
    // Model (Vinga & Almeida 2003; Zielezinski 2017): an alignment-free distance maps each sequence
    // to its k-mer frequency vector f (count / (L−k+1)) and measures d = √(Σ_w (f1[w]−f2[w])²) over
    // the union of k-mers. It is a metric: d(x,x)=0 (identity of indiscernibles), d(x,y)=d(y,x)
    // (symmetry), d ≥ 0 (non-negativity), and d(x,y) > 0 for compositionally different x ≠ y.
    //
    // Axis mapping (documented — one Euclidean metric is implemented): metric → {Identity regime
    // d(x,x), Discrimination regime d(x,y)}. The combinatorial point: across every k and length the
    // identity regime gives exactly 0, the discrimination regime gives a strictly positive value that
    // equals an INDEPENDENT Euclidean recomputation and is symmetric.
    // ═══════════════════════════════════════════════════════════════════════

    public enum DistanceRegime { Identity, Discrimination }

    /// <summary>Independent ground truth: Euclidean distance between normalized k-mer frequency vectors.</summary>
    private static double EuclideanFreqDistance(string a, string b, int k)
    {
        Dictionary<string, double> Freq(string s)
        {
            var c = BruteCount(s, k);
            int total = c.Values.Sum();
            return total == 0 ? new() : c.ToDictionary(kv => kv.Key, kv => (double)kv.Value / total);
        }
        var f1 = Freq(a);
        var f2 = Freq(b);
        var keys = new HashSet<string>(f1.Keys);
        keys.UnionWith(f2.Keys);
        double sum = 0;
        foreach (var key in keys)
        {
            double d = f1.GetValueOrDefault(key) - f2.GetValueOrDefault(key);
            sum += d * d;
        }
        return Math.Sqrt(sum);
    }

    [Test, Combinatorial]
    public void KmerDistance_MetricAxioms_AcrossKMetricLength(
        [Values(3, 5, 7, 11)] int k,
        [Values(DistanceRegime.Identity, DistanceRegime.Discrimination)] DistanceRegime metric,
        [Values(50, 120, 240)] int seqLen)
    {
        string a = DiverseDna(seqLen, 0xA11CEu);

        if (metric == DistanceRegime.Identity)
        {
            KmerDistanceMetricFor(a, a, k).Should().Be(0.0, "a sequence is at zero distance from itself");
        }
        else
        {
            string b = DiverseDna(seqLen, 0xB0B0u);
            double d = KmerAnalyzer.KmerDistance(a, b, k);

            d.Should().BeGreaterThan(0, "compositionally different sequences are at positive distance");
            d.Should().BeApproximately(EuclideanFreqDistance(a, b, k), 1e-12, "matches the independent Euclidean recomputation");
            KmerAnalyzer.KmerDistance(b, a, k).Should().BeApproximately(d, 1e-12, "the distance is symmetric");
        }

        // Non-negativity holds in both regimes.
        KmerAnalyzer.KmerDistance(a, DiverseDna(seqLen, 0xCAFEu), k).Should().BeGreaterThanOrEqualTo(0);
    }

    private static double KmerDistanceMetricFor(string a, string b, int k) => KmerAnalyzer.KmerDistance(a, b, k);

    /// <summary>
    /// Interaction witness — the distance grows with compositional divergence: a sequence is closer
    /// to a lightly-edited copy of itself than to an unrelated sequence.
    /// </summary>
    [Test]
    public void KmerDistance_IncreasesWithDivergence()
    {
        const int k = 4;
        string a = DiverseDna(200, 0xD1u);
        string near = string.Concat(a.AsSpan(0, 180), DiverseDna(20, 0xD2u)); // 90% shared prefix
        string far = DiverseDna(200, 0xFFFFu);                      // unrelated

        double dNear = KmerAnalyzer.KmerDistance(a, near, k);
        double dFar = KmerAnalyzer.KmerDistance(a, far, k);

        dNear.Should().BeGreaterThan(0);
        dFar.Should().BeGreaterThan(dNear, "a more divergent sequence is farther away");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: KMER-UNIQUE-001 — Unique vs recurrent k-mers (K-mer)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 162.
    // Spec: canonical FindUniqueKmers (count == 1) / FindKmersWithMinCount (count ≥ t). ADVANCED §10.
    // Dimensions: k(4) × minCount(3) × seqLen(3). Grid 4×3×3 = 36 (full, exhaustive ⊇ pairwise).
    //
    // Model (Clavijo 2018 BioInfoLogics; Compeau & Pevzner): a "unique" k-mer occurs EXACTLY once;
    // a recurrent k-mer occurs ≥ t times. The two partition the distinct k-mers at t = 2:
    // distinct = unique ⊎ recurrent(≥2), and the recurrent sets are nested in t.
    //
    // The combinatorial point: across k, the count floor (minCount) and length, the unique set is
    // exactly the count-1 k-mers (independent recount), the min-count set is exactly the count-≥t
    // k-mers, and the two are complementary at t=2 — verified cell-by-cell against the brute histogram.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void KmerUnique_UniqueAndRecurrentSets_AcrossKMinCountLength(
        [Values(3, 5, 7, 11)] int k,
        [Values(1, 2, 3)] int minCount,
        [Values(60, 140, 220)] int seqLen)
    {
        string seq = RecurrentAndUniqueSeq((seqLen - 64) / 2 + 32); // length grows with seqLen
        var brute = BruteCount(seq, k);

        var unique = KmerAnalyzer.FindUniqueKmers(seq, k).ToHashSet();
        unique.Should().BeEquivalentTo(brute.Where(kv => kv.Value == 1).Select(kv => kv.Key),
            "unique k-mers are exactly those occurring once");

        var minCountSet = KmerAnalyzer.FindKmersWithMinCount(seq, k, minCount).ToList();
        minCountSet.Should().OnlyContain(r => brute[r.Kmer] >= minCount && r.Kmer.Length == k,
            "every reported k-mer clears the count floor");
        minCountSet.Select(r => r.Kmer).Should().BeEquivalentTo(brute.Where(kv => kv.Value >= minCount).Select(kv => kv.Key),
            "the min-count set is exactly the count-≥t k-mers");
        minCountSet.Select(r => r.Count).Should().BeInDescendingOrder("results are ordered by count");

        // Complementarity at t=2: distinct = unique ⊎ recurrent(≥2), disjointly.
        var recurrent = KmerAnalyzer.FindKmersWithMinCount(seq, k, 2).Select(r => r.Kmer).ToHashSet();
        unique.Should().NotIntersectWith(recurrent, "a k-mer cannot be both unique and recurrent");
        unique.Union(recurrent).Should().BeEquivalentTo(brute.Keys, "unique ⊎ recurrent(≥2) = all distinct k-mers");
    }

    /// <summary>
    /// Interaction witness — minCount is a monotone floor over the recurrent sets, and FindUniqueKmers
    /// is precisely the distinct k-mers that do NOT survive the ≥2 floor.
    /// </summary>
    [Test]
    public void KmerUnique_MinCountNesting_AndUniqueComplement()
    {
        const int k = 6;
        string seq = RecurrentAndUniqueSeq(80);

        var ge1 = KmerAnalyzer.FindKmersWithMinCount(seq, k, 1).Select(r => r.Kmer).ToHashSet();
        var ge2 = KmerAnalyzer.FindKmersWithMinCount(seq, k, 2).Select(r => r.Kmer).ToHashSet();
        var ge3 = KmerAnalyzer.FindKmersWithMinCount(seq, k, 3).Select(r => r.Kmer).ToHashSet();

        ge3.Should().BeSubsetOf(ge2);
        ge2.Should().BeSubsetOf(ge1);
        ge2.Should().NotBeEmpty("the repeated block guarantees recurrent k-mers");

        var unique = KmerAnalyzer.FindUniqueKmers(seq, k).ToHashSet();
        unique.Should().BeEquivalentTo(ge1.Except(ge2), "unique = distinct minus recurrent(≥2)");
        unique.Should().NotBeEmpty("the diverse padding guarantees unique k-mers");
    }
}
