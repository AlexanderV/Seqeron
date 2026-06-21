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

        KmerAnalyzer.FindMostFrequentKmers(text, 4).Should().BeEquivalentTo(new[] { "CATG", "GCAT" });
        var counts = KmerAnalyzer.CountKmers(text, 4);
        counts["CATG"].Should().Be(3);
        counts["GCAT"].Should().Be(3);
    }
}
