// 08_DIFFERENTIAL_TESTING rows 156-162 (K-mer extras). Independent oracles: async==sync counting,
// forward+revcomp both-strand counts, manual Euclidean k-mer distance, base-4 odometer k-mer generation,
// naive position scan, manual composition stats, and the count==1 unique set. Counts use a sorted-array
// run count (different data structure than the production hash map).

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class KmerExtraDifferentialTests
{
    private const double Tol = 1e-12;

    private static Dictionary<string, int> SortedCount(string seq, int k)
    {
        var result = new Dictionary<string, int>();
        if (string.IsNullOrEmpty(seq) || k <= 0 || k > seq.Length) return result;
        var s = seq.ToUpperInvariant();
        var w = new List<string>();
        for (int i = 0; i + k <= s.Length; i++) w.Add(s.Substring(i, k));
        w.Sort(StringComparer.Ordinal);
        int run = 0;
        for (int i = 0; i < w.Count; i++)
        {
            run++;
            if (i == w.Count - 1 || w[i] != w[i + 1]) { result[w[i]] = run; run = 0; }
        }
        return result;
    }

    private static readonly Dictionary<char, char> Comp = new() { ['A'] = 'T', ['T'] = 'A', ['G'] = 'C', ['C'] = 'G' };
    private static string RevComp(string s)
    {
        var a = s.ToUpperInvariant().Select(c => Comp[c]).ToArray();
        Array.Reverse(a);
        return new string(a);
    }

    // ---- Row 156: KMER-ASYNC-001 — async == sync ----

    [Test]
    [Category("KMER-ASYNC-001")]
    public void CountKmersAsync_EqualsSync()
    {
        const string seq = "ACGTACGTACGT";
        var sync = KmerAnalyzer.CountKmers(seq, 3);
        var asyncResult = KmerAnalyzer.CountKmersAsync(seq, 3).GetAwaiter().GetResult();
        Assert.That(asyncResult, Is.EquivalentTo(sync));
    }

    // ---- Row 157: KMER-BOTH-001 — both-strand = forward[w] + forward[RC(w)] ----

    [Test]
    [Category("KMER-BOTH-001")]
    [TestCase("ACGTACGT", 3)]
    [TestCase("AAGGCC", 2)]
    public void CountKmersBothStrands_MatchesForwardPlusRevComp(string seq, int k)
    {
        var fwd = SortedCount(seq, k);
        var rc = SortedCount(RevComp(seq), k);
        var expected = new Dictionary<string, int>(fwd);
        foreach (var kv in rc) expected[kv.Key] = expected.GetValueOrDefault(kv.Key) + kv.Value;

        Assert.That(KmerAnalyzer.CountKmersBothStrands(seq, k), Is.EquivalentTo(expected));
    }

    // ---- Row 158: KMER-DIST-001 — Euclidean distance of frequency vectors ----

    [Test]
    [Category("KMER-DIST-001")]
    [TestCase("ACGTACGT", "ACGTTTTT", 2)]
    [TestCase("AAAA", "TTTT", 1)]
    public void KmerDistance_MatchesManualEuclidean(string s1, string s2, int k)
    {
        Dictionary<string, double> Freq(string s)
        {
            var c = SortedCount(s, k);
            double total = c.Values.Sum();
            return total == 0 ? new() : c.ToDictionary(kv => kv.Key, kv => kv.Value / total);
        }
        var f1 = Freq(s1); var f2 = Freq(s2);
        var all = new HashSet<string>(f1.Keys); all.UnionWith(f2.Keys);
        double sum = all.Sum(km => Math.Pow(f1.GetValueOrDefault(km, 0) - f2.GetValueOrDefault(km, 0), 2));
        Assert.That(KmerAnalyzer.KmerDistance(s1, s2, k), Is.EqualTo(Math.Sqrt(sum)).Within(Tol));
    }

    // ---- Row 159: KMER-GENERATE-001 — Cartesian product, base-4 odometer order ----

    [Test]
    [Category("KMER-GENERATE-001")]
    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    public void GenerateAllKmers_MatchesBase4Odometer(int k)
    {
        const string alpha = "ACGT";
        int total = (int)Math.Pow(4, k);
        var expected = new List<string>();
        for (int idx = 0; idx < total; idx++)
        {
            var chars = new char[k];
            int x = idx;
            for (int p = k - 1; p >= 0; p--) { chars[p] = alpha[x % 4]; x /= 4; }
            expected.Add(new string(chars));
        }
        Assert.That(KmerAnalyzer.GenerateAllKmers(k).ToList(), Is.EqualTo(expected));
    }

    // ---- Row 160: KMER-POSITIONS-001 — naive overlapping position scan ----

    [Test]
    [Category("KMER-POSITIONS-001")]
    [TestCase("AAAA", "AA")]
    [TestCase("ACGTACGT", "CGT")]
    [TestCase("GATTACA", "T")]
    public void FindKmerPositions_MatchesNaiveScan(string seq, string kmer)
    {
        var s = seq.ToUpperInvariant(); var km = kmer.ToUpperInvariant();
        var expected = new List<int>();
        for (int i = 0; i + km.Length <= s.Length; i++)
            if (s.Substring(i, km.Length) == km) expected.Add(i);
        Assert.That(KmerAnalyzer.FindKmerPositions(seq, kmer).ToList(), Is.EqualTo(expected));
    }

    // ---- Row 161: KMER-STATS-001 — composition stats vs manual ----

    [Test]
    [Category("KMER-STATS-001")]
    [TestCase("ACGTACGTAC", 2)]
    [TestCase("AAAAA", 1)]
    public void AnalyzeKmers_MatchesManualStats(string seq, int k)
    {
        var counts = SortedCount(seq, k);
        var vals = counts.Values.ToList();
        var st = KmerAnalyzer.AnalyzeKmers(seq, k);
        Assert.That(st.TotalKmers, Is.EqualTo(vals.Sum()));
        Assert.That(st.UniqueKmers, Is.EqualTo(counts.Count));
        Assert.That(st.MaxCount, Is.EqualTo(vals.Max()));
        Assert.That(st.MinCount, Is.EqualTo(vals.Min()));
        Assert.That(st.AverageCount, Is.EqualTo(Math.Round(vals.Average(), 2)).Within(Tol));
    }

    // ---- Row 162: KMER-UNIQUE-001 — k-mers occurring exactly once ----

    [Test]
    [Category("KMER-UNIQUE-001")]
    [TestCase("ACGTACGTAC", 3)]
    [TestCase("AAAA", 2)]
    public void FindUniqueKmers_MatchesCountEqualsOne(string seq, int k)
    {
        var expected = SortedCount(seq, k).Where(kv => kv.Value == 1).Select(kv => kv.Key).ToHashSet();
        Assert.That(KmerAnalyzer.FindUniqueKmers(seq, k).ToHashSet(), Is.EquivalentTo(expected));
    }
}
