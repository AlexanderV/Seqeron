// 08_DIFFERENTIAL_TESTING rows 32-34 (K-mer). The hash-map k-mer routines are checked against an
// INDEPENDENT sorted-array counting oracle (sort all windows, count consecutive runs) and a manual
// count/total frequency computation — a different data structure than the production Dictionary.

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class KmerDifferentialTests
{
    private const double Tol = 1e-12;

    // Sorted-array counting: enumerate windows, sort, count runs.
    private static Dictionary<string, int> SortedCountOracle(string seq, int k)
    {
        var result = new Dictionary<string, int>();
        if (string.IsNullOrEmpty(seq) || k <= 0 || k > seq.Length) return result;
        var s = seq.ToUpperInvariant();
        var windows = new List<string>();
        for (int i = 0; i + k <= s.Length; i++) windows.Add(s.Substring(i, k));
        windows.Sort(StringComparer.Ordinal);

        int run = 0;
        for (int i = 0; i < windows.Count; i++)
        {
            run++;
            if (i == windows.Count - 1 || windows[i] != windows[i + 1])
            {
                result[windows[i]] = run;
                run = 0;
            }
        }
        return result;
    }

    // ---- Row 32: KMER-COUNT-001 — CountKmers (hash map) vs sorted-array counting ----

    [Test]
    [Category("KMER-COUNT-001")]
    [TestCase("AAAA", 1)]
    [TestCase("AAAA", 2)]
    [TestCase("ACGTACGT", 3)]
    [TestCase("ACGTACGTACGT", 4)]
    [TestCase("AC", 5)]          // k > length -> empty
    [TestCase("GATTACA", 1)]
    public void CountKmers_MatchesSortedArrayOracle(string seq, int k)
    {
        var actual = KmerAnalyzer.CountKmers(seq, k);
        var expected = SortedCountOracle(seq, k);
        Assert.That(actual, Is.EquivalentTo(expected));
    }

    // ---- Row 33: KMER-FREQ-001 — GetKmerFrequencies vs manual count/total ----

    [Test]
    [Category("KMER-FREQ-001")]
    [TestCase("ACGTACGT", 2)]
    [TestCase("AAAA", 1)]
    [TestCase("ACGTACGTACGT", 3)]
    public void GetKmerFrequencies_MatchesCountOverTotalOracle(string seq, int k)
    {
        var actual = KmerAnalyzer.GetKmerFrequencies(seq, k);
        var counts = SortedCountOracle(seq, k);
        int total = counts.Values.Sum(); // = number of windows = n-k+1

        Assert.That(actual.Count, Is.EqualTo(counts.Count));
        foreach (var kvp in counts)
            Assert.That(actual[kvp.Key], Is.EqualTo((double)kvp.Value / total).Within(Tol), kvp.Key);

        // Frequencies form a probability distribution.
        Assert.That(actual.Values.Sum(), Is.EqualTo(1.0).Within(1e-9));
    }

    // ---- Row 34: KMER-FIND-001 — FindMostFrequentKmers vs sorted-array arg-max ----

    [Test]
    [Category("KMER-FIND-001")]
    [TestCase("ACGTACGTAC", 2)]   // "AC" appears most
    [TestCase("AAAATTTT", 1)]
    [TestCase("ACGTACGT", 4)]     // all unique -> tie, every k-mer returned
    public void FindMostFrequentKmers_MatchesSortedArgMax(string seq, int k)
    {
        var actual = KmerAnalyzer.FindMostFrequentKmers(seq, k).ToHashSet();
        var counts = SortedCountOracle(seq, k);
        int max = counts.Values.Max();
        var expected = counts.Where(c => c.Value == max).Select(c => c.Key).ToHashSet();
        Assert.That(actual, Is.EquivalentTo(expected));
    }
}
