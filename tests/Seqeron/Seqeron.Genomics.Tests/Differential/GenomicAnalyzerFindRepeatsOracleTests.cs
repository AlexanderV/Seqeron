// ST-S0 WP 0.4 — GenomicAnalyzer.FindRepeats: the O(n^2) GetAllSuffixes implementation was replaced by a
// suffix-tree traversal. The previous implementation is kept here, unchanged, as the REFERENCE ORACLE: the
// new method must return the same RepeatInfo sequence (same strings, same ascending positions, same order)
// on random and adversarial inputs. A separate allocation test pins the complexity (the old method
// materialised every suffix: ~n^2 bytes).

using System.Text;

namespace Seqeron.Genomics.Tests.Differential;

[TestFixture]
public class GenomicAnalyzerFindRepeatsOracleTests
{
    // ───────────────────────── reference oracle (verbatim copy of the pre-WP-0.4 implementation) ─────────────────────────

    private static IEnumerable<RepeatInfo> ReferenceFindRepeats(DnaSequence sequence, int minLength)
    {
        var tree = sequence.SuffixTree;
        var found = new HashSet<string>();

        int effectiveMinLength = Math.Max(1, minLength);

        var suffixes = tree.GetAllSuffixes().OrderBy(s => s).ToList();

        for (int i = 0; i < suffixes.Count - 1; i++)
        {
            string s1 = suffixes[i];
            string s2 = suffixes[i + 1];

            int lcpLen = 0;
            int maxLen = Math.Min(s1.Length, s2.Length);
            while (lcpLen < maxLen && s1[lcpLen] == s2[lcpLen])
            {
                lcpLen++;
            }

            for (int prefixLen = effectiveMinLength; prefixLen <= lcpLen; prefixLen++)
            {
                string repeat = s1.Substring(0, prefixLen);
                if (found.Add(repeat))
                {
                    var positions = tree.FindAllOccurrences(repeat).OrderBy(p => p).ToList();
                    if (positions.Count >= 2)
                    {
                        yield return new RepeatInfo(repeat, positions);
                    }
                }
            }
        }
    }

    // ───────────────────────── inputs ─────────────────────────

    private static string RandomDna(int length, int seed, string alphabet = "ACGT")
    {
        var rng = new Random(seed);
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(alphabet[rng.Next(alphabet.Length)]);
        return sb.ToString();
    }

    private static string Fibonacci(int steps)
    {
        string a = "A", b = "AC";
        for (int i = 0; i < steps; i++)
            (a, b) = (b, b + a);
        return b;
    }

    private static string ThueMorse(int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(int.PopCount(i) % 2 == 0 ? 'A' : 'C');
        return sb.ToString();
    }

    /// <summary>De Bruijn sequence B(k, n) over <paramref name="alphabet"/> (cyclic form, then linearised by appending n-1 symbols).</summary>
    private static string DeBruijn(string alphabet, int n)
    {
        int k = alphabet.Length;
        var a = new int[k * n];
        var seq = new List<int>();
        void Db(int t, int p)
        {
            if (t > n)
            {
                if (n % p == 0)
                    for (int j = 1; j <= p; j++) seq.Add(a[j]);
            }
            else
            {
                a[t] = a[t - p];
                Db(t + 1, p);
                for (int j = a[t - p] + 1; j < k; j++)
                {
                    a[t] = j;
                    Db(t + 1, t);
                }
            }
        }
        Db(1, 1);
        var sb = new StringBuilder();
        foreach (int i in seq) sb.Append(alphabet[i]);
        for (int i = 0; i < n - 1; i++) sb.Append(alphabet[seq[i]]);
        return sb.ToString();
    }

    private static IEnumerable<string> AdversarialTexts()
    {
        yield return "A";
        yield return "AC";
        yield return "AA";
        yield return "AAA";
        yield return new string('A', 12);
        yield return new string('A', 64);
        yield return "AAAAAAAAACAAAAAAAAA";
        yield return "CAAAAAAAAAAAAAAAAAAAC";
        yield return string.Concat(Enumerable.Repeat("AC", 20));
        yield return string.Concat(Enumerable.Repeat("ACG", 15));
        yield return string.Concat(Enumerable.Repeat("ACGT", 12));
        yield return string.Concat(Enumerable.Repeat("ACGTTGCA", 6));
        yield return Fibonacci(8);            // 89 chars over {A,C}
        yield return Fibonacci(10);           // 233 chars
        yield return ThueMorse(64);
        yield return ThueMorse(257);
        yield return DeBruijn("ACGT", 3);
        yield return DeBruijn("AC", 6);
        yield return DeBruijn("ACGT", 4);
        yield return "ACGTACGTAC" + "TTTT" + "ACGTACGTAC";
        yield return "GATTACA";
        yield return "ATCGATCGA";
        yield return "ACGTAGCTAGCTAGCTTACGT";
    }

    private static readonly int[] MinLengths = { int.MinValue, -5, 0, 1, 2, 3, 4, 7, 12, 100 };

    // ───────────────────────── differential tests ─────────────────────────

    private static void AssertSameAsOracle(string text, int minLength)
    {
        var dna = new DnaSequence(text);
        var expected = ReferenceFindRepeats(dna, minLength).ToList();
        var actual = GenomicAnalyzer.FindRepeats(dna, minLength).ToList();

        string context = $"text='{(text.Length > 60 ? text[..60] + "..." : text)}' (n={text.Length}), minLength={minLength}";
        Assert.That(actual.Select(r => r.Sequence).ToList(), Is.EqualTo(expected.Select(r => r.Sequence).ToList()),
            "sequences and order differ: " + context);
        for (int i = 0; i < expected.Count; i++)
        {
            Assert.That(actual[i].Positions.ToArray(), Is.EqualTo(expected[i].Positions.ToArray()),
                $"positions of repeat #{i} '{expected[i].Sequence}' differ: " + context);
        }
    }

    [Test]
    public void FindRepeats_EmptySequence_YieldsNothing()
    {
        var dna = new DnaSequence("");

        Assert.That(GenomicAnalyzer.FindRepeats(dna, 1), Is.Empty);
        Assert.That(ReferenceFindRepeats(dna, 1), Is.Empty);
    }

    [Test]
    public void FindRepeats_AdversarialTexts_MatchTheReferenceOracle()
    {
        int compared = 0;
        foreach (string text in AdversarialTexts())
        {
            foreach (int minLength in MinLengths)
            {
                AssertSameAsOracle(text, minLength);
                compared++;
            }
        }
        Assert.That(compared, Is.GreaterThan(200));
    }

    [Test]
    public void FindRepeats_RandomDna_MatchesTheReferenceOracle([Values(2, 3, 5, 9, 17, 40, 120, 400, 1200)] int length)
    {
        // Several seeds per length, several alphabets (binary alphabets create many nested repeats).
        for (int seed = 1; seed <= 6; seed++)
        {
            foreach (string alphabet in new[] { "ACGT", "AC", "AT", "ACG" })
            {
                string text = RandomDna(length, seed * 1000 + length, alphabet);
                foreach (int minLength in new[] { -1, 1, 2, 3, 5, 8 })
                    AssertSameAsOracle(text, minLength);
            }
        }
    }

    [Test]
    public void FindRepeats_RandomDna_WithPlantedLongRepeat_MatchesTheReferenceOracle()
    {
        for (int seed = 1; seed <= 5; seed++)
        {
            string unit = RandomDna(60, seed);
            string text = RandomDna(150, seed + 100) + unit + RandomDna(90, seed + 200) + unit + RandomDna(30, seed + 300) + unit[10..50];
            foreach (int minLength in new[] { 1, 4, 10, 30, 59, 60, 61 })
                AssertSameAsOracle(text, minLength);
        }
    }

    // ───────────────────────── complexity ─────────────────────────

    [Test]
    public void FindRepeats_RandomDna_AllocatesLinearlyNotQuadratically()
    {
        // GetAllSuffixes materialised every suffix: for n = 20 000 that is ~n^2 bytes = 400 MB of strings.
        // The tree traversal needs a few MB (visitor lists, one key list per internal node, the few results).
        const int n = 20_000;
        var dna = new DnaSequence(RandomDna(n, seed: 7));
        _ = dna.SuffixTree; // build the tree outside the measurement (construction is the same for both methods)

        long before = GC.GetAllocatedBytesForCurrentThread();
        int count = GenomicAnalyzer.FindRepeats(dna, minLength: 12).Count();
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.That(count, Is.GreaterThan(0).And.LessThan(1_000), "a few 12-mer repeats are expected in 20 kbp of random DNA");
        Assert.That(allocated, Is.LessThan(2_000L * n),
            $"FindRepeats allocated {allocated:N0} bytes for n={n}: that is not O(n); the old O(n^2) method needed ~{(long)n * n:N0}");
    }

    [Test]
    public void FindRepeats_LargeRandomDna_ProducesVerifiedRepeats()
    {
        // 100 kbp: the old method would have needed ~10 GB for the suffix list. Every reported repeat is
        // verified directly against the text; the count of repeats of a fixed length is cross-checked with
        // a k-mer dictionary built independently of the suffix tree.
        const int n = 100_000;
        const int minLength = 14;
        string text = RandomDna(n, seed: 11);
        var dna = new DnaSequence(text);

        var repeats = GenomicAnalyzer.FindRepeats(dna, minLength).ToList();

        var occurrences = new Dictionary<string, List<int>>();
        for (int i = 0; i + minLength <= n; i++)
        {
            string kmer = text.Substring(i, minLength);
            if (!occurrences.TryGetValue(kmer, out var list))
                occurrences[kmer] = list = new List<int>(1);
            list.Add(i);
        }
        int repeatedKmers = occurrences.Values.Count(l => l.Count >= 2);

        Assert.Multiple(() =>
        {
            Assert.That(repeats.Count(r => r.Length == minLength), Is.EqualTo(repeatedKmers),
                "number of repeated 14-mers must equal the independent k-mer count");
            Assert.That(repeats.All(r => r.Length >= minLength && r.Count >= 2), Is.True);
            Assert.That(repeats.Select(r => r.Sequence).Distinct().Count(), Is.EqualTo(repeats.Count), "every repeat is reported once");
            foreach (var r in repeats.Where(r => r.Length == minLength))
                Assert.That(r.Positions.ToArray(), Is.EqualTo(occurrences[r.Sequence].ToArray()), $"positions of {r.Sequence}");
            foreach (var r in repeats.Take(200))
            {
                Assert.That(r.Positions, Is.Ordered.Ascending.And.Unique);
                foreach (int p in r.Positions)
                    Assert.That(string.CompareOrdinal(text, p, r.Sequence, 0, r.Length), Is.EqualTo(0));
            }
        });
    }
}
