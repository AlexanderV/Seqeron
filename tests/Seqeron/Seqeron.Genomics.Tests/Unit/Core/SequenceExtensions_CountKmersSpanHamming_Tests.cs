using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// B01 finisher A1-7 (F20) and A1-8 (F21): span <c>CountKmersSpan</c> (one ASCII upper-casing pass +
/// alternate-lookup counting) and span <c>HammingDistance</c> (ASCII-only case folding).
/// </summary>
[TestFixture]
public class SequenceExtensions_CountKmersSpanHamming_Tests
{
    /// <summary>Straightforward reference: one substring per window, upper-cased (ASCII input).</summary>
    private static List<KeyValuePair<string, int>> ReferenceCounts(string s, int k)
    {
        var order = new List<string>();
        var counts = new Dictionary<string, int>();
        for (int i = 0; i + k <= s.Length; i++)
        {
            string key = s.Substring(i, k).ToUpperInvariant();
            if (counts.TryGetValue(key, out int c))
                counts[key] = c + 1;
            else
            {
                counts[key] = 1;
                order.Add(key);
            }
        }
        return order.Select(key => new KeyValuePair<string, int>(key, counts[key])).ToList();
    }

    [Test]
    [Description("A1-7: equals a substring-based reference count (same keys, counts and first-occurrence order) on random ASCII input incl. lower/mixed case and N; lengths cross the stackalloc/ArrayPool threshold")]
    public void CountKmersSpan_RandomInputs_EqualReferenceCountAndOrder()
    {
        var rng = new Random(20261009);
        const string alphabet = "ACGTacgtNn";
        for (int run = 0; run < 400; run++)
        {
            int length = rng.Next(0, run % 4 == 0 ? 2000 : 300);
            var chars = new char[length];
            int mode = run % 3; // 0 upper, 1 lower, 2 mixed
            for (int i = 0; i < length; i++)
            {
                char c = alphabet[rng.Next(alphabet.Length)];
                chars[i] = mode == 0 ? char.ToUpperInvariant(c) : mode == 1 ? char.ToLowerInvariant(c) : c;
            }
            string s = new(chars);
            int k = rng.Next(1, 13);

            var observed = s.AsSpan().CountKmersSpan(k).ToList();
            Assert.That(observed, Is.EqualTo(ReferenceCounts(s, k)), $"run {run}: L={length}, k={k}");
        }
    }

    [Test]
    [Description("A1-7: Python Counter(s.upper()[i:i+3]) = skbio DNA('ACGTACGTNA').kmer_frequencies(3, overlap=True) = {ACG:2, CGT:2, GTA:1, TAC:1, GTN:1, TNA:1}")]
    public void CountKmersSpan_MixedCase_MatchesScikitBio()
    {
        var counts = "acgTAcgtNa".AsSpan().CountKmersSpan(3);
        Assert.That(counts, Is.EqualTo(new Dictionary<string, int>
        {
            ["ACG"] = 2, ["CGT"] = 2, ["GTA"] = 1, ["TAC"] = 1, ["GTN"] = 1, ["TNA"] = 1,
        }));
        Assert.That(counts.Keys, Is.EqualTo(new[] { "ACG", "CGT", "GTA", "TAC", "GTN", "TNA" }), "first-occurrence order");
    }

    [Test]
    [Description("A1-7: an input without lower-case letters is counted in place; the input span is never modified")]
    public void CountKmersSpan_DoesNotModifyInput()
    {
        char[] data = "acgtACGT".ToCharArray();
        var counts = new ReadOnlySpan<char>(data).CountKmersSpan(4);
        Assert.That(new string(data), Is.EqualTo("acgtACGT"));
        Assert.That(counts, Is.EqualTo(new Dictionary<string, int>
        {
            ["ACGT"] = 2, ["CGTA"] = 1, ["GTAC"] = 1, ["TACG"] = 1,
        }));
    }

    [Test]
    [Description("A1-7 lock (F3/F18 consistency): U+017F 'ſ' is not folded to 'S' (char.ToUpperInvariant would); U+212A Kelvin is kept as is")]
    public void CountKmersSpan_NonAscii_IsNotFoldedIntoAscii()
    {
        var counts = "GſGSK".AsSpan().CountKmersSpan(1);
        Assert.That(counts, Is.EqualTo(new Dictionary<string, int>
        {
            ["G"] = 2, ["ſ"] = 1, ["S"] = 1, ["K"] = 1,
        }));
        Assert.That("gſ".AsSpan().CountKmersSpan(2).Keys, Is.EqualTo(new[] { "Gſ" }));
    }

    [Test]
    [Description("A1-7: k > L → empty; k ≤ 0 throws (unchanged contract)")]
    public void CountKmersSpan_Contract_Unchanged()
    {
        Assert.That("ACG".AsSpan().CountKmersSpan(4), Is.Empty);
        Assert.That(ReadOnlySpan<char>.Empty.CountKmersSpan(1), Is.Empty);
        Assert.Throws<ArgumentOutOfRangeException>(() => "ACG".AsSpan().CountKmersSpan(0));
    }

    [Test]
    [Description("A1-7: a repeated k-mer costs no allocation: 100,000 windows of one k-mer allocate far less than one string per window (old code: ≥ 100,000 × 28 B)")]
    public void CountKmersSpan_RepeatedKmer_AllocatesPerDistinctKmerOnly()
    {
        string s = new('A', 100_008);
        s.AsSpan().CountKmersSpan(8); // warm-up (JIT)
        long before = GC.GetAllocatedBytesForCurrentThread();
        var counts = s.AsSpan().CountKmersSpan(8);
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.That(counts["AAAAAAAA"], Is.EqualTo(100_001));
        Assert.That(allocated, Is.LessThan(100_000), $"allocated {allocated} B");
    }

    [Test]
    [Description("A1-8: scipy hamming(list('ACGTACGTNN'), list('ACGAACGTNR')) × 10 = 2 (after str.upper())")]
    public void HammingDistance_MixedCase_MatchesScipy()
    {
        Assert.That("ACGTacgtNN".AsSpan().HammingDistance("acgaACGTnR".AsSpan()), Is.EqualTo(2));
    }

    [Test]
    [Description("A1-8: equals a char-by-char reference after ASCII upper-casing on random mixed-case inputs")]
    public void HammingDistance_RandomInputs_EqualReference()
    {
        var rng = new Random(8);
        const string alphabet = "ACGTNacgtn";
        for (int run = 0; run < 300; run++)
        {
            int length = rng.Next(0, 100);
            string a = new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
            string b = new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
            int expected = Enumerable.Range(0, length).Count(i => char.ToUpperInvariant(a[i]) != char.ToUpperInvariant(b[i]));
            Assert.That(a.AsSpan().HammingDistance(b.AsSpan()), Is.EqualTo(expected), $"'{a}' vs '{b}'");
        }
    }

    [Test]
    [Description("A1-8 lock (F3/F18 consistency): 'ſ' vs 'S'/'s' and Kelvin vs 'k' are mismatches; 'ſ' vs 'ſ' matches")]
    public void HammingDistance_NonAscii_IsNotFoldedIntoAscii()
    {
        Assert.That("AſT".AsSpan().HammingDistance("AST".AsSpan()), Is.EqualTo(1));
        Assert.That("AſT".AsSpan().HammingDistance("ast".AsSpan()), Is.EqualTo(1));
        Assert.That("K".AsSpan().HammingDistance("k".AsSpan()), Is.EqualTo(1));
        Assert.That("aſt".AsSpan().HammingDistance("AſT".AsSpan()), Is.EqualTo(0));
    }

    [Test]
    public void HammingDistance_UnequalLength_Throws()
    {
        Assert.Throws<ArgumentException>(() => "ACG".AsSpan().HammingDistance("AC".AsSpan()));
    }
}
