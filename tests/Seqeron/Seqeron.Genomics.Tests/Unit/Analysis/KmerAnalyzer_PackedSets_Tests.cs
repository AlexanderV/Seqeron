using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// B09 WP78 R14-18 (F360): the packed k-mer sets (2 bits per base, sorted ulong arrays, k &lt;= 32) behind <see cref="KmerAnalyzer.JaccardSimilarity(string, string, int, KmerCountingOptions)"/>,
/// <see cref="KmerAnalyzer.MashDistance"/> and the containment indices give the numbers of the string sets (<see cref="KmerAnalyzer.DistinctKmers(string, int, KmerCountingOptions)"/>, which the packed engine
/// does not touch) on random sequences: literal, ACGT-only and canonical counting, lower case, N / IUPAC symbols, k = 1 ... 40 (above 32 the string route is used).
/// </summary>
[TestFixture]
[Category("GENOMIC-SIMILARITY-001")]
public class KmerAnalyzer_PackedSets_Tests
{
    private static string Random(Random rng, string alphabet, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    private static (int Shared, int Union) FromStringSets(string a, string b, int k, KmerCountingOptions options)
    {
        var setA = KmerAnalyzer.DistinctKmers(a, k, options);
        var setB = KmerAnalyzer.DistinctKmers(b, k, options);
        int shared = setA.Count(setB.Contains);
        return (shared, setA.Count + setB.Count - shared);
    }

    [Test]
    public void Jaccard_EqualsTheStringSets_OnRandomInputs()
    {
        var rng = new Random(7805);
        var optionSets = new[] { default(KmerCountingOptions), new KmerCountingOptions(AcgtOnly: true), new KmerCountingOptions(Canonical: true), new KmerCountingOptions(Canonical: true, AcgtOnly: true) };
        for (int round = 0; round < 800; round++)
        {
            string alphabet = rng.Next(3) == 0 ? "ACGT" : rng.Next(2) == 0 ? "ACGTacgt" : "ACGTNRYacgtn-U";
            string a = Random(rng, alphabet, rng.Next(0, 400));
            string b = rng.Next(3) == 0 ? Random(rng, alphabet, rng.Next(0, 400)) : a.Length > 20 ? string.Concat(a.AsSpan(rng.Next(0, 10)), Random(rng, "ACGT", rng.Next(0, 30))) : a;
            int k = rng.Next(1, 41);
            var options = optionSets[rng.Next(optionSets.Length)];
            var (shared, union) = FromStringSets(a, b, k, options);
            double expected = union == 0 ? 0.0 : (double)shared / union;
            Assert.That(KmerAnalyzer.JaccardSimilarity(a, b, k, options), Is.EqualTo(expected), $"k={k} {options} a={a} b={b}");
            double mash = shared == union ? 0.0 : KmerAnalyzer.MashDistanceFromJaccard((double)shared / union, k);
            Assert.That(KmerAnalyzer.MashDistance(a, b, k, options), Is.EqualTo(mash), $"mash k={k} {options}");
        }
    }

    [Test]
    public void Containment_EqualsTheStringSets_OnRandomInputs()
    {
        var rng = new Random(7806);
        for (int round = 0; round < 300; round++)
        {
            string a = Random(rng, "ACGT", rng.Next(0, 300));
            string b = Random(rng, "ACGT", rng.Next(0, 300)) + a.Take(rng.Next(0, 100));
            int k = rng.Next(1, 33);
            var options = new KmerCountingOptions(Canonical: rng.Next(2) == 0);
            var setA = KmerAnalyzer.DistinctKmers(a, k, options);
            var setB = KmerAnalyzer.DistinctKmers(b, k, options);
            double expected = setA.Count == 0 ? 0.0 : (double)setA.Count(setB.Contains) / setA.Count;
            Assert.That(KmerAnalyzer.ContainmentIndex(a, b, k, options), Is.EqualTo(expected), $"k={k} {options}");
        }
    }

    [Test]
    public void ReverseComplementScoresOneHundred_WithCanonicalKmers_AtTheWordLimits()
    {
        var rng = new Random(7807);
        foreach (int k in new[] { 1, 15, 16, 31, 32 })
        {
            string a = Random(rng, "ACGT", 500);
            string rc = Seqeron.Genomics.Core.DnaSequence.GetReverseComplementString(a);
            Assert.That(KmerAnalyzer.JaccardSimilarity(a, rc, k, new KmerCountingOptions(Canonical: true)), Is.EqualTo(1.0), $"k={k}");
        }
    }
}
