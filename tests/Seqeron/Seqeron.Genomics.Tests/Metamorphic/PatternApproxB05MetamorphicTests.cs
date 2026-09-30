namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic tests for the review 2026-09 B05 follow-up (PAT-APPROX-002/003).
///
/// Relations:
///   • BA1J-RC      FindFrequentKmersWithMismatchesAndReverseComplements(T) ==
///                  FindFrequentKmersWithMismatchesAndReverseComplements(revcomp(T)) as a (k-mer, score) set,
///                  since Count_d(rc(T), P) = Count_d(T, rc(P)) (Hamming distance is preserved by rc).
///   • BA1J-CLOSED  the BA1J result set is closed under reverse complement (score(P) = score(rc(P))).
///   • BA1J≥BA1I    the BA1J maximum is ≥ the BA1I maximum (score(P) = Count_d(P) + Count_d(rc(P)) ≥ Count_d(P)).
///   • ED-REVERSE   reversing both strings preserves Levenshtein distance (Myers engine) and the
///                  alignment distance (the reversed CIGAR is also optimal).
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class PatternApproxB05MetamorphicTests
{
    private static string RandomDna(Random rng, int length) =>
        new(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]).ToArray());

    private static HashSet<(string, int)> Ba1j(string text, int k, int d) =>
        ApproximateMatcher.FindFrequentKmersWithMismatchesAndReverseComplements(text, k, d).ToHashSet();

    [Test]
    public void Ba1j_TextAndReverseComplement_GiveSameSet()
    {
        var rng = new Random(20260930);
        for (int trial = 0; trial < 150; trial++)
        {
            string text = RandomDna(rng, rng.Next(1, 80));
            int k = rng.Next(1, 6);
            int d = rng.Next(0, 3);
            string rc = DnaSequence.GetReverseComplementString(text);

            Assert.That(Ba1j(rc, k, d).SetEquals(Ba1j(text, k, d)), Is.True, $"text={text} k={k} d={d}");
        }
    }

    [Test]
    public void Ba1j_ResultClosedUnderReverseComplement_AndDominatesBa1i()
    {
        var rng = new Random(20260931);
        for (int trial = 0; trial < 150; trial++)
        {
            string text = RandomDna(rng, rng.Next(4, 80));
            int k = rng.Next(1, 6);
            int d = rng.Next(0, 3);
            var set = Ba1j(text, k, d);

            foreach (var (kmer, score) in set)
                Assert.That(set.Contains((DnaSequence.GetReverseComplementString(kmer), score)), Is.True, kmer);

            var ba1i = ApproximateMatcher.FindFrequentKmersWithMismatches(text, k, d).ToList();
            if (ba1i.Count > 0)
                Assert.That(set.First().Item2, Is.GreaterThanOrEqualTo(ba1i[0].Count));
        }
    }

    [Test]
    public void EditDistance_ReversingBothStrings_PreservesDistance()
    {
        var rng = new Random(20260932);
        for (int trial = 0; trial < 500; trial++)
        {
            string a = RandomDna(rng, rng.Next(0, 150));
            string b = RandomDna(rng, rng.Next(0, 150));
            string ra = new(a.Reverse().ToArray());
            string rb = new(b.Reverse().ToArray());

            int d = ApproximateMatcher.EditDistance(a, b);
            Assert.That(ApproximateMatcher.EditDistance(ra, rb), Is.EqualTo(d));
            if (a.Length <= 40 && b.Length <= 40)
                Assert.That(ApproximateMatcher.GetEditAlignment(ra, rb).Distance, Is.EqualTo(d));
        }
    }
}
