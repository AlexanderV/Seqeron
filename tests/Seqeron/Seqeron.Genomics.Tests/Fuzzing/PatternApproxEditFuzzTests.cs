namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the review 2026-09 B05 follow-up on PAT-APPROX-002 (Myers bit-parallel engine,
/// traceback, Damerau variants). Boundary exploitation: pattern lengths at the 64-bit word
/// boundaries (63/64/65/127/128/129), symbols outside ASCII (the Peq dictionary path), empty
/// strings, maxEdits up to int.MaxValue. Every input must give the DP-reference result and never an
/// unhandled runtime exception.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PatternApproxEditFuzzTests
{
    private static string RandomString(Random rng, int length, string alphabet) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    [TestCase(63)]
    [TestCase(64)]
    [TestCase(65)]
    [TestCase(127)]
    [TestCase(128)]
    [TestCase(129)]
    public void Myers_WordBoundaryLengths_EqualDp(int m)
    {
        var rng = new Random(m);
        foreach (string alphabet in new[] { "A", "AC", "ACGT", "ACGTé中" })
        {
            for (int trial = 0; trial < 20; trial++)
            {
                string a = RandomString(rng, m, alphabet);
                string b = RandomString(rng, rng.Next(0, 2 * m), alphabet);
                Assert.That(ApproximateMatcher.EditDistance(a, b), Is.EqualTo(ApproximateMatcher.EditDistanceDp(a, b)));
                Assert.That(ApproximateMatcher.EditDistance(b, a), Is.EqualTo(ApproximateMatcher.EditDistanceDp(b, a)));

                if (b.Length > 0)
                {
                    int k = rng.Next(0, 3) == 0 ? int.MaxValue : rng.Next(0, m + 1);
                    Assert.That(ApproximateMatcher.FindEditEndPositions(b, a, k).ToList(),
                        Is.EqualTo(ApproximateMatcher.FindEditEndPositionsDp(b, a, k).ToList()));
                }
            }
        }
    }

    [Test]
    public void DamerauAndAlignment_RandomInputs_NoException()
    {
        var rng = new Random(424242);
        for (int trial = 0; trial < 500; trial++)
        {
            string a = RandomString(rng, rng.Next(0, 20), "ABé￿");
            string b = RandomString(rng, rng.Next(0, 20), "ABé￿");
            int lev = ApproximateMatcher.EditDistance(a, b);
            Assert.That(ApproximateMatcher.GetEditAlignment(a, b).Distance, Is.EqualTo(lev));
            Assert.That(ApproximateMatcher.OptimalStringAlignmentDistance(a, b), Is.InRange(0, lev));
            Assert.That(ApproximateMatcher.DamerauLevenshteinDistance(a, b), Is.InRange(0, lev));
        }
    }

    [Test]
    public void FindWithEdits_ExtremeMaxEdits_ReportsEveryWindowWithAlignment()
    {
        string text = "ACGTNACGT";
        var hits = ApproximateMatcher.FindWithEdits(text, "GT", int.MaxValue).ToList();

        Assert.That(hits, Has.Count.EqualTo(text.Length * (text.Length + 1) / 2), "every non-empty window");
        Assert.That(hits.All(h => h.Alignment is not null && h.Alignment.Distance == h.Distance), Is.True);
    }

    [Test]
    public void Ba1j_BoundaryInputs_NoException()
    {
        foreach (var (text, k, d) in new[] { ("N", 1, 0), ("ACGT", 4, 4), ("ACGT", 5, 1), ("nnnnACGTnnnn", 3, 1), ("A", 1, 3) })
            Assert.DoesNotThrow(() => ApproximateMatcher.FindFrequentKmersWithMismatchesAndReverseComplements(text, k, d).ToList());
    }
}
