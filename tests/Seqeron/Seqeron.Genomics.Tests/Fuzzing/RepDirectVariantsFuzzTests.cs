namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the REP-DIRECT-001 enumeration variants (review B04, audit WP2): reverse-complement maximal
/// pairs, maximal k-mismatch repeats and supermaximal repeats never crash or hang on arbitrary characters and
/// extreme parameters, results satisfy their contracts, and invalid parameters fail eagerly.
/// Test Unit: REP-DIRECT-001.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Repeats")]
public class RepDirectVariantsFuzzTests
{
    private static string RandomText(Random rng, int n) =>
        new(Enumerable.Range(0, n).Select(_ => rng.Next(4) == 0 ? (char)rng.Next(0, 0x3000) : "ACGTNUacgtu-"[rng.Next(12)]).ToArray());

    [Test, CancelAfter(60000)]
    public void ReverseComplement_ArbitraryInput_NoCrash_ContractHolds()
    {
        var rng = new Random(11);
        int[] spacings = { int.MinValue, -50, -1, 0, 1, 7, int.MaxValue };
        for (int t = 0; t < 400; t++)
        {
            string s = RandomText(rng, rng.Next(0, 400));
            int minL = rng.Next(2, 8);
            int maxL = rng.Next(3) == 0 ? int.MaxValue : minL + rng.Next(0, 20);
            int minSp = spacings[rng.Next(spacings.Length)];
            foreach (var r in RepeatFinder.FindReverseComplementRepeats(s, minL, maxL, minSp))
            {
                Assert.That(r.Length, Is.InRange(minL, maxL));
                Assert.That(r.FirstPosition, Is.LessThanOrEqualTo(r.SecondPosition));
                Assert.That(r.Spacing, Is.GreaterThanOrEqualTo(minSp));
                Assert.That(r.SecondSequence, Is.EqualTo(new string(r.RepeatSequence.Reverse().Select(SequenceExtensions.GetComplementBase).ToArray())));
            }
        }
    }

    [Test, CancelAfter(60000)]
    public void KMismatch_ArbitraryInput_NoCrash_ContractHolds()
    {
        var rng = new Random(22);
        for (int t = 0; t < 300; t++)
        {
            string s = RandomText(rng, rng.Next(0, 300));
            int k = rng.Next(0, 5), minL = k + rng.Next(1, 12);
            if (minL < 2) minL = 2;
            int maxL = rng.Next(3) == 0 ? int.MaxValue : minL + rng.Next(0, 30);
            int minSp = rng.Next(3) == 0 ? int.MinValue : rng.Next(-20, 20);
            foreach (var r in RepeatFinder.FindApproximateDirectRepeats(s, minL, k, maxL, minSp, rng.Next(2) == 0))
            {
                Assert.That(r.Length, Is.InRange(minL, maxL));
                Assert.That(r.Mismatches, Is.InRange(0, k));
                Assert.That(r.FirstPosition, Is.LessThan(r.SecondPosition));
                Assert.That(r.Spacing, Is.GreaterThanOrEqualTo(minSp));
                Assert.That(r.FirstCopy, Has.Length.EqualTo(r.Length));
            }
        }
    }

    [Test, CancelAfter(60000)]
    public void Supermaximal_ArbitraryInput_NoCrash_ContractHolds()
    {
        var rng = new Random(33);
        for (int t = 0; t < 400; t++)
        {
            string s = RandomText(rng, rng.Next(0, 500));
            int minL = rng.Next(1, 10);
            foreach (var r in RepeatFinder.FindSupermaximalRepeats(s, minL))
            {
                Assert.That(r.Length, Is.GreaterThanOrEqualTo(minL));
                Assert.That(r.Positions, Has.Count.GreaterThanOrEqualTo(2));
                Assert.That(r.Positions, Is.Ordered.And.Unique);
                Assert.That(r.Sequence.All(c => "ACGT".Contains(c)));
            }
        }
    }

    [Test]
    public void Homopolymers_And_ExtremeParameters()
    {
        string a = new('A', 300);
        // Homopolymer A300: reverse-complement pairs need T, none; supermaximal = A299 at 0,1; k-mismatch k=0 = direct.
        Assert.That(RepeatFinder.FindReverseComplementRepeats(a, 2, int.MaxValue, int.MinValue), Is.Empty);
        var sup = RepeatFinder.FindSupermaximalRepeats(a, 1).Single();
        Assert.That((sup.Length, string.Join(",", sup.Positions)), Is.EqualTo((299, "0,1")));
        Assert.That(RepeatFinder.FindApproximateDirectRepeats(a, 5, 0, int.MaxValue, int.MinValue).Count(),
            Is.EqualTo(RepeatFinder.FindDirectRepeats(a, 5, int.MaxValue, int.MinValue).Count()));
        // (AT)150 is its own reverse complement: whole-sequence pair (0,0,300) present.
        string at = string.Concat(Enumerable.Repeat("AT", 150));
        Assert.That(RepeatFinder.FindReverseComplementRepeats(at, 2, int.MaxValue, int.MinValue)
            .Any(r => r.FirstPosition == 0 && r.SecondPosition == 0 && r.Length == 300));
    }

    [Test]
    public void InvalidParameters_ThrowEagerly()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindReverseComplementRepeats("ACGT", int.MinValue));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindReverseComplementRepeats("ACGT", 10, 9));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateDirectRepeats("ACGT", 3, int.MaxValue));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindApproximateDirectRepeats("ACGT", 3, int.MinValue));
        Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindSupermaximalRepeats("ACGT", int.MinValue));
    }
}
