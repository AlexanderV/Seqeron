namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the REP-INV-001 options (maxMismatches / maxArmLength / allowWobble) and the einverted-equivalent
/// scored search (review B04, F21–F24): arbitrary characters and extreme parameters never crash or hang; invalid
/// parameters fail eagerly with <see cref="ArgumentOutOfRangeException"/>.
/// Test Unit: REP-INV-001.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Repeats")]
public class RepInvOptionsFuzzTests
{
    private static string RandomText(Random rng, int n) =>
        new(Enumerable.Range(0, n).Select(_ => rng.Next(4) == 0 ? (char)rng.Next(0, 0x3000) : "ACGTNUacgtu-"[rng.Next(12)]).ToArray());

    [Test, CancelAfter(60000)]
    public void Options_ArbitraryInput_NoCrash()
    {
        var rng = new Random(314);
        int[] loops = { 0, 1, 5, 50, int.MaxValue };
        for (int t = 0; t < 400; t++)
        {
            string s = RandomText(rng, rng.Next(0, 300));
            int minArm = rng.Next(2, 6), minLoop = rng.Next(0, 4);
            int maxLoop = Math.Max(minLoop, loops[rng.Next(loops.Length)]);
            int k = rng.Next(4) == 0 ? int.MaxValue : rng.Next(0, 5);
            int maxArm = rng.Next(3) == 0 ? int.MaxValue : minArm + rng.Next(0, 10);
            var res = RepeatFinder.FindInvertedRepeats(s, minArm, maxLoop, minLoop, k, maxArm, rng.Next(2) == 0).ToList();
            Assert.That(res.All(r => r.ArmLength >= minArm && r.LoopLength >= minLoop && r.Mismatches <= k), s);
        }
    }

    [Test, CancelAfter(60000)]
    public void Scored_ArbitraryInputAndParameters_NoCrash()
    {
        var rng = new Random(2718);
        for (int t = 0; t < 300; t++)
        {
            string s = RandomText(rng, rng.Next(0, 400));
            int gap = rng.Next(0, 30), threshold = rng.Next(0, 80), match = rng.Next(0, 10), mismatch = -rng.Next(0, 10);
            int maxRep = rng.Next(3) == 0 ? int.MaxValue : rng.Next(2, 500);
            var res = RepeatFinder.FindInvertedRepeatsScored(s, gap, threshold, match, mismatch, maxRep).ToList();
            Assert.That(res.All(r => r.LeftArmStart >= 0 && r.RightArmEnd < s.Length && r.LeftArmEnd < r.RightArmStart), s);
        }
    }

    [Test]
    public void InvalidParameters_ThrowEagerly()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeats("ACGT", maxMismatches: int.MinValue));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeats("ACGT", 10, maxArmLength: 9));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeatsScored("ACGT", gapPenalty: int.MaxValue));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeatsScored("ACGT", mismatchScore: int.MinValue));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeatsScored("ACGT", maxRepeatLength: int.MinValue));
        });
    }
}
