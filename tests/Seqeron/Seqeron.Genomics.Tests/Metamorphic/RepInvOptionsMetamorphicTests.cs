namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic relations for the REP-INV-001 options (review B04, F21–F24): strand symmetry of
/// mismatch / arm-bounded stems, reverse symmetry of wobble stems, and case / N-prefix invariance of the
/// einverted-equivalent scored search. Fixed seeds; every relation is exact.
/// Test Unit: REP-INV-001.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Repeats")]
public class RepInvOptionsMetamorphicTests
{
    private static string Random(Random rng, string alphabet, int length) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static string ReverseComplement(string s) =>
        new(s.Reverse().Select(c => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c }).ToArray());

    private static List<(int, int, int, int)> Mirror(IEnumerable<InvertedRepeatResult> r, int n) =>
        r.Select(x => (n - x.RightArmStart - x.ArmLength, n - x.LeftArmStart - x.ArmLength, x.ArmLength, x.Mismatches))
         .OrderBy(t => t.Item1).ThenBy(t => t.Item2).ToList();

    private static List<(int, int, int, int)> Plain(IEnumerable<InvertedRepeatResult> r) =>
        r.Select(x => (x.LeftArmStart, x.RightArmStart, x.ArmLength, x.Mismatches)).ToList();

    /// <summary>MR: stems of revcomp(S) are the mirror images of the stems of S (Watson–Crick, any k / maxArm).</summary>
    [Test]
    public void MismatchStems_ReverseComplementMirror()
    {
        var rng = new Random(20260930);
        for (int t = 0; t < 300; t++)
        {
            string s = Random(rng, rng.Next(3) == 0 ? "AT" : "ACGT", rng.Next(0, 150));
            int k = rng.Next(0, 4), minLoop = rng.Next(0, 4);
            int maxArm = rng.Next(2) == 0 ? int.MaxValue : rng.Next(4, 12);
            var a = RepeatFinder.FindInvertedRepeats(s, 4, minLoop + 20, minLoop, k, maxArm);
            var b = RepeatFinder.FindInvertedRepeats(ReverseComplement(s), 4, minLoop + 20, minLoop, k, maxArm);
            Assert.That(Plain(b), Is.EqualTo(Mirror(a, s.Length)), $"{s} k={k} maxArm={maxArm}");
        }
    }

    /// <summary>MR: with wobble (G·U is not complement-symmetric) the stems of reverse(S) mirror those of S.</summary>
    [Test]
    public void WobbleStems_ReverseMirror()
    {
        var rng = new Random(7);
        for (int t = 0; t < 300; t++)
        {
            string s = Random(rng, rng.Next(2) == 0 ? "GTAC" : "GUA", rng.Next(0, 120));
            int k = rng.Next(0, 3);
            var a = RepeatFinder.FindInvertedRepeats(s, 3, 15, 2, k, allowWobble: true);
            var b = RepeatFinder.FindInvertedRepeats(new string(s.Reverse().ToArray()), 3, 15, 2, k, allowWobble: true);
            Assert.That(Plain(b), Is.EqualTo(Mirror(a, s.Length)), s);
        }
    }

    /// <summary>MR: scored repeats are case-invariant and shift by p under an N×p prefix (N never pairs).</summary>
    [Test]
    public void Scored_CaseInvariant_NPrefixShifts()
    {
        var rng = new Random(11);
        for (int t = 0; t < 120; t++)
        {
            string s = Unit.Analysis.RepeatFinder_InvertedRepeatOptions_Tests.PlantRepeat(rng, Random(rng, "ACGT", rng.Next(40, 400)));
            int maxRep = rng.Next(2) == 0 ? 2000 : rng.Next(20, 200), p = rng.Next(1, 50);
            var a = RepeatFinder.FindInvertedRepeatsScored(s, threshold: 30, maxRepeatLength: maxRep).ToList();
            var lower = RepeatFinder.FindInvertedRepeatsScored(s.ToLowerInvariant(), threshold: 30, maxRepeatLength: maxRep).ToList();
            var shifted = RepeatFinder.FindInvertedRepeatsScored(new string('N', p) + s, threshold: 30, maxRepeatLength: maxRep)
                .Select(r => r with
                {
                    LeftArmStart = r.LeftArmStart - p, LeftArmEnd = r.LeftArmEnd - p,
                    RightArmStart = r.RightArmStart - p, RightArmEnd = r.RightArmEnd - p
                }).ToList();
            Assert.Multiple(() =>
            {
                Assert.That(lower, Is.EqualTo(a), s);
                Assert.That(shifted, Is.EqualTo(a), $"p={p} maxRep={maxRep} {s}");
            });
        }
    }
}
