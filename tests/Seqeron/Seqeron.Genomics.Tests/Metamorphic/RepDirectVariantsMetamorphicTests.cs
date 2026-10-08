namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic relations for the REP-DIRECT-001 enumeration variants (review B04, audit WP2):
/// strand symmetry of reverse-complement maximal pairs, reversal / complement symmetry of k-mismatch repeats,
/// reversal symmetry of supermaximal repeats, and relabelling / case invariance. Fixed seeds; every relation is exact.
/// Test Unit: REP-DIRECT-001.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Repeats")]
public class RepDirectVariantsMetamorphicTests
{
    private static string Random(Random rng, string alphabet, int length) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static string Map(string s, string from, string to) =>
        new(s.Select(c => from.IndexOf(c) is var i and >= 0 ? to[i] : c).ToArray());

    private static string ReverseComplement(string s) => Map(new string(s.Reverse().ToArray()), "ACGT", "TGCA");

    /// <summary>MR: pairs of revcomp(S) are the mirror images (n−k−L, n−i−L, L) of the pairs of S.</summary>
    [Test]
    public void ReverseComplementRepeats_StrandMirror()
    {
        var rng = new Random(101);
        for (int t = 0; t < 60; t++)
        {
            string s = Random(rng, t % 3 == 0 ? "ACGTN" : "ACGT", rng.Next(20, 600));
            int n = s.Length;
            var mirrored = RepeatFinder.FindReverseComplementRepeats(s, 4, int.MaxValue, int.MinValue)
                .Select(r => (n - r.SecondPosition - r.Length, n - r.FirstPosition - r.Length, r.Length))
                .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item3).ToList();
            var direct = RepeatFinder.FindReverseComplementRepeats(ReverseComplement(s), 4, int.MaxValue, int.MinValue)
                .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
            Assert.That(direct, Is.EqualTo(mirrored), s);
        }
    }

    /// <summary>MR: a complement-preserving relabelling (A↔C, T↔G) and lowercasing leave the pair set unchanged.</summary>
    [Test]
    public void ReverseComplementRepeats_ComplementPreservingRelabel_AndCase_Invariant()
    {
        var rng = new Random(202);
        for (int t = 0; t < 40; t++)
        {
            string s = Random(rng, "ACGT", rng.Next(10, 500));
            var baseline = RepeatFinder.FindReverseComplementRepeats(s, 3, 60, -10).Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
            Assert.That(RepeatFinder.FindReverseComplementRepeats(Map(s, "ACGT", "CATG"), 3, 60, -10)
                .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)), Is.EqualTo(baseline));
            Assert.That(RepeatFinder.FindReverseComplementRepeats(s.ToLowerInvariant(), 3, 60, -10)
                .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)), Is.EqualTo(baseline));
        }
    }

    /// <summary>MR: reversing S maps each k-mismatch repeat (i, j, L, m) to (n−j−L, n−i−L, L, m); complementing S keeps it.</summary>
    [Test]
    public void KMismatchRepeats_ReverseMirror_ComplementInvariant()
    {
        var rng = new Random(303);
        for (int t = 0; t < 50; t++)
        {
            string s = Random(rng, t % 4 == 0 ? "ACGTN" : "ACGT", rng.Next(20, 400));
            int n = s.Length, k = rng.Next(0, 3), minL = rng.Next(k + 4, k + 12);
            bool literal = t % 2 == 0;
            var forward = RepeatFinder.FindApproximateDirectRepeats(s, minL, k, int.MaxValue, int.MinValue, literal).ToList();
            var mirrored = forward.Select(r => (n - r.SecondPosition - r.Length, n - r.FirstPosition - r.Length, r.Length, r.Mismatches))
                .OrderBy(x => x.Item1).ThenBy(x => x.Item2).ThenBy(x => x.Item3).ToList();
            var reversed = RepeatFinder.FindApproximateDirectRepeats(new string(s.Reverse().ToArray()), minL, k, int.MaxValue, int.MinValue, literal)
                .Select(r => (r.FirstPosition, r.SecondPosition, r.Length, r.Mismatches)).ToList();
            Assert.That(reversed, Is.EqualTo(mirrored), s);
            Assert.That(RepeatFinder.FindApproximateDirectRepeats(Map(s, "ACGT", "TGCA"), minL, k, int.MaxValue, int.MinValue, literal)
                .Select(r => (r.FirstPosition, r.SecondPosition, r.Length, r.Mismatches)),
                Is.EqualTo(forward.Select(r => (r.FirstPosition, r.SecondPosition, r.Length, r.Mismatches))));
        }
    }

    /// <summary>MR: increasing k never loses coverage — every (k)-repeat lies inside some (k+1)-repeat on its diagonal.</summary>
    [Test]
    public void KMismatchRepeats_LargerK_CoversSmallerK()
    {
        var rng = new Random(404);
        for (int t = 0; t < 30; t++)
        {
            string s = Random(rng, "ACGT", rng.Next(50, 400));
            int k = rng.Next(0, 3), minL = k + 6;
            var small = RepeatFinder.FindApproximateDirectRepeats(s, minL, k, int.MaxValue, int.MinValue).ToList();
            var large = RepeatFinder.FindApproximateDirectRepeats(s, minL, k + 1, int.MaxValue, int.MinValue).ToList();
            foreach (var r in small)
            {
                int d = r.SecondPosition - r.FirstPosition;
                Assert.That(large.Any(q => q.SecondPosition - q.FirstPosition == d && q.FirstPosition <= r.FirstPosition
                                           && q.FirstPosition + q.Length >= r.FirstPosition + r.Length), $"{s} {r}");
            }
        }
    }

    /// <summary>MR: supermaximal repeats of reverse(S) are the reversed strings at mirrored positions; relabelling/case invariant.</summary>
    [Test]
    public void SupermaximalRepeats_ReverseMirror_RelabelInvariant()
    {
        var rng = new Random(505);
        for (int t = 0; t < 60; t++)
        {
            string s = Random(rng, t % 3 == 0 ? "ACGTN" : "ACGT", rng.Next(10, 800));
            int n = s.Length, minL = rng.Next(1, 8);
            var forward = RepeatFinder.FindSupermaximalRepeats(s, minL).ToList();
            var mirrored = forward.Select(r => (new string(r.Sequence.Reverse().ToArray()),
                    string.Join(",", r.Positions.Select(p => n - p - r.Length).OrderBy(p => p))))
                .OrderBy(x => x.Item2).ThenBy(x => x.Item1).ToList();
            var reversed = RepeatFinder.FindSupermaximalRepeats(new string(s.Reverse().ToArray()), minL)
                .Select(r => (r.Sequence, string.Join(",", r.Positions))).OrderBy(x => x.Item2).ThenBy(x => x.Item1).ToList();
            Assert.That(reversed, Is.EqualTo(mirrored), s);
            Assert.That(RepeatFinder.FindSupermaximalRepeats(Map(s, "ACGT", "GTAC").ToLowerInvariant(), minL)
                .Select(r => (r.Length, string.Join(",", r.Positions))),
                Is.EqualTo(forward.Select(r => (r.Length, string.Join(",", r.Positions)))));
        }
    }
}
