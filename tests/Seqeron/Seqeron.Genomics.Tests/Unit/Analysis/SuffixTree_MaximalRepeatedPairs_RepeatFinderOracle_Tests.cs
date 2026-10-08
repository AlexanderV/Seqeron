using System;
using System.Linq;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Cross-check of <c>SuffixTree.FindMaximalRepeatedPairs</c> (B05, suffix tree, Gusfield 1997 §7.12.3)
/// against the independent B04 enumeration behind <see cref="RepeatFinder.FindDirectRepeats(string, int, int, int)"/>
/// (suffix array + LCP intervals, non-ACGT symbols unique). With maxLength = int.MaxValue and
/// minSpacing = int.MinValue, FindDirectRepeats returns every forward maximal pair of an upper-cased
/// sequence ordered by (FirstPosition, SecondPosition) — the same set and order as the suffix tree with
/// non-ACGT characters declared unique.
/// </summary>
[TestFixture]
public class SuffixTree_MaximalRepeatedPairs_RepeatFinderOracle_Tests
{
    private static (int, int, int)[] Tree(string seq, int minLength)
        => global::SuffixTree.SuffixTree.Build(seq)
            .FindMaximalRepeatedPairs(minLength, static c => c is not ('A' or 'C' or 'G' or 'T'))
            .Select(p => (p.FirstPosition, p.SecondPosition, p.Length)).ToArray();

    private static (int, int, int)[] Oracle(string seq, int minLength)
        => RepeatFinder.FindDirectRepeats(seq, minLength, int.MaxValue, int.MinValue)
            .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToArray();

    [Test]
    public void LiteralWithN_EqualsRepeatFinder()
    {
        const string seq = "ACGTNACGTNACGT";
        Assert.That(Tree(seq, 3), Is.EqualTo(new[] { (0, 5, 4), (0, 10, 4), (5, 10, 4) }));
        Assert.That(Oracle(seq, 3), Is.EqualTo(Tree(seq, 3)));
    }

    [Test]
    public void RandomSequencesWithAmbiguityCodes_EqualRepeatFinder([Values(0, 1, 2)] int seed)
    {
        var rng = new Random(20261004 + seed);
        for (int iter = 0; iter < 40; iter++)
        {
            int n = rng.Next(0, 600);
            var chars = new char[n];
            for (int i = 0; i < n; i++)
                chars[i] = rng.Next(25) == 0 ? "NRY"[rng.Next(3)] : "ACGT"[rng.Next(4)];
            string seq = new(chars);
            int minLength = rng.Next(2, 8);
            Assert.That(Tree(seq, minLength), Is.EqualTo(Oracle(seq, minLength)), $"seed {seed} iter {iter} L={minLength}");
        }
    }
}
