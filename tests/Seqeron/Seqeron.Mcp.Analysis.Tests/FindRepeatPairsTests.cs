using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// R14-17 (B09 F298): the <c>find_repeat_pairs</c> MCP tool (GenomicAnalyzer.FindRepeatPairs = MUMmer repeat-match list). Oracle: the installed MUMmer 3.23 <c>repeat-match -n 3</c> on ACGTTTACGAAAT prints
/// "1 7 3" (the direct pair, 1-based) and three reverse pairs of lengths 4, 3, 3 (<c>12r</c>, <c>9r</c>, <c>4r</c>).
/// </summary>
[TestFixture]
public class FindRepeatPairsTests
{
    private const string Seq = "ACGTTTACGAAAT";

    [Test]
    public void Forward_IsTheDirectPairOfRepeatMatch()
    {
        var r = AnalysisTools.FindRepeatPairs(Seq, 3, "Forward");
        Assert.That(r.Items, Has.Length.EqualTo(1));
        var p = r.Items[0];
        Assert.That((p.FirstPosition, p.SecondPosition, p.Length, p.Spacing, p.IsReverseComplement), Is.EqualTo((0, 6, 3, 3, false)));
    }

    [Test]
    public void Both_HasThreeReversePairsOfLengths_4_3_3_AndEqualsTheLibrary()
    {
        var r = AnalysisTools.FindRepeatPairs(Seq, 3);
        Assert.That(r.Items.Count(i => !i.IsReverseComplement), Is.EqualTo(1));
        Assert.That(r.Items.Where(i => i.IsReverseComplement).Select(i => i.Length).OrderByDescending(x => x), Is.EqualTo(new[] { 4, 3, 3 }));
        var lib = GenomicAnalyzer.FindRepeatPairs(Seq, 3);
        Assert.That(r.Items.Select(i => (i.FirstPosition, i.SecondPosition, i.Length, i.IsReverseComplement)), Is.EqualTo(lib.Select(p => (p.FirstPosition, p.SecondPosition, p.Length, p.IsReverseComplement))));
        var reverse = AnalysisTools.FindRepeatPairs(Seq, 3, "ReverseComplement");
        Assert.That(reverse.Items.All(i => i.IsReverseComplement), Is.True);
        Assert.That(reverse.Items, Has.Length.EqualTo(3));
    }

    [Test]
    public void IupacAndLowerCase_NeverMatch_AndEndAPair()
    {
        // an N is never part of a copy: ACN... cannot pair, so min 3 finds nothing and min 2 finds the AC pair that ends at the N; lower case is folded
        Assert.That(AnalysisTools.FindRepeatPairs("acntttacnaaat", 3, "Forward").Items, Is.Empty);
        var two = AnalysisTools.FindRepeatPairs("acntttacnaaat", 2, "Forward");
        // (0, 6, 2) = the AC copies that stop at the N (positions 2 and 8); the other pairs are the TT and AA runs; no copy covers an N
        var spans = two.Items.Select(i => (i.FirstPosition, i.SecondPosition, i.Length)).ToList();
        Assert.That(spans, Does.Contain((0, 6, 2)));
        Assert.That(spans, Is.EquivalentTo(new[] { (0, 6, 2), (3, 4, 2), (9, 10, 2) }));
        Assert.That(two.Items.Any(i => i.FirstPosition <= 2 && 2 < i.FirstPosition + i.Length || i.SecondPosition <= 8 && 8 < i.SecondPosition + i.Length), Is.False);
        Assert.That(AnalysisTools.FindRepeatPairs("acgtttacgaaat", 3, "Forward").Items, Has.Length.EqualTo(1));
    }

    [Test]
    public void Validation()
    {
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRepeatPairs("", 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindRepeatPairs(Seq, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRepeatPairs(Seq, 3, "Sideways"));
    }
}
