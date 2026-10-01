using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_supermaximal_repeats</c> MCP tool. Expected values = Vmatch 2.3.1 -supermax, locked in
/// RepeatFinder_RepeatVariants_Tests. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindSupermaximalRepeatsTests
{
    private static string Supers(SupermaximalRepeatItem[] r) =>
        string.Join(";", r.Select(x => $"{x.Length}:{string.Join(",", x.Positions)}"));

    [Test]
    public void FindSupermaximalRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindSupermaximalRepeats("GATTACAGATTACA", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindSupermaximalRepeats(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindSupermaximalRepeats(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindSupermaximalRepeats("ACGT", 0));
    }

    [TestCase("AAAAAAAACGTTGCAACGTAAAA", 2, "7:0,1;5:6,14")]
    [TestCase("GATTACAGATTACA", 3, "7:0,7")]
    [TestCase("ACGTACGTTTTTTTTTACGTACGT", 3, "8:0,16;8:7,8")]
    public void FindSupermaximalRepeats_Binding_MatchesVmatchSupermax(string seq, int minLength, string expected)
    {
        Assert.That(Supers(AnalysisTools.FindSupermaximalRepeats(seq, minLength).Items), Is.EqualTo(expected));
    }

    [Test]
    public void FindSupermaximalRepeats_SequenceAndOccurrences()
    {
        var r = AnalysisTools.FindSupermaximalRepeats("CAGCAGCAGTTTCAGCAG", 3).Items.Single();
        Assert.That((r.Sequence, r.Length, r.Positions), Is.EqualTo(("CAGCAG", 6, new[] { 0, 3, 12 })));
    }
}
