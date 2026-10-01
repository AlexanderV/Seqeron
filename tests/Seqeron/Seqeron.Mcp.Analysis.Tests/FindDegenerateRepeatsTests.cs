using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_degenerate_repeats</c> MCP tool. Expected values = Vmatch 2.3.1 (-e / -h, -p, -allmax),
/// locked in RepeatFinder_DegenerateRepeats_Tests (Evidence REP-DIRECT-001 §WP8). NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindDegenerateRepeatsTests
{
    private static string Tuples(DegenerateRepeatItem[] r) =>
        string.Join(";", r.Select(x => $"{x.FirstPosition},{x.SecondPosition},{x.FirstLength},{x.SecondLength},{x.Distance}"));

    [Test]
    public void FindDegenerateRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindDegenerateRepeats("GATTACANGATTACA", 6));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindDegenerateRepeats(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindDegenerateRepeats(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindDegenerateRepeats("ACGTACGT", 4, 1, "levenshtein"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindDegenerateRepeats("ACGTACGT", 4, 0));
    }

    [Test]
    public void FindDegenerateRepeats_Binding_EditDirect_MatchesVmatchAllmax()
    {
        // vmatch -l 8 -e 1 -allmax
        var r = AnalysisTools.FindDegenerateRepeats("ACGTTGCATGCAAACGTAGCATGCAGGGTTTACGTTGCTTGCAAACG", 8, 1, "edit", false, int.MaxValue, int.MinValue).Items;
        Assert.That(Tuples(r), Is.EqualTo("0,13,12,12,1;0,31,16,16,1;5,18,8,8,1;8,39,9,8,1"));
        Assert.That((r[3].FirstCopy, r[3].SecondCopy, r[3].Spacing), Is.EqualTo(("TGCAAACGT", "TGCAAACG", 22)));
        Assert.That(r.All(x => !x.IsReverseComplement));
    }

    [Test]
    public void FindDegenerateRepeats_Palindromic_EditAndHamming_MatchVmatch()
    {
        // vmatch -p -l 8 -e 1 vs -p -l 8 -h 1 (-allmax)
        const string seq = "TTGACCGTAACCCCCGTTACGGTCAACC";
        Assert.That(Tuples(AnalysisTools.FindDegenerateRepeats(seq, 8, 1, "edit", true, int.MaxValue, int.MinValue).Items),
            Is.EqualTo("0,14,12,12,1;0,15,11,12,1;13,13,9,9,1"));
        Assert.That(Tuples(AnalysisTools.FindDegenerateRepeats(seq, 8, 1, "hamming", true, int.MaxValue, int.MinValue).Items),
            Is.EqualTo("0,14,12,12,1;13,13,9,9,1"));
    }
}
