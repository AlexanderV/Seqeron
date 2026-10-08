using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>standardize_repeat_motif</c> MCP tool. Expected = misa.pl .statistics row names and Krait
/// motif.py StandardMotif(level).standard(), locked in RepeatFinder_MisaCompound_Tests. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class StandardizeRepeatMotifTests
{
    [Test]
    public void StandardizeRepeatMotif_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.StandardizeRepeatMotif("CA"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.StandardizeRepeatMotif(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.StandardizeRepeatMotif(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.StandardizeRepeatMotif("ACN"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.StandardizeRepeatMotif("AC", 5));
    }

    [TestCase("CTG", new[] { "CTG", "TGC", "AGC", "ACG", "ACG" }, "AGC/CTG")]
    [TestCase("ACAT", new[] { "ACAT", "ATAC", "ATAC", "ATAC", "ATAC" }, "ACAT/ATGT")]
    [TestCase("CTGGG", new[] { "CTGGG", "TGGGC", "AGCCC", "ACCCG", "ACCCG" }, "AGCCC/CTGGG")]
    public void StandardizeRepeatMotif_Binding_MatchesKraitAndMisa(string motif, string[] krait, string misa)
    {
        for (int level = 0; level <= 4; level++)
        {
            var r = AnalysisTools.StandardizeRepeatMotif(motif, level);
            Assert.That((r.StandardMotif, r.CanonicalClass, r.Level), Is.EqualTo((krait[level], misa, level)), $"level {level}");
        }
    }
}
