using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_protein_low_complexity_regions</c> MCP tool.
/// Expected values from ProteinMotifFinder's own SEG unit test
/// (ProteinMotifFinder_FindLowComplexityRegions_Tests, flank+Q20+flank -> region 6..37,
/// complexity 0.0; 20 distinct residues -> none), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindProteinLowComplexityRegionsTests
{
    [Test]
    public void FindProteinLowComplexityRegions_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindProteinLowComplexityRegions("ACDEFGHIKLMN" + new string('Q', 20) + "NMLKIHGFEDCA"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindProteinLowComplexityRegions(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindProteinLowComplexityRegions(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindProteinLowComplexityRegions("QQQQ", 0));
    }

    [Test]
    public void FindProteinLowComplexityRegions_Binding_InvokesSuccessfully()
    {
        // flank[12] + Q20 + flank[12]: SEG extends the poly-Q region to residues 6..37.
        string sequence = "ACDEFGHIKLMN" + new string('Q', 20) + "NMLKIHGFEDCA";
        var regions = AnalysisTools.FindProteinLowComplexityRegions(sequence).Items;
        Assert.Multiple(() =>
        {
            Assert.That(regions, Has.Length.EqualTo(1));
            Assert.That(regions[0].Start, Is.EqualTo(6));
            Assert.That(regions[0].End, Is.EqualTo(37));
            Assert.That(regions[0].Complexity, Is.EqualTo(0.0).Within(1e-10));
        });

        // All 20 distinct residues -> maximal complexity, no low-complexity region.
        var none = AnalysisTools.FindProteinLowComplexityRegions("ACDEFGHIKLMNPQRSTVWY").Items;
        Assert.That(none, Is.Empty);
    }
}
