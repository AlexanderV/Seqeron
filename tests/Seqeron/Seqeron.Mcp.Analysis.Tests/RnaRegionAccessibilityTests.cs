using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>rna_region_accessibility</c> MCP tool (RNA-ACCESS-001). Expected values: ViennaRNA 2.7.2 exact constrained-ensemble ratio
/// exp(-(G_constrained - G_ensemble)/RT) with the region hard-constrained unpaired (NNDB-precision parameter file; equal to
/// RNA.pfl_fold_up(seq, u, n, n) at 37 C), NOT the wrapper output.
/// </summary>
[TestFixture]
public class RnaRegionAccessibilityTests
{
    private const string Seq = "GUCUCUGUGCGCAGCACAUGCUCGAGUCGUCGUGA";

    [Test]
    public void RnaRegionAccessibility_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.RnaRegionAccessibility(Seq, 34, 14));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaRegionAccessibility("", 0, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaRegionAccessibility(null!, 0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.RnaRegionAccessibility(Seq, 34, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.RnaRegionAccessibility(Seq, 35, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.RnaRegionAccessibility(Seq, 5, 20));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.RnaRegionAccessibility(Seq, 20, 5, temperature: 0));
    }

    [Test]
    public void RnaRegionAccessibility_Binding_MatchesViennaRna()
    {
        var r14 = AnalysisTools.RnaRegionAccessibility(Seq.ToLowerInvariant(), 34, 14);
        var r1 = AnalysisTools.RnaRegionAccessibility(Seq, 17, 1);
        Assert.Multiple(() =>
        {
            Assert.That(r14.UnpairedProbability, Is.EqualTo(0.030532).Within(0.0015));     // ViennaRNA constrained ensemble / pfl_fold_up
            Assert.That(r14.OpeningEnergy, Is.EqualTo(2.1503).Within(0.03));               // -RT ln P
            Assert.That(r14.WindowStart, Is.EqualTo(21));
            Assert.That(r14.WindowEnd, Is.EqualTo(34));
            Assert.That(r14.WindowLength, Is.EqualTo(14));
            Assert.That(r14.Temperature, Is.EqualTo(310.15));
            Assert.That(r1.UnpairedProbability, Is.EqualTo(0.143155).Within(0.0015));
        });
    }

    [Test]
    public void RnaRegionAccessibility_DelegatesToCanonicalMethods()
    {
        var r = AnalysisTools.RnaRegionAccessibility(Seq, 26, 14, minLoopSize: 4, temperature: 293.15);
        Assert.Multiple(() =>
        {
            Assert.That(r.UnpairedProbability, Is.EqualTo(RnaSecondaryStructure.CalculateRegionUnpairedProbability(Seq, 26, 14, 4, 293.15)));
            Assert.That(r.OpeningEnergy, Is.EqualTo(RnaSecondaryStructure.CalculateRegionOpeningEnergy(Seq, 26, 14, 4, 293.15)));
            Assert.That(r.Temperature, Is.EqualTo(293.15));
        });
    }

    [Test]
    public void RnaRegionAccessibility_NoPairPossible_ProbabilityOneEnergyZero()
    {
        var r = AnalysisTools.RnaRegionAccessibility("NNNNNNNN", 7, 4);
        Assert.Multiple(() =>
        {
            Assert.That(r.UnpairedProbability, Is.EqualTo(1.0));
            Assert.That(r.OpeningEnergy, Is.EqualTo(0.0).Within(1e-12));
        });
    }
}
