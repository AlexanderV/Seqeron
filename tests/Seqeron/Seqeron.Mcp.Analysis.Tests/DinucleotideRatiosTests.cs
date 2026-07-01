using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>dinucleotide_ratios</c> MCP tool.
/// Expected values computed by hand from SequenceStatistics.CalculateDinucleotideRatios:
/// rho_XY = f_XY / (f_X * f_Y). NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class DinucleotideRatiosTests
{
    [Test]
    public void DinucleotideRatios_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.DinucleotideRatios("ATAT"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DinucleotideRatios(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DinucleotideRatios(null!));
    }

    [Test]
    public void DinucleotideRatios_Binding_InvokesSuccessfully()
    {
        // "ATAT": comp A=2,T=2,total=4 -> f_A=f_T=0.5. dinucFreq AT=2/3, TA=1/3.
        //   rho_AT = (2/3) / (0.5*0.5) = (2/3)/0.25 = 8/3 = 2.6667
        //   rho_TA = (1/3) / 0.25 = 4/3 = 1.3333
        var r = AnalysisTools.DinucleotideRatios("ATAT").Ratios;
        Assert.Multiple(() =>
        {
            Assert.That(r["AT"], Is.EqualTo(8.0 / 3.0).Within(1e-12));
            Assert.That(r["TA"], Is.EqualTo(4.0 / 3.0).Within(1e-12));
        });
    }
}
