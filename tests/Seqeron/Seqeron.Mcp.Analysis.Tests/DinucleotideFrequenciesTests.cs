using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>dinucleotide_frequencies</c> MCP tool.
/// Expected values computed by hand from SequenceStatistics.CalculateDinucleotideFrequencies:
/// f_XY = count(XY) / (N-1 counted dinucleotides). NOT from the wrapper's output.
/// </summary>
[TestFixture]
public class DinucleotideFrequenciesTests
{
    [Test]
    public void DinucleotideFrequencies_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.DinucleotideFrequencies("ATAT"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DinucleotideFrequencies(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.DinucleotideFrequencies(null!));
    }

    [Test]
    public void DinucleotideFrequencies_Binding_InvokesSuccessfully()
    {
        // "AAAA": dinucs AA,AA,AA -> total 3 -> AA = 1.0.
        var homo = AnalysisTools.DinucleotideFrequencies("AAAA").Frequencies;
        Assert.Multiple(() =>
        {
            Assert.That(homo, Has.Count.EqualTo(1));
            Assert.That(homo["AA"], Is.EqualTo(1.0).Within(1e-12));
        });

        // "ATAT": AT,TA,AT -> AT=2/3, TA=1/3.
        var alt = AnalysisTools.DinucleotideFrequencies("ATAT").Frequencies;
        Assert.Multiple(() =>
        {
            Assert.That(alt["AT"], Is.EqualTo(2.0 / 3.0).Within(1e-12));
            Assert.That(alt["TA"], Is.EqualTo(1.0 / 3.0).Within(1e-12));
            Assert.That(alt.Values.Sum(), Is.EqualTo(1.0).Within(1e-12));
        });
    }
}
