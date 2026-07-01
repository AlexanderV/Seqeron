using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>hydrophobicity_profile</c> MCP tool.
/// Expected values computed from the Kyte-Doolittle scale (I=4.5, A=1.8, V=4.2) and
/// the window-mean definition, NOT the wrapper output.
/// </summary>
[TestFixture]
public class HydrophobicityProfileTests
{
    [Test]
    public void HydrophobicityProfile_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.HydrophobicityProfile("IIIII", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.HydrophobicityProfile("", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.HydrophobicityProfile(null!, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.HydrophobicityProfile("IIIII", 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.HydrophobicityProfile("IIIII", -1));
    }

    [Test]
    public void HydrophobicityProfile_Binding_InvokesSuccessfully()
    {
        // "IIIII" window 3 -> I=4.5 -> 3 windows all 4.5.
        var iso = AnalysisTools.HydrophobicityProfile("IIIII", 3).Values;
        Assert.Multiple(() =>
        {
            Assert.That(iso, Has.Length.EqualTo(3));
            Assert.That(iso, Is.All.EqualTo(4.5).Within(1e-9));
        });

        // "AV" window 2 -> (1.8 + 4.2)/2 = 3.0.
        var av = AnalysisTools.HydrophobicityProfile("AV", 2).Values;
        Assert.Multiple(() =>
        {
            Assert.That(av, Has.Length.EqualTo(1));
            Assert.That(av[0], Is.EqualTo(3.0).Within(1e-9));
        });

        // Window larger than sequence -> empty.
        var empty = AnalysisTools.HydrophobicityProfile("AV", 5).Values;
        Assert.That(empty, Is.Empty);
    }
}
