using NUnit.Framework;
using Seqeron.Mcp.Population.Tools;

namespace Seqeron.Mcp.Population.Tests;

[TestFixture]
public class ScanSelectionSignalsTests
{
    private static RegionItem R(string id, double tajD, double fst, double ihs) =>
        new(id, 0, 10_000, tajD, fst, ihs);

    [Test]
    public void ScanSelectionSignals_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PopulationTools.ScanSelectionSignals(new[] { R("Region1", -2.5, 0.1, 0.5) }));
        Assert.DoesNotThrow(() => PopulationTools.ScanSelectionSignals(Array.Empty<RegionItem>()));
    }

    [Test]
    public void ScanSelectionSignals_Binding_InvokesSuccessfully()
    {
        // Tajima's D = -2.5 below threshold -2.0 → emits a TajimasD signal.
        var tajD = PopulationTools.ScanSelectionSignals(
            new[] { R("Region1", -2.5, 0.1, 0.5) }, tajimaDThreshold: -2.0);
        Assert.Multiple(() =>
        {
            Assert.That(tajD.Items, Has.Count.GreaterThanOrEqualTo(1));
            Assert.That(tajD.Items.Any(s => s.TestType == "TajimasD"), Is.True);
            Assert.That(tajD.Items.First(s => s.TestType == "TajimasD").Score, Is.EqualTo(-2.5).Within(1e-10));
        });

        // High Fst = 0.5 above threshold 0.25 → emits an Fst signal.
        var highFst = PopulationTools.ScanSelectionSignals(
            new[] { R("Region1", 0.0, 0.5, 0.0) }, fstThreshold: 0.25);
        Assert.That(highFst.Items.Any(s => s.TestType == "Fst"), Is.True);
    }

    [Test]
    public void ScanSelectionSignals_Binding_NeutralRegion_EmitsNothing()
    {
        var neutral = PopulationTools.ScanSelectionSignals(
            new[] { R("Region1", 0.0, 0.1, 0.5) });
        Assert.That(neutral.Items, Is.Empty);
    }
}
