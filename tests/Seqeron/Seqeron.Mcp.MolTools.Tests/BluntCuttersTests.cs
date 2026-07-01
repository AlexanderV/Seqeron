using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class BluntCuttersTests
{
    [Test]
    public void BluntCutters_Schema_ValidatesCorrectly()
    {
        // Parameterless query tool: must never throw and must return a non-null result.
        Assert.DoesNotThrow(() => MolToolsTools.blunt_cutters());
        var result = MolToolsTools.blunt_cutters();
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Enzymes, Is.Not.Null);
    }

    [Test]
    public void BluntCutters_Binding_InvokesSuccessfully()
    {
        var result = MolToolsTools.blunt_cutters();
        var names = result.Enzymes.Select(e => e.Name).ToHashSet();

        Assert.Multiple(() =>
        {
            // The built-in database has exactly 10 blunt cutters (CutFwd == CutRev).
            Assert.That(result.Enzymes, Has.Count.EqualTo(10));

            // Every returned enzyme must actually be blunt.
            Assert.That(result.Enzymes.All(e => e.IsBluntEnd), Is.True);
            Assert.That(result.Enzymes.All(e => e.CutPositionForward == e.CutPositionReverse), Is.True);

            // Known blunt cutters present.
            foreach (var expected in new[]
                     { "AluI", "RsaI", "HaeIII", "DpnI", "EcoRV", "SmaI", "HincII", "ScaI", "StuI", "SwaI" })
                Assert.That(names, Does.Contain(expected), $"expected blunt cutter {expected}");

            // Sticky cutters must NOT appear.
            Assert.That(names, Does.Not.Contain("EcoRI"));
            Assert.That(names, Does.Not.Contain("BamHI"));

            // Spot-check a known record.
            var ecoRV = result.Enzymes.Single(e => e.Name == "EcoRV");
            Assert.That(ecoRV.RecognitionSequence, Is.EqualTo("GATATC"));
            Assert.That(ecoRV.CutPositionForward, Is.EqualTo(3));
            Assert.That(ecoRV.CutPositionReverse, Is.EqualTo(3));
        });
    }
}
