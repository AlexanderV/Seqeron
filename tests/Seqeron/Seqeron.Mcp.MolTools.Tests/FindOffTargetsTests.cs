using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class FindOffTargetsTests
{
    private const string Guide = "ATATATATATATATATATAT"; // 20 nt (SpCas9)

    [Test]
    public void FindOffTargets_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.find_off_targets(Guide, "GTATATATATATATATATATTGG"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_off_targets("", "GTATATATATATATATATATTGG"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_off_targets(Guide, ""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_off_targets(Guide, "GTATATATATATATATATATTGG", 6));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_off_targets(Guide, "GTATATATATATATATATATTGG", -1));
        // Wrong guide length -> underlying method throws.
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_off_targets("ACGT", "GTATATATATATATATATATTGG"));
    }

    [Test]
    public void FindOffTargets_Binding_InvokesSuccessfully()
    {
        // Genome target "GTATATATATATATATATAT" differs from the guide at position 0 (A->G),
        // followed by a TGG PAM -> exactly one 1-mismatch off-target.
        var offTargets = MolToolsTools.find_off_targets(Guide, "GTATATATATATATATATATTGG", 3).OffTargets;

        Assert.Multiple(() =>
        {
            Assert.That(offTargets, Has.Count.EqualTo(1));
            Assert.That(offTargets[0].Mismatches, Is.EqualTo(1));
            Assert.That(offTargets[0].MismatchPositions, Is.EqualTo(new[] { 0 }));
            Assert.That(offTargets[0].Sequence, Is.EqualTo("GTATATATATATATATATAT"));
        });

        // The exact on-target (0 mismatches) is NOT reported as an off-target.
        var onlyOnTarget = MolToolsTools.find_off_targets(Guide, "ATATATATATATATATATATTGG", 3).OffTargets;
        Assert.That(onlyOnTarget, Is.Empty);
    }
}
