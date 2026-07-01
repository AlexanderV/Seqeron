using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class RestrictionMapTests
{
    [Test]
    public void RestrictionMap_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.restriction_map("AAAGAATTCAAA", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_map("", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_map(null!, new[] { "EcoRI" }));
        // Empty enzyme list is valid (means "all built-in enzymes").
        Assert.DoesNotThrow(() => MolToolsTools.restriction_map("AAAGAATTCAAA", Array.Empty<string>()));
    }

    [Test]
    public void RestrictionMap_Binding_InvokesSuccessfully()
    {
        // "AAAGAATTCAAA" mapped for EcoRI (present, single forward site) + BamHI (absent).
        var map = MolToolsTools.restriction_map("AAAGAATTCAAA", new[] { "EcoRI", "BamHI" });

        Assert.Multiple(() =>
        {
            Assert.That(map.SequenceLength, Is.EqualTo(12));
            // Both strands reported for the palindromic EcoRI site.
            Assert.That(map.SitesByEnzyme["EcoRI"], Is.EqualTo(new[] { 3, 3 }));
            // TotalSites counts forward strand only.
            Assert.That(map.TotalSites, Is.EqualTo(1));
            // EcoRI is a unique cutter (one forward-strand site).
            Assert.That(map.UniqueCutters, Does.Contain("EcoRI"));
            // BamHI does not cut this sequence.
            Assert.That(map.NonCutters, Does.Contain("BamHI"));
            Assert.That(map.NonCutters, Does.Not.Contain("EcoRI"));
        });
    }
}
