using NUnit.Framework;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class RemoveRestrictionSitesTests
{
    private static CodonUsageTableInput EColi() => new(Preset: "EColiK12");

    [Test]
    public void RemoveRestrictionSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.remove_restriction_sites("GAATTC", new[] { "GAATTC" }, EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.remove_restriction_sites("", new[] { "GAATTC" }, EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.remove_restriction_sites("GAATTC", Array.Empty<string>(), EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.remove_restriction_sites("GAATTC", null!, EColi()));
    }

    [Test]
    public void RemoveRestrictionSites_Binding_InvokesSuccessfully()
    {
        // "GAATTC" encodes Glu-Phe (GAA=E, UUC=F). Removing the EcoRI site synonymously
        // (GAA -> GAG etc.) must eliminate the RNA site "GAAUUC" while keeping the length.
        var result = MolToolsTools.remove_restriction_sites("GAATTC", new[] { "GAATTC" }, EColi());

        Assert.Multiple(() =>
        {
            Assert.That(result.OptimizedSequence, Does.Not.Contain("GAAUUC"));
            Assert.That(result.OptimizedSequence.Length, Is.EqualTo(6));
            // Output is RNA-alphabet (no T).
            Assert.That(result.OptimizedSequence, Does.Not.Contain("T"));
        });

        // A sequence without the site is returned unchanged (as RNA).
        var noSite = MolToolsTools.remove_restriction_sites("AUGAUG", new[] { "GAATTC" }, EColi());
        Assert.That(noSite.OptimizedSequence, Is.EqualTo("AUGAUG"));
    }
}
