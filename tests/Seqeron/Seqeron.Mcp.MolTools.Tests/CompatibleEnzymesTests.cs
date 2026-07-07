using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CompatibleEnzymesTests
{
    [Test]
    public void CompatibleEnzymes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.compatible_enzymes());
        var result = MolToolsTools.compatible_enzymes();
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Pairs, Is.Not.Null);
    }

    [Test]
    public void CompatibleEnzymes_Binding_InvokesSuccessfully()
    {
        var result = MolToolsTools.compatible_enzymes();

        Assert.Multiple(() =>
        {
            // The 10 blunt cutters are pairwise compatible ("blunt"): C(10,2) = 45 blunt pairs.
            var bluntPairs = result.Pairs.Where(p => p.CompatibleEnd == "blunt").ToList();
            Assert.That(bluntPairs, Has.Count.EqualTo(45));

            // Every reported pair carries a non-empty compatible end and two distinct enzymes.
            Assert.That(result.Pairs.All(p => !string.IsNullOrEmpty(p.CompatibleEnd)), Is.True);
            Assert.That(result.Pairs.All(p => p.Enzyme1 != p.Enzyme2), Is.True);

            // BamHI (GGATCC,1,5) and BglII (AGATCT,1,5) both leave a 5' "GATC" overhang -> compatible.
            bool bamBgl = result.Pairs.Any(p =>
                p.CompatibleEnd == "GATC" &&
                new[] { p.Enzyme1, p.Enzyme2 }.Contains("BamHI") &&
                new[] { p.Enzyme1, p.Enzyme2 }.Contains("BglII"));
            Assert.That(bamBgl, Is.True, "BamHI/BglII should be GATC-compatible");
        });
    }
}
