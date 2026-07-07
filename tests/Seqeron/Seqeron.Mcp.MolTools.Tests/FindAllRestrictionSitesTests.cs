using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class FindAllRestrictionSitesTests
{
    [Test]
    public void FindAllRestrictionSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.find_all_restriction_sites("AAAGAATTCAAA"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_all_restriction_sites(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_all_restriction_sites(null!));
    }

    [Test]
    public void FindAllRestrictionSites_Binding_InvokesSuccessfully()
    {
        // Scanning every enzyme must still find the EcoRI GAATTC site at position 3
        // (forward + reverse, since EcoRI is palindromic).
        var sites = MolToolsTools.find_all_restriction_sites("AAAGAATTCAAA").Sites;

        var ecoRIforward = sites.Where(s => s.Enzyme.Name == "EcoRI" && s.IsForwardStrand).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(ecoRIforward, Has.Count.EqualTo(1));
            Assert.That(ecoRIforward[0].Position, Is.EqualTo(3));
            Assert.That(ecoRIforward[0].CutPosition, Is.EqualTo(4));
            // Both strands reported for the palindromic site.
            Assert.That(sites.Count(s => s.Enzyme.Name == "EcoRI"), Is.EqualTo(2));
        });

        // A poly-A sequence matches no built-in recognition sequence -> no sites.
        Assert.That(MolToolsTools.find_all_restriction_sites("AAAAAAAA").Sites, Is.Empty);
    }
}
