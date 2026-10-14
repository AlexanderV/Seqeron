using NUnit.Framework;
using Seqeron.Genomics.MolTools;
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

    /// <summary>
    /// topology reaches the all-enzyme scan. Reference: Biopython 1.88
    /// EcoRI.search(Seq("TTCAAAAAAAGAA"), linear=False) == [12] (1-based), linear == [].
    /// </summary>
    [Test]
    public void FindAllRestrictionSites_Circular_FindsOriginSpanningSite()
    {
        var circular = MolToolsTools
            .find_all_restriction_sites("TTCAAAAAAAGAA", MoleculeTopology.Circular)
            .Sites
            .Where(s => s.IsForwardStrand && s.Enzyme.Name == "EcoRI")
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.find_all_restriction_sites("TTCAAAAAAAGAA").Sites
                .Any(s => s.Enzyme.Name == "EcoRI"), Is.False, "the linear scan finds no EcoRI site");
            Assert.That(circular, Has.Count.EqualTo(1));
            Assert.That(circular[0].CutPosition, Is.EqualTo(11));
        });
    }
}
