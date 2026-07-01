using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class FindRestrictionSitesTests
{
    [Test]
    public void FindRestrictionSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.find_restriction_sites("AAAGAATTCAAA", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_restriction_sites("", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_restriction_sites(null!, new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_restriction_sites("AAAGAATTCAAA", Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_restriction_sites("AAAGAATTCAAA", null!));
    }

    [Test]
    public void FindRestrictionSites_Binding_InvokesSuccessfully()
    {
        // "AAAGAATTCAAA": EcoRI (GAATTC) occurs once at index 3.
        // Palindromic enzyme -> one forward + one reverse site, both at position 3.
        var sites = MolToolsTools.find_restriction_sites("AAAGAATTCAAA", new[] { "EcoRI" }).Sites;

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(2));

            var fwd = sites.Single(s => s.IsForwardStrand);
            Assert.That(fwd.Position, Is.EqualTo(3));
            Assert.That(fwd.RecognizedSequence, Is.EqualTo("GAATTC"));
            // Cut position = position + CutPositionForward(1) = 4.
            Assert.That(fwd.CutPosition, Is.EqualTo(4));

            var rev = sites.Single(s => !s.IsForwardStrand);
            Assert.That(rev.Position, Is.EqualTo(3));
            // Reverse cut = forwardPos + CutPositionReverse(5) = 8.
            Assert.That(rev.CutPosition, Is.EqualTo(8));
        });

        // A sequence with no EcoRI site returns no sites.
        Assert.That(MolToolsTools.find_restriction_sites("AAAAAAAA", new[] { "EcoRI" }).Sites, Is.Empty);
    }
}
