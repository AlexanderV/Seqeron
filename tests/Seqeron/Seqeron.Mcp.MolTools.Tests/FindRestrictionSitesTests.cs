using NUnit.Framework;
using Seqeron.Genomics.MolTools;
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

    /// <summary>
    /// The topology parameter reaches the search. Reference (Biopython 1.88 Bio.Restriction):
    /// EcoRI.search(Seq("TTCAAAAAAAGAA"), linear=False) == [12] (1-based) while the linear
    /// search returns [] — the GAATTC site spans the origin starting at position 10.
    /// </summary>
    [Test]
    public void FindRestrictionSites_Circular_FindsOriginSpanningSite()
    {
        var linear = MolToolsTools.find_restriction_sites("TTCAAAAAAAGAA", new[] { "EcoRI" }).Sites;
        var circular = MolToolsTools
            .find_restriction_sites("TTCAAAAAAAGAA", new[] { "EcoRI" }, MoleculeTopology.Circular)
            .Sites;

        Assert.Multiple(() =>
        {
            Assert.That(linear, Is.Empty, "no GAATTC lies wholly inside the linear molecule");
            Assert.That(circular, Has.Count.EqualTo(2), "forward + reverse record of the one site");
            var fwd = circular.Single(s => s.IsForwardStrand);
            Assert.That(fwd.Position, Is.EqualTo(10));
            Assert.That(fwd.CutPosition, Is.EqualTo(11), "Biopython reports 12 (1-based)");
            Assert.That(fwd.RecognizedSequence, Is.EqualTo("GAATTC"));
        });
    }

    /// <summary>
    /// Linear cut feasibility is applied by default. Reference:
    /// Bio.Restriction.Sau3AI.search(Seq("GATCGATCGATC")) == [5] — only the middle GATC has both
    /// of its cuts inside the molecule.
    /// </summary>
    [Test]
    public void FindRestrictionSites_Linear_EdgeCutterSitesAreNotReported()
    {
        var sites = MolToolsTools.find_restriction_sites("GATCGATCGATC", new[] { "Sau3AI" }).Sites;

        Assert.Multiple(() =>
        {
            Assert.That(sites.Count(s => s.IsForwardStrand), Is.EqualTo(1));
            Assert.That(sites.Single(s => s.IsForwardStrand).CutPosition, Is.EqualTo(4));
        });
    }
}
