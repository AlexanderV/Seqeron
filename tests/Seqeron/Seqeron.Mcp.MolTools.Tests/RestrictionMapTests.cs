using NUnit.Framework;
using Seqeron.Genomics.MolTools;
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

    /// <summary>
    /// NonCutters is taken from the considered set, resolved case-insensitively, so a cutter is
    /// never also a non-cutter; with an empty set the whole built-in table is considered.
    /// </summary>
    [Test]
    public void RestrictionMap_NonCutters_PartitionTheConsideredSet()
    {
        var lower = MolToolsTools.restriction_map("AAAGAATTCAAA", new[] { "ecori", "bamhi" });
        var all = MolToolsTools.restriction_map("AAAGAATTCAAA", Array.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(lower.SitesByEnzyme.ContainsKey("EcoRI"), Is.True);
            Assert.That(lower.NonCutters, Is.EquivalentTo(new[] { "BamHI" }));
            Assert.That(all.NonCutters, Is.Not.Empty,
                "with no names the considered set is the whole table, so non-cutters are reported");
            Assert.That(all.NonCutters, Does.Not.Contain("EcoRI"));
            Assert.That(all.SitesByEnzyme.Count + all.NonCutters.Count,
                Is.EqualTo(RestrictionAnalyzer.Enzymes.Count));
        });
    }

    /// <summary>
    /// topology reaches CreateMap. Reference: Biopython 1.88
    /// EcoRI.search(Seq("TTCAAAAAAAGAA"), linear=False) == [12] (1-based), linear == [].
    /// </summary>
    [Test]
    public void RestrictionMap_Circular_FindsOriginSpanningSite()
    {
        var map = MolToolsTools.restriction_map(
            "TTCAAAAAAAGAA", new[] { "EcoRI" }, MoleculeTopology.Circular);

        Assert.Multiple(() =>
        {
            Assert.That(map.TotalSites, Is.EqualTo(1));
            Assert.That(map.NonCutters, Is.Empty);
            Assert.That(map.SitesByEnzyme["EcoRI"], Is.EqualTo(new[] { 10, 10 }));
        });
    }
}
