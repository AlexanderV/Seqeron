using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Optional <c>methylation</c> parameter (none|dam|dcm|dam+dcm|emboss) of the site / digest / map tools.
/// Expected values are EMBOSS <c>restrict -methylation</c> output (see
/// RestrictionAnalyzer_Methylation_Tests): on
/// GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA EcoRI sits at 0 and 18, MboI at 10 and 36
/// unless Dam blocks it.
/// </summary>
[TestFixture]
public class RestrictionMethylationToolsTests
{
    private const string Seq = "GAATTCAAAAGATCAAAAGAATTCAAACCCGGGAAAGATCGAAATCGAAAA";

    [Test]
    public void FindRestrictionSites_DefaultAndNone_AreUnchanged()
    {
        var plain = MolToolsTools.find_restriction_sites(Seq, new[] { "MboI" }).Sites;
        var none = MolToolsTools.find_restriction_sites(Seq, new[] { "MboI" }, MoleculeTopology.Linear, "none").Sites;
        Assert.That(plain.Where(s => s.IsForwardStrand).Select(s => s.Position), Is.EqualTo(new[] { 10, 36 }));
        Assert.That(none, Is.EqualTo(plain));
    }

    [TestCase("dam")]
    [TestCase("DAM+dcm")]
    [TestCase("emboss")]
    public void FindRestrictionSites_DamBlocksMboIButNotEcoRI(string methylation)
    {
        Assert.That(MolToolsTools.find_restriction_sites(Seq, new[] { "MboI" }, MoleculeTopology.Linear, methylation).Sites, Is.Empty);
        var eco = MolToolsTools.find_restriction_sites(Seq, new[] { "EcoRI" }, MoleculeTopology.Linear, methylation).Sites;
        Assert.That(eco.Where(s => s.IsForwardStrand).Select(s => s.Position), Is.EqualTo(new[] { 0, 18 }));
    }

    [Test]
    public void FindRestrictionSites_DcmDoesNotBlockMboI() =>
        Assert.That(MolToolsTools.find_restriction_sites(Seq, new[] { "MboI" }, MoleculeTopology.Linear, "dcm").Sites, Has.Count.EqualTo(4));

    [Test]
    public void UnknownMethylationSet_IsAnError()
    {
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_restriction_sites(Seq, new[] { "MboI" }, MoleculeTopology.Linear, "cpg"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest(Seq, new[] { "MboI" }, MoleculeTopology.Linear, "cpg"));
    }

    [Test]
    public void RestrictionDigest_AndSummary_CarryMethylation()
    {
        var plain = MolToolsTools.restriction_digest(Seq, new[] { "EcoRI", "MboI" }).Fragments;
        var dam = MolToolsTools.restriction_digest(Seq, new[] { "EcoRI", "MboI" }, MoleculeTopology.Linear, "dam").Fragments;
        Assert.That(plain, Has.Count.GreaterThan(dam.Count));
        Assert.That(dam.Select(f => f.Length), Is.EqualTo(new[] { 1, 18, 32 }));

        var summary = MolToolsTools.digest_summary(Seq, new[] { "EcoRI", "MboI" }, MoleculeTopology.Linear, "dam");
        Assert.That(summary.FragmentSizes, Is.EqualTo(new[] { 32, 18, 1 }));
    }

    [Test]
    public void RestrictionMap_CarriesMethylation_AlsoForTheRebaseSet()
    {
        var map = MolToolsTools.restriction_map(Seq, new[] { "EcoRI", "MboI" }, MoleculeTopology.Linear, "common", "dam");
        Assert.That(map.NonCutters, Does.Contain("MboI"));
        Assert.That(map.SitesByEnzyme.ContainsKey("EcoRI"), Is.True);

        var plainAll = MolToolsTools.restriction_map(Seq, Array.Empty<string>(), MoleculeTopology.Linear, "rebase");
        var damAll = MolToolsTools.restriction_map(Seq, Array.Empty<string>(), MoleculeTopology.Linear, "rebase", "dam");
        Assert.That(plainAll.SitesByEnzyme.ContainsKey("MboI"), Is.True);
        Assert.That(damAll.SitesByEnzyme.ContainsKey("MboI"), Is.False);
        Assert.That(damAll.SitesByEnzyme.ContainsKey("EcoRI"), Is.True);
    }
}
