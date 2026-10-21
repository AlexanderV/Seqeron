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

    // ---- B08 round 3 G3 (F83): partial_digest. Oracle: constructed set of F53 over the cuts EMBOSS
    // `restrict -methylation` keeps (L2: EcoRI 1,19; MboI 10,36 blocked; SmaI 30; TaqI 39,45 -> X = {0,1,19,30,39,45,51}).
    private static readonly string[] PartialEnzymes = { "EcoRI", "MboI", "SmaI", "TaqI" };

    [Test]
    public void PartialDigest_DefaultAndNone_AreUnchanged()
    {
        var plain = MolToolsTools.partial_digest(Seq, PartialEnzymes).Fragments;
        var none = MolToolsTools.partial_digest(Seq, PartialEnzymes, MoleculeTopology.Linear, null, "none").Fragments;
        Assert.That(plain, Has.Count.EqualTo(36));      // 7 cuts -> 8 * 9 / 2
        Assert.That(none, Is.EqualTo(plain));
    }

    [TestCase("dam")]
    [TestCase("DAM+dcm")]
    [TestCase("emboss")]
    public void PartialDigest_Methylation_UsesTheSurvivingCutSet(string methylation)
    {
        var fragments = MolToolsTools.partial_digest(Seq, PartialEnzymes, MoleculeTopology.Linear, null, methylation).Fragments;
        Assert.That(fragments.Select(f => f.Fragment.Length).OrderBy(x => x),
            Is.EqualTo(new[] { 1, 6, 6, 9, 11, 12, 15, 18, 19, 20, 21, 26, 29, 30, 32, 38, 39, 44, 45, 50, 51 }));
        Assert.That(fragments.Max(f => f.MissedCuts), Is.EqualTo(5));
        // missedCuts = 0 prefix is restriction_digest with the same methylation, member for member
        Assert.That(fragments.Where(f => f.MissedCuts == 0).Select(f => f.Fragment),
            Is.EqualTo(MolToolsTools.restriction_digest(Seq, PartialEnzymes, MoleculeTopology.Linear, methylation).Fragments));
    }

    [Test]
    public void PartialDigest_Methylation_RespectsTopologyAndLimit()
    {
        const string l3 = "TTTTAACGGTACCGTGCTTTTGAATTCTTTTGCACGGTACCGTTTTTGGTACCTTT";
        var linear = MolToolsTools.partial_digest(l3, new[] { "KpnI", "Acc65I", "EcoRI" }, MoleculeTopology.Linear, null, "emboss").Fragments;
        Assert.That(linear.Select(f => (f.Fragment.Length, f.MissedCuts)).OrderBy(t => t),
            Is.EqualTo(new[] { (4, 0), (4, 0), (8, 1), (22, 0), (26, 0), (30, 1), (34, 2), (48, 1), (52, 2), (56, 3) }));

        var circular = MolToolsTools.partial_digest(l3, new[] { "KpnI", "Acc65I", "EcoRI" }, MoleculeTopology.Circular, 1, "emboss").Fragments;
        Assert.That(circular, Has.Count.EqualTo(6));    // 3 cuts, <= 1 missed: 3 * 2
        Assert.That(circular.Max(f => f.MissedCuts), Is.EqualTo(1));
    }

    [Test]
    public void PartialDigest_UnknownMethylationSet_IsAnError() =>
        Assert.Throws<ArgumentException>(() => MolToolsTools.partial_digest(Seq, PartialEnzymes, MoleculeTopology.Linear, null, "cpg"));
}
