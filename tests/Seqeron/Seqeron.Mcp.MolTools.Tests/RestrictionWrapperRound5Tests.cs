using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Findings of audit round 5 on the restriction MCP wrappers (B08 A5R-1, A5R-4, A5R-5 and the surviving wrapper mutants
/// A5R-T6, A5R-T7). The expected numbers come from the documented geometry (BaeI cuts 10 before the site and 23 after it,
/// REBASE emboss 404 via Biopython 1.88), from EMBOSS' Dam semantics, and from the table itself.
/// </summary>
[TestFixture]
public class RestrictionWrapperRound5Tests
{
    // ------------------------------------------------------------------------------------------------
    // A5R-1: the partner may be spelled the way GetEnzyme accepts (trimmed, REBASE hyphen)
    // ------------------------------------------------------------------------------------------------

    [TestCase("MspI", " HpaII", "CG", "CG")]
    [TestCase("MspI", "hpaii", "CG", "CG")]
    [TestCase("AluI", "CviKI-1", "blunt", "blunt")]
    [TestCase("AluI", "  cviki_1 ", "blunt", "blunt")]
    [TestCase("CviKI-1", "AluI", "blunt", "blunt")]
    public void EnzymesCompatible_AcceptsEverySpellingTheLibraryResolves(string first, string second, string end, string end2)
    {
        var result = MolToolsTools.enzymes_compatible(first, second);

        Assert.Multiple(() =>
        {
            Assert.That(result.Compatible, Is.True);
            Assert.That(result.CompatibleEnd, Is.EqualTo(end));
            Assert.That(result.End2, Is.EqualTo(end2));
        });
    }

    [Test]
    public void EnzymesCompatible_IncompatiblePairAndUnknownName_AreUnchanged()
    {
        var incompatible = MolToolsTools.enzymes_compatible("EcoRI", "EcoRV");

        Assert.Multiple(() =>
        {
            Assert.That(incompatible.Compatible, Is.False);
            Assert.That(incompatible.CompatibleEnd, Is.Null);
            Assert.That(incompatible.End2, Is.Null);
            Assert.Throws<ArgumentException>(() => MolToolsTools.enzymes_compatible("EcoRI", "NoSuchEnzyme1"));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-4: the custom-enzyme tools refuse an empty sequence like the six named ones
    // ------------------------------------------------------------------------------------------------

    [TestCase(null)]
    [TestCase("")]
    public void CustomEnzymeTools_RefuseANullOrEmptySequence(string? sequence)
    {
        var find = Assert.Throws<ArgumentException>(() => MolToolsTools.find_sites_custom_enzyme(sequence!, "GAATTC", 1, 5));
        var digest = Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest_custom_enzyme(sequence!, "GAATTC", 1, 5));
        Assert.Multiple(() =>
        {
            Assert.That(find!.ParamName, Is.EqualTo("sequence"));
            Assert.That(digest!.ParamName, Is.EqualTo("sequence"));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-5: an undefined topology (any number deserializes into the enum) is refused by every restriction tool
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void RestrictionTools_RefuseAnUndefinedTopology()
    {
        const string seq = "AAAGAATTCAAA";
        const MoleculeTopology bad = (MoleculeTopology)7;
        var names = new[] { "EcoRI" };

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.find_restriction_sites(seq, names, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.find_all_restriction_sites(seq, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.restriction_digest(seq, names, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.partial_digest(seq, names, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.digest_summary(seq, names, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.restriction_map(seq, names, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.find_sites_custom_enzyme(seq, "GAATTC", 1, 5, null, bad));
            Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.restriction_digest_custom_enzyme(seq, "GAATTC", 1, 5, null, bad));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-T6: every argument of the two custom-enzyme tools reaches the algorithm
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void FindSitesCustomEnzyme_ForwardsMethylation_AndTopology()
    {
        // A custom GATC enzyme: Dam methylates the A of GATC, which lies inside the window, so the site is blocked.
        const string gatc = "AAAGATCAAA";
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.find_sites_custom_enzyme(gatc, "GATC", 0, 4).Sites, Has.Count.EqualTo(2), "palindromic: forward and reverse record");
            Assert.That(MolToolsTools.find_sites_custom_enzyme(gatc, "GATC", 0, 4, null, MoleculeTopology.Linear, null, null, "dam").Sites, Is.Empty);
        });

        // EcoRI written by hand across the origin of a circle: GAA at 9..11, TTC at 0..2.
        const string wrapped = "TTCAAAAAAGAA";
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.find_sites_custom_enzyme(wrapped, "GAATTC", 1, 5).Sites, Is.Empty, "a linear molecule has no such site");
            var circular = MolToolsTools.find_sites_custom_enzyme(wrapped, "GAATTC", 1, 5, null, MoleculeTopology.Circular).Sites;
            Assert.That(circular, Has.Count.EqualTo(2));
            Assert.That(circular.Single(s => s.IsForwardStrand).Position, Is.EqualTo(9));
        });
    }

    [Test]
    public void RestrictionDigestCustomEnzyme_ForwardsSecondCut_AndMethylation()
    {
        // BaeI = ACNNNNGTAYC, first cut 10 before the site (-10, -15), second cut pair (23, 18): the site at 20..30 of a
        // 51-mer cuts the top strand at 10 and at 43, giving 10 / 33 / 8; with the first cut alone 10 / 41.
        string seq = new string('A', 20) + "ACAAAAGTACC" + new string('A', 20);
        var both = MolToolsTools.restriction_digest_custom_enzyme(seq, "ACNNNNGTAYC", -10, -15, "BaeI", MoleculeTopology.Linear, 23, 18).Fragments;
        var first = MolToolsTools.restriction_digest_custom_enzyme(seq, "ACNNNNGTAYC", -10, -15, "BaeI").Fragments;

        Assert.Multiple(() =>
        {
            Assert.That(both.Select(f => f.Length), Is.EqualTo(new[] { 10, 33, 8 }));
            Assert.That(first.Select(f => f.Length), Is.EqualTo(new[] { 10, 41 }));
            Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest_custom_enzyme(seq, "ACNNNNGTAYC", -10, -15, "BaeI", MoleculeTopology.Linear, 23, null));
        });

        // methylation: two GATC sites, 3 fragments; Dam blocks both, the molecule stays whole.
        const string twoSites = "AAAGATCAAAAGATCAAA";
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.restriction_digest_custom_enzyme(twoSites, "GATC", 0, 4).Fragments, Has.Count.EqualTo(3));
            Assert.That(MolToolsTools.restriction_digest_custom_enzyme(twoSites, "GATC", 0, 4, null, MoleculeTopology.Linear, null, null, "dam").Fragments, Has.Count.EqualTo(1));
        });
    }

    // ------------------------------------------------------------------------------------------------
    // A5R-T7: end2 is the second enzyme's own window, in all three forms of compatible_enzymes
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void CompatibleEnzymes_ReportsBothEnds_ForARevcompOnlyPair_InThePerEnzymeAndTheAllPairsForm()
    {
        // BsmI GAATGCN^ / BsrI ACTGGN^ anneal only through the reverse complement: the queried enzyme's window is CN, the partner's GN.
        var perEnzyme = MolToolsTools.compatible_enzymes("BsmI", "rebase").Pairs.Single(p => p.Enzyme2 == "BsrI");
        Assert.Multiple(() =>
        {
            Assert.That(perEnzyme.CompatibleEnd, Is.EqualTo("CN"));
            Assert.That(perEnzyme.End2, Is.EqualTo("GN"));
        });

        // The all-pairs listing carries the pair in table order with the FIRST enzyme's end first: Bse1I / BsmI is GN / CN.
        EnzymeCompatibilityPair? row = null;
        int total = MolToolsTools.compatible_enzymes(null, "rebase", 1, 0).TotalPairs!.Value;
        for (int offset = 0; offset < total && row is null; offset += 10000)
        {
            row = MolToolsTools.compatible_enzymes(null, "rebase", 10000, offset).Pairs
                .FirstOrDefault(p => p.Enzyme1 == "Bse1I" && p.Enzyme2 == "BsmI");
        }

        Assert.That(row, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(row!.CompatibleEnd, Is.EqualTo("GN"));
            Assert.That(row.End2, Is.EqualTo("CN"));
        });
    }
}
