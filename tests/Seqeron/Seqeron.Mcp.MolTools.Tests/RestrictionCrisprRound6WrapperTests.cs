using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// The MCP-visible side of audit round 6 (B08 A6-1 .. A6-3, U3): what a caller of the tools sees for the library defects
/// the confirmation audit found. Each test failed on the code the auditor read (bec761ac).
/// </summary>
[TestFixture]
public class RestrictionCrisprRound6WrapperTests
{
    private const string Guide = "GACGCATAAAGATGAGACGC";

    private static readonly string RepeatGenome = string.Join("TTTTT", Enumerable.Repeat(Guide + "TGG", 3));

    // A6-1: the alternative-PAM cutoff is refused whatever the genome holds (a repeat guide used to return 0 instead).
    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    public void CrisprMitSpecificityScore_InvalidAlternativePamCutoff_IsRefused_ForARepeatGuideToo(double cutoff)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.crispr_mit_specificity_score(
            Guide, RepeatGenome, 4, CrisprSystemType.SpCas9, include_alternative_pams: true, min_alternative_pam_hit_score: cutoff, max_occurrences: 2));
    }

    [TestCase(-1.0)]
    [TestCase(double.NaN)]
    public void CrisprCfdSpecificityScore_InvalidAlternativePamCutoff_IsRefused_ForARepeatGuideToo(double cutoff)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.crispr_cfd_specificity_score(
            Guide, RepeatGenome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, include_alternative_pams: true, min_alternative_pam_hit_score: cutoff, max_occurrences: 2));
    }

    [Test]
    public void CrisprSpecificityScores_ARepeatGuide_StillScoresZero_WithAValidCutoff()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.crispr_mit_specificity_score(Guide, RepeatGenome, 4, include_alternative_pams: true, max_occurrences: 2).Specificity, Is.EqualTo(0.0));
            Assert.That(MolToolsTools.crispr_cfd_specificity_score(Guide, RepeatGenome, 4, include_alternative_pams: true, max_occurrences: 2).Specificity, Is.EqualTo(0.0));
        });
    }

    // A6-2: a U is not one of the 15 IUPAC codes, so it is never part of a site on either strand.
    [TestCase("AAAGAAUUCAAA")]
    [TestCase("AAAGAATUCAAA")]
    public void FindRestrictionSites_AUracilInTheSite_IsNoSite(string sequence)
    {
        foreach (var topology in new[] { MoleculeTopology.Linear, MoleculeTopology.Circular })
        {
            Assert.Multiple(() =>
            {
                Assert.That(MolToolsTools.find_restriction_sites(sequence, ["EcoRI"], topology).Sites, Is.Empty, $"{topology} sites");
                Assert.That(MolToolsTools.restriction_digest(sequence, ["EcoRI"], topology).Fragments, Has.Count.EqualTo(1), $"{topology} digest");
            });
        }
    }

    [Test]
    public void FindRestrictionSites_TheSameSiteWithT_IsFound()
    {
        Assert.That(MolToolsTools.find_restriction_sites("AAAGAATTCAAA", ["EcoRI"]).Sites, Is.Not.Empty);
    }

    // A6-3: the stagger of a custom enzyme is bounded; the call is refused at once instead of exhausting memory.
    [TestCase(0, int.MinValue)]
    [TestCase(int.MaxValue, int.MinValue)]
    [TestCase(int.MaxValue, 5)]
    [TestCase(50_000_000, 5)]
    [TestCase(0, 65)]
    public void RestrictionDigestCustomEnzyme_StaggerBeyondTheLimit_IsRefused(int forward, int reverse)
    {
        var e = Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest_custom_enzyme(
            "AAAGAATTCAAAGAATTCAAAGGGGGGGGGGGG", "GAATTC", forward, reverse, topology: MoleculeTopology.Circular));

        Assert.That(e!.Message, Does.Contain("64"));
    }

    [Test]
    public void RestrictionDigestCustomEnzyme_StaggerAtTheLimit_Digests()
    {
        var result = MolToolsTools.restriction_digest_custom_enzyme(
            "AAAGAATTCAAAGAATTCAAAGGGGGGGGGGGG", "GAATTC", 0, 64, topology: MoleculeTopology.Circular);

        Assert.That(result.Fragments, Is.Not.Empty);
    }

    // U3: the per-enzyme form of compatible_enzymes reports enzyme1 under the table's spelling, as enzymes_compatible does.
    [TestCase(" bamhi ")]
    [TestCase("BAMHI")]
    [TestCase("bamhi")]
    public void CompatibleEnzymes_PerEnzymeForm_EchoesTheCanonicalName(string spelling)
    {
        var pairs = MolToolsTools.compatible_enzymes(spelling).Pairs;

        Assert.That(pairs, Is.Not.Empty);
        Assert.That(pairs.Select(p => p.Enzyme1).Distinct(), Is.EqualTo(new[] { "BamHI" }));
    }

    [Test]
    public void CompatibleEnzymes_PerEnzymeForm_ResolvesAHyphenSpelling_ToTheTableName()
    {
        string canonical = MolToolsTools.get_enzyme("CviKI-1").Enzyme!.Name;

        var pairs = MolToolsTools.compatible_enzymes("  cviki-1 ", "rebase").Pairs;

        Assert.That(pairs.Select(p => p.Enzyme1).Distinct(), Is.EqualTo(new[] { canonical }));
    }
}
