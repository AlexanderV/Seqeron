using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Pins the two MCP wrappers that expose the PUBLISHED guide-level off-target aggregates
/// (CRISPR-OFF-001): <c>crispr_mit_specificity_score</c> (MIT/Hsu 2013, CRISPOR
/// <c>calcMitGuideScore</c>) and <c>crispr_cfd_specificity_score</c> (Doench 2016 CFD, CRISPOR
/// <c>annotateOfftargets</c>). Both must delegate to <see cref="CrisprDesigner"/> — the expected
/// values here are the same reference numbers the library-level tests lock, so a wrapper that
/// re-implemented or rescaled anything fails.
/// </summary>
[TestFixture]
public class CrisprPublishedSpecificityScoreTests
{
    private const string Guide = "GACGCATAAAGATGAGACGC";

    // One forward NGG off-target: protospacer differs at index 0 (PAM-distal, Hsu weight 0), so the
    // MIT single-hit score is 100 and the aggregate is 100/(100+100)*100 = 50.
    private const string GenomeOneOffTarget = "TACGCATAAAGATGAGACGCAGG";

    // Only the exact on-target: no off-targets.
    private const string GenomeNoOffTargets = "GACGCATAAAGATGAGACGCAGG";

    // revComp("TACGCATAAAGATGAGACGC" + "AAG"): one REVERSE-strand off-target whose own PAM is "AAG",
    // used with SpCas9-NAG so the CFD PAM score is the non-trivial "AG" entry 0.259259259.
    private const string GenomeReverseNagOffTarget = "CTTGCGTCTCATCTTTATGCGTA";

    [Test]
    public void CrisprMitSpecificityScore_Schema_ValidatesCorrectly()
    {
        Assert.Multiple(() =>
        {
            Assert.DoesNotThrow(() => MolToolsTools.crispr_mit_specificity_score(Guide, GenomeOneOffTarget));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_mit_specificity_score("", GenomeOneOffTarget));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_mit_specificity_score(null!, GenomeOneOffTarget));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_mit_specificity_score(Guide, ""));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_mit_specificity_score(Guide, null!));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_mit_specificity_score("ACGT", GenomeOneOffTarget));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_mit_specificity_score(Guide, GenomeOneOffTarget, -1));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_mit_specificity_score(Guide, GenomeOneOffTarget, 6));
        });
    }

    [Test]
    public void CrisprMitSpecificityScore_Binding_MatchesLibraryAndReference()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.crispr_mit_specificity_score(Guide, GenomeOneOffTarget).Specificity,
                Is.EqualTo(50.0).Within(1e-9));
            Assert.That(MolToolsTools.crispr_mit_specificity_score(Guide, GenomeNoOffTargets).Specificity,
                Is.EqualTo(100.0).Within(1e-12));
        });
    }

    /// <summary>
    /// A 21-nt SaCas9 guide with one real off-target must NOT report a perfect 100 (defect F9), and it
    /// must be scored with the reference's SaCas9 model (Tycko 2018, F14): single-hit score
    /// 4.2613097364785695 → aggregate 95.91285612347565.
    /// </summary>
    [Test]
    public void CrisprMitSpecificityScore_SaCas9_RealOffTarget_IsNot100()
    {
        var result = MolToolsTools.crispr_mit_specificity_score(
            "ACGTACGTACGTACGTACGTA", "ACGTACGTACGTACGTACGTCAAGAAT", 3, CrisprSystemType.SaCas9);

        Assert.Multiple(() =>
        {
            Assert.That(result.Specificity, Is.Not.EqualTo(100.0));
            Assert.That(result.Specificity, Is.EqualTo(95.91285612347565).Within(1e-9));
        });
    }

    /// <summary>
    /// Cas12a-class / CasX (PAM-before-target) systems are refused by both published aggregates:
    /// the reference reports no MIT or CFD off-target score for them.
    /// </summary>
    [TestCase(CrisprSystemType.Cas12a, 23)]
    [TestCase(CrisprSystemType.AsCas12a, 23)]
    [TestCase(CrisprSystemType.LbCas12a, 24)]
    [TestCase(CrisprSystemType.CasX, 20)]
    public void PublishedAggregates_PamBeforeTargetSystems_Throw(CrisprSystemType system, int guideLength)
    {
        string guide = new string('A', guideLength);
        const string genome = "TTTAACGTACGTACGTACGTACGTACGT";

        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() =>
                MolToolsTools.crispr_mit_specificity_score(guide, genome, 3, system));
            Assert.Throws<ArgumentException>(() =>
                MolToolsTools.crispr_cfd_specificity_score(guide, genome, 3, system));
        });
    }

    [Test]
    public void CrisprCfdSpecificityScore_Schema_ValidatesCorrectly()
    {
        Assert.Multiple(() =>
        {
            Assert.DoesNotThrow(() => MolToolsTools.crispr_cfd_specificity_score(Guide, GenomeOneOffTarget));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_cfd_specificity_score("", GenomeOneOffTarget));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_cfd_specificity_score(null!, GenomeOneOffTarget));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_cfd_specificity_score(Guide, ""));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_cfd_specificity_score(Guide, null!));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_cfd_specificity_score("ACGT", GenomeOneOffTarget));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_cfd_specificity_score(Guide, GenomeOneOffTarget, -1));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_cfd_specificity_score(Guide, GenomeOneOffTarget, 6));
            // CFD needs a 20-nt protospacer: SaCas9 (21 nt) is refused, as in the reference.
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_cfd_specificity_score(
                "ACGTACGTACGTACGTACGTA", "ACGTACGTACGTACGTACGTCAAGAAT", 3, CrisprSystemType.SaCas9));
        });
    }

    /// <summary>
    /// Both published CFD aggregations are reachable over MCP and carry the reference numbers.
    /// The single reverse-strand off-target has CFD = rG:dA@1 (1.0) x PAM "AG" (0.259259259), so
    /// MitStyle = 100/(100+0.259259259)*100 = 99.74141115651946 and ParkinsonNormalized = 100/Σ.
    /// </summary>
    [Test]
    public void CrisprCfdSpecificityScore_BothAggregations_MatchReference()
    {
        var mitStyle = MolToolsTools.crispr_cfd_specificity_score(
            Guide, GenomeReverseNagOffTarget, 3, CrisprSystemType.SpCas9_NAG);
        var parkinson = MolToolsTools.crispr_cfd_specificity_score(
            Guide, GenomeReverseNagOffTarget, 3, CrisprSystemType.SpCas9_NAG, CfdAggregation.ParkinsonNormalized);

        Assert.Multiple(() =>
        {
            Assert.That(mitStyle.Specificity, Is.EqualTo(99.74141115651946).Within(1e-9));
            Assert.That(parkinson.Specificity, Is.EqualTo(385.71428610000004).Within(1e-6));
        });
    }

    /// <summary>ParkinsonNormalized divides by the CFD sum, so a zero sum must throw, not return 100.</summary>
    [Test]
    public void CrisprCfdSpecificityScore_ParkinsonNormalized_NoOffTargets_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => MolToolsTools.crispr_cfd_specificity_score(
            Guide, GenomeNoOffTargets, 3, CrisprSystemType.SpCas9, CfdAggregation.ParkinsonNormalized));
    }

    /// <summary>
    /// The three specificity tools are genuinely different models on the same input, so the honest
    /// descriptions matter: the unsourced heuristic gives 98, the MIT aggregate 50, and the CFD
    /// aggregate its own value. A caller cannot be reading the same number under three names.
    /// </summary>
    [Test]
    public void SpecificityTools_HeuristicAndPublishedModels_AreDistinct()
    {
        double heuristic = MolToolsTools.crispr_specificity_score(Guide, GenomeOneOffTarget).Specificity;
        double mit = MolToolsTools.crispr_mit_specificity_score(Guide, GenomeOneOffTarget).Specificity;
        double cfd = MolToolsTools.crispr_cfd_specificity_score(Guide, GenomeOneOffTarget).Specificity;

        Assert.Multiple(() =>
        {
            Assert.That(heuristic, Is.EqualTo(98.0).Within(1e-12));
            Assert.That(mit, Is.EqualTo(50.0).Within(1e-9));
            Assert.That(cfd, Is.Not.EqualTo(heuristic));
            Assert.That(cfd, Is.Not.EqualTo(mit));
        });
    }
}
