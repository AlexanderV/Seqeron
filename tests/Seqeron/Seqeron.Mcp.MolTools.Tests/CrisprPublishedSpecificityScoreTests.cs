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
    #region F50 — include_alternative_pams (CRISPOR's x0.2 alternative-PAM MIT downweight)

    // Genome with two off-targets of Guide: one on the main AGG PAM (calcHitScore 60.5) and one on an
    // alternative TAG PAM (31.4999...). Reference values from a verbatim Python copy of crispor.py's
    // calcHitScore / calcMitGuideScore plus annotateOfftargets' "if pam==\"NGG\" and otSeq[-2:]!=\"GG\":
    // mitScore = mitScore * 0.2" (L1947-1958):
    //   main PAM only              -> 100/(100+60.5)*100            = 62.30529595015576
    //   with the alt-PAM hit x0.2  -> 100/(100+60.5+6.2999...)*100   = 59.95203836930455
    //   (at FULL weight it would be 52.083333333333336, so the 0.2 is visible in the value)
    private const string GenomeMainAndAltPam =
        "GACGCTTAAAGATGAGACGC" + "AGG" + "TTTT" + "GACGCATAAAGATGAGACTC" + "TAG";

    /// <summary>
    /// F50: the new optional parameter is OFF by default — the tool's existing behaviour is unchanged —
    /// and ON it folds in the alternative-PAM (NAG) hit with the reference's 0.2 weight.
    /// </summary>
    [Test]
    public void CrisprMitSpecificityScore_IncludeAlternativePams_DefaultsOffAndDownweightsWhenOn()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.crispr_mit_specificity_score(Guide, GenomeMainAndAltPam, 4).Specificity,
                Is.EqualTo(62.30529595015576).Within(1e-9), "default must stay main-PAM only");
            Assert.That(
                MolToolsTools.crispr_mit_specificity_score(
                    Guide, GenomeMainAndAltPam, 4, CrisprSystemType.SpCas9, include_alternative_pams: false).Specificity,
                Is.EqualTo(62.30529595015576).Within(1e-9));
            Assert.That(
                MolToolsTools.crispr_mit_specificity_score(
                    Guide, GenomeMainAndAltPam, 4, CrisprSystemType.SpCas9, include_alternative_pams: true).Specificity,
                Is.EqualTo(59.95203836930455).Within(1e-9));
        });
    }

    /// <summary>
    /// F50: the flag changes nothing where no alternative-PAM hit exists, it keeps refusing
    /// PAM-before-target systems, and it keeps validating max_mismatches.
    /// </summary>
    [Test]
    public void CrisprMitSpecificityScore_IncludeAlternativePams_InertWithoutAltHits_AndStillValidates()
    {
        Assert.Multiple(() =>
        {
            Assert.That(
                MolToolsTools.crispr_mit_specificity_score(
                    Guide, GenomeOneOffTarget, 4, CrisprSystemType.SpCas9, include_alternative_pams: true).Specificity,
                Is.EqualTo(50.0).Within(1e-9));
            Assert.That(
                MolToolsTools.crispr_mit_specificity_score(
                    Guide, GenomeNoOffTargets, 4, CrisprSystemType.SpCas9, include_alternative_pams: true).Specificity,
                Is.EqualTo(100.0).Within(1e-12));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_mit_specificity_score(
                "ACGTACGTACGTACGTACGTACG", "TTTAACGTACGTACGTACGTACGTACC", 3, CrisprSystemType.Cas12a,
                include_alternative_pams: true));
            Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_mit_specificity_score(
                Guide, GenomeMainAndAltPam, 6, CrisprSystemType.SpCas9, include_alternative_pams: true));
        });
    }

    /// <summary>
    /// B08 round 2 E1/F66-F67: the tool now scans the NGA alternative PAM (was silently dropped, 100.0) and
    /// forwards <c>min_alternative_pam_hit_score</c>. Guide GACGCATAAAGATGAGACGC, genome = 1-mismatch hit
    /// (calcHitScore 60.5, crispor.py code re-run) + CGA + filler: 100/(100+0.2x60.5)x100 = 89.20606601248885;
    /// a cutoff of 61 (> 60.5) drops that hit -> 100.
    /// </summary>
    [Test]
    public void CrisprMitSpecificityScore_IncludeAlternativePams_ScansNga_AndForwardsTheCutoff()
    {
        const string genome = "GACGCTTAAAGATGAGACGC" + "CGA" + "TTTTTT";
        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.crispr_mit_specificity_score(Guide, genome, 4).Specificity,
                Is.EqualTo(100.0));
            Assert.That(
                MolToolsTools.crispr_mit_specificity_score(
                    Guide, genome, 4, CrisprSystemType.SpCas9, include_alternative_pams: true).Specificity,
                Is.EqualTo(89.20606601248885).Within(1e-9));
            Assert.That(
                MolToolsTools.crispr_mit_specificity_score(
                    Guide, genome, 4, CrisprSystemType.SpCas9, include_alternative_pams: true,
                    min_alternative_pam_hit_score: 61.0).Specificity,
                Is.EqualTo(100.0));
        });
    }

    /// <summary>
    /// F50: only the MIT aggregate is downweighted — the CFD tool is unaffected by the rule and still
    /// sees only its own system's PAM, so its value on this genome does not move.
    /// </summary>
    [Test]
    public void CrisprCfdSpecificityScore_NotAffectedByTheAlternativePamRule()
    {
        double cfd = MolToolsTools.crispr_cfd_specificity_score(Guide, GenomeMainAndAltPam, 4).Specificity;
        double mitMain = MolToolsTools.crispr_mit_specificity_score(Guide, GenomeMainAndAltPam, 4).Specificity;
        double mitAlt = MolToolsTools.crispr_mit_specificity_score(
            Guide, GenomeMainAndAltPam, 4, CrisprSystemType.SpCas9, include_alternative_pams: true).Specificity;

        Assert.Multiple(() =>
        {
            Assert.That(mitAlt, Is.LessThan(mitMain), "the alt-PAM hit lowers the MIT aggregate");
            Assert.That(cfd, Is.Not.EqualTo(mitMain));
            Assert.That(cfd, Is.Not.EqualTo(mitAlt));
        });
    }

    #endregion

    #region B08 round 3 G1 — crispr_cfd_specificity_score include_alternative_pams

    /// <summary>
    /// G1/F80: default OFF (a NAG-only guide scores 100, today's behaviour); ON it reproduces CRISPOR's wider
    /// search: calcCfdScore of the TAG hit = 0.18518518492592592 -> 100/(100+0.185...)*100 = 99.81515711670932
    /// (crispor.py calcCfdScore exec'd in Python, see CrisprDesigner_OffTarget_Tests region G1). The cutoff
    /// parameter is forwarded: the mixed genome (NGG 60.5, NAG 4-mm 0.3129, NGA 31.5) gives 99.24515917350747
    /// by default (ALTPAMMINSCORE 1.0), 99.20727233740708 at 0.
    /// </summary>
    [Test]
    public void CrisprCfdSpecificityScore_IncludeAlternativePams_DefaultsOff_AndScansNagNga_AndForwardsTheCutoff()
    {
        const string nagOnly = "GACGCTTAAAGATGAGACGC" + "TAG" + "TTTTTT";
        const string mixed = "GACGCTTAAAGATGAGACGC" + "AGG" + "TTTT" + "GTCGCTTCAAGATGAGACTC" + "TAG" + "TTTT"
            + "GACGCATAAAGATGAGACTC" + "CGA";

        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.crispr_cfd_specificity_score(Guide, nagOnly, 4).Specificity, Is.EqualTo(100.0));
            Assert.That(
                MolToolsTools.crispr_cfd_specificity_score(
                    Guide, nagOnly, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, include_alternative_pams: true).Specificity,
                Is.EqualTo(99.81515711670932).Within(1e-9));
            Assert.That(
                MolToolsTools.crispr_cfd_specificity_score(
                    Guide, mixed, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, include_alternative_pams: true).Specificity,
                Is.EqualTo(99.24515917350747).Within(1e-9));
            Assert.That(
                MolToolsTools.crispr_cfd_specificity_score(
                    Guide, mixed, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, include_alternative_pams: true,
                    min_alternative_pam_hit_score: 0.0).Specificity,
                Is.EqualTo(99.20727233740708).Within(1e-9));
        });
    }

    #endregion
}
