using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Every optional argument of the CRISPR MCP wrappers reaches the algorithm (B08 round 5 A5-11). The auditor replaced
/// <c>system_type</c> / <c>parameters</c> / <c>max_mismatches</c> by their defaults in five wrappers at once and the whole
/// MCP test dll still passed, because every test used the SpCas9 defaults. Each test below fails for such a wrapper: the
/// Cas12a scenario is invalid for SpCas9 (a 23-nt guide, a TTTV PAM before the target), and the mismatch budget is
/// probed on both sides of the one site that decides it.
/// </summary>
[TestFixture]
public class CrisprWrapperArgumentForwardingTests
{
    // 23-nt Cas12a spacer (47.8 % GC, no TTTT).
    private const string Cas12aGuide = "ATCGGATCTAGCATGCAAGTCCA";
    // The same spacer with one mismatch at index 20 (C -> G), PAM-distal for a PAM-before system.
    private const string Cas12aOffTarget = "ATCGGATCTAGCATGCAAGTGCA";

    private static string Cas12aLocus(string spacer) => "GCGC" + "TTTA" + spacer + "GCGCGC";

    [Test]
    public void FindPamSites_ForwardsSystemType()
    {
        var sites = MolToolsTools.find_pam_sites(Cas12aLocus(Cas12aGuide), CrisprSystemType.Cas12a).Sites;

        var forward = sites.Single(s => s.IsForwardStrand && s.Position == 4);
        Assert.Multiple(() =>
        {
            Assert.That(forward.PamSequence, Is.EqualTo("TTTA"));
            Assert.That(forward.TargetSequence, Is.EqualTo(Cas12aGuide));
            Assert.That(forward.TargetStart, Is.EqualTo(8));
            Assert.That(forward.System, Is.EqualTo(CrisprDesigner.GetSystem(CrisprSystemType.Cas12a)));
        });
    }

    [Test]
    public void EvaluateGuideRna_ForwardsSystemType()
    {
        var candidate = MolToolsTools.evaluate_guide_rna(Cas12aGuide, CrisprSystemType.Cas12a);

        Assert.Multiple(() =>
        {
            Assert.That(candidate.Sequence, Is.EqualTo(Cas12aGuide));
            Assert.That(candidate.System, Is.EqualTo(CrisprDesigner.GetSystem(CrisprSystemType.Cas12a)));
            Assert.That(candidate.FullGuideRna, Is.Null, "no scaffold is appended for a Cas12a spacer (H9)");
        });
        // With the SpCas9 default the same 23-nt guide is the wrong length.
        Assert.Throws<ArgumentException>(() => MolToolsTools.evaluate_guide_rna(Cas12aGuide));
    }

    [Test]
    public void DesignGuideRnas_ForwardsSystemType_AndParameters()
    {
        // Cas12a cut of the forward site at PAM position 4: Position + pamLength + 18 = 26.
        string locus = Cas12aLocus(Cas12aGuide);

        var designed = MolToolsTools.design_guide_rnas(locus, 26, 26, CrisprSystemType.Cas12a).Guides;
        Assert.That(designed, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(designed[0].Sequence, Is.EqualTo(Cas12aGuide));
            Assert.That(designed[0].System, Is.EqualTo(CrisprDesigner.GetSystem(CrisprSystemType.Cas12a)));
            Assert.That(designed[0].Position, Is.EqualTo(8));
        });

        // parameters: a GC window the guide (47.8 %) does not satisfy removes it.
        var strict = GuideRnaParameters.Default with { MinGcContent = 60, MaxGcContent = 100 };
        Assert.That(MolToolsTools.design_guide_rnas(locus, 26, 26, CrisprSystemType.Cas12a, strict).Guides, Is.Empty);

        // ...and the same call with the SpCas9 default designs nothing at that cut.
        Assert.That(MolToolsTools.design_guide_rnas(locus, 26, 26).Guides.Where(g => g.System.Name == "Cas12a/Cpf1"), Is.Empty);
    }

    [Test]
    public void FindOffTargets_ForwardsSystemType()
    {
        var hits = MolToolsTools.find_off_targets(Cas12aGuide, Cas12aLocus(Cas12aOffTarget), 3, CrisprSystemType.Cas12a).OffTargets;

        var forward = hits.Single(h => h.IsForwardStrand);
        Assert.Multiple(() =>
        {
            Assert.That(forward.Mismatches, Is.EqualTo(1));
            Assert.That(forward.MismatchPositions, Is.EqualTo(new[] { 20 }));
            Assert.That(forward.Sequence, Is.EqualTo(Cas12aOffTarget));
        });
        Assert.Throws<ArgumentException>(() => MolToolsTools.find_off_targets(Cas12aGuide, Cas12aLocus(Cas12aOffTarget), 3),
            "the SpCas9 default rejects a 23-nt guide");
    }

    [Test]
    public void FindOffTargets_ForwardsMaxMismatches_OnBothSidesOfTheSiteThatDecidesIt()
    {
        // One SpCas9 site with exactly four mismatches (indices 0..3) and a TGG PAM.
        const string guide = "ATATATATATATATATATAT";
        const string fourMismatches = "CACAATATATATATATATAT";
        string genome = fourMismatches + "TGG";

        Assert.Multiple(() =>
        {
            Assert.That(MolToolsTools.find_off_targets(guide, genome, 3).OffTargets, Is.Empty, "budget 3 < 4 mismatches");
            Assert.That(MolToolsTools.find_off_targets(guide, genome, 4).OffTargets, Has.Count.EqualTo(1), "budget 4 reaches it");
            Assert.That(MolToolsTools.find_off_targets(guide, genome, 5).OffTargets, Has.Count.EqualTo(1));
            Assert.That(MolToolsTools.find_off_targets(guide, genome, 0).OffTargets, Is.Empty);
        });
    }

    [Test]
    public void CrisprSpecificityScore_ForwardsSystemType()
    {
        // One Cas12a off-target with a single PAM-distal mismatch: the documented heuristic charges 2 points
        // (5 per seed mismatch, 2 per distal one; the PAM-before seed is the first 12 positions).
        double score = MolToolsTools.crispr_specificity_score(Cas12aGuide, Cas12aLocus(Cas12aOffTarget), CrisprSystemType.Cas12a).Specificity;

        Assert.That(score, Is.EqualTo(98.0).Within(1e-9));
        Assert.Throws<ArgumentException>(() => MolToolsTools.crispr_specificity_score(Cas12aGuide, Cas12aLocus(Cas12aOffTarget)),
            "the SpCas9 default rejects a 23-nt guide");
    }
}
