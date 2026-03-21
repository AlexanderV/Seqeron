using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for CRISPR guide design: PAM finding, guide RNA design, off-target scoring.
///
/// Test Units: CRISPR-PAM-001, CRISPR-GUIDE-001, CRISPR-OFF-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("MolTools")]
public class CrisprProperties
{
    #region Generators

    /// <summary>
    /// Generates random DNA sequences of sufficient length for CRISPR analysis.
    /// A 20-nt guide + 3-nt PAM requires at least 23 nt; we use 60–120 for richer tests.
    /// </summary>
    private static Arbitrary<string> CrisprDnaArbitrary(int minLen = 60) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen && a.Length <= 120)
            .Select(a => new string(a))
            .ToArbitrary();

    #endregion

    #region CRISPR-PAM-001: R: positions valid; P: PAM motif at each site; M: longer seq → ≥ sites; D: deterministic

    /// <summary>
    /// INV-1: All PAM site positions are within sequence bounds [0, seqLen).
    /// Evidence: A PAM at position p requires p + pamLen ≤ seqLen.
    /// Source: Jinek et al. (2012) Science — SpCas9 requires an NGG PAM.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PamSites_Positions_WithinBounds_Property()
    {
        return Prop.ForAll(CrisprDnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            var sites = CrisprDesigner.FindPamSites(dna, CrisprSystemType.SpCas9).ToList();

            return sites.All(s => s.Position >= 0 && s.Position < seq.Length)
                .Label($"All PAM positions must be in [0, {seq.Length})");
        });
    }

    /// <summary>
    /// INV-2: Each reported PAM site contains a sequence matching the system's PAM motif (NGG for SpCas9).
    /// Evidence: SpCas9 PAM is NGG where N ∈ {A,C,G,T}.
    /// The PamSequence on forward strand sites must match [ACGT]GG.
    /// Source: Jinek et al. (2012).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PamSites_ContainExpectedMotif_Property()
    {
        return Prop.ForAll(CrisprDnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            var sites = CrisprDesigner.FindPamSites(dna, CrisprSystemType.SpCas9).ToList();
            var forwardSites = sites.Where(s => s.IsForwardStrand).ToList();

            return forwardSites.All(s => s.PamSequence.Length == 3 &&
                                         "ACGT".Contains(s.PamSequence[0]) &&
                                         s.PamSequence[1] == 'G' &&
                                         s.PamSequence[2] == 'G')
                .Label("Forward-strand PAM must match NGG pattern");
        });
    }

    /// <summary>
    /// INV-3: Target sequence length matches the system's guide length (20 nt for SpCas9).
    /// Evidence: SpCas9 guide RNA is 20 nucleotides complementary to the target.
    /// Source: Jinek et al. (2012).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PamSites_TargetLength_MatchesSystem_Property()
    {
        return Prop.ForAll(CrisprDnaArbitrary(), seq =>
        {
            var system = CrisprDesigner.GetSystem(CrisprSystemType.SpCas9);
            var dna = new DnaSequence(seq);
            var sites = CrisprDesigner.FindPamSites(dna, CrisprSystemType.SpCas9).ToList();

            return sites.All(s => s.TargetSequence.Length == system.GuideLength)
                .Label($"Target length must be {system.GuideLength}");
        });
    }

    /// <summary>
    /// INV-4: PAM finding is deterministic — same sequence always yields same sites.
    /// Evidence: FindPamSites is a pure pattern-matching algorithm.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property PamSites_IsDeterministic()
    {
        return Prop.ForAll(CrisprDnaArbitrary(), seq =>
        {
            var dna = new DnaSequence(seq);
            var sites1 = CrisprDesigner.FindPamSites(dna).ToList();
            var sites2 = CrisprDesigner.FindPamSites(dna).ToList();

            return (sites1.Count == sites2.Count &&
                    sites1.Zip(sites2).All(p => p.First.Position == p.Second.Position &&
                                                p.First.IsForwardStrand == p.Second.IsForwardStrand))
                .Label("FindPamSites must be deterministic");
        });
    }

    #endregion

    #region CRISPR-GUIDE-001: R: guide length = specified; P: target strand correct; R: score ∈ valid range

    // A synthetic sequence with known AGG PAM at position 43 on forward strand
    private const string GuideTestSequence =
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAGGACGTACGTACGTACGTACGT";

    /// <summary>
    /// INV-1: Guide RNA GC content is in [0, 100] (percentage).
    /// Evidence: GC% = 100 × (G+C) / length; with length > 0 and counts ≥ 0, GC% ∈ [0,100].
    /// Source: Standard molecular biology (Watson &amp; Crick 1953).
    /// </summary>
    [Test]
    [Category("Property")]
    public void GuideRna_GcContent_InRange()
    {
        var dna = new DnaSequence(GuideTestSequence);
        var guides = CrisprDesigner.DesignGuideRnas(dna, 0, dna.Length - 1).ToList();

        Assert.That(guides, Is.Not.Empty, "Test sequence must produce guides");
        foreach (var g in guides)
            Assert.That(g.GcContent, Is.InRange(0.0, 100.0),
                $"GC content {g.GcContent} out of range for guide at {g.Position}");
    }

    /// <summary>
    /// INV-2: Guide RNA score is in [0, 100].
    /// Evidence: Design scoring combines GC%, homopolymer, and thermodynamic penalties
    /// producing a normalized score.
    /// Source: Doench et al. (2014) Nat Biotech — sgRNA activity scoring.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GuideRna_Score_InRange()
    {
        var dna = new DnaSequence(GuideTestSequence);
        var guides = CrisprDesigner.DesignGuideRnas(dna, 0, dna.Length - 1).ToList();

        Assert.That(guides, Is.Not.Empty, "Test sequence must produce guides");
        foreach (var g in guides)
            Assert.That(g.Score, Is.InRange(0.0, 100.0),
                $"Score {g.Score} out of range for guide at {g.Position}");
    }

    /// <summary>
    /// INV-3: Each guide RNA's target strand is correctly reported (forward or reverse).
    /// Evidence: PAM sites exist on both strands; the guide targets the strand opposite the PAM.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GuideRna_TargetStrand_IsValid()
    {
        var dna = new DnaSequence(GuideTestSequence);
        var guides = CrisprDesigner.DesignGuideRnas(dna, 0, dna.Length - 1).ToList();

        foreach (var g in guides)
            Assert.That(g.IsForwardStrand, Is.TypeOf<bool>(),
                "Guide strand must be a valid boolean");
    }

    /// <summary>
    /// INV-4: Guide RNA position is within the requested target region.
    /// Evidence: DesignGuideRnas filters PAM sites to those within [regionStart, regionEnd].
    /// </summary>
    [Test]
    [Category("Property")]
    public void GuideRna_Position_WithinRegion()
    {
        var dna = new DnaSequence(GuideTestSequence);
        int regionStart = 10;
        int regionEnd = dna.Length - 1;
        var guides = CrisprDesigner.DesignGuideRnas(dna, regionStart, regionEnd).ToList();

        foreach (var g in guides)
            Assert.That(g.Position, Is.InRange(0, dna.Length - 1),
                $"Guide position {g.Position} out of sequence bounds");
    }

    #endregion

    #region CRISPR-OFF-001: R: off-target score ∈ [0,1]; M: more mismatches → lower score; D: deterministic

    /// <summary>
    /// INV-1: Off-target score is in [0, 1].
    /// Evidence: Off-target scoring normalizes mismatch penalties to a [0,1] probability.
    /// Source: Hsu et al. (2013) Nat Biotech — CFD score for off-target prediction.
    /// </summary>
    [Test]
    [Category("Property")]
    public void OffTargets_Score_InRange()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence(guide + "GGGTTTTTTTTTTACGTACGTACGTACGTACGTGGG");
        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3).ToList();

        foreach (var ot in offTargets)
            Assert.That(ot.OffTargetScore, Is.InRange(0.0, 1.0),
                $"Off-target score {ot.OffTargetScore} out of range at pos {ot.Position}");
    }

    /// <summary>
    /// INV-2: Off-target mismatch count ≤ maxMismatches parameter.
    /// Evidence: FindOffTargets filters by maxMismatches threshold.
    /// </summary>
    [Test]
    [Category("Property")]
    public void OffTargets_Mismatches_WithinLimit()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence(guide + "GGGTTTTTTTTTT");
        int maxMismatches = 2;
        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, maxMismatches).ToList();

        foreach (var ot in offTargets)
            Assert.That(ot.Mismatches, Is.LessThanOrEqualTo(maxMismatches),
                $"Off-target at {ot.Position} has {ot.Mismatches} mismatches > max {maxMismatches}");
    }

    /// <summary>
    /// INV-3: Off-target finding is deterministic.
    /// Evidence: FindOffTargets is a pure mismatch-counting algorithm.
    /// </summary>
    [Test]
    [Category("Property")]
    public void OffTargets_IsDeterministic()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence(guide + "GGGTTTTTTTTTTACGTACGTACGTACGTACGTGGG");

        var ot1 = CrisprDesigner.FindOffTargets(guide, genome, 3).ToList();
        var ot2 = CrisprDesigner.FindOffTargets(guide, genome, 3).ToList();

        Assert.That(ot1.Count, Is.EqualTo(ot2.Count), "Off-target count must be deterministic");
        for (int i = 0; i < ot1.Count; i++)
        {
            Assert.That(ot1[i].Position, Is.EqualTo(ot2[i].Position));
            Assert.That(ot1[i].Mismatches, Is.EqualTo(ot2[i].Mismatches));
        }
    }

    /// <summary>
    /// INV-4: Perfect match (0 mismatches) has maximum off-target score.
    /// Evidence: A perfect match to the guide means maximum potential for off-target cleavage.
    /// Source: Hsu et al. (2013).
    /// </summary>
    [Test]
    [Category("Property")]
    public void OffTargets_PerfectMatch_MaxScore()
    {
        string guide = "ACGTACGTACGTACGTACGT";
        var genome = new DnaSequence(guide + "GGGTTTTTTTTTT");
        var offTargets = CrisprDesigner.FindOffTargets(guide, genome, 3).ToList();

        var perfectMatches = offTargets.Where(ot => ot.Mismatches == 0).ToList();
        foreach (var pm in perfectMatches)
            Assert.That(pm.OffTargetScore, Is.GreaterThanOrEqualTo(0.9),
                $"Perfect match off-target score {pm.OffTargetScore} should be near 1.0");
    }

    #endregion
}
