using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for splice site prediction: donor/acceptor sites, gene structure.
/// Uses FsCheck to verify invariants with random DNA sequences containing canonical splice signals.
///
/// Test Units: SPLICE-DONOR-001, SPLICE-ACCEPTOR-001, SPLICE-PREDICT-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Annotation")]
public class SplicingProperties
{
    #region Generators

    /// <summary>
    /// Generates random DNA sequences from valid bases (A, C, G, T).
    /// </summary>
    private static Arbitrary<string> DnaArbitrary(int minLen = 50) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// Canonical test sequence with known GT donor and AG acceptor.
    /// Exon1(21bp) + GT...intron...AG + Exon2(21bp).
    /// </summary>
    private const string CanonicalTestSequence =
        "ATGATGAAAGCCGCCATGGCG" +       // exon1 (21bp)
        "GTAAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCAG" + // intron (~85bp with GT...AG)
        "GCGATGAAAGCCGCCATGGCG";        // exon2 (21bp)

    #endregion

    #region SPLICE-DONOR-001: R: score ∈ [0,1]; P: canonical GT at donor site; D: deterministic

    /// <summary>
    /// INV-1: Donor site scores are in [0, 1].
    /// Evidence: Scores are computed as normalized PWM log-likelihood ratios,
    /// bounded by min/max possible scores.
    /// Source: Shapiro &amp; Senapathy (1987) — consensus-based scoring.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DonorSites_Score_InRange()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var donors = SpliceSitePredictor.FindDonorSites(seq, minScore: 0.0).ToList();
            return donors.All(d => d.Score >= -1e-9 && d.Score <= 1.0 + 1e-9)
                .Label("All donor site scores must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-2 (P): EVERY canonical donor site reported sits exactly on a GT dinucleotide — the bases
    /// at [Position, Position+1] of the input are G then T. The default predictor only emits sites
    /// whose intron starts with GT (GU after T→U), so this must hold for any sequence, not just a
    /// hand-picked one. Source: Breathnach &amp; Chambon (1981); Mount (1982) — the GT-AG rule.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DonorSites_SitOnCanonicalGtDinucleotide()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var donors = SpliceSitePredictor.FindDonorSites(seq, minScore: 0.0).ToList();
            return donors.All(d => seq.Substring(d.Position, 2).ToUpperInvariant() == "GT")
                .Label("every canonical donor site must sit on a GT dinucleotide");
        });
    }

    /// <summary>
    /// INV-3: Donor site positions are within sequence bounds.
    /// Evidence: Position references a valid index in the input DNA string.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DonorSites_Positions_WithinBounds()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var donors = SpliceSitePredictor.FindDonorSites(seq, minScore: 0.0).ToList();
            return donors.All(d => d.Position >= 0 && d.Position < seq.Length)
                .Label("All donor positions must be in [0, seqLen)");
        });
    }

    /// <summary>
    /// INV-4: FindDonorSites is deterministic.
    /// Evidence: PWM scanning is a deterministic computation.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DonorSites_AreDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var d1 = SpliceSitePredictor.FindDonorSites(seq, minScore: 0.0).ToList();
            var d2 = SpliceSitePredictor.FindDonorSites(seq, minScore: 0.0).ToList();
            bool same = d1.Count == d2.Count &&
                        d1.Zip(d2).All(p => p.First.Position == p.Second.Position &&
                                             Math.Abs(p.First.Score - p.Second.Score) < 1e-10);
            return same.Label("FindDonorSites must be deterministic");
        });
    }

    #endregion

    #region SPLICE-ACCEPTOR-001: R: score ∈ [0,1]; P: canonical AG at acceptor site; D: deterministic

    /// <summary>
    /// INV-5: Acceptor site scores are in [0, 1].
    /// Evidence: Scores are normalized PWM log-likelihood ratios.
    /// Source: Shapiro &amp; Senapathy (1987).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AcceptorSites_Score_InRange()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var acceptors = SpliceSitePredictor.FindAcceptorSites(seq, minScore: 0.0).ToList();
            return acceptors.All(a => a.Score >= -1e-9 && a.Score <= 1.0 + 1e-9)
                .Label("All acceptor site scores must be in [0, 1]");
        });
    }

    /// <summary>
    /// INV-6 (P): EVERY canonical acceptor site reported ends the intron on an AG dinucleotide.
    /// The predictor reports Position as the G of the 3' AG (last intron nucleotide), so the input
    /// bases at [Position-1, Position] are A then G. Holds for any sequence. Source: Breathnach &amp;
    /// Chambon (1981) — the GT-AG rule.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AcceptorSites_EndOnCanonicalAgDinucleotide()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var acceptors = SpliceSitePredictor.FindAcceptorSites(seq, minScore: 0.0).ToList();
            return acceptors.All(a => a.Position >= 1 &&
                                      seq.Substring(a.Position - 1, 2).ToUpperInvariant() == "AG")
                .Label("every canonical acceptor site must end on an AG dinucleotide");
        });
    }

    /// <summary>
    /// INV-7: Acceptor site positions are within sequence bounds.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AcceptorSites_Positions_WithinBounds()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var acceptors = SpliceSitePredictor.FindAcceptorSites(seq, minScore: 0.0).ToList();
            return acceptors.All(a => a.Position >= 0 && a.Position < seq.Length)
                .Label("All acceptor positions must be in [0, seqLen)");
        });
    }

    /// <summary>
    /// INV-8: FindAcceptorSites is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property AcceptorSites_AreDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(50), seq =>
        {
            var a1 = SpliceSitePredictor.FindAcceptorSites(seq, minScore: 0.0).ToList();
            var a2 = SpliceSitePredictor.FindAcceptorSites(seq, minScore: 0.0).ToList();
            bool same = a1.Count == a2.Count &&
                        a1.Zip(a2).All(p => p.First.Position == p.Second.Position &&
                                             Math.Abs(p.First.Score - p.Second.Score) < 1e-10);
            return same.Label("FindAcceptorSites must be deterministic");
        });
    }

    #endregion

    #region SPLICE-PREDICT-001: R: exon start < end; P: introns flanked by GT…AG; D: deterministic

    /// <summary>
    /// INV-9: All exons have start &lt; end.
    /// Evidence: An exon spans at least one nucleotide.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_Exons_StartLessThanEnd()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(
            CanonicalTestSequence, minExonLength: 10, minIntronLength: 30);

        foreach (var exon in structure.Exons)
            Assert.That(exon.Start, Is.LessThan(exon.End),
                $"Exon start={exon.Start} must be < end={exon.End}");
    }

    /// <summary>
    /// INV-10: Gene structure has at least one exon for non-empty input.
    /// Evidence: Any DNA region is at minimum a single-exon gene.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_HasAtLeastOneExon()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(
            CanonicalTestSequence, minExonLength: 10, minIntronLength: 30);

        Assert.That(structure.Exons.Count, Is.GreaterThanOrEqualTo(1),
            "Gene structure should have at least one exon");
    }

    /// <summary>
    /// INV-11: Intron count = exon count − 1 (when introns exist).
    /// Evidence: Introns are defined as the gaps between consecutive exons.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_IntronCount_IsExonsMinusOne()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(
            CanonicalTestSequence, minExonLength: 10, minIntronLength: 30);

        if (structure.Introns.Count > 0)
            Assert.That(structure.Introns.Count, Is.EqualTo(structure.Exons.Count - 1),
                $"Expected {structure.Exons.Count - 1} introns for {structure.Exons.Count} exons");
    }

    /// <summary>
    /// INV-12: Spliced sequence is ≤ original length (introns removed).
    /// Evidence: Splicing removes intronic sequences, reducing total length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_SplicedSequence_ShorterOrEqual()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(
            CanonicalTestSequence, minExonLength: 10, minIntronLength: 30);

        Assert.That(structure.SplicedSequence.Length, Is.LessThanOrEqualTo(CanonicalTestSequence.Length),
            "Spliced sequence should be ≤ original length");
    }

    /// <summary>
    /// INV-13: Gene structure overall score is in [0, 1].
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_OverallScore_InRange()
    {
        var structure = SpliceSitePredictor.PredictGeneStructure(
            CanonicalTestSequence, minExonLength: 10, minIntronLength: 30);

        Assert.That(structure.OverallScore, Is.InRange(0.0, 1.0),
            $"Overall score {structure.OverallScore} out of range");
    }

    /// <summary>
    /// INV-14: PredictGeneStructure is deterministic.
    /// Evidence: Greedy non-overlapping intron selection from deterministic PWM scanning.
    /// </summary>
    [Test]
    [Category("Property")]
    public void GeneStructure_IsDeterministic()
    {
        var s1 = SpliceSitePredictor.PredictGeneStructure(CanonicalTestSequence, minExonLength: 10, minIntronLength: 30);
        var s2 = SpliceSitePredictor.PredictGeneStructure(CanonicalTestSequence, minExonLength: 10, minIntronLength: 30);

        Assert.Multiple(() =>
        {
            Assert.That(s1.Exons.Count, Is.EqualTo(s2.Exons.Count));
            Assert.That(s1.Introns.Count, Is.EqualTo(s2.Introns.Count));
            Assert.That(s1.SplicedSequence, Is.EqualTo(s2.SplicedSequence));
            Assert.That(Math.Abs(s1.OverallScore - s2.OverallScore), Is.LessThan(1e-10));
        });
    }

    #endregion

    #region SPLICE-MAXENT3-001: R: score finite; D: deterministic; requires 23-nt acceptor window

    // ScoreAcceptorMaxEnt — MaxEntScan score3ss (Yeo & Burge 2004): a maximum-entropy log-odds score over a
    // fixed 23-nt 3' acceptor window (20 intron + 3 exon nt). The window length is a hard contract.

    private static Arbitrary<string> FixedDnaArbitrary(int length) =>
        (from cs in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(length) select new string(cs)).ToArbitrary();

    /// <summary>INV-1 (R): the MaxEnt 3' acceptor score is finite for any valid 23-nt window.</summary>
    [FsCheck.NUnit.Property]
    public Property AcceptorMaxEnt_Score_IsFinite()
    {
        return Prop.ForAll(FixedDnaArbitrary(23), window =>
            double.IsFinite(SpliceSitePredictor.ScoreAcceptorMaxEnt(window)).Label("acceptor MaxEnt score must be finite"));
    }

    /// <summary>INV-2 (D): the MaxEnt 3' acceptor score is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property AcceptorMaxEnt_IsDeterministic()
    {
        return Prop.ForAll(FixedDnaArbitrary(23), window =>
            (SpliceSitePredictor.ScoreAcceptorMaxEnt(window) == SpliceSitePredictor.ScoreAcceptorMaxEnt(window))
                .Label("ScoreAcceptorMaxEnt must be deterministic"));
    }

    /// <summary>INV-3 (contract): a window whose length is not 23 nt is rejected.</summary>
    [FsCheck.NUnit.Property]
    public Property AcceptorMaxEnt_WrongLength_Throws()
    {
        var gen = (from len in Gen.Choose(0, 40).Where(l => l != 23)
                   from cs in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
                   select new string(cs)).ToArbitrary();
        return Prop.ForAll(gen, window =>
        {
            try { SpliceSitePredictor.ScoreAcceptorMaxEnt(window); return false.Label($"len {window.Length} accepted"); }
            catch (ArgumentException) { return true.ToProperty(); }
        });
    }

    #endregion

    #region SPLICE-MAXENT5-001: R: score finite; D: deterministic; requires 9-nt donor window

    // ScoreDonorMaxEnt — MaxEntScan score5ss (Yeo & Burge 2004): a maximum-entropy log-odds score over a
    // fixed 9-nt 5' donor window (3 exon + 6 intron nt).

    /// <summary>INV-1 (R): the MaxEnt 5' donor score is finite for any valid 9-nt window.</summary>
    [FsCheck.NUnit.Property]
    public Property DonorMaxEnt_Score_IsFinite()
    {
        return Prop.ForAll(FixedDnaArbitrary(9), window =>
            double.IsFinite(SpliceSitePredictor.ScoreDonorMaxEnt(window)).Label("donor MaxEnt score must be finite"));
    }

    /// <summary>INV-2 (D): the MaxEnt 5' donor score is deterministic.</summary>
    [FsCheck.NUnit.Property]
    public Property DonorMaxEnt_IsDeterministic()
    {
        return Prop.ForAll(FixedDnaArbitrary(9), window =>
            (SpliceSitePredictor.ScoreDonorMaxEnt(window) == SpliceSitePredictor.ScoreDonorMaxEnt(window))
                .Label("ScoreDonorMaxEnt must be deterministic"));
    }

    /// <summary>INV-3 (contract): a window whose length is not 9 nt is rejected.</summary>
    [FsCheck.NUnit.Property]
    public Property DonorMaxEnt_WrongLength_Throws()
    {
        var gen = (from len in Gen.Choose(0, 20).Where(l => l != 9)
                   from cs in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(len)
                   select new string(cs)).ToArbitrary();
        return Prop.ForAll(gen, window =>
        {
            try { SpliceSitePredictor.ScoreDonorMaxEnt(window); return false.Label($"len {window.Length} accepted"); }
            catch (ArgumentException) { return true.ToProperty(); }
        });
    }

    #endregion
}
