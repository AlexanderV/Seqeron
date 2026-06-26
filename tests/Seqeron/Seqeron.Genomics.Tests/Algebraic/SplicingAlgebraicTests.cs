namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the Splicing area — the MaxEntScan maximum-entropy
/// splice-site scores (Yeo &amp; Burge 2004).
///
/// Algebraic testing pins the reference-window identity values (the published
/// MaxEntScan score3/score5 anchors), the T≡U spelling invariance of the model,
/// and the determinism of rescoring.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 250, 251.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Splicing")]
public class SplicingAlgebraicTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SPLICE-MAXENT3-001 — MaxEntScan 3' acceptor score (Splicing), row 250.
    // ID — the canonical 23-nt reference window reproduces the published score3 = 2.89 bits.
    // IDEMP — the score is a pure, deterministic function.
    //   — SpliceSitePredictor.ScoreAcceptorMaxEnt; Yeo & Burge (2004) J Comput Biol 11:377.
    // ═══════════════════════════════════════════════════════════════════════

    // 23 nt: 20 intron + 3 exon, the documented maxentpy score3 example.
    private const string AcceptorReferenceWindow = "ttccaaacgaacttttgtAGgga";

    [Test]
    public void AcceptorMaxEnt_Identity_ReferenceWindowReproduces289()
    {
        double score = SpliceSitePredictor.ScoreAcceptorMaxEnt(AcceptorReferenceWindow);
        System.Math.Round(score, 2).Should().Be(2.89,
            "MaxEntScan score3 of the canonical documented window must be 2.89 bits (Yeo & Burge 2004)");
    }

    [Test]
    public void AcceptorMaxEnt_Invariant_DnaAndRnaSpellingAgree()
    {
        // The maximum-entropy model reads T as U: the DNA and RNA spellings of the same window
        // must score identically.
        double dna = SpliceSitePredictor.ScoreAcceptorMaxEnt(AcceptorReferenceWindow);
        double rna = SpliceSitePredictor.ScoreAcceptorMaxEnt(AcceptorReferenceWindow.Replace('t', 'u').Replace('T', 'U'));
        rna.Should().BeApproximately(dna, 1e-12);
    }

    [Test]
    public void AcceptorMaxEnt_Idempotent_Deterministic()
    {
        SpliceSitePredictor.ScoreAcceptorMaxEnt(AcceptorReferenceWindow)
            .Should().Be(SpliceSitePredictor.ScoreAcceptorMaxEnt(AcceptorReferenceWindow));
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SPLICE-MAXENT5-001 — MaxEntScan 5' donor score (Splicing), row 251.
    // ID — the canonical 9-nt reference window cagGTAAGT reproduces score5 = 10.86 bits.
    // IDEMP — the score is a pure, deterministic function.
    //   — SpliceSitePredictor.ScoreDonorMaxEnt; Yeo & Burge (2004).
    // ═══════════════════════════════════════════════════════════════════════

    // 9 nt: 3 exon + 6 intron, invariant GT at 0-based 3-4 — the documented maxentpy score5 example.
    private const string DonorReferenceWindow = "cagGTAAGT";

    [Test]
    public void DonorMaxEnt_Identity_ReferenceWindowReproduces1086()
    {
        double score = SpliceSitePredictor.ScoreDonorMaxEnt(DonorReferenceWindow);
        System.Math.Round(score, 2).Should().Be(10.86,
            "MaxEntScan score5 of the canonical documented window must be 10.86 bits (Yeo & Burge 2004)");
    }

    [Test]
    public void DonorMaxEnt_Invariant_DnaAndRnaSpellingAgree()
    {
        double dna = SpliceSitePredictor.ScoreDonorMaxEnt(DonorReferenceWindow);
        double rna = SpliceSitePredictor.ScoreDonorMaxEnt(DonorReferenceWindow.Replace('T', 'U'));
        rna.Should().BeApproximately(dna, 1e-12);
    }

    [Test]
    public void DonorMaxEnt_Idempotent_Deterministic()
    {
        SpliceSitePredictor.ScoreDonorMaxEnt(DonorReferenceWindow)
            .Should().Be(SpliceSitePredictor.ScoreDonorMaxEnt(DonorReferenceWindow));
    }
}
