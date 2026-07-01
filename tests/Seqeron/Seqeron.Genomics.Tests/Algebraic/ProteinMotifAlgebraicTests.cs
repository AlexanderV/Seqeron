using System;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the ProteinMotif area — the Plan7 profile-HMM domain
/// search (HMMER3-style).
///
/// Algebraic testing pins the zero-reference of the log-odds bit score (a null /
/// non-homologous sequence carries no positive signal), the exact HMMER pipeline
/// score identity (search bits = pre-null2 bits − null2 bias), the sign of the
/// biased-composition correction, and the determinism of rescoring.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, row 239.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("ProteinMotif")]
public class ProteinMotifAlgebraicTests
{
    // SRC_HUMAN SH3 core (UniProt P12931, residues 84–137) — a genuine SH3 domain.
    private const string Sh3TruePositive =
        "TFVALYDYESRTETDLSFKKGERLQIVNNTEGDWWLAHSLSTGQTGYIPSNYVAP";

    // Low-complexity, non-domain sequence: drawn from a biased background, homologous
    // to no profile — the null-model negative control.
    private const string NullModelSequence =
        "AAAAAAAAAAAAAAEEEEEEEEEEEEEEKKKKKKKKKKKK";

    private const string Sh3Accession = "PF00018";

    private static Plan7ProfileHmm Sh3Profile() =>
        Plan7ProfileHmm.Parse(ReadEmbedded("PF00018_SH3_1.hmm"));

    private static string ReadEmbedded(string fileName)
    {
        var asm = typeof(Plan7ProfileHmm).Assembly;
        string resourceName = $"Seqeron.Genomics.Analysis.Resources.{fileName}";
        using var stream = asm.GetManifestResourceStream(resourceName)
            ?? throw new InvalidOperationException($"Missing embedded resource: {resourceName}");
        using var reader = new System.IO.StreamReader(stream);
        return reader.ReadToEnd();
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PROTMOTIF-HMM-001 — Plan7 profile-HMM domain search (ProteinMotif), row 239.
    //
    // Model: the reported score is a log-odds bit score, bits = log2 P(seq|profile)/P(seq|null).
    //        For a sequence that follows the null model the odds ratio is ≈ 1, so its bit
    //        score sits at the zero reference (no positive evidence), whereas a genuine
    //        homolog scores far above. The HMMER pipeline reports
    //        search-bits = pre-null2-bits − null2-bias, with null2-bias = logsumexp(0,…) ≥ 0.
    //   — Eddy (2011) PLoS Comput Biol 7:e1002195; Durbin et al. (1998) §5.4;
    //     Plan7ProfileHmm.{HmmSearchBitScore,LocalForwardBitScore,Null2BiasBits};
    //     ProteinMotifFinder.ScoreDomainHmm. TestSpec tests/TestSpecs/PROTMOTIF-HMM-001.md.
    //
    // Laws (row 239): ID — a null-model sequence → bit score ≈ 0 (no positive signal).
    //                 IDEMP — rescoring the same sequence is stable (deterministic).
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Hmm_Identity_NullModelSequenceScoresAtZeroReference()
    {
        // ID: the log-odds of a non-homologous (null-composition) sequence carries no positive
        // evidence — its bit score sits at or below the zero reference, while a genuine SH3
        // domain scores far above it. (bit = log2 P(seq|model)/P(seq|null); null ⇒ odds ≈ 1.)
        double nullBits = ProteinMotifFinder.ScoreDomainHmm(NullModelSequence, Sh3Accession);
        double homologBits = ProteinMotifFinder.ScoreDomainHmm(Sh3TruePositive, Sh3Accession);

        nullBits.Should().BeLessThanOrEqualTo(0.0,
            "a null-model sequence must not produce a positive log-odds bit score");
        homologBits.Should().BeGreaterThan(10.0,
            "a genuine SH3 domain must score far above the null zero reference");
    }

    [Test]
    public void Hmm_Distributive_PipelineScoreIdentity()
    {
        // The exact HMMER pipeline identity (p7_pipeline.c):
        //   HmmSearchBitScore = LocalForwardBitScore − Null2BiasBits.
        var hmm = Sh3Profile();
        foreach (string seq in new[] { Sh3TruePositive, NullModelSequence })
        {
            double search = hmm.HmmSearchBitScore(seq);
            double pre = hmm.LocalForwardBitScore(seq);
            double bias = hmm.Null2BiasBits(seq);
            search.Should().BeApproximately(pre - bias, 1e-6,
                $"pipeline identity must hold for \"{seq[..Math.Min(8, seq.Length)]}…\"");
        }
    }

    [Test]
    public void Hmm_Null2BiasNeverRaisesScore()
    {
        // null2 bias = logsumexp(0, …) ≥ 0 ⇒ the biased-composition correction can only
        // lower the reported score: HmmSearchBitScore ≤ LocalForwardBitScore.
        var hmm = Sh3Profile();
        foreach (string seq in new[] { Sh3TruePositive, NullModelSequence })
        {
            hmm.Null2BiasBits(seq).Should().BeGreaterThanOrEqualTo(0.0, "null2 bias is a log-sum-exp with a 0 term");
            hmm.HmmSearchBitScore(seq).Should().BeLessThanOrEqualTo(hmm.LocalForwardBitScore(seq) + 1e-9,
                "the bias correction never raises the score");
        }
    }

    [Test]
    public void Hmm_Idempotent_RescoringIsStable()
    {
        // IDEMP: the per-region ensemble is reseeded to a fixed pipeline seed, so rescoring the
        // same (profile, sequence) is bit-for-bit stable.
        var hmm = Sh3Profile();
        hmm.HmmSearchBitScore(Sh3TruePositive).Should().Be(hmm.HmmSearchBitScore(Sh3TruePositive));
        hmm.LocalForwardBitScore(Sh3TruePositive).Should().Be(hmm.LocalForwardBitScore(Sh3TruePositive));
        ProteinMotifFinder.ScoreDomainHmm(Sh3TruePositive, Sh3Accession)
            .Should().Be(ProteinMotifFinder.ScoreDomainHmm(Sh3TruePositive, Sh3Accession));
    }
}
