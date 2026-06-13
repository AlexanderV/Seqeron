// PROTMOTIF-SP-001 — Signal Peptide Cleavage-Site Prediction (von Heijne 1986 weight matrix)
// Evidence: docs/Evidence/PROTMOTIF-SP-001-Evidence.md
// TestSpec: tests/TestSpecs/PROTMOTIF-SP-001.md
// Source:   von Heijne G (1986) Nucleic Acids Res 14(11):4683-4690; EMBOSS 6.6.0 sigcleave (data/Esig.euk, sigcleave.c).

using System;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for PROTMOTIF-SP-001: von Heijne (1986) signal-peptide cleavage-site prediction.
/// Expected values are derived from the EMBOSS <c>sigcleave</c> reference output and an independent
/// re-derivation of the log-odds score (natural log of observed/expected residue frequency).
/// </summary>
[TestFixture]
[Category("PROTMOTIF-SP-001")]
public class ProteinMotifFinder_PredictSignalPeptide_Tests
{
    #region Evidence-sourced constants

    // UniProt P17644 (ACH2_DROME), 576 aa — the EMBOSS sigcleave worked example.
    // Source: https://rest.uniprot.org/uniprotkb/P17644.fasta (accessed 2026-06-14)
    private const string Ach2Drome =
        "MAPGCCTTRPRPIALLAHIWRHCKPLCLLLVLLLLCETVQANPDAKRLYDDLLSNYNRLI" +
        "RPVSNNTDTVLVKLGLRLSQLIDLNLKDQILTTNVWLEHEWQDHKFKWDPSEYGGVTELY" +
        "VPSEHIWLPDIVLYNNADGEYVVTTMTKAILHYTGKVVWTPPAIFKSSCEIDVRYFPFDQ" +
        "QTCFMKFGSWTYDGDQIDLKHISQKNDKDNKVEIGIDLREYYPSVEWDILGVPAERHEKY" +
        "YPCCAEPYPDIFFNITLRRKTLFYTVNLIIPCVGISYLSVLVFYLPADSGEKIALCISIL" +
        "LSQTMFFLLISEIIPSTSLALPLLGKYLLFTMLLVGLSVVITIIILNIHYRKPSTHKMRP" +
        "WIRSFFIKRLPKLLLMRVPKDLLRDLAANKINYGLKFSKTKFGQALMDEMQMNSGGSSPD" +
        "SLRRMQGRVGAGGCNGMHVTTATNRFSGLVGALGGGLSTLSGYNGLPSVLSGLDDSLSDV" +
        "AARKKYPFELEKAIHNVMFIQHHMQRQDEFNAEDQDWGFVAMVMDRLFLWLFMIASLVGT" +
        "FVILGEAPSLYDDTKAIDVQLSDVAKQIYNLTEKKN";

    // EMBOSS sigcleave: "Maximum score 13.739"; mature peptide starts at residue 42.
    // Independent re-derivation: Score = 13.7390400704164 at mature-start index 42.
    private const int Ach2ExpectedCleavagePosition = 42;
    private const double Ach2ExpectedScore = 13.7390400704164;

    // Runner-up site reported by EMBOSS ("(2) 26->38", mature start 39): score 12.135281025149874.
    private const double Ach2RunnerUpScore = 12.135281025149874;

    // The 15-residue scoring window (positions -13..+2 inclusive) ending at -1 = A (residue 41),
    // with +1 = N (residue 42, mature start) and +2 = P (residue 43). EMBOSS shows -13..-1 = "LLVLLLLCETVQA".
    private const string Ach2ExpectedWindow = "LLVLLLLCETVQANP";

    #endregion

    #region MUST: ACH2_DROME reference example

    // M1 — EMBOSS sigcleave worked example: mature protein starts at residue 42.
    [Test]
    public void PredictSignalPeptide_Ach2Drome_ReturnsCleavagePosition42()
    {
        var sp = PredictSignalPeptide(Ach2Drome);

        Assert.That(sp, Is.Not.Null, "ACH2_DROME must yield a best site");
        Assert.That(sp!.Value.CleavagePosition, Is.EqualTo(Ach2ExpectedCleavagePosition),
            "EMBOSS sigcleave reports the mature protein starting at residue 42 (cleavage 41|42)");
    }

    // M2 — von Heijne weight at the best site equals EMBOSS's "Maximum score 13.739".
    [Test]
    public void PredictSignalPeptide_Ach2Drome_ScoreMatchesEmboss()
    {
        var sp = PredictSignalPeptide(Ach2Drome);

        Assert.That(sp, Is.Not.Null, "ACH2_DROME must yield a best site");
        Assert.That(sp!.Value.Score, Is.EqualTo(Ach2ExpectedScore).Within(1e-3),
            "Sum of log-odds weights over positions -13..+2 must equal EMBOSS maximum score 13.739");
    }

    // M3 — SignalSequence is residues 1..41; WindowSequence is the 15-residue -13..+2 window.
    [Test]
    public void PredictSignalPeptide_Ach2Drome_SignalAndWindowSequences()
    {
        var sp = PredictSignalPeptide(Ach2Drome);

        Assert.That(sp, Is.Not.Null, "ACH2_DROME must yield a best site");
        Assert.Multiple(() =>
        {
            Assert.That(sp!.Value.SignalSequence, Is.EqualTo(Ach2Drome.Substring(0, 41)),
                "Signal peptide = residues 1..CleavagePosition-1 (the 41 residues before the mature start)");
            Assert.That(sp.Value.WindowSequence, Is.EqualTo(Ach2ExpectedWindow),
                "Scoring window -13..+2 ends at -1=A (res 41); EMBOSS shows -13..-1 = LLVLLLLCETVQA");
        });
    }

    // M4 — the returned site is the global argmax: strictly higher than the known runner-up site.
    [Test]
    public void PredictSignalPeptide_Ach2Drome_IsGlobalArgmax()
    {
        var sp = PredictSignalPeptide(Ach2Drome);

        Assert.That(sp, Is.Not.Null, "ACH2_DROME must yield a best site");
        Assert.That(sp!.Value.Score, Is.GreaterThan(Ach2RunnerUpScore),
            "Best score 13.739 must strictly exceed the EMBOSS runner-up (mature start 39) score 12.135");
    }

    // M5 — Score >= 3.5 default threshold => flagged as a likely signal peptide.
    [Test]
    public void PredictSignalPeptide_Ach2Drome_IsLikelySignalPeptide()
    {
        var sp = PredictSignalPeptide(Ach2Drome);

        Assert.That(sp, Is.Not.Null, "ACH2_DROME must yield a best site");
        Assert.That(sp!.Value.IsLikelySignalPeptide, Is.True,
            "Score 13.739 >= default minWeight 3.5 → likely signal peptide (EMBOSS acceptance level)");
    }

    #endregion

    #region MUST: log-odds formula and preconditions

    // M10 — Score equals the hand-summed log-odds for a controlled 15-residue window.
    // The single candidate site that produces a full window for a 15-residue input is index 13
    // (positions -13..-1 map to indices 0..12, +1 to index 13, +2 out of range and skipped).
    [Test]
    public void PredictSignalPeptide_LogOddsFormula_MatchesHandComputation()
    {
        // 15 residues: AAAAAAAAAAAGAN  -> deliberately uses A (rich at -1, col 12) and G.
        // Independently computed via S(i)=Σ ln(count/expect) using EMBOSS data/Esig.euk.
        const string seq = "AAAAAAAAAAAAGAN";
        var sp = PredictSignalPeptide(seq);

        Assert.That(sp, Is.Not.Null, "15-residue sequence must yield exactly one full-window site");

        double expected = HandScoreEukaryotic(seq, sp!.Value.CleavagePosition - 1);
        Assert.That(sp.Value.Score, Is.EqualTo(expected).Within(1e-9),
            "Score must equal the sum of ln(count/expect) log-odds weights over the window");
    }

    // M6 — null input returns null.
    [Test]
    public void PredictSignalPeptide_NullInput_ReturnsNull()
    {
        Assert.That(PredictSignalPeptide(null!), Is.Null, "Null sequence must return null");
    }

    // M7 — empty input returns null.
    [Test]
    public void PredictSignalPeptide_EmptyInput_ReturnsNull()
    {
        Assert.That(PredictSignalPeptide(""), Is.Null, "Empty sequence must return null");
    }

    // M8 — sequence shorter than one full 15-residue window returns null.
    [Test]
    public void PredictSignalPeptide_ShorterThanWindow_ReturnsNull()
    {
        Assert.That(PredictSignalPeptide("MKTLLLTLVVVTLV"), Is.Null,
            "A 14-residue sequence is below the 15-residue window minimum → null");
    }

    // M9 — exactly 15 residues yields a (single-window) result.
    [Test]
    public void PredictSignalPeptide_ExactlyWindowLength_ReturnsResult()
    {
        var sp = PredictSignalPeptide("AAAAAAAAAAAAGAN");

        Assert.That(sp, Is.Not.Null, "A 15-residue sequence has exactly one full window and must yield a result");
        Assert.That(sp!.Value.CleavagePosition, Is.InRange(1, 15),
            "INV-02: cleavage position must be within [1, length]");
    }

    #endregion

    #region SHOULD

    // S1 — case-insensitive: lower-case input produces identical site and score.
    [Test]
    public void PredictSignalPeptide_CaseInsensitive()
    {
        var upper = PredictSignalPeptide(Ach2Drome);
        var lower = PredictSignalPeptide(Ach2Drome.ToLowerInvariant());

        Assert.That(upper, Is.Not.Null);
        Assert.That(lower, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(lower!.Value.CleavagePosition, Is.EqualTo(upper!.Value.CleavagePosition),
                "Cleavage position must be identical regardless of case");
            Assert.That(lower.Value.Score, Is.EqualTo(upper!.Value.Score).Within(1e-12),
                "Score must be identical regardless of case (input is upper-cased)");
        });
    }

    // S2 — raising minWeight above the score flips the likelihood flag but not the selection.
    [Test]
    public void PredictSignalPeptide_MinWeightAboveScore_NotLikelyButSameSite()
    {
        var sp = PredictSignalPeptide(Ach2Drome, minWeight: 14.0);

        Assert.That(sp, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(sp!.Value.CleavagePosition, Is.EqualTo(Ach2ExpectedCleavagePosition),
                "Threshold does not change the argmax site");
            Assert.That(sp.Value.Score, Is.EqualTo(Ach2ExpectedScore).Within(1e-3),
                "Threshold does not change the score");
            Assert.That(sp.Value.IsLikelySignalPeptide, Is.False,
                "Score 13.739 < minWeight 14.0 → not flagged as likely");
        });
    }

    // S3 — prokaryotic matrix gives a distinct score (different matrix, data/Esig.pro).
    [Test]
    public void PredictSignalPeptide_ProkaryoticMatrix_DiffersFromEukaryotic()
    {
        var euk = PredictSignalPeptide(Ach2Drome, prokaryote: false);
        var pro = PredictSignalPeptide(Ach2Drome, prokaryote: true);

        Assert.That(euk, Is.Not.Null);
        Assert.That(pro, Is.Not.Null);
        Assert.That(pro!.Value.Score, Is.Not.EqualTo(euk!.Value.Score).Within(1e-6),
            "Prokaryotic (Esig.pro) and eukaryotic (Esig.euk) matrices produce different scores");
    }

    #endregion

    #region COULD

    // C1 — non-standard residues (X) contribute 0 and must not throw.
    [Test]
    public void PredictSignalPeptide_NonStandardResidue_NoThrow()
    {
        // X is not in the 20-residue matrix; it must be skipped (contributes 0), not crash.
        string seq = "AAAAAXAAAAAAGAN";
        SignalPeptide? sp = null;

        Assert.DoesNotThrow(() => sp = PredictSignalPeptide(seq),
            "Non-standard residue X must be tolerated (skipped), not throw");
        Assert.That(sp, Is.Not.Null, "A 15-residue sequence still yields a windowed result");
    }

    #endregion

    #region Helpers (independent re-derivation of the EMBOSS eukaryotic score)

    // Eukaryotic count matrix and Expect column, verbatim from EMBOSS 6.6.0 data/Esig.euk.
    private static readonly string Residues = "ACDEFGHIKLMNPQRSTVWY";

    private static readonly int[][] Counts =
    {
        new[] { 16, 13, 14, 15, 20, 18, 18, 17, 25, 15, 47,  6, 80, 18,  6 }, // A
        new[] {  3,  6,  9,  7,  9, 14,  6,  8,  5,  6, 19,  3,  9,  8,  3 }, // C
        new[] {  0,  0,  0,  0,  0,  0,  0,  0,  5,  3,  0,  5,  0, 10, 11 }, // D
        new[] {  0,  0,  0,  1,  0,  0,  0,  0,  3,  7,  0,  7,  0, 13, 14 }, // E
        new[] { 13,  9, 11, 11,  6,  7, 18, 13,  4,  5,  0, 13,  0,  6,  4 }, // F
        new[] {  4,  4,  3,  6,  3, 13,  3,  2, 19, 34,  5,  7, 39, 10,  7 }, // G
        new[] {  0,  0,  0,  0,  0,  1,  1,  0,  5,  0,  0,  6,  0,  4,  2 }, // H
        new[] { 15, 15,  8,  6, 11,  5,  4,  8,  5,  1, 10,  5,  0,  8,  7 }, // I
        new[] {  0,  0,  0,  1,  0,  0,  1,  0,  0,  4,  0,  2,  0, 11,  9 }, // K
        new[] { 71, 68, 72, 79, 78, 45, 64, 49, 10, 23,  8, 20,  1,  8,  4 }, // L
        new[] {  0,  3,  7,  4,  1,  6,  2,  2,  0,  0,  0,  1,  0,  1,  2 }, // M
        new[] {  0,  1,  0,  1,  1,  0,  0,  0,  3,  3,  0, 10,  0,  4,  7 }, // N
        new[] {  2,  0,  2,  0,  0,  4,  1,  8, 20, 14,  0,  1,  3,  0, 22 }, // P
        new[] {  0,  0,  0,  1,  0,  6,  1,  0, 10,  8,  0, 18,  3, 19, 10 }, // Q
        new[] {  2,  0,  0,  0,  0,  1,  0,  0,  7,  4,  0, 15,  0, 12,  9 }, // R
        new[] {  9,  3,  8,  6, 13, 10, 15, 16, 26, 11, 23, 17, 20, 15, 10 }, // S
        new[] {  2, 10,  5,  4,  5, 13,  7,  7, 12,  6, 17,  8,  6,  3, 10 }, // T
        new[] { 20, 25, 15, 18, 13, 15, 11, 27,  0, 12, 32,  3,  0,  8, 17 }, // V
        new[] {  4,  3,  3,  1,  1,  2,  6,  3,  1,  3,  0,  9,  0,  2,  0 }, // W
        new[] {  0,  1,  4,  0,  0,  1,  3,  1,  1,  2,  0,  5,  0,  1,  7 }, // Y
    };

    private static readonly double[] Expect =
    {
        14.5, 4.5, 8.9, 10.0, 5.6, 12.1, 3.4, 7.4, 11.3, 12.1,
        2.7, 7.1, 7.4, 6.3, 7.6, 11.4, 9.7, 11.1, 1.8, 5.6,
    };

    /// <summary>
    /// Independent re-implementation of the EMBOSS log-odds score for candidate +1 index <paramref name="i"/>,
    /// summing ln(count/expect) over the 15-column window with the -3/-1 zero-count penalty (1e-10) and
    /// 1.0 elsewhere. Mirrors data/Esig.euk + sigcleave.c without sharing the production code path.
    /// </summary>
    private static double HandScoreEukaryotic(string seq, int i)
    {
        double total = 0.0;
        for (int col = 0; col < 15; col++)
        {
            int j = i - 13 + col;
            if (j < 0 || j >= seq.Length)
                continue;

            int row = Residues.IndexOf(char.ToUpperInvariant(seq[j]));
            if (row < 0)
                continue;

            double c = Counts[row][col];
            if (c == 0.0)
                c = (col == 10 || col == 12) ? 1.0e-10 : 1.0;

            total += Math.Log(c / Expect[row]);
        }

        return total;
    }

    #endregion
}
