// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.PairwiseAlignment;
using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 (audit A2-28, B09 F122): <see cref="EmbossPairwisePrograms.NeedleAlign"/> = EMBOSS 6.6.0 <c>needle</c> with its complete option set (<c>-datafile -gapopen -gapextend -endweight -endopen -endextend</c>): the port of
/// <c>embAlignPathCalcWithEndGapPenalties</c> / <c>embAlignWalkNWMatrixUsingCompass</c> in single precision, against the installed binary: TestData/Alignment/needle_oracle.tsv (scripts/alignment/generate_needle_oracle.py), 887 random
/// DNA / IUPAC / protein pairs (1 - 200 residues; integer, dyadic and non-dyadic gap penalties; with and without end-gap weights; 10 matrices): both rows and the printed score are identical on every pair.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class SequenceAligner_Needle_Tests
{
    private sealed record Row(string Id, string Matrix, double GapOpen, double GapExtend, bool EndWeight, double EndOpen, double EndExtend, string Seq1, string Seq2, string Row1, string Row2, string Score);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.needle_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        var inv = CultureInfo.InvariantCulture;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], double.Parse(f[2], inv), double.Parse(f[3], inv), f[4] == "Y", double.Parse(f[5], inv), double.Parse(f[6], inv), f[7], f[8], f[9], f[10], f[11]));
        }

        return list;
    });

    // ajAlignSetScoreR: "%.3f" of the float score, up to two trailing zeros removed
    private static string PrintedScore(double exact)
    {
        string t = CFormat.Fixed((float)exact, 3);
        for (int i = 1; i < 3 && t[^1] == '0'; i++)
            t = t[..^1];
        return t;
    }

    [Test]
    public void NeedleAlign_RowsAndScore_AreIdenticalToEmbossOnAllOracleCases()
    {
        var bad = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var a = EmbossPairwisePrograms.NeedleAlign(r.Seq1, r.Seq2, SubstitutionMatrix.Get(r.Matrix), r.GapOpen, r.GapExtend, r.EndWeight, r.EndOpen, r.EndExtend);
            if (a.AlignedSequence1 != r.Row1 || a.AlignedSequence2 != r.Row2 || PrintedScore(a.ExactScore) != r.Score)
                bad.Add($"{r.Id} {r.Matrix} {r.GapOpen}/{r.GapExtend} ew={r.EndWeight}");
        }

        Assert.That(Oracle.Value, Has.Count.GreaterThanOrEqualTo(100));
        Assert.That(Oracle.Value.Count(r => r.EndWeight), Is.GreaterThanOrEqualTo(100));
        Assert.That(Oracle.Value.Count(r => !r.EndWeight), Is.GreaterThanOrEqualTo(100));
        Assert.That(bad, Is.Empty, $"{bad.Count} of {Oracle.Value.Count}: {string.Join(", ", bad.Take(15))}");
    }

    private static bool Dyadic(double v) => Math.Abs((v * 8) - Math.Round(v * 8)) < 1e-9;

    /// <summary>
    /// The score of the program is the score of ITS recurrence (<c>embAlignPathCalcWithEndGapPenalties</c>), not of the textbook one: an independent score-only dynamic program written from the C text
    /// (<see cref="EmbossReferenceDp.Needle"/>) gives the printed score of every oracle row whose penalties are multiples of 1/8 (exact in plain doubles) - the rows are equal to the oracle in
    /// <see cref="NeedleAlign_RowsAndScore_AreIdenticalToEmbossOnAllOracleCases"/>, this test says WHY the plain affine optimum is sometimes higher (the next two tests).
    /// </summary>
    [Test]
    public void NeedleAlign_Score_EqualsTheIndependentRecurrenceOfTheProgram_OnEveryDyadicOracleRow()
    {
        var bad = new List<string>();
        int tested = 0;
        foreach (var r in Oracle.Value.Where(x => x.Seq1.Length > 1 && x.Seq2.Length > 1 && Dyadic(x.GapOpen) && Dyadic(x.GapExtend) && Dyadic(x.EndOpen) && Dyadic(x.EndExtend)))
        {
            tested++;
            double reference = EmbossReferenceDp.Needle(r.Seq1, r.Seq2, SubstitutionMatrix.Get(r.Matrix), r.GapOpen, r.GapExtend, r.EndWeight, r.EndOpen, r.EndExtend);
            if (Math.Abs(reference - EmbossPairwisePrograms.NeedleAlign(r.Seq1, r.Seq2, SubstitutionMatrix.Get(r.Matrix), r.GapOpen, r.GapExtend, r.EndWeight, r.EndOpen, r.EndExtend).ExactScore) > 1e-9)
                bad.Add(r.Id);
        }

        Assert.That(tested, Is.GreaterThanOrEqualTo(400));
        Assert.That(bad, Is.Empty, $"{bad.Count} of {tested}: {string.Join(", ", bad.Take(15))}");
    }

    /// <summary>Rows (of the first 250 / 200 sample rows) on which needle's score is BELOW the plain affine optimum: exactly the committed ids, each reproduced by the program's recurrence (previous test), none by a tolerance.</summary>
    private static readonly string[] BelowPlainOptimumFreeEnds = ["234"];

    private static readonly string[] BelowPlainOptimumEndWeights = ["234"];

    [Test]
    public void NeedleAlign_DefaultEndGaps_AreFree_ScoreIsTheOverlapAffineOptimum_ExceptTheEnumeratedRows()
    {
        var below = new List<string>();
        int total = 0;
        foreach (var r in Oracle.Value.Where(x => !x.EndWeight && x.GapExtend > 0 && x.GapOpen >= x.GapExtend && x.Seq1.Length > 1 && x.Seq2.Length > 1 && Dyadic(x.GapOpen) && Dyadic(x.GapExtend)).Take(250))
        {
            total++;
            var m = SubstitutionMatrix.Get(r.Matrix);
            var n = EmbossPairwisePrograms.NeedleAlign(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend);
            var g = SequenceAligner.SemiGlobalAlignAffine(r.Seq1, r.Seq2, MatrixScoring.Emboss(m, r.GapOpen, r.GapExtend), EndGaps.All);
            Assert.That(n.ExactScore, Is.LessThanOrEqualTo(g.ExactScore + 1e-6), r.Id);
            if (n.ExactScore < g.ExactScore - 1e-6)
                below.Add(r.Id);
        }

        Assert.That(total, Is.EqualTo(235));
        Assert.That(below, Is.EqualTo(BelowPlainOptimumFreeEnds));
    }

    [Test]
    public void NeedleAlign_EndWeightWithTheGapPenalties_IsTheGlobalAffineOptimum_ExceptTheEnumeratedRows()
    {
        // end gaps charged exactly like inner gaps = the plain global Gotoh optimum (EMBOSS needle mis-initialises the boundary of a one-residue sequence: score 0, excluded)
        var below = new List<string>();
        int total = 0;
        foreach (var r in Oracle.Value.Where(x => x.GapExtend > 0 && x.GapOpen >= x.GapExtend && x.Seq1.Length > 1 && x.Seq2.Length > 1 && Dyadic(x.GapOpen) && Dyadic(x.GapExtend)).Take(200))
        {
            total++;
            var m = SubstitutionMatrix.Get(r.Matrix);
            var n = EmbossPairwisePrograms.NeedleAlign(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend, endWeight: true, endOpen: r.GapOpen, endExtend: r.GapExtend);
            var g = SequenceAligner.GlobalAlignAffine(r.Seq1, r.Seq2, MatrixScoring.Emboss(m, r.GapOpen, r.GapExtend));
            Assert.That(n.ExactScore, Is.LessThanOrEqualTo(g.ExactScore + 1e-6), r.Id);
            if (n.ExactScore < g.ExactScore - 1e-6)
                below.Add(r.Id);
        }

        Assert.That(total, Is.EqualTo(200));
        Assert.That(below, Is.EqualTo(BelowPlainOptimumEndWeights));
    }

    [Test]
    public void RoundFloat_ReproducesAjRoundFloat_ToSixteenMantissaBits()
    {
        Assert.That(EmbossPairwisePrograms.RoundFloat(10f, 8), Is.EqualTo(10f));
        Assert.That(EmbossPairwisePrograms.RoundFloat(0.5f, 8), Is.EqualTo(0.5f));
        Assert.That(EmbossPairwisePrograms.RoundFloat(0.1f, 8), Is.EqualTo(0.10000038f));
        Assert.That(EmbossPairwisePrograms.RoundFloat(0f, 8), Is.EqualTo(0f));
    }

    [Test]
    public void NeedleAlign_Contract_RangesNullsAndEmpty()
    {
        Assert.Throws<ArgumentNullException>(() => EmbossPairwisePrograms.NeedleAlign(null!, "ACGT"));
        Assert.Throws<ArgumentNullException>(() => EmbossPairwisePrograms.NeedleAlign("ACGT", null!));
        Assert.Throws<ArgumentException>(() => EmbossPairwisePrograms.NeedleAlign("", "ACGT"));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAlign("ACGT", "ACGT", null, 0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAlign("ACGT", "ACGT", null, 10, 11));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAlign("ACGT", "ACGT", null, 10, 0.5, true, 0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.NeedleAlign("ACGT", "ACGT", null, 10, 0.5, true, 10, -1));
    }

    [Test]
    public void NeedleAlign_DocumentedExample_HbaVersusHbb_IsTheNeedleOutput()
    {
        // the first alignment of the needle documentation / the oracle file: 34 columns, score 146.0
        var a = EmbossPairwisePrograms.NeedleAlign("ACGTACGTAGCTAGCTAGCATCGATCGATCGAT", "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT");
        Assert.That(a.AlignedSequence1, Is.EqualTo("ACGTACGTAGCTAGCTAGCATCG-ATCGATCGAT"));
        Assert.That(a.AlignedSequence2, Is.EqualTo("ACGTACGTAGCTAGATAGCATCGAATCGATCGAT"));
        Assert.That(a.ExactScore, Is.EqualTo(146.0));
    }
}
