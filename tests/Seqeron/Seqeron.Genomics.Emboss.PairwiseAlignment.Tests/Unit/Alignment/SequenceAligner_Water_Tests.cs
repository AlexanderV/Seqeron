// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-LOCAL-001 (audit A2-28, B09 F122): <see cref="EmbossPairwisePrograms.WaterAlign"/> = EMBOSS 6.6.0 <c>water</c> (<c>embAlignPathCalcSW</c> + <c>embAlignWalkSWMatrix</c> in the program's arithmetic) against the installed binary:
/// TestData/Alignment/water_oracle.tsv (scripts/alignment/generate_water_oracle.py), 860 random DNA / IUPAC / protein pairs (2 - 200 residues, integer / dyadic / non-dyadic gap penalties, 10 matrices): rows, start coordinates and
/// printed score are identical on every pair.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-LOCAL-001")]
public class SequenceAligner_Water_Tests
{
    private sealed record Row(string Id, string Matrix, double GapOpen, double GapExtend, string Seq1, string Seq2, string Row1, string Row2, int Start1, int Start2, string Score);


    [Test]
    public void WaterAlign_RowsStartsAndScore_AreIdenticalToEmbossOnAllOracleCases()
    {
        var bad = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var a = EmbossPairwisePrograms.WaterAlign(r.Seq1, r.Seq2, SubstitutionMatrix.Get(r.Matrix), r.GapOpen, r.GapExtend);
            if (a.AlignedSequence1 != r.Row1 || a.AlignedSequence2 != r.Row2 || a.StartPosition1 != r.Start1 || a.StartPosition2 != r.Start2 || PrintedScore(a.ExactScore) != r.Score)
                bad.Add($"{r.Id} {r.Matrix} {r.GapOpen}/{r.GapExtend}");
        }

        Assert.That(Oracle.Value, Has.Count.GreaterThanOrEqualTo(100));
        Assert.That(bad, Is.Empty, $"{bad.Count} of {Oracle.Value.Count}: {string.Join(", ", bad.Take(15))}");
    }

    /// <summary>
    /// The score of the program is the maximum of ITS recurrence (<c>embAlignPathCalcSW</c>): the running maximum <c>ret</c> is taken only over the cells whose row AND column index are at least 1 (the first
    /// row and column are initialised in separate loops and never enter the maximum), before the floor at zero. <see cref="EmbossReferenceDp.Water"/> is a score-only program written from that C text: it gives
    /// the printed score of every oracle row with penalties that are multiples of 1/8, and - with the border cells added to the maximum - the textbook local optimum (the next test).
    /// </summary>
    [Test]
    public void WaterAlign_Score_EqualsTheIndependentRecurrenceOfTheProgram_OnEveryDyadicOracleRow()
    {
        var bad = new List<string>();
        int tested = 0;
        foreach (var r in Oracle.Value.Where(x => x.Seq1.Length > 1 && x.Seq2.Length > 1 && Dyadic(x.GapOpen) && Dyadic(x.GapExtend)))
        {
            tested++;
            var m = SubstitutionMatrix.Get(r.Matrix);
            if (Math.Abs(EmbossReferenceDp.Water(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend) - EmbossPairwisePrograms.WaterAlign(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend).ExactScore) > 1e-9)
                bad.Add(r.Id);
        }

        Assert.That(tested, Is.GreaterThanOrEqualTo(600));
        Assert.That(bad, Is.Empty, $"{bad.Count} of {tested}: {string.Join(", ", bad.Take(15))}");
    }

    [Test]
    public void WaterAlign_Score_IsTheTextbookLocalOptimum_ExceptTheEnumeratedBorderRows()
    {
        var differ = new List<string>();
        int total = 0;
        foreach (var r in Oracle.Value.Where(x => x.GapExtend > 0 && x.GapOpen >= x.GapExtend && Dyadic(x.GapOpen) && Dyadic(x.GapExtend)).Take(250))
        {
            total++;
            var m = SubstitutionMatrix.Get(r.Matrix);
            var w = EmbossPairwisePrograms.WaterAlign(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend);
            var l = SequenceAligner.LocalAlignAffine(r.Seq1, r.Seq2, MatrixScoring.Emboss(m, r.GapOpen, r.GapExtend));
            // the textbook optimum is the program's recurrence with the border cells in the maximum
            Assert.That(EmbossReferenceDp.Water(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend, includeBorderCells: true), Is.EqualTo(l.ExactScore).Within(1e-9), r.Id);
            if (Math.Abs(w.ExactScore - l.ExactScore) > 1e-6)
            {
                differ.Add(r.Id);
                Assert.That(w.ExactScore, Is.LessThan(l.ExactScore), r.Id);
            }
        }

        Assert.That(total, Is.EqualTo(250));
        Assert.That(differ, Is.EqualTo(DifferFromTextbookLocalOptimum));
    }

    [Test]
    public void WaterAlign_Contract_RangesNullsAndShortInput()
    {
        Assert.Throws<ArgumentNullException>(() => EmbossPairwisePrograms.WaterAlign(null!, "ACGT"));
        Assert.Throws<ArgumentNullException>(() => EmbossPairwisePrograms.WaterAlign("ACGT", null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.WaterAlign("ACGT", "ACGT", null, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossPairwisePrograms.WaterAlign("ACGT", "ACGT", null, 10, 11));
        // a one-residue sequence (water 6.6.0 `water a.fa b.fa -aformat3 pair`): the program reports the best single column of the one-row path matrix with the score -FLT_MAX (A5-7 measured: ACGT / G -> G at 3, ACGT / T -> T at 4,
        // TTTT / A -> T,A at 1; the earlier "shorter input is empty" reading was wrong)
        var one = EmbossPairwisePrograms.WaterAlign("A", "ACGT");
        Assert.That((one.AlignedSequence1, one.AlignedSequence2, one.StartPosition1, one.StartPosition2), Is.EqualTo(("A", "A", 0, 0)));
        Assert.That(one.ExactScore, Is.EqualTo(-float.MaxValue));
        var g = EmbossPairwisePrograms.WaterAlign("ACGT", "G");
        Assert.That((g.AlignedSequence1, g.AlignedSequence2, g.StartPosition1), Is.EqualTo(("G", "G", 2)));
        Assert.That(EmbossPairwisePrograms.WaterAlign("", "ACGT").AlignedSequence1, Is.Empty);
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.water_oracle.tsv";
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
            list.Add(new Row(f[0], f[1], double.Parse(f[2], inv), double.Parse(f[3], inv), f[4], f[5], f[6], f[7], int.Parse(f[8], inv), int.Parse(f[9], inv), f[10]));
        }

        return list;
    });

    private static string PrintedScore(double exact)
    {
        string t = CFormat.Fixed((float)exact, 3);
        for (int i = 1; i < 3 && t[^1] == '0'; i++)
            t = t[..^1];
        return t;
    }

    /// <summary>
    /// Rows (of the first 250 sample rows) on which water's score differs from the textbook affine local optimum (<c>LocalAlignAffine</c>): exactly the committed ids. In each of them EVERY optimal textbook alignment
    /// ends in the first row or column (the textbook maximum is the program's with the border cells added, <c>includeBorderCells</c>) or the matrix has no positive cell (the program reports the negative
    /// maximum of the inner cells, the textbook optimum is 0).
    /// </summary>
    private static readonly string[] DifferFromTextbookLocalOptimum = ["23", "51", "69", "240", "292", "305", "405"];

    private static bool Dyadic(double v) => Math.Abs((v * 8) - Math.Round(v * 8)) < 1e-9;
}
