using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 (audit A2-28, B09 F122): <see cref="SequenceAligner.NeedleAlign"/> = EMBOSS 6.6.0 <c>needle</c> with its complete option set (<c>-datafile -gapopen -gapextend -endweight -endopen -endextend</c>): the port of
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
            var a = SequenceAligner.NeedleAlign(r.Seq1, r.Seq2, SubstitutionMatrix.Get(r.Matrix), r.GapOpen, r.GapExtend, r.EndWeight, r.EndOpen, r.EndExtend);
            if (a.AlignedSequence1 != r.Row1 || a.AlignedSequence2 != r.Row2 || PrintedScore(a.ExactScore) != r.Score)
                bad.Add($"{r.Id} {r.Matrix} {r.GapOpen}/{r.GapExtend} ew={r.EndWeight}");
        }

        Assert.That(Oracle.Value, Has.Count.GreaterThanOrEqualTo(100));
        Assert.That(Oracle.Value.Count(r => r.EndWeight), Is.GreaterThanOrEqualTo(100));
        Assert.That(Oracle.Value.Count(r => !r.EndWeight), Is.GreaterThanOrEqualTo(100));
        Assert.That(bad, Is.Empty, $"{bad.Count} of {Oracle.Value.Count}: {string.Join(", ", bad.Take(15))}");
    }

    [Test]
    public void NeedleAlign_DefaultEndGaps_AreFree_ScoreNeverExceedsTheOverlapAffineOptimum()
    {
        int same = 0;
        int total = 0;
        foreach (var r in Oracle.Value.Where(x => !x.EndWeight && x.GapExtend > 0 && x.GapOpen >= x.GapExtend && x.Seq1.Length > 1 && x.Seq2.Length > 1 && Math.Abs((x.GapOpen * 8) - Math.Round(x.GapOpen * 8)) < 1e-9 && Math.Abs((x.GapExtend * 8) - Math.Round(x.GapExtend * 8)) < 1e-9).Take(250))
        {
            total++;
            var m = SubstitutionMatrix.Get(r.Matrix);
            var n = SequenceAligner.NeedleAlign(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend);
            var g = SequenceAligner.SemiGlobalAlignAffine(r.Seq1, r.Seq2, MatrixScoring.Emboss(m, r.GapOpen, r.GapExtend), EndGaps.All);
            Assert.That(n.ExactScore, Is.LessThanOrEqualTo(g.ExactScore + 1e-6), r.Id);
            if (Math.Abs(n.ExactScore - g.ExactScore) < 1e-6)
                same++;
        }
        Assert.That(same, Is.GreaterThanOrEqualTo(total * 9 / 10), "EMBOSS reaches the plain affine optimum on at least 90 % of the pairs (its recurrences forbid a few gap-state transitions)");
    }

    [Test]
    public void NeedleAlign_EndWeightWithTheGapPenalties_NeverExceedsTheGlobalAffineOptimum()
    {
        // end gaps charged exactly like inner gaps = the plain global Gotoh optimum (EMBOSS needle mis-initialises the boundary of a one-residue sequence: score 0, excluded)
        int checkedCases = 0;
        int same = 0;
        int total = 0;
        foreach (var r in Oracle.Value.Where(x => x.GapExtend > 0 && x.GapOpen >= x.GapExtend && x.Seq1.Length > 1 && x.Seq2.Length > 1 && Math.Abs((x.GapOpen * 8) - Math.Round(x.GapOpen * 8)) < 1e-9 && Math.Abs((x.GapExtend * 8) - Math.Round(x.GapExtend * 8)) < 1e-9).Take(200))
        {
            total++;
            var m = SubstitutionMatrix.Get(r.Matrix);
            var n = SequenceAligner.NeedleAlign(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend, endWeight: true, endOpen: r.GapOpen, endExtend: r.GapExtend);
            var g = SequenceAligner.GlobalAlignAffine(r.Seq1, r.Seq2, MatrixScoring.Emboss(m, r.GapOpen, r.GapExtend));
            Assert.That(n.ExactScore, Is.LessThanOrEqualTo(g.ExactScore + 1e-6), r.Id);
            if (Math.Abs(n.ExactScore - g.ExactScore) < 1e-6)
                same++;
            checkedCases++;
        }
        Assert.That(same, Is.GreaterThanOrEqualTo(total * 9 / 10), "EMBOSS reaches the plain affine optimum on at least 90 % of the pairs (its recurrences forbid a few gap-state transitions)");

        Assert.That(checkedCases, Is.GreaterThan(100));
    }

    [Test]
    public void RoundFloat_ReproducesAjRoundFloat_ToSixteenMantissaBits()
    {
        Assert.That(SequenceAligner.RoundFloat(10f, 8), Is.EqualTo(10f));
        Assert.That(SequenceAligner.RoundFloat(0.5f, 8), Is.EqualTo(0.5f));
        Assert.That(SequenceAligner.RoundFloat(0.1f, 8), Is.EqualTo(0.10000038f));
        Assert.That(SequenceAligner.RoundFloat(0f, 8), Is.EqualTo(0f));
    }

    [Test]
    public void NeedleAlign_Contract_RangesNullsAndEmpty()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.NeedleAlign(null!, "ACGT"));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.NeedleAlign("ACGT", null!));
        Assert.Throws<ArgumentException>(() => SequenceAligner.NeedleAlign("", "ACGT"));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAlign("ACGT", "ACGT", null, 0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAlign("ACGT", "ACGT", null, 10, 11));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAlign("ACGT", "ACGT", null, 10, 0.5, true, 0.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.NeedleAlign("ACGT", "ACGT", null, 10, 0.5, true, 10, -1));
    }

    [Test]
    public void NeedleAlign_DocumentedExample_HbaVersusHbb_IsTheNeedleOutput()
    {
        // the first alignment of the needle documentation / the oracle file: 34 columns, score 146.0
        var a = SequenceAligner.NeedleAlign("ACGTACGTAGCTAGCTAGCATCGATCGATCGAT", "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT");
        Assert.That(a.AlignedSequence1, Is.EqualTo("ACGTACGTAGCTAGCTAGCATCG-ATCGATCGAT"));
        Assert.That(a.AlignedSequence2, Is.EqualTo("ACGTACGTAGCTAGATAGCATCGAATCGATCGAT"));
        Assert.That(a.ExactScore, Is.EqualTo(146.0));
    }
}
