using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-LOCAL-001 (audit A2-28, B09 F122): <see cref="SequenceAligner.WaterAlign"/> = EMBOSS 6.6.0 <c>water</c> (<c>embAlignPathCalcSW</c> + <c>embAlignWalkSWMatrix</c> in the program's arithmetic) against the installed binary:
/// TestData/Alignment/water_oracle.tsv (scripts/alignment/generate_water_oracle.py), 860 random DNA / IUPAC / protein pairs (2 - 200 residues, integer / dyadic / non-dyadic gap penalties, 10 matrices): rows, start coordinates and
/// printed score are identical on every pair.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-LOCAL-001")]
public class SequenceAligner_Water_Tests
{
    private sealed record Row(string Id, string Matrix, double GapOpen, double GapExtend, string Seq1, string Seq2, string Row1, string Row2, int Start1, int Start2, string Score);

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

    [Test]
    public void WaterAlign_RowsStartsAndScore_AreIdenticalToEmbossOnAllOracleCases()
    {
        var bad = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var a = SequenceAligner.WaterAlign(r.Seq1, r.Seq2, SubstitutionMatrix.Get(r.Matrix), r.GapOpen, r.GapExtend);
            if (a.AlignedSequence1 != r.Row1 || a.AlignedSequence2 != r.Row2 || a.StartPosition1 != r.Start1 || a.StartPosition2 != r.Start2 || PrintedScore(a.ExactScore) != r.Score)
                bad.Add($"{r.Id} {r.Matrix} {r.GapOpen}/{r.GapExtend}");
        }

        Assert.That(Oracle.Value, Has.Count.GreaterThanOrEqualTo(100));
        Assert.That(bad, Is.Empty, $"{bad.Count} of {Oracle.Value.Count}: {string.Join(", ", bad.Take(15))}");
    }

    [Test]
    public void WaterAlign_Score_NeverExceedsTheEmbossAffineLocalOptimum()
    {
        int same = 0;
        int total = 0;
        foreach (var r in Oracle.Value.Where(x => x.GapExtend > 0 && x.GapOpen >= x.GapExtend && Math.Abs((x.GapOpen * 8) - Math.Round(x.GapOpen * 8)) < 1e-9 && Math.Abs((x.GapExtend * 8) - Math.Round(x.GapExtend * 8)) < 1e-9).Take(250))
        {
            total++;
            var m = SubstitutionMatrix.Get(r.Matrix);
            var w = SequenceAligner.WaterAlign(r.Seq1, r.Seq2, m, r.GapOpen, r.GapExtend);
            var l = SequenceAligner.LocalAlignAffine(r.Seq1, r.Seq2, MatrixScoring.Emboss(m, r.GapOpen, r.GapExtend));
            Assert.That(w.ExactScore, Is.LessThanOrEqualTo(l.ExactScore + 1e-6), r.Id);
            if (Math.Abs(w.ExactScore - l.ExactScore) < 1e-6)
                same++;
        }
        Assert.That(same, Is.GreaterThanOrEqualTo(total * 9 / 10), "EMBOSS reaches the plain affine optimum on at least 90 % of the pairs (its recurrences forbid a few gap-state transitions)");
    }

    [Test]
    public void WaterAlign_Contract_RangesNullsAndShortInput()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.WaterAlign(null!, "ACGT"));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.WaterAlign("ACGT", null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.WaterAlign("ACGT", "ACGT", null, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.WaterAlign("ACGT", "ACGT", null, 10, 11));
        var one = SequenceAligner.WaterAlign("A", "ACGT");
        Assert.That(one.AlignedSequence1, Is.Empty);
        Assert.That(one.Score, Is.EqualTo(0));
    }
}
