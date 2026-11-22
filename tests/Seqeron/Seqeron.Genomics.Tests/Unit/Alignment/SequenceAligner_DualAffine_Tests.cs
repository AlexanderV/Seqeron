using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit A2-29, B09 F123): ksw2-style two-piece affine gap alignment (<see cref="SequenceAligner.GlobalAlignDualAffine"/>, <see cref="SequenceAligner.LocalAlignDualAffine"/>, optional band).
/// TestData/Alignment/dual_affine_oracle.tsv (scripts/alignment/generate_dual_affine_oracle.py): 1 600 cases (400 global, 400 local, 400 banded global, 400 banded local; the full-matrix scores are Biopython 1.88 PairwiseAligner with
/// the callable gap score max(o1 + k e1, o2 + k e2), asserted equal to an independent Waterman-Smith-Beyer program that also provides the banded scores); every score is equal and every returned alignment re-scores to it.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class SequenceAligner_DualAffine_Tests
{
    private sealed record Row(string Id, string Mode, int Band, string Seq1, string Seq2, DualAffineScoring Scoring, int Score);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.dual_affine_oracle.tsv";
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
            list.Add(new Row(f[0], f[1], int.Parse(f[2], inv), f[3], f[4], new DualAffineScoring(int.Parse(f[5], inv), int.Parse(f[6], inv), int.Parse(f[7], inv), int.Parse(f[8], inv), int.Parse(f[9], inv), int.Parse(f[10], inv)), int.Parse(f[11], inv)));
        }

        return list;
    });

    private static long Rescore(AlignmentResult r, DualAffineScoring sc)
    {
        long total = 0;
        int run = 0;
        char kind = ' ';
        for (int i = 0; i < r.AlignedSequence1.Length; i++)
        {
            char x = r.AlignedSequence1[i];
            char y = r.AlignedSequence2[i];
            if (x == '-' || y == '-')
            {
                char t = x == '-' ? 'h' : 'v';
                if (t != kind && run > 0)
                {
                    total += sc.GapScore(run);
                    run = 0;
                }

                kind = t;
                run++;
            }
            else
            {
                if (run > 0)
                    total += sc.GapScore(run);
                run = 0;
                kind = ' ';
                total += x == y ? sc.Match : sc.Mismatch;
            }
        }

        return run > 0 ? total + sc.GapScore(run) : total;
    }

    private static AlignmentResult Run(Row r) => r.Mode == "global"
        ? SequenceAligner.GlobalAlignDualAffine(r.Seq1, r.Seq2, r.Scoring, r.Band)
        : SequenceAligner.LocalAlignDualAffine(r.Seq1, r.Seq2, r.Scoring, r.Band);

    [Test]
    public void DualAffine_Scores_EqualTheOracleOnAllCases_AndEveryAlignmentRescoresToIt()
    {
        var bad = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var a = Run(r);
            string ungapped1 = a.AlignedSequence1.Replace("-", "");
            string ungapped2 = a.AlignedSequence2.Replace("-", "");
            // R14-1: End is 0-based INCLUSIVE; an empty local result is -1 / -1
            bool coordsOk = ungapped1.Length == 0 && ungapped2.Length == 0 && a.StartPosition1 == (r.Mode == "local" ? -1 : 0) && a.EndPosition1 == (r.Mode == "local" ? -1 : r.Seq1.Length - 1)
                || (ungapped1.Length > 0 || ungapped2.Length > 0) && r.Seq1.ToUpperInvariant().Substring(a.StartPosition1, a.EndPosition1 - a.StartPosition1 + 1) == ungapped1
                && r.Seq2.ToUpperInvariant().Substring(a.StartPosition2, a.EndPosition2 - a.StartPosition2 + 1) == ungapped2;
            if (a.Score != r.Score || Rescore(a, r.Scoring) != r.Score || !coordsOk)
                bad.Add($"{r.Id} {r.Mode} band {r.Band}");
        }

        foreach (string mode in new[] { "global", "local" })
        {
            Assert.That(Oracle.Value.Count(r => r.Mode == mode && r.Band < 0), Is.GreaterThanOrEqualTo(300), mode + " full");
            Assert.That(Oracle.Value.Count(r => r.Mode == mode && r.Band >= 0), Is.GreaterThanOrEqualTo(300), mode + " banded");
        }

        Assert.That(bad, Is.Empty, $"{bad.Count} of {Oracle.Value.Count}: {string.Join(", ", bad.Take(15))}");
    }

    [Test]
    public void DualAffine_WithOneRegime_EqualsTheSingleAffineAligner()
    {
        var sc = new DualAffineScoring(2, -3, -5, -2, -5, -2);
        var m = new ScoringMatrix(2, -3, -5, -2);
        foreach (var (a, b) in new[] { ("ACGTACGTAGCTAGCTAGC", "ACGTAGTAGCTAGGCTAGC"), ("AAAAACCCCCGGGGG", "AAAACCCCGGGGG"), ("ACGT", "TGCA") })
        {
            Assert.That(SequenceAligner.GlobalAlignDualAffine(a, b, sc).Score, Is.EqualTo(SequenceAligner.GlobalAlignAffine(a, b, m).Score), a);
            Assert.That(SequenceAligner.LocalAlignDualAffine(a, b, sc).Score, Is.EqualTo(SequenceAligner.LocalAlignAffine(a, b, m).Score), a);
        }
    }

    [Test]
    public void DualAffine_LongGap_UsesTheCheaperLongRegime()
    {
        // a 30-residue deletion: min(4 + 2k, 24 + k) = 54 with the long regime, 64 with the short one
        string left = "ACGTTGCATGCCGTAAGCTTGACCAGT";
        string right = "TTGACCGGATCCATGCAAGTCGTACGT";
        string middle = "GATTACAGATTACAGATTACAGATTACAGA";
        var r = SequenceAligner.GlobalAlignDualAffine(left + middle + right, left + right, DualAffineScoring.Minimap2);
        Assert.That(r.Score, Is.EqualTo((2 * (left.Length + right.Length)) - 54));
    }

    [Test]
    public void DualAffine_Contract()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlignDualAffine(null!, "A", DualAffineScoring.Minimap2));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignDualAffine("A", null!, DualAffineScoring.Minimap2));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlignDualAffine("A", "A", null!));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignDualAffine("A", "A", new DualAffineScoring(1, -1, 1, -1, -2, -1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.GlobalAlignDualAffine("A", "A", DualAffineScoring.Minimap2, -2));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignDualAffine("AAAAAAAA", "AA", DualAffineScoring.Minimap2, 2), "the end cell is outside the band");
        var e = SequenceAligner.GlobalAlignDualAffine("", "ACG", DualAffineScoring.Minimap2);
        Assert.That(e.AlignedSequence1, Is.EqualTo("---"));
        Assert.That(e.Score, Is.EqualTo(-10), "min(4 + 2*3, 24 + 3) = 10");
        Assert.That(SequenceAligner.LocalAlignDualAffine("", "ACG", DualAffineScoring.Minimap2).Score, Is.EqualTo(0));
    }
}
