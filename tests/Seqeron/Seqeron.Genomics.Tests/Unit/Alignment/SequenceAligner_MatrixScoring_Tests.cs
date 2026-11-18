using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 (audit A1-2, B09 F41): pairwise alignment scored with a
/// <see cref="SubstitutionMatrix"/> (BLOSUM62 / PAM250 / EDNAFULL) and real-valued gap scores through
/// <see cref="MatrixScoring"/> — linear and affine, global / local / semi-global (fitting).
///
/// Oracle (TestData/Alignment/matrix_align_oracle.tsv, scripts/alignment/generate_matrix_align_oracle.py): for random protein,
/// DNA and IUPAC-DNA pairs and fractional gap scores (0.5, 0.1, 0.25, 0.3, ...): Biopython 1.88
/// PairwiseAligner(substitution_matrix=...), EMBOSS needle (-endweight Y) / water 6.6.0, parasail 1.3.4 (integer gaps).
/// Gap convention: GapOpen + k·GapExtend = Biopython open_gap_score + (k-1)·extend_gap_score = EMBOSS gapopen G + (k-1)·E
/// with G = -(GapOpen + GapExtend), E = -GapExtend (MatrixScoring.Emboss).
/// </summary>
[TestFixture]
public class SequenceAligner_MatrixScoring_Tests
{
    private sealed record Row(
        string Mode, string Model, string Matrix, double GapOpen, double GapExtend, string Seq1, string Seq2,
        double Biopython, double? Emboss, double? Parasail);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.matrix_align_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine(); // header
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            double? Opt(string s) => s == "NA" ? null : double.Parse(s, CultureInfo.InvariantCulture);
            rows.Add(new Row(f[0], f[1], f[2], D(f[3]), D(f[4]), f[5], f[6], D(f[7]), Opt(f[8]), Opt(f[9])));
        }
        return rows;
    }

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static SubstitutionMatrix MatrixOf(string name) => name switch
    {
        "BLOSUM62" => SubstitutionMatrix.Blosum62,
        "PAM250" => SubstitutionMatrix.Pam250,
        "EDNAFULL" => SubstitutionMatrix.EdnaFull,
        _ => throw new ArgumentException(name),
    };

    private static AlignmentResult Align(Row r)
    {
        var scoring = new MatrixScoring(MatrixOf(r.Matrix), r.GapOpen, r.GapExtend);
        bool affine = r.Model == "affine";
        return (r.Mode, affine) switch
        {
            ("global", false) => SequenceAligner.GlobalAlign(r.Seq1, r.Seq2, scoring),
            ("global", true) => SequenceAligner.GlobalAlignAffine(r.Seq1, r.Seq2, scoring),
            ("local", false) => SequenceAligner.LocalAlign(r.Seq1, r.Seq2, scoring),
            ("local", true) => SequenceAligner.LocalAlignAffine(r.Seq1, r.Seq2, scoring),
            ("semi", false) => SequenceAligner.SemiGlobalAlign(r.Seq1, r.Seq2, scoring),
            _ => SequenceAligner.SemiGlobalAlignAffine(r.Seq1, r.Seq2, scoring),
        };
    }

    /// <summary>Independent re-scoring of aligned rows: substitution matrix per column, gap runs per row
    /// (affine: GapOpen + k·GapExtend per maximal run; linear: GapExtend per gap position).</summary>
    private static double Rescore(string r1, string r2, SubstitutionMatrix m, double open, double extend, bool affine)
    {
        double total = 0;
        for (int c = 0; c < r1.Length; c++)
        {
            if (r1[c] != '-' && r2[c] != '-')
            {
                total += m.GetScore(r1[c], r2[c]);
                continue;
            }
            bool newRun = c == 0 || (r1[c] == '-' ? r1[c - 1] != '-' : r2[c - 1] != '-');
            total += extend + (affine && newRun ? open : 0);
        }
        return total;
    }

    // --- oracle differential ------------------------------------------------------------------------------

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    [Test]
    public void Oracle_HasThousandsOfRows_AllSixModes_AndFractionalGaps()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(1900));
            Assert.That(rows.Select(r => r.Mode + r.Model).Distinct().Count(), Is.EqualTo(6));
            Assert.That(rows.Select(r => r.Matrix).Distinct().Order().ToArray(), Is.EqualTo(new[] { "BLOSUM62", "EDNAFULL", "PAM250" }));
            Assert.That(rows.Count(r => r.GapExtend != Math.Round(r.GapExtend) || r.GapOpen != Math.Round(r.GapOpen)), Is.GreaterThan(800),
                "fractional gap scores");
        });
    }

    /// <summary>Every row: Seqeron score = Biopython 1.88 score (tolerance 1e-6, Biopython sums doubles) and = parasail
    /// 1.3.4 when the gaps are integers.</summary>
    [Test]
    public void Oracle_Score_EqualsBiopythonAndParasail_AllRows()
    {
        var failures = new List<string>();
        int parasailRows = 0;
        foreach (var r in Oracle.Value)
        {
            var a = Align(r);
            if (Math.Abs(a.ExactScore - r.Biopython) > 1e-6)
                failures.Add($"{r.Mode}/{r.Model}/{r.Matrix} go={r.GapOpen} ge={r.GapExtend} {r.Seq1} {r.Seq2}: {a.ExactScore} vs Biopython {r.Biopython}");
            if (r.Parasail is { } p)
            {
                parasailRows++;
                if (Math.Abs(a.ExactScore - p) > 1e-9)
                    failures.Add($"{r.Mode}/{r.Model}/{r.Matrix} {r.Seq1} {r.Seq2}: {a.ExactScore} vs parasail {p}");
            }
        }
        Assert.That(parasailRows, Is.GreaterThan(500));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    /// <summary>EMBOSS needle (end gaps penalised: -endweight Y) / water 6.6.0 agree with Seqeron on all but a small set of rows;
    /// the disagreements are EMBOSS-internal (embaln.c embAlignPathCalcWithEndGapPenalties opens the last-row/column
    /// end gap only from the match matrix, so needle can return a sub-optimal score; verified: every disagreeing needle row is
    /// lower than the Biopython/parasail/Seqeron optimum) — never a higher global score.</summary>
    [Test]
    public void Oracle_Emboss_NeedleAndWater_AgreeExceptKnownEndGapQuirk()
    {
        int rows = 0, equal = 0, needleHigher = 0;
        foreach (var r in Oracle.Value.Where(r => r.Emboss is not null))
        {
            rows++;
            double ours = Align(r).ExactScore;
            double emboss = r.Emboss!.Value;
            if (Math.Abs(ours - emboss) <= 0.051) equal++;
            else if (r.Mode == "global" && emboss > ours) needleHigher++;
        }
        Assert.Multiple(() =>
        {
            Assert.That(rows, Is.GreaterThan(1000));
            Assert.That((double)equal / rows, Is.GreaterThan(0.98), $"{equal}/{rows}");
            Assert.That(needleHigher, Is.EqualTo(0), "needle never beats the optimum");
        });
    }

    /// <summary>The returned alignment really has the reported score (independent re-scoring), reproduces the inputs and
    /// keeps the coordinate contracts of the integer aligners.</summary>
    [Test]
    public void Oracle_Alignment_RescoresToScore_AndReproducesInputs()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var a = Align(r);
            var m = MatrixOf(r.Matrix);
            bool affine = r.Model == "affine";
            string r1 = a.AlignedSequence1, r2 = a.AlignedSequence2;
            string tag = $"{r.Mode}/{r.Model}/{r.Matrix} {r.Seq1} {r.Seq2}";
            if (r1.Length != r2.Length)
            {
                failures.Add(tag + ": row length");
                continue;
            }

            if (r.Mode == "semi")
            {
                if (r1.Replace("-", "") != r.Seq1 || r2.Replace("-", "") != r.Seq2) failures.Add(tag + ": rows do not reproduce inputs");
                int lead = a.StartPosition2, trail = r.Seq2.Length - 1 - a.EndPosition2;
                string c1 = r1.Substring(lead, r1.Length - lead - trail), c2 = r2.Substring(lead, r2.Length - lead - trail);
                if (Math.Abs(Rescore(c1, c2, m, r.GapOpen, r.GapExtend, affine) - a.ExactScore) > 1e-6) failures.Add(tag + ": rescore");
            }
            else if (r.Mode == "global")
            {
                if (r1.Replace("-", "") != r.Seq1 || r2.Replace("-", "") != r.Seq2) failures.Add(tag + ": rows do not reproduce inputs");
                if (Math.Abs(Rescore(r1, r2, m, r.GapOpen, r.GapExtend, affine) - a.ExactScore) > 1e-6) failures.Add(tag + ": rescore");
            }
            else
            {
                if (a.ExactScore == 0)
                {
                    if (r1.Length != 0) failures.Add(tag + ": non-empty zero local");
                    continue;
                }
                if (r.Seq1.Substring(a.StartPosition1, a.EndPosition1 - a.StartPosition1 + 1) != r1.Replace("-", "")) failures.Add(tag + ": local seq1 interval");
                if (r.Seq2.Substring(a.StartPosition2, a.EndPosition2 - a.StartPosition2 + 1) != r2.Replace("-", "")) failures.Add(tag + ": local seq2 interval");
                if (Math.Abs(Rescore(r1, r2, m, r.GapOpen, r.GapExtend, affine) - a.ExactScore) > 1e-6) failures.Add(tag + ": rescore");
            }
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    // --- sourced fixed values -----------------------------------------------------------------------------

    private const string Hba = "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
    private const string Hbb = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";

    /// <summary>EMBOSS needle 6.6.0 -gapopen 10 -gapextend 0.5 -datafile EBLOSUM62 -endweight Y -endopen 10 -endextend 0.5
    /// of HBA_HUMAN vs HBB_HUMAN: Score 292.5 (= Biopython open_gap_score -10, extend_gap_score -0.5).</summary>
    [Test]
    public void HbaVsHbb_EmbossDefaults_Blosum62_Score292_5()
    {
        var r = SequenceAligner.GlobalAlignAffine(Hba, Hbb, MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.5));
        Assert.Multiple(() =>
        {
            Assert.That(r.ExactScore, Is.EqualTo(292.5));
            Assert.That(r.Score, Is.EqualTo(293), "rounded half away from zero");
            var stats = SequenceAligner.CalculateStatistics(r, SubstitutionMatrix.Blosum62);
            Assert.That(stats.AlignmentLength, Is.EqualTo(r.AlignedSequence1.Length));
        });
    }

    /// <summary>EMBOSS water 6.6.0 -gapopen 10 -gapextend 0.5 EBLOSUM62 of the same pair: Score 293.5.</summary>
    [Test]
    public void HbaVsHbb_Water_Blosum62_Score293_5()
    {
        var r = SequenceAligner.LocalAlignAffine(Hba, Hbb, MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.5));
        Assert.That(r.ExactScore, Is.EqualTo(293.5));
    }

    /// <summary>Hand-checkable: C/C = 9, W/W = 11 (BLOSUM62); "CW" vs "CAW": the only gap has length 1, so the
    /// affine penalty is GapOpen + GapExtend = EMBOSS gapopen (10): 9 + 11 - 10 = 10; linear (0.5 per gap position): 20 - 0.5 = 19.5;
    /// gap of length 2: 20 - (10 + 0.5) = 9.5.</summary>
    [Test]
    public void SingleGap_ConventionIsEmbossGapopenForLengthOne()
    {
        var emboss = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62, 10, 0.5);
        Assert.Multiple(() =>
        {
            Assert.That(emboss.GapOpen, Is.EqualTo(-9.5));
            Assert.That(emboss.GapExtend, Is.EqualTo(-0.5));
            var affine = SequenceAligner.GlobalAlignAffine("CW", "CAW", emboss);
            Assert.That(affine.ExactScore, Is.EqualTo(10.0));
            Assert.That(affine.AlignedSequence1, Is.EqualTo("C-W"));
            Assert.That(affine.AlignedSequence2, Is.EqualTo("CAW"));
            var linear = SequenceAligner.GlobalAlign("CW", "CAW", emboss);
            Assert.That(linear.ExactScore, Is.EqualTo(19.5));
            var two = SequenceAligner.GlobalAlignAffine("CW", "CAAW", emboss);
            Assert.That(two.ExactScore, Is.EqualTo(9.5));
        });
    }

    // --- one DP: identical to the integer aligners --------------------------------------------------------

    /// <summary>A match/mismatch SubstitutionMatrix with integer gaps must reproduce the integer aligners row for row
    /// (same cores, same tie order diagonal &gt; up &gt; left / M &gt; X &gt; Y, same coordinates).</summary>
    [Test]
    public void IntegerMatrix_IntegerGaps_IsIdenticalToScalarAligners_RowForRow()
    {
        var rnd = new Random(1118);
        const string alpha = "ACGT";
        var table = new int[4, 4];
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                table[i, j] = i == j ? 2 : -3;
        var matrix = new SubstitutionMatrix("M2-3", alpha, table);
        var scalar = new ScoringMatrix(2, -3, -5, -2);
        var ms = new MatrixScoring(matrix, -5, -2);

        for (int t = 0; t < 400; t++)
        {
            string a = Random(rnd, alpha, rnd.Next(0, 25)), b = Random(rnd, alpha, rnd.Next(0, 30));
            Assert.Multiple(() =>
            {
                AssertSame(SequenceAligner.GlobalAlign(a, b, scalar), SequenceAligner.GlobalAlign(a, b, ms), "NW " + a + " " + b);
                AssertSame(SequenceAligner.GlobalAlignAffine(a, b, scalar), SequenceAligner.GlobalAlignAffine(a, b, ms), "NWG " + a + " " + b);
                AssertSame(SequenceAligner.LocalAlign(a, b, scalar), SequenceAligner.LocalAlign(a, b, ms), "SW " + a + " " + b);
                AssertSame(SequenceAligner.LocalAlignAffine(a, b, scalar), SequenceAligner.LocalAlignAffine(a, b, ms), "SWG " + a + " " + b);
                AssertSame(SequenceAligner.SemiGlobalAlign(a, b, scalar), SequenceAligner.SemiGlobalAlign(a, b, ms), "SG " + a + " " + b);
                AssertSame(SequenceAligner.SemiGlobalAlignAffine(a, b, scalar), SequenceAligner.SemiGlobalAlignAffine(a, b, ms), "SGG " + a + " " + b);
            });
        }
    }

    private static string Random(Random rnd, string alpha, int n) =>
        new(Enumerable.Range(0, n).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray());

    private static void AssertSame(AlignmentResult legacy, AlignmentResult matrix, string tag)
    {
        Assert.That(matrix, Is.EqualTo(legacy with { ExactScore = matrix.ExactScore }), tag);
        Assert.That(matrix.ExactScore, Is.EqualTo((double)legacy.Score), tag);
    }

    /// <summary>The integer aligners are untouched: ExactScore of their results equals Score.</summary>
    [Test]
    public void IntegerAligners_ExactScore_EqualsScore()
    {
        var r = SequenceAligner.GlobalAlign("ACGTACGT", "ACGTTCGT", SequenceAligner.SimpleDna);
        Assert.That(r.ExactScore, Is.EqualTo((double)r.Score));
        Assert.That(AlignmentResult.Empty.ExactScore, Is.EqualTo(0.0));
    }

    // --- overloads ----------------------------------------------------------------------------------------

    [Test]
    public void ProteinAndDnaSequenceOverloads_EqualStringOverloads()
    {
        var ms = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62);
        var p1 = new ProteinSequence("MVLSPADKTNVKAAWGKVG");
        var p2 = new ProteinSequence("MVHLTPEEKSAVTALWGKV");
        var dm = MatrixScoring.Emboss(SubstitutionMatrix.EdnaFull, 10, 0.5);
        var d1 = new DnaSequence("ACGTACGTTAGC");
        var d2 = new DnaSequence("ACGTTCGTAGC");
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.GlobalAlign(p1, p2, ms), Is.EqualTo(SequenceAligner.GlobalAlign(p1.Sequence, p2.Sequence, ms)));
            Assert.That(SequenceAligner.GlobalAlignAffine(p1, p2, ms), Is.EqualTo(SequenceAligner.GlobalAlignAffine(p1.Sequence, p2.Sequence, ms)));
            Assert.That(SequenceAligner.LocalAlign(p1, p2, ms), Is.EqualTo(SequenceAligner.LocalAlign(p1.Sequence, p2.Sequence, ms)));
            Assert.That(SequenceAligner.LocalAlignAffine(p1, p2, ms), Is.EqualTo(SequenceAligner.LocalAlignAffine(p1.Sequence, p2.Sequence, ms)));
            Assert.That(SequenceAligner.SemiGlobalAlign(p1, p2, ms), Is.EqualTo(SequenceAligner.SemiGlobalAlign(p1.Sequence, p2.Sequence, ms)));
            Assert.That(SequenceAligner.SemiGlobalAlignAffine(p1, p2, ms), Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine(p1.Sequence, p2.Sequence, ms)));
            Assert.That(SequenceAligner.GlobalAlign(d1, d2, dm), Is.EqualTo(SequenceAligner.GlobalAlign(d1.Sequence, d2.Sequence, dm)));
            Assert.That(SequenceAligner.GlobalAlignAffine(d1, d2, dm), Is.EqualTo(SequenceAligner.GlobalAlignAffine(d1.Sequence, d2.Sequence, dm)));
            Assert.That(SequenceAligner.LocalAlign(d1, d2, dm), Is.EqualTo(SequenceAligner.LocalAlign(d1.Sequence, d2.Sequence, dm)));
            Assert.That(SequenceAligner.LocalAlignAffine(d1, d2, dm), Is.EqualTo(SequenceAligner.LocalAlignAffine(d1.Sequence, d2.Sequence, dm)));
            Assert.That(SequenceAligner.SemiGlobalAlign(d1, d2, dm), Is.EqualTo(SequenceAligner.SemiGlobalAlign(d1.Sequence, d2.Sequence, dm)));
            Assert.That(SequenceAligner.SemiGlobalAlignAffine(d1, d2, dm), Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine(d1.Sequence, d2.Sequence, dm)));
        });
    }

    [Test]
    public void LowerCaseInput_IsCaseInsensitive_AndFractionalScaleIsExact()
    {
        var ms = new MatrixScoring(SubstitutionMatrix.Blosum62, -0.3, -0.1);
        var upper = SequenceAligner.GlobalAlignAffine("MVLSPADK", "MVHLTPEK", ms);
        var lower = SequenceAligner.GlobalAlignAffine("mvlspadk", "mvhltpek", ms);
        Assert.Multiple(() =>
        {
            Assert.That(lower, Is.EqualTo(upper));
            // 0.1-resolution gaps: exact decimal, no floating-point drift (0.3 + 0.1·k)
            Assert.That(upper.ExactScore * 10 - Math.Round(upper.ExactScore * 10), Is.EqualTo(0).Within(1e-9));
        });
    }

    // --- validation ---------------------------------------------------------------------------------------

    [Test]
    public void InvalidInput_Throws()
    {
        var ms = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlign((string)null!, "ACD", ms));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlign("ACD", (string)null!, ms));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlign("ACD", "ACD", (MatrixScoring)null!));
            var unknown = Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("AC1D", "ACD", ms));
            Assert.That(unknown!.Message, Does.Contain("'1'").And.Contain("BLOSUM62"));
            Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlign("ACD", "AC-D", ms), "gap symbol in input");
            Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACD", "ACD", new MatrixScoring(SubstitutionMatrix.Blosum62, 1, -1)), "positive gap open");
            Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACD", "ACD", new MatrixScoring(SubstitutionMatrix.Blosum62, 0, double.NaN)));
            Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACD", "ACD", new MatrixScoring(SubstitutionMatrix.Blosum62, 0, -0.0000001)), "finer than 1e-6");
            Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign(new string('A', 3000), new string('A', 3000),
                new MatrixScoring(SubstitutionMatrix.Blosum62, -0.000001, -0.000001)), "integer range");
        });
    }

    [Test]
    public void EmptyInputs_FollowTheIntegerAlignerContracts()
    {
        var ms = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62);
        var global = SequenceAligner.GlobalAlign("", "ACD", ms);
        var local = SequenceAligner.LocalAlign("", "ACD", ms);
        var semi = SequenceAligner.SemiGlobalAlign("", "ACD", ms);
        Assert.Multiple(() =>
        {
            Assert.That(global.AlignedSequence1, Is.EqualTo("---"));
            Assert.That(global.ExactScore, Is.EqualTo(-1.5));
            Assert.That(local.AlignedSequence1, Is.Empty);
            Assert.That(local.AlignmentType, Is.EqualTo(AlignmentType.Local));
            Assert.That(semi.ExactScore, Is.EqualTo(0.0));
            Assert.That(SequenceAligner.GlobalAlign("", "", ms).ExactScore, Is.EqualTo(0.0));
        });
    }

    [Test]
    public void ZeroGapOpen_Affine_EqualsLinear_WithMatrixScoring()
    {
        var ms = new MatrixScoring(SubstitutionMatrix.Pam250, 0, -0.5);
        var rnd = new Random(77);
        for (int i = 0; i < 100; i++)
        {
            string a = Random(rnd, "ACDEFGHIKLMNPQRSTVWY", rnd.Next(5, 40)), b = Random(rnd, "ACDEFGHIKLMNPQRSTVWY", rnd.Next(5, 40));
            Assert.That(SequenceAligner.GlobalAlignAffine(a, b, ms).ExactScore, Is.EqualTo(SequenceAligner.GlobalAlign(a, b, ms).ExactScore));
            Assert.That(SequenceAligner.LocalAlignAffine(a, b, ms).ExactScore, Is.EqualTo(SequenceAligner.LocalAlign(a, b, ms).ExactScore));
            Assert.That(SequenceAligner.SemiGlobalAlignAffine(a, b, ms).ExactScore, Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b, ms).ExactScore));
        }
    }
}
