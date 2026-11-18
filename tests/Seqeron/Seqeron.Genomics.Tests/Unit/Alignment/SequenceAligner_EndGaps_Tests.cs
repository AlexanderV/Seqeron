using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-SEMI-001 (audit A1-3, B09 F44): ends-free (generalised semi-global) alignment — each of the four end gaps
/// (leading/trailing residues of sequence 1/2) free or penalised, <see cref="EndGaps"/> (all 16 policies: global, fitting,
/// overlap in both orientations, all-free = EMBOSS needle's default), linear + affine, scalar + substitution-matrix scoring.
///
/// Oracle (TestData/Alignment/endgaps_oracle.tsv, scripts/alignment/generate_endgaps_oracle.py): Biopython 1.88
/// <c>PairwiseAligner</c> (<c>left/right_insertion/deletion_score = 0</c>; <c>align().score</c> and the enumerated co-optimal
/// alignments), parasail 1.3.4 <c>nw</c>/<c>sg_*</c>, EMBOSS needle 6.6.0 default (<c>-endweight N</c>).
/// </summary>
[TestFixture]
public class SequenceAligner_EndGaps_Tests
{
    private sealed record Row(
        string Kind, string Model, int Flags, string Matrix, int Match, int Mismatch, double GapOpen, double GapExtend,
        string Seq1, string Seq2, double Biopython, double BiopythonScoreApi, double? Parasail, double? Emboss, string Alignments);

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.endgaps_oracle.tsv";
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
            double? Opt(string s) => s == "NA" ? null : D(s);
            rows.Add(new Row(f[0], f[1], int.Parse(f[2], CultureInfo.InvariantCulture), f[3],
                int.Parse(f[4], CultureInfo.InvariantCulture), int.Parse(f[5], CultureInfo.InvariantCulture),
                D(f[6]), D(f[7]), f[8], f[9], D(f[10]), D(f[11]), Opt(f[12]), Opt(f[13]), f[14]));
        }
        return rows;
    }

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static SubstitutionMatrix MatrixOf(string name) => name switch
    {
        "BLOSUM62" => SubstitutionMatrix.Blosum62,
        "EDNAFULL" => SubstitutionMatrix.EdnaFull,
        _ => throw new ArgumentException(name),
    };

    private static AlignmentResult Align(Row r)
    {
        var ends = (EndGaps)r.Flags;
        bool affine = r.Model == "affine";
        if (r.Kind == "scalar")
        {
            var sc = new ScoringMatrix(r.Match, r.Mismatch, affine ? (int)r.GapOpen : 0, (int)r.GapExtend);
            return affine
                ? SequenceAligner.SemiGlobalAlignAffine(r.Seq1, r.Seq2, sc, ends)
                : SequenceAligner.SemiGlobalAlign(r.Seq1, r.Seq2, sc, ends);
        }
        var ms = new MatrixScoring(MatrixOf(r.Matrix), r.GapOpen, r.GapExtend);
        return affine
            ? SequenceAligner.SemiGlobalAlignAffine(r.Seq1, r.Seq2, ms, ends)
            : SequenceAligner.SemiGlobalAlign(r.Seq1, r.Seq2, ms, ends);
    }

    private static string Describe(Row r) =>
        $"{r.Kind}/{r.Model}/flags={r.Flags}/{r.Matrix} sc=({r.Match},{r.Mismatch},{r.GapOpen},{r.GapExtend}) {r.Seq1} {r.Seq2}";

    /// <summary>
    /// Independent re-scoring of aligned rows under an end-gap policy. A maximal run of columns that are gaps in
    /// sequence 2 (residue of sequence 1 against '-') is free when it starts the alignment and
    /// <see cref="EndGaps.FreeStart1"/> is set, or ends it and <see cref="EndGaps.FreeEnd1"/> is set; the same for runs of
    /// gaps in sequence 1 with the sequence-2 flags; every other gap run costs GapOpen + k·GapExtend (affine) or
    /// k·GapExtend (linear); aligned columns score by <paramref name="sub"/>.
    /// </summary>
    private static double Rescore(string r1, string r2, EndGaps ends, Func<char, char, double> sub, double open, double extend, bool affine)
    {
        double total = 0;
        int n = r1.Length;
        int c = 0;
        while (c < n)
        {
            if (r1[c] != '-' && r2[c] != '-')
            {
                total += sub(r1[c], r2[c]);
                c++;
                continue;
            }
            bool gapInSeq2 = r2[c] == '-'; // residue of sequence 1 against a gap
            int start = c;
            while (c < n && (gapInSeq2 ? r2[c] == '-' && r1[c] != '-' : r1[c] == '-' && r2[c] != '-')) c++;
            bool atStart = start == 0, atEnd = c == n;
            bool free = gapInSeq2
                ? (atStart && ends.HasFlag(EndGaps.FreeStart1)) || (atEnd && ends.HasFlag(EndGaps.FreeEnd1))
                : (atStart && ends.HasFlag(EndGaps.FreeStart2)) || (atEnd && ends.HasFlag(EndGaps.FreeEnd2));
            if (!free)
                total += (c - start) * extend + (affine ? open : 0);
        }
        return total;
    }

    private static Func<char, char, double> ScalarSub(int match, int mismatch) => (a, b) => a == b ? match : mismatch;

    private static Func<char, char, double> MatrixSub(SubstitutionMatrix m) => (a, b) => m.GetScore(a, b);

    /// <summary>Within each maximal run of gap-only columns the order of "residue of seq 1 vs gap" and "gap vs residue of
    /// seq 2" columns carries no information (Biopython lists one canonical order); normalise to all seq-1 columns first.</summary>
    private static string Canonical(string r1, string r2)
    {
        var o1 = new System.Text.StringBuilder();
        var o2 = new System.Text.StringBuilder();
        int c = 0;
        while (c < r1.Length)
        {
            if (r1[c] != '-' && r2[c] != '-')
            {
                o1.Append(r1[c]); o2.Append(r2[c]); c++;
                continue;
            }
            int s = c;
            while (c < r1.Length && (r1[c] == '-' || r2[c] == '-') && !(r1[c] == '-' && r2[c] == '-')) c++;
            for (int k = s; k < c; k++) if (r2[k] == '-') { o1.Append(r1[k]); o2.Append('-'); }
            for (int k = s; k < c; k++) if (r1[k] == '-') { o1.Append('-'); o2.Append(r2[k]); }
        }
        return o1 + "|" + o2;
    }

    // --- oracle differential ------------------------------------------------------------------------------

    [Test]
    public void Oracle_Covers_All16Policies_LinearAndAffine_ScalarAndMatrix()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(3000));
            Assert.That(rows.Select(r => r.Flags).Distinct().Count(), Is.EqualTo(16));
            Assert.That(rows.Select(r => r.Kind + r.Model).Distinct().Count(), Is.EqualTo(4));
            Assert.That(rows.Count(r => r.Parasail is not null), Is.GreaterThan(1500));
            Assert.That(rows.Count(r => r.Alignments != "NA"), Is.GreaterThan(1500));
            Assert.That(rows.Count(r => r.Emboss is not null), Is.GreaterThan(250));
        });
    }

    /// <summary>Every row: Seqeron score = Biopython 1.88 optimum (align().score), within 1e-6.</summary>
    [Test]
    public void Oracle_Score_EqualsBiopython_AllRows()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var a = Align(r);
            if (Math.Abs(a.ExactScore - r.Biopython) > 1e-6)
                failures.Add($"{Describe(r)}: {a.ExactScore} vs Biopython {r.Biopython}");
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    /// <summary>
    /// parasail 1.3.4: equal on every row of the policies whose parasail variant is exact (nw, sg_qb, sg_db, sg_qx, sg_dx,
    /// sg_qb_db). The variants with a free END flag (sg, sg_qe, sg_de, sg_qb_de, sg_qe_db, sg_qe_de) never exceed the optimum but
    /// can be lower: parasail does not consider the border end cells (m,0)/(0,n) (the alignment that ends with one sequence
    /// entirely against gaps) — verified against an exhaustive reference where "skip those two end cells" reproduces parasail on
    /// 300/300 random rows.
    /// </summary>
    [Test]
    public void Oracle_Score_EqualsParasail_ExactVariants_AndNeverAboveOptimumForTheOthers()
    {
        int[] exact = { 0, 1, 4, 3, 12, 5 };
        int exactRows = 0, deviating = 0;
        var failures = new List<string>();
        foreach (var r in Oracle.Value.Where(r => r.Parasail is not null))
        {
            double ours = Align(r).ExactScore;
            double p = r.Parasail!.Value;
            if (exact.Contains(r.Flags))
            {
                exactRows++;
                if (Math.Abs(ours - p) > 1e-9) failures.Add($"{Describe(r)}: {ours} vs parasail {p}");
            }
            else if (Math.Abs(ours - p) > 1e-9)
            {
                deviating++;
                if (p > ours) failures.Add($"{Describe(r)}: parasail {p} above optimum {ours}");
            }
        }
        Assert.Multiple(() =>
        {
            Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
            Assert.That(exactRows, Is.GreaterThan(800));
            Assert.That(deviating, Is.LessThan(60), "parasail end-cell quirk stays rare");
        });
    }

    /// <summary>EMBOSS needle 6.6.0 with its DEFAULT <c>-endweight N</c> (end gaps free on both sequences) = EndGaps.All: equal on
    /// almost every row; needle's last-row/column end gap is opened only from the match matrix (embaln.c), so the few disagreements
    /// are lower than the optimum, never higher.</summary>
    [Test]
    public void Oracle_Emboss_NeedleDefaultEndweight_IsEndGapsAll()
    {
        int rows = 0, equal = 0, higher = 0;
        foreach (var r in Oracle.Value.Where(r => r.Emboss is not null))
        {
            Assert.That(r.Flags, Is.EqualTo(15), "EMBOSS rows are the all-free policy");
            rows++;
            double ours = Align(r).ExactScore;
            if (Math.Abs(ours - r.Emboss!.Value) <= 0.051) equal++;
            else if (r.Emboss.Value > ours) higher++;
        }
        Assert.Multiple(() =>
        {
            Assert.That(rows, Is.GreaterThan(250));
            Assert.That((double)equal / rows, Is.GreaterThan(0.95), $"{equal}/{rows}");
            Assert.That(higher, Is.EqualTo(0), "needle never beats the optimum");
        });
    }

    /// <summary>The returned alignment is one of Biopython's co-optimal alignments (up to the order of adjacent gap columns of
    /// opposite type) whenever Biopython enumerates at most 40 of them.</summary>
    [Test]
    public void Oracle_Alignment_IsOneOfBiopythonCoOptimalAlignments()
    {
        int checkedRows = 0;
        var failures = new List<string>();
        foreach (var r in Oracle.Value.Where(r => r.Alignments != "NA"))
        {
            var a = Align(r);
            var set = r.Alignments.Split(';').Select(x => { var p = x.Split('|'); return Canonical(p[0], p[1]); }).ToHashSet();
            checkedRows++;
            if (!set.Contains(Canonical(a.AlignedSequence1, a.AlignedSequence2)))
                failures.Add($"{Describe(r)}: {a.AlignedSequence1}|{a.AlignedSequence2} not in {r.Alignments}");
        }
        Assert.That(checkedRows, Is.GreaterThan(1500));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    /// <summary>The alignment reproduces both inputs, has the reported score under independent re-scoring with the policy, and its
    /// coordinates are the consumed intervals (free flanks excluded).</summary>
    [Test]
    public void Oracle_Alignment_RescoresToScore_AndCoordinatesAreTheConsumedIntervals()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var a = Align(r);
            var ends = (EndGaps)r.Flags;
            bool affine = r.Model == "affine";
            Func<char, char, double> sub = r.Kind == "scalar" ? ScalarSub(r.Match, r.Mismatch) : MatrixSub(MatrixOf(r.Matrix));
            string tag = Describe(r);
            string r1 = a.AlignedSequence1, r2 = a.AlignedSequence2;
            if (r1.Length != r2.Length || r1.Replace("-", "") != r.Seq1 || r2.Replace("-", "") != r.Seq2)
            {
                failures.Add(tag + ": rows do not reproduce the inputs");
                continue;
            }
            if (Math.Abs(Rescore(r1, r2, ends, sub, r.GapOpen, r.GapExtend, affine) - a.ExactScore) > 1e-6)
                failures.Add(tag + $": rescore {Rescore(r1, r2, ends, sub, r.GapOpen, r.GapExtend, affine)} vs {a.ExactScore}");

            // coordinates: the consumed interval of each sequence = residues outside the leading/trailing free flanks
            if (a.AlignmentType != (ends == EndGaps.None ? AlignmentType.Global : AlignmentType.SemiGlobal))
                failures.Add(tag + ": type");
            int lead1 = a.StartPosition1, lead2 = a.StartPosition2;
            int trail1 = r.Seq1.Length - 1 - a.EndPosition1, trail2 = r.Seq2.Length - 1 - a.EndPosition2;
            if (lead1 < 0 || lead2 < 0 || trail1 < 0 || trail2 < 0 || (lead1 > 0 && lead2 > 0) || (trail1 > 0 && trail2 > 0))
            {
                failures.Add(tag + $": coordinates {lead1},{lead2},{trail1},{trail2}");
                continue;
            }
            if ((lead1 > 0 && !ends.HasFlag(EndGaps.FreeStart1)) || (lead2 > 0 && !ends.HasFlag(EndGaps.FreeStart2))
                || (trail1 > 0 && !ends.HasFlag(EndGaps.FreeEnd1)) || (trail2 > 0 && !ends.HasFlag(EndGaps.FreeEnd2)))
                failures.Add(tag + ": overhang on a penalised end");
            // flank columns: first lead1 columns are residues of seq1 against gaps (or lead2 of seq2), same at the end
            int lead = lead1 + lead2, trail = trail1 + trail2;
            for (int c = 0; c < lead; c++)
                if ((lead1 > 0 ? r2[c] != '-' || r1[c] == '-' : r1[c] != '-' || r2[c] == '-')) { failures.Add(tag + ": leading flank columns"); break; }
            for (int c = r1.Length - trail; c < r1.Length; c++)
                if ((trail1 > 0 ? r2[c] != '-' || r1[c] == '-' : r1[c] != '-' || r2[c] == '-')) { failures.Add(tag + ": trailing flank columns"); break; }
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    // --- identities with the existing aligners ---------------------------------------------------------------

    private static IEnumerable<(string A, string B)> RandomPairs(int count, int seed, int maxLen = 30)
    {
        var rnd = new Random(seed);
        string Rand(int len, string alpha) => new(Enumerable.Range(0, len).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray());
        for (int i = 0; i < count; i++)
        {
            string alpha = i % 3 == 0 ? "AC" : "ACGT";
            string a = Rand(rnd.Next(0, maxLen), alpha);
            string b = rnd.Next(3) == 0 ? Rand(rnd.Next(0, maxLen), alpha) : MutateCopy(a, rnd, alpha) + Rand(rnd.Next(0, 6), alpha);
            yield return (a, b);
        }
    }

    private static string MutateCopy(string s, Random rnd, string alpha)
    {
        var sb = new System.Text.StringBuilder();
        foreach (char ch in s)
        {
            double x = rnd.NextDouble();
            if (x < 0.1) sb.Append(alpha[rnd.Next(alpha.Length)]);
            else if (x < 0.15) { }
            else if (x < 0.2) { sb.Append(ch); sb.Append(alpha[rnd.Next(alpha.Length)]); }
            else sb.Append(ch);
        }
        return sb.ToString();
    }

    private static readonly ScoringMatrix[] Schemes =
    {
        new(2, -3, 0, -2), new(2, -3, -3, -1), new(1, -1, -2, -1), new(5, -4, -10, -1), new(3, -2, 0, -4), new(2, -3, -2, 0),
    };

    private static void AssertSame(AlignmentResult x, AlignmentResult y, string tag)
    {
        Assert.Multiple(() =>
        {
            Assert.That(x.AlignedSequence1, Is.EqualTo(y.AlignedSequence1), tag);
            Assert.That(x.AlignedSequence2, Is.EqualTo(y.AlignedSequence2), tag);
            Assert.That(x.Score, Is.EqualTo(y.Score), tag);
            Assert.That(x.AlignmentType, Is.EqualTo(y.AlignmentType), tag);
            Assert.That((x.StartPosition1, x.StartPosition2, x.EndPosition1, x.EndPosition2),
                Is.EqualTo((y.StartPosition1, y.StartPosition2, y.EndPosition1, y.EndPosition2)), tag);
        });
    }

    /// <summary>EndGaps.Fitting is exactly the classic SemiGlobalAlign (rows, score, coordinates, tie order) for the linear and affine
    /// aligners, string and DnaSequence overloads; EndGaps.None is exactly GlobalAlign / GlobalAlignAffine (tagged Global).</summary>
    [Test]
    public void Fitting_IsTheClassicSemiGlobal_AndNone_IsGlobal_RowForRow()
    {
        foreach (var (a, b) in RandomPairs(300, 11))
            foreach (var sc in Schemes)
            {
                string tag = $"{a} {b} {sc}";
                AssertSame(SequenceAligner.SemiGlobalAlign(a, b, sc, EndGaps.Fitting), SequenceAligner.SemiGlobalAlign(a, b, sc), "linear fitting " + tag);
                AssertSame(SequenceAligner.SemiGlobalAlignAffine(a, b, sc, EndGaps.Fitting), SequenceAligner.SemiGlobalAlignAffine(a, b, sc), "affine fitting " + tag);
                AssertSame(SequenceAligner.SemiGlobalAlign(a, b, sc, EndGaps.None), SequenceAligner.GlobalAlign(a, b, sc), "linear none " + tag);
                AssertSame(SequenceAligner.SemiGlobalAlignAffine(a, b, sc, EndGaps.None), SequenceAligner.GlobalAlignAffine(a, b, sc), "affine none " + tag);
                if (a.Length > 0 && b.Length > 0)
                {
                    var da = new DnaSequence(a);
                    var db = new DnaSequence(b);
                    foreach (var ends in new[] { EndGaps.Fitting, EndGaps.None, EndGaps.All, EndGaps.OverlapSuffix1Prefix2 })
                    {
                        AssertSame(SequenceAligner.SemiGlobalAlign(da, db, sc, ends), SequenceAligner.SemiGlobalAlign(a, b, sc, ends), "dna linear " + tag);
                        AssertSame(SequenceAligner.SemiGlobalAlignAffine(da, db, sc, ends), SequenceAligner.SemiGlobalAlignAffine(a, b, sc, ends), "dna affine " + tag);
                    }
                }
            }
    }

    [Test]
    public void MatrixScoring_FittingAndNone_EqualTheClassicMatrixAligners()
    {
        var ms = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62);
        var rnd = new Random(5);
        const string alpha = "ACDEFGHIKLMNPQRSTVWY";
        for (int i = 0; i < 200; i++)
        {
            string a = new(Enumerable.Range(0, rnd.Next(1, 40)).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray());
            string b = MutateCopy(a, rnd, alpha);
            AssertSame(SequenceAligner.SemiGlobalAlign(a, b, ms, EndGaps.Fitting), SequenceAligner.SemiGlobalAlign(a, b, ms), "linear " + a + b);
            AssertSame(SequenceAligner.SemiGlobalAlignAffine(a, b, ms, EndGaps.Fitting), SequenceAligner.SemiGlobalAlignAffine(a, b, ms), "affine " + a + b);
            AssertSame(SequenceAligner.SemiGlobalAlign(a, b, ms, EndGaps.None), SequenceAligner.GlobalAlign(a, b, ms), "linear none " + a + b);
            AssertSame(SequenceAligner.SemiGlobalAlignAffine(new ProteinSequence(a), new ProteinSequence(b), ms, EndGaps.None),
                SequenceAligner.GlobalAlignAffine(a, b, ms), "affine none protein " + a + b);
            Assert.That(SequenceAligner.SemiGlobalAlign(new ProteinSequence(a), new ProteinSequence(b), ms, EndGaps.All).ExactScore,
                Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b, ms, EndGaps.All).ExactScore));
        }
    }

    // --- semantics ----------------------------------------------------------------------------------------

    /// <summary>Overlap (dovetail) alignment, values from Biopython 1.88 (left_deletion_score = right_insertion_score = 0,
    /// match 1, mismatch -1, linear gap -1; alignment unique): a suffix of sequence 1 against a prefix of sequence 2.</summary>
    [Test]
    public void Overlap_SuffixOfSeq1_PrefixOfSeq2_RowsAndCoordinates()
    {
        var sc = new ScoringMatrix(1, -1, 0, -1);
        var r = SequenceAligner.SemiGlobalAlign("GGGGACGTAC", "ACGTACTTTT", sc, EndGaps.OverlapSuffix1Prefix2);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequence1, Is.EqualTo("GGGGACGTAC----"));
            Assert.That(r.AlignedSequence2, Is.EqualTo("----ACGTACTTTT"));
            Assert.That(r.Score, Is.EqualTo(6));
            Assert.That(r.AlignmentType, Is.EqualTo(AlignmentType.SemiGlobal));
            Assert.That((r.StartPosition1, r.EndPosition1, r.StartPosition2, r.EndPosition2), Is.EqualTo((4, 9, 0, 5)));
        });

        // the other orientation finds the same overlap when the sequences are swapped
        var s = SequenceAligner.SemiGlobalAlign("ACGTACTTTT", "GGGGACGTAC", sc, EndGaps.OverlapSuffix2Prefix1);
        Assert.Multiple(() =>
        {
            Assert.That(s.AlignedSequence1, Is.EqualTo("----ACGTACTTTT"));
            Assert.That(s.AlignedSequence2, Is.EqualTo("GGGGACGTAC----"));
            Assert.That(s.Score, Is.EqualTo(6));
            Assert.That((s.StartPosition1, s.EndPosition1, s.StartPosition2, s.EndPosition2), Is.EqualTo((0, 5, 4, 9)));
        });
    }

    /// <summary>Containment: sequence 2 fitted inside sequence 1 is the mirror image of fitting (FittingSequence2 = parasail sg_qx).</summary>
    [Test]
    public void FittingSequence2_IsFittingOfTheSwappedPair()
    {
        foreach (var (a, b) in RandomPairs(200, 21))
            foreach (var sc in Schemes)
            {
                Assert.That(SequenceAligner.SemiGlobalAlign(a, b, sc, EndGaps.FittingSequence2).Score,
                    Is.EqualTo(SequenceAligner.SemiGlobalAlign(b, a, sc, EndGaps.Fitting).Score), $"{a} {b}");
                Assert.That(SequenceAligner.SemiGlobalAlignAffine(a, b, sc, EndGaps.OverlapSuffix2Prefix1).Score,
                    Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine(b, a, sc, EndGaps.OverlapSuffix1Prefix2).Score), $"{a} {b}");
            }
    }

    /// <summary>Freeing an end gap can only raise the score (monotone in the flag set), and All is the maximum.</summary>
    [Test]
    public void Score_IsMonotoneInTheFreeFlags()
    {
        foreach (var (a, b) in RandomPairs(120, 31, 18))
            foreach (var sc in Schemes.Take(3))
                for (int f = 0; f < 16; f++)
                    for (int bit = 1; bit < 16; bit <<= 1)
                    {
                        if ((f & bit) != 0) continue;
                        Assert.That(SequenceAligner.SemiGlobalAlignAffine(a, b, sc, (EndGaps)(f | bit)).Score,
                            Is.GreaterThanOrEqualTo(SequenceAligner.SemiGlobalAlignAffine(a, b, sc, (EndGaps)f).Score), $"{a} {b} {f}+{bit}");
                    }
    }

    /// <summary>Independent exhaustive reference: enumerate EVERY alignment of short sequences (all column sequences), score each
    /// with the rescorer of the policy, take the maximum — equals the dynamic program for all 16 policies, linear and affine.</summary>
    [Test]
    public void BruteForce_EveryAlignment_MaximumEqualsDynamicProgram_All16Policies()
    {
        var rnd = new Random(41);
        int pairs = 0;
        for (int t = 0; t < 120; t++)
        {
            string Rand(int len) => new(Enumerable.Range(0, len).Select(_ => "ACG"[rnd.Next(3)]).ToArray());
            string a = Rand(rnd.Next(0, 6)), b = Rand(rnd.Next(0, 6));
            var sc = Schemes[t % Schemes.Length];
            var all = new List<(string, string)>();
            Enumerate(a, b, 0, 0, "", "", all);
            pairs++;
            for (int f = 0; f < 16; f++)
            {
                var ends = (EndGaps)f;
                foreach (bool affine in new[] { false, true })
                {
                    double best = all.Max(x => Rescore(x.Item1, x.Item2, ends, ScalarSub(sc.Match, sc.Mismatch), sc.GapOpen, sc.GapExtend, affine));
                    var r = affine ? SequenceAligner.SemiGlobalAlignAffine(a, b, sc, ends) : SequenceAligner.SemiGlobalAlign(a, b, sc, ends);
                    Assert.That(r.Score, Is.EqualTo(best), $"{a} {b} flags={f} affine={affine} {sc}");
                }
            }
        }
        Assert.That(pairs, Is.EqualTo(120));
    }

    private static void Enumerate(string a, string b, int i, int j, string r1, string r2, List<(string, string)> all)
    {
        if (i == a.Length && j == b.Length) { all.Add((r1, r2)); return; }
        if (i < a.Length && j < b.Length) Enumerate(a, b, i + 1, j + 1, r1 + a[i], r2 + b[j], all);
        if (i < a.Length) Enumerate(a, b, i + 1, j, r1 + a[i], r2 + "-", all);
        if (j < b.Length) Enumerate(a, b, i, j + 1, r1 + "-", r2 + b[j], all);
    }

    [Test]
    public void EmptyAndNull_Inputs_FollowTheAlignerContract()
    {
        var sc = new ScoringMatrix(2, -3, -2, -1);
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.SemiGlobalAlign((string)null!, "A", sc, EndGaps.All));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.SemiGlobalAlignAffine("A", (string)null!, sc, EndGaps.All));
        // everything free: an empty side costs nothing; penalised end gaps charge the all-gap alignment
        var free = SequenceAligner.SemiGlobalAlign("", "ACGT", sc, EndGaps.All);
        Assert.Multiple(() =>
        {
            Assert.That((free.AlignedSequence1, free.AlignedSequence2, free.Score), Is.EqualTo(("----", "ACGT", 0)));
            Assert.That((free.StartPosition1, free.EndPosition1, free.StartPosition2, free.EndPosition2), Is.EqualTo((0, -1, 0, -1)));
        });
        var penalised = SequenceAligner.SemiGlobalAlign("ACGT", "", sc, EndGaps.FreeStart2 | EndGaps.FreeEnd2);
        Assert.That(penalised.Score, Is.EqualTo(-4));
        Assert.That(SequenceAligner.SemiGlobalAlign("ACGT", "", sc, EndGaps.FreeEnd1).Score, Is.EqualTo(0));
        Assert.That(SequenceAligner.SemiGlobalAlignAffine("", "", sc, EndGaps.All).Score, Is.EqualTo(0));
    }

    /// <summary>The documented degenerate case on which Biopython's PairwiseAligner.score() disagrees with its own align()
    /// (score('AC','T') = -2 with left_deletion_score = right_insertion_score = 0, the alignment AC-/--T it returns scores 0): the
    /// optimum, Seqeron's value and Biopython align().score are 0 (zero-length overlap).</summary>
    [Test]
    public void ZeroLengthOverlap_ScoresZero_BiopythonAlignScore()
    {
        var sc = new ScoringMatrix(2, -3, 0, -1);
        var r = SequenceAligner.SemiGlobalAlign("AC", "T", sc, EndGaps.OverlapSuffix1Prefix2);
        Assert.Multiple(() =>
        {
            Assert.That(r.Score, Is.EqualTo(0));
            Assert.That(r.AlignedSequence1, Is.EqualTo("AC-"));
            Assert.That(r.AlignedSequence2, Is.EqualTo("--T"));
            Assert.That((r.StartPosition1, r.EndPosition1, r.StartPosition2, r.EndPosition2), Is.EqualTo((2, 1, 0, -1)));
        });
        Assert.That(Oracle.Value.Count(x => Math.Abs(x.Biopython - x.BiopythonScoreApi) > 1e-9), Is.GreaterThan(0),
            "the Biopython score() quirk is present in the oracle data");
    }
}
