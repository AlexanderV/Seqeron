using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 (audit A1-5, B09 F45): alignment with an arbitrary gap-cost function w(k)
/// (Waterman-Smith-Beyer 1976) — <c>GlobalAlignGeneralGap</c>, <c>LocalAlignGeneralGap</c>, <c>SemiGlobalAlignGeneralGap</c>.
///
/// Oracle (TestData/Alignment/general_gap_oracle.tsv, scripts/alignment/generate_general_gap_oracle.py): 960 random global/local pairs
/// with linear, affine, logarithmic, square-root (concave), quadratic (convex), capped, step and two-slope gap functions, scalar and
/// BLOSUM62/EDNAFULL scoring — Biopython 1.88 <c>PairwiseAligner(gap_score=callable)</c> (Waterman-Smith-Beyer algorithm) = an independent
/// three-state reference on every row; 480 ends-free rows with the independent reference only (Biopython ignores the end-gap
/// properties for a callable gap score).
/// </summary>
[TestFixture]
public class SequenceAligner_GeneralGap_Tests
{
    private sealed record Row(
        string Mode, int Flags, string Kind, double A, double B, double C, double D, string Scoring, int Match, int Mismatch,
        string Seq1, string Seq2, double? Biopython, double PyRef);

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static double Dbl(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.general_gap_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], int.Parse(f[1], CultureInfo.InvariantCulture), f[2], Dbl(f[3]), Dbl(f[4]), Dbl(f[5]), Dbl(f[6]), f[7],
                int.Parse(f[8], CultureInfo.InvariantCulture), int.Parse(f[9], CultureInfo.InvariantCulture), f[10], f[11],
                f[12] == "NA" ? null : Dbl(f[12]), Dbl(f[13])));
        }
        return rows;
    }

    private static Func<int, double> GapFunction(string kind, double a, double b, double c, double d) => kind switch
    {
        "linear" => k => -(a * k),
        "affine" => k => -(a + b * k),
        "log" => k => -(a + b * Math.Log(k)),
        "sqrt" => k => -(a + b * Math.Sqrt(k)),
        "quad" => k => -(a + b * k * k),
        "cap" => k => -Math.Min(a + b * k, c),
        "step" => k => -(k <= c ? a : b),
        "twoslope" => k => k <= c ? -(a + b * k) : -(a + b * c + d * (k - c)),
        _ => throw new ArgumentException(kind),
    };

    private static Func<int, double> GapOf(Row r) => GapFunction(r.Kind, r.A, r.B, r.C, r.D);

    private static SubstitutionMatrix? MatrixOf(string scoring) => scoring switch
    {
        "BLOSUM62" => SubstitutionMatrix.Blosum62,
        "EDNAFULL" => SubstitutionMatrix.EdnaFull,
        _ => null,
    };

    private static AlignmentType TypeOf(string mode) => mode switch
    {
        "global" => AlignmentType.Global,
        "local" => AlignmentType.Local,
        _ => AlignmentType.SemiGlobal,
    };

    /// <summary>Public API (fast path for linear/affine functions, general program otherwise).</summary>
    private static AlignmentResult Align(Row r)
    {
        var w = GapOf(r);
        var m = MatrixOf(r.Scoring);
        var sc = new ScoringMatrix(r.Match, r.Mismatch, 0, 0);
        var ends = (EndGaps)r.Flags;
        return r.Mode switch
        {
            "global" => m is null ? SequenceAligner.GlobalAlignGeneralGap(r.Seq1, r.Seq2, w, sc) : SequenceAligner.GlobalAlignGeneralGap(r.Seq1, r.Seq2, w, m),
            "local" => m is null ? SequenceAligner.LocalAlignGeneralGap(r.Seq1, r.Seq2, w, sc) : SequenceAligner.LocalAlignGeneralGap(r.Seq1, r.Seq2, w, m),
            _ => SequenceAligner.SemiGlobalAlignGeneralGap(r.Seq1, r.Seq2, w, sc, ends),
        };
    }

    /// <summary>The general program itself (never the linear/affine fast path).</summary>
    private static AlignmentResult AlignGeneric(Row r)
    {
        var w = GapOf(r);
        var m = MatrixOf(r.Scoring);
        var type = TypeOf(r.Mode);
        var ends = (EndGaps)r.Flags;
        return m is null
            ? SequenceAligner.AlignGeneralGapGeneric(r.Seq1, r.Seq2, w, new ScoringMatrix(r.Match, r.Mismatch, 0, 0), type, ends)
            : SequenceAligner.AlignGeneralGapGeneric(r.Seq1, r.Seq2, w, m, type, ends);
    }

    private static string Describe(Row r) => $"{r.Mode}/{r.Flags}/{r.Kind}({r.A},{r.B},{r.C},{r.D})/{r.Scoring} {r.Seq1} {r.Seq2}";

    /// <summary>Independent re-scoring: every maximal run of residues of one sequence against gaps costs w(run length) — except
    /// the free end runs of the policy; aligned columns score by <paramref name="sub"/>.</summary>
    private static double Rescore(string r1, string r2, EndGaps ends, Func<char, char, double> sub, Func<int, double> w)
    {
        double total = 0;
        int c = 0, n = r1.Length;
        while (c < n)
        {
            if (r1[c] != '-' && r2[c] != '-')
            {
                total += sub(r1[c], r2[c]);
                c++;
                continue;
            }
            bool gapInSeq2 = r2[c] == '-';
            int start = c;
            while (c < n && (gapInSeq2 ? r2[c] == '-' && r1[c] != '-' : r1[c] == '-' && r2[c] != '-')) c++;
            bool atStart = start == 0, atEnd = c == n;
            bool free = gapInSeq2
                ? (atStart && ends.HasFlag(EndGaps.FreeStart1)) || (atEnd && ends.HasFlag(EndGaps.FreeEnd1))
                : (atStart && ends.HasFlag(EndGaps.FreeStart2)) || (atEnd && ends.HasFlag(EndGaps.FreeEnd2));
            if (!free) total += w(c - start);
        }
        return total;
    }

    private static Func<char, char, double> SubOf(Row r)
    {
        var m = MatrixOf(r.Scoring);
        return m is null ? (a, b) => a == b ? r.Match : r.Mismatch : (a, b) => m.GetScore(a, b);
    }

    // --- oracle --------------------------------------------------------------------------------------------

    [Test]
    public void Oracle_Covers_AllGapFamilies_GlobalLocalSemi_AndScalarAndMatrixScoring()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows.Count(r => r.Mode != "semi"), Is.GreaterThanOrEqualTo(900));
            Assert.That(rows.Select(r => r.Kind).Distinct().Count(), Is.EqualTo(8));
            Assert.That(rows.Select(r => r.Mode).Distinct().Count(), Is.EqualTo(3));
            Assert.That(rows.Select(r => r.Scoring).Distinct().Order().ToArray(), Is.EqualTo(new[] { "BLOSUM62", "EDNAFULL", "scalar" }));
            Assert.That(rows.Count(r => r.Mode == "semi"), Is.GreaterThanOrEqualTo(400));
            Assert.That(rows.Where(r => r.Biopython is not null).All(r => Math.Abs(r.Biopython!.Value - r.PyRef) < 1e-9), Is.True,
                "Biopython 1.88 callable gap_score = the independent reference on every row");
        });
    }

    /// <summary>Every row: the general program = Biopython 1.88 (callable gap_score) where available and = the independent reference on
    /// every row, global, local and ends-free; the public API (with its linear/affine fast path) returns the same score.</summary>
    [Test]
    public void Oracle_Score_GeneralProgramAndPublicApi_EqualBiopythonAndReference()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            double expected = r.Biopython ?? r.PyRef;
            double generic = AlignGeneric(r).ExactScore;
            double api = Align(r).ExactScore;
            double tol = 1e-8 * Math.Max(1, Math.Abs(expected));
            if (Math.Abs(generic - expected) > tol) failures.Add($"{Describe(r)}: general {generic} vs {expected}");
            if (Math.Abs(api - expected) > tol) failures.Add($"{Describe(r)}: API {api} vs {expected}");
            if (Math.Abs(generic - r.PyRef) > tol) failures.Add($"{Describe(r)}: general {generic} vs pyref {r.PyRef}");
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    /// <summary>The returned alignment reproduces the inputs, re-scores (gap runs charged w(k), free ends excluded) to its score, and its
    /// coordinates are consistent.</summary>
    [Test]
    public void Oracle_Alignment_RescoresToScore_AndCoordinatesAreConsistent()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            foreach (var res in new[] { Align(r), AlignGeneric(r) })
            {
                string tag = Describe(r);
                string r1 = res.AlignedSequence1, r2 = res.AlignedSequence2;
                var w = GapOf(r);
                var sub = SubOf(r);
                var ends = (EndGaps)r.Flags;
                if (r.Mode == "local")
                {
                    if (res.ExactScore == 0)
                    {
                        if (r1.Length != 0) failures.Add(tag + ": non-empty zero local");
                        continue;
                    }
                    if (r.Seq1.Substring(res.StartPosition1, res.EndPosition1 - res.StartPosition1 + 1) != r1.Replace("-", "")
                        || r.Seq2.Substring(res.StartPosition2, res.EndPosition2 - res.StartPosition2 + 1) != r2.Replace("-", ""))
                        failures.Add(tag + ": local intervals");
                    if (r1.Length > 0 && (r1[0] == '-' || r2[0] == '-' || r1[^1] == '-' || r2[^1] == '-')) failures.Add(tag + ": local alignment starts/ends with a gap");
                }
                else
                {
                    if (r1.Replace("-", "") != r.Seq1 || r2.Replace("-", "") != r.Seq2) { failures.Add(tag + ": rows do not reproduce the inputs"); continue; }
                    if (r.Mode == "global" && (res.StartPosition1, res.StartPosition2, res.EndPosition1, res.EndPosition2) != (0, 0, r.Seq1.Length - 1, r.Seq2.Length - 1))
                        failures.Add(tag + ": global coordinates");
                    if (res.AlignmentType != (r.Mode == "global" ? AlignmentType.Global : AlignmentType.SemiGlobal)) failures.Add(tag + ": type");
                }
                double rescored = Rescore(r1, r2, r.Mode == "semi" ? ends : EndGaps.None, sub, w);
                if (Math.Abs(rescored - res.ExactScore) > 1e-8 * Math.Max(1, Math.Abs(rescored))) failures.Add($"{tag}: rescore {rescored} vs {res.ExactScore}");
            }
        }
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    // --- fast path: linear / affine functions --------------------------------------------------------------

    private static IEnumerable<(string A, string B)> Pairs(int count, int seed, int maxLen)
    {
        var rnd = new Random(seed);
        string Rand(int len) => new(Enumerable.Range(0, len).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
        for (int i = 0; i < count; i++)
        {
            string a = Rand(rnd.Next(0, maxLen));
            var sb = new System.Text.StringBuilder();
            foreach (char ch in a)
            {
                double x = rnd.NextDouble();
                if (x < 0.1) sb.Append("ACGT"[rnd.Next(4)]);
                else if (x < 0.2) { }
                else sb.Append(ch);
            }
            yield return (a, sb + Rand(rnd.Next(0, 5)));
        }
    }

    /// <summary>A linear or affine gap function is aligned by the O(m·n) dynamic programs: rows identical to the affine / linear aligners
    /// (so the fast path is what ran), and the general program gives the same score.</summary>
    [Test]
    public void AffineAndLinearFunctions_UseTheAffineAndLinearAligners_RowForRow()
    {
        foreach (var (a, b) in Pairs(150, 3, 28))
            foreach (var (open, ext) in new[] { (-3, -1), (-5, -2), (0, -2), (-10, -1) })
            {
                var sc = new ScoringMatrix(2, -3, open, ext);
                Func<int, double> w = k => open + ext * k;
                string tag = $"{a} {b} {open} {ext}";
                AssertSame(SequenceAligner.GlobalAlignGeneralGap(a, b, w, sc), SequenceAligner.GlobalAlignAffine(a, b, sc), "global " + tag);
                AssertSame(SequenceAligner.LocalAlignGeneralGap(a, b, w, sc), SequenceAligner.LocalAlignAffine(a, b, sc), "local " + tag);
                AssertSame(SequenceAligner.SemiGlobalAlignGeneralGap(a, b, w, sc), SequenceAligner.SemiGlobalAlignAffine(a, b, sc), "semi " + tag);
                foreach (var ends in new[] { EndGaps.None, EndGaps.All, EndGaps.OverlapSuffix1Prefix2, EndGaps.FittingSequence2 })
                {
                    AssertSame(SequenceAligner.SemiGlobalAlignGeneralGap(a, b, w, sc, ends), SequenceAligner.SemiGlobalAlignAffine(a, b, sc, ends), $"ends {ends} " + tag);
                    var generic = SequenceAligner.AlignGeneralGapGeneric(a, b, w, sc, AlignmentType.SemiGlobal, ends);
                    Assert.That(generic.Score, Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine(a, b, sc, ends).Score), $"generic ends {ends} " + tag);
                }
                if (open == 0)
                    AssertSame(SequenceAligner.GlobalAlignGeneralGap(a, b, w, sc), SequenceAligner.GlobalAlign(a, b, sc), "linear " + tag);
            }
    }

    private static void AssertSame(AlignmentResult x, AlignmentResult y, string tag)
    {
        Assert.Multiple(() =>
        {
            Assert.That(x.AlignedSequence1, Is.EqualTo(y.AlignedSequence1), tag);
            Assert.That(x.AlignedSequence2, Is.EqualTo(y.AlignedSequence2), tag);
            Assert.That(x.Score, Is.EqualTo(y.Score), tag);
            Assert.That(x.AlignmentType, Is.EqualTo(y.AlignmentType), tag);
            Assert.That((x.StartPosition1, x.StartPosition2, x.EndPosition1, x.EndPosition2), Is.EqualTo((y.StartPosition1, y.StartPosition2, y.EndPosition1, y.EndPosition2)), tag);
        });
    }

    /// <summary>EMBOSS-style fractional affine function with a substitution matrix is delegated with exact (fractional) scores; the
    /// general program agrees.</summary>
    [Test]
    public void FractionalAffineFunction_WithBlosum62_EqualsTheMatrixAffineAligner_AndTheGeneralProgram()
    {
        var ms = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62); // gapopen 10, gapextend 0.5 => w(k) = -(9.5 + 0.5 k)
        Func<int, double> w = k => -(9.5 + 0.5 * k);
        var rnd = new Random(8);
        const string alpha = "ACDEFGHIKLMNPQRSTVWY";
        for (int i = 0; i < 60; i++)
        {
            string a = new(Enumerable.Range(0, rnd.Next(1, 30)).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray());
            string b = new(Enumerable.Range(0, rnd.Next(1, 30)).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray());
            var fast = SequenceAligner.GlobalAlignGeneralGap(a, b, w, SubstitutionMatrix.Blosum62);
            var affine = SequenceAligner.GlobalAlignAffine(a, b, ms);
            var generic = SequenceAligner.AlignGeneralGapGeneric(a, b, w, SubstitutionMatrix.Blosum62, AlignmentType.Global, EndGaps.None);
            Assert.Multiple(() =>
            {
                Assert.That(fast.ExactScore, Is.EqualTo(affine.ExactScore).Within(1e-9), a + b);
                Assert.That(fast.AlignedSequence1, Is.EqualTo(affine.AlignedSequence1));
                Assert.That(generic.ExactScore, Is.EqualTo(affine.ExactScore).Within(1e-9), a + b);
            });
        }
    }

    // --- semantics ----------------------------------------------------------------------------------------

    /// <summary>Exhaustive enumeration of every alignment of short sequences, scored by the independent rescorer: the maximum equals the
    /// general program for concave, convex, step and capped gap functions, every end-gap policy (and local).</summary>
    [Test]
    public void BruteForce_EveryAlignment_MaximumEqualsTheGeneralProgram()
    {
        var functions = new (string Name, Func<int, double> W)[]
        {
            ("log", k => -(1 + 2 * Math.Log(k))), ("quad", k => -(0.5 * k * k)), ("step", k => k <= 2 ? -1 : -6),
            ("cap", k => -Math.Min(1 + k, 4)), ("twoslope", k => k <= 2 ? -(1 + 0.5 * k) : -(2 + 2 * (k - 2))),
            ("bonus", k => k == 3 ? -0.5 : -2.0 * k), // non-monotone: a gap of length 3 is cheap
        };
        var rnd = new Random(77);
        int checkedPairs = 0;
        for (int t = 0; t < 60; t++)
        {
            string Rand(int len) => new(Enumerable.Range(0, len).Select(_ => "ACG"[rnd.Next(3)]).ToArray());
            string a = Rand(rnd.Next(0, t < 15 ? 5 : 6)), b = Rand(rnd.Next(0, t < 15 ? 5 : 6));
            var all = new List<(string, string)>();
            Enumerate(a, b, 0, 0, "", "", all);
            var sc = new ScoringMatrix(2, -3, 0, 0);
            Func<char, char, double> sub = (x, y) => x == y ? 2 : -3;
            foreach (var (name, w) in functions)
            {
                for (int f = 0; f < 16; f++)
                {
                    double best = all.Max(x => Rescore(x.Item1, x.Item2, (EndGaps)f, sub, w));
                    var r = SequenceAligner.SemiGlobalAlignGeneralGap(a, b, w, sc, (EndGaps)f);
                    Assert.That(r.ExactScore, Is.EqualTo(best).Within(1e-9), $"{name} {a} {b} flags={f}");
                }
                if (t >= 15) continue; // the substring-pair brute force is the expensive part
                double localBest = Math.Max(0, LocalBrute(a, b, w, sub));
                Assert.That(SequenceAligner.LocalAlignGeneralGap(a, b, w, sc).ExactScore, Is.EqualTo(localBest).Within(1e-9), $"local {name} {a} {b}");
            }
            checkedPairs++;
        }
        Assert.That(checkedPairs, Is.EqualTo(60));
    }

    /// <summary>Best local alignment by brute force: the best global alignment of any pair of substrings (a local alignment is a global
    /// alignment of the two substrings it covers).</summary>
    private static double LocalBrute(string a, string b, Func<int, double> w, Func<char, char, double> sub)
    {
        double best = 0;
        for (int i = 0; i <= a.Length; i++)
            for (int j = i + 1; j <= a.Length; j++)
                for (int p = 0; p <= b.Length; p++)
                    for (int q = p + 1; q <= b.Length; q++)
                    {
                        var all = new List<(string, string)>();
                        Enumerate(a[i..j], b[p..q], 0, 0, "", "", all);
                        // a local alignment neither starts nor ends with a gap column (it would only lower the score)
                        foreach (var (x, y) in all)
                            best = Math.Max(best, Rescore(x, y, EndGaps.None, sub, w));
                    }
        return best;
    }

    private static void Enumerate(string a, string b, int i, int j, string r1, string r2, List<(string, string)> all)
    {
        if (i == a.Length && j == b.Length) { all.Add((r1, r2)); return; }
        if (i < a.Length && j < b.Length) Enumerate(a, b, i + 1, j + 1, r1 + a[i], r2 + b[j], all);
        if (i < a.Length) Enumerate(a, b, i + 1, j, r1 + a[i], r2 + "-", all);
        if (j < b.Length) Enumerate(a, b, i, j + 1, r1 + "-", r2 + b[j], all);
    }

    /// <summary>A gap run is charged once as w(k), never split (Biopython 1.88 PairwiseAligner(gap_score=callable); values computed with it,
    /// match 2 / mismatch -3): with the convex w(k) = -k^2 the run of four C's would cost -16 (total 0), so the optimum (12 co-optimal
    /// alignments in Biopython) spreads over short runs and a mismatch: 5; with the concave w(k) = -(2 + 3 ln k) the single long run is
    /// optimal and unique: 16 - (2 + 3 ln 4) = 9.841116916640328.</summary>
    [Test]
    public void ConvexAndConcaveGaps_BiopythonValues()
    {
        var sc = new ScoringMatrix(2, -3, 0, 0);
        var convex = SequenceAligner.GlobalAlignGeneralGap("AAAATTTT", "AAAACCCCTTTT", k => -k * k, sc);
        Assert.That(convex.ExactScore, Is.EqualTo(5));
        var concave = SequenceAligner.GlobalAlignGeneralGap("AAAATTTT", "AAAACCCCTTTT", k => -(2 + 3 * Math.Log(k)), sc);
        Assert.Multiple(() =>
        {
            Assert.That(concave.ExactScore, Is.EqualTo(9.841116916640328).Within(1e-12));
            Assert.That(concave.AlignedSequence1, Is.EqualTo("AAAA----TTTT"));
            Assert.That(concave.AlignedSequence2, Is.EqualTo("AAAACCCCTTTT"));
        });
    }

    [Test]
    public void GapFunction_IsEvaluatedOncePerLength_AndMustBeValid()
    {
        var calls = new List<int>();
        var r = SequenceAligner.GlobalAlignGeneralGap("ACGTACGT", "ACGGT", k => { calls.Add(k); return -(1 + Math.Sqrt(k)); });
        Assert.That(calls, Is.EqualTo(Enumerable.Range(1, 8).ToList()), "w(k) evaluated exactly once for k = 1 .. max(m, n)");
        Assert.That(r.ExactScore, Is.Not.EqualTo(0));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignGeneralGap("ACGT", "ACG", k => double.NaN));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignGeneralGap("ACGT", "ACG", k => double.PositiveInfinity));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlignGeneralGap("ACGT", "ACG", null!));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignGeneralGap((string)null!, "ACG", k => -k));
        // -infinity forbids a gap length: equal lengths still align without gaps, unequal ones need an allowed length
        var forbidden = SequenceAligner.GlobalAlignGeneralGap("ACGT", "ACGT", k => double.NegativeInfinity);
        Assert.That(forbidden.Score, Is.EqualTo(4)); // SimpleDna: 4 matches
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignGeneralGap("ACGT", "ACG", k => double.NegativeInfinity));
        // a residue outside the substitution matrix alphabet
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignGeneralGap("AC1T", "ACGT", k => -k, SubstitutionMatrix.Blosum62));
    }

    [Test]
    public void EmptyInputs_FollowTheAlignerContract()
    {
        Func<int, double> w = k => -(1 + Math.Log(k));
        var g = SequenceAligner.GlobalAlignGeneralGap("", "ACGT", w);
        Assert.Multiple(() =>
        {
            Assert.That((g.AlignedSequence1, g.AlignedSequence2), Is.EqualTo(("----", "ACGT")));
            Assert.That(g.ExactScore, Is.EqualTo(-(1 + Math.Log(4))).Within(1e-12));
            Assert.That(SequenceAligner.GlobalAlignGeneralGap("", "", w).ExactScore, Is.EqualTo(0));
            var l = SequenceAligner.LocalAlignGeneralGap("", "ACGT", w);
            Assert.That((l.AlignedSequence1, l.Score, l.AlignmentType, l.StartPosition1, l.EndPosition2), Is.EqualTo(("", 0, AlignmentType.Local, -1, -1)));
        });
    }

    [Test]
    public void DnaAndProteinOverloads_EqualTheStringOverloads()
    {
        Func<int, double> w = k => -(2 + 1.5 * Math.Sqrt(k));
        var sc = new ScoringMatrix(2, -3, 0, 0);
        const string a = "ACGTTGCAAGCT", b = "ACGTGCAGCTT";
        AssertSame(SequenceAligner.GlobalAlignGeneralGap(new DnaSequence(a), new DnaSequence(b), w, sc), SequenceAligner.GlobalAlignGeneralGap(a, b, w, sc), "global dna");
        AssertSame(SequenceAligner.LocalAlignGeneralGap(new DnaSequence(a), new DnaSequence(b), w, sc), SequenceAligner.LocalAlignGeneralGap(a, b, w, sc), "local dna");
        AssertSame(SequenceAligner.SemiGlobalAlignGeneralGap(new DnaSequence(a), new DnaSequence(b), w, sc, EndGaps.All), SequenceAligner.SemiGlobalAlignGeneralGap(a, b, w, sc, EndGaps.All), "semi dna");
        const string p = "MKVLAAGIVGLLLAQ", q = "MKVLAGIVGLLLAPQ";
        var blosum = SubstitutionMatrix.Blosum62;
        AssertSame(SequenceAligner.GlobalAlignGeneralGap(new ProteinSequence(p), new ProteinSequence(q), w, blosum), SequenceAligner.GlobalAlignGeneralGap(p, q, w, blosum), "global protein");
        AssertSame(SequenceAligner.LocalAlignGeneralGap(new ProteinSequence(p), new ProteinSequence(q), w, blosum), SequenceAligner.LocalAlignGeneralGap(p, q, w, blosum), "local protein");
        AssertSame(SequenceAligner.SemiGlobalAlignGeneralGap(new ProteinSequence(p), new ProteinSequence(q), w, blosum, EndGaps.OverlapSuffix1Prefix2),
            SequenceAligner.SemiGlobalAlignGeneralGap(p, q, w, blosum, EndGaps.OverlapSuffix1Prefix2), "semi protein");
        AssertSame(SequenceAligner.GlobalAlignGeneralGap(new DnaSequence(a), new DnaSequence(b), w, SubstitutionMatrix.EdnaFull), SequenceAligner.GlobalAlignGeneralGap(a, b, w, SubstitutionMatrix.EdnaFull), "global dna matrix");
    }
}
