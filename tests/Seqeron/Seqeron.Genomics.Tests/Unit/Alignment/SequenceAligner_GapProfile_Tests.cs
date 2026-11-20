using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-SEMI-001 / ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit A3-10, B09 F137): per-direction and per-end gap scores — Biopython 1.88 <c>PairwiseAligner</c>'s twelve
/// <c>{insertion,deletion}_{internal,left,right}_{open,extend}</c> scores — on the existing linear / affine cores (<see cref="GapScoreProfile"/>).
///
/// Oracle (TestData/Alignment/gap_profile_oracle.tsv, scripts/alignment/generate_gap_profile_oracle.py): 1000 random pairs x global | local x linear | affine x scalar |
/// EDNAFULL | BLOSUM62 (fractional scores) with random values for all twelve properties (free / equal-to-interior / independently priced ends), <c>PairwiseAligner.score</c>.
/// </summary>
[TestFixture]
public class SequenceAligner_GapProfile_Tests
{
    private sealed record Row(string Seq1, string Seq2, string Mode, bool Affine, string Scoring, double[] Scores, double Biopython, double Model);

    // order of the twelve columns of the oracle
    private static readonly string[] Props =
    {
        "target_internal_open_gap_score", "target_internal_extend_gap_score", "target_left_open_gap_score", "target_left_extend_gap_score",
        "target_right_open_gap_score", "target_right_extend_gap_score", "query_internal_open_gap_score", "query_internal_extend_gap_score",
        "query_left_open_gap_score", "query_left_extend_gap_score", "query_right_open_gap_score", "query_right_extend_gap_score",
    };

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.gap_profile_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], f[2], f[3] == "affine", f[4],
                f[5].Split(',').Select(x => double.Parse(x, CultureInfo.InvariantCulture)).ToArray(), double.Parse(f[6], CultureInfo.InvariantCulture), double.Parse(f[7], CultureInfo.InvariantCulture)));
        }
        return rows;
    }

    private static GapScoreProfile Profile(Row r) =>
        GapScoreProfile.FromBiopython(Props.Zip(r.Scores).ToDictionary(p => p.First, p => p.Second));

    private static AlignmentResult Align(Row r)
    {
        var profile = Profile(r);
        if (r.Scoring.StartsWith("scalar", StringComparison.Ordinal))
        {
            var p = r.Scoring.Split(':');
            var sm = new ScoringMatrix(int.Parse(p[1], CultureInfo.InvariantCulture), int.Parse(p[2], CultureInfo.InvariantCulture), 0, 0);
            return (r.Mode, r.Affine) switch
            {
                ("global", false) => SequenceAligner.GlobalAlign(r.Seq1, r.Seq2, sm, profile),
                ("global", true) => SequenceAligner.GlobalAlignAffine(r.Seq1, r.Seq2, sm, profile),
                (_, false) => SequenceAligner.LocalAlign(r.Seq1, r.Seq2, sm, profile),
                _ => SequenceAligner.LocalAlignAffine(r.Seq1, r.Seq2, sm, profile),
            };
        }

        var ms = new MatrixScoring(r.Scoring == "EDNAFULL" ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62, 0, 0);
        return (r.Mode, r.Affine) switch
        {
            ("global", false) => SequenceAligner.GlobalAlign(r.Seq1, r.Seq2, ms, profile),
            ("global", true) => SequenceAligner.GlobalAlignAffine(r.Seq1, r.Seq2, ms, profile),
            (_, false) => SequenceAligner.LocalAlign(r.Seq1, r.Seq2, ms, profile),
            _ => SequenceAligner.LocalAlignAffine(r.Seq1, r.Seq2, ms, profile),
        };
    }

    // Independent scorer: columns of the rows, a gap run is a start run when it begins at the first column, an end run when it finishes at the last one (global only).
    private static double Rescore(Row r, AlignmentResult res)
    {
        var d = Props.Zip(r.Scores).ToDictionary(p => p.First, p => p.Second);
        double Run(string dir, string pos, int k, bool affine)
        {
            double o = d[$"{dir}_{pos}_open_gap_score"], e = d[$"{dir}_{pos}_extend_gap_score"];
            return affine ? o + (k - 1) * e : k * e; // Biopython convention: first position open, further positions extend (linear rows have open == extend)
        }

        Func<char, char, double> sub;
        if (r.Scoring.StartsWith("scalar", StringComparison.Ordinal))
        {
            var p = r.Scoring.Split(':');
            double match = int.Parse(p[1], CultureInfo.InvariantCulture), mis = int.Parse(p[2], CultureInfo.InvariantCulture);
            sub = (x, y) => x == y ? match : mis;
        }
        else
        {
            var m = r.Scoring == "EDNAFULL" ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62;
            sub = (x, y) => m.ScoreAt(m.IndexOf(x), m.IndexOf(y));
        }

        string a = res.AlignedSequence1, b = res.AlignedSequence2;
        double total = 0;
        int cols = a.Length;
        for (int c = 0; c < cols;)
        {
            if (a[c] != '-' && b[c] != '-') { total += sub(a[c], b[c]); c++; continue; }
            bool gapInSeq2 = b[c] == '-'; // vertical run, query / deletion (seq 1 residues against gaps)
            int start = c;
            while (c < cols && (gapInSeq2 ? b[c] == '-' : a[c] == '-')) c++;
            string pos = r.Mode == "global" && start == 0 ? "left" : r.Mode == "global" && c == cols ? "right" : "internal";
            total += Run(gapInSeq2 ? "query" : "target", pos, c - start, r.Affine);
        }

        return total;
    }

    [Test]
    public void GapProfile_Score_EqualsBiopython_AllRows()
    {
        var rows = Oracle.Value;
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(500));
        int global = 0, local = 0, affine = 0, matrix = 0, freeEnds = 0, corner = 0;
        foreach (var r in rows)
        {
            var res = Align(r);
            double actual = res.ExactScore;
            string label = $"{r.Mode} {(r.Affine ? "affine" : "linear")} {r.Scoring} {r.Seq1} / {r.Seq2} [{string.Join(",", r.Scores)}]";
            Assert.That(actual, Is.EqualTo(r.Model).Within(1e-9), label);
            if (Math.Abs(r.Model - r.Biopython) < 1e-9)
                Assert.That(actual, Is.EqualTo(r.Biopython).Within(1e-9), label);
            else
                corner++; // the documented Biopython 1.88 linear Needleman-Wunsch corner (see the generator header)
            if (r.Mode == "global") global++; else local++;
            if (r.Affine) affine++;
            if (!r.Scoring.StartsWith("scalar", StringComparison.Ordinal)) matrix++;
            if (r.Scores.Where((_, i) => i is >= 2 and <= 5 or >= 8).Any(v => v == 0)) freeEnds++;
        }

        // 998 of the 1000 Biopython scores are equal; the other two are the linear Needleman-Wunsch corner of Biopython 1.88: first column of the last row charged `nA * right_gap_extend_B`
        // (an all-gap alignment whose leading run consumes ALL of sequence 1), where its own Gotoh routine and Seqeron use the start score; the generator reproduces them with that rule.
        Assert.That(corner, Is.EqualTo(2));
        Assert.That(rows.Count - corner, Is.GreaterThanOrEqualTo(990));
        Assert.That(global, Is.GreaterThan(700));
        Assert.That(local, Is.GreaterThan(150));
        Assert.That(affine, Is.GreaterThan(500));
        Assert.That(matrix, Is.GreaterThan(250));
        Assert.That(freeEnds, Is.GreaterThan(300));
    }

    [Test]
    public void GapProfile_Rows_AreAnAlignmentOfTheInputsWithTheReportedScore()
    {
        foreach (var r in Oracle.Value)
        {
            var res = Align(r);
            double reported = res.ExactScore;
            Assert.That(Rescore(r, res), Is.EqualTo(reported).Within(1e-9), $"{r.Mode} {r.Scoring} {r.Seq1} / {r.Seq2} rows {res.AlignedSequence1} / {res.AlignedSequence2}");
            Assert.That(res.AlignedSequence1.Length, Is.EqualTo(res.AlignedSequence2.Length));
            if (r.Mode == "global")
            {
                Assert.That(res.AlignedSequence1.Replace("-", ""), Is.EqualTo(r.Seq1.ToUpperInvariant()));
                Assert.That(res.AlignedSequence2.Replace("-", ""), Is.EqualTo(r.Seq2.ToUpperInvariant()));
            }
            else
            {
                Assert.That(r.Seq1.ToUpperInvariant().Contains(res.AlignedSequence1.Replace("-", ""), StringComparison.Ordinal), Is.True);
                Assert.That(r.Seq2.ToUpperInvariant().Contains(res.AlignedSequence2.Replace("-", ""), StringComparison.Ordinal), Is.True);
            }
        }
    }

    [Test]
    public void EndGapsPolicies_AreFreeMembersOfAProfile_AllSixteen()
    {
        // GapScoreProfile.FromEndGaps(flags) reproduces the score of every EndGaps policy of the legacy ends-free cores (linear + affine, scalar + matrix).
        var rnd = new Random(1234);
        const string alpha = "ACGT";
        string Make(int n) => new(Enumerable.Range(0, n).Select(_ => alpha[rnd.Next(4)]).ToArray());
        var scalar = new ScoringMatrix(2, -3, -4, -2);
        var matrix = new MatrixScoring(SubstitutionMatrix.EdnaFull, -2.5, -0.5);
        int compared = 0;
        for (int flags = 0; flags < 16; flags++)
        {
            var ends = (EndGaps)flags;
            var profile = GapScoreProfile.FromEndGaps(ends);
            for (int k = 0; k < 14; k++)
            {
                string a = Make(rnd.Next(1, 18)), b = Make(rnd.Next(1, 18));
                Assert.That(SequenceAligner.SemiGlobalAlign(a, b, scalar, profile).Score, Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b, scalar, ends).Score), $"linear {ends} {a} {b}");
                Assert.That(SequenceAligner.SemiGlobalAlignAffine(a, b, scalar, profile).Score, Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine(a, b, scalar, ends).Score), $"affine {ends} {a} {b}");
                var lm = SequenceAligner.SemiGlobalAlign(a, b, matrix, profile);
                Assert.That(lm.ExactScore, Is.EqualTo(SequenceAligner.SemiGlobalAlign(a, b, matrix, ends).ExactScore).Within(1e-9), $"matrix linear {ends} {a} {b}");
                var am = SequenceAligner.SemiGlobalAlignAffine(a, b, matrix, profile);
                Assert.That(am.ExactScore, Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine(a, b, matrix, ends).ExactScore).Within(1e-9), $"matrix affine {ends} {a} {b}");
                compared += 4;
            }
        }

        Assert.That(compared, Is.EqualTo(16 * 14 * 4));
    }

    [Test]
    public void DefaultProfile_IsTheOrdinaryAlignment_RowForRow()
    {
        var rnd = new Random(77);
        var sc = new ScoringMatrix(1, -1, -2, -1);
        for (int k = 0; k < 60; k++)
        {
            string a = new(Enumerable.Range(0, rnd.Next(1, 25)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            string b = new(Enumerable.Range(0, rnd.Next(1, 25)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var plain = SequenceAligner.GlobalAlign(a, b, sc);
            var viaProfile = SequenceAligner.GlobalAlign(a, b, sc, new GapScoreProfile());
            Assert.That((viaProfile.AlignedSequence1, viaProfile.AlignedSequence2, viaProfile.Score), Is.EqualTo((plain.AlignedSequence1, plain.AlignedSequence2, plain.Score)));
            var plainAffine = SequenceAligner.GlobalAlignAffine(a, b, sc);
            var viaAffine = SequenceAligner.GlobalAlignAffine(a, b, sc, new GapScoreProfile());
            Assert.That((viaAffine.AlignedSequence1, viaAffine.AlignedSequence2, viaAffine.Score), Is.EqualTo((plainAffine.AlignedSequence1, plainAffine.AlignedSequence2, plainAffine.Score)));
            var plainLocal = SequenceAligner.LocalAlignAffine(a, b, sc);
            var viaLocal = SequenceAligner.LocalAlignAffine(a, b, sc, new GapScoreProfile());
            Assert.That((viaLocal.AlignedSequence1, viaLocal.AlignedSequence2, viaLocal.Score), Is.EqualTo((plainLocal.AlignedSequence1, plainLocal.AlignedSequence2, plainLocal.Score)));
        }
    }

    [Test]
    public void Coordinates_FreeOverhangIsExcluded_PricedEndStaysInside()
    {
        // Biopython end_insertion_score = 0 on GATTACA / CCGATTTACAGG: the fitted interval is ref[2..9] (F40); the profile with free Start2 / End2 gives the same interval.
        var free = GapScoreProfile.FromEndGaps(EndGaps.Fitting);
        var fit = SequenceAligner.SemiGlobalAlign("GATTACA", "CCGATTTACAGG", new ScoringMatrix(2, -3, -2, -2), free);
        Assert.That((fit.AlignmentType, fit.StartPosition1, fit.EndPosition1, fit.StartPosition2, fit.EndPosition2), Is.EqualTo((AlignmentType.SemiGlobal, 0, 6, 2, 9)));
        // all ends priced: Global, full spans
        var priced = SequenceAligner.GlobalAlign("GATTACA", "CCGATTTACAGG", new ScoringMatrix(2, -3, -2, -2), new GapScoreProfile());
        Assert.That((priced.AlignmentType, priced.StartPosition2, priced.EndPosition2), Is.EqualTo((AlignmentType.Global, 0, 11)));
    }

    [Test]
    public void WorkedExample_LeftGapInTarget_ChargedWithTheLeftScore()
    {
        // run in Biopython 1.88: mode global, match 1, mismatch -1, open/extend -5, left_insertion = 0: score("AC", "GGAC") = 2; priced -5 -5 +... = -8 (GG gap: -5 + -5 + 2)
        var sc = new ScoringMatrix(1, -1, 0, 0);
        var freeLeft = new GapScoreProfile(Interior1: GapScore.FromBiopython(-5, -5), Interior2: GapScore.FromBiopython(-5, -5), Start1: GapScore.FromBiopython(-5, -5), End1: GapScore.FromBiopython(-5, -5),
            Start2: GapScore.Free, End2: GapScore.FromBiopython(-5, -5));
        Assert.That(SequenceAligner.GlobalAlignAffine("AC", "GGAC", sc, freeLeft).Score, Is.EqualTo(2));
        Assert.That(SequenceAligner.GlobalAlignAffine("GGAC", "AC", sc, freeLeft).Score, Is.EqualTo(-8));
    }

    [Test]
    public void FromBiopython_AcceptsBothNameFamilies_AndAggregates()
    {
        // new names: insertion = gap in the target (sequence 1), deletion = gap in the query
        var p = GapScoreProfile.FromBiopython(new Dictionary<string, double>
        {
            ["open_gap_score"] = -4,
            ["extend_gap_score"] = -1,
            ["open_left_insertion_score"] = 0,
            ["extend_left_insertion_score"] = 0,
            ["end_deletion_score"] = -2,
        });
        Assert.That(p.Start2, Is.EqualTo(GapScore.Free));
        Assert.That(p.Start1, Is.EqualTo(GapScore.FromBiopython(-2, -2)));
        Assert.That(p.End1, Is.EqualTo(GapScore.FromBiopython(-2, -2)));
        Assert.That(p.Interior1, Is.EqualTo(GapScore.FromBiopython(-4, -1)));
        Assert.That(p.End2, Is.EqualTo(GapScore.FromBiopython(-4, -1)));
        // deprecated names resolve to the same slots
        var q = GapScoreProfile.FromBiopython(new Dictionary<string, double>
        {
            ["target_left_open_gap_score"] = -3, ["target_left_extend_gap_score"] = -1, ["query_right_open_gap_score"] = -6, ["query_right_extend_gap_score"] = -2,
        });
        Assert.That(q.Start2, Is.EqualTo(GapScore.FromBiopython(-3, -1)));
        Assert.That(q.End1, Is.EqualTo(GapScore.FromBiopython(-6, -2)));
        Assert.That(q.Interior1, Is.Null);
        // an unknown / non-finite property is refused
        Assert.Throws<ArgumentException>(() => GapScoreProfile.FromBiopython(new Dictionary<string, double> { ["left_gap_scorex"] = 1 }));
        Assert.Throws<ArgumentException>(() => GapScoreProfile.FromBiopython(new Dictionary<string, double> { ["open_gap_score"] = double.NaN }));
        Assert.Throws<ArgumentException>(() => GapScoreProfile.FromBiopython(new Dictionary<string, double> { ["insertion_deletion_score"] = 1 }));
    }

    [Test]
    public void Validation_RefusesPositiveGapScores_AndTheLinearSpaceRange()
    {
        var sc = new ScoringMatrix(1, -1, -2, -1);
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACGT", "ACG", sc, new GapScoreProfile(Start1: new GapScore(0, 1))));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACGT", "ACG", sc, new GapScoreProfile(End2: new GapScore(2, -1))));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlign("ACGT", "ACG", sc, new GapScoreProfile(End2: new GapScore(double.NaN, -1))));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlign("ACGT", "ACG", sc, (GapScoreProfile)null!));
        // an input beyond the full-matrix limit is no longer refused (A4-8): the exact linear-space program takes it (tests below)
        string big = new('A', 5001);
        Assert.That(SequenceAligner.GlobalAlignAffine(big, big, sc, new GapScoreProfile(End1: GapScore.Free)).AlignedSequence1, Is.EqualTo(big));
    }

    [Test]
    public void FractionalProfile_ScalesExactly_AndReportsExactScore()
    {
        var ms = new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -0.5);
        var profile = new GapScoreProfile(Start2: GapScore.Free, End2: GapScore.FromBiopython(-1.5, -0.25));
        var r = SequenceAligner.GlobalAlignAffine("HEAGAWGHEE", "PAWHEAEHEE", ms, profile);
        Assert.That(r.ExactScore, Is.EqualTo(21.5).Within(1e-9), "Biopython 1.88 PairwiseAligner: BLOSUM62, open -10.5 / extend -0.5, left insertion 0, right insertion open -1.5 / extend -0.25"); // a multiple of 1/4
    }

    // ------------------------------------------------------------------ A4-8: the exact linear-space programs with a gap profile

    private static GapScore? RandomGap(Random rnd, bool affine)
    {
        int k = rnd.Next(10);
        if (k < 3) return null;                       // the interior score of the scoring
        if (k < 5) return GapScore.Free;
        double open = affine ? -rnd.Next(0, 7) * 0.5 : 0;
        double extend = -rnd.Next(0, 5) * 0.25;
        return new GapScore(open, extend);
    }

    private static string RandomResidues(Random rnd, string alphabet, int n)
    {
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = alphabet[rnd.Next(alphabet.Length)];
        return new string(c);
    }

    private static string Mutated(Random rnd, string s, string alphabet)
    {
        var sb = new System.Text.StringBuilder();
        foreach (char c in s)
        {
            int r = rnd.Next(100);
            if (r < 8) sb.Append(alphabet[rnd.Next(alphabet.Length)]);
            else if (r < 13) { /* deletion */ }
            else if (r < 18) { sb.Append(c).Append(alphabet[rnd.Next(alphabet.Length)]); }
            else sb.Append(c);
        }

        return sb.ToString();
    }

    /// <summary>The full-matrix result and the forced linear-space result of one method family.</summary>
    private static (AlignmentResult Full, AlignmentResult Linear) Both(int family, string a, string b, bool matrix, GapScoreProfile profile, long baseCells)
    {
        var sm = new ScoringMatrix(2, -3, -4, -1);
        var ms = new MatrixScoring(SubstitutionMatrix.Blosum62, -3.5, -0.5);
        return (family, matrix) switch
        {
            (0, false) => (SequenceAligner.GlobalAlign(a, b, sm, profile), SequenceAligner.GlobalAlignLinearSpace(a, b, sm, profile, baseCells)),
            (1, false) => (SequenceAligner.GlobalAlignAffine(a, b, sm, profile), SequenceAligner.GlobalAlignAffineLinearSpace(a, b, sm, profile, baseCells)),
            (2, false) => (SequenceAligner.SemiGlobalAlign(a, b, sm, profile), SequenceAligner.SemiGlobalAlignLinearSpace(a, b, sm, profile, baseCells)),
            (3, false) => (SequenceAligner.SemiGlobalAlignAffine(a, b, sm, profile), SequenceAligner.SemiGlobalAlignAffineLinearSpace(a, b, sm, profile, baseCells)),
            (4, false) => (SequenceAligner.LocalAlign(a, b, sm, profile), SequenceAligner.LocalAlignLinearSpace(a, b, sm, profile, baseCells)),
            (5, false) => (SequenceAligner.LocalAlignAffine(a, b, sm, profile), SequenceAligner.LocalAlignAffineLinearSpace(a, b, sm, profile, baseCells)),
            (0, true) => (SequenceAligner.GlobalAlign(a, b, ms, profile), SequenceAligner.GlobalAlignLinearSpace(a, b, ms, profile, baseCells)),
            (1, true) => (SequenceAligner.GlobalAlignAffine(a, b, ms, profile), SequenceAligner.GlobalAlignAffineLinearSpace(a, b, ms, profile, baseCells)),
            (2, true) => (SequenceAligner.SemiGlobalAlign(a, b, ms, profile), SequenceAligner.SemiGlobalAlignLinearSpace(a, b, ms, profile, baseCells)),
            (3, true) => (SequenceAligner.SemiGlobalAlignAffine(a, b, ms, profile), SequenceAligner.SemiGlobalAlignAffineLinearSpace(a, b, ms, profile, baseCells)),
            (4, true) => (SequenceAligner.LocalAlign(a, b, ms, profile), SequenceAligner.LocalAlignLinearSpace(a, b, ms, profile, baseCells)),
            _ => (SequenceAligner.LocalAlignAffine(a, b, ms, profile), SequenceAligner.LocalAlignAffineLinearSpace(a, b, ms, profile, baseCells)),
        };
    }

    [Test]
    public void LinearSpace_ForcedRecursion_EqualsTheFullMatrix_RowForRow()
    {
        var rnd = new Random(20261201);
        int cases = 0, freeEnds = 0, pricedEnds = 0, withGaps = 0;
        long[] bases = { 1, 2, 9, 40, 300 };
        for (int t = 0; t < 600; t++)
        {
            int family = t % 6;
            bool matrix = (t / 6) % 2 == 1;
            bool affine = family is 1 or 3 or 5;
            string alphabet = matrix ? "ACDEFGHIKLMNPQRSTVWY" : "ACGT";
            string a = RandomResidues(rnd, alphabet, rnd.Next(1, 90));
            string b = rnd.Next(4) == 0 ? RandomResidues(rnd, alphabet, rnd.Next(1, 90)) : Mutated(rnd, a, alphabet);
            if (b.Length == 0) b = "A";
            var profile = new GapScoreProfile(RandomGap(rnd, affine), RandomGap(rnd, affine), RandomGap(rnd, affine), RandomGap(rnd, affine), RandomGap(rnd, affine), RandomGap(rnd, affine));
            var (full, linear) = Both(family, a, b, matrix, profile, bases[t % bases.Length]);
            Assert.That(linear, Is.EqualTo(full), $"family {family} matrix {matrix} base {bases[t % bases.Length]}: {a} / {b} {profile}");
            cases++;
            if (profile.HasFreeEnd) freeEnds++; else pricedEnds++;
            if (full.AlignedSequence1.Contains('-') || full.AlignedSequence2.Contains('-')) withGaps++;
        }

        Assert.That(cases, Is.EqualTo(600));
        Assert.That(freeEnds, Is.GreaterThan(300));
        Assert.That(pricedEnds, Is.GreaterThan(30));
        Assert.That(withGaps, Is.GreaterThan(400));
    }

    [Test]
    public void LinearSpace_ForcedRecursion_Oracle_EqualsBiopython_AllRows()
    {
        // the 1000 Biopython rows with the base-case size forced to ONE cell (a sub-problem is solved by the full-matrix core only when one side has <= 1 residue)
        var rows = Oracle.Value;
        int checkedRows = 0;
        foreach (var r in rows)
        {
            var full = Align(r);
            var profile = Profile(r);
            AlignmentResult forced;
            if (r.Scoring.StartsWith("scalar", StringComparison.Ordinal))
            {
                var p = r.Scoring.Split(':');
                var sm = new ScoringMatrix(int.Parse(p[1], CultureInfo.InvariantCulture), int.Parse(p[2], CultureInfo.InvariantCulture), 0, 0);
                forced = (r.Mode, r.Affine) switch
                {
                    ("global", false) => SequenceAligner.GlobalAlignLinearSpace(r.Seq1, r.Seq2, sm, profile, 1),
                    ("global", true) => SequenceAligner.GlobalAlignAffineLinearSpace(r.Seq1, r.Seq2, sm, profile, 1),
                    (_, false) => SequenceAligner.LocalAlignLinearSpace(r.Seq1, r.Seq2, sm, profile, 1),
                    _ => SequenceAligner.LocalAlignAffineLinearSpace(r.Seq1, r.Seq2, sm, profile, 1),
                };
            }
            else
            {
                var ms = new MatrixScoring(r.Scoring == "EDNAFULL" ? SubstitutionMatrix.EdnaFull : SubstitutionMatrix.Blosum62, 0, 0);
                forced = (r.Mode, r.Affine) switch
                {
                    ("global", false) => SequenceAligner.GlobalAlignLinearSpace(r.Seq1, r.Seq2, ms, profile, 1),
                    ("global", true) => SequenceAligner.GlobalAlignAffineLinearSpace(r.Seq1, r.Seq2, ms, profile, 1),
                    (_, false) => SequenceAligner.LocalAlignLinearSpace(r.Seq1, r.Seq2, ms, profile, 1),
                    _ => SequenceAligner.LocalAlignAffineLinearSpace(r.Seq1, r.Seq2, ms, profile, 1),
                };
            }

            string label = $"{r.Mode} {(r.Affine ? "affine" : "linear")} {r.Scoring} {r.Seq1} / {r.Seq2}";
            Assert.That(forced, Is.EqualTo(full), label);
            Assert.That(forced.ExactScore, Is.EqualTo(r.Model).Within(1e-9), label);
            checkedRows++;
        }

        Assert.That(checkedRows, Is.GreaterThanOrEqualTo(500));
    }

    // deterministic xorshift32 pair, mirrored in Python for the Biopython oracle of the large smoke test below
    private static (string A, string B) XorshiftPair()
    {
        uint x = 12345;
        uint Next() { x ^= x << 13; x ^= x >> 17; x ^= x << 5; return x; }
        string Gen(int n, uint seed, out string mutated)
        {
            x = seed;
            var sbA = new System.Text.StringBuilder();
            for (int i = 0; i < n; i++) sbA.Append("ACGT"[(int)(Next() & 3)]);
            string a = sbA.ToString();
            var sb = new System.Text.StringBuilder();
            foreach (char c in a)
            {
                uint r = Next() % 100;
                if (r < 4) sb.Append("ACGT"[(int)(Next() & 3)]);
                else if (r < 6) { }
                else if (r < 8) { sb.Append(c); sb.Append("ACGT"[(int)(Next() & 3)]); }
                else sb.Append(c);
            }

            mutated = sb.ToString();
            return a;
        }

        string a0 = Gen(8000, 12345, out string b);
        string flank = Gen(700, 777, out _);
        return (flank + a0, b);
    }

    [Test]
    public void LargeInput_AboveTheFullMatrixLimit_IsAlignedByTheLinearSpaceProgram_WithBiopythonScores_AndLinearMemory()
    {
        var (a, b) = XorshiftPair();
        long cells = (long)(a.Length + 1) * (b.Length + 1);
        Assert.That(cells, Is.GreaterThan(SequenceAligner.LinearSpaceCellThreshold * 2)); // 69.6 M cells: the automatic switch applies
        var sm = new ScoringMatrix(2, -3, -5, -2);
        var affineProfile = new GapScoreProfile(Start1: GapScore.Free, End1: GapScore.FromBiopython(-3, -1), End2: GapScore.Free);
        var linearScoring = new ScoringMatrix(2, -3, 0, -2);
        var linearProfile = new GapScoreProfile(Start1: GapScore.Free, End1: GapScore.FromBiopython(-1, -1), End2: GapScore.Free);

        SequenceAligner.GlobalAlignAffine(a.Substring(0, 300), b.Substring(0, 300), sm, affineProfile); // JIT warm-up
        long before = GC.GetAllocatedBytesForCurrentThread();
        var globalAffine = SequenceAligner.GlobalAlignAffine(a, b, sm, affineProfile);
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        long fullBytes = cells * 12; // three int matrices of the full-matrix affine program: 835 MB
        TestContext.Out.WriteLine(string.Create(CultureInfo.InvariantCulture, $"A4-8 memory: {cells} cells, allocated {allocated} bytes (the full affine matrices: {fullBytes} bytes)"));
        Assert.That(allocated, Is.LessThan(fullBytes / 4), () => $"allocated {allocated:N0} bytes, full matrices {fullBytes:N0}");

        // Biopython 1.88 PairwiseAligner on the same pair (scripts: see the test header of this region): open -7 / extend -2, query_left 0, query_right -3 / -1, target_right 0
        Assert.That(globalAffine.ExactScore, Is.EqualTo(12242.0));
        Assert.That(SequenceAligner.LocalAlignAffine(a, b, sm, affineProfile).ExactScore, Is.EqualTo(12242.0));
        var globalLinear = SequenceAligner.GlobalAlign(a, b, linearScoring, linearProfile);
        // global linear: 13824 = an independent numpy dynamic program of the same scheme (prefix-maximum row update) = Biopython's own LOCAL score of the pair; Biopython 1.88's linear Needleman-Wunsch
        // routine itself prints 12499 here: the documented corner of the first column of the last row (see the oracle above), not a difference of the programs
        Assert.That(globalLinear.ExactScore, Is.EqualTo(13824.0));
        Assert.That(SequenceAligner.LocalAlign(a, b, linearScoring, linearProfile).ExactScore, Is.EqualTo(13824.0));

        // the rows are a legal alignment of the two sequences and independent of the size of the full-matrix leaves
        Assert.That(globalAffine.AlignedSequence1.Replace("-", ""), Is.EqualTo(a));
        Assert.That(globalAffine.AlignedSequence2.Replace("-", ""), Is.EqualTo(b));
        Assert.That(SequenceAligner.GlobalAlignAffineLinearSpace(a, b, sm, affineProfile, 1L << 18), Is.EqualTo(globalAffine));
    }
}
