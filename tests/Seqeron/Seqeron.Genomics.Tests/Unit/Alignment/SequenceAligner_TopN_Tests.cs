using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-LOCAL-001 (audit A1-7, B09 F47): Waterman-Eggert (1987) non-intersecting local alignments
/// (<see cref="SequenceAligner.LocalAlignTopN(string, string, ScoringMatrix?, int, int)"/> and affine/matrix variants) and
/// the enumeration of all co-optimal local alignments (<see cref="SequenceAligner.LocalAlignCoOptimal(string, string, ScoringMatrix?, int)"/>).
///
/// Oracle (TestData/Alignment/topn_oracle.tsv, scripts/alignment/generate_topn_oracle.py): an independent full-recompute
/// Waterman-Eggert reference written from the definition, exhaustive enumeration of all local alignments of tiny inputs,
/// EMBOSS matcher 6.6.0 (every reported alignment verified optimal for the definition, round by round), Biopython 1.88
/// <c>PairwiseAligner(mode='local')</c> co-optimal set.
/// Definition: an aligned pair is a diagonal column (a_i opposite b_j); alignment k is a maximum-score local alignment among those
/// that contain none of the pairs of alignments 1..k-1.
/// </summary>
[TestFixture]
public class SequenceAligner_TopN_Tests
{
    private sealed record Row(
        string Kind, string Model, string Matrix, int Match, int Mismatch, double GapOpen, double GapExtend,
        string Seq1, string Seq2, int MaxAlignments, double Biopython, string TopN, string Emboss, string CoOptimal);

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.topn_oracle.tsv";
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
            rows.Add(new Row(f[0], f[1], f[2], int.Parse(f[3], CultureInfo.InvariantCulture), int.Parse(f[4], CultureInfo.InvariantCulture),
                D(f[5]), D(f[6]), f[7], f[8], int.Parse(f[9], CultureInfo.InvariantCulture), D(f[10]), f[11], f[12], f[13]));
        }
        return rows;
    }

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    /// <summary>Canonical text of an alignment: score|start1|start2|end1|end2|row1|row2 (the oracle format).</summary>
    private static string Key(AlignmentResult r) =>
        string.Create(CultureInfo.InvariantCulture, $"{r.ExactScore:R}|{r.StartPosition1}|{r.StartPosition2}|{r.EndPosition1}|{r.EndPosition2}|{r.AlignedSequence1}|{r.AlignedSequence2}");

    private static (double Score, int S1, int S2, int E1, int E2, string R1, string R2) Parse(string entry)
    {
        var f = entry.Split('|');
        return (D(f[0]), int.Parse(f[1], CultureInfo.InvariantCulture), int.Parse(f[2], CultureInfo.InvariantCulture),
            int.Parse(f[3], CultureInfo.InvariantCulture), int.Parse(f[4], CultureInfo.InvariantCulture), f[5], f[6]);
    }

    /// <summary>Oracle entries with the score normalised to the C# round-trip format (the generator prints Python floats, "6.0").</summary>
    private static string[] Entries(string s) => s.Length == 0 ? [] : s.Split(';').Select(e =>
    {
        int bar = e.IndexOf('|');
        return D(e[..bar]).ToString("R", CultureInfo.InvariantCulture) + e[bar..];
    }).ToArray();

    private static SubstitutionMatrix MatrixOf(string name) => name == "BLOSUM62" ? SubstitutionMatrix.Blosum62 : throw new ArgumentException(name);

    private static IReadOnlyList<AlignmentResult> TopN(Row r, int k)
    {
        bool affine = r.Model == "affine";
        if (r.Kind == "scalar")
        {
            var sc = new ScoringMatrix(r.Match, r.Mismatch, affine ? (int)r.GapOpen : 0, (int)r.GapExtend);
            return affine ? SequenceAligner.LocalAlignAffineTopN(r.Seq1, r.Seq2, sc, k) : SequenceAligner.LocalAlignTopN(r.Seq1, r.Seq2, sc, k);
        }
        var ms = new MatrixScoring(MatrixOf(r.Matrix), r.GapOpen, r.GapExtend);
        return affine ? SequenceAligner.LocalAlignAffineTopN(r.Seq1, r.Seq2, ms, k) : SequenceAligner.LocalAlignTopN(r.Seq1, r.Seq2, ms, k);
    }

    private static IReadOnlyList<AlignmentResult> CoOptimal(Row r, int cap)
    {
        bool affine = r.Model == "affine";
        if (r.Kind == "scalar")
        {
            var sc = new ScoringMatrix(r.Match, r.Mismatch, affine ? (int)r.GapOpen : 0, (int)r.GapExtend);
            return affine ? SequenceAligner.LocalAlignAffineCoOptimal(r.Seq1, r.Seq2, sc, cap) : SequenceAligner.LocalAlignCoOptimal(r.Seq1, r.Seq2, sc, cap);
        }
        var ms = new MatrixScoring(MatrixOf(r.Matrix), r.GapOpen, r.GapExtend);
        return affine ? SequenceAligner.LocalAlignAffineCoOptimal(r.Seq1, r.Seq2, ms, cap) : SequenceAligner.LocalAlignCoOptimal(r.Seq1, r.Seq2, ms, cap);
    }

    private static string Describe(Row r) => $"{r.Kind}/{r.Model}/{r.Matrix} sc=({r.Match},{r.Mismatch},{r.GapOpen},{r.GapExtend}) {r.Seq1} {r.Seq2} k={r.MaxAlignments}";

    // ------------------------------------------------------------------------------------------------ oracle tests

    [Test]
    public void TopN_EqualsIndependentWatermanEggertReference_OnEveryOracleRow()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            var expected = Entries(r.TopN);
            var actual = TopN(r, r.MaxAlignments).Select(Key).ToArray();
            if (!expected.SequenceEqual(actual))
                failures.Add($"{Describe(r)}\n   expected {string.Join(';', expected)}\n   actual   {string.Join(';', actual)}");
        }
        Assert.That(Oracle.Value.Count, Is.GreaterThanOrEqualTo(1000));
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void TopN_FirstAlignmentOfEveryRow_IsTheBiopythonOptimum()
    {
        foreach (var r in Oracle.Value)
        {
            var first = TopN(r, 1);
            if (r.Biopython <= 0) { Assert.That(first, Is.Empty, Describe(r)); continue; }
            Assert.That(first, Has.Count.EqualTo(1), Describe(r));
            Assert.That(first[0].ExactScore, Is.EqualTo(r.Biopython).Within(1e-9), Describe(r));
        }
    }

    [Test]
    public void TopN_AgreesWithEmbossMatcherAlternatives()
    {
        // matcher rows are verified in the generator to be a valid Waterman-Eggert sequence (each round optimal, no pair reused).
        // Tie handling differs between tools (a later round depends on WHICH co-optimal alignment was taken), so the whole score list
        // is compared and a divergence is only acceptable on a minority of tie-rich rows; the optimum (round 1) always agrees.
        int rows = 0, equal = 0;
        foreach (var r in Oracle.Value.Where(x => x.Emboss != "NA"))
        {
            var expected = r.Emboss.Split(',').Select(s => D(s)).ToArray();
            var actual = TopN(r, expected.Length).Select(a => a.ExactScore).ToArray();
            Assert.That(actual.Length, Is.GreaterThanOrEqualTo(1), Describe(r));
            Assert.That(actual[0], Is.EqualTo(expected[0]).Within(1e-9), Describe(r));
            rows++;
            if (actual.SequenceEqual(expected)) equal++;
        }
        Assert.That(rows, Is.GreaterThanOrEqualTo(300));
        Assert.That(equal, Is.GreaterThanOrEqualTo((int)(rows * 0.9)), $"{equal} of {rows} matcher score lists equal");
    }

    [Test]
    public void CoOptimal_EqualsBiopythonCoOptimalSet_OnEveryOracleRow()
    {
        int checkedRows = 0;
        var failures = new List<string>();
        foreach (var r in Oracle.Value.Where(x => x.CoOptimal != "NA"))
        {
            var expected = Entries(r.CoOptimal).ToHashSet();
            var actual = CoOptimal(r, 100).Select(Key).ToList();
            checkedRows++;
            if (!expected.SetEquals(actual) || actual.Count != actual.Distinct().Count())
                failures.Add($"{Describe(r)}\n   only Biopython: {string.Join(';', expected.Except(actual))}\n   only Seqeron: {string.Join(';', actual.Except(expected))}");
            else if (expected.Count > 0 && r.Kind == "scalar" && r.Match > 0 && r.Mismatch < 0)
                Assert.That(actual[0], Is.EqualTo(Key(r.Model == "affine"
                    ? SequenceAligner.LocalAlignAffine(r.Seq1, r.Seq2, new ScoringMatrix(r.Match, r.Mismatch, (int)r.GapOpen, (int)r.GapExtend))
                    : SequenceAligner.LocalAlign(r.Seq1, r.Seq2, new ScoringMatrix(r.Match, r.Mismatch, 0, (int)r.GapExtend)))), "first co-optimal = LocalAlign result");
        }
        Assert.That(checkedRows, Is.GreaterThanOrEqualTo(1000));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {checkedRows}: " + string.Join("\n", failures.Take(5)));
    }

    // ------------------------------------------------------------------------------------------------ definition / properties

    private static int Sub(ScoringMatrix s, char x, char y) => x == y ? s.Match : s.Mismatch;

    /// <summary>Score of an alignment re-derived from its rows (linear: GapExtend per gap; affine: GapOpen + k·GapExtend per gap run).</summary>
    private static int Rescore(AlignmentResult r, ScoringMatrix s, bool affine)
    {
        int total = 0;
        char prevKind = 'P';
        for (int c = 0; c < r.AlignedSequence1.Length; c++)
        {
            char x = r.AlignedSequence1[c], y = r.AlignedSequence2[c];
            if (x != '-' && y != '-') { total += Sub(s, x, y); prevKind = 'P'; }
            else
            {
                char kind = x == '-' ? 'Y' : 'X';
                total += s.GapExtend + (affine && kind != prevKind ? s.GapOpen : 0);
                prevKind = kind;
            }
        }
        return total;
    }

    private static string Residues(string aligned) => aligned.Replace("-", "");

    private static void AssertValid(IReadOnlyList<AlignmentResult> list, string a, string b, ScoringMatrix s, bool affine, string what)
    {
        var used = new HashSet<(int, int)>();
        int previousScore = int.MaxValue;
        foreach (var r in list)
        {
            Assert.That(r.AlignmentType, Is.EqualTo(AlignmentType.Local), what);
            Assert.That(r.Score, Is.GreaterThan(0).And.LessThanOrEqualTo(previousScore), what + " non-increasing positive scores");
            previousScore = r.Score;
            Assert.That(Rescore(r, s, affine), Is.EqualTo(r.Score), what + " rescore");
            Assert.That(a.Substring(r.StartPosition1, r.EndPosition1 - r.StartPosition1 + 1), Is.EqualTo(Residues(r.AlignedSequence1)), what + " rows spell a");
            Assert.That(b.Substring(r.StartPosition2, r.EndPosition2 - r.StartPosition2 + 1), Is.EqualTo(Residues(r.AlignedSequence2)), what + " rows spell b");
            foreach (var p in SequenceAligner.AlignedPairs(r))
                Assert.That(used.Add(p), Is.True, $"{what}: aligned pair {p} shared with an earlier alignment");
        }
    }

    /// <summary>Exhaustive maximum over ALL local alignments of tiny inputs that use none of the forbidden pairs (independent of every DP).</summary>
    private static int Brute(string a, string b, ScoringMatrix s, bool affine, HashSet<(int, int)> forbidden)
    {
        int best = 0;
        int m = a.Length, n = b.Length;
        void Rec(int i, int j, int score, char last)
        {
            if (last == 'P' && score > best) best = score;
            if (i < m && j < n && !forbidden.Contains((i + 1, j + 1)))
                Rec(i + 1, j + 1, score + Sub(s, a[i], b[j]), 'P');
            if (i < m) Rec(i + 1, j, score + s.GapExtend + (affine && last != 'X' ? s.GapOpen : 0), 'X');
            if (j < n) Rec(i, j + 1, score + s.GapExtend + (affine && last != 'Y' ? s.GapOpen : 0), 'Y');
        }
        for (int i0 = 0; i0 < m; i0++)
            for (int j0 = 0; j0 < n; j0++)
                if (!forbidden.Contains((i0 + 1, j0 + 1)))
                    Rec(i0 + 1, j0 + 1, Sub(s, a[i0], b[j0]), 'P');
        return best;
    }

    private static readonly ScoringMatrix[] Schemes =
    [
        new(2, -3, -5, -2), new(1, -1, -2, -1), new(5, -4, -10, -1), new(3, -2, 0, -4), new(1, -2, -5, -2), new(2, -1, -4, -3), new(4, -5, 0, -6),
    ];

    [Test]
    public void TopN_EveryRoundIsTheExhaustiveOptimum_OnTinyRandomPairs()
    {
        var rnd = new Random(1987);
        for (int t = 0; t < 600; t++)
        {
            string alpha = t % 3 == 0 ? "AC" : "ACGT";
            string a = string.Concat(Enumerable.Range(0, rnd.Next(2, 6)).Select(_ => alpha[rnd.Next(alpha.Length)]));
            string b = string.Concat(Enumerable.Range(0, rnd.Next(2, 6)).Select(_ => alpha[rnd.Next(alpha.Length)]));
            var s = Schemes[rnd.Next(Schemes.Length)];
            bool affine = rnd.Next(2) == 1;
            var list = affine ? SequenceAligner.LocalAlignAffineTopN(a, b, s, 50) : SequenceAligner.LocalAlignTopN(a, b, s, 50);
            AssertValid(list, a, b, s, affine, $"{a}/{b}/{s}/{affine}");

            var forbidden = new HashSet<(int, int)>();
            foreach (var r in list)
            {
                Assert.That(Brute(a, b, s, affine, forbidden), Is.EqualTo(r.Score), $"{a}/{b}/{s}/{affine}: round is not the optimum");
                foreach (var p in SequenceAligner.AlignedPairs(r)) forbidden.Add(p);
            }
            Assert.That(Brute(a, b, s, affine, forbidden), Is.EqualTo(0), $"{a}/{b}/{s}/{affine}: list ended early");
        }
    }

    [Test]
    public void TopN_FirstAlignmentIsLocalAlign_AndListIsPrefixStable()
    {
        var rnd = new Random(4242);
        for (int t = 0; t < 300; t++)
        {
            string unit = string.Concat(Enumerable.Range(0, rnd.Next(2, 6)).Select(_ => "ACGT"[rnd.Next(4)]));
            string a = Random(rnd, 0, 8) + string.Concat(Enumerable.Repeat(unit, rnd.Next(2, 5))) + Random(rnd, 0, 8);
            string b = Random(rnd, 0, 8) + string.Concat(Enumerable.Repeat(unit, rnd.Next(1, 4))) + Random(rnd, 0, 8);
            var s = Schemes[rnd.Next(Schemes.Length)];
            bool affine = rnd.Next(2) == 1;
            var all = affine ? SequenceAligner.LocalAlignAffineTopN(a, b, s, 25) : SequenceAligner.LocalAlignTopN(a, b, s, 25);
            var one = affine ? SequenceAligner.LocalAlignAffine(a, b, s) : SequenceAligner.LocalAlign(a, b, s);
            if (one.Score == 0) { Assert.That(all, Is.Empty); continue; }
            Assert.That(all[0], Is.EqualTo(one));
            AssertValid(all, a, b, s, affine, $"{a}/{b}");
            var firstThree = affine ? SequenceAligner.LocalAlignAffineTopN(a, b, s, 3) : SequenceAligner.LocalAlignTopN(a, b, s, 3);
            Assert.That(firstThree, Is.EqualTo(all.Take(3).ToList()), "top-3 is the prefix of top-25");
        }
    }

    private static string Random(Random rnd, int min, int max) =>
        string.Concat(Enumerable.Range(0, rnd.Next(min, max)).Select(_ => "ACGT"[rnd.Next(4)]));

    [Test]
    public void TopN_ReportsEveryRepeatCopy_OfATandemRepeat()
    {
        // 3 copies of ACGTAC in sequence 1 against one copy in sequence 2: Smith-Waterman reports one copy, Waterman-Eggert all
        // three (each of 6 matches = 6 under +1/-1/-2 ... here +2/-3): 12, 12, 12 and then nothing that shares no pair.
        var s = new ScoringMatrix(2, -3, -5, -2);
        var list = SequenceAligner.LocalAlignTopN("ACGTACTTTACGTACGGGACGTAC", "ACGTAC", s, 10);
        Assert.That(list.Select(r => r.Score).Take(3), Is.EqualTo(new[] { 12, 12, 12 }));
        Assert.That(list.Take(3).Select(r => r.StartPosition1), Is.EqualTo(new[] { 0, 9, 18 }));
        Assert.That(list.Take(3).Select(r => r.AlignedSequence1), Is.All.EqualTo("ACGTAC"));
        Assert.That(list.Skip(3).All(r => r.Score < 12), Is.True);
    }

    [Test]
    public void TopN_MinScoreAndMaxAlignments_Limit_TheList()
    {
        var s = new ScoringMatrix(2, -3, -5, -2);
        const string a = "ACGTACTTTACGTACGGGACGTAC";
        var all = SequenceAligner.LocalAlignTopN(a, "ACGTAC", s, 50);
        Assert.That(all.Count, Is.GreaterThan(3));
        Assert.That(SequenceAligner.LocalAlignTopN(a, "ACGTAC", s, 2), Has.Count.EqualTo(2));
        var strong = SequenceAligner.LocalAlignTopN(a, "ACGTAC", s, 50, minScore: 12);
        Assert.That(strong, Has.Count.EqualTo(3));
        Assert.That(strong.All(r => r.Score >= 12));
        Assert.That(SequenceAligner.LocalAlignTopN(a, "ACGTAC", s, 50, minScore: 13), Is.Empty);
        Assert.That(SequenceAligner.LocalAlignTopN(a, "ACGTAC", s, 50, minScore: -100).Select(Key), Is.EqualTo(all.Select(Key)), "scores <= 0 are never reported");
    }

    [Test]
    public void TopN_EmptyAndNoMatch_ReturnEmptyList()
    {
        Assert.That(SequenceAligner.LocalAlignTopN("", "ACGT"), Is.Empty);
        Assert.That(SequenceAligner.LocalAlignTopN("ACGT", ""), Is.Empty);
        Assert.That(SequenceAligner.LocalAlignAffineTopN("AAAA", "CCCC"), Is.Empty);
        Assert.That(SequenceAligner.LocalAlignCoOptimal("", "ACGT"), Is.Empty);
        Assert.That(SequenceAligner.LocalAlignAffineCoOptimal("AAAA", "CCCC"), Is.Empty);
    }

    [Test]
    public void Api_Validates_Arguments()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignTopN((string)null!, "A"));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignAffineTopN("A", (string)null!));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignTopN((DnaSequence)null!, new DnaSequence("A")));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignCoOptimal("A", (string)null!));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignTopN("A", "A", (MatrixScoring)null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.LocalAlignTopN("A", "A", (ScoringMatrix?)null, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.LocalAlignAffineTopN("A", "A", (ScoringMatrix?)null, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.LocalAlignCoOptimal("A", "A", (ScoringMatrix?)null, 0));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignTopN("ACGT", "ACGZ", new MatrixScoring(SubstitutionMatrix.EdnaFull, -10, -1)));
    }

    [Test]
    public void TopN_AboveTheCellLimit_UsesTheLinearSpaceProgram_AndCoOptimalStillRejectsIt()
    {
        // A1-7c (F77): above TopNMaxCells the top-N alignments are computed in linear space (no longer rejected);
        // the co-optimal enumerator keeps the matrices and rejects the input with a clear error.
        string big = new('A', 5001);
        var top = SequenceAligner.LocalAlignTopN(big, big, (ScoringMatrix?)null, 2);
        Assert.That(top, Has.Count.EqualTo(2));
        Assert.That(top[0].Score, Is.EqualTo(5001));
        Assert.That(top[1].Score, Is.GreaterThan(0).And.LessThan(5001));
        var ex = Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignCoOptimal(big, big));
        Assert.That(ex!.Message, Does.Contain("25,000,000"));
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignAffineCoOptimal(big, big));
    }

    [Test]
    public void TopN_DnaProteinAndMatrixOverloads_AgreeWithTheStringOverloads()
    {
        var s = new ScoringMatrix(2, -3, -5, -2);
        const string a = "ACGTACTTTACGTACGGGACGTAC", b = "ACGTAC";
        var viaString = SequenceAligner.LocalAlignTopN(a, b, s, 5).Select(Key).ToArray();
        Assert.That(SequenceAligner.LocalAlignTopN(new DnaSequence(a), new DnaSequence(b), s, 5).Select(Key), Is.EqualTo(viaString));
        var viaStringAffine = SequenceAligner.LocalAlignAffineTopN(a, b, s, 5).Select(Key).ToArray();
        Assert.That(SequenceAligner.LocalAlignAffineTopN(new DnaSequence(a), new DnaSequence(b), s, 5).Select(Key), Is.EqualTo(viaStringAffine));

        // lower case is accepted (upper-cased), matrix overloads on DNA and protein
        var ms = new MatrixScoring(SubstitutionMatrix.EdnaFull, -10, -1);
        var m1 = SequenceAligner.LocalAlignAffineTopN(a.ToLowerInvariant(), b, ms, 4).Select(Key).ToArray();
        Assert.That(SequenceAligner.LocalAlignAffineTopN(new DnaSequence(a), new DnaSequence(b), ms, 4).Select(Key), Is.EqualTo(m1));
        Assert.That(SequenceAligner.LocalAlignTopN(new DnaSequence(a), new DnaSequence(b), ms, 4), Has.Count.EqualTo(4));

        var pm = new MatrixScoring(SubstitutionMatrix.Blosum62, -9.5, -0.5);
        var p1 = new ProteinSequence("MKVLAAGIVGLLLAQPSMKVLAAGIVGLLLA");
        var p2 = new ProteinSequence("MKVLAAGIVGLLL");
        var viaProt = SequenceAligner.LocalAlignAffineTopN(p1, p2, pm, 5);
        Assert.That(viaProt.Select(Key), Is.EqualTo(SequenceAligner.LocalAlignAffineTopN(p1.Sequence, p2.Sequence, pm, 5).Select(Key)));
        Assert.That(viaProt[0].ExactScore, Is.EqualTo(SequenceAligner.LocalAlignAffine(p1, p2, pm).ExactScore));
        Assert.That(viaProt[1].ExactScore, Is.EqualTo(viaProt[0].ExactScore), "the second copy of the repeat scores like the first");
        Assert.That(SequenceAligner.LocalAlignTopN(p1, p2, pm, 3)[0].ExactScore, Is.EqualTo(SequenceAligner.LocalAlign(p1, p2, pm).ExactScore));
        Assert.That(SequenceAligner.LocalAlignAffineTopN(p1, p2, pm, 5, minScore: viaProt[0].ExactScore).Count, Is.EqualTo(2));
    }

    // ------------------------------------------------------------------------------------------------ co-optimal enumeration

    [Test]
    public void CoOptimal_HandComputedExample_AllThreeAlignmentsOfAnAmbiguousGap()
    {
        var s = new ScoringMatrix(1, -1, 0, -1);
        // AAAT vs AAT: the unique optimum is the ungapped AAT|AAT (3); every gapped alternative scores 2.
        var list = SequenceAligner.LocalAlignCoOptimal("AAAT", "AAT", s);
        Assert.That(list.Select(r => r.Score), Is.All.EqualTo(3));
        Assert.That(list.Select(r => r.AlignedSequence1 + "/" + r.AlignedSequence2), Is.EquivalentTo(new[] { "AAT/AAT" }));
        // two equally good places: "ACGTACGT" vs "ACGT" -> two copies, each 4 matches
        var two = SequenceAligner.LocalAlignCoOptimal("ACGTACGT", "ACGT", s);
        Assert.That(two.Select(r => r.StartPosition1), Is.EqualTo(new[] { 0, 4 }));
        Assert.That(two.Select(r => r.EndPosition1), Is.EqualTo(new[] { 3, 7 }));
        // a gap that can sit in three places: GAAC vs GAAAC (4 matches x 5 - 1 gap = 19; the ungapped alternatives only reach 15)
        var gap = SequenceAligner.LocalAlignCoOptimal("GAAC", "GAAAC", new ScoringMatrix(5, -3, 0, -1));
        Assert.That(gap.Select(r => r.AlignedSequence1 + "/" + r.AlignedSequence2), Is.EquivalentTo(new[] { "G-AAC/GAAAC", "GA-AC/GAAAC", "GAA-C/GAAAC" }));
        Assert.That(gap.Select(r => r.Score), Is.All.EqualTo(19));
    }

    [Test]
    public void CoOptimal_CapBoundsTheNumberOfAlignments_AndEveryOneIsValidAndDistinct()
    {
        var s = new ScoringMatrix(5, -3, 0, -1);
        string a = "GAAAAAAAAAAC", b = "GAAAAAAAAAAAC";   // G + 10 A + C against G + 11 A + C: 11 placements of the single gap
        var all = SequenceAligner.LocalAlignCoOptimal(a, b, s, 1000);
        Assert.That(all, Has.Count.EqualTo(11));
        Assert.That(all.Select(Key).Distinct().Count(), Is.EqualTo(all.Count));
        var capped = SequenceAligner.LocalAlignCoOptimal(a, b, s, 4);
        Assert.That(capped, Has.Count.EqualTo(4));
        Assert.That(capped.Select(Key), Is.EqualTo(all.Take(4).Select(Key)));
        foreach (var r in all) Assert.That(Rescore(r, s, affine: false), Is.EqualTo(r.Score));
    }

    [Test]
    public void CoOptimal_MatrixOverloads_ReturnExactScores()
    {
        var pm = new MatrixScoring(SubstitutionMatrix.Blosum62, -9.5, -0.5);
        var p1 = new ProteinSequence("MKVLAAGIVGLLLAQPSMKVLAAGIVGLLLA");
        var p2 = new ProteinSequence("MKVLAAGIVGLLL");
        var list = SequenceAligner.LocalAlignAffineCoOptimal(p1, p2, pm);
        Assert.That(list, Has.Count.EqualTo(2));
        double best = SequenceAligner.LocalAlignAffine(p1, p2, pm).ExactScore;
        Assert.That(list.Select(r => r.ExactScore), Is.All.EqualTo(best));
        Assert.That(SequenceAligner.LocalAlignCoOptimal(p1, p2, pm), Has.Count.EqualTo(2));
        Assert.That(SequenceAligner.LocalAlignCoOptimal(new DnaSequence("ACGTACGT"), new DnaSequence("ACGT"), new MatrixScoring(SubstitutionMatrix.EdnaFull, -10, -1)), Has.Count.EqualTo(2));
        Assert.That(SequenceAligner.LocalAlignAffineCoOptimal(new DnaSequence("ACGTACGT"), new DnaSequence("ACGT"), new ScoringMatrix(1, -1, -2, -1)), Has.Count.EqualTo(2));
        Assert.That(SequenceAligner.LocalAlignCoOptimal(new DnaSequence("ACGTACGT"), new DnaSequence("ACGT")), Has.Count.EqualTo(2));
    }

    [Test]
    public void AlignedPairs_AreTheDiagonalColumns_InOneBasedCoordinates()
    {
        var r = new AlignmentResult("AC-GT", "ACTG-", 0, AlignmentType.Local, 2, 5, 5, 8);
        Assert.That(SequenceAligner.AlignedPairs(r), Is.EqualTo(new[] { (3, 6), (4, 7), (5, 9) }.Select(p => (p.Item1, p.Item2))));
        Assert.That(SequenceAligner.AlignedPairs(new AlignmentResult("", "", 0, AlignmentType.Local, -1, -1, -1, -1)), Is.Empty);
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.AlignedPairs(null!));
    }
}
