using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 (audit A1-7b, B09 F59): enumeration of ALL co-optimal global / ends-free alignments
/// (<see cref="SequenceAligner.GlobalAlignCoOptimal(string, string, ScoringMatrix?, int)"/>,
/// <see cref="SequenceAligner.GlobalAlignAffineCoOptimal(string, string, ScoringMatrix?, int)"/>,
/// <see cref="SequenceAligner.SemiGlobalAlignCoOptimal(string, string, ScoringMatrix?, EndGaps, int)"/> and the matrix overloads).
///
/// Oracle (TestData/Alignment/global_cooptimal_oracle.tsv, scripts/alignment/generate_global_cooptimal_oracle.py): Biopython 1.88
/// <c>PairwiseAligner(mode='global').align(a, b)</c> — the full set of co-optimal alignments for all 16 end-gap policies, linear and
/// affine gaps, scalar scores and BLOSUM62 with fractional gaps — and an independent exhaustive enumeration of every column sequence of
/// tiny inputs (the generator asserts that Biopython's set equals the exhaustive one on every tiny row). The tests below additionally
/// enumerate every column sequence themselves (<see cref="Exhaustive"/>) and re-score every returned alignment independently.
/// </summary>
[TestFixture]
public class SequenceAligner_GlobalCoOptimal_Tests
{
    private sealed record Row(
        string Kind, string Model, int Flags, string Matrix, int Match, int Mismatch, double GapOpen, double GapExtend,
        string Seq1, string Seq2, double Biopython, string Count, string Alignments);

    private static readonly Lazy<List<Row>> Oracle = new(ReadOracle);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.global_cooptimal_oracle.tsv";
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
            rows.Add(new Row(f[0], f[1], int.Parse(f[2], CultureInfo.InvariantCulture), f[3], int.Parse(f[4], CultureInfo.InvariantCulture),
                int.Parse(f[5], CultureInfo.InvariantCulture), D(f[6]), D(f[7]), f[8], f[9], D(f[10]), f[11], f[12]));
        }
        return rows;
    }

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static string RowKey(AlignmentResult r) => r.AlignedSequence1 + "|" + r.AlignedSequence2;

    private static string[] Rows(string alignments) => alignments.Length == 0 || alignments == "NA" ? [] : alignments.Split(';');

    private static SubstitutionMatrix MatrixOf(string name) => name == "BLOSUM62" ? SubstitutionMatrix.Blosum62 : throw new ArgumentException(name);

    private static bool IsAffine(Row r) => r.Model == "affine";

    private static EndGaps Ends(Row r) => (EndGaps)r.Flags;

    private static IReadOnlyList<AlignmentResult> CoOptimal(Row r, int cap)
    {
        bool affine = IsAffine(r);
        var ends = Ends(r);
        if (r.Kind == "scalar")
        {
            var sc = new ScoringMatrix(r.Match, r.Mismatch, affine ? (int)r.GapOpen : 0, (int)r.GapExtend);
            return affine
                ? SequenceAligner.SemiGlobalAlignAffineCoOptimal(r.Seq1, r.Seq2, sc, ends, cap)
                : SequenceAligner.SemiGlobalAlignCoOptimal(r.Seq1, r.Seq2, sc, ends, cap);
        }
        var ms = new MatrixScoring(MatrixOf(r.Matrix), r.GapOpen, r.GapExtend);
        return affine
            ? SequenceAligner.SemiGlobalAlignAffineCoOptimal(r.Seq1, r.Seq2, ms, ends, cap)
            : SequenceAligner.SemiGlobalAlignCoOptimal(r.Seq1, r.Seq2, ms, ends, cap);
    }

    /// <summary>The single (deterministic) aligner of the same row: GlobalAlign* for no free end, SemiGlobalAlign* otherwise.</summary>
    private static AlignmentResult Single(Row r)
    {
        bool affine = IsAffine(r);
        var ends = Ends(r);
        if (r.Kind == "scalar")
        {
            var sc = new ScoringMatrix(r.Match, r.Mismatch, affine ? (int)r.GapOpen : 0, (int)r.GapExtend);
            return affine ? SequenceAligner.SemiGlobalAlignAffine(r.Seq1, r.Seq2, sc, ends) : SequenceAligner.SemiGlobalAlign(r.Seq1, r.Seq2, sc, ends);
        }
        var ms = new MatrixScoring(MatrixOf(r.Matrix), r.GapOpen, r.GapExtend);
        return affine ? SequenceAligner.SemiGlobalAlignAffine(r.Seq1, r.Seq2, ms, ends) : SequenceAligner.SemiGlobalAlign(r.Seq1, r.Seq2, ms, ends);
    }

    private static string Describe(Row r) => $"{r.Kind}/{r.Model}/flags={r.Flags}/{r.Matrix} sc=({r.Match},{r.Mismatch},{r.GapOpen},{r.GapExtend}) {r.Seq1} {r.Seq2}";

    /// <summary>
    /// Independent score + consumed intervals of a gapped row pair under an end-gap policy (written from the definition): the first column
    /// decides which side's leading flank is free (a gap in row 2 = leading residues of sequence 1, FreeStart1; a gap in row 1 = FreeStart2),
    /// the last column decides the trailing flank; internal gaps cost open + e per run of one row (affine) or e per position (linear).
    /// </summary>
    private static (double Score, int Start1, int Start2, int End1, int End2) Rescore(
        string r1, string r2, Func<char, char, double> sub, double go, double ge, bool affine, EndGaps ends, int m, int n)
    {
        int len = r1.Length;
        int lead = 0, trail = 0;
        int start1 = 0, start2 = 0, end1 = m - 1, end2 = n - 1;
        if (len > 0 && r2[0] == '-' && ends.HasFlag(EndGaps.FreeStart1))
        {
            while (lead < len && r2[lead] == '-') lead++;
            start1 = lead;
        }
        else if (len > 0 && r1[0] == '-' && ends.HasFlag(EndGaps.FreeStart2))
        {
            while (lead < len && r1[lead] == '-') lead++;
            start2 = lead;
        }
        if (len > 0 && r2[len - 1] == '-' && ends.HasFlag(EndGaps.FreeEnd1))
        {
            while (trail < len - lead && r2[len - 1 - trail] == '-') trail++;
            end1 = m - 1 - trail;
        }
        else if (len > 0 && r1[len - 1] == '-' && ends.HasFlag(EndGaps.FreeEnd2))
        {
            while (trail < len - lead && r1[len - 1 - trail] == '-') trail++;
            end2 = n - 1 - trail;
        }

        double score = 0;
        char last = ' ';
        for (int k = lead; k < len - trail; k++)
        {
            char x = r1[k], y = r2[k];
            if (x != '-' && y != '-') { score += sub(x, y); last = 'M'; }
            else if (y == '-') { score += !affine || last == 'X' ? ge : go + ge; last = 'X'; }
            else { score += !affine || last == 'Y' ? ge : go + ge; last = 'Y'; }
        }
        return (score, start1, start2, end1, end2);
    }

    private static Func<char, char, double> SubOf(Row r) =>
        r.Kind == "scalar"
            ? (x, y) => x == y ? r.Match : r.Mismatch
            : (x, y) => SubstitutionMatrix.Blosum62.GetScore(x, y);

    // ------------------------------------------------------------------------------------------------ oracle tests

    [Test]
    public void CoOptimalSet_EqualsBiopythonSet_OnEveryOracleRow()
    {
        var failures = new List<string>();
        int compared = 0, multi = 0;
        foreach (var r in Oracle.Value)
        {
            if (r.Alignments == "NA") continue;
            var expected = Rows(r.Alignments).ToHashSet();
            var actual = CoOptimal(r, 10_000).ToList();
            compared++;
            if (expected.Count > 1) multi++;
            var keys = actual.Select(RowKey).ToList();
            if (keys.Count != expected.Count || !keys.ToHashSet().SetEquals(expected) || keys.Distinct().Count() != keys.Count)
                failures.Add($"{Describe(r)}\n   expected {expected.Count}: {string.Join(';', expected.OrderBy(x => x).Take(4))}\n   actual   {keys.Count} ({keys.Distinct().Count()} distinct): {string.Join(';', keys.OrderBy(x => x).Take(4))}");
        }
        Assert.That(compared, Is.GreaterThanOrEqualTo(1800));
        Assert.That(multi, Is.GreaterThanOrEqualTo(800), "the oracle must contain many rows with several co-optimal alignments");
        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(5)));
    }

    [Test]
    public void EveryAlignment_RescoresToTheBiopythonOptimum_WithConsumedIntervals()
    {
        foreach (var r in Oracle.Value)
        {
            if (r.Alignments == "NA") continue;
            var sub = SubOf(r);
            foreach (var a in CoOptimal(r, 10_000))
            {
                var (score, s1, s2, e1, e2) = Rescore(a.AlignedSequence1, a.AlignedSequence2, sub, r.GapOpen, r.GapExtend, IsAffine(r), Ends(r), r.Seq1.Length, r.Seq2.Length);
                Assert.That(score, Is.EqualTo(r.Biopython).Within(1e-9), $"{Describe(r)} {RowKey(a)}");
                Assert.That(a.ExactScore, Is.EqualTo(r.Biopython).Within(1e-9), Describe(r));
                Assert.That((a.StartPosition1, a.StartPosition2, a.EndPosition1, a.EndPosition2), Is.EqualTo((s1, s2, e1, e2)), $"{Describe(r)} {RowKey(a)}");
                Assert.That(a.AlignedSequence1.Replace("-", ""), Is.EqualTo(r.Seq1.ToUpperInvariant()), Describe(r));
                Assert.That(a.AlignedSequence2.Replace("-", ""), Is.EqualTo(r.Seq2.ToUpperInvariant()), Describe(r));
                Assert.That(a.AlignmentType, Is.EqualTo(Ends(r) == EndGaps.None ? AlignmentType.Global : AlignmentType.SemiGlobal), Describe(r));
            }
        }
    }

    [Test]
    public void FirstAlignment_IsTheSingleAligners_OnEveryOracleRow()
    {
        foreach (var r in Oracle.Value)
        {
            var first = CoOptimal(r, 1);
            Assert.That(first, Has.Count.EqualTo(1), Describe(r));
            var single = Single(r);
            // Zero gap costs with a free flank: a flank is also a run of zero-cost gap columns of the path, the single aligner reports
            // the cells its traceback ends in, the enumerator the consumed interval of the rows (a function of the rows alone, checked
            // in EveryAlignment_Rescores...); rows, score and type are identical.
            bool zeroCostFlank = Ends(r) != EndGaps.None && r.GapExtend == 0 && r.GapOpen == 0;
            if (zeroCostFlank)
                Assert.That(first[0] with { StartPosition1 = single.StartPosition1, StartPosition2 = single.StartPosition2, EndPosition1 = single.EndPosition1, EndPosition2 = single.EndPosition2 }, Is.EqualTo(single), Describe(r));
            else
                Assert.That(first[0], Is.EqualTo(single), Describe(r));
        }
    }

    [Test]
    public void Cap_OnRowsWithMoreThanTheOracleLimit_ReturnsDistinctOptimalPrefix()
    {
        int capped = 0;
        foreach (var r in Oracle.Value.Where(x => x.Alignments == "NA"))
        {
            var sub = SubOf(r);
            var a = CoOptimal(r, 50);
            var b = CoOptimal(r, 20);
            capped++;
            Assert.That(a, Has.Count.EqualTo(50), Describe(r));
            Assert.That(a.Select(RowKey).Distinct().Count(), Is.EqualTo(50), Describe(r));
            Assert.That(b.Select(RowKey), Is.EqualTo(a.Take(20).Select(RowKey)), "a smaller cap is a prefix of a larger one (deterministic order)");
            Assert.That(a[0], Is.EqualTo(Single(r)), Describe(r));
            foreach (var x in a)
                Assert.That(Rescore(x.AlignedSequence1, x.AlignedSequence2, sub, r.GapOpen, r.GapExtend, IsAffine(r), Ends(r), r.Seq1.Length, r.Seq2.Length).Score,
                    Is.EqualTo(r.Biopython).Within(1e-9), Describe(r));
            if (int.TryParse(r.Count, out int count))
                Assert.That(count, Is.GreaterThan(150), "NA rows have more alignments than the oracle limit");
        }
        Assert.That(capped, Is.GreaterThanOrEqualTo(30));
    }

    // ------------------------------------------------------------------------------------------------ exhaustive enumeration

    /// <summary>Every gapped column sequence of (a, b) with its score under the model; returns the maximum and the optimal rows.</summary>
    private static (int Best, HashSet<string> Rows) Exhaustive(string a, string b, ScoringMatrix sc, bool affine, EndGaps ends)
    {
        var all = new Dictionary<string, int>();
        var cols = new List<(char X, char Y)>();

        void Rec(int i, int j)
        {
            if (i == a.Length && j == b.Length)
            {
                var r1 = new string(cols.Select(c => c.X).ToArray());
                var r2 = new string(cols.Select(c => c.Y).ToArray());
                all[r1 + "|" + r2] = (int)Rescore(r1, r2, (x, y) => x == y ? sc.Match : sc.Mismatch, sc.GapOpen, sc.GapExtend, affine, ends, a.Length, b.Length).Score;
                return;
            }
            if (i < a.Length && j < b.Length) { cols.Add((a[i], b[j])); Rec(i + 1, j + 1); cols.RemoveAt(cols.Count - 1); }
            if (i < a.Length) { cols.Add((a[i], '-')); Rec(i + 1, j); cols.RemoveAt(cols.Count - 1); }
            if (j < b.Length) { cols.Add(('-', b[j])); Rec(i, j + 1); cols.RemoveAt(cols.Count - 1); }
        }

        Rec(0, 0);
        int best = all.Values.Max();
        return (best, all.Where(kv => kv.Value == best).Select(kv => kv.Key).ToHashSet());
    }

    [Test]
    public void CoOptimalSet_EqualsExhaustiveEnumeration_AllSixteenPolicies_LinearAndAffine()
    {
        var rnd = new Random(20261120);
        var schemes = new[]
        {
            new ScoringMatrix(1, -1, -1, -1), new ScoringMatrix(2, -3, 0, -2), new ScoringMatrix(1, -4, -1, -1), new ScoringMatrix(0, 0, 0, 0),
            new ScoringMatrix(1, 0, 0, -1), new ScoringMatrix(3, -2, -1, 0), new ScoringMatrix(2, -1, -3, -1), new ScoringMatrix(1, -1, 0, 0),
            new ScoringMatrix(2, -3, -3, -1),
        };
        int checkedSets = 0, multi = 0;
        for (int t = 0; t < 140; t++)
        {
            string alpha = new[] { "AC", "ACG", "A" }[rnd.Next(3)];
            string a = new(Enumerable.Range(0, rnd.Next(1, 6)).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray());
            string b = new(Enumerable.Range(0, rnd.Next(1, 6)).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray());
            var sc = schemes[rnd.Next(schemes.Length)];
            for (int flags = 0; flags < 16; flags++)
                foreach (bool affine in new[] { false, true })
                {
                    var scoring = affine ? sc : sc with { GapOpen = 0 };
                    var (best, expected) = Exhaustive(a, b, scoring, affine, (EndGaps)flags);
                    var actual = affine
                        ? SequenceAligner.SemiGlobalAlignAffineCoOptimal(a, b, scoring, (EndGaps)flags, 100_000)
                        : SequenceAligner.SemiGlobalAlignCoOptimal(a, b, scoring, (EndGaps)flags, 100_000);
                    var keys = actual.Select(RowKey).ToList();
                    string what = $"{a} {b} flags={flags} affine={affine} sc=({scoring.Match},{scoring.Mismatch},{scoring.GapOpen},{scoring.GapExtend})";
                    Assert.That(keys.Count, Is.EqualTo(expected.Count), what);
                    Assert.That(keys.ToHashSet().SetEquals(expected), Is.True, what);
                    Assert.That(actual.Select(x => x.Score), Is.All.EqualTo(best), what);
                    checkedSets++;
                    if (expected.Count > 1) multi++;
                }
        }
        Assert.That(checkedSets, Is.EqualTo(140 * 32));
        Assert.That(multi, Is.GreaterThan(1500));
    }

    // ------------------------------------------------------------------------------------------------ examples, overloads, contract

    [Test]
    public void GlobalExample_GaacVsGaaac_ThreeGapPlacements()
    {
        // match 5 / mismatch -3 / gap -1: the single gap can sit at any of the three A's of AAA; order = diagonal > up > left from (m, n)
        var sc = new ScoringMatrix(5, -3, -2, -1);
        var all = SequenceAligner.GlobalAlignCoOptimal("GAAC", "GAAAC", sc);
        Assert.That(all.Select(a => a.AlignedSequence1 + "/" + a.AlignedSequence2), Is.EqualTo(new[] { "G-AAC/GAAAC", "GA-AC/GAAAC", "GAA-C/GAAAC" }));
        Assert.That(all.Select(a => a.Score), Is.All.EqualTo(19));
        Assert.That(all.Select(a => a.AlignmentType), Is.All.EqualTo(AlignmentType.Global));
        Assert.That(all[0], Is.EqualTo(SequenceAligner.GlobalAlign("GAAC", "GAAAC", sc)));

        var affine = SequenceAligner.GlobalAlignAffineCoOptimal("GAAC", "GAAAC", sc);
        Assert.That(affine, Has.Count.EqualTo(3));
        Assert.That(affine.Select(a => a.Score), Is.All.EqualTo(17));   // 4 x 5 - (GapOpen + 1 x GapExtend = -3)
        Assert.That(affine[0], Is.EqualTo(SequenceAligner.GlobalAlignAffine("GAAC", "GAAAC", sc)));
    }

    [Test]
    public void BiopythonTutorialExample_GaactVsGat_TwoAlignments()
    {
        // Biopython 1.88 PairwiseAligner() default (match 1, mismatch 0, all gaps 0): target GAACT, query GAT -> 2 alignments
        // (verified with the installed Biopython: ['GAACT/GA--T', 'GAACT/G-A-T'] as sets of rows, score 3).
        var sc = new ScoringMatrix(1, 0, 0, 0);
        var all = SequenceAligner.GlobalAlignCoOptimal("GAACT", "GAT", sc);
        Assert.That(all.Select(RowKey).OrderBy(x => x), Is.EqualTo(new[] { "GAACT|G-A-T", "GAACT|GA--T" }));
        Assert.That(all.Select(a => a.Score), Is.All.EqualTo(3));
    }

    [Test]
    public void StringDnaProteinAndMatrixOverloads_ReturnTheSameAlignments()
    {
        var sc = new ScoringMatrix(2, -3, -4, -1);
        var fromString = SequenceAligner.GlobalAlignCoOptimal("acgtacgt", "acgacgt", sc);
        var fromDna = SequenceAligner.GlobalAlignCoOptimal(new DnaSequence("ACGTACGT"), new DnaSequence("ACGACGT"), sc);
        Assert.That(fromDna.Select(RowKey), Is.EqualTo(fromString.Select(RowKey)), "strings are upper-cased");

        var ms = new MatrixScoring(SubstitutionMatrix.Blosum62, -9.5, -0.5);
        var s = SequenceAligner.GlobalAlignAffineCoOptimal("MKVLAAGIVG", "MKVLAGIVG", ms);
        var p = SequenceAligner.GlobalAlignAffineCoOptimal(new ProteinSequence("MKVLAAGIVG"), new ProteinSequence("MKVLAGIVG"), ms);
        Assert.That(p.Select(RowKey), Is.EqualTo(s.Select(RowKey)));
        Assert.That(s.Count, Is.GreaterThanOrEqualTo(2));   // the gap can sit at either A
        Assert.That(s[0], Is.EqualTo(SequenceAligner.GlobalAlignAffine("MKVLAAGIVG", "MKVLAGIVG", ms)));
        Assert.That(s.Select(a => a.ExactScore), Is.All.EqualTo(s[0].ExactScore));

        var dnaMatrix = new MatrixScoring(SubstitutionMatrix.EdnaFull, -2.5, -0.5);
        var d = SequenceAligner.SemiGlobalAlignAffineCoOptimal(new DnaSequence("ACGTACGT"), new DnaSequence("TTACGTTT"), dnaMatrix, EndGaps.All);
        Assert.That(d[0], Is.EqualTo(SequenceAligner.SemiGlobalAlignAffine("ACGTACGT", "TTACGTTT", dnaMatrix, EndGaps.All)));
        Assert.That(d.Select(a => a.AlignmentType), Is.All.EqualTo(AlignmentType.SemiGlobal));
    }

    [Test]
    public void EndGapsNone_IsTheGlobalProblem_AndDefaultPolicyOfTheGlobalNames()
    {
        var sc = new ScoringMatrix(1, -1, -2, -1);
        var viaSemi = SequenceAligner.SemiGlobalAlignCoOptimal("ACGTTGCA", "ACGGCA", sc, EndGaps.None);
        var viaGlobal = SequenceAligner.GlobalAlignCoOptimal("ACGTTGCA", "ACGGCA", sc);
        Assert.That(viaSemi, Is.EqualTo(viaGlobal));
    }

    [Test]
    public void EmptySide_HasTheSingleAllGapAlignment()
    {
        var sc = new ScoringMatrix(1, -1, -3, -2);
        var a = SequenceAligner.GlobalAlignCoOptimal("", "ACGT", sc);
        Assert.That(a, Has.Count.EqualTo(1));
        Assert.That((a[0].AlignedSequence1, a[0].AlignedSequence2, a[0].Score), Is.EqualTo(("----", "ACGT", -8)));
        Assert.That(a[0], Is.EqualTo(SequenceAligner.GlobalAlign("", "ACGT", sc)));
        var b = SequenceAligner.GlobalAlignAffineCoOptimal("ACGT", "", sc);
        Assert.That(b, Has.Count.EqualTo(1));
        Assert.That(b[0], Is.EqualTo(SequenceAligner.GlobalAlignAffine("ACGT", "", sc)));
        var both = SequenceAligner.GlobalAlignCoOptimal("", "", sc);
        Assert.That(both, Has.Count.EqualTo(1));
        Assert.That(both[0].Score, Is.EqualTo(0));
        var free = SequenceAligner.SemiGlobalAlignCoOptimal("", "ACGT", sc, EndGaps.All);
        Assert.That(free, Has.Count.EqualTo(1));
        Assert.That(free[0].Score, Is.EqualTo(0));
    }

    [Test]
    public void ZeroGapCosts_FreeFlankIsListedOnce()
    {
        // With free ends and gap costs 0 the same rows would be reachable through several cells; every row pair must be listed once.
        var sc = new ScoringMatrix(1, 0, 0, 0);
        for (int flags = 0; flags < 16; flags++)
        {
            var list = SequenceAligner.SemiGlobalAlignAffineCoOptimal("ACA", "AC", sc, (EndGaps)flags, 100_000);
            Assert.That(list.Select(RowKey).Distinct().Count(), Is.EqualTo(list.Count), $"flags={flags}");
            var lin = SequenceAligner.SemiGlobalAlignCoOptimal("ACA", "AC", sc, (EndGaps)flags, 100_000);
            Assert.That(lin.Select(RowKey).Distinct().Count(), Is.EqualTo(lin.Count), $"flags={flags}");
        }
    }

    [Test]
    public void Validation_NullCapEndGapsAndMemoryGuard()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlignCoOptimal((string)null!, "ACGT"));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlignAffineCoOptimal("ACGT", (string)null!));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlignCoOptimal((DnaSequence)null!, new DnaSequence("ACGT")));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.GlobalAlignCoOptimal("ACGT", "ACGT", (MatrixScoring)null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.GlobalAlignCoOptimal("ACGT", "ACGT", maxAlignments: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.GlobalAlignAffineCoOptimal("ACGT", "ACGT", maxAlignments: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.SemiGlobalAlignCoOptimal("ACGT", "ACGT", (ScoringMatrix?)null, (EndGaps)99));
        // one residue outside the matrix alphabet
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignCoOptimal("MKV1", "MKV", new MatrixScoring(SubstitutionMatrix.Blosum62, -1, -1)));
        // (m + 1)(n + 1) above TopNMaxCells: the O(m n) matrices are never allocated
        string big = new('A', 5_100);
        Assert.That((long)(big.Length + 1) * (big.Length + 1), Is.GreaterThan(SequenceAligner.TopNMaxCells));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignCoOptimal(big, big));
        Assert.Throws<ArgumentException>(() => SequenceAligner.GlobalAlignAffineCoOptimal(big, big));
    }

    [Test]
    public void ManyTies_CapStopsTheEnumeration_Quickly()
    {
        // poly-A against poly-A with a long insertion: C(30, 10) gap placements; the cap must bound the work
        string a = new('A', 30), b = new('A', 20);
        long before = GC.GetAllocatedBytesForCurrentThread();
        var list = SequenceAligner.GlobalAlignCoOptimal(a, b, new ScoringMatrix(1, -1, -2, -1), maxAlignments: 500);
        Assert.That(list, Has.Count.EqualTo(500));
        Assert.That(list.Select(RowKey).Distinct().Count(), Is.EqualTo(500));
        Assert.That(list.Select(x => x.Score), Is.All.EqualTo(20 - 10));
        // A5-4: deterministic proxy of "the cap bounds the work": enumerating all C(30, 10) = 30 045 015 placements would allocate gigabytes
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(256L << 20));
    }
}
