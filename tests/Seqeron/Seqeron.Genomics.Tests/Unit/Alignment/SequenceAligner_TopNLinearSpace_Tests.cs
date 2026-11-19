using System.Diagnostics;
using System.Globalization;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-LOCAL-001 (audit A1-7c, B09 F77): Waterman-Eggert top-N local alignments in LINEAR space
/// (<see cref="SequenceAligner.LocalAlignTopNLinearSpace(string, string, ScoringMatrix?, int, int, long)"/>, the forbidden-pair predicate threaded through
/// the locate passes and the interval solves of the exact linear-space local program; automatic above
/// <see cref="SequenceAligner.TopNMaxCells"/> cells). The in-memory declumper (the 1987 matrices) and the independent oracle of
/// <c>SequenceAligner_TopN_Tests</c> (full-recompute Waterman-Eggert reference, EMBOSS matcher 6.6.0, Biopython 1.88) are the references:
/// the linear-space program must reproduce every alignment (rows, score, coordinates, order, ties).
/// </summary>
[TestFixture]
public class SequenceAligner_TopNLinearSpace_Tests
{
    private static string Key(AlignmentResult r) =>
        string.Create(CultureInfo.InvariantCulture, $"{r.Score}|{r.ExactScore:R}|{r.StartPosition1}|{r.StartPosition2}|{r.EndPosition1}|{r.EndPosition2}|{r.AlignedSequence1}|{r.AlignedSequence2}");

    private static string[] Keys(IEnumerable<AlignmentResult> list) => list.Select(Key).ToArray();

    private static string Random(Random rng, string alphabet, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++) chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    /// <summary>A pair with shared (mutated) segments and repeats so that several non-intersecting alignments exist.</summary>
    private static (string A, string B) Pair(Random rng, string alphabet, int maxLen)
    {
        string unit = Random(rng, alphabet, rng.Next(1, 6));
        var a = new System.Text.StringBuilder();
        var b = new System.Text.StringBuilder();
        int blocks = rng.Next(1, 5);
        for (int k = 0; k < blocks; k++)
        {
            string seg = Random(rng, alphabet, rng.Next(2, Math.Max(3, maxLen / 4)));
            a.Append(Random(rng, alphabet, rng.Next(0, 6))).Append(seg);
            var mutated = seg.ToCharArray();
            for (int i = 0; i < mutated.Length; i++)
                if (rng.Next(6) == 0) mutated[i] = alphabet[rng.Next(alphabet.Length)];
            b.Append(Random(rng, alphabet, rng.Next(0, 6))).Append(new string(mutated));
            if (rng.Next(3) == 0) { a.Append(unit).Append(unit); b.Append(unit); }
        }
        a.Append(Random(rng, alphabet, rng.Next(0, 8)));
        b.Append(Random(rng, alphabet, rng.Next(0, 8)));
        return (a.ToString(), b.ToString());
    }

    private static readonly long[] BaseCellSizes = [1, 9, 40, 200];

    // ----------------------------------------------------------------------- 300+ random pairs, linear path forced

    [Test]
    public void LinearSpace_ScalarLinearAndAffine_EqualTheInMemoryDeclumper_OnRandomPairs()
    {
        var rng = new Random(20261119);
        var models = new (string Name, ScoringMatrix Sc)[]
        {
            ("1/-1/-2", new ScoringMatrix(1, -1, 0, -2)),
            ("2/-3/-4", new ScoringMatrix(2, -3, 0, -4)),
            ("5/-4/-10/-1", new ScoringMatrix(5, -4, -10, -1)),
            ("1/0/-1/-1", new ScoringMatrix(1, 0, -1, -1)),
            ("3/-2/0/0", new ScoringMatrix(3, -2, 0, 0)),
        };
        string[] alphabets = ["AC", "ACG", "ACGT"];
        int pairs = 0;
        var failures = new List<string>();
        for (int iter = 0; iter < 160; iter++)
        {
            var (a, b) = Pair(rng, alphabets[iter % alphabets.Length], 36);
            var (name, sc) = models[iter % models.Length];
            int k = rng.Next(1, 12);
            int min = iter % 4 == 0 ? rng.Next(0, 4) : 0;
            long baseCells = BaseCellSizes[iter % BaseCellSizes.Length];
            foreach (bool affine in new[] { false, true })
            {
                var expected = Keys(affine ? SequenceAligner.LocalAlignAffineTopN(a, b, sc, k, min) : SequenceAligner.LocalAlignTopN(a, b, sc, k, min));
                var actual = Keys(affine
                    ? SequenceAligner.LocalAlignAffineTopNLinearSpace(a, b, sc, k, min, baseCells)
                    : SequenceAligner.LocalAlignTopNLinearSpace(a, b, sc, k, min, baseCells));
                pairs++;
                if (!expected.SequenceEqual(actual))
                    failures.Add($"{(affine ? "affine" : "linear")} {name} k={k} min={min} base={baseCells} {a} {b}\n  expected {string.Join(';', expected)}\n  actual   {string.Join(';', actual)}");
            }
        }
        Assert.That(pairs, Is.GreaterThanOrEqualTo(300));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {pairs}: " + string.Join("\n", failures.Take(3)));
    }

    [Test]
    public void LinearSpace_MatrixScoring_EqualTheInMemoryDeclumper_OnRandomPairs()
    {
        var rng = new Random(77);
        var protein = new (MatrixScoring Sc, string Name)[]
        {
            (new MatrixScoring(SubstitutionMatrix.Blosum62, -4, -4), "BLOSUM62 linear -4"),
            (new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -0.5), "BLOSUM62 -10/-0.5"),
            (new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -1), "BLOSUM62 -10/-1"),
        };
        var dna = new MatrixScoring(SubstitutionMatrix.EdnaFull, -10, -0.5);
        int pairs = 0;
        var failures = new List<string>();
        for (int iter = 0; iter < 80; iter++)
        {
            bool isProtein = iter % 2 == 0;
            var (a, b) = Pair(rng, isProtein ? "ACDEFGHIKLMNPQRSTVWY" : "ACGTN", 40);
            var sc = isProtein ? protein[iter / 2 % protein.Length].Sc : dna;
            int k = rng.Next(1, 9);
            double min = iter % 5 == 0 ? rng.Next(0, 12) / 2.0 : 0;
            long baseCells = BaseCellSizes[iter % BaseCellSizes.Length];
            foreach (bool affine in new[] { false, true })
            {
                var expected = Keys(affine ? SequenceAligner.LocalAlignAffineTopN(a, b, sc, k, min) : SequenceAligner.LocalAlignTopN(a, b, sc, k, min));
                var actual = Keys(affine
                    ? SequenceAligner.LocalAlignAffineTopNLinearSpace(a, b, sc, k, min, baseCells)
                    : SequenceAligner.LocalAlignTopNLinearSpace(a, b, sc, k, min, baseCells));
                pairs++;
                if (!expected.SequenceEqual(actual))
                    failures.Add($"{(affine ? "affine" : "linear")} {(isProtein ? "protein" : "dna")} k={k} min={min} base={baseCells} {a} {b}\n  expected {string.Join(';', expected)}\n  actual   {string.Join(';', actual)}");
            }
        }
        Assert.That(pairs, Is.GreaterThanOrEqualTo(160));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {pairs}: " + string.Join("\n", failures.Take(3)));
    }

    [Test]
    public void LinearSpace_ProteinAndDnaMatrixOverloads_AgreeWithTheStringOverloads()
    {
        var ms = new MatrixScoring(SubstitutionMatrix.Blosum62, -10, -1);
        const string a = "MKVLAAGIVGLLLAQPSMKVLAAGIVG", b = "MKVLAAGIVGPPPPPMKVLAGIVGLL";
        var viaString = Keys(SequenceAligner.LocalAlignAffineTopNLinearSpace(a, b, ms, 6, 0, 16));
        Assert.That(viaString, Is.Not.Empty);
        Assert.That(Keys(SequenceAligner.LocalAlignAffineTopNLinearSpace(new ProteinSequence(a), new ProteinSequence(b), ms, 6, 0, 16)), Is.EqualTo(viaString));
        Assert.That(Keys(SequenceAligner.LocalAlignAffineTopN(a, b, ms, 6)), Is.EqualTo(viaString));
        var dms = new MatrixScoring(SubstitutionMatrix.EdnaFull, -10, -1);
        const string d1 = "ACGTACGTTTACGTACGGACGTAC", d2 = "ACGTACGTAAGTACGGACGTAC";
        var dnaVia = Keys(SequenceAligner.LocalAlignTopNLinearSpace(d1, d2, dms, 5, 0, 16));
        Assert.That(Keys(SequenceAligner.LocalAlignTopNLinearSpace(new DnaSequence(d1), new DnaSequence(d2), dms, 5, 0, 16)), Is.EqualTo(dnaVia));
        Assert.That(Keys(SequenceAligner.LocalAlignTopNLinearSpace(new DnaSequence(d1), new DnaSequence(d2), new ScoringMatrix(2, -3, 0, -2), 5, 0, 16)),
            Is.EqualTo(Keys(SequenceAligner.LocalAlignTopN(d1, d2, new ScoringMatrix(2, -3, 0, -2), 5))));
        Assert.That(Keys(SequenceAligner.LocalAlignAffineTopNLinearSpace(new DnaSequence(d1), new DnaSequence(d2), new ScoringMatrix(2, -3, -4, -2), 5, 0, 16)),
            Is.EqualTo(Keys(SequenceAligner.LocalAlignAffineTopN(d1, d2, new ScoringMatrix(2, -3, -4, -2), 5))));
        Assert.That(Keys(SequenceAligner.LocalAlignAffineTopNLinearSpace(new DnaSequence(d1), new DnaSequence(d2), dms, 5, 0, 16)),
            Is.EqualTo(Keys(SequenceAligner.LocalAlignAffineTopN(d1, d2, dms, 5))));
        Assert.That(Keys(SequenceAligner.LocalAlignTopNLinearSpace(new ProteinSequence(a), new ProteinSequence(b), ms, 4, 0, 16)),
            Is.EqualTo(Keys(SequenceAligner.LocalAlignTopN(a, b, ms, 4))));
    }

    // ----------------------------------------------------------------------- the automatic switch (threshold lowered)

    [Test]
    public void AutomaticSwitch_AboveTheCellLimit_GivesTheInMemoryResult()
    {
        var rng = new Random(5);
        var sc = new ScoringMatrix(2, -3, -5, -2);
        for (int iter = 0; iter < 40; iter++)
        {
            var (a, b) = Pair(rng, "ACGT", 40);
            var ps = new PairScoring(sc);
            foreach (bool affine in new[] { false, true })
            {
                var memory = Keys(SequenceAligner.LocalTopNCore(a, b, ps, affine, 8, 0));
                // fullMatrixCells = 10 forces the linear path for every non-trivial input; baseCells = 30 forces the divide and conquer.
                var forced = Keys(SequenceAligner.LocalTopNCore(a, b, ps, affine, 8, 0, fullMatrixCells: 10, baseCells: 30));
                Assert.That(forced, Is.EqualTo(memory), $"{(affine ? "affine" : "linear")} {a} {b}");
            }
        }
    }

    [Test]
    public void LinearSpace_PositiveGapScore_IsRejectedAboveTheLimitWithAClearError()
    {
        var ps = new PairScoring(new ScoringMatrix(2, -3, 1, 1));
        var ex = Assert.Throws<ArgumentException>(() => SequenceAligner.LocalTopNCore("ACGTACGT", "ACGTACGT", ps, affine: true, 3, 0, fullMatrixCells: 10));
        Assert.That(ex!.Message, Does.Contain("non-positive gap scores"));
        Assert.That(SequenceAligner.LocalAlignAffineTopN("ACGTACGT", "ACGTACGT", new ScoringMatrix(2, -3, 1, 1), 3), Is.Not.Empty);
    }

    [Test]
    public void LinearSpace_EmptyAndNoPositivePair_GiveEmptyLists_AndArgumentsAreValidated()
    {
        Assert.That(SequenceAligner.LocalAlignTopNLinearSpace("", "ACGT"), Is.Empty);
        Assert.That(SequenceAligner.LocalAlignAffineTopNLinearSpace("AAAA", "CCCC"), Is.Empty);
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignTopNLinearSpace((string)null!, "A"));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignAffineTopNLinearSpace("A", (string)null!));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignTopNLinearSpace((DnaSequence)null!, new DnaSequence("A")));
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.LocalAlignTopNLinearSpace("A", "A", (MatrixScoring)null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.LocalAlignTopNLinearSpace("A", "A", (ScoringMatrix?)null, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAligner.LocalAlignTopNLinearSpace("A", "A", new MatrixScoring(SubstitutionMatrix.EdnaFull, -10, -1), 1, double.NaN));
    }

    // ----------------------------------------------------------------------- the independent oracle (matcher / full-recompute reference)

    [Test]
    public void LinearSpace_EqualsTheIndependentOracle_OnEveryRow()
    {
        // TestData/Alignment/topn_oracle.tsv: independent full-recompute Waterman-Eggert reference (column topN), EMBOSS matcher 6.6.0 (column emboss).
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.topn_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)!;
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        int rows = 0, matcherRows = 0, matcherEqual = 0;
        var failures = new List<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#') continue;
            var f = line.Split('\t');
            bool affine = f[1] == "affine";
            double go = double.Parse(f[5], CultureInfo.InvariantCulture), ge = double.Parse(f[6], CultureInfo.InvariantCulture);
            int k = int.Parse(f[9], CultureInfo.InvariantCulture);
            long baseCells = BaseCellSizes[rows % BaseCellSizes.Length];
            IReadOnlyList<AlignmentResult> actual;
            if (f[0] == "scalar")
            {
                var sc = new ScoringMatrix(int.Parse(f[3], CultureInfo.InvariantCulture), int.Parse(f[4], CultureInfo.InvariantCulture), affine ? (int)go : 0, (int)ge);
                actual = affine ? SequenceAligner.LocalAlignAffineTopNLinearSpace(f[7], f[8], sc, k, 0, baseCells) : SequenceAligner.LocalAlignTopNLinearSpace(f[7], f[8], sc, k, 0, baseCells);
            }
            else
            {
                var ms = new MatrixScoring(SubstitutionMatrix.Blosum62, go, ge);
                actual = affine ? SequenceAligner.LocalAlignAffineTopNLinearSpace(f[7], f[8], ms, k, 0, baseCells) : SequenceAligner.LocalAlignTopNLinearSpace(f[7], f[8], ms, k, 0, baseCells);
            }

            var expected = f[11].Length == 0 ? [] : f[11].Split(';');
            var got = actual.Select(r => string.Create(CultureInfo.InvariantCulture, $"{r.ExactScore:R}|{r.StartPosition1}|{r.StartPosition2}|{r.EndPosition1}|{r.EndPosition2}|{r.AlignedSequence1}|{r.AlignedSequence2}")).ToArray();
            expected = expected.Select(e => { int bar = e.IndexOf('|'); return double.Parse(e[..bar], CultureInfo.InvariantCulture).ToString("R", CultureInfo.InvariantCulture) + e[bar..]; }).ToArray();
            if (!expected.SequenceEqual(got))
                failures.Add($"{f[0]}/{f[1]} {f[7]} {f[8]} k={k}\n  expected {string.Join(';', expected)}\n  actual   {string.Join(';', got)}");
            if (f[12] != "NA")
            {
                var em = f[12].Split(',').Select(s => double.Parse(s, CultureInfo.InvariantCulture)).ToArray();
                matcherRows++;
                if (actual.Select(r => r.ExactScore).Take(em.Length).SequenceEqual(em) && actual.Count >= em.Length) matcherEqual++;
            }
            rows++;
        }
        Assert.That(rows, Is.GreaterThanOrEqualTo(1000));
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}: " + string.Join("\n", failures.Take(3)));
        Assert.That(matcherRows, Is.GreaterThanOrEqualTo(300));
        Assert.That(matcherEqual, Is.GreaterThanOrEqualTo((int)(matcherRows * 0.9)), $"{matcherEqual} of {matcherRows} matcher score lists equal");
    }

    // ----------------------------------------------------------------------- large input (above TopNMaxCells), automatic path

    private static (string A, string B) LargePair(int length, int seed)
    {
        var rng = new Random(seed);
        var a = Random(rng, "ACGT", length).ToCharArray();
        var b = Random(rng, "ACGT", length).ToCharArray();
        // Planted similar segments (several non-intersecting alignments), with substitutions and indels.
        for (int s = 0; s < 6; s++)
        {
            int len = 150 + 60 * s;
            int pa = rng.Next(0, length - len), pb = rng.Next(0, length - len);
            for (int i = 0; i < len; i++)
                b[pb + i] = rng.Next(12) == 0 ? "ACGT"[rng.Next(4)] : a[pa + i];
        }
        return (new string(a), new string(b));
    }

    [Test]
    public void LargeInput_AboveTheCellLimit_LinearModel_EqualsTheInMemoryDeclumper_AndStaysInLinearMemory()
    {
        var (a, b) = LargePair(5300, 1);   // (5301)^2 = 28.1 M cells > TopNMaxCells
        Assert.That((long)(a.Length + 1) * (b.Length + 1), Is.GreaterThan(SequenceAligner.TopNMaxCells));
        var sc = new ScoringMatrix(2, -3, 0, -4);

        long before = GC.GetAllocatedBytesForCurrentThread();
        var linear = SequenceAligner.LocalAlignTopN(a, b, sc, 4);
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;

        Assert.That(linear, Has.Count.EqualTo(4));
        var memory = SequenceAligner.LocalTopNCore(a, b, new PairScoring(sc), affine: false, 4, 0, fullMatrixCells: long.MaxValue);
        Assert.That(Keys(linear), Is.EqualTo(Keys(memory)));
        Assert.That(Key(linear[0]), Is.EqualTo(Key(SequenceAligner.LocalAlign(a, b, sc))), "alignment 1 is LocalAlign");

        var used = new HashSet<(int, int)>();
        foreach (var r in linear)
            foreach (var p in SequenceAligner.AlignedPairs(r))
                Assert.That(used.Add(p), Is.True, "no aligned pair is reported twice");
        // The matrix of this input needs 113 MB (28.1 M cells x 4 B); the linear program allocates far less in total.
        TestContext.Out.WriteLine($"allocated by the linear-space top-4: {allocated / 1_000_000} MB (matrix would be {(long)(a.Length + 1) * (b.Length + 1) * 4 / 1_000_000} MB)");
        Assert.That(allocated, Is.LessThan((long)(a.Length + 1) * (b.Length + 1) * 4));
    }

    [Test]
    public void LargeInput_AboveTheCellLimit_AffineModel_EqualsTheInMemoryDeclumper()
    {
        var (a, b) = LargePair(5100, 2);   // (5101)^2 = 26.0 M cells
        var sc = new ScoringMatrix(2, -3, -5, -2);
        var linear = SequenceAligner.LocalAlignAffineTopN(a, b, sc, 3);
        var memory = SequenceAligner.LocalTopNCore(a, b, new PairScoring(sc), affine: true, 3, 0, fullMatrixCells: long.MaxValue);
        Assert.That(linear, Has.Count.EqualTo(3));
        Assert.That(Keys(linear), Is.EqualTo(Keys(memory)));
        Assert.That(Key(linear[0]), Is.EqualTo(Key(SequenceAligner.LocalAlignAffine(a, b, sc))), "alignment 1 is LocalAlignAffine");
    }

    [Test]
    public void CoOptimalEnumeration_StillRejectsInputsAboveTheCellLimit()
    {
        string big = new('A', 5001);
        Assert.Throws<ArgumentException>(() => SequenceAligner.LocalAlignAffineCoOptimal(big, big));
    }
}
