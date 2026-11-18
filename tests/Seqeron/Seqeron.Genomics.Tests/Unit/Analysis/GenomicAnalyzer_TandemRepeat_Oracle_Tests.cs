// GENOMIC-TANDEM-001 (B09 F35-F38) — GenomicAnalyzer.FindTandemRepeats against independent oracles.
// Evidence: docs/Evidence/GENOMIC-TANDEM-001-Evidence.md, TestSpec: tests/TestSpecs/GENOMIC-TANDEM-001.md
//
// Definition under test (Kolpakov & Kucherov 1999 maximal repetitions; Gusfield 1997 sec. 9; Benson 1999 "two or more
// contiguous copies"; microsatellite = primitive unit repeated, MISA "reject false type motifs"): for unit length p,
// each maximal interval with period p and >= minReps complete copies whose unit S[a..a+p) is primitive, reported once
// at its left end a with floor((e-a)/p) copies; order = unit length, then start.
//
// Oracle 1: Python, 1156 (sequence, minUnit, minReps, maxUnit) rows, three independent methods asserted equal
//           (brute-force Kolpakov-Kucherov runs, boolean-array runs, regex) — scripts/analysis/generate_tandem_oracle.py.
// Oracle 2: in-test brute force (every interval, its smallest period by definition) on random strings.

using System.Diagnostics;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_TandemRepeat_Oracle_Tests
{
    private sealed record Row(string Sequence, int MinUnit, int MinReps, int MaxUnit, int Count, string Hits);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.tandem_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], int.Parse(f[1]), int.Parse(f[2]), int.Parse(f[3]), int.Parse(f[4]), f[5]));
        }
        return rows;
    }

    private static string Format(IEnumerable<TandemRepeat> items)
    {
        string s = string.Join(";", items.Select(t => $"{t.Unit}:{t.Position}:{t.Repetitions}"));
        return s.Length == 0 ? "-" : s;
    }

    private static IEnumerable<TandemRepeat> Run(string seq, int minUnit, int minReps, int maxUnit = 0)
        => maxUnit == 0
            ? GenomicAnalyzer.FindTandemRepeats(new DnaSequence(seq), minUnit, minReps)
            : GenomicAnalyzer.FindTandemRepeats(new DnaSequence(seq), minUnit, minReps, maxUnit);

    #region Python oracle

    [Test]
    public void FindTandemRepeats_MatchesPythonOracle_AllRows()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(1156));
        int nonEmpty = 0;
        foreach (var row in rows)
        {
            var actual = Run(row.Sequence, row.MinUnit, row.MinReps, row.MaxUnit).ToList();
            Assert.That(actual, Has.Count.EqualTo(row.Count), $"{row.Sequence} min={row.MinUnit} reps={row.MinReps} max={row.MaxUnit}");
            Assert.That(Format(actual), Is.EqualTo(row.Hits), $"{row.Sequence} min={row.MinUnit} reps={row.MinReps} max={row.MaxUnit}");
            if (actual.Count > 0) nonEmpty++;
        }
        Assert.That(nonEmpty, Is.GreaterThan(500), "the corpus must exercise non-trivial results");
    }

    [Test]
    public void FindTandemRepeats_AgreesWithRepeatFinderMicrosatellites_OnOracleCorpus()
    {
        // The facade must be the canonical engine's result mapped 1:1 (no re-implementation).
        foreach (var row in ReadOracle())
        {
            int max = row.MaxUnit == 0 ? int.MaxValue : row.MaxUnit;
            var ga = Run(row.Sequence, row.MinUnit, row.MinReps, row.MaxUnit).Select(t => (t.Unit, t.Position, t.Repetitions));
            var rf = RepeatFinder.FindMicrosatellites(new DnaSequence(row.Sequence), row.MinUnit, max, row.MinReps)
                .Select(m => (m.RepeatUnit, m.Position, m.RepeatCount));
            Assert.That(ga, Is.EqualTo(rf), $"{row.Sequence} min={row.MinUnit} reps={row.MinReps} max={row.MaxUnit}");
        }
    }

    #endregion

    #region In-test brute force (Kolpakov-Kucherov runs by definition)

    private static bool IsPrimitive(string u) => (u + u).IndexOf(u, 1, StringComparison.Ordinal) == u.Length;

    private static List<(string Unit, int Position, int Copies)> BruteForce(string s, int minUnit, int minReps, int maxUnit)
    {
        int n = s.Length;
        var found = new HashSet<(int p, int a, int c)>();
        for (int a = 0; a < n; a++)
        for (int e = a + 2; e <= n; e++)
        {
            int m = e - a;
            int q = 1;
            while (q <= m)
            {
                bool ok = true;
                for (int x = 0; x + q < m && ok; x++)
                    ok = s[a + x] == s[a + x + q];
                if (ok) break;
                q++;
            }

            if (2 * q > m) continue;
            if (a > 0 && s[a - 1] == s[a - 1 + q]) continue; // extendable to the left
            if (e < n && s[e] == s[e - q]) continue;          // extendable to the right
            int copies = m / q;
            if (q >= minUnit && q <= maxUnit && copies >= minReps && IsPrimitive(s.Substring(a, q)))
                found.Add((q, a, copies));
        }

        return found.OrderBy(t => t.p).ThenBy(t => t.a).Select(t => (s.Substring(t.a, t.p), t.a, t.c)).ToList();
    }

    [Test]
    public void FindTandemRepeats_MatchesInTestBruteForce_RandomStrings()
    {
        var rnd = new Random(20261118);
        string[] alphabets = ["AC", "ACG", "ACGT", "AT"];
        for (int iter = 0; iter < 1500; iter++)
        {
            string alpha = alphabets[rnd.Next(alphabets.Length)];
            int n = rnd.Next(0, 40);
            var sb = new System.Text.StringBuilder();
            for (int i = 0; i < n; i++) sb.Append(alpha[rnd.Next(alpha.Length)]);
            if (n >= 6 && rnd.Next(3) > 0)
            {
                var unit = new System.Text.StringBuilder();
                for (int i = rnd.Next(1, 5); i > 0; i--) unit.Append(alpha[rnd.Next(alpha.Length)]);
                string block = string.Concat(Enumerable.Repeat(unit.ToString(), rnd.Next(2, 6)));
                sb.Insert(rnd.Next(sb.Length + 1), block);
            }

            string s = sb.ToString();
            int minUnit = rnd.Next(1, 4), minReps = rnd.Next(2, 5);
            int maxUnit = rnd.Next(minUnit, 8);
            var expected = BruteForce(s, minUnit, minReps, maxUnit);
            var actual = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(s), minUnit, minReps, maxUnit)
                .Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
            Assert.That(actual, Is.EqualTo(expected), $"{s} min={minUnit} reps={minReps} max={maxUnit}");

            var unbounded = BruteForce(s, minUnit, minReps, int.MaxValue);
            var actualUnbounded = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(s), minUnit, minReps)
                .Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
            Assert.That(actualUnbounded, Is.EqualTo(unbounded), $"{s} min={minUnit} reps={minReps} unbounded");
        }
    }

    #endregion

    #region Run-enumeration engine (Lyndon runs, Bannai et al. 2017) vs per-unit-length engine

    [Test]
    public void RunsEngine_MatchesPythonOracle_AllRows()
    {
        foreach (var row in ReadOracle())
        {
            int max = row.MaxUnit == 0 ? int.MaxValue : row.MaxUnit;
            var actual = GenomicAnalyzer.FindTandemRepeatsByRuns(new DnaSequence(row.Sequence), row.MinUnit, row.MinReps, max).ToList();
            Assert.That(Format(actual), Is.EqualTo(row.Hits), $"{row.Sequence} min={row.MinUnit} reps={row.MinReps} max={row.MaxUnit}");
        }
    }

    [Test]
    public void RunsEngine_MatchesBruteForceRuns_RandomAndStructuredStrings()
    {
        var rnd = new Random(77);
        var texts = new List<string>();
        for (int iter = 0; iter < 1200; iter++)
        {
            string alpha = new[] { "AC", "ACG", "ACGT", "AT" }[rnd.Next(4)];
            int n = rnd.Next(0, 60);
            texts.Add(new string(Enumerable.Range(0, n).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray()));
        }

        // highly periodic families: Fibonacci / Thue-Morse words, powers, nested repeats
        string fa = "A", fb = "AC";
        for (int k = 0; k < 7; k++) { (fa, fb) = (fb, fb + fa); texts.Add(fb[..Math.Min(fb.Length, 70)]); }
        string tm = "A";
        for (int k = 0; k < 6; k++) tm += new string(tm.Select(c => c == 'A' ? 'C' : 'A').ToArray());
        texts.Add(tm[..60]);
        texts.Add(string.Concat(Enumerable.Repeat("ACGTACGTAC", 6)));
        texts.Add(string.Concat(Enumerable.Repeat("AAAAAAAAAAAAAAAAAAAA", 2)));
        texts.Add(string.Concat(Enumerable.Repeat("ACACGACACGACACT", 4)));

        foreach (string s in texts)
        {
            var expected = BruteForce(s, 1, 2, int.MaxValue);
            var actual = GenomicAnalyzer.FindTandemRepeatsByRuns(new DnaSequence(s), 1, 2, int.MaxValue)
                .Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
            Assert.That(actual, Is.EqualTo(expected), s);

            var viaScan = RepeatFinder.FindMicrosatellites(new DnaSequence(s), 1, int.MaxValue, 2)
                .Select(m => (m.RepeatUnit, m.Position, m.RepeatCount)).ToList();
            Assert.That(actual, Is.EqualTo(viaScan), s);
        }
    }

    [Test]
    public void RunsEngine_ExhaustiveSmallStrings_AgreeWithBruteForceRuns()
    {
        // Every string over {A,C} up to length 14 and over {A,C,G} up to length 9 (completeness of the Lyndon-root candidates).
        int checkedStrings = 0;
        foreach (var (alphabet, maxLen) in new[] { ("AC", 14), ("ACG", 9) })
        {
            for (int len = 2; len <= maxLen; len++)
            {
                long total = (long)Math.Pow(alphabet.Length, len);
                var chars = new char[len];
                for (long code = 0; code < total; code++)
                {
                    long c = code;
                    for (int i = 0; i < len; i++) { chars[i] = alphabet[(int)(c % alphabet.Length)]; c /= alphabet.Length; }
                    string s = new(chars);
                    var runs = MaximalRepetitions.FindRuns(s);
                    var expected = BruteForce(s, 1, 2, int.MaxValue);
                    var actual = runs.Select(r => (s.Substring(r.Start, r.Period), r.Start, (r.End - r.Start) / r.Period)).ToList();
                    if (!actual.SequenceEqual(expected))
                        Assert.Fail($"{s}: runs {string.Join(";", actual)} != brute force {string.Join(";", expected)}");
                    checkedStrings++;
                }
            }
        }

        Assert.That(checkedStrings, Is.GreaterThan(60_000));
    }

    [Test]
    public void RunsEngine_AgreesWithPerUnitLengthScan_OnLongLowEntropyStrings()
    {
        var rnd = new Random(4242);
        for (int iter = 0; iter < 60; iter++)
        {
            string alpha = new[] { "AC", "ACG", "AT" }[iter % 3];
            var sb = new System.Text.StringBuilder();
            while (sb.Length < 2500)
            {
                // mixture of random stretches and planted repeats of unit length 1..40
                if (rnd.Next(3) == 0)
                    for (int i = rnd.Next(1, 30); i > 0; i--) sb.Append(alpha[rnd.Next(alpha.Length)]);
                else
                {
                    var unit = new System.Text.StringBuilder();
                    for (int i = rnd.Next(1, 40); i > 0; i--) unit.Append(alpha[rnd.Next(alpha.Length)]);
                    for (int c = rnd.Next(2, 7); c > 0; c--) sb.Append(unit);
                }
            }

            var dna = new DnaSequence(sb.ToString());
            int minReps = rnd.Next(2, 4);
            var viaRuns = GenomicAnalyzer.FindTandemRepeatsByRuns(dna, 1, minReps, int.MaxValue)
                .Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
            var viaScan = RepeatFinder.FindMicrosatellites(dna, 1, int.MaxValue, minReps)
                .Select(m => (m.RepeatUnit, m.Position, m.RepeatCount)).ToList();
            Assert.That(viaRuns, Is.EqualTo(viaScan), $"iteration {iter} reps={minReps}");
            Assert.That(viaRuns, Is.Not.Empty);
        }
    }

    [Test]
    public void WideUnitRange_SwitchesToRunsEngine_WithIdenticalResult()
    {
        // 700 bp, minUnit 1 .. unbounded -> > RunsEngineMinUnitLengths admissible unit lengths -> run enumeration.
        string unit = "ACGTTGCAAGCTTAGG"; // 16 bp
        string seq = RandomDna(200, 5) + string.Concat(Enumerable.Repeat(unit, 6)) + RandomDna(300, 6);
        var dna = new DnaSequence(seq);
        Assert.That(seq.Length / 2 - 1 + 1, Is.GreaterThan(GenomicAnalyzer.RunsEngineMinUnitLengths));

        var viaRuns = GenomicAnalyzer.FindTandemRepeats(dna, 1, 2).Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
        var viaScan = RepeatFinder.FindMicrosatellites(dna, 1, int.MaxValue, 2).Select(m => (m.RepeatUnit, m.Position, m.RepeatCount)).ToList();
        Assert.That(viaRuns, Is.EqualTo(viaScan));
        Assert.That(viaRuns, Does.Contain((unit, 200, 6)));
    }

    #endregion

    #region Contract (primitive units, maximal runs, order) — B04/DUP_MAP §17 request

    [Test]
    public void NonPrimitiveUnit_IsNotReportedAlongsideItsPrimitiveRoot()
    {
        // Defect repro: ATATATAT, minUnit 1, minReps 2 returned AT x4 @0 AND ATAT x2 @0 (a non-primitive unit).
        var hits = Run("ATATATAT", 1, 2).ToList();
        Assert.That(Format(hits), Is.EqualTo("AT:0:4"));

        // AAAAAA: only the primitive A x6 — neither AA x3 nor AAA x2.
        Assert.That(Format(Run("AAAAAA", 1, 2)), Is.EqualTo("A:0:6"));

        // minUnit above the primitive period: a non-primitive multiple is not a microsatellite/tandem unit either.
        Assert.That(Run("ATATATAT", 3, 2), Is.Empty);
        Assert.That(Run("ATATATAT", 4, 2), Is.Empty);
    }

    [Test]
    public void MaximalRun_IsReportedOnce_AtItsLeftEnd_WithCompleteCopiesOnly()
    {
        // ATATATA: run [0,7) period 2 -> AT x3 (the trailing A is a partial copy); TA @1 is a rotation, not reported.
        Assert.That(Format(Run("ATATATA", 2, 2)), Is.EqualTo("AT:0:3"));
        // Flanks do not change the run: T + (CAG x4) + T (left/right extension fails at both ends).
        Assert.That(Format(Run("TCAGCAGCAGCAGT", 3, 2)), Is.EqualTo("CAG:1:4"));
        // The unit is the motif phase at the run's LEFT end: GCAGCAGCAGCAGG starts the period-3 run at the G
        // (GCA x4 @0; MISA / pytrf convention), not at CAG @1.
        Assert.That(Format(Run("GCAGCAGCAGCAGG", 3, 2)), Is.EqualTo("GCA:0:4"));
        Assert.That(Format(Run("GACACACT", 2, 2)), Is.EqualTo("AC:1:3"));
    }

    [Test]
    public void Result_IsOrderedByUnitLengthThenPosition()
    {
        // A x3 @0 and @9 (unit 1) precede the unit-2 run AC x3 @2 (the run starts at the A of AAA|CACAC).
        var hits = Run("AAACACACGAAA", 1, 3).ToList();
        Assert.That(Format(hits), Is.EqualTo("A:0:3;A:9:3;AC:2:3"));
        Assert.That(hits.Select(h => h.Unit.Length), Is.Ordered);
    }

    [Test]
    public void RunsOfDifferentUnitLengths_OverlapAndAreReportedIndependently()
    {
        // AAAGAAAGAAAG: AAAG x3 @0 plus the three A x3 homopolymer runs @0, @4, @8.
        Assert.That(Format(Run("AAAGAAAGAAAG", 1, 2)), Is.EqualTo("A:0:3;A:4:3;A:8:3;AAAG:0:3"));
    }

    [Test]
    public void MaxUnitLength_RestrictsUnitLength_AndValidates()
    {
        // Oracle (scripts/analysis/generate_tandem_oracle.py): ACG x4 @0 (p=3), GT x2 @11 (p=2), TTAGC x3 @14 (p=5),
        // T x2 @14/@19/@24 (the TT of TTAGC; p=1).
        string seq = "ACGACGACGACGTGTTAGCTTAGCTTAGC";
        Assert.That(Format(Run(seq, 1, 2)), Is.EqualTo("T:14:2;T:19:2;T:24:2;GT:11:2;ACG:0:4;TTAGC:14:3"));
        Assert.That(Format(Run(seq, 1, 2, 3)), Is.EqualTo("T:14:2;T:19:2;T:24:2;GT:11:2;ACG:0:4"));
        Assert.That(Format(Run(seq, 1, 2, 1)), Is.EqualTo("T:14:2;T:19:2;T:24:2"));
        Assert.That(Format(Run(seq, 5, 2, 5)), Is.EqualTo("TTAGC:14:3"));
        Assert.That(Run(seq, 4, 2, 4), Is.Empty);

        var dna = new DnaSequence(seq);
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindTandemRepeats(dna, 3, 2, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindTandemRepeats(dna, 0, 2, 5));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindTandemRepeats(dna, 1, 1, 5));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindTandemRepeats(null!, 1, 2, 5));
        // eager validation (the exception surfaces at the call, not on enumeration)
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindTandemRepeats(dna, 5, 2, 4));
    }

    [Test]
    public void MaxUnitLength_LargerThanSequence_IsHarmless_EmptyAndShortSequences()
    {
        Assert.That(Run("ACACAC", 1, 2, int.MaxValue), Is.Not.Empty);
        Assert.That(Format(Run("ACACAC", 1, 2, int.MaxValue)), Is.EqualTo("AC:0:3"));
        Assert.That(Run("", 2, 2, 5), Is.Empty);
        Assert.That(Run("", 2, 2), Is.Empty);
        Assert.That(Run("A", 1, 2, 5), Is.Empty);
        Assert.That(Run("AC", 5, 2, 9), Is.Empty);
    }

    [Test]
    public void LowerCaseInput_IsNormalised_ByDnaSequence()
    {
        Assert.That(Format(Run("atgatgatg", 3, 3)), Is.EqualTo("ATG:0:3"));
    }

    #endregion

    #region Complexity

    private static string RandomDna(int n, int seed)
    {
        var rnd = new Random(seed);
        var chars = new char[n];
        for (int i = 0; i < n; i++) chars[i] = "ACGT"[rnd.Next(4)];
        return new string(chars);
    }

    [Test]
    public void OneMegabase_BoundedUnitLength_RunsInLinearTimePerUnitLength()
    {
        // 1 Mbp random sequence with planted repeats; the old detector was O(n^2 * m) substring scans (hours).
        var chars = RandomDna(1_000_000, 1).ToCharArray();
        string planted = string.Concat(Enumerable.Repeat("CAG", 40));
        for (int i = 0; i < planted.Length; i++) chars[500_000 + i] = planted[i];
        var dna = new DnaSequence(new string(chars));

        var sw = Stopwatch.StartNew();
        var hits = GenomicAnalyzer.FindTandemRepeats(dna, 2, 5, 60).ToList();
        sw.Stop();

        var cag = hits.Single(h => h.Unit == "CAG" && h.Position == 500_000);
        Assert.That(cag.Repetitions, Is.EqualTo(40));
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(60), $"took {sw.Elapsed}");
    }


    [Test]
    public void OneMegabase_HighlyRepetitive_Unbounded_StillFast()
    {
        // Worst case for the suffix-array doubling (every prefix doubling round is needed): a 1 Mbp homopolymer and a
        // 1 Mbp period-8 string. The only run of the homopolymer is A x 1 000 000.
        var poly = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(new string('A', 1_000_000))).ToList();
        Assert.That(poly, Has.Count.EqualTo(0), "default minUnitLength 2: the primitive unit is A (length 1)");
        var poly1 = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(new string('A', 1_000_000)), 1, 2).ToList();
        Assert.That(Format(poly1), Is.EqualTo("A:0:1000000"));

        var sw = Stopwatch.StartNew();
        var periodic = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(string.Concat(Enumerable.Repeat("ACGTTGCA", 125_000)))).ToList();
        sw.Stop();
        Assert.That(periodic.Single(h => h.Unit.Length == 8).Repetitions, Is.EqualTo(125_000));
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(120), $"took {sw.Elapsed}");
    }

    [Test]
    public void OneMegabase_UnboundedUnitLength_UsesRunEnumeration()
    {
        // Unbounded unit length on 1 Mbp: a per-unit-length scan would need up to 500 000 passes over the sequence
        // (about 5e11 comparisons); the run enumeration is O(n log^2 n).
        var chars = RandomDna(1_000_000, 2).ToCharArray();
        string planted = string.Concat(Enumerable.Repeat("ACGTTGCAAGCTTAGGCATTG", 30)); // 21 bp unit x 30
        for (int i = 0; i < planted.Length; i++) chars[300_000 + i] = planted[i];
        var dna = new DnaSequence(new string(chars));

        var sw = Stopwatch.StartNew();
        var hits = GenomicAnalyzer.FindTandemRepeats(dna).ToList();
        sw.Stop();

        var big = hits.Single(h => h.Unit == "ACGTTGCAAGCTTAGGCATTG" && h.Position == 300_000);
        Assert.That(big.Repetitions, Is.EqualTo(30));
        Assert.That(hits.Select(h => h.Unit.Length), Is.Ordered);
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(120), $"took {sw.Elapsed}");
    }
    #endregion
}
