// GENOMIC-MOTIFS-001 (B09 A1-20) — mismatch-tolerant known-motif search: GenomicAnalyzer.FindKnownMotifHits(maxMismatches)
// and ApproximateMatcher.FindWithMismatches(iupacPattern: true).
// Evidence: docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md; TestSpec: tests/TestSpecs/GENOMIC-MOTIFS-001.md
// Oracle: EMBOSS fuzznuc 6.6.0 (-pmismatch 1..3 -complement Y, hit positions and the "*mismatch n" notes) == independent Python
//         Hamming scan, 6576 (sequence, pattern, k) rows, scripts/analysis/generate_known_motif_oracle.py --mismatches
//         (TestData/Analysis/known_motif_mismatch_oracle.tsv). fuzznuc's overhanging hits on sequences shorter than the
//         pattern are not part of the oracle (documented quirk).

using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_KnownMotifMismatch_Tests
{
    private sealed record Row(string Sequence, string Pattern, int K, KnownMotifHit[] Forward, KnownMotifHit[] Reverse);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.known_motif_mismatch_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], int.Parse(f[2]), Parse(f[3], false), Parse(f[4], true)));
        }

        return rows;

        static KnownMotifHit[] Parse(string s, bool reverse) => s == "-"
            ? Array.Empty<KnownMotifHit>()
            : s.Split(',').Select(x => x.Split(':')).Select(x => new KnownMotifHit(int.Parse(x[0]), reverse, int.Parse(x[1]))).ToArray();
    }

    private static bool IsAcgt(string pattern) => pattern.All(c => "ACGT".Contains(c));

    #region fuzznuc oracle

    [Test]
    public void FindKnownMotifHits_Mismatches_MatchFuzznucOracle_PositionsStrandsAndCounts()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(6576));
        int withHits = 0, withMismatchingHits = 0;
        foreach (var row in rows)
        {
            var hits = GenomicAnalyzer.FindKnownMotifHits(new DnaSequence(row.Sequence), new[] { row.Pattern },
                bothStrands: true, allowIupac: true, maxMismatches: row.K);

            var expected = row.Forward.Concat(row.Reverse)
                .OrderBy(h => h.Position).ThenBy(h => h.IsReverseStrand).ToArray();
            if (expected.Length == 0)
            {
                Assert.That(hits, Is.Empty, $"{row.Sequence} / {row.Pattern} / k={row.K}");
                continue;
            }

            Assert.That(hits.Keys, Is.EqualTo(new[] { row.Pattern }));
            Assert.That(hits[row.Pattern], Is.EqualTo(expected), $"{row.Sequence} / {row.Pattern} / k={row.K}");
            withHits++;
            if (expected.Any(h => h.Mismatches > 0)) withMismatchingHits++;
        }

        Assert.That(withHits, Is.GreaterThan(4000));
        Assert.That(withMismatchingHits, Is.GreaterThan(4000));
    }

    [Test]
    public void FindKnownMotifHits_Mismatches_ForwardOnly_And_LiteralMode_AgreeWithOracleForAcgtPatterns()
    {
        int checkedRows = 0;
        foreach (var row in ReadOracle().Where(r => IsAcgt(r.Pattern)))
        {
            var dna = new DnaSequence(row.Sequence);
            // Literal mode (IUPAC off) and IUPAC mode coincide for a pattern without ambiguity codes.
            var literal = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { row.Pattern }, bothStrands: true, maxMismatches: row.K);
            var iupac = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { row.Pattern }, bothStrands: true, allowIupac: true, maxMismatches: row.K);
            Assert.That(literal.Keys, Is.EqualTo(iupac.Keys));
            foreach (var key in literal.Keys)
                Assert.That(literal[key], Is.EqualTo(iupac[key]), $"{row.Sequence} / {row.Pattern} / k={row.K}");

            // Forward strand only = the forward hits of the oracle.
            var fwdOnly = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { row.Pattern }, bothStrands: false, maxMismatches: row.K);
            Assert.That(fwdOnly.TryGetValue(row.Pattern, out var f) ? f : Array.Empty<KnownMotifHit>(),
                Is.EqualTo(row.Forward), $"{row.Sequence} / {row.Pattern} / k={row.K}");
            checkedRows++;
        }

        Assert.That(checkedRows, Is.GreaterThan(1500));
    }

    #endregion

    #region fuzznuc worked examples (positions 0-based, forward-strand window starts, "*mismatch n" notes)

    [Test]
    public void FuzznucWorkedExample_EcoRI_OneMismatch()
    {
        // fuzznuc -pattern GAATTC -pmismatch 1 -complement Y on GAATTCAAGATTCGGATTTTGAATTTT:
        // + 1-6 (exact), + 21-26 (*mismatch 1), - 1-6, - 21-26 (*mismatch 1)  [1-based]
        var dna = new DnaSequence("GAATTCAAGATTCGGATTTTGAATTTT");
        var hits = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "GAATTC" }, bothStrands: true, maxMismatches: 1)["GAATTC"];
        Assert.That(hits, Is.EqualTo(new[]
        {
            new KnownMotifHit(0, false, 0), new KnownMotifHit(0, true, 0),
            new KnownMotifHit(20, false, 1), new KnownMotifHit(20, true, 1),
        }));
    }

    [Test]
    public void FuzznucWorkedExample_Degenerate_TwoMismatches_BothStrands()
    {
        // fuzznuc -pattern RAATTC -pmismatch 2 -complement Y on the same sequence ('-' = windows matching GAATTY).
        var dna = new DnaSequence("GAATTCAAGATTCGGATTTTGAATTTT");
        var hits = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "RAATTC" }, bothStrands: true, allowIupac: true, maxMismatches: 2)["RAATTC"];
        var forward = hits.Where(h => !h.IsReverseStrand).Select(h => (h.Position, h.Mismatches)).ToArray();
        var reverse = hits.Where(h => h.IsReverseStrand).Select(h => (h.Position, h.Mismatches)).ToArray();
        Assert.That(forward, Is.EqualTo(new[] { (0, 0), (7, 1), (13, 2), (14, 2), (20, 1), (21, 2) }));
        Assert.That(reverse, Is.EqualTo(new[] { (0, 0), (7, 2), (13, 1), (14, 1), (20, 0), (21, 2) }));
    }

    [Test]
    public void FuzznucWorkedExample_NonPalindrome_ReverseHitHasItsOwnCount()
    {
        // fuzznuc -pattern AACCGT -pmismatch 1 -complement Y on AACCGTACGGTA: + 1-6 exact; - 7-12 (*mismatch 1): ACGGTA vs ACGGTT
        var dna = new DnaSequence("AACCGTACGGTA");
        var hits = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "AACCGT" }, bothStrands: true, maxMismatches: 1)["AACCGT"];
        Assert.That(hits, Is.EqualTo(new[] { new KnownMotifHit(0, false, 0), new KnownMotifHit(6, true, 1) }));
    }

    #endregion

    #region Contract

    [Test]
    public void ZeroMismatches_IsTheExactSearch_WithZeroCounts()
    {
        var dna = new DnaSequence("GAATTCAAGATTCGGATTTTGAATTTTTATAAAG");
        var exact = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "GAATTC", "TATAAA", "AA" }, bothStrands: true);
        var zero = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "GAATTC", "TATAAA", "AA" }, bothStrands: true, maxMismatches: 0);
        Assert.That(zero.Keys, Is.EqualTo(exact.Keys));
        foreach (var key in exact.Keys)
        {
            Assert.That(zero[key], Is.EqualTo(exact[key]));
            Assert.That(exact[key], Has.All.Matches<KnownMotifHit>(h => h.Mismatches == 0));
        }
    }

    [Test]
    public void IncreasingTheBound_OnlyAddsHits_AndKeepsTheCounts()
    {
        var rnd = new Random(2026);
        for (int trial = 0; trial < 60; trial++)
        {
            var dna = new DnaSequence(new string(Enumerable.Range(0, 120).Select(_ => "ACGT"[rnd.Next(4)]).ToArray()));
            string motif = new string(Enumerable.Range(0, 6).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var prev = new HashSet<KnownMotifHit>();
            for (int k = 0; k <= 3; k++)
            {
                var hits = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { motif }, bothStrands: true, maxMismatches: k);
                var current = hits.TryGetValue(motif, out var l) ? l.ToHashSet() : new HashSet<KnownMotifHit>();
                Assert.That(current.IsSupersetOf(prev), $"{motif} k={k}");
                Assert.That(current, Has.All.Matches<KnownMotifHit>(h => h.Mismatches <= k));
                prev = current;
            }
        }
    }

    [Test]
    public void LiteralMode_AmbiguityCodeIsOneMismatch_IupacModeAgreesWithItsBaseSet()
    {
        var dna = new DnaSequence("ACGTA");
        // windows ACG, CGT, GTA against ANG
        var literal1 = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "ANG" }, maxMismatches: 1)["ANG"];
        Assert.That(literal1, Is.EqualTo(new[] { new KnownMotifHit(0, false, 1) }), "N is a literal symbol: ACG has 1 mismatch (C vs N)");
        Assert.That(GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "ANG" }, maxMismatches: 0), Is.Empty);

        var iupac0 = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "ANG" }, allowIupac: true)["ANG"];
        Assert.That(iupac0, Is.EqualTo(new[] { new KnownMotifHit(0, false, 0) }));
        var iupac2 = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "ANG" }, allowIupac: true, maxMismatches: 2)["ANG"];
        Assert.That(iupac2, Is.EqualTo(new[] { new KnownMotifHit(0, false, 0), new KnownMotifHit(1, false, 2), new KnownMotifHit(2, false, 2) }));
    }

    [Test]
    public void BoundAtLeastTheMotifLength_MakesEveryWindowAHit_ShorterSequenceHasNone()
    {
        var dna = new DnaSequence("ACGTAC");
        var all = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "TTT" }, maxMismatches: 3)["TTT"];
        Assert.That(all.Select(h => h.Position), Is.EqualTo(new[] { 0, 1, 2, 3 }));

        // Sequence shorter than the motif: no window exists (fuzznuc would report an overhanging hit; documented quirk).
        Assert.That(GenomicAnalyzer.FindKnownMotifHits(new DnaSequence("AC"), new[] { "AAA" }, bothStrands: true, maxMismatches: 2), Is.Empty);
    }

    [Test]
    public void Contract_NegativeBound_NullArguments_InvalidIupac()
    {
        var dna = new DnaSequence("ACGTACGT");
        Assert.That(() => GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "ACG" }, maxMismatches: -1),
            NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
        Assert.That(() => GenomicAnalyzer.FindKnownMotifHits((DnaSequence)null!, new[] { "A" }, maxMismatches: 1), NUnit.Framework.Throws.ArgumentNullException);
        Assert.That(() => GenomicAnalyzer.FindKnownMotifHits(dna, null!, maxMismatches: 1), NUnit.Framework.Throws.ArgumentNullException);
        Assert.That(() => GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "ACX" }, allowIupac: true, maxMismatches: 1),
            NUnit.Framework.Throws.ArgumentException);
        Assert.That(GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "", "  " }, maxMismatches: 1), Is.Empty);

        // Keys are upper-cased distinct motifs, as in the exact search.
        var hits = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "acg", "ACG" }, maxMismatches: 1);
        Assert.That(hits.Keys, Is.EqualTo(new[] { "ACG" }));
    }

    #endregion

    #region ApproximateMatcher IUPAC-pattern overload (canonical Hamming matcher, B05)

    [Test]
    public void ApproximateMatcher_IupacPattern_CountsBaseSetMembership()
    {
        // pattern RAATTC vs GGATTT: R/G ok, A/G x, A/A ok, T, T ok, C/T x -> 2 mismatches at pattern positions 1 and 5
        var hit = ApproximateMatcher.FindWithMismatches("CCGGATTTCC", "RAATTC", 2, iupacPattern: true).Single(h => h.Position == 2);
        Assert.That(hit.Position, Is.EqualTo(2));
        Assert.That(hit.Distance, Is.EqualTo(2));
        Assert.That(hit.MismatchPositions, Is.EqualTo(new[] { 1, 5 }));
        Assert.That(hit.MatchedSequence, Is.EqualTo("GGATTT"));
        Assert.That(ApproximateMatcher.FindWithMismatches("CCGGATTTCC", "RAATTC", 1, iupacPattern: true).Select(h => h.Position), Does.Not.Contain(2));
    }

    private static List<string> Project(IEnumerable<ApproximateMatchResult> hits)
        => hits.Select(h => $"{h.Position}:{h.MatchedSequence}:{h.Distance}:{string.Join(",", h.MismatchPositions)}").ToList();

    [Test]
    public void ApproximateMatcher_IupacFlagFalse_AndAcgtPatterns_EqualTheLiteralOverloads()
    {
        var rnd = new Random(77);
        for (int t = 0; t < 100; t++)
        {
            string s = new string(Enumerable.Range(0, rnd.Next(5, 60)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            string p = new string(Enumerable.Range(0, rnd.Next(1, 6)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            int k = rnd.Next(0, 3);
            var literal = Project(ApproximateMatcher.FindWithMismatches(s, p, k));
            Assert.That(Project(ApproximateMatcher.FindWithMismatches(s, p, k, false)), Is.EqualTo(literal));
            Assert.That(Project(ApproximateMatcher.FindWithMismatches(s, p, k, true)), Is.EqualTo(literal));
        }
    }

    [Test]
    public void ApproximateMatcher_IupacPattern_SequenceAmbiguityIsAlwaysAMismatch_InvalidPatternThrows()
    {
        Assert.That(ApproximateMatcher.FindWithMismatches("ACNGT", "ACNGT", 0, true), Is.Empty,
            "a sequence N matches no pattern code (literal sequence side, = MotifFinder.FindDegenerateMotif)");
        var one = ApproximateMatcher.FindWithMismatches("ACNGT", "ACNGT", 1, true).Single();
        Assert.That(one.MismatchPositions, Is.EqualTo(new[] { 2 }));
        Assert.That(() => ApproximateMatcher.FindWithMismatches("ACGT", "AXG", 1, true), NUnit.Framework.Throws.ArgumentException);
        Assert.That(() => ApproximateMatcher.FindWithMismatches("ACGT", "ACG", -1, true), NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
        Assert.That(ApproximateMatcher.FindWithMismatches("", "A", 1, true), Is.Empty);
        Assert.That(ApproximateMatcher.FindWithMismatches("ACGT", "", 1, true), Is.Empty);
        Assert.That(ApproximateMatcher.FindWithMismatches("AC", "ACG", 1, true), Is.Empty);
    }

    #endregion
}
