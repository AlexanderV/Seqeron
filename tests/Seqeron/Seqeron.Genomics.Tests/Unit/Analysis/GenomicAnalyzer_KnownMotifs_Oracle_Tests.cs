// GENOMIC-MOTIFS-001 (B09 F18-F20) — FindMotif / FindKnownMotifs / FindKnownMotifHits against an independent oracle.
// Evidence: docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-MOTIFS-001.md
// Oracle: EMBOSS fuzznuc 6.6.0 (-pmismatch 0 -complement Y) == independent Python `re` lookahead, 2448
//         (sequence, pattern) pairs, scripts/analysis/generate_known_motif_oracle.py
//         (TestData/Analysis/known_motif_oracle.tsv). Reverse hits are forward-strand window starts.

using System.Diagnostics;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_KnownMotifs_Oracle_Tests
{
    private static readonly char[] Acgt = { 'A', 'C', 'G', 'T' };

    private sealed record OracleRow(string Sequence, string Pattern, int[] Forward, int[] Reverse);

    // One parsed copy per fixture (audit A3-15): the oracle was re-read and re-parsed on every call, i.e. in every test.
    private static readonly Lazy<List<OracleRow>> OracleRows = new(ReadOracleRows);

    private static List<OracleRow> ReadOracle() => OracleRows.Value;

    private static List<OracleRow> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.known_motif_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<OracleRow>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new OracleRow(f[0], f[1], Parse(f[2]), Parse(f[3])));
        }
        return rows;

        static int[] Parse(string s) => s == "-" ? Array.Empty<int>() : s.Split(',').Select(int.Parse).ToArray();
    }

    private static bool IsExact(string pattern) => pattern.All(c => Array.IndexOf(Acgt, c) >= 0);

    #region fuzznuc oracle

    [Test]
    public void FindKnownMotifHits_BothStrandsIupac_MatchesFuzznucOracle()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(2448));
        int withReverse = 0;
        foreach (var row in rows)
        {
            var hits = GenomicAnalyzer.FindKnownMotifHits(new DnaSequence(row.Sequence), new[] { row.Pattern },
                bothStrands: true, allowIupac: true);

            var expected = row.Forward.Select(p => new KnownMotifHit(p, false))
                .Concat(row.Reverse.Select(p => new KnownMotifHit(p, true)))
                .OrderBy(h => h.Position).ThenBy(h => h.IsReverseStrand).ToArray();
            if (expected.Length == 0)
            {
                Assert.That(hits, Is.Empty, $"{row.Sequence} / {row.Pattern}");
                continue;
            }

            withReverse += row.Reverse.Length > 0 ? 1 : 0;
            Assert.That(hits[row.Pattern], Is.EqualTo(expected), $"{row.Sequence} / {row.Pattern}");
        }

        Assert.That(withReverse, Is.GreaterThan(700), "the oracle must exercise the reverse strand");
    }

    [Test]
    public void FindKnownMotifHits_ForwardOnlyIupac_MatchesFuzznucForwardHits()
    {
        foreach (var row in ReadOracle())
        {
            var hits = GenomicAnalyzer.FindKnownMotifHits(new DnaSequence(row.Sequence), new[] { row.Pattern },
                bothStrands: false, allowIupac: true);
            var actual = hits.TryGetValue(row.Pattern, out var h) ? h.Select(x => x.Position).ToArray() : Array.Empty<int>();
            Assert.That(actual, Is.EqualTo(row.Forward), $"{row.Sequence} / {row.Pattern}");
            Assert.That(hits.Values.SelectMany(x => x), Has.None.Matches<KnownMotifHit>(x => x.IsReverseStrand));
        }
    }

    // Exact (pure ACGT) patterns: FindMotif, FindKnownMotifs and FindKnownMotifHits (exact mode, both strands)
    // all equal the oracle; IUPAC symbols are literal in exact mode and match nothing in a DNA sequence.
    [Test]
    public void ExactMode_AllEntryPoints_MatchOracle()
    {
        foreach (var row in ReadOracle())
        {
            var dna = new DnaSequence(row.Sequence);
            var viaFindMotif = GenomicAnalyzer.FindMotif(dna, row.Pattern);
            var viaKnown = GenomicAnalyzer.FindKnownMotifs(dna, new[] { row.Pattern });
            var viaHits = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { row.Pattern }, bothStrands: true);

            if (IsExact(row.Pattern))
            {
                Assert.That(viaFindMotif, Is.EqualTo(row.Forward), $"FindMotif {row.Sequence} / {row.Pattern}");
                Assert.That(viaKnown.TryGetValue(row.Pattern, out var k) ? k : Array.Empty<int>(), Is.EqualTo(row.Forward));
                var rev = viaHits.TryGetValue(row.Pattern, out var h)
                    ? h.Where(x => x.IsReverseStrand).Select(x => x.Position).ToArray()
                    : Array.Empty<int>();
                Assert.That(rev, Is.EqualTo(row.Reverse), $"reverse {row.Sequence} / {row.Pattern}");
            }
            else
            {
                Assert.That(viaFindMotif, Is.Empty, $"literal IUPAC symbol {row.Pattern}");
                Assert.That(viaKnown, Is.Empty);
                Assert.That(viaHits, Is.Empty);
            }
        }
    }

    // Worked fuzznuc example: GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC.
    //  GAATTC (palindrome) -> forward 0, 9 and reverse 0, 9 (fuzznuc prints both strands at 1-6 and 10-15);
    //  TATAWAW             -> forward 17 (fuzznuc 18-24 '+'), reverse 16 (17-23 '-').
    [Test]
    public void FindKnownMotifHits_FuzznucWorkedExample()
    {
        var dna = new DnaSequence("GAATTCAAAGAATTCTTTATATAAAGTTCAAGAAC");

        var eco = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "GAATTC" }, bothStrands: true)["GAATTC"];
        Assert.That(eco, Is.EqualTo(new[]
        {
            new KnownMotifHit(0, false), new KnownMotifHit(0, true),
            new KnownMotifHit(9, false), new KnownMotifHit(9, true),
        }), "a palindromic site is a hit on both strands at the same forward position");

        var tata = GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "tataWaW" }, bothStrands: true, allowIupac: true);
        Assert.That(tata["TATAWAW"], Is.EqualTo(new[] { new KnownMotifHit(16, true), new KnownMotifHit(17, false) }));
    }

    #endregion

    #region FindMotif contract (B05 + B06 cross-batch requests)

    // B05 request: GenomicAnalyzer.FindMotif returned suffix-tree (DFS) order; it is now ascending like
    // MotifFinder.FindExactMotif. The sequence below puts the suffix-tree order out of sequence.
    [Test]
    public void FindMotif_PositionsAscending_EqualToMotifFinderAndStringIndexOf()
    {
        var rnd = new Random(18);
        for (int trial = 0; trial < 200; trial++)
        {
            string s = new string(Enumerable.Range(0, rnd.Next(1, 120)).Select(_ => "AC"[rnd.Next(2)]).ToArray());
            var dna = new DnaSequence(s);
            string motif = new string(Enumerable.Range(0, rnd.Next(1, 5)).Select(_ => "AC"[rnd.Next(2)]).ToArray());

            var expected = new List<int>();
            for (int i = s.IndexOf(motif, StringComparison.Ordinal); i >= 0; i = s.IndexOf(motif, i + 1, StringComparison.Ordinal))
                expected.Add(i);

            var actual = GenomicAnalyzer.FindMotif(dna, motif);
            Assert.That(actual, Is.EqualTo(expected), $"{s} / {motif}");
            Assert.That(actual, Is.EqualTo(MotifFinder.FindExactMotif(dna, motif).ToList()));
            Assert.That(GenomicAnalyzer.FindMotif(dna, motif.ToLowerInvariant()), Is.EqualTo(expected));
        }
    }

    // Rosalind SUBS sample (1-based 2 4 10 = 0-based 1 3 9); the suffix tree enumerates it as 1, 9, 3.
    [Test]
    public void FindMotif_RosalindSubsSample_AscendingNotTreeOrder()
    {
        var dna = new DnaSequence("GATATATGCATATACTT");
        Assert.That(dna.SuffixTree.FindAllOccurrences("ATAT"), Is.EquivalentTo(new[] { 1, 3, 9 }));
        Assert.That(GenomicAnalyzer.FindMotif(dna, "ATAT"), Is.EqualTo(new[] { 1, 3, 9 }));
    }

    [Test]
    public void FindMotif_Contract_NullEmptyAndOverlap()
    {
        var dna = new DnaSequence("AAAAA");
        Assert.That(GenomicAnalyzer.FindMotif(dna, "AAA"), Is.EqualTo(new[] { 0, 1, 2 }));
        Assert.That(GenomicAnalyzer.FindMotif(dna, ""), Is.Empty);
        Assert.That(GenomicAnalyzer.FindMotif(dna, null!), Is.Empty);
        Assert.That(GenomicAnalyzer.FindMotif(dna, "AAAAAA"), Is.Empty);
        Assert.That(GenomicAnalyzer.FindMotif(dna, "N"), Is.Empty, "IUPAC symbols are literal in the exact search");
        Assert.That(() => GenomicAnalyzer.FindMotif((DnaSequence)null!, "A"), NUnit.Framework.Throws.ArgumentNullException);
        Assert.That(() => GenomicAnalyzer.FindKnownMotifs((DnaSequence)null!, new[] { "A" }), NUnit.Framework.Throws.ArgumentNullException);
        Assert.That(() => GenomicAnalyzer.FindKnownMotifHits((DnaSequence)null!, new[] { "A" }), NUnit.Framework.Throws.ArgumentNullException);
        Assert.That(() => GenomicAnalyzer.FindKnownMotifHits(dna, null!), NUnit.Framework.Throws.ArgumentNullException);
    }

    #endregion

    #region Engine choice (KMP scan vs suffix-tree index) and IUPAC validation

    // Fewer than 128 distinct queries are KMP scans, 128 or more build the suffix tree: both must return the
    // same positions as a brute-force scan, for the same motifs.
    [Test]
    public void FindKnownMotifs_ScanAndIndexRoutes_ReturnIdenticalPositions()
    {
        var rnd = new Random(19);
        string s = new string(Enumerable.Range(0, 600).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
        s = s.Insert(100, "GAATTCGAATTC").Insert(400, "AAAAAAAA");
        var dna = new DnaSequence(s);

        var core = new[] { "GAATTC", "AAA", "ACGT", "A", "CCGG", "TTTTTTTTTT" };
        var filler = new List<string>();
        while (filler.Count < 200)
        {
            int len = 4 + filler.Count / 40;
            filler.Add(new string(Enumerable.Range(0, len).Select(_ => "ACGT"[rnd.Next(4)]).ToArray()));
        }
        var all = core.Concat(filler).ToList();
        Assert.That(all.Distinct().Count(), Is.GreaterThan(128));

        var scan = GenomicAnalyzer.FindKnownMotifs(dna, core);
        var indexed = GenomicAnalyzer.FindKnownMotifs(dna, all);

        foreach (string motif in all.Distinct())
        {
            var expected = new List<int>();
            for (int i = s.IndexOf(motif, StringComparison.Ordinal); i >= 0; i = s.IndexOf(motif, i + 1, StringComparison.Ordinal))
                expected.Add(i);
            if (expected.Count == 0)
                Assert.That(indexed.ContainsKey(motif), Is.False, motif);
            else
                Assert.That(indexed[motif], Is.EqualTo(expected), $"indexed route {motif}");
            if (core.Contains(motif))
            {
                if (expected.Count == 0)
                    Assert.That(scan.ContainsKey(motif), Is.False);
                else
                    Assert.That(scan[motif], Is.EqualTo(expected), $"scan route {motif}");
            }
        }

        // both strands double the query count: 70 distinct motifs x 2 strands crosses the index threshold too
        var seventy = all.Take(70).ToList();
        var asHits = GenomicAnalyzer.FindKnownMotifHits(dna, seventy, bothStrands: true);
        foreach (string motif in seventy)
        {
            string rc = DnaSequence.GetReverseComplementString(motif);
            var fwd = new List<int>();
            for (int i = s.IndexOf(motif, StringComparison.Ordinal); i >= 0; i = s.IndexOf(motif, i + 1, StringComparison.Ordinal))
                fwd.Add(i);
            var rev = new List<int>();
            for (int i = s.IndexOf(rc, StringComparison.Ordinal); i >= 0; i = s.IndexOf(rc, i + 1, StringComparison.Ordinal))
                rev.Add(i);
            var got = asHits.TryGetValue(motif, out var h) ? h : Array.Empty<KnownMotifHit>();
            Assert.That(got.Where(x => !x.IsReverseStrand).Select(x => x.Position), Is.EqualTo(fwd), motif);
            Assert.That(got.Where(x => x.IsReverseStrand).Select(x => x.Position), Is.EqualTo(rev), motif);
        }
    }

    [Test]
    public void FindKnownMotifHits_InvalidIupacSymbol_Throws_AndExactModeKeepsItLiteral()
    {
        var dna = new DnaSequence("ACGTACGT");
        Assert.That(() => GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "AXG" }, allowIupac: true),
            NUnit.Framework.Throws.ArgumentException);
        Assert.That(GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "AXG" }), Is.Empty);
        Assert.That(GenomicAnalyzer.FindKnownMotifHits(dna, new[] { "", "  ", null! }, bothStrands: true), Is.Empty);
    }

    // Complexity: one motif over 3 Mbp must not build a suffix tree (23 s / 1.7 GB at 10 Mbp, B06 measurement)
    // and still equal the string.IndexOf scan.
    [Test]
    public void FindKnownMotifs_SingleMotifOnLargeSequence_IsLinearScan()
    {
        var rnd = new Random(20);
        var chars = new char[3_000_000];
        for (int i = 0; i < chars.Length; i++) chars[i] = "ACGT"[rnd.Next(4)];
        string s = new string(chars);
        var dna = new DnaSequence(s);

        var sw = Stopwatch.StartNew();
        var result = GenomicAnalyzer.FindKnownMotifs(dna, new[] { "GAATTC" });
        sw.Stop();

        var expected = new List<int>();
        for (int i = s.IndexOf("GAATTC", StringComparison.Ordinal); i >= 0; i = s.IndexOf("GAATTC", i + 1, StringComparison.Ordinal))
            expected.Add(i);
        Assert.That(result["GAATTC"], Is.EqualTo(expected));
        // A5-4 complexity guard, not a latency bound: the single-motif path is a linear scan, the suffix-tree build of 3 Mbp takes minutes in Debug; the bound is 10x the earlier Debug budget so a loaded shared machine cannot fail it while a regression to the quadratic code still does
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(200), "a linear scan; the suffix-tree build of 3 Mbp takes far longer in Debug");
    }

    #endregion
}
