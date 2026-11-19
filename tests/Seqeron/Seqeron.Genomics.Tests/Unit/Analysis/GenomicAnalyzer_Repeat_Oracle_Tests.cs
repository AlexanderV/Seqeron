// GENOMIC-REPEAT-001 (B09 F28-F30) — FindRepeats / FindMaximalRepeats / FindLongestRepeat / FindAllLongestRepeats
// against independent oracles. Evidence: docs/Evidence/GENOMIC-REPEAT-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-REPEAT-001.md
// Oracle 1: Python (brute-force substring dictionary == suffix array + Kasai LCP + lcp-intervals, asserted equal),
//           628 (sequence, minLength) rows, scripts/analysis/generate_repeat_oracle.py (TestData/Analysis/repeat_oracle.tsv).
// Oracle 2: in-test brute force (every substring, left/right diversity straight from the text) on random strings.

using System.Diagnostics;
using System.Security.Cryptography;
using System.Text;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_Repeat_Oracle_Tests
{
    private sealed record Row(string Sequence, int MinLength, int Count, string Digest, string Maximal, string Longest);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.repeat_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], int.Parse(f[1]), int.Parse(f[2]), f[3], f[4], f[5]));
        }
        return rows;
    }

    private static string Format(IEnumerable<RepeatInfo> items)
    {
        string s = string.Join(";", items.Select(r => $"{r.Sequence}:{string.Join(",", r.Positions)}"));
        return s.Length == 0 ? "-" : s;
    }

    private static string Digest(IEnumerable<RepeatInfo> items)
    {
        var sb = new StringBuilder();
        foreach (var r in items)
            sb.Append(r.Sequence).Append(':').Append(string.Join(",", r.Positions)).Append('\n');
        return Convert.ToHexString(SHA256.HashData(Encoding.ASCII.GetBytes(sb.ToString()))).ToLowerInvariant()[..16];
    }

    #region Python oracle

    [Test]
    public void FindRepeats_MatchesPythonOracle_SetPositionsAndOrder()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(628));
        foreach (var row in rows)
        {
            var actual = GenomicAnalyzer.FindRepeats(new DnaSequence(row.Sequence), row.MinLength).ToList();
            Assert.That(actual, Has.Count.EqualTo(row.Count), $"{row.Sequence} minLength={row.MinLength}");
            Assert.That(Digest(actual), Is.EqualTo(row.Digest), $"{row.Sequence} minLength={row.MinLength}");
        }
    }

    [Test]
    public void FindMaximalRepeats_MatchesPythonOracle_GusfieldDefinition()
    {
        int nonEmpty = 0;
        foreach (var row in ReadOracle())
        {
            var actual = GenomicAnalyzer.FindMaximalRepeats(new DnaSequence(row.Sequence), row.MinLength).ToList();
            Assert.That(Format(actual), Is.EqualTo(row.Maximal), $"{row.Sequence} minLength={row.MinLength}");
            if (actual.Count > 0) nonEmpty++;
        }
        Assert.That(nonEmpty, Is.GreaterThan(300), "the corpus must exercise non-trivial maximal-repeat sets");
    }

    [Test]
    public void FindAllLongestRepeats_AndFindLongestRepeat_MatchPythonOracle()
    {
        foreach (var row in ReadOracle().Where(r => r.MinLength == 1))
        {
            var dna = new DnaSequence(row.Sequence);
            var all = GenomicAnalyzer.FindAllLongestRepeats(dna);
            Assert.That(Format(all), Is.EqualTo(row.Longest), row.Sequence);

            var lrs = GenomicAnalyzer.FindLongestRepeat(dna);
            if (row.Longest == "-")
                Assert.That(lrs.IsEmpty, Is.True, row.Sequence);
            else
                Assert.That(Format(new[] { lrs }), Is.EqualTo(row.Longest.Split(';')[0]),
                    $"{row.Sequence}: FindLongestRepeat must be the tied repeat with the leftmost first occurrence");
        }
    }

    #endregion

    #region In-test brute force on random strings

    private static (string Sequence, List<int> Positions)[] BruteRepeats(string s, int minLength)
    {
        var occ = new Dictionary<string, List<int>>(StringComparer.Ordinal);
        for (int i = 0; i < s.Length; i++)
            for (int j = i + 1; j <= s.Length; j++)
            {
                string sub = s.Substring(i, j - i);
                if (!occ.TryGetValue(sub, out var list))
                    occ[sub] = list = new List<int>();
                list.Add(i);
            }
        var rank = new int[s.Length];
        var order = Enumerable.Range(0, s.Length).OrderBy(p => s.Substring(p), StringComparer.Ordinal).ToArray();
        for (int r = 0; r < order.Length; r++) rank[order[r]] = r;
        return occ.Where(kv => kv.Value.Count >= 2 && kv.Key.Length >= Math.Max(1, minLength))
            .OrderBy(kv => kv.Value.Min(p => rank[p])).ThenBy(kv => kv.Key.Length)
            .Select(kv => (kv.Key, kv.Value)).ToArray();
    }

    private static bool IsMaximal(string s, string sub, List<int> positions)
    {
        var left = positions.Select(p => p == 0 ? '^' : s[p - 1]).Distinct().Count();
        var right = positions.Select(p => p + sub.Length == s.Length ? '$' : s[p + sub.Length]).Distinct().Count();
        return left >= 2 && right >= 2;
    }

    [Test]
    public void FindRepeats_And_FindMaximalRepeats_MatchBruteForce_OnRandomStrings(
        [Values("ACGT", "AC", "A", "ACG")] string alphabet)
    {
        var rng = new Random(alphabet.Length * 7919);
        for (int trial = 0; trial < 120; trial++)
        {
            int n = rng.Next(1, 70);
            string s = new(Enumerable.Range(0, n).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
            var dna = new DnaSequence(s);
            foreach (int minLength in new[] { -3, 1, 2, 4 })
            {
                var brute = BruteRepeats(s, minLength);
                var actual = GenomicAnalyzer.FindRepeats(dna, minLength).ToList();
                Assert.That(actual.Select(r => r.Sequence).ToArray(), Is.EqualTo(brute.Select(b => b.Sequence).ToArray()),
                    $"{s} minLength={minLength}");
                for (int i = 0; i < brute.Length; i++)
                    Assert.That(actual[i].Positions.ToArray(), Is.EqualTo(brute[i].Positions.ToArray()), $"{s} '{brute[i].Sequence}'");

                var maximal = brute.Where(b => IsMaximal(s, b.Sequence, b.Positions)).ToArray();
                var actualMax = GenomicAnalyzer.FindMaximalRepeats(dna, minLength).ToList();
                Assert.That(actualMax.Select(r => r.Sequence).ToArray(), Is.EqualTo(maximal.Select(b => b.Sequence).ToArray()),
                    $"maximal: {s} minLength={minLength}");
                Assert.That(actualMax.Count, Is.LessThanOrEqualTo(Math.Max(1, n)), "at most n maximal repeats (Gusfield §7.12)");
            }
        }
    }

    #endregion

    #region Semantics and contract

    [Test]
    public void FindMaximalRepeats_ACGTACGT_IsOnlyACGT_WhereFindRepeatsListsEveryPrefixAndSuffix()
    {
        var dna = new DnaSequence("ACGTACGT");
        Assert.That(Format(GenomicAnalyzer.FindMaximalRepeats(dna, 1)), Is.EqualTo("ACGT:0,4"));
        Assert.That(GenomicAnalyzer.FindRepeats(dna, 1).Select(r => r.Sequence),
            Is.EqualTo(new[] { "A", "AC", "ACG", "ACGT", "C", "CG", "CGT", "G", "GT", "T" }));
    }

    [Test]
    public void FindMaximalRepeats_Homopolymer_EveryProperRunIsMaximal_BecauseTheTextEndsAreDiverse()
    {
        // A^k (k <= n-1) occurs at 0 (text start on the left) and 1 (A on the left): left-diverse;
        // the occurrence at n-k ends at the text end: right-diverse.
        var result = GenomicAnalyzer.FindMaximalRepeats(new DnaSequence("AAAAA"), 1).ToList();
        Assert.That(result.Select(r => r.Sequence), Is.EqualTo(new[] { "A", "AA", "AAA", "AAAA" }));
        Assert.That(result[3].Positions, Is.EqualTo(new[] { 0, 1 }));
    }

    [Test]
    public void FindMaximalRepeats_IsSubsetOfFindRepeats_AndContainsEveryLongestRepeat()
    {
        var dna = new DnaSequence("GATTACAGATTACATTACAGATCAGATTACA");
        var all = GenomicAnalyzer.FindRepeats(dna, 2).Select(r => r.Sequence).ToHashSet();
        var maximal = GenomicAnalyzer.FindMaximalRepeats(dna, 2).ToList();
        Assert.That(maximal.Select(r => r.Sequence).All(all.Contains), Is.True);
        foreach (var lrs in GenomicAnalyzer.FindAllLongestRepeats(dna))
            Assert.That(maximal.Any(m => m.Sequence == lrs.Sequence), Is.True, $"a longest repeat is always maximal: {lrs.Sequence}");
    }

    [Test]
    public void FindMaximalRepeats_Strings_EqualTheStringsOfTheMaximalPairsOfRepeatFinder()
    {
        // Gusfield: a string is a maximal repeat iff some pair of its occurrences is left- and right-maximal
        // (REPuter / repeat-match -f pairs = RepeatFinder.FindDirectRepeats without spacing or length limits).
        var rng = new Random(77);
        foreach (string alphabet in new[] { "ACGT", "AC", "ACG" })
            for (int trial = 0; trial < 80; trial++)
            {
                string s = new(Enumerable.Range(0, rng.Next(2, 90)).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
                var dna = new DnaSequence(s);
                var fromPairs = RepeatFinder.FindDirectRepeats(dna, minLength: 2, maxLength: int.MaxValue, minSpacing: -1_000_000)
                    .Select(d => d.RepeatSequence).Distinct().OrderBy(x => x, StringComparer.Ordinal).ToArray();
                var maximal = GenomicAnalyzer.FindMaximalRepeats(dna, 2).Select(r => r.Sequence).OrderBy(x => x, StringComparer.Ordinal).ToArray();
                Assert.That(maximal, Is.EqualTo(fromPairs), s);
            }
    }

    [Test]
    public void NullSequence_ThrowsArgumentNullException_AtTheCall()
    {
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindRepeats((DnaSequence)null!, 3));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindMaximalRepeats((DnaSequence)null!, 3));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindLongestRepeat((DnaSequence)null!));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindAllLongestRepeats((DnaSequence)null!));
    }

    [Test]
    public void EmptyAndRepeatFreeSequences_YieldNothing()
    {
        foreach (var text in new[] { "", "A", "ACGT", "ACGTC".Substring(0, 3) })
        {
            var dna = new DnaSequence(text);
            Assert.That(GenomicAnalyzer.FindRepeats(dna, 1), Is.Empty, text);
            Assert.That(GenomicAnalyzer.FindMaximalRepeats(dna, 1), Is.Empty, text);
            Assert.That(GenomicAnalyzer.FindAllLongestRepeats(dna), Is.Empty, text);
            Assert.That(GenomicAnalyzer.FindLongestRepeat(dna).IsEmpty, Is.True, text);
        }
    }

    [Test]
    public void FindLongestRepeat_TieBreak_IsTheLeftmostFirstOccurrence_NotTreeOrder()
    {
        // Brute force (generate_repeat_oracle.py): the longest repeats of TTGGCACGTAACGTCTTGG are TTGG@{0,15} and ACGT@{5,10}
        // (length 4, no longer substring repeats). The suffix tree visits the 'A' branch before the 'T' branch, so a
        // tree-order choice would return ACGT; the contract is the leftmost first occurrence, TTGG.
        var dna = new DnaSequence("TTGGCACGTAACGTCTTGG");
        Assert.That(Format(GenomicAnalyzer.FindAllLongestRepeats(dna)), Is.EqualTo("TTGG:0,15;ACGT:5,10"));
        var lrs = GenomicAnalyzer.FindLongestRepeat(dna);
        Assert.That(lrs.Sequence, Is.EqualTo("TTGG"));
        Assert.That(lrs.Positions, Is.EqualTo(new[] { 0, 15 }));
    }

    #endregion

    #region Complexity (1 Mbp)

    [Test]
    public void FindRepeats_And_FindMaximalRepeats_OneMegabase_RunInSeconds_AndAreVerified()
    {
        const int n = 1_000_000;
        var rng = new Random(20261118);
        var chars = new char[n];
        for (int i = 0; i < n; i++) chars[i] = "ACGT"[rng.Next(4)];
        // plant a 200-nt repeat three times, the second copy with a different left flank and the third with a different right flank
        string unit = new(chars, 500_000, 200);
        unit.CopyTo(0, chars, 100_000, 200);
        unit.CopyTo(0, chars, 800_000, 200);
        string text = new(chars);
        var dna = new DnaSequence(text);

        var sw = Stopwatch.StartNew();
        var repeats = GenomicAnalyzer.FindRepeats(dna, 20).ToList();
        var maximal = GenomicAnalyzer.FindMaximalRepeats(dna, 20).ToList();
        sw.Stop();

        Assert.Multiple(() =>
        {
            Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(60), $"1 Mbp FindRepeats + FindMaximalRepeats took {sw.Elapsed.TotalSeconds:F1} s");
            Assert.That(repeats.All(r => r.Length >= 20 && r.Count >= 2), Is.True);
            Assert.That(maximal.Count, Is.LessThanOrEqualTo(repeats.Count));
            var planted = maximal.Where(r => r.Length >= 190).ToList();
            Assert.That(planted.Count, Is.GreaterThanOrEqualTo(1), "the planted 200-mer (or its maximal core) must be reported");
            foreach (var r in maximal.Take(300))
            {
                Assert.That(r.Positions, Is.Ordered.Ascending.And.Unique);
                foreach (int p in r.Positions)
                    Assert.That(string.CompareOrdinal(text, p, r.Sequence, 0, r.Length), Is.EqualTo(0));
                Assert.That(IsMaximal(text, r.Sequence, r.Positions.ToList()), Is.True, "reported maximal repeat must be left- and right-diverse");
            }
        });
        TestContext.Out.WriteLine($"1 Mbp: {repeats.Count} repeats, {maximal.Count} maximal, {sw.Elapsed.TotalSeconds:F1} s");
    }

    #endregion
}
