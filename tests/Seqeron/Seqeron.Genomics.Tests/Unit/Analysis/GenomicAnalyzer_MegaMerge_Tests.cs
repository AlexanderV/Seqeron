// GENOMIC-COMMON-001 (B09 audit A3-13b, F149) — EMBOSS megamerger: GenomicAnalyzer.MegaMerge.
// Oracle: scripts/analysis/generate_megamerger_oracle.py — 600 related / unrelated DNA pairs (substitutions, indels, N, IUPAC letters, lower-case stretches, flanks giving every start / end case),
// word size 2-16, -prefer N / Y; the merged sequence and the report text of megamerger 6.6.0. Evidence: docs/Evidence/GENOMIC-COMMON-001-Evidence.md.

using System.Text.Json;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_MegaMerge_Tests
{
    private sealed record Case(string Seq1, string Seq2, int WordSize, bool Prefer, string? Merged, string? Report, string Status);

    private static readonly Lazy<List<Case>> OracleCases = new(ReadOracle);

    private static List<Case> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.megamerger_oracle.jsonl";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var cases = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            using var doc = JsonDocument.Parse(line);
            var r = doc.RootElement;
            cases.Add(new Case(r.GetProperty("seq1").GetString()!, r.GetProperty("seq2").GetString()!, r.GetProperty("wordsize").GetInt32(), r.GetProperty("prefer").GetBoolean(),
                r.GetProperty("merged").GetString(), r.GetProperty("report").GetString(), r.GetProperty("status").GetString()!));
        }
        return cases;
    }

    [Test]
    public void MegaMerge_EqualsMegamerger_AllCases()
    {
        var cases = OracleCases.Value;
        Assert.That(cases, Has.Count.EqualTo(600));
        var failures = new List<string>();
        int merged = 0, mismatchRegions = 0, neitherStart = 0, neitherEnd = 0, prefer = 0, finalSeq1 = 0, finalSeq2 = 0, initialSeq1 = 0, initialSeq2 = 0, closer = 0, nomerge = 0, lower = 0, upper = 0;
        foreach (var c in cases)
        {
            var actual = GenomicAnalyzer.MegaMerge(c.Seq1, c.Seq2, c.WordSize, c.Prefer, "seqA", "seqB");
            if (c.Status == "nomerge")
            {
                nomerge++;
                if (actual.MergedSequence.Length != 0 || actual.Report.Length != 0)
                    failures.Add($"w={c.WordSize}: the program writes nothing but the port merged {actual.MergedSequence}");
                continue;
            }

            if (actual.MergedSequence != c.Merged || actual.Report != c.Report)
                failures.Add($"w={c.WordSize} prefer={c.Prefer}\n{c.Seq1}\n{c.Seq2}\n--- expected\n{c.Merged}\n{c.Report}--- actual\n{actual.MergedSequence}\n{actual.Report}");
            merged++;
            if (c.Report!.Contains("Mismatch region found")) mismatchRegions++;
            if (c.Report.Contains("overlap is at the start")) neitherStart++;
            if (c.Report.Contains("overlap is at the end")) neitherEnd++;
            if (c.Prefer) prefer++;
            if (c.Report.Contains("Using seqA") && c.Report.Contains("as the final sequence")) finalSeq1++;
            if (c.Report.Contains("Using seqB") && c.Report.Contains("as the final sequence")) finalSeq2++;
            if (c.Report.Contains("as the initial sequence") && c.Report.Contains("Using seqA 1-")) initialSeq1++;
            if (c.Report.Contains("as the initial sequence") && c.Report.Contains("Using seqB 1-")) initialSeq2++;
            if (c.Report.Contains("Mismatch is closer")) closer++;
            if (c.Merged!.Any(char.IsLower)) lower++;
            if (c.Merged!.Any(char.IsUpper)) upper++;
        }

        Assert.That(failures, Is.Empty, string.Join("\n=====\n", failures.Take(2)));
        Assert.That(merged, Is.GreaterThan(500));
        Assert.That(nomerge, Is.GreaterThan(30));
        Assert.That(mismatchRegions, Is.GreaterThan(200));
        Assert.That(neitherStart, Is.GreaterThan(20));
        Assert.That(neitherEnd, Is.GreaterThan(20));
        Assert.That(prefer, Is.GreaterThan(100));
        Assert.That(finalSeq1, Is.GreaterThan(50));
        Assert.That(finalSeq2, Is.GreaterThan(50));
        Assert.That(initialSeq1, Is.GreaterThan(30));
        Assert.That(initialSeq2, Is.GreaterThan(30));
        Assert.That(closer, Is.GreaterThan(100));
        Assert.That(lower, Is.GreaterThan(500));
        Assert.That(upper, Is.GreaterThan(300));
    }

    [Test]
    public void WorkedValues_RunOnTheBinary()
    {
        // megamerger a=ACGTACGTAC b=ACGTACGTACGGGGGG -wordsize 4: one match 1-10, then b 11-16 as the final sequence
        var r = GenomicAnalyzer.MegaMerge("ACGTACGTAC", "ACGTACGTACGGGGGG", 4, false, "a", "b");
        Assert.That(r.MergedSequence, Is.EqualTo("acgtacgtacgggggg"));
        Assert.That(r.Report, Is.EqualTo("# Report of megamerger of: a and b\n\na overlap starts at 1\nb overlap starts at 1\n\nMatching region a 1-10 : b 1-10\nLength of match: 10\n\na overlap ends at 10\nb overlap ends at 10\n\nUsing b 11-16 as the final sequence\n\n"));
        Assert.That(r.Matches.Single(), Is.EqualTo(new WordMatch(0, 1, 0, 0, 10)));
        // the second sequence shorter than the word: the program only warns and writes nothing; the first shorter than the word: "No match found" (exit 1)
        Assert.That(GenomicAnalyzer.MegaMerge("ACGTACGTACGGGGGG", "ACGTACGTAC", 12).MergedSequence, Is.Empty);
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.MegaMerge("ACGTACGTAC", "ACGTACGTACGGGGGG", 20));
    }

    [Test]
    public void Validation()
    {
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.MegaMerge(null!, "ACGT"));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.MegaMerge("ACGTACGT", "ACGTACGT", 1));
        Assert.Throws<OperationCanceledException>(() => GenomicAnalyzer.MegaMerge("ACGTACGT", "ACGTACGT", 4, cancellationToken: new CancellationToken(true)));
    }

    // ------------------------------------------------------------------ A4-6: MinimalMatchSet without the O(k) list shifts and whole-list re-sorts

    /// <summary>embWordMatchMin as it was before A4-6 (List.RemoveAt(0) / RemoveAt(i), the whole list sorted again after a truncation): the reference the new loop must equal.</summary>
    private static List<(int S1, int S2, int Len)> ReferenceMinimalMatchSet(IEnumerable<(int S1, int S2, int Len)> matches, int wordLength, out string? exception)
    {
        exception = null;
        var matchlist = matches.Select(m => new GenomicAnalyzer.MegaMatch(m.S1, m.S2, m.Len)).ToList();
        matchlist.Sort(GenomicAnalyzer.MegaOrder);
        var minlist = new List<GenomicAnalyzer.MegaMatch>();
        try
        {
            while (matchlist.Count > 0)
            {
                var thismatch = matchlist[0];
                matchlist.RemoveAt(0);
                minlist.Add(thismatch);
                int deadx1 = thismatch.Start1, deady1 = thismatch.Start2;
                int deadx2 = thismatch.Start1 + thismatch.Length - 1, deady2 = thismatch.Start2 + thismatch.Length - 1;
                bool truncated = false;
                for (int i = 0; i < matchlist.Count;)
                {
                    int result = GenomicAnalyzer.DeadZone(matchlist[i], deadx1, deady1, deadx2, deady2, wordLength);
                    if (result == 1)
                    {
                        matchlist.RemoveAt(i);
                        continue;
                    }

                    if (result == 2)
                        truncated = true;
                    i++;
                }

                if (truncated)
                    matchlist.Sort(GenomicAnalyzer.MegaOrder);
            }
        }
        catch (InvalidOperationException e)
        {
            exception = e.Message;
            return [];
        }

        minlist.Sort((a, b) => a.Start1 != b.Start1 ? a.Start1.CompareTo(b.Start1) : a.Start2.CompareTo(b.Start2));
        return minlist.Select(m => (m.Start1, m.Start2, m.Length)).ToList();
    }

    private static List<(int S1, int S2, int Len)> NewMinimalMatchSet(IEnumerable<(int S1, int S2, int Len)> matches, int wordLength, out string? exception)
    {
        exception = null;
        var list = matches.Select(m => new GenomicAnalyzer.MegaMatch(m.S1, m.S2, m.Len)).ToList();
        try
        {
            GenomicAnalyzer.MinimalMatchSet(list, wordLength);
        }
        catch (InvalidOperationException e)
        {
            exception = e.Message;
            return [];
        }

        return list.Select(m => (m.Start1, m.Start2, m.Length)).ToList();
    }

    [Test]
    public void MinimalMatchSet_EqualsTheListShiftingOriginal_OnRandomRepetitiveMatchLists_ExceptionsIncluded()
    {
        var rnd = new Random(20261127);
        int truncating = 0, threw = 0, total = 0;
        for (int round = 0; round < 400; round++)
        {
            // repetitive sequences give crossing dead zones (shortened matches) and, now and then, the program's own "same diagonal" abort
            string unit = new(Enumerable.Range(0, rnd.Next(1, 7)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            string a = string.Concat(Enumerable.Repeat(unit, rnd.Next(5, 40))) + new string(Enumerable.Range(0, rnd.Next(0, 30)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            string b = new string(Enumerable.Range(0, rnd.Next(0, 20)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray()) + string.Concat(Enumerable.Repeat(unit, rnd.Next(5, 40)));
            int w = rnd.Next(2, 9);
            var words = GenomicAnalyzer.SeqMatchAll([a, b], w).Select(m => (m.Position1, m.Position2, m.Length)).ToList();
            if (words.Count == 0)
                continue;
            total++;
            var expected = ReferenceMinimalMatchSet(words, w, out string? expectedException);
            var actual = NewMinimalMatchSet(words, w, out string? actualException);
            Assert.That(actualException, Is.EqualTo(expectedException), $"round {round}");
            Assert.That(actual, Is.EqualTo(expected), $"round {round}: {a} / {b} w={w}");
            if (expectedException is not null) threw++;
            if (expected.Count > 0 && words.Count > expected.Count) truncating++;
        }

        Assert.That(total, Is.GreaterThan(300));
        Assert.That(truncating, Is.GreaterThan(100), "rounds in which the minimal set is smaller than the match list");
        Assert.That(threw, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void MinimalMatchSet_OfSixtyThousandDiagonalMatches_IsFast()
    {
        // the matches of two poly-A sequences of 30 000: one on every diagonal (60 000 of them); the longest one kills all the others (60 000 removals: List.RemoveAt shifted the tail every time)
        int n = 30_000, w = 20;
        var list = new List<GenomicAnalyzer.MegaMatch>();
        for (int d = -(n - w); d <= n - w; d++)
            list.Add(new GenomicAnalyzer.MegaMatch(Math.Max(0, d), Math.Max(0, -d), n - Math.Abs(d)));
        var sw = System.Diagnostics.Stopwatch.StartNew();
        GenomicAnalyzer.MinimalMatchSet(list, w);
        sw.Stop();
        Assert.That(list, Has.Count.EqualTo(1));
        Assert.That((list[0].Start1, list[0].Start2, list[0].Length), Is.EqualTo((0, 0, n)));
        Assert.That(sw.Elapsed.TotalSeconds, Is.LessThan(10), $"MinimalMatchSet of 60 000 matches took {sw.Elapsed.TotalSeconds:F2} s");
    }

    [Test]
    public void PolyA_IsTheProgramsResult_ForTheLengthsTheProgramFinishes()
    {
        // megamerger 6.6.0 on two poly-A sequences of 600 / 1 500 residues (1 s / 20 s: its own minimal-match loop is quadratic too): the whole sequence in lower case, "Matching region a 1-n : b 1-n", "Length of match: n"
        foreach (int n in new[] { 600, 1500 })
        {
            var r = GenomicAnalyzer.MegaMerge(new string('A', n), new string('A', n), 20);
            Assert.That(r.MergedSequence, Is.EqualTo(new string('a', n)));
            Assert.That(r.Report, Does.Contain($"Matching region seq1 1-{n} : seq2 1-{n}").And.Contain($"Length of match: {n}"));
            Assert.That(r.Matches, Is.EqualTo(new[] { new WordMatch(0, 1, 0, 0, n) }));
        }
    }

    [Test]
    public void TandemRepeat_MegaMerge_IsTheProgramsResult_AndCancellable()
    {
        // a tandem repeat of a 7-mer against itself (megamerger 6.6.0 on the 420-base version: one region 1-420)
        string tandem = string.Concat(Enumerable.Repeat("ACGTTGA", 570));
        var r = GenomicAnalyzer.MegaMerge(tandem, tandem, 12);
        Assert.That(r.MergedSequence, Is.EqualTo(tandem.ToLowerInvariant()));
        Assert.That(r.Matches, Is.EqualTo(new[] { new WordMatch(0, 1, 0, 0, tandem.Length) }));
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() => GenomicAnalyzer.MegaMerge(tandem, tandem, 12, cancellationToken: cts.Token));
    }
}
