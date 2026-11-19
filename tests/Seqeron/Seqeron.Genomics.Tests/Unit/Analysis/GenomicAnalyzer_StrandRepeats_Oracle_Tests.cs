// GENOMIC-REPEAT-001 (B09 A2-33 / A2-35) — both-strand repeats, repeat pairs and the string overloads of the repeat finders.
// Oracle 1: MUMmer 3.23 repeat-match (default = both strands, -f = forward only) on 1900 (sequence, L) rows plus an independent
//           Python brute force of the both-strand repeat classes over S # revcomp(S) — scripts/analysis/generate_strand_repeat_oracle.py.
// Oracle 2: in-test brute force of the single-strand repeats / maximal repeats with N runs (a non-ACGT symbol never matches and is
//           its own left / right context).
// Evidence: docs/Evidence/GENOMIC-REPEAT-001-Evidence.md.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_StrandRepeats_Oracle_Tests
{
    private sealed record Row(string Sequence, int MinLength, string Maximal, string All, string Pairs, string ForwardPairs);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.strand_repeat_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], int.Parse(f[1]), f[2], f[3], f[4], f[5]));
        }
        return rows;
    }

    private static string Format(IEnumerable<StrandedRepeatInfo> items)
    {
        string s = string.Join(";", items.OrderBy(i => i.Sequence, StringComparer.Ordinal).Select(i =>
            $"{i.Sequence}:{(i.Positions.Count == 0 ? "-" : string.Join(",", i.Positions))}:{(i.ReversePositions.Count == 0 ? "-" : string.Join(",", i.ReversePositions))}"));
        return s.Length == 0 ? "-" : s;
    }

    private static string Format(IEnumerable<RepeatPair> pairs)
    {
        string s = string.Join(";", pairs.OrderBy(p => p.FirstPosition).ThenBy(p => p.SecondPosition).ThenBy(p => p.Length).ThenBy(p => p.IsReverseComplement).Select(p => $"{p.FirstPosition}:{p.SecondPosition}:{p.Length}:{(p.IsReverseComplement ? 'r' : 'f')}"));
        return s.Length == 0 ? "-" : s;
    }

    private static string FormatSorted(string s)
        => s == "-" ? s : string.Join(";", s.Split(';').Select(h => h.Split(':')).OrderBy(p => int.Parse(p[0])).ThenBy(p => int.Parse(p[1])).ThenBy(p => int.Parse(p[2]))
            .ThenBy(p => p[3]).Select(p => string.Join(":", p)));

    [Test]
    public void FindRepeatsBothStrands_Maximal_MatchesOracle_AllRows()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(1900));
        int nonEmpty = 0;
        foreach (var row in rows)
        {
            string actual = Format(GenomicAnalyzer.FindRepeatsBothStrands(row.Sequence, row.MinLength, maximalOnly: true));
            Assert.That(actual, Is.EqualTo(row.Maximal), $"{row.Sequence} L={row.MinLength}");
            if (actual != "-") nonEmpty++;
        }
        Assert.That(nonEmpty, Is.GreaterThan(1500));
    }

    [Test]
    public void FindRepeatsBothStrands_AllSubstrings_MatchesOracle()
    {
        int checkedRows = 0;
        foreach (var row in ReadOracle().Where(r => r.All != "?"))
        {
            Assert.That(Format(GenomicAnalyzer.FindRepeatsBothStrands(row.Sequence, row.MinLength)), Is.EqualTo(row.All), $"{row.Sequence} L={row.MinLength}");
            checkedRows++;
        }
        Assert.That(checkedRows, Is.GreaterThan(400));
    }

    [Test]
    public void DnaSequenceOverload_EqualsStringOverload_OnAcgtRows()
    {
        foreach (var row in ReadOracle().Where(r => r.Sequence.All(c => "ACGT".Contains(c))))
        {
            var dna = new DnaSequence(row.Sequence);
            Assert.That(Format(GenomicAnalyzer.FindRepeatsBothStrands(dna, row.MinLength, maximalOnly: true)), Is.EqualTo(row.Maximal), row.Sequence);
        }
    }

    [Test]
    public void FindRepeatPairs_MatchesRepeatMatch_BothStrandsAndForwardOnly()
    {
        int both = 0, forward = 0;
        foreach (var row in ReadOracle())
        {
            if (row.Pairs != "-")
            {
                Assert.That(Format(GenomicAnalyzer.FindRepeatPairs(row.Sequence, row.MinLength)), Is.EqualTo(FormatSorted(row.Pairs)), $"{row.Sequence} L={row.MinLength}");
                both++;
            }

            if (row.ForwardPairs != "-")
            {
                Assert.That(Format(GenomicAnalyzer.FindRepeatPairs(row.Sequence, row.MinLength, StrandSelection.Forward)), Is.EqualTo(FormatSorted(row.ForwardPairs)), $"{row.Sequence} L={row.MinLength} -f");
                forward++;
            }
        }
        Assert.That(both, Is.GreaterThan(400));
        Assert.That(forward, Is.GreaterThan(400));
    }

    [Test]
    public void MaximalClasses_AreTheStringsOfTheRepeatMatchPairs()
    {
        // the class set of repeat-match's pairs == the maximal both-strand repeat classes (the generator asserts it on the oracle side too)
        foreach (var row in ReadOracle().Where(r => r.Pairs != "-"))
        {
            var fromPairs = GenomicAnalyzer.FindRepeatPairs(row.Sequence, row.MinLength)
                .Select(p => row.Sequence.Substring(p.FirstPosition, p.Length))
                .Select(a => string.CompareOrdinal(a, DnaSequence.GetReverseComplementString(a)) <= 0 ? a : DnaSequence.GetReverseComplementString(a))
                .Distinct().OrderBy(a => a, StringComparer.Ordinal).ToList();
            var classes = GenomicAnalyzer.FindRepeatsBothStrands(row.Sequence, row.MinLength, maximalOnly: true).Select(c => c.Sequence).OrderBy(a => a, StringComparer.Ordinal).ToList();
            Assert.That(classes, Is.EqualTo(fromPairs), $"{row.Sequence} L={row.MinLength}");
        }
    }

    [Test]
    public void FindAllLongestRepeatsBothStrands_AreTheLongestMaximalClasses()
    {
        foreach (var row in ReadOracle().Where(r => r.MinLength == 1))
        {
            var maximal = GenomicAnalyzer.FindRepeatsBothStrands(row.Sequence, 1, maximalOnly: true);
            int max = maximal.Count == 0 ? 0 : maximal.Max(m => m.Length);
            var expected = Format(maximal.Where(m => m.Length == max));
            var all = GenomicAnalyzer.FindAllLongestRepeatsBothStrands(row.Sequence);
            Assert.That(Format(all), Is.EqualTo(expected), row.Sequence);
            if (row.Sequence.All(c => "ACGT".Contains(c)))
                Assert.That(GenomicAnalyzer.FindLongestRepeatBothStrands(new DnaSequence(row.Sequence)).Sequence, Is.EqualTo(all.Count == 0 ? string.Empty : all[0].Sequence), row.Sequence);
        }
    }

    [Test]
    public void Palindrome_OneSite_IsARepeatOnBothStrands()
    {
        var r = GenomicAnalyzer.FindRepeatsBothStrands(new DnaSequence("ACGAATTCCA"), 6, maximalOnly: true).ToList();
        Assert.That(r.Select(x => (x.Sequence, x.Count)), Has.Member(("GAATTC", 2)));
        var g = r.Single(x => x.Sequence == "GAATTC");
        Assert.That(g.Positions, Is.EqualTo(new[] { 2 }));
        Assert.That(g.ReversePositions, Is.EqualTo(new[] { 2 }));
        // ... but not a forward repeat
        Assert.That(GenomicAnalyzer.FindRepeats(new DnaSequence("ACGAATTCCA"), 6), Is.Empty);
    }

    [Test]
    public void ReverseComplementCopy_IsFoundWithBothPositionLists()
    {
        // ACGGTTC ... GAACCGT : the second copy is the reverse complement
        var r = GenomicAnalyzer.FindRepeatsBothStrands(new DnaSequence("ACGGTTCAAAAGAACCGT"), 7, maximalOnly: true).Single();
        Assert.That(r.Sequence, Is.EqualTo("ACGGTTC"));
        Assert.That(r.Positions, Is.EqualTo(new[] { 0 }));
        Assert.That(r.ReversePositions, Is.EqualTo(new[] { 11 }));
    }

    #region Single-strand string overloads (A2-35): brute force with N runs

    private static bool IsAcgt(char c) => c is 'A' or 'C' or 'G' or 'T';

    // Brute force: distinct all-ACGT substrings of length >= L occurring >= 2 times; maximal = left- and right-diverse (a non-ACGT neighbour is its own value).
    private static Dictionary<string, List<int>> Brute(string s, int minLength, bool maximalOnly)
    {
        var result = new Dictionary<string, List<int>>();
        for (int l = Math.Max(1, minLength); l <= s.Length; l++)
        {
            var occ = new Dictionary<string, List<int>>();
            for (int i = 0; i + l <= s.Length; i++)
            {
                string w = s.Substring(i, l);
                if (!w.All(IsAcgt)) continue;
                if (!occ.TryGetValue(w, out var list)) occ[w] = list = [];
                list.Add(i);
            }

            foreach (var (w, list) in occ)
            {
                if (list.Count < 2) continue;
                if (maximalOnly)
                {
                    var left = list.Select(i => i > 0 && IsAcgt(s[i - 1]) ? s[i - 1].ToString() : "#L" + i).Distinct().Count();
                    var right = list.Select(i => i + l < s.Length && IsAcgt(s[i + l]) ? s[i + l].ToString() : "#R" + i).Distinct().Count();
                    if (left < 2 || right < 2) continue;
                }

                result[w] = list;
            }
        }

        return result;
    }

    private static string Canon(IEnumerable<RepeatInfo> repeats)
        => string.Join(";", repeats.OrderBy(r => r.Sequence, StringComparer.Ordinal).Select(r => $"{r.Sequence}:{string.Join(",", r.Positions)}"));

    private static string Canon(Dictionary<string, List<int>> d)
        => string.Join(";", d.OrderBy(kv => kv.Key, StringComparer.Ordinal).Select(kv => $"{kv.Key}:{string.Join(",", kv.Value)}"));

    private static IEnumerable<string> Corpus()
    {
        var rnd = new Random(20261120);
        for (int iter = 0; iter < 600; iter++)
        {
            string alpha = new[] { "ACGT", "AC", "ACGTN", "ACGTNR", "acgtn" }[rnd.Next(5)];
            int n = rnd.Next(0, 34);
            yield return new string(Enumerable.Range(0, n).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray());
        }

        yield return "ACGTNACGTNACGT";
        yield return "NNNNNNNN";
        yield return "ACNNNACNNNAC";
        yield return "ACGTACGTNNNNACGTACGT";
        yield return "ACGTUACGTU";
    }

    [Test]
    public void StringOverloads_MatchBruteForce_WithAmbiguitySymbols()
    {
        int withN = 0, nonEmpty = 0;
        foreach (string s in Corpus())
        {
            string upper = s.ToUpperInvariant();
            if (upper.Any(c => !IsAcgt(c))) withN++;
            foreach (int min in new[] { 1, 2, 3 })
            {
                var all = Brute(upper, min, maximalOnly: false);
                var maximal = Brute(upper, min, maximalOnly: true);
                Assert.That(Canon(GenomicAnalyzer.FindRepeats(s, min)), Is.EqualTo(Canon(all)), $"{s} min={min} all");
                Assert.That(Canon(GenomicAnalyzer.FindMaximalRepeats(s, min)), Is.EqualTo(Canon(maximal)), $"{s} min={min} maximal");
                if (all.Count > 0) nonEmpty++;
            }

            var maxAll = Brute(upper, 1, maximalOnly: false);
            int longest = maxAll.Count == 0 ? 0 : maxAll.Keys.Max(k => k.Length);
            var expectedLongest = maxAll.Where(kv => kv.Key.Length == longest).OrderBy(kv => kv.Value[0]).Select(kv => $"{kv.Key}:{string.Join(",", kv.Value)}");
            var actualLongest = GenomicAnalyzer.FindAllLongestRepeats(s).Select(r => $"{r.Sequence}:{string.Join(",", r.Positions)}");
            Assert.That(actualLongest, Is.EqualTo(expectedLongest), $"{s} longest");
            Assert.That(GenomicAnalyzer.FindLongestRepeat(s).Sequence, Is.EqualTo(longest == 0 ? string.Empty : expectedLongest.First().Split(':')[0]), $"{s} first longest");
        }

        Assert.That(withN, Is.GreaterThan(300));
        Assert.That(nonEmpty, Is.GreaterThan(600));
    }

    [Test]
    public void StringOverloads_EqualDnaSequenceOverloads_OnAcgt()
    {
        var rnd = new Random(11);
        for (int iter = 0; iter < 200; iter++)
        {
            string s = new string(Enumerable.Range(0, rnd.Next(1, 60)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var dna = new DnaSequence(s);
            foreach (int min in new[] { 1, 3 })
            {
                Assert.That(GenomicAnalyzer.FindRepeats(s, min).Select(r => (r.Sequence, string.Join(",", r.Positions))),
                    Is.EquivalentTo(GenomicAnalyzer.FindRepeats(dna, min).Select(r => (r.Sequence, string.Join(",", r.Positions)))), s);
                Assert.That(GenomicAnalyzer.FindMaximalRepeats(s.ToLowerInvariant(), min).Select(r => (r.Sequence, string.Join(",", r.Positions))),
                    Is.EquivalentTo(GenomicAnalyzer.FindMaximalRepeats(dna, min).Select(r => (r.Sequence, string.Join(",", r.Positions)))), s);
            }

            Assert.That(GenomicAnalyzer.FindAllLongestRepeats(s).Select(r => (r.Sequence, string.Join(",", r.Positions))),
                Is.EqualTo(GenomicAnalyzer.FindAllLongestRepeats(dna).Select(r => (r.Sequence, string.Join(",", r.Positions)))), s);
        }
    }

    [Test]
    public void NRun_IsNeverARepeat_AndAnAdjacentRepeatIsUnchanged()
    {
        Assert.That(GenomicAnalyzer.FindRepeats("NNNNNNNN", 1), Is.Empty);
        Assert.That(GenomicAnalyzer.FindLongestRepeat("NNNNNNNN").IsEmpty, Is.True);
        Assert.That(GenomicAnalyzer.FindMaximalRepeats("ACGTNACGTN", 1).Select(r => r.Sequence), Is.EqualTo(new[] { "ACGT" }));
        // a protein: an explicit alphabet makes the engine alphabet-agnostic
        var protein = GenomicAnalyzer.FindMaximalRepeats("MKVLAAGMKVLTT", 3, "ACDEFGHIKLMNPQRSTVWY").Select(r => (r.Sequence, string.Join(",", r.Positions))).ToList();
        Assert.That(protein, Is.EqualTo(new[] { ("MKVL", "0,7") }));
    }

    [Test]
    public void StringOverloads_Validation()
    {
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindRepeats((string)null!, 2));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindLongestRepeat((string)null!));
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindMaximalRepeats("ACGT", 2, "AC1"));
        Assert.That(GenomicAnalyzer.FindRepeats(string.Empty, 1), Is.Empty);
        Assert.That(GenomicAnalyzer.FindAllLongestRepeats(string.Empty), Is.Empty);
        Assert.That(GenomicAnalyzer.FindRepeatsBothStrands(string.Empty, 1), Is.Empty);
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindRepeatPairs("ACGT", 1));
    }

    #endregion
}
