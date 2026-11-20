// GENOMIC-COMMON-001 / GENOMIC-REPEAT-001 (B09 A3-13) — EMBOSS seqmatchall (GenomicAnalyzer.SeqMatchAll) and EMBOSS palindrome (gaplimit 0) vs GenomicAnalyzer.FindPalindromes.
// Oracles: scripts/analysis/generate_seqmatchall_oracle.py (330 sets, seqmatchall 6.6.0) and scripts/analysis/generate_palindrome_gap0_oracle.py (330 sequences, palindrome 6.6.0).
// Evidence: docs/Evidence/GENOMIC-COMMON-001-Evidence.md, docs/Evidence/GENOMIC-REPEAT-001-Evidence.md.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_EmbossPrograms_Oracle_Tests
{
    private static List<string[]> ReadRows(string file)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Analysis." + file;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<string[]>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0) rows.Add(line.Split('\t'));
        return rows;
    }

    [Test]
    public void SeqMatchAll_EqualsSeqmatchall_AllSets()
    {
        var rows = ReadRows("seqmatchall_oracle.tsv");
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(300));
        int nonEmpty = 0, protein = 0, lower = 0, withN = 0, pairs = 0;
        foreach (var f in rows)
        {
            var seqs = f[0].Split(',');
            int k = int.Parse(f[1]);
            var hits = GenomicAnalyzer.SeqMatchAll(seqs, k)
                .OrderBy(m => m.Sequence1).ThenBy(m => m.Sequence2).ThenBy(m => m.Position1).ThenBy(m => m.Position2).ThenBy(m => m.Length)
                .Select(m => $"{m.Sequence1}:{m.Sequence2}:{m.Position1}:{m.Position2}:{m.Length}");
            var expected = f[2] == "-" ? Array.Empty<string>() : f[2].Split(' ').OrderBy(x => x, StringComparer.Ordinal).ToArray();
            Assert.That(hits.OrderBy(x => x, StringComparer.Ordinal).ToArray(), Is.EqualTo(expected), $"k={k} {f[0]}");
            if (expected.Length > 0) nonEmpty++;
            if (f[0].Any(c => "EFHIKLMPQRSVWY".Contains(char.ToUpperInvariant(c)))) protein++;
            if (f[0].Any(char.IsLower)) lower++;
            if (f[0].Contains('N')) withN++;
            pairs += seqs.Length * (seqs.Length - 1) / 2;
        }

        Assert.That(nonEmpty, Is.GreaterThan(280));
        Assert.That(protein, Is.GreaterThan(50));
        Assert.That(lower, Is.GreaterThan(50));
        Assert.That(withN, Is.GreaterThan(20));
        Assert.That(pairs, Is.GreaterThan(600));
    }

    [Test]
    public void SeqMatchAll_OnNucleotides_IsTheMaximalExactMatchModeOfFindCommonRegions()
    {
        // one engine on N-free nucleotide sets: the canonical suffix-tree MEM stream (mummer -maxmatch), per ordered pair, length(seq1) > k <= length(seq2); the sets with N follow the program's table scan
        foreach (var f in ReadRows("seqmatchall_oracle.tsv").Where(r => !r[0].Any(c => "EFHIKLMPQRSVWYN".Contains(char.ToUpperInvariant(c)))).Take(200))
        {
            var seqs = f[0].Split(',');
            int k = int.Parse(f[1]);
            var viaCommon = new List<string>();
            for (int i = 0; i < seqs.Length; i++)
                for (int j = i + 1; j < seqs.Length; j++)
                {
                    if (seqs[i].Length <= k || seqs[j].Length < k) continue;
                    foreach (var r in GenomicAnalyzer.FindCommonRegions(seqs[i], seqs[j], k, CommonRegionMode.MaximalExactMatches))
                        viaCommon.Add($"{i}:{j}:{r.PositionInFirst}:{r.PositionInSecond}:{r.Length}");
                }

            var viaApi = GenomicAnalyzer.SeqMatchAll(seqs, k).Select(m => $"{m.Sequence1}:{m.Sequence2}:{m.Position1}:{m.Position2}:{m.Length}");
            Assert.That(viaApi.OrderBy(x => x, StringComparer.Ordinal).ToArray(), Is.EqualTo(viaCommon.OrderBy(x => x, StringComparer.Ordinal).ToArray()));
        }
    }

    [Test]
    public void SeqMatchAll_WorkedExample_AndArguments()
    {
        // seqmatchall -wordsize 6 (run by hand): "10 a 1..10 b 3..12", "6 a 1..6 b 7..12", "6 a 5..10 b 3..8" (the third sequence is shorter than 6 and compares with nothing)
        var r = GenomicAnalyzer.SeqMatchAll(["ACGTACGTAC", "GGACGTACGTACTT", "TTTT"], 6);
        Assert.That(r.Select(m => (m.Sequence1, m.Sequence2, m.Position1, m.Position2, m.Length)).OrderBy(x => x.Length).ThenBy(x => x.Position1).ToArray(),
            Is.EqualTo(new[] { (0, 1, 0, 6, 6), (0, 1, 4, 2, 6), (0, 1, 0, 2, 10) }));
        // N breaks every match of a nucleotide set, R / U match literally; X breaks a protein match, B / Z / U match (seqmatchall 6.6.0 runs)
        Assert.That(GenomicAnalyzer.SeqMatchAll(["ACGTACGNACGTACGT", "ACGTACGNACGTACGT"], 4).Select(m => m.Length).OrderDescending().ToArray(), Is.EqualTo(new[] { 8, 7, 7, 7, 4, 4, 4, 4 }));
        // ... but not always: after a skipped N the program's scan restarts right behind it and may keep a word with an N inside (oracle rows with 15-25 % N)
        Assert.That(GenomicAnalyzer.SeqMatchAll(["ACGTACGRACGTACGT", "ACGTACGRACGTACGT"], 4).Max(m => m.Length), Is.EqualTo(16));
        Assert.That(GenomicAnalyzer.SeqMatchAll(["MKVLAAXLLGGAMKVLAAXLLGGA", "MKVLAAXLLGGAMKVLAAXLLGGA"], 4).Max(m => m.Length), Is.EqualTo(11));
        Assert.That(GenomicAnalyzer.SeqMatchAll(["ACGT"], 2), Is.Empty);
        Assert.That(GenomicAnalyzer.SeqMatchAll(["ACGT", "ACGT"], 4), Is.Empty, "length(seq1) must exceed the word size");
        Assert.That(GenomicAnalyzer.SeqMatchAll(["ACGTA", "TACGT"], 4).Select(m => m.Length), Is.EqualTo(new[] { 4 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.SeqMatchAll(["ACGT", "ACGT"], 0));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.SeqMatchAll(null!, 4));
    }

    [Test]
    public void FindPalindromes_PerCentreLongest_EqualsEmbossPalindrome_NoLoopNoMismatch()
    {
        // EMBOSS palindrome (-gaplimit 0 -nummismatches 0 -overlap) reports, per centre, the longest reverse-complement palindrome of arm >= minpallen, capped to arm maxpallen
        // (the inner 2 * maxpallen bases when the palindrome is longer); FindPalindromes(2m, 2M) lists every even palindromic substring in that length range, nested ones included.
        var rows = ReadRows("palindrome_gap0_oracle.tsv");
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(300));
        int hitsTotal = 0, extras = 0, capped = 0;
        foreach (var f in rows)
        {
            string seq = f[0].ToUpperInvariant();
            int m = int.Parse(f[1]), bigM = int.Parse(f[2]);
            var all = GenomicAnalyzer.FindPalindromes(new DnaSequence(seq), 2 * m, 2 * bigM).ToList();
            var longest = all.GroupBy(p => p.Position + (p.Length / 2)).Select(g => g.OrderByDescending(p => p.Length).First())
                .Select(p => $"{p.Position + 1}:{p.Position + p.Length}:{p.Length / 2}").OrderBy(x => x, StringComparer.Ordinal).ToArray();
            var expected = f[3] == "-" ? Array.Empty<string>() : f[3].Split(' ').OrderBy(x => x, StringComparer.Ordinal).ToArray();
            Assert.That(longest, Is.EqualTo(expected), $"{seq} m={m} M={bigM}");
            hitsTotal += expected.Length;
            extras += all.Count - longest.Length;
            if (expected.Any(h => int.Parse(h.Split(':')[2]) == bigM)) capped++;
        }

        Assert.That(hitsTotal, Is.GreaterThan(200));
        Assert.That(extras, Is.GreaterThan(50), "FindPalindromes also lists the nested palindromes of a centre, which EMBOSS palindrome never reports");
        Assert.That(capped, Is.GreaterThan(10));
    }
}
