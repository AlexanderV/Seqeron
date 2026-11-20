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
}
