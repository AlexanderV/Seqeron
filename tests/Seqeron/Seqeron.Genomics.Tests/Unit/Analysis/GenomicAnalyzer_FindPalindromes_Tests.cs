// GENOMIC-REPEAT-001 (B09 F27) — GenomicAnalyzer.FindPalindromes is a facade over the canonical
// RepeatFinder.FindPalindromes (REP-PALIN-001). Evidence: docs/Evidence/GENOMIC-REPEAT-001-Evidence.md
// TestSpec: tests/TestSpecs/GENOMIC-REPEAT-001.md (section 8)
// Oracle: Biopython 1.88 Seq.reverse_complement over every window (== pure-str translate/reverse), 910 rows,
//         scripts/analysis/generate_palindrome_oracle.py (TestData/Analysis/palindrome_oracle.tsv); Rosalind REVP sample.

using System.Diagnostics;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_FindPalindromes_Tests
{
    private sealed record Row(string Sequence, int MinLength, int MaxLength, (int Position, int Length)[] Hits);

    private static List<Row> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.palindrome_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            var hits = f[3] == "-"
                ? Array.Empty<(int, int)>()
                : f[3].Split(',').Select(h => { var a = h.Split(':'); return (int.Parse(a[0]), int.Parse(a[1])); }).ToArray();
            rows.Add(new Row(f[0], int.Parse(f[1]), int.Parse(f[2]), hits));
        }
        return rows;
    }

    [Test]
    public void FindPalindromes_MatchesBiopythonOracle_PositionThenLengthOrder()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(910));
        int withHits = 0;
        foreach (var row in rows)
        {
            var actual = GenomicAnalyzer.FindPalindromes(new DnaSequence(row.Sequence), row.MinLength, row.MaxLength)
                .Select(p => (p.Position, p.Length)).ToArray();
            Assert.That(actual, Is.EqualTo(row.Hits), $"{row.Sequence} [{row.MinLength},{row.MaxLength}]");
            if (actual.Length > 0) withHits++;
        }
        Assert.That(withHits, Is.GreaterThan(500));
    }

    [Test]
    public void FindPalindromes_RosalindRevpSample_PositionThenLength()
    {
        // Rosalind REVP sample, 0-based: ATGCAT@3, TGCA@4, GCATGC@5, CATG@6, TATA@16, ATAT@17, ATGCAT@19, TGCA@20.
        var result = GenomicAnalyzer.FindPalindromes(new DnaSequence("TCAATGCATGCGGGTCTATATGCAT"), 4, 12).ToList();
        Assert.That(result.Select(p => (p.Position, p.Length)),
            Is.EqualTo(new[] { (3, 6), (4, 4), (5, 6), (6, 4), (16, 4), (17, 4), (19, 6), (20, 4) }));
        Assert.That(result.Select(p => p.Sequence),
            Is.EqualTo(new[] { "ATGCAT", "TGCA", "GCATGC", "CATG", "TATA", "ATAT", "ATGCAT", "TGCA" }));
    }

    [Test]
    public void FindPalindromes_EqualsTheCanonicalRepeatFinder()
    {
        var rng = new Random(11);
        for (int trial = 0; trial < 60; trial++)
        {
            string s = new(Enumerable.Range(0, rng.Next(0, 120)).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
            var dna = new DnaSequence(s);
            var canonical = RepeatFinder.FindPalindromes(dna, 4, 12).Select(p => (p.Sequence, p.Position)).ToArray();
            var facade = GenomicAnalyzer.FindPalindromes(dna).Select(p => (p.Sequence, p.Position)).ToArray();
            Assert.That(facade, Is.EqualTo(canonical), s);
        }
    }

    [Test]
    public void FindPalindromes_EcoRiSite_SmokeTest()
    {
        var palindromes = GenomicAnalyzer.FindPalindromes(new DnaSequence("AAAGAATTCAAA"), minLength: 6, maxLength: 6).ToList();
        Assert.That(palindromes, Has.Count.EqualTo(1));
        Assert.That(palindromes[0].Sequence, Is.EqualTo("GAATTC"));
        Assert.That(palindromes[0].Position, Is.EqualTo(3));
        Assert.That(palindromes[0].Length, Is.EqualTo(6));
    }

    #region Contract (validation was missing before the facade)

    [Test]
    public void FindPalindromes_NullSequence_ThrowsArgumentNullException_AtTheCall()
    {
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindPalindromes((DnaSequence)null!));
    }

    [TestCase(2)]
    [TestCase(3)]
    [TestCase(5)]
    [TestCase(0)]
    [TestCase(-4)]
    public void FindPalindromes_MinLengthOddOrBelowFour_ThrowsArgumentOutOfRange_AtTheCall(int minLength)
    {
        // Before: minLength 2 reported every dinucleotide AT/CG/GC/TA; odd minLength silently scanned only odd (impossible) lengths.
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindPalindromes(new DnaSequence("GAATTC"), minLength, 12));
    }

    [Test]
    public void FindPalindromes_MaxLengthBelowMinLength_ThrowsArgumentOutOfRange()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindPalindromes(new DnaSequence("GAATTC"), 8, 6));
    }

    [Test]
    public void FindPalindromes_OddMaxLength_IsRoundedDown_AndEmptyOrShortSequencesYieldNothing()
    {
        var dna = new DnaSequence("GCGCGC");
        Assert.That(GenomicAnalyzer.FindPalindromes(dna, 4, 5).Select(p => p.Length), Is.All.EqualTo(4));
        Assert.That(GenomicAnalyzer.FindPalindromes(new DnaSequence(""), 4, 12), Is.Empty);
        Assert.That(GenomicAnalyzer.FindPalindromes(new DnaSequence("GAA"), 4, 12), Is.Empty);
    }

    #endregion

    #region Complexity

    [Test]
    public void FindPalindromes_OneMegabase_IsLinearInTheLength_NotQuadraticInTheWindow()
    {
        // The old scan did one substring + reverse-complement string per window and length (O(n * L^2) character work);
        // the centre-expansion core is O(n * maxLength).
        var rng = new Random(5);
        var chars = new char[1_000_000];
        for (int i = 0; i < chars.Length; i++) chars[i] = "ACGT"[rng.Next(4)];
        var dna = new DnaSequence(new string(chars));

        long before = GC.GetAllocatedBytesForCurrentThread();
        var result = GenomicAnalyzer.FindPalindromes(dna, 4, 12).ToList();
        // A5-4: deterministic proxy of "linear, not one substring + reverse complement per window and length": the old scan made about 10 million strings (hundreds of MB), the result list is ~83k entries
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(128L << 20));
        // expected count of palindromic windows of length 2h in random DNA: n * 4^-h (4: 62.5k, 6: 15.6k, ... ) ~ 83k
        Assert.That(result.Count, Is.InRange(60_000, 110_000));
        Assert.That(result.Select(p => (p.Position, p.Length)), Is.Ordered);
    }

    #endregion
}
