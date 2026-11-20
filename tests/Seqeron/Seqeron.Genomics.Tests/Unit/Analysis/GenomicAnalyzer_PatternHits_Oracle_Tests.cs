// GENOMIC-MOTIFS-001 (B09 A2-32 / A2-35) — the EMBOSS fuzznuc pattern language: FuzznucPattern (parser / expansion) and GenomicAnalyzer.FindPatternHits.
// Oracle: scripts/analysis/generate_fuzznuc_pattern_oracle.py — 3300 (sequence, pattern, k) rows, EMBOSS fuzznuc 6.6.0 (-pmismatch k -complement Y) == an independent Python
// reference (asserted by the generator; includes the two fuzznuc range-search classes read in embpat.c). Worked values below were run on the binary by hand.
// Evidence: docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_PatternHits_Oracle_Tests
{
    private sealed record Row(string Sequence, string Pattern, int K, string Forward, string Reverse);

    // One parsed copy per fixture (audit A3-15): the oracle was re-read and re-parsed on every call, i.e. in every test.
    private static readonly Lazy<List<Row>> OracleRows = new(ReadOracleRows);

    private static List<Row> ReadOracle() => OracleRows.Value;

    private static List<Row> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.fuzznuc_pattern_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            rows.Add(new Row(f[0], f[1], int.Parse(f[2]), f[3], f[4]));
        }
        return rows;
    }

    private static string Fmt(IEnumerable<PatternHit> hits)
    {
        string s = string.Join(" ", hits.OrderBy(h => h.Position).ThenBy(h => h.Length).Select(h => $"{h.Position}:{h.Length}:{h.Mismatches}"));
        return s.Length == 0 ? "-" : s;
    }

    #region fuzznuc oracle

    [Test]
    public void FindPatternHits_BothStrands_MatchesFuzznuc_AllRows()
    {
        var rows = ReadOracle();
        Assert.That(rows, Has.Count.EqualTo(3300));
        int nonEmpty = 0, ranges = 0, mismatches = 0, anchored = 0, reverse = 0;
        foreach (var r in rows)
        {
            var hits = GenomicAnalyzer.FindPatternHits(r.Sequence, r.Pattern, r.K, bothStrands: true);
            Assert.That(Fmt(hits.Where(h => !h.IsReverseStrand)), Is.EqualTo(r.Forward), $"{r.Sequence} {r.Pattern} k={r.K} forward");
            Assert.That(Fmt(hits.Where(h => h.IsReverseStrand)), Is.EqualTo(r.Reverse), $"{r.Sequence} {r.Pattern} k={r.K} reverse");
            if (hits.Count > 0) nonEmpty++;
            if (r.Pattern.Contains(',')) ranges++;
            if (r.K > 0) mismatches++;
            if (r.Pattern.Contains('<') || r.Pattern.Contains('>')) anchored++;
            if (r.Reverse != "-") reverse++;
        }

        Assert.That(nonEmpty, Is.GreaterThan(2700));
        Assert.That(ranges, Is.GreaterThan(1500));
        Assert.That(mismatches, Is.GreaterThan(1000));
        Assert.That(anchored, Is.GreaterThan(800));
        Assert.That(reverse, Is.GreaterThan(2000));
    }

    [Test]
    public void FindPatternHits_ForwardOnly_AndDnaSequenceOverload()
    {
        foreach (var r in ReadOracle().Take(1200))
        {
            var dna = new DnaSequence(r.Sequence);
            var forward = GenomicAnalyzer.FindPatternHits(dna, r.Pattern, r.K);
            Assert.That(forward.All(h => !h.IsReverseStrand), Is.True);
            Assert.That(Fmt(forward), Is.EqualTo(r.Forward), $"{r.Sequence} {r.Pattern} k={r.K}");
            Assert.That(GenomicAnalyzer.FindPatternHits(dna, FuzznucPattern.Parse(r.Pattern), r.K, true).Count(h => h.IsReverseStrand),
                Is.EqualTo(r.Reverse == "-" ? 0 : r.Reverse.Split(' ').Length));
        }
    }

    #endregion

    #region Values run on the fuzznuc 6.6.0 binary

    private static string[] Hits(string seq, string pattern, int k = 0, bool both = false)
        => GenomicAnalyzer.FindPatternHits(seq, pattern, k, both).Select(h => $"{(h.IsReverseStrand ? '-' : '+')}{h.Position}:{h.Length}:{h.Mismatches}").ToArray();

    [Test]
    public void Range_ReportsEveryMatchingWindow()
    {
        // fuzznuc -pattern "A(1,3)C" on AAAACAAAAT -> 1-based 4..5, 3..5, 2..5
        Assert.That(Hits("AAAACAAAAT", "A(1,3)C"), Is.EquivalentTo(new[] { "+3:2:0", "+2:3:0", "+1:4:0" }));
        // fuzznuc -pattern "CA(0,2)G" on ACGTACGT -> 2..3 and 6..7 (A is optional)
        Assert.That(Hits("ACGTACGT", "CA(0,2)G"), Is.EquivalentTo(new[] { "+1:2:0", "+5:2:0" }));
    }

    [Test]
    public void RangeWithMismatches_ExtraRepeatsMustMatchExactly()
    {
        // fuzznuc -pattern "A(1,3)C" -pmismatch 1 on AAGACCAAC (1-based: 1-2:1 1-3:1 2-3:1 3-5:1 4-5:0 5-6:1 6-9:1 7-8:1 7-9:0 8-9:0);
        // note: 4..6 ('ACC' against the A-A-C branch) is absent: the second A is an extra repeat and may not mismatch.
        Assert.That(Hits("AAGACCAAC", "A(1,3)C", 1), Is.EquivalentTo(new[]
        {
            "+0:2:1", "+0:3:1", "+1:2:1", "+2:3:1", "+3:2:0", "+4:2:1", "+5:4:1", "+6:2:1", "+6:3:0", "+7:2:0",
        }));
    }

    [Test]
    public void Anchors_FollowTheStrand()
    {
        // fuzznuc -complement Y: "<ACG" on AAACGT -> '-' 4..6 ; on CGTAAA nothing; "ACG>" on CGTAAA -> '-' 1..3
        Assert.That(Hits("AAACGT", "<ACG", 0, both: true), Is.EqualTo(new[] { "-3:3:0" }));
        Assert.That(Hits("CGTAAA", "<ACG", 0, both: true), Is.Empty);
        Assert.That(Hits("CGTAAA", "ACG>", 0, both: true), Is.EqualTo(new[] { "-0:3:0" }));
        Assert.That(Hits("ACGTTTACGT", "<ACG"), Is.EqualTo(new[] { "+0:3:0" }));
        Assert.That(Hits("ACGTTTACGT", "CGT>"), Is.EqualTo(new[] { "+7:3:0" }));
    }

    [Test]
    public void StartAnchoredExactRange_ReturnsOnlyTheLongestMatchAtTheAnchor()
    {
        // fuzznuc -pattern "<A(1,3)" on AAAACAAAAT reports only 1..3 (regex search), whereas "A(1,3)" lists every window at position 0
        Assert.That(Hits("AAAACAAAAT", "<A(1,3)"), Is.EqualTo(new[] { "+0:3:0" }));
        Assert.That(Hits("AAAACAAAAT", "A(1,3)").Where(h => h.StartsWith("+0:")), Is.EquivalentTo(new[] { "+0:1:0", "+0:2:0", "+0:3:0" }));
        // with mismatches the recursive matcher lists all lengths at the anchor
        Assert.That(Hits("AAAACAAAAT", "<A(1,3)", 0).Length, Is.EqualTo(1));
    }

    [Test]
    public void ClassesExclusionsAndRepeats()
    {
        Assert.That(Hits("ACGTACGT", "N(3)"), Is.EqualTo(new[] { "+0:3:0", "+1:3:0", "+2:3:0", "+3:3:0", "+4:3:0", "+5:3:0" }));
        Assert.That(Hits("ACGTACGT", "{A}(2)"), Is.EqualTo(new[] { "+1:2:0", "+2:2:0", "+5:2:0", "+6:2:0" }));
        Assert.That(Hits("ACGTACGT", "[AC](2)"), Is.EqualTo(new[] { "+0:2:0", "+4:2:0" }));
        Assert.That(Hits("ACGTACGT", "{ACGT}"), Is.Empty);
        Assert.That(Hits("ACGTACGT", "a-c-g."), Is.EqualTo(new[] { "+0:3:0", "+4:3:0" }));
        Assert.That(Hits("ACGTACGT", "RYK"), Is.EqualTo(new[] { "+0:3:0", "+4:3:0" }));
    }

    #endregion

    #region Parser

    [TestCase("[CG](5)TG{A}N(1,5)C")]
    [TestCase("<ACG>")]
    [TestCase("A-C-G.")]
    [TestCase("AC(2,3)G")]
    [TestCase("  [ac]{gt}(2)  ")]
    [TestCase("N(1,1)")]
    public void Parse_AcceptsTheDocumentedLanguage(string pattern)
        => Assert.That(FuzznucPattern.TryParse(pattern, out var parsed, out var error), Is.True, error);

    [TestCase("")]
    [TestCase("[]")]
    [TestCase("{}")]
    [TestCase("[AC")]
    [TestCase("AC]")]
    [TestCase("{A")]
    [TestCase("A(2,1)")]
    [TestCase("A(0)")]
    [TestCase("A(3,)")]
    [TestCase("A()")]
    [TestCase("A(1,2,3)")]
    [TestCase("A(x)")]
    [TestCase("A(2")]
    [TestCase("A<C")]
    [TestCase("<A>C")]
    [TestCase("[RT]")]
    [TestCase("{N}")]
    [TestCase("U")]
    [TestCase("ACGU")]
    [TestCase("A(0,0)")]
    [TestCase("<")]
    [TestCase(">")]
    [TestCase(".")]
    [TestCase("A(1)(2)")]
    public void Parse_RejectsWhatTheLanguageDoesNotDefine(string pattern)
    {
        Assert.That(FuzznucPattern.TryParse(pattern, out var parsed, out var error), Is.False);
        Assert.That(parsed, Is.Null);
        Assert.That(error, Is.Not.Empty);
        Assert.Throws<ArgumentException>(() => FuzznucPattern.Parse(pattern));
    }

    [Test]
    public void Parse_NormalisesAndExpands()
    {
        var p = FuzznucPattern.Parse("<ac-[TG](2)-{A}N(1,2)>.");
        Assert.That(p.ToString(), Is.EqualTo("<AC[GT](2)[CGT]N(1,2)>"));
        Assert.That((p.AnchoredAtStart, p.AnchoredAtEnd, p.MinLength, p.MaxLength, p.IsFixedLength), Is.EqualTo((true, true, 6, 7, false)));
        Assert.That(FuzznucPattern.Parse("AC(2,3)G").Expand(), Is.EqualTo(new[] { "ACCG", "ACCCG" }));
        Assert.That(FuzznucPattern.Parse("CA(0,2)G").Expand(), Is.EqualTo(new[] { "CG", "CAG", "CAAG" }));
        Assert.That(FuzznucPattern.Parse("[AC]{AG}").Expand(), Is.EqualTo(new[] { "MY" }));
        Assert.That(FuzznucPattern.Parse("{ACGT}A").Expand(), Is.Empty);
        // A(1,2)A(1,2): two distinct paths to AAA with different extras
        var e = FuzznucPattern.Parse("A(1,2)A(1,2)").ExpandWithExtras();
        Assert.That(e.Select(x => (x.Iupac, string.Join(",", x.ExtraPositions))), Is.EqualTo(new[] { ("AA", ""), ("AAA", "1"), ("AAA", "2"), ("AAAA", "1,3") }));
        Assert.Throws<ArgumentException>(() => FuzznucPattern.Parse("N(0,400)N(0,400)").Expand());
    }

    [Test]
    public void ReverseComplement_ReversesElementsComplementsSetsAndSwapsAnchors()
    {
        var p = FuzznucPattern.Parse("<[AC]G(2)T(1,3)>");
        Assert.That(p.ReverseComplement().ToString(), Is.EqualTo("<A(1,3)C(2)[GT]>"));
        Assert.That(p.ReverseComplement().ReverseComplement().ToString(), Is.EqualTo(p.ToString()));
        // every reverse expansion is the reverse complement of a forward one
        foreach (string text in new[] { "AC(2,3)G", "[AC]{AG}N(0,2)T", "RYSWKM(1,2)BDHV" })
        {
            var q = FuzznucPattern.Parse(text);
            Assert.That(q.ReverseComplement().Expand().ToHashSet(),
                Is.EquivalentTo(q.Expand().Select(DnaSequence.GetReverseComplementString).ToHashSet()), text);
        }
    }

    [Test]
    public void ParsePatternFile_PlainAndFasta()
    {
        var plain = FuzznucPattern.ParsePatternFile("ACGT\n\nTTGA\nA(2,3)C\n");
        Assert.That(plain.Select(x => (x.Name, x.Pattern.ToString())), Is.EqualTo(new[] { ("pattern1", "ACGT"), ("pattern2", "TTGA"), ("pattern3", "AC(2,3)".Replace("AC(2,3)", "A(2,3)C")) }));
        var fasta = FuzznucPattern.ParsePatternFile(">p1 first\nACGT\n>p2\nTT\nGA\n>\nCC\n");
        Assert.That(fasta.Select(x => (x.Name, x.Pattern.ToString())), Is.EqualTo(new[] { ("p1", "ACGT"), ("p2", "TTGA"), ("pattern3", "CC") }));
        Assert.That(FuzznucPattern.ParsePatternFile("ACGT", "xx").Single().Name, Is.EqualTo("xx1"));
        Assert.That(FuzznucPattern.ParsePatternFile(" \n"), Is.Empty);
        Assert.Throws<ArgumentException>(() => FuzznucPattern.ParsePatternFile("ACGT\n[AC"));
    }

    [Test]
    public void NamedPatterns_AndValidation()
    {
        var named = FuzznucPattern.ParsePatternFile(">p1\nACGT\n>p2\nTTGA\n>p3\nCCCC\n");
        var result = GenomicAnalyzer.FindPatternHits("ACGTACGTTTGAACC", named);
        Assert.That(result.Keys, Is.EqualTo(new[] { "p1", "p2", "p3" }));
        Assert.That(result["p1"].Select(h => h.Position), Is.EqualTo(new[] { 0, 4 }));
        Assert.That(result["p2"].Select(h => h.Position), Is.EqualTo(new[] { 8 }));
        Assert.That(result["p3"], Is.Empty);
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindPatternHits("ACGT", new[] { named[0], named[0] }));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindPatternHits("ACGT", "AC", -1));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindPatternHits((string)null!, "AC"));
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindPatternHits("ACGT", "A(2,1)"));
    }

    [Test]
    public void StringSequences_WithN_NeverMatchNonAcgt()
    {
        // the canonical contract: a sequence symbol other than A/C/G/T matches no pattern symbol, N included, and costs a mismatch with k > 0
        Assert.That(Hits("ACNTAAAT", "ACNT"), Is.Empty);
        Assert.That(Hits("ACNTAAAT", "ACNT", 1), Is.EquivalentTo(new[] { "+0:4:1", "+4:4:1" })); // the N of the sequence is a mismatch; AAAT differs from ACNT at one place
        Assert.That(Hits("ACGTNNACGT", "NN"), Is.EqualTo(new[] { "+0:2:0", "+1:2:0", "+2:2:0", "+6:2:0", "+7:2:0", "+8:2:0" }));
        Assert.That(Hits("acgtnnacgt", "[AC]G"), Is.EqualTo(new[] { "+1:2:0", "+7:2:0" }));
    }

    #endregion
}
