// GENOMIC-MOTIFS-001 (B09 audit A3-12a, F146) — fuzznuc with an IUPAC ambiguity code IN THE SEQUENCE (GenomicAnalyzer.FindPatternHits(..., sequenceAmbiguity: true)).
// Oracle: scripts/analysis/generate_fuzznuc_iupac_sequence_oracle.py — 1350 (sequence, pattern, k) rows, EMBOSS fuzznuc 6.6.0 (-pmismatch k -complement Y) == an independent Python reference of the
// measured rules (asserted by the generator). Worked values below were run on the binary by hand. Evidence: docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_PatternHits_IupacSequence_Tests
{
    private sealed record Row(string Sequence, string Pattern, int K, string Forward, string Reverse);

    private static readonly Lazy<List<Row>> OracleRows = new(ReadOracleRows);

    private static List<Row> ReadOracleRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.fuzznuc_iupac_sequence_oracle.tsv";
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

    [Test]
    public void FindPatternHits_SequenceAmbiguity_MatchesFuzznuc_AllRows()
    {
        var rows = OracleRows.Value;
        Assert.That(rows, Has.Count.EqualTo(1350));
        int nonEmpty = 0, withAmbiguity = 0, mismatches = 0, ranges = 0, reverse = 0, classes = 0, exclusions = 0, anchored = 0, x = 0, u = 0;
        foreach (var r in rows)
        {
            var hits = GenomicAnalyzer.FindPatternHits(r.Sequence, r.Pattern, r.K, bothStrands: true, sequenceAmbiguity: true);
            Assert.That(Fmt(hits.Where(h => !h.IsReverseStrand)), Is.EqualTo(r.Forward), $"{r.Sequence} {r.Pattern} k={r.K} forward");
            Assert.That(Fmt(hits.Where(h => h.IsReverseStrand)), Is.EqualTo(r.Reverse), $"{r.Sequence} {r.Pattern} k={r.K} reverse");
            if (hits.Count > 0) nonEmpty++;
            if (r.Sequence.IndexOfAny("RYSWKMBDHVN".ToCharArray()) >= 0) withAmbiguity++;
            if (r.K > 0) mismatches++;
            if (r.Pattern.Contains(',')) ranges++;
            if (r.Reverse != "-") reverse++;
            if (r.Pattern.Contains('[')) classes++;
            if (r.Pattern.Contains('{')) exclusions++;
            if (r.Pattern.Contains('<') || r.Pattern.Contains('>')) anchored++;
            if (r.Sequence.Contains('X')) x++;
            if (r.Sequence.Contains('U')) u++;
        }

        Assert.That(nonEmpty, Is.GreaterThan(1100));
        Assert.That(withAmbiguity, Is.GreaterThan(1300));
        Assert.That(mismatches, Is.GreaterThan(400));
        Assert.That(ranges, Is.GreaterThan(400));
        Assert.That(reverse, Is.GreaterThan(800));
        Assert.That(classes, Is.GreaterThan(400));
        Assert.That(exclusions, Is.GreaterThan(300));
        Assert.That(anchored, Is.GreaterThan(250));
        Assert.That(x, Is.GreaterThan(30));
        Assert.That(u, Is.GreaterThan(30));
    }

    private static string[] Hits(string seq, string pattern, int k = 0, bool both = false, bool ambiguity = true)
        => GenomicAnalyzer.FindPatternHits(seq, pattern, k, both, sequenceAmbiguity: ambiguity).Select(h => $"{(h.IsReverseStrand ? '-' : '+')}{h.Position}:{h.Length}:{h.Mismatches}").ToArray();

    [Test]
    public void WorkedValues_RunOnTheBinary()
    {
        // fuzznuc -pattern X -complement Y on ARGRYN (1-based GFF starts of single bases)
        Assert.That(Hits("ARGRYN", "A", both: true), Is.EqualTo(new[] { "+0:1:0" }));
        Assert.That(Hits("ARGRYN", "R", both: true), Is.EquivalentTo(new[] { "+0:1:0", "+1:1:0", "+2:1:0", "+3:1:0", "-4:1:0" }));
        Assert.That(Hits("ARGRYN", "[AG]", both: true), Is.EquivalentTo(new[] { "+0:1:0", "+2:1:0" }));
        Assert.That(Hits("ARGRYN", "{A}"), Is.EquivalentTo(new[] { "+1:1:0", "+2:1:0", "+3:1:0", "+4:1:0", "+5:1:0" }));
        Assert.That(Hits("ARGRYN", "N"), Has.Length.EqualTo(6));
        Assert.That(Hits("ARGRYN", "M", both: true), Is.EquivalentTo(new[] { "+0:1:0", "-2:1:0" }));
        Assert.That(Hits("ARGRYN", "V", both: true), Is.EquivalentTo(new[] { "+0:1:0", "+1:1:0", "+2:1:0", "+3:1:0", "-2:1:0", "-4:1:0" }));
        // an exclusion {ACGT} matches exactly the symbols that are not plain bases; U reads as T
        Assert.That(Hits("ACXGT", "{ACGT}", both: true), Is.EquivalentTo(new[] { "+2:1:0", "-2:1:0" }));
        Assert.That(Hits("ACUGT", "T"), Is.EquivalentTo(new[] { "+2:1:0", "+4:1:0" }));
    }

    [Test]
    public void Default_KeepsTheLiteralContract_AmbiguityMatchesNothing()
    {
        Assert.That(Hits("ARGRYN", "R", ambiguity: false), Is.EqualTo(new[] { "+0:1:0", "+2:1:0" }));   // only the plain A and G
        Assert.That(Hits("ARGRYN", "N", ambiguity: false), Is.EqualTo(new[] { "+0:1:0", "+2:1:0" }));
        // an A/C/G/T string gives the same hits in both readings
        foreach (var r in OracleRows.Value.Where(r => r.Sequence.All(c => "ACGT".Contains(c))).Take(100))
            Assert.That(GenomicAnalyzer.FindPatternHits(r.Sequence, r.Pattern, r.K, true, sequenceAmbiguity: true),
                Is.EqualTo(GenomicAnalyzer.FindPatternHits(r.Sequence, r.Pattern, r.K, true)), $"{r.Sequence} {r.Pattern}");
    }

    [Test]
    public void NamedPatterns_AndParsedPattern_TakeTheOption()
    {
        var parsed = FuzznucPattern.Parse("R");
        Assert.That(GenomicAnalyzer.FindPatternHits("ARGRYN", parsed, 0, false, default, true), Has.Count.EqualTo(4));
        var named = GenomicAnalyzer.FindPatternHits("ARGRYN", FuzznucPattern.ParsePatternFile(">p1\nR\n>p2\nY\n"), 0, false, default, true);
        Assert.That(named["p1"], Has.Count.EqualTo(4));
        Assert.That(named["p2"].Select(h => h.Position), Is.EqualTo(new[] { 4 }));
    }

    [Test]
    public void ExpandPaths_KeepTheKindOfEachElement()
    {
        var p = FuzznucPattern.Parse("R[AG]{AG}N(0,2)");
        Assert.That(p.Elements.Select(e => e.Kind), Is.EqualTo(new[] { PatternElementKind.Letter, PatternElementKind.Class, PatternElementKind.Exclusion, PatternElementKind.Letter }));
        var paths = p.ExpandPaths();
        Assert.That(paths.Select(x => string.Join(",", x.Elements)), Is.EqualTo(new[] { "0,1,2", "0,1,2,3", "0,1,2,3,3" }));
        Assert.That(paths.Select(x => x.EndsWithEmptyElement), Is.EqualTo(new[] { true, false, false }));
        Assert.That(string.Join(";", paths[2].ExtraPositions), Is.EqualTo("3;4"));
        // plain expansion unchanged: a class [ACGT] and the letter N have the same base set
        Assert.That(FuzznucPattern.Parse("[ACGT]N").Expand(), Is.EqualTo(new[] { "NN" }));
        Assert.That(p.Matches(0, 'R'), Is.True);
        Assert.That(p.Matches(1, 'R'), Is.False);
        Assert.That(p.Matches(2, 'R'), Is.True);
        Assert.That(p.Matches(3, 'X'), Is.True);
    }
}
