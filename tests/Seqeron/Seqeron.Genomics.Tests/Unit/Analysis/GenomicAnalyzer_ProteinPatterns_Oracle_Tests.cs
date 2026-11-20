// GENOMIC-MOTIFS-001 (B09 A3-12) — EMBOSS fuzzpro / fuzztran: ProteinPattern, GenomicAnalyzer.FindProteinPatternHits and FindTranslatedPatternHits.
// Oracles: scripts/analysis/generate_fuzzpro_oracle.py (900 random syntax-tree patterns, fuzzpro 6.6.0) and scripts/analysis/generate_fuzztran_oracle.py (760 random DNA / pattern /
// frame / table cases, fuzztran 6.6.0). Evidence: docs/Evidence/GENOMIC-MOTIFS-001-Evidence.md.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_ProteinPatterns_Oracle_Tests
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
    public void FindProteinPatternHits_EqualsFuzzpro_AllRows()
    {
        var rows = ReadRows("fuzzpro_oracle.tsv");
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(800));
        int hitsRows = 0, rejected = 0, ranges = 0, anchored = 0, ambiguity = 0, lower = 0, wildcard = 0;
        foreach (var f in rows)
        {
            string label = $"{f[0]} / {f[1]}";
            if (f[2] == "ERR")
            {
                Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindProteinPatternHits(f[0], f[1]), label);
                rejected++;
                continue;
            }

            var hits = GenomicAnalyzer.FindProteinPatternHits(f[0], f[1]).Select(h => $"{h.Position}:{h.Length}").ToArray();
            var expected = f[2] == "-" ? Array.Empty<string>() : f[2].Split(' ');
            Assert.That(hits, Is.EqualTo(expected), label);
            if (expected.Length > 0) hitsRows++;
            if (f[1].Contains(',')) ranges++;
            if (f[1].Contains('<') || f[1].Contains('>')) anchored++;
            if (f[1].ToUpperInvariant().Any(c => c is 'B' or 'Z' or 'J')) ambiguity++;
            if (f[0].Any(char.IsLower)) lower++;
            if (f[1].Contains('x')) wildcard++;
        }

        Assert.That(hitsRows, Is.GreaterThan(200));
        Assert.That(rejected, Is.GreaterThan(30));
        Assert.That(ranges, Is.GreaterThan(300));
        Assert.That(anchored, Is.GreaterThan(150));
        Assert.That(ambiguity, Is.GreaterThan(40));
        Assert.That(lower, Is.GreaterThan(150));
        Assert.That(wildcard, Is.GreaterThan(60));
    }

    [Test]
    public void FindTranslatedPatternHits_EqualsFuzztran_AllRows()
    {
        var rows = ReadRows("fuzztran_oracle.tsv");
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(700));
        int hitRows = 0, reverse = 0, forward = 0, tables = 0, rejected = 0;
        foreach (var f in rows)
        {
            string label = $"{f[0]} / {f[1]} frame {f[2]} table {f[3]}";
            int table = int.Parse(f[3]);
            if (f[4] == "ERR")
            {
                Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindTranslatedPatternHits(f[0], f[1], f[2], table), label);
                rejected++;
                continue;
            }

            var hits = GenomicAnalyzer.FindTranslatedPatternHits(f[0], f[1], f[2], table);
            var expected = f[4] == "-" ? Array.Empty<string[]>() : f[4].Split(' ').Select(h => h.Split(':')).ToArray();
            string Proj(string frame, string ps, string pe, string tl) => $"{frame}:{ps}:{pe}:{tl}";
            Assert.That(hits.Select(h => Proj(h.Frame.ToString(), h.ProteinStart.ToString(), h.ProteinEnd.ToString(), h.Translation)).OrderBy(x => x, StringComparer.Ordinal).ToArray(),
                Is.EqualTo(expected.Select(e => Proj(e[0], e[1], e[2], e[3])).OrderBy(x => x, StringComparer.Ordinal).ToArray()), label);
            // the printed coordinates of the forward frames equal the program's GFF start / end exactly (the reverse frames print S - |frame| - 3 start + 4 / S - |frame| - 3 end + 2, A5-7)
            var printed = hits.Where(h => !h.IsReverseStrand).Select(h => $"{h.Frame}:{h.ProteinStart}:{h.EmbossStart}:{h.EmbossEnd}").OrderBy(x => x, StringComparer.Ordinal).ToArray();
            var gff = expected.Where(e => int.Parse(e[0]) > 0).Select(e => $"{e[0]}:{e[1]}:{e[4]}:{e[5]}").OrderBy(x => x, StringComparer.Ordinal).ToArray();
            Assert.That(printed, Is.EqualTo(gff), label + " coordinates");
            // the program passes the END of the sequence instead of its begin to fuzztran_SourceFeature for frame -3 of -frame R / 6 (fuzztran.c)
            int Anchor(TranslatedPatternHit h) => h.Frame == -3 && f[2].Trim().ToUpperInvariant() is "R" or "6" ? Seqeron.Genomics.Infrastructure.EmbossSequenceReader.Clean(f[0]).Length : 1;
            Assert.That(hits.Where(h => h.IsReverseStrand).All(h => h.EmbossStart == Anchor(h) - Math.Abs(h.Frame) - 3 * h.ProteinStart + 4 && h.EmbossEnd == Anchor(h) - Math.Abs(h.Frame) - 3 * h.ProteinEnd + 2), Is.True);
            // the forward-strand location is consistent with the protein hit
            foreach (var h in hits)
            {
                Assert.That(h.Length, Is.EqualTo(3 * h.Translation.Length), label);
                Assert.That(h.Position, Is.InRange(0, f[0].Length - h.Length), label);
                if (h.Frame > 0) Assert.That((h.Position + 1, h.Position + h.Length), Is.EqualTo((h.EmbossStart, h.EmbossEnd)), label);
            }

            if (expected.Length > 0) hitRows++;
            if (hits.Any(h => h.IsReverseStrand)) reverse++;
            if (hits.Any(h => !h.IsReverseStrand)) forward++;
            if (table != 0) tables++;
        }

        Assert.That(hitRows, Is.GreaterThan(300));
        Assert.That(reverse, Is.GreaterThan(100));
        Assert.That(forward, Is.GreaterThan(150));
        Assert.That(tables, Is.GreaterThan(200));
        Assert.That(rejected, Is.GreaterThanOrEqualTo(0));
    }

    [Test]
    public void ReverseFrameHit_LocatesTheCodonsOnTheForwardStrand()
    {
        // the reverse complement of GGATTA is TAATCC: frame -1 translates it (TAA = stop, TCC = S); the residue S covers forward bases 0..2 (GGA, whose reverse complement is TCC)
        var hits = GenomicAnalyzer.FindTranslatedPatternHits("GGATTA", "S", "-1");
        Assert.That(hits.Select(h => (h.Frame, h.Translation, h.Position, h.Length)).ToArray(), Is.EqualTo(new[] { (-1, "S", 0, 3) }));
        var forward = GenomicAnalyzer.FindTranslatedPatternHits("AAAGGATTA", "G", "1");
        Assert.That(forward.Select(h => (h.Frame, h.Translation, h.Position, h.Length)).ToArray(), Is.EqualTo(new[] { (1, "G", 3, 3) }));
    }

    [Test]
    public void ProteinPattern_Parser_Letters_Anchors_AndErrors()
    {
        // letters: B = DN, Z = EQ, J = IL, x = any (fuzzpro 6.6.0 runs); a single trailing '*' is not part of the sequence
        Assert.That(GenomicAnalyzer.FindProteinPatternHits("ABDNEQIL", "B").Select(h => h.Position), Is.EqualTo(new[] { 2, 3 }));
        Assert.That(GenomicAnalyzer.FindProteinPatternHits("ABDNEQIL", "Z").Select(h => h.Position), Is.EqualTo(new[] { 4, 5 }));
        Assert.That(GenomicAnalyzer.FindProteinPatternHits("ABDNEQIL", "J").Select(h => h.Position), Is.EqualTo(new[] { 6, 7 }));
        Assert.That(GenomicAnalyzer.FindProteinPatternHits("ACDEF*", "x(2)>").Select(h => (h.Position, h.Length)), Is.EqualTo(new[] { (3, 2) }));
        // a start-anchored range reports only the longest match; an end-anchored one every window
        Assert.That(GenomicAnalyzer.FindProteinPatternHits("AAAC", "<A(1,3)").Select(h => h.Length), Is.EqualTo(new[] { 3 }));
        Assert.That(GenomicAnalyzer.FindProteinPatternHits("CAAA", "A(1,3)>").Select(h => (h.Position, h.Length)), Is.EqualTo(new[] { (1, 3), (2, 2), (3, 1) }));
        Assert.That(ProteinPattern.Parse("[DE](2)HS{P}x(2)Px(2,4)C").ToString(), Does.StartWith("[DE](2)HS{"));
        foreach (string bad in new[] { "", "[BA]", "{Z}", "[X]", "[JA]", "*", "1", "A(0)", "A(2,1)", "A(,2)", "[]", "{}", "A>C", "A<C", "[A", "A)", "A(" })
            Assert.That(ProteinPattern.TryParse(bad, out _, out string? error), Is.False, bad + " -> " + error);
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindProteinPatternHits(null!, "A"));
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindTranslatedPatternHits("ACGT", "A", "5"));
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindTranslatedPatternHits("ACGT", "A", "1", 99));
    }
}
