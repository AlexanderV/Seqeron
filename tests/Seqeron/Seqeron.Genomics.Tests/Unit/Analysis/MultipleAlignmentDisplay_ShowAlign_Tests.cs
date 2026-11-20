// ALIGN-STATS-001 (B09 audit A3-13c, F148) — EMBOSS showalign: MultipleAlignmentDisplay.ShowAlign.
// Oracle: scripts/analysis/generate_showalign_oracle.py — 600 random alignments (nucleotide / protein, case, gaps, ragged rows) with random -refseq -bottom -show -order -similarcase
// -consensus -number -ruler -width -margin -uppercase -matrix -plurality -setcase -identity -gaps; the text written by showalign 6.6.0. Evidence: docs/Evidence/ALIGN-STATS-001-Evidence.md.

using System.Globalization;
using System.Text.Json;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class MultipleAlignmentDisplay_ShowAlign_Tests
{
    private sealed record Case(string[] Names, string[] Rows, Dictionary<string, string> Args, string Out);

    private static readonly Lazy<List<Case>> OracleCases = new(() => ReadOracle());

    private static List<Case> ReadOracle(string file = "showalign_oracle.jsonl")
    {
        string name = "Seqeron.Genomics.Tests.TestData.Analysis." + file;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var cases = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            using var doc = JsonDocument.Parse(line);
            var root = doc.RootElement;
            var args = root.GetProperty("args").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.ToString());
            cases.Add(new Case(
                root.GetProperty("names").EnumerateArray().Select(e => e.GetString()!).ToArray(),
                root.GetProperty("rows").EnumerateArray().Select(e => e.GetString()!).ToArray(),
                args, root.GetProperty("out").GetString()!));
        }
        return cases;
    }

    private static bool Flag(Dictionary<string, string> a, string key, bool dflt) => a.TryGetValue(key, out var v) ? v == "Y" : dflt;

    private static ShowAlignOptions OptionsOf(Dictionary<string, string> a)
    {
        var uppercase = new List<(int, int)>();
        if (a.TryGetValue("uppercase", out var u))
        {
            var p = u.Split('-');
            uppercase.Add((int.Parse(p[0], CultureInfo.InvariantCulture), int.Parse(p[1], CultureInfo.InvariantCulture)));
        }

        return new ShowAlignOptions
        {
            Matrix = a.TryGetValue("matrix", out var m) ? SubstitutionMatrix.Get(m) : null,
            RefSeq = a["refseq"],
            Bottom = Flag(a, "bottom", true),
            Show = (ShowAlignContent)a["show"][0],
            Order = (ShowAlignOrder)a["order"][0],
            SimilarCase = Flag(a, "similarcase", true),
            Consensus = Flag(a, "consensus", true),
            Number = Flag(a, "number", true),
            Ruler = Flag(a, "ruler", true),
            Width = int.Parse(a["width"], CultureInfo.InvariantCulture),
            Margin = int.Parse(a["margin"], CultureInfo.InvariantCulture),
            Uppercase = uppercase,
            Plurality = a.TryGetValue("plurality", out var pl) ? float.Parse(pl, CultureInfo.InvariantCulture) : 50f,
            SetCase = a.TryGetValue("setcase", out var sc) ? float.Parse(sc, CultureInfo.InvariantCulture) : null,
            Identity = a.TryGetValue("identity", out var id) ? float.Parse(id, CultureInfo.InvariantCulture) : 0f,
            Gaps = Flag(a, "gaps", true),
            Begin = a.TryGetValue("sbegin", out var sb) ? int.Parse(sb, CultureInfo.InvariantCulture) : 0,
            End = a.TryGetValue("send", out var se) ? int.Parse(se, CultureInfo.InvariantCulture) : 0,
        };
    }

    [Test]
    public void ShowAlign_EqualsShowalign_WithSbeginAndSend()
    {
        // A4-9: 400 rows of the installed 6.6.0 program with random -sbegin / -send (0, 1-based, negative = from the end, beyond the alignment, an end before the begin); the displayed columns are cut,
        // the consensus is that of the whole alignment and the ticks / numbers keep the positions of the whole alignment (showalign.c: begin = ajSeqsetGetBegin - 1, end = ajSeqsetGetEnd - 1)
        var cases = ReadOracle("showalign_range_oracle.jsonl");
        Assert.That(cases, Has.Count.EqualTo(400));
        var failures = new List<string>();
        int begins = 0, ends = 0, negative = 0, beyond = 0, cut = 0;
        foreach (var c in cases)
        {
            string actual = MultipleAlignmentDisplay.ShowAlign(c.Names, c.Rows, OptionsOf(c.Args));
            if (actual != c.Out)
                failures.Add($"{string.Join(" ", c.Args.Select(kv => $"-{kv.Key} {kv.Value}"))}\n{string.Join("\n", c.Names.Zip(c.Rows, (n, r) => $">{n} {r}"))}\n--- expected\n{c.Out}--- actual\n{actual}");
            int mlen = c.Rows.Max(r => r.Length);
            int sb = c.Args.TryGetValue("sbegin", out var b0) ? int.Parse(b0, CultureInfo.InvariantCulture) : 0, se = c.Args.TryGetValue("send", out var e0) ? int.Parse(e0, CultureInfo.InvariantCulture) : 0;
            if (sb != 0) begins++;
            if (se != 0) ends++;
            if (sb < 0 || se < 0) negative++;
            if (sb > mlen || se > mlen || sb < -mlen || se < -mlen) beyond++;
            if (sb != 0 || se != 0) cut++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count}: " + string.Join("\n=====\n", failures.Take(2)));
        Assert.That(begins, Is.GreaterThan(200));
        Assert.That(ends, Is.GreaterThan(200));
        Assert.That(negative, Is.GreaterThan(100));
        Assert.That(beyond, Is.GreaterThan(50));
        Assert.That(cut, Is.GreaterThan(300));
    }

    [Test]
    public void ShowAlign_SequenceRange_WorkedValues()
    {
        // measured: three 18-column rows, -sbegin 3 -send 12 shows the columns 3..12 with the numbers of the whole alignment; -sbegin 30 = the last column; -sbegin 8 -send 4 = column 8 only
        string[] names = ["s1", "s2", "s3"];
        string[] rows = ["ACGT--ACGTAAGGCCTT", "ACGTTTACGTA-GGCCTA", "ACGTTTACCTAAGG-CTT"];
        string text = MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Begin = 3, End = 12 });
        Assert.That(text, Is.EqualTo("                  10  \n          --:----|----:---\ns1        ..--......\ns2        .........-\ns3        ......C...\nConsensus GTTTACGTAA\n\n"));
        Assert.That(MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Begin = 8, End = 4 }), Does.Contain("Consensus C\n"));
        Assert.That(MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Begin = 30 }), Does.Contain("Consensus T\n"));
    }

    [Test]
    public void ShowAlign_EqualsShowalign_AllCases()
    {
        var cases = OracleCases.Value;
        Assert.That(cases, Has.Count.EqualTo(600));
        var failures = new List<string>();
        int refSeq = 0, bottom = 0, similar = 0, dissimilar = 0, identities = 0, all = 0, alpha = 0, bySimilarity = 0, narrow = 0, protein = 0, ranges = 0, ragged = 0, noMargin = 0;
        foreach (var c in cases)
        {
            string actual = MultipleAlignmentDisplay.ShowAlign(c.Names, c.Rows, OptionsOf(c.Args));
            if (actual != c.Out)
                failures.Add($"{string.Join(" ", c.Args.Select(kv => $"-{kv.Key} {kv.Value}"))}\n{string.Join("\n", c.Names.Zip(c.Rows, (n, r) => $">{n} {r}"))}\n--- expected\n{c.Out}--- actual\n{actual}");
            if (c.Args["refseq"] != "0") refSeq++;
            if (c.Args.TryGetValue("bottom", out var b) && b == "N") bottom++;
            if (c.Args["show"] == "S") similar++;
            if (c.Args["show"] == "D") dissimilar++;
            if (c.Args["show"] == "I") identities++;
            if (c.Args["show"] == "A") all++;
            if (c.Args["order"] == "A") alpha++;
            if (c.Args["order"] == "S") bySimilarity++;
            if (int.Parse(c.Args["width"]) < 20) narrow++;
            if (c.Rows[0].ToUpperInvariant().Any(ch => "EFHIKLMPQ".Contains(ch))) protein++;
            if (c.Args.ContainsKey("uppercase")) ranges++;
            if (c.Rows.Select(r => r.Length).Distinct().Count() > 1) ragged++;
            if (c.Args["margin"] == "0") noMargin++;
        }

        Assert.That(failures, Is.Empty, string.Join("\n=====\n", failures.Take(2)));
        Assert.That(refSeq, Is.GreaterThan(150));
        Assert.That(bottom, Is.GreaterThan(30));
        Assert.That(similar, Is.GreaterThan(50));
        Assert.That(dissimilar, Is.GreaterThan(30));
        Assert.That(identities, Is.GreaterThan(50));
        Assert.That(all, Is.GreaterThan(30));
        Assert.That(alpha, Is.GreaterThan(50));
        Assert.That(bySimilarity, Is.GreaterThan(30));
        Assert.That(narrow, Is.GreaterThan(100));
        Assert.That(protein, Is.GreaterThan(150));
        Assert.That(ranges, Is.GreaterThan(100));
        Assert.That(ragged, Is.GreaterThan(30));
        Assert.That(noMargin, Is.GreaterThan(30));
    }

    [Test]
    public void WorkedValues_RunOnTheBinary()
    {
        // showalign on s1 acgtac-gtaGT / s2 acgtaccgtaGA / s3 acctaccgtTGA (default options)
        string text = MultipleAlignmentDisplay.ShowAlign(["s1", "s2", "s3"], ["acgtac-gtaGT", "acgtaccgtaGA", "acctaccgtTGA"]);
        Assert.That(text, Is.EqualTo(
            "                  10  \n          ----:----|--\ns1        ......-....T\ns2        ............\ns3        ..C......T..\nConsensus acgtaccgtaGA\n\n"));
        // -show I -refseq 2: the reference first in full, then the others as identities, the consensus last
        text = MultipleAlignmentDisplay.ShowAlign(["s1", "s2", "s3"], ["acgtac-gtaGT", "acgtaccgtaGA", "acctaccgtTGA"],
            new ShowAlignOptions { Show = ShowAlignContent.Identities, RefSeq = "2" });
        Assert.That(text, Does.Contain("s2        acgtaccgtaGA\ns1        acgtac-gtaG.\ns3        ac.taccgt.GA\nConsensus acgtaccgtaGA\n"));
    }

    [Test]
    public void Validation()
    {
        Assert.Throws<ArgumentException>(() => MultipleAlignmentDisplay.ShowAlign([], []));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentDisplay.ShowAlign(["a"], ["AC", "AC"]));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentDisplay.ShowAlign(["a", "b"], ["AC", "AC"], new ShowAlignOptions { RefSeq = "3" }));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentDisplay.ShowAlign(["a", "b"], ["AC", "AC"], new ShowAlignOptions { RefSeq = "nope" }));
        Assert.Throws<ArgumentOutOfRangeException>(() => MultipleAlignmentDisplay.ShowAlign(["a", "b"], ["AC", "AC"], new ShowAlignOptions { Width = 0 }));
    }

    // ------------------------------------------------------------------ A4-2 (b) / A4-5: extremes of width and margin, allocation linear in the output

    [Test]
    public void Width_BeyondTheAlignment_ShowsOneBlock_TheSameTextAsTheAlignmentWidth_AndInstantly()
    {
        string[] rows = ["acgtac-gtaGT", "acgtaccgtaGA", "acctaccgtTGA"];
        string[] names = ["s1", "s2", "s3"];
        long before = GC.GetAllocatedBytesForCurrentThread();
        string huge = MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Width = int.MaxValue });
        // A5-4: the block loop is bounded by the alignment length, not by Width (it ran to int.MaxValue and took 9.4 s for a 9-column alignment): the allocation of the call is the text and a few lists
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(8L << 20));
        Assert.That(huge, Is.EqualTo(MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Width = 12 })));
        Assert.That(MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Width = int.MaxValue - 1, Margin = 3 }),
            Is.EqualTo(MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Width = 12, Margin = 3 })));
    }

    [Test]
    public void Margin_OrWidth_AboveTheOutputLimit_IsRefusedBeforeAnythingIsBuilt()
    {
        string[] rows = ["acgtac-gtaGT", "acgtaccgtaGA"];
        string[] names = ["s1", "s2"];
        long before = GC.GetAllocatedBytesForCurrentThread();
        var ex = Assert.Throws<ArgumentException>(() => MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Margin = int.MaxValue }))!;
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(64L << 20), "the limit is checked before the display is built");
        Assert.That(ex.Message, Does.StartWith("EMBOSS showalign builds the whole display text").And.Contain("MultipleAlignmentDisplay.ShowAlignMaxOutputCharacters").And.Contain("margin"));
        Assert.That(MultipleAlignmentDisplay.ShowAlignMaxOutputCharacters, Is.EqualTo(134_217_728L));
        Assert.Throws<ArgumentException>(() => MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Margin = 200_000_000 }));
        // width 1 over a long alignment is the same quantity from the other side: 2 x 40 million columns would be 40 million blocks
        string longRow = new('A', 40_000_000);
        Assert.Throws<ArgumentException>(() => MultipleAlignmentDisplay.ShowAlign(names, [longRow, longRow], new ShowAlignOptions { Width = 1 }));
        // the last margin that fits is accepted
        Assert.DoesNotThrow(() => MultipleAlignmentDisplay.ShowAlign(names, rows, new ShowAlignOptions { Margin = 1000 }));
    }

    [Test]
    public void Output_AllocationIsLinearInTheText_NotQuadratic()
    {
        // width 1 over two 20 000-column rows: 20 000 blocks; copying the whole row for every block allocated about 1.6 GB (60 GB for two 100 kb rows), the slice allocates the text
        string a = new string('A', 20_000), b = new string('A', 19_990) + "CCCCCCCCCC";
        long before = GC.GetAllocatedBytesForCurrentThread();
        string text = MultipleAlignmentDisplay.ShowAlign(["x", "y"], [a, b], new ShowAlignOptions { Width = 1 });
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(text.Length, Is.GreaterThan(20_000 * 5));
        Assert.That(allocated, Is.LessThan(120L << 20), $"allocated {allocated} bytes for a {text.Length}-character text");
        // the blocks still carry the right residue: the last block of row y is the C that differs
        Assert.That(text, Does.Contain("y         C\n"));
    }
}
