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

    private static readonly Lazy<List<Case>> OracleCases = new(ReadOracle);

    private static List<Case> ReadOracle()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Analysis.showalign_oracle.jsonl";
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
        };
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
}
