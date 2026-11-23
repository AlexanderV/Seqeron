using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-26, B09 F128): <see cref="PoaProgram.Run"/> (the POA 2.0 command line on top of <see cref="SequenceAligner.PoaRun"/>) against the installed Debian <c>poa</c>
/// 2.0+20060928. Every committed run (scripts/alignment/generate_poa_oracle.py, poa_oracle.tsv) must give the program's exit code, standard output and the PO, PIR and CLUSTAL files byte for byte.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_PoaOracle_Tests
{
    internal sealed record Row(string Id, string Kind, string[] Args, Dictionary<string, string> Files, int Rc, string Stdout, Dictionary<string, string> OutFiles);

    internal static readonly Lazy<List<Row>> Oracle = new(Load);

    private static Dictionary<string, string> Files(string field)
    {
        var d = new Dictionary<string, string>();
        string text = SequenceAligner_ProbConsOracle_Tests.Unescape(field);
        if (text.Length == 0) return d;
        foreach (var part in text.Split('\x1f'))
        {
            int eq = part.IndexOf('=');
            d[part[..eq]] = part[(eq + 1)..];
        }

        return d;
    }

    private static List<Row> Load()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.poa_oracle.tsv";
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#' || line.StartsWith("id\t", StringComparison.Ordinal))
                continue;
            var f = line.Split('\t');
            string args = SequenceAligner_ProbConsOracle_Tests.Unescape(f[2]);
            list.Add(new Row(f[0], f[1], args == "-" ? [] : args.Split('\x1f'), Files(f[3]), int.Parse(f[4]), SequenceAligner_ProbConsOracle_Tests.Unescape(f[5]), Files(f[6])));
        }

        return list;
    }

    public static IEnumerable<TestCaseData> Runs() => Oracle.Value.Select(r => new TestCaseData(r.Id).SetArgDisplayNames(r.Id));

    [Test]
    public void Oracle_Has350RunsOfEveryKind()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(350));
            Assert.That(rows.Select(r => r.Kind).Distinct(), Is.SupersetOf(new[] { "default", "global", "progressive", "fuse", "hb", "order", "nuc", "trunc", "case", "unknown", "format", "pairscores", "verbose", "matrix", "long", "error" }));
            Assert.That(rows.Count(r => r.Rc != 0), Is.GreaterThanOrEqualTo(5));
        });
    }

    [TestCaseSource(nameof(Runs))]
    public void Run_EqualsPoa(string id)
    {
        var r = Oracle.Value.First(x => x.Id == id);
        var res = PoaProgram.Run(r.Args, r.Files);
        Assert.Multiple(() =>
        {
            Assert.That(res.ExitCode, Is.EqualTo(r.Rc), id + " exit code: " + string.Join(' ', r.Args));
            Assert.That(res.Stdout, Is.EqualTo(r.Stdout), id + " stdout: " + string.Join(' ', r.Args));
            Assert.That(res.Files.Keys.OrderBy(k => k, StringComparer.Ordinal), Is.EqualTo(r.OutFiles.Keys.OrderBy(k => k, StringComparer.Ordinal)), id + " file names");
            foreach (var kv in r.OutFiles)
                Assert.That(res.Files[kv.Key], Is.EqualTo(kv.Value), id + " file " + kv.Key + ": " + string.Join(' ', r.Args));
        });
    }
}
