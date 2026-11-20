using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-23, B09 F127): <see cref="ProbConsProgram.Run"/> (the ProbCons 1.12 command line on top of <see cref="SequenceAligner.ProbConsRun"/>) against the installed
/// Debian <c>probcons</c> 1.12-14. Every committed run (scripts/alignment/generate_probcons_oracle.py, probcons_oracle.tsv) must give the program's exit code, standard output,
/// standard error (heading, parameter set, tree, -v progress lines, error text) and written files (pairwise alignments, -annot values, trained parameter file) byte for byte.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ProbConsOracle_Tests
{
    internal sealed record Row(string Id, string Kind, string[] Args, Dictionary<string, string> Files, int Rc, string Stdout, string Stderr, Dictionary<string, string> OutFiles);

    internal static readonly Lazy<List<Row>> Oracle = new(Load);

    internal static string Unescape(string s)
    {
        var sb = new System.Text.StringBuilder();
        int i = 0;
        while (i < s.Length)
        {
            if (s[i] != '\\' || i + 1 >= s.Length)
            {
                sb.Append(s[i++]);
                continue;
            }

            char n = s[i + 1];
            sb.Append(n switch { 'n' => '\n', 't' => '\t', 'r' => '\r', _ => n });
            i += 2;
        }

        return sb.ToString();
    }

    private static Dictionary<string, string> Files(string field)
    {
        var d = new Dictionary<string, string>();
        string text = Unescape(field);
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
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.probcons_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#' || line.StartsWith("id\t", StringComparison.Ordinal))
                continue;
            var f = line.Split('\t');
            string args = Unescape(f[2]);
            list.Add(new Row(f[0], f[1], args == "-" ? [] : args.Split('\x1f'), Files(f[3]), int.Parse(f[4]), Unescape(f[5]), Unescape(f[6]), Files(f[7])));
        }

        return list;
    }

    public static IEnumerable<TestCaseData> Runs() => Oracle.Value.Select(r => new TestCaseData(r.Id).SetArgDisplayNames(r.Id));

    [Test]
    public void Oracle_Has400RunsOfEveryKind()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(400));
            Assert.That(rows.Select(r => r.Kind).Distinct(), Is.SupersetOf(new[] { "default", "c", "ir", "clustalw", "order", "annot", "pairs", "viterbi", "penalties", "cutoff", "pre", "train", "verbose", "format", "files", "long", "tiny", "error", "paramfile" }));
            Assert.That(rows.Count(r => r.Rc != 0), Is.GreaterThanOrEqualTo(15));
        });
    }

    [TestCaseSource(nameof(Runs))]
    public void Run_EqualsProbCons(string id)
    {
        var r = Oracle.Value.First(x => x.Id == id);
        var res = ProbConsProgram.Run(r.Args, r.Files);
        Assert.Multiple(() =>
        {
            Assert.That(res.ExitCode, Is.EqualTo(r.Rc), id + " exit code: " + string.Join(' ', r.Args));
            Assert.That(res.Stdout, Is.EqualTo(r.Stdout), id + " stdout: " + string.Join(' ', r.Args));
            Assert.That(res.Stderr, Is.EqualTo(r.Stderr), id + " stderr: " + string.Join(' ', r.Args));
            Assert.That(res.Files.Keys.OrderBy(k => k, StringComparer.Ordinal), Is.EqualTo(r.OutFiles.Keys.OrderBy(k => k, StringComparer.Ordinal)), id + " file names");
            foreach (var kv in r.OutFiles)
                Assert.That(res.Files[kv.Key], Is.EqualTo(kv.Value), id + " file " + kv.Key);
        });
    }
}
