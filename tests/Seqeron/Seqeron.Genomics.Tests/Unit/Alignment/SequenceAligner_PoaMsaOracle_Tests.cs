using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-26b, B09 F223): MSA / PO file input of <see cref="PoaProgram.Run"/> (<c>-read_msa</c>, <c>-read_msa2</c>, <c>-read_msa_list</c>, <c>-subset</c>, <c>-subset2</c>, <c>-remove</c>, <c>-remove2</c>) against
/// the installed Debian <c>poa</c> 2.0+20060928 (scripts/alignment/generate_poa_msa_oracle.py, poa_msa_oracle.tsv): exit code, standard output, the <c>...Read / Wrote</c> lines of standard error and the PO / PIR / CLUSTAL files
/// byte for byte. Rows on which the program loops forever (status -999) or crashes (status above 128 / negative) must stop with a non-zero status here.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_PoaMsaOracle_Tests
{
    internal sealed record Row(string Id, string Kind, string[] Args, Dictionary<string, string> Files, int Rc, string Stdout, string Stderr, Dictionary<string, string> OutFiles);

    private static readonly Lazy<List<Row>> Oracle = new(Load);

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
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.poa_msa_oracle.tsv";
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
            list.Add(new Row(f[0], f[1], SequenceAligner_ProbConsOracle_Tests.Unescape(f[2]).Split('\x1f'), Files(f[3]), int.Parse(f[4]), SequenceAligner_ProbConsOracle_Tests.Unescape(f[5]),
                SequenceAligner_ProbConsOracle_Tests.Unescape(f[6]), Files(f[7])));
        }

        return list;
    }

    [Test]
    public void Oracle_HasTheRecordedRuns()
    {
        var rows = Oracle.Value;
        Assert.Multiple(() =>
        {
            Assert.That(rows, Has.Count.GreaterThanOrEqualTo(360));
            Assert.That(rows.Select(r => r.Kind).Distinct(), Is.SupersetOf(new[] { "single", "two", "msa+fasta", "subset", "remove", "list", "subset2", "remove2", "err", "hand-pir", "hand-clustal", "hand-po-two" }));
            Assert.That(rows.Count(r => r.Rc != 0), Is.GreaterThanOrEqualTo(20));
        });
    }

    private static string Diff(string actual, string expected)
    {
        var a = actual.Split('\n');
        var e = expected.Split('\n');
        for (int i = 0; i < Math.Max(a.Length, e.Length); i++)
        {
            string x = i < a.Length ? a[i] : "<end>", y = i < e.Length ? e[i] : "<end>";
            if (x != y) return $"line {i + 1}: actual '{x}' expected '{y}' (lines {a.Length} / {e.Length})";
        }

        return "equal";
    }

    [Test]
    public void Run_EveryCommandLine_EqualsPoa()
    {
        var failures = new List<string>();
        foreach (var r in Oracle.Value)
        {
            PoaProgramResult res;
            try
            {
                res = PoaProgram.Run(r.Args, r.Files);
            }
            catch (Exception e)
            {
                failures.Add($"{r.Id} {r.Kind} threw {e.GetType().Name}: {e.Message}");
                continue;
            }

            bool bad = r.Rc != 0;
            if (bad != (res.ExitCode != 0))
            {
                failures.Add($"{r.Id} {r.Kind} exit {res.ExitCode} (program {r.Rc}): {string.Join(' ', r.Args)} | {res.Stderr}");
                continue;
            }

            if (bad && (r.Rc > 128 || r.Rc < 0)) continue; // a crash / endless loop of the program: refused here
            string stderr = string.Concat(res.Stderr.Split('\n').Where(l => l.StartsWith("...", StringComparison.Ordinal)).Select(l => l + "\n"));
            if (res.Stdout != r.Stdout) failures.Add($"{r.Id} {r.Kind} stdout: {string.Join(' ', r.Args)}\n{Diff(res.Stdout, r.Stdout)}");
            else if (stderr != r.Stderr) failures.Add($"{r.Id} {r.Kind} stderr: {string.Join(' ', r.Args)}\n{Diff(stderr, r.Stderr)}");
            else if (!res.Files.Keys.OrderBy(k => k, StringComparer.Ordinal).SequenceEqual(r.OutFiles.Keys.OrderBy(k => k, StringComparer.Ordinal)))
                failures.Add($"{r.Id} {r.Kind} file names: {string.Join(' ', r.Args)}");
            else
            {
                foreach (var kv in r.OutFiles)
                {
                    if (res.Files[kv.Key] != kv.Value)
                    {
                        failures.Add($"{r.Id} {r.Kind} file {kv.Key}: {string.Join(' ', r.Args)}\n{Diff(res.Files[kv.Key], kv.Value)}");
                        break;
                    }
                }
            }
        }

        Assert.That(failures, Is.Empty, string.Join("\n=====\n", failures.Take(4)) + $"\n({failures.Count} of {Oracle.Value.Count}: {string.Join(' ', failures.Select(f => f.Split(' ')[0] + ':' + f.Split(' ')[1]))})");
    }
}
