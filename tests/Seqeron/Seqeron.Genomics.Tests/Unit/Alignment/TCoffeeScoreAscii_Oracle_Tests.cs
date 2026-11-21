using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-27d, B09 F237): the T-Coffee 13.41 <c>score_ascii</c> evaluation (<c>t_coffee -infile aln -mode evaluate -output score_ascii</c>, evaluate mode <c>triplet</c>) written byte for byte
/// by <see cref="SequenceAligner.TCoffeeEvaluate"/> + <see cref="TCoffeeEvaluation.ToScoreAscii"/> (scripts/alignment/generate_tcoffee_score_oracle.py, tcoffee_score_oracle.tsv).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class TCoffeeScoreAscii_Oracle_Tests
{
    private sealed record Row(string Id, string Kind, string Fasta, string Output);

    private static string Unesc(string s) => Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "t" => "\t", "n" => "\n", "r" => "\r", var c => c });

    private static readonly Lazy<List<Row>> Rows = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.tcoffee_score_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[2], Unesc(f[3]), Unesc(f[4])));
        }

        return list;
    });

    private static (List<string> Names, List<string> Rows) Parse(string fasta)
    {
        var names = new List<string>();
        var rows = new List<string>();
        foreach (var block in fasta.Split('>', StringSplitOptions.RemoveEmptyEntries))
        {
            var lines = block.Split('\n', StringSplitOptions.RemoveEmptyEntries);
            names.Add(lines[0].Split(' ')[0]);
            rows.Add(string.Concat(lines.Skip(1)));
        }

        return (names, rows);
    }

    private static string FirstDifference(string expected, string actual)
    {
        var e = expected.Split('\n');
        var a = actual.Split('\n');
        for (int i = 0; i < Math.Max(e.Length, a.Length); i++)
        {
            string x = i < e.Length ? e[i] : "<end>", y = i < a.Length ? a[i] : "<end>";
            if (x != y)
                return $"line {i}: expected '{x}' actual '{y}'";
        }

        return "";
    }

    [Test]
    public void Oracle_Has200Evaluations()
    {
        Assert.That(Rows.Value.Count(r => !r.Output.StartsWith("ERROR:", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(200));
    }

    [Test]
    public void ScoreAscii_EveryRow_IsByteEqualToTCoffee()
    {
        var failures = new List<string>();
        foreach (var r in Rows.Value)
        {
            var (names, rows) = Parse(r.Fasta);
            if (r.Output.StartsWith("ERROR:", StringComparison.Ordinal))
            {
                Assert.That(() => SequenceAligner.TCoffeeEvaluate(names, rows), NUnit.Framework.Throws.InstanceOf<Exception>(), r.Id);
                continue;
            }

            int cpu = int.Parse(Regex.Match(r.Output, @"CPU TIME:(\d+) sec").Groups[1].Value);
            string actual = SequenceAligner.TCoffeeEvaluate(names, rows).ToScoreAscii(50, cpu);
            if (actual != r.Output)
                failures.Add($"{r.Id} ({r.Kind}, {names.Count} rows): {FirstDifference(r.Output, actual)}");
        }

        Assert.That(failures, Is.Empty, string.Join("\n", failures.Take(8)) + $"\n({failures.Count})");
    }
}
