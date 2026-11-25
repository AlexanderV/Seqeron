// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.FeatureTools;
using Seqeron.Genomics.Emboss.SequenceFormats;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-27d, B09 F233): the EMBOSS 6.6.0 alignment formats of three or more rows (<c>unknown</c>, <c>multiple</c>, <c>simple</c>, <c>srs</c>: <c>alignConsStats</c> consensus markup) written byte for byte as the installed
/// <c>aligncopy</c> does, and the refusal of the pair-markup formats (scripts/alignment/generate_emboss_multirow_oracle.py, emboss_multirow_oracle.tsv).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class EmbossMultiRowAlignment_Oracle_Tests
{
    private sealed record Row(string Id, string Format, int Width, string Fasta, string Output);

    private static string Unesc(string s) => System.Text.RegularExpressions.Regex.Replace(s, @"\\(.)", m => m.Groups[1].Value switch { "t" => "\t", "n" => "\n", "r" => "\r", var c => c });

    private static readonly Lazy<List<Row>> Rows = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_multirow_oracle.tsv";
        using var stream = TestDataResource.Open(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], int.Parse(f[2]), Unesc(f[3]), Unesc(f[4])));
        }

        return list;
    });

    private static List<EmbossAlignedSequence> Parse(string fasta)
    {
        var rows = new List<EmbossAlignedSequence>();
        foreach (var block in fasta.Split('>', StringSplitOptions.RemoveEmptyEntries))
        {
            var lines = block.Split('\n', StringSplitOptions.RemoveEmptyEntries);
            rows.Add(new EmbossAlignedSequence(lines[0].Split(' ')[0].Split('|')[^1], string.Concat(lines.Skip(1))));
        }

        return rows;
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

        return "same";
    }

    [Test]
    public void Oracle_Has200RowsPerMultiRowFormat()
    {
        Assert.Multiple(() =>
        {
            foreach (var f in new[] { "unknown", "multiple", "simple", "srs", "tcoffee" })
                Assert.That(Rows.Value.Count(r => r.Format == f && !r.Output.StartsWith("ERROR:", StringComparison.Ordinal)), Is.GreaterThanOrEqualTo(200), f);
        });
    }

    [Test]
    public void Write_EveryRow_EqualsAligncopy()
    {
        var failures = new List<string>();
        foreach (var r in Rows.Value)
        {
            Assert.That(EmbossAlignmentWriter.TryParseFormat(r.Format, out var format), Is.True, r.Format);
            var cmd = new List<string> { "-sequences a.fa", "-aformat " + r.Format };
            cmd.Add("-outfile stdout");
            cmd.Add("-auto");
            if (r.Width != 0)
                cmd.Add("-awidth " + r.Width);
            var options = new EmbossAlignOptions { Program = "aligncopy", CommandLine = cmd, RunDate = new DateTime(2026, 11, 21, 19, 0, 0, DateTimeKind.Unspecified), ReportFile = "stdout", Width = r.Width };
            var rows = Parse(r.Fasta);
            try
            {
                var w = new EmbossAlignmentWriter(format, options);
                w.Write(new EmbossMultiAlignment(rows));
                string actual = w.Close();
                if (r.Output.StartsWith("ERROR:", StringComparison.Ordinal))
                    failures.Add($"{r.Id} {r.Format} wrote but the program refused: {r.Output}");
                else if (actual != r.Output)
                    failures.Add($"{r.Id} {r.Format} w{r.Width}\n{FirstDifference(r.Output, actual)}");
            }
            catch (ArgumentException ex)
            {
                if (!r.Output.StartsWith("ERROR:Error: " + ex.Message.Split(" (Parameter")[0] + "\n", StringComparison.Ordinal))
                    failures.Add($"{r.Id} {r.Format} threw '{ex.Message}' but the program: {r.Output[..Math.Min(r.Output.Length, 200)]}");
            }
        }

        Assert.That(failures, Is.Empty, string.Join("\n=====\n", failures.Take(3)) + $"\n({failures.Count}: {string.Join(", ", failures.Take(40).Select(f => f.Split(' ')[0]))})");
    }
}
