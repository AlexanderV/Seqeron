using System.IO.Compression;
using System.Text;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3e (B09 F241): the feature readers of EMBOSS 6.6.0 (<c>ajfeatread.c</c>: embl, genbank, refseq, genpept, refseqp, gff2, gff3, the gff3 of EMBOSS before 6.4, swiss, pir, dasgff) read every feature table of a file as the
/// installed <c>featcopy -features FILE -outfeat out -offormat debug</c> does (the debug trace prints type, source, score, strand, positions, flags, group, exon, tags and sub-features of every feature):
/// TestData/Alignment/featread_oracle.jsonl.gz, scripts/alignment/generate_featread_oracle.py; no format (all the formats the library tries for an unknown input), the <c>format::file</c> prefix (a format that does not read the
/// file, an unknown name), <c>-fbegin1</c> / <c>-fend1</c> (the trimmed first table), damaged and mixed files, the group field the GFF2 reader loops on forever.
/// </summary>
[TestFixture]
public class EmbossFeatureReaders_Oracle_Tests
{
    // ajFeattableWrite of the debug format
    private static string WriteDebug(EmbossSeqFeatureTable table)
    {
        if (table.Rev || table.GetBegin() > 1 || table.GetEnd() < table.Length)
        {
            var copy = table.NewFtable();
            copy.Trim();
            return copy.Debug();
        }

        return table.Debug();
    }

    internal static (int Rc, string Output, string Error) RunDebug(string text, string ufo, IReadOnlyList<string> extra)
    {
        int fbegin = 1;
        int? fend = null;
        for (int i = 0; i + 1 < extra.Count; i += 2)
        {
            if (extra[i] == "-fbegin1")
                fbegin = int.Parse(extra[i + 1], System.Globalization.CultureInfo.InvariantCulture);
            else if (extra[i] == "-fend1")
                fend = int.Parse(extra[i + 1], System.Globalization.CultureInfo.InvariantCulture);
        }

        var context = new EmbossFeatureReadContext();
        EmbossFeatureInput input;
        try
        {
            input = new EmbossFeatureInput(ufo, name => name == "in.dat" ? text : null, fbegin, fend, context);
        }
        catch (EmbossFeatureDiedException e)
        {
            return (1, "", string.Join("\n", [.. e.Errors, e.Message]) + "\n");
        }

        var sb = new StringBuilder();
        EmbossSeqFeatureTable? table;
        while ((table = input.Next()) is not null)
            sb.Append(WriteDebug(table));
        return (0, sb.ToString(), string.Join("\n", input.Errors) + (input.Errors.Count > 0 ? "\n" : ""));
    }

    [Test]
    public void FeatureReaders_EqualFeatcopyDebug_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.featread_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var perKind = new Dictionary<string, int>();
        int rows = 0, hang = 0, died = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string kind = j.GetProperty("kind").GetString()!;
            string text = j.GetProperty("text").GetString()!;
            string ufo = j.GetProperty("ufo").GetString()!;
            var extra = j.GetProperty("extra").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            string expected = j.GetProperty("debug").GetString()!;
            string err = j.GetProperty("err").GetString()!;
            string id = $"row {j.GetProperty("id")} ({kind} {ufo} {string.Join(' ', extra)})";
            if (j.GetProperty("hang").GetBoolean())
            {
                hang++;
                try
                {
                    RunDebug(text, ufo, extra);
                    failures.Add($"{id}: the program loops forever, the port returned");
                }
                catch (InvalidOperationException)
                {
                }

                continue;
            }

            try
            {
                var (grc, output, gerr) = RunDebug(text, ufo, extra);
                string[] eLines = [.. err.Split('\n').Where(l => l.Length > 0)];
                string[] gLines = [.. gerr.Split('\n').Where(l => l.Length > 0)];
                if (grc != rc)
                    failures.Add($"{id}: exit code {grc} expected {rc}");
                else if (output != expected)
                    failures.Add($"{id}: output differs at {FirstDifference(output, expected)}");
                else if (!eLines.SequenceEqual(gLines))
                    failures.Add($"{id}: errors '{string.Join('|', gLines)}' expected '{string.Join('|', eLines)}'");
                else if (rc == 0 && expected.Length > 0)
                    perKind[kind] = perKind.GetValueOrDefault(kind) + 1;
                if (rc != 0)
                    died++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(8)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(1000));
    }

    private static string FirstDifference(string got, string expected)
    {
        int i = 0;
        while (i < got.Length && i < expected.Length && got[i] == expected[i])
            i++;
        string Cut(string s) => s[Math.Max(0, i - 40)..Math.Min(s.Length, i + 60)].Replace("\n", "\\n", StringComparison.Ordinal);
        return $"{i}: got '{Cut(got)}' expected '{Cut(expected)}'";
    }
}
