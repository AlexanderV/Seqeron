using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F257, WP58): <c>entret</c> of EMBOSS 6.6.0 (<c>emboss/entret.c</c>): the text of every entry the sequence readers of <c>ajseqread.c</c> read (<c>ajSeqGetEntryS</c>, stored by <c>ajTextinStoreReadline</c>) through the installed program:
/// <c>entret -auto -sequence SPEC [-sformat1 F] [-sbegin1 N] [-firstonly] -outfile o.txt</c>. TestData/Alignment/entret_oracle.jsonl.gz, scripts/alignment/generate_entret_oracle.py: the inputs of the reader oracle in 62 formats
/// (writer outputs of seqret, entries with feature tables, hand-made texts, damaged files), with a format named, the <c>format::</c> prefix, no format, queries on the file, ranges. The output file and the messages of the library are
/// compared byte for byte; the runs the reference cannot run (a fault, a loop, memory errors) are refused by the port.
/// </summary>
[TestFixture]
public class EmbossEntret_Oracle_Tests
{
    private static string FirstDifference(string a, string b)
    {
        int n = Math.Min(a.Length, b.Length);
        int i = 0;
        while (i < n && a[i] == b[i])
            i++;
        return $"offset {i} (got '{(i < a.Length ? a[i..Math.Min(a.Length, i + 30)] : "<end>").Replace("\n", "\\n")}', expected '{(i < b.Length ? b[i..Math.Min(b.Length, i + 30)] : "<end>").Replace("\n", "\\n")}')";
    }

    [Test]
    public void Entret_EqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.entret_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var formats = new HashSet<string>();
        int rows = 0, compared = 0, refused = 0, withOutput = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string fmt = j.GetProperty("fmt").GetString()!;
            var files = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            string spec = j.GetProperty("spec").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            bool memerr = j.TryGetProperty("memerr", out var me) && me.GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25 or -6 or -24;
            string expected = j.GetProperty("out").GetString()!, err = j.GetProperty("err").GetString()!;
            string id = $"row {j.GetProperty("id")} [{fmt} {j.GetProperty("how").GetString()} {spec} {string.Join(' ', args)}]";
            try
            {
                var q = EmbossSequenceQualifiers.ParseFile(string.Join(' ', args.Where(a => a != "-firstonly")));
                var o = EmbossEntretProgram.Entret(new EmbossEntretRequest(spec, files, q, args.Contains("-firstonly")));
                if (memerr)
                    continue;
                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                string messages = o.Messages.Count > 0 ? string.Join("\n", o.Messages) + "\n" : "";
                if ((o.Died ? 1 : 0) != rc)
                    failures.Add($"{id}: exit code {(o.Died ? 1 : 0)} expected {rc}");
                else if ((o.Died ? "" : o.Text) != expected)
                    failures.Add($"{id}: text differs at {FirstDifference(o.Text, expected)}");
                else if (messages != err)
                    failures.Add($"{id}: messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else
                {
                    compared++;
                    if (expected.Length > 0)
                    {
                        withOutput++;
                        formats.Add(fmt);
                    }
                }
            }
            catch (Exception e) when (e is InvalidOperationException or NotSupportedException)
            {
                if (memerr || refusal)
                    refused++;
                else
                    failures.Add($"{id}: refused ({e.Message[..Math.Min(70, e.Message.Length)]}) but the reference ran (exit code {rc})");
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(700));
        Assert.That(compared, Is.GreaterThanOrEqualTo(600), $"compared {compared}, refused {refused}");
        Assert.That(withOutput, Is.GreaterThanOrEqualTo(500));
        Assert.That(formats.Count, Is.GreaterThanOrEqualTo(50));
    }
}
