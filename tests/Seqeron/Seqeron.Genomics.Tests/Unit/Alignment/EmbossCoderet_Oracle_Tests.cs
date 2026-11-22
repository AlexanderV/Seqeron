#pragma warning disable S127
using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F259, WP58): <c>coderet</c> of EMBOSS 6.6.0 (<c>emboss/coderet.c</c>) through the installed program:
/// <c>coderet -auto -seqall in.dat [-sformat1 F] -outfile o.cod [-cdsoutseq [F::]f | ''] [-mrnaoutseq ..] [-translationoutseq ..] [-restoutseq ..] [-osformat3 F .. -osformat6 F]</c>.
/// TestData/Alignment/coderet_oracle.jsonl.gz, scripts/alignment/generate_coderet_oracle.py: EMBL entries with CDS / mRNA / exon features (joins, complements, orders, partial ends, <c>/translation</c> tags and notes, several entries per file)
/// and the entries with feature tables of the reader oracle (embl, genbank, gff3, gff2, nbrf, swiss, refseqp, FASTA, damaged files); the output files named by the first sequence or by the user, with a format prefix or <c>-osformatN</c>,
/// one to four of them. The log file, every sequence file the program wrote (names and texts) and the messages of the library are compared byte for byte; the runs the reference cannot run (a fault, a loop, memory errors) are refused by the port.
/// </summary>
[TestFixture]
public class EmbossCoderet_Oracle_Tests
{
    private static string FirstDifference(string a, string b)
    {
        int n = Math.Min(a.Length, b.Length);
        int i = 0;
        while (i < n && a[i] == b[i])
            i++;
        return $"offset {i} (got '{(i < a.Length ? a[i..Math.Min(a.Length, i + 30)] : "<end>").Replace("\n", "\\n")}', expected '{(i < b.Length ? b[i..Math.Min(b.Length, i + 30)] : "<end>").Replace("\n", "\\n")}')";
    }

    private static EmbossCoderetRequest Request(string text, List<string> a, DateTime date)
    {
        var q = new List<string>();
        string? cds = null, mrna = null, trn = null, rest = null, f3 = null, f4 = null, f5 = null, f6 = null;
        for (int i = 0; i < a.Count; i++)
        {
            switch (a[i])
            {
                case "-cdsoutseq": cds = a[++i]; break;
                case "-mrnaoutseq": mrna = a[++i]; break;
                case "-translationoutseq": trn = a[++i]; break;
                case "-restoutseq": rest = a[++i]; break;
                case "-osformat3": f3 = a[++i]; break;
                case "-osformat4": f4 = a[++i]; break;
                case "-osformat5": f5 = a[++i]; break;
                case "-osformat6": f6 = a[++i]; break;
                default:
                    q.Add(a[i]);
                    if (a[i] is "-sformat1" or "-sbegin1")
                        q.Add(a[++i]);
                    break;
            }
        }

        var quals = EmbossSequenceQualifiers.ParseFile(string.Join(' ', q));
        return new EmbossCoderetRequest("in.dat", new Dictionary<string, string> { ["in.dat"] = text }, quals, cds, mrna, trn, rest, f3, f4, f5, f6) { RunDate = date };
    }

    [Test]
    public void Coderet_EqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.coderet_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var formats = new HashSet<string>();
        int rows = 0, compared = 0, refused = 0, withFiles = 0, sequencesWritten = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string fmt = j.GetProperty("fmt").GetString()!;
            string text = j.GetProperty("text").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            bool memerr = j.GetProperty("memerr").GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25 or -6 or -24;
            string expectedLog = j.GetProperty("log").GetString()!, err = j.GetProperty("err").GetString()!;
            var expectedFiles = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
            string id = $"row {j.GetProperty("id")} [{fmt} {j.GetProperty("how").GetString()} {string.Join(' ', args)}]";
            try
            {
                var o = EmbossCoderetProgram.Coderet(Request(text, args, date));
                if (memerr)
                    continue;
                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                string messages = o.Messages.Count > 0 ? string.Join("\n", o.Messages) + "\n" : "";
                string? why = null;
                if ((o.Died ? 1 : 0) != rc)
                    why = $"exit code {(o.Died ? 1 : 0)} expected {rc}";
                else if (!o.Died && o.Report != expectedLog)
                    why = $"log differs at {FirstDifference(o.Report, expectedLog)}";
                else if (!o.Died && !o.Files.Keys.Order(StringComparer.Ordinal).SequenceEqual(expectedFiles.Keys.Order(StringComparer.Ordinal)))
                    why = $"files {string.Join(",", o.Files.Keys)} expected {string.Join(",", expectedFiles.Keys)}";
                else if (!o.Died && o.Files.Where(kv => kv.Value != expectedFiles[kv.Key]).Select(kv => $"file {kv.Key} differs at {FirstDifference(kv.Value, expectedFiles[kv.Key])}").FirstOrDefault() is { } fileDiff)
                    why = fileDiff;
                else if (messages != err)
                    why = $"messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'";
                if (why is not null)
                {
                    failures.Add($"{id}: {why}");
                    continue;
                }

                compared++;
                if (expectedFiles.Values.Any(v => v.Length > 0))
                {
                    withFiles++;
                    formats.Add(fmt);
                    sequencesWritten += expectedFiles.Values.Sum(v => v.Count(c => c == '\n'));
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
        Assert.That(rows, Is.GreaterThanOrEqualTo(1000));
        Assert.That(compared, Is.GreaterThanOrEqualTo(900), $"compared {compared}, refused {refused}");
        Assert.That(withFiles, Is.GreaterThanOrEqualTo(430));
        Assert.That(formats.Count, Is.GreaterThanOrEqualTo(3));
        Assert.That(sequencesWritten, Is.GreaterThan(0));
    }
}
