#pragma warning disable S107

using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F255, WP58): the USA forms of the sequence input of EMBOSS 6.6.0 (<c>seqinUsaProcess</c> / <c>seqinListProcess</c> of ajseqread.c, <c>ajQuerystrParse*</c> of ajquery.c, <c>acdSetSeqall</c> of ajacd.c) through the installed
/// <c>seqret -auto [-feature] -sequence SPEC [-sformat1 F] [-sbegin1 N] [-send1 N] [-sreverse1] [-sid1 ID] [-sdbname1 DB] [-squick1] [-slower1] [-supper1] [-scircular1] [-iquery1 Q] [-ioffset1 N] [-firstonly]</c>.
/// TestData/Alignment/usa_oracle.jsonl.gz, scripts/alignment/generate_usa_oracle.py: 1 500 runs. SPEC is a file name with a <c>[begin:end:r]</c> range or the <c>:begin:end</c> / <c>:begin..end:r</c> forms, a <c>format::</c> prefix, a list file
/// (<c>@list</c>, <c>list::list</c>, <c>LIST:list</c>; comments, bad records, missing files, nested lists, ranges and formats of the entries) or a query on a file (<c>file:id</c>, <c>file:acc:value</c>, <c>file:des:value</c>, wildcards, <c>file{a|b}</c>, <c>%offset</c>);
/// the files are made from the inputs of the reader oracle (writer outputs of seqret, entries with feature tables), so the range and the reversal reach the sequences and their features. The FASTA file, the GFF3 file and the messages of the
/// library are compared byte for byte; the runs on which the reference program never stops (a list that names itself, a reader that never consumes an entry), dies by a fault (a treecon file with no sequence) or reads memory it does not own are refused
/// by the port with <see cref="InvalidOperationException"/>.
/// </summary>
[TestFixture]
public class EmbossUsa_Oracle_Tests
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
    public void UsaForms_EqualSeqret_AllRows() => CompareRows("Seqeron.Genomics.Tests.TestData.Alignment.usa_oracle.jsonl.gz", 1000, 900, 500, 200, 200, 100, 100, 0);

    private static void CompareRows(string name, int minRows, int minCompared, int minOutput, int minLists, int minRanges, int minQueries, int minReversed, int unused)
    {
        _ = unused;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, compared = 0, refused = 0, withOutput = 0, lists = 0, ranges = 0, queries = 0, reversed = 0, skipped = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var files = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            string[]? definitions = j.TryGetProperty("definitions", out var dd) ? [dd.GetString()!] : null;
            string spec = j.GetProperty("spec").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            bool feature = j.GetProperty("feature").GetBoolean();
            int rc = j.GetProperty("rc").GetInt32();
            bool memerr = j.TryGetProperty("memerr", out var me) && me.GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25 or -6 or -24;
            string fa = j.GetProperty("fa").GetString()!, gff = j.GetProperty("gff").GetString()!, err = j.GetProperty("err").GetString()!;
            string id = $"row {j.GetProperty("id")} [{spec} {string.Join(' ', args)}{(feature ? " -feature" : "")}]";
            try
            {
                var q = EmbossSequenceQualifiers.ParseFile(string.Join(' ', args.Where(a => a != "-firstonly")));
                var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
                var o = EmbossSeqretPrograms.Seqret(new EmbossSeqretRequest(spec, files, q, feature, args.Contains("-firstonly"), "o.fa", "fasta", "o.gff", "gff3", date) { Definitions = definitions });
                if (memerr)
                {
                    skipped++;
                    continue;
                }

                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                string messages = o.Messages.Count > 0 ? string.Join("\n", o.Messages) + "\n" : "";
                int ourRc = o.Died ? 1 : 0;
                if (ourRc != rc)
                    failures.Add($"{id}: exit code {ourRc} expected {rc}");
                else if (o.SequenceText != fa)
                    failures.Add($"{id}: sequences differ at {FirstDifference(o.SequenceText, fa)}");
                else if (feature && !o.Died && o.FeatureText != gff)
                    failures.Add($"{id}: features differ at {FirstDifference(o.FeatureText, gff)}");
                else if (messages != err)
                    failures.Add($"{id}: messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else
                {
                    compared++;
                    if (fa.Length > 0)
                        withOutput++;
                    if (spec.Contains('@') || spec.Contains("list:", StringComparison.OrdinalIgnoreCase))
                        lists++;
                    if (spec.Contains('[') || args.Contains("-sbegin1") || args.Contains("-send1"))
                        ranges++;
                    if (spec.Contains(':') && !spec.Contains("::") && !spec.StartsWith('@'))
                        queries++;
                    if (args.Contains("-sreverse1") || spec.Contains(":r]", StringComparison.OrdinalIgnoreCase))
                        reversed++;
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
        Assert.That(rows, Is.GreaterThanOrEqualTo(minRows));
        Assert.That(compared, Is.GreaterThanOrEqualTo(minCompared), $"compared {compared}, refused {refused}, skipped {skipped}");
        Assert.That(withOutput, Is.GreaterThanOrEqualTo(minOutput));
        Assert.That(lists, Is.GreaterThanOrEqualTo(minLists));
        Assert.That(ranges, Is.GreaterThanOrEqualTo(minRanges));
        Assert.That(queries, Is.GreaterThanOrEqualTo(minQueries));
        Assert.That(reversed, Is.GreaterThanOrEqualTo(minReversed));
    }
}
