#pragma warning disable S127

using System.IO.Compression;
using System.Text;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 WP60, F275 and following): the sequence DATABASES whose access method reads an INDEX made by the EMBOSS 6.6.0 indexers (<c>dbxfasta</c>, <c>dbxflat</c>, <c>dbxgcg</c> for the method <c>emboss</c>; <c>dbifasta</c>, <c>dbiflat</c>, <c>dbigcg</c> for
/// <c>emblcd</c> and <c>gcg</c>) queried through the installed <c>seqret -auto [-feature] -sequence SPEC ...</c> and <c>entret -auto -sequence SPEC -outfile o.fa</c> (<c>EMBOSSRC</c> names the work directory; the files of a row are the database files and the index
/// files the indexer made, as bytes in a latin-1 text). TestData/Alignment/dbindex_*_oracle.jsonl.gz, scripts/alignment/generate_dbindex_oracle.py. The FASTA file, the GFF3 file, the entries and the messages of the library are compared byte for byte; the runs on
/// which the reference never stops or dies by a fault are refused with <see cref="InvalidOperationException"/>.
/// </summary>
[TestFixture]
public class EmbossDbIndex_Oracle_Tests
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
    public void MethodEmboss_EqualsSeqretAndEntret_AllRows() => CompareRows("Seqeron.Genomics.Tests.TestData.Alignment.dbindex_emboss_oracle.jsonl.gz", 1000, 1000, 380);

    [Test]
    public void MethodEmblcd_EqualsSeqretAndEntret_AllRows() => CompareRows("Seqeron.Genomics.Tests.TestData.Alignment.dbindex_emblcd_oracle.jsonl.gz", 1000, 1000, 300);

    private static void CompareRows(string name, int minRows, int minCompared, int minOutput)
    {
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, compared = 0, refused = 0, withOutput = 0, queries = 0, wild = 0, entrets = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            var files = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            string[] definitions = [j.GetProperty("definitions").GetString()!];
            string spec = j.GetProperty("spec").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            bool feature = j.GetProperty("feature").GetBoolean();
            string program = j.GetProperty("program").GetString()!;
            int rc = j.GetProperty("rc").GetInt32();
            bool memerr = j.TryGetProperty("memerr", out var me) && me.GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25 or -6 or -24;
            string fa = j.GetProperty("fa").GetString()!, gff = j.GetProperty("gff").GetString()!, err = j.GetProperty("err").GetString()!;
            string id = $"row {j.GetProperty("id")} [{j.GetProperty("how").GetString()} {program} {spec} {string.Join(' ', args)}{(feature ? " -feature" : "")}]";
            try
            {
                var q = EmbossSequenceQualifiers.ParseFile(string.Join(' ', args.Where(a => a != "-firstonly")));
                var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
                string sequenceText;
                string featureText = "";
                IReadOnlyList<string> messages;
                bool died;
                if (program == "entret")
                {
                    var e = EmbossEntretProgram.Entret(new EmbossEntretRequest(spec, files, q, args.Contains("-firstonly")) { Definitions = definitions });
                    sequenceText = e.Text;
                    messages = e.Messages;
                    died = e.Died;
                }
                else
                {
                    var o = EmbossSeqretPrograms.Seqret(new EmbossSeqretRequest(spec, files, q, feature, args.Contains("-firstonly"), "o.fa", "fasta", "o.gff", "gff3", date) { Definitions = definitions });
                    sequenceText = o.SequenceText;
                    featureText = o.FeatureText;
                    messages = o.Messages;
                    died = o.Died;
                }

                if (memerr)
                    continue;
                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                string text = messages.Count > 0 ? string.Join("\n", messages) + "\n" : "";
                int ourRc = died ? 1 : 0;
                if (ourRc != rc)
                    failures.Add($"{id}: exit code {ourRc} expected {rc} ({text.Replace("\n", "|")}, expected '{err.Replace("\n", "|")}')");
                else if (sequenceText != fa)
                    failures.Add($"{id}: sequences differ at {FirstDifference(sequenceText, fa)}");
                else if (feature && !died && featureText != gff)
                    failures.Add($"{id}: features differ at {FirstDifference(featureText, gff)}");
                else if (text != err)
                    failures.Add($"{id}: messages '{text.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else
                {
                    compared++;
                    if (fa.Length > 0)
                        withOutput++;
                    if (spec.Contains('-') || spec.Contains("::"))
                        queries++;
                    if (spec.Contains('*') || spec.Contains('?'))
                        wild++;
                    if (program == "entret")
                        entrets++;
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
        Assert.That(compared, Is.GreaterThanOrEqualTo(minCompared), $"compared {compared}, refused {refused}");
        Assert.That(withOutput, Is.GreaterThanOrEqualTo(minOutput));
        Assert.That(queries, Is.GreaterThan(0));
        Assert.That(wild, Is.GreaterThan(0));
        Assert.That(entrets, Is.GreaterThan(0));
        _ = Encoding.UTF8;
    }
}
