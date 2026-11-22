#pragma warning disable S127

using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F262, WP58): <c>featreport</c> of EMBOSS 6.6.0 (<c>emboss/featreport.c</c>, <c>ajReportWrite</c> of a feature table) through the installed program:
/// <c>featreport -auto -sequence in.dat -features f.dat -outfile o.rep [-sformat1 F] [-sbegin1 N] [-send1 N] [-sreverse1] [-fbegin2 N] [-fend2 N] [-rformat3 FORMAT] [-raccshow3] [-rdesshow3] [-rscoreshow3 N] [-rstrandshow3 N] [-rusashow3] [-rmaxall3 N] [-rmaxseq3 N]</c>.
/// TestData/Alignment/featreport_oracle.jsonl.gz, scripts/alignment/generate_featreport_oracle.py: feature files in embl, genbank, swiss, nbrf, gff3, refseqp, gff2 (the feature library readers) and EMBL entries with CDS / mRNA / exon features, as the sequence file or beside FASTA
/// sequences, damaged files, every report format of the library. The report file (the feature formats gff, gff3, embl, genbank, pir, swiss and dasgff are written by the feature writers of the library) and the messages are compared byte for
/// byte; the runs the reference cannot run are refused by the port.
/// </summary>
[TestFixture]
public class EmbossFeatreport_Oracle_Tests
{
    private static string FirstDifference(string a, string b)
    {
        int n = Math.Min(a.Length, b.Length);
        int i = 0;
        while (i < n && a[i] == b[i])
            i++;
        return $"offset {i} (got '{(i < a.Length ? a[i..Math.Min(a.Length, i + 40)] : "<end>").Replace("\n", "\\n")}', expected '{(i < b.Length ? b[i..Math.Min(b.Length, i + 40)] : "<end>").Replace("\n", "\\n")}')";
    }

    [Test]
    public void Featreport_EqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.featreport_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var formats = new HashSet<string>();
        int rows = 0, compared = 0, refused = 0, withReport = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string fmt = j.GetProperty("fmt").GetString()!;
            string text = j.GetProperty("text").GetString()!, ftext = j.GetProperty("ftext").GetString()!;
            var a = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            var cmd = j.GetProperty("cmd").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            bool memerr = j.GetProperty("memerr").GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25 or -6 or -24;
            string expected = j.GetProperty("out").GetString()!, err = j.GetProperty("err").GetString()!;
            var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
            string id = $"row {j.GetProperty("id")} [{fmt} {j.GetProperty("how").GetString()} {string.Join(' ', a)}]";
            try
            {
                var q = new List<string>();
                int featureBegin = 1;
                int? featureEnd = null;
                for (int i = 0; i < a.Count; i++)
                {
                    switch (a[i])
                    {
                        case "-auto": break;
                        case "-sequence": case "-features": case "-outfile": i++; break;
                        case "-fbegin2": featureBegin = int.Parse(a[++i]); break;
                        case "-fend2": featureEnd = int.Parse(a[++i]); break;
                        default:
                            q.Add(a[i]);
                            if (i + 1 < a.Count && !a[i + 1].StartsWith('-'))
                                q.Add(a[++i]);
                            break;
                    }
                }

                var (inputs, rq) = EmbossSequenceQualifiers.ParseInputsAndReport("featreport", string.Join(' ', q));
                var o = EmbossFeatreportProgram.Featreport(new EmbossFeatreportRequest("in.dat", "f.dat", new Dictionary<string, string> { ["in.dat"] = text, ["f.dat"] = ftext }, inputs[0], rq, featureBegin, featureEnd, cmd, date, "o.rep"));
                if (memerr)
                    continue;
                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                var got = o.Messages.SelectMany(m => m.Split('\n')).Where(l => l.Length > 0).ToList();
                var want = err.Split('\n').Where(l => l.Length > 0).ToList();
                string? why = null;
                if ((o.Died ? 1 : 0) != rc)
                    why = $"exit code {(o.Died ? 1 : 0)} expected {rc}";
                else if (!o.Died && o.Report != expected)
                    why = $"report differs at {FirstDifference(o.Report, expected)}";
                else if (!got.SequenceEqual(want))
                    why = $"messages '{string.Join("|", got)}' expected '{string.Join("|", want)}'";
                if (why is not null)
                {
                    failures.Add($"{id}: {why}");
                    continue;
                }

                compared++;
                if (expected.Length > 0)
                {
                    withReport++;
                    formats.Add(a.Contains("-rformat3") ? a[a.IndexOf("-rformat3") + 1] : "gff");
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
        Assert.That(compared, Is.GreaterThanOrEqualTo(990), $"compared {compared}, refused {refused}");
        Assert.That(withReport, Is.GreaterThanOrEqualTo(380));
        Assert.That(formats.Count, Is.GreaterThanOrEqualTo(20));
    }
}
