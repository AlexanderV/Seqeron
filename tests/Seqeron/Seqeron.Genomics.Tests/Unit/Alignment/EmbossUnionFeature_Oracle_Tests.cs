#pragma warning disable S127

using System.IO.Compression;
using System.Text.Json;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F261, WP58): <c>union -feature</c> of EMBOSS 6.6.0 (<c>emboss/union.c</c>: <c>union_CopyFeatures</c>, the source features, the feature file) through the installed program:
/// <c>union -auto -sequence in.dat [-sformat1 F] -feature [-source] [-findoverlap] [-overlapfile ov.txt] -outseq OUT [-osformat2 F] [-offormat2 FF] [-ofname2 NAME] [-sbegin1 N] [-send1 N]</c>.
/// TestData/Alignment/union_feature_oracle.jsonl.gz, scripts/alignment/generate_union_feature_oracle.py: EMBL entries made for the program (fragments of one sequence that overlap by random lengths, with features: joins, complements, orders,
/// between positions), the entries with feature tables of the reader oracle (embl, genbank, gff3, gff2, swiss, nbrf, refseqp), damaged files; the feature formats gff3, gff, gff2, embl, genbank, swiss, pir, tagseq (unknown, written as gff3), debug,
/// default and explicit file names, ranges that trim the first sequence (the features are shifted by twice the offset, unsigned). Every file the program wrote (the sequence file, the feature file, the overlap file) and the messages of the
/// library are compared byte for byte; the runs the reference cannot run (a fault, memory errors) are refused by the port.
/// </summary>
[TestFixture]
public class EmbossUnionFeature_Oracle_Tests
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
    public void UnionFeature_EqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.union_feature_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var formats = new HashSet<string>();
        int rows = 0, compared = 0, refused = 0, withFeatures = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string fmt = j.GetProperty("fmt").GetString()!;
            string text = j.GetProperty("text").GetString()!;
            var a = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            bool memerr = j.GetProperty("memerr").GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25 or -6 or -24;
            string err = j.GetProperty("err").GetString()!;
            var expFiles = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
            string id = $"row {j.GetProperty("id")} [{fmt} {j.GetProperty("how").GetString()} {string.Join(' ', a)}]";
            try
            {
                var q = new List<string>();
                string outseq = "stdout";
                string? osf = null, off = null, ofn = null, odb = null;
                bool feature = false, source = false, findOverlap = false, overlapFile = false;
                for (int i = 0; i < a.Count; i++)
                {
                    switch (a[i])
                    {
                        case "-feature": feature = true; break;
                        case "-source": source = true; break;
                        case "-findoverlap": findOverlap = true; break;
                        case "-overlapfile": overlapFile = true; i++; break;
                        case "-outseq": outseq = a[++i]; break;
                        case "-osformat2": osf = a[++i]; break;
                        case "-offormat2": off = a[++i]; break;
                        case "-ofname2": ofn = a[++i]; break;
                        case "-osdbname2": odb = a[++i]; break;
                        default:
                            q.Add(a[i]);
                            if (a[i] is "-sformat1" or "-sbegin1" or "-send1")
                                q.Add(a[++i]);
                            break;
                    }
                }

                var quals = EmbossSequenceQualifiers.ParseFile(string.Join(' ', q));
                var o = EmbossUnionProgram.Union(new EmbossUnionRequest("in.dat", new Dictionary<string, string> { ["in.dat"] = text }, quals, findOverlap, overlapFile, outseq, osf, date) { Feature = feature, Source = source, OutFeatFormat = off, OutFeatName = ofn, OutDbName = odb, CommandLine = "-auto -sequence in.dat " + string.Join(' ', a.Select(x => x.Length == 0 || x.Contains(' ', StringComparison.Ordinal) ? "\"" + x + "\"" : x)) });
                if (memerr)
                    continue;
                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                string messages = o.Messages.Count > 0 ? string.Join("\n", o.Messages) + "\n" : "";
                string outName = outseq;
                var um = Regex.Match(outName, "^([A-Za-z0-9-]*)::?(.*)$");
                if (um.Success)
                    outName = um.Groups[2].Value;
                var got = new Dictionary<string, string>();
                if (!o.Died)
                    got[outName] = o.SequenceText;
                if (o.FeatureText.Length > 0 || o.FeatureFileName.Length > 0)
                    got[o.FeatureFileName] = o.FeatureText;
                if (overlapFile)
                    got["ov.txt"] = o.OverlapText;
                string? why = null;
                if ((o.Died ? 1 : 0) != rc)
                    why = $"exit code {(o.Died ? 1 : 0)} expected {rc}";
                else if (!o.Died)
                {
                    foreach (var kv in expFiles)
                    {
                        if (!got.TryGetValue(kv.Key, out string? g))
                            why = $"missing file {kv.Key} (got {string.Join(",", got.Keys)})";
                        else if (g != kv.Value)
                            why = $"file {kv.Key} differs at {FirstDifference(g, kv.Value)}";
                        if (why is not null)
                            break;
                    }

                    why ??= got.Where(kv => !expFiles.ContainsKey(kv.Key) && kv.Value.Length > 0).Select(kv => $"extra file {kv.Key}").FirstOrDefault();
                }

                if (why is null && messages != err)
                    why = $"messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'";
                if (why is not null)
                {
                    failures.Add($"{id}: {why}");
                    continue;
                }

                compared++;
                if (expFiles.Any(kv => kv.Key != outName && kv.Key != "ov.txt" && kv.Value.Length > 0))
                {
                    withFeatures++;
                    formats.Add(fmt);
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
        Assert.That(compared, Is.GreaterThanOrEqualTo(950), $"compared {compared}, refused {refused}");
        Assert.That(withFeatures, Is.GreaterThanOrEqualTo(700));
        Assert.That(formats.Count, Is.GreaterThanOrEqualTo(6));
    }
}
