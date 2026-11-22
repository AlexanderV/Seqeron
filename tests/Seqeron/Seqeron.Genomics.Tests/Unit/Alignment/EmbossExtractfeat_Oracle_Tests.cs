#pragma warning disable S127
using System.Globalization;
using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F258, WP58): <c>extractfeat</c> of EMBOSS 6.6.0 (<c>emboss/extractfeat.c</c>) through the installed program:
/// <c>extractfeat -auto -sequence in.dat [-sformat1 F] [-before N] [-after N] [-join] [-featinname] [-describe T] [-source S] [-type T] [-sense N] [-minscore X] [-maxscore X] [-tag T] [-value V] -outseq o.out -osformat2 OF</c>.
/// TestData/Alignment/extractfeat_oracle.jsonl.gz, scripts/alignment/generate_extractfeat_oracle.py: entries with feature tables in the formats whose readers carry features (embl, genbank, swiss, nbrf, gff3, refseqp, gff2), FASTA files without
/// features, damaged files; every filter, the regions before / after (padding, negative values), join, the names and descriptions, 19 output formats, the dates of the headers (faked). The output file and the messages of the library are compared
/// byte for byte; the runs the reference cannot run (a fault, a loop, memory errors) are refused by the port.
/// </summary>
[TestFixture]
public class EmbossExtractfeat_Oracle_Tests
{
    private static string FirstDifference(string a, string b)
    {
        int n = Math.Min(a.Length, b.Length);
        int i = 0;
        while (i < n && a[i] == b[i])
            i++;
        return $"offset {i} (got '{(i < a.Length ? a[i..Math.Min(a.Length, i + 30)] : "<end>").Replace("\n", "\\n")}', expected '{(i < b.Length ? b[i..Math.Min(b.Length, i + 30)] : "<end>").Replace("\n", "\\n")}')";
    }

    // acdArgSave (ajUtilGetCmdline): an argument that is empty, has a blank or one of the characters *?[]{}|!&^ is in double quotes
    private static string Quote(string x) => x.Length == 0 || x.Any(char.IsWhiteSpace) || x.IndexOfAny(['*', '?', '[', ']', '{', '}', '|', '!', '&', '^']) >= 0 ? "\"" + x + "\"" : x;

    private static EmbossExtractfeatRequest Request(string spec, string text, List<string> a, DateTime date)
    {
        var q = new List<string>();
        int before = 0, after = 0, sense = 0;
        bool join = false, featInName = false;
        string dbname = "", describe = "", source = "*", type = "*", tag = "*", value = "*", osf = "fasta";
        string? osdir = null, osext = null;
        bool single = false;
        float minscore = 0, maxscore = 0;
        for (int i = 0; i < a.Count; i++)
        {
            switch (a[i])
            {
                case "-before": before = int.Parse(a[++i], CultureInfo.InvariantCulture); break;
                case "-after": after = int.Parse(a[++i], CultureInfo.InvariantCulture); break;
                case "-join": join = true; break;
                case "-featinname": featInName = true; break;
                case "-describe": describe = a[++i]; break;
                case "-source": source = a[++i]; break;
                case "-type": type = a[++i]; break;
                case "-sense": sense = int.Parse(a[++i], CultureInfo.InvariantCulture); break;
                case "-minscore": minscore = float.Parse(a[++i], CultureInfo.InvariantCulture); break;
                case "-maxscore": maxscore = float.Parse(a[++i], CultureInfo.InvariantCulture); break;
                case "-tag": tag = a[++i]; break;
                case "-value": value = a[++i]; break;
                case "-osformat2": osf = a[++i]; break;
                case "-osdbname2": dbname = a[++i]; break;
                case "-ossingle2": single = true; break;
                case "-osextension2": osext = a[++i]; break;
                case "-osdirectory2": osdir = a[++i]; break;
                case "-osname2": i++; break;
                default:
                    q.Add(a[i]);
                    if (a[i] is "-sformat1" or "-sbegin1")
                        q.Add(a[++i]);
                    break;
            }
        }

        var quals = EmbossSequenceQualifiers.ParseFile(string.Join(' ', q));
        return new EmbossExtractfeatRequest(spec, new Dictionary<string, string> { ["in.dat"] = text }, quals, before, after, join, featInName, describe, source, type, sense, minscore, maxscore, tag, value, "o.out", osf) { RunDate = date, OutDbName = dbname.Length > 0 ? dbname : null, OutSingle = single, OutExtension = osext, OutDirectory = osdir, ExistingDirectories = ["sub"], CommandLine = "-auto -sequence " + spec + " " + string.Join(' ', a.Take(a.Count - 2).Select(Quote)) + " -outseq o.out -osformat2 " + osf };
    }

    // the files the program wrote: the output file (created when the run starts, unless -ossingle) and the files of -ossingle
    private static Dictionary<string, string> Written(EmbossSeqretOutput o, List<string> args)
    {
        var files = new Dictionary<string, string>(o.SingleFiles);
        if (!args.Contains("-ossingle2"))
            files[o.SequenceFileName] = o.SequenceText;
        return files;
    }

    [Test]
    public void Extractfeat_EqualsTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.extractfeat_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var formats = new HashSet<string>();
        var outFormats = new HashSet<string>();
        int rows = 0, compared = 0, refused = 0, withOutput = 0, singles = 0, directories = 0, singlesWritten = 0, directoriesWritten = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string fmt = j.GetProperty("fmt").GetString()!;
            string text = j.GetProperty("text").GetString()!;
            string spec = j.GetProperty("spec").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            bool memerr = j.GetProperty("memerr").GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25 or -6 or -24;
            string err = j.GetProperty("err").GetString()!;
            var expectedFiles = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
            string expected = string.Concat(expectedFiles.Values);
            var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
            string id = $"row {j.GetProperty("id")} [{fmt} {j.GetProperty("how").GetString()} {spec} {string.Join(' ', args)}]";
            try
            {
                var o = EmbossExtractfeatProgram.Extractfeat(Request(spec, text, args, date));
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
                else if (!o.Died && !Written(o, args).Keys.Order(StringComparer.Ordinal).SequenceEqual(expectedFiles.Keys.Order(StringComparer.Ordinal)))
                    failures.Add($"{id}: files [{string.Join(',', Written(o, args).Keys.Order(StringComparer.Ordinal))}] expected [{string.Join(',', expectedFiles.Keys.Order(StringComparer.Ordinal))}]");
                else if (!o.Died && Written(o, args).Where(kv => kv.Value != expectedFiles[kv.Key]).Select(kv => $"file {kv.Key} differs at {FirstDifference(kv.Value, expectedFiles[kv.Key])}").FirstOrDefault() is { } fileDiff)
                    failures.Add($"{id}: {fileDiff}");
                else if (messages != err)
                    failures.Add($"{id}: messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else
                {
                    compared++;
                    if (args.Contains("-ossingle2"))
                    {
                        singles++;
                        if (expected.Length > 0)
                            singlesWritten++;
                    }

                    if (args.Contains("-osdirectory2"))
                    {
                        directories++;
                        if (expected.Length > 0)
                            directoriesWritten++;
                    }
                    if (expected.Length > 0)
                    {
                        withOutput++;
                        formats.Add(fmt);
                        outFormats.Add(args[args.IndexOf("-osformat2") + 1]);
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
        Assert.That(rows, Is.GreaterThanOrEqualTo(1000));
        Assert.That(compared, Is.GreaterThanOrEqualTo(900), $"compared {compared}, refused {refused}");
        Assert.That(withOutput, Is.GreaterThanOrEqualTo(350));
        Assert.That(singles, Is.GreaterThanOrEqualTo(100), "the -ossingle runs");
        Assert.That(directories, Is.GreaterThanOrEqualTo(100), "the -osdirectory runs");
        Assert.That(singlesWritten, Is.GreaterThanOrEqualTo(30), "the -ossingle runs that wrote files");
        Assert.That(directoriesWritten, Is.GreaterThanOrEqualTo(40), "the -osdirectory runs that wrote files");
        Assert.That(formats.Count, Is.GreaterThanOrEqualTo(6));
        Assert.That(outFormats.Count, Is.GreaterThanOrEqualTo(15));
    }
}
