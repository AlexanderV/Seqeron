using System.IO.Compression;
using System.Text.Json;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F253, F254): EMBOSS 6.6.0 <c>union</c> (<c>emboss/union.c</c>: the sequences of a file joined into one, with the overlap of consecutive sequences found and counted) and <c>seqretsplit</c>
/// (<c>emboss/seqretsplit.c</c>: seqret that writes one file per sequence and renames duplicate names) through the installed programs. TestData/Alignment/seqmerge_oracle.jsonl.gz,
/// scripts/alignment/generate_seqmerge_oracle.py: 200 runs of each program on fragments of one sequence that overlap by random lengths (union) or sequences with repeated names in any case (seqretsplit) and on the
/// inputs of the reader oracle (15 % damaged), with <c>-sformat1</c>, <c>-slower</c> / <c>-supper</c> / <c>-snucleotide</c> / <c>-sprotein</c>, any of the 14 output formats (<c>-osformat2</c> or the <c>format::</c>
/// prefix of the output name), <c>-findoverlap</c>, <c>-overlapfile</c>, <c>-firstonly</c>, <c>-ossingle2</c> / <c>-noossingle2</c>, <c>-osextension2</c>, <c>-osdirectory2</c>. Every file and the messages of the library are
/// compared byte for byte.
/// </summary>
[TestFixture]
public class EmbossSeqmerge_Oracle_Tests
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
    public void UnionAndSeqretsplit_EqualTheInstalledPrograms_AllRows()
    {
        var asm = System.Reflection.Assembly.GetExecutingAssembly();
        var lines = asm.GetManifestResourceNames().Where(n => n.StartsWith("Seqeron.Genomics.Tests.TestData.Alignment.seqmerge_oracle", StringComparison.Ordinal)).Order().SelectMany(n =>
        {
            using var stream = asm.GetManifestResourceStream(n)!;
            using var gz = new GZipStream(stream, CompressionMode.Decompress);
            using var reader = new StreamReader(gz);
            return reader.ReadToEnd().Split('\n', StringSplitOptions.RemoveEmptyEntries);
        }).ToList();
        var failures = new List<string>();
        var perProgram = new Dictionary<string, int>();
        var withFiles = new Dictionary<string, int>();
        int overlaps = 0, renamed = 0;
        foreach (string line in lines)
        {
            var j = JsonDocument.Parse(line).RootElement;
            string prog = j.GetProperty("prog").GetString()!;
            perProgram[prog] = perProgram.GetValueOrDefault(prog) + 1;
            string text = j.GetProperty("text").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            string? osformat = j.GetProperty("osformat").GetString();
            string outseq = j.GetProperty("outseq").GetString()!;
            int rc = j.GetProperty("rc").GetInt32();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25;
            string err = j.GetProperty("err").GetString()!;
            string? sformat = null;
            bool lower = false, upper = false, protein = false, nucleotide = false;
            for (int i = 0; i < args.Count; i++)
            {
                switch (args[i])
                {
                    case "-sformat1":
                        sformat = args[i + 1];
                        break;
                    case "-slower1":
                        lower = true;
                        break;
                    case "-supper1":
                        upper = true;
                        break;
                    case "-sprotein1":
                        protein = true;
                        break;
                    case "-snucleotide1":
                        nucleotide = true;
                        break;
                }
            }

            var q = new EmbossSequenceQualifiers { Format = sformat, Lower = lower, Upper = upper, Protein = protein, Nucleotide = nucleotide };
            var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
            string id = $"row {j.GetProperty("id")} [{prog} {j.GetProperty("how").GetString()} {string.Join(' ', args)} -> {osformat} {outseq}]";
            var um = Regex.Match(outseq, "^([A-Za-z0-9-]*)::?(.*)$");
            string mainName = um.Success ? um.Groups[2].Value : outseq;
            try
            {
                var ours = new Dictionary<string, string>();
                IReadOnlyList<string> messages;
                bool died;
                if (prog == "union")
                {
                    bool overlapFile = j.GetProperty("overlapfile").ValueKind != JsonValueKind.Null;
                    var o = EmbossUnionProgram.Union(new EmbossUnionRequest("in.dat", new Dictionary<string, string> { ["in.dat"] = text }, q, j.GetProperty("findoverlap").GetBoolean(), overlapFile, outseq, osformat, date));
                    died = o.Died;
                    messages = o.Messages;
                    if (!died)
                    {
                        ours[mainName] = o.SequenceText;
                        if (overlapFile)
                        {
                            ours["ov.txt"] = o.OverlapText;
                            if (o.OverlapText.Split('\n').Any(l => l.Length > 0 && l != "0"))
                                overlaps++;
                        }
                    }
                }
                else
                {
                    bool? single = j.GetProperty("single").ValueKind == JsonValueKind.Null ? null : j.GetProperty("single").GetBoolean();
                    string? osext = j.GetProperty("osext").GetString();
                    string? osdir = j.GetProperty("osdir").GetString();
                    var o = EmbossSeqretPrograms.Seqretsplit(new EmbossSeqretRequest("in.dat", new Dictionary<string, string> { ["in.dat"] = text }, q, false, j.GetProperty("firstonly").GetBoolean(), outseq, osformat, null, null, date)
                    { OutSingle = single, OutDirectory = osdir, OutExtension = osext });
                    died = o.Died;
                    messages = o.Messages;
                    if (!died)
                    {
                        if (!(single ?? true))
                            ours[osdir is null ? mainName : Path.Combine(osdir, mainName)] = o.SequenceText;
                        foreach (var kv in o.SingleFiles)
                            ours[kv.Key] = kv.Value;
                    }

                    if (o.Messages.Any(m => m.StartsWith("Warning: Duplicate name", StringComparison.Ordinal)))
                        renamed++;
                }

                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                var expected = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
                string text2 = messages.Count > 0 ? string.Join("\n", messages) + "\n" : "";
                if ((died ? 1 : 0) != rc)
                    failures.Add($"{id}: exit code {(died ? 1 : 0)} expected {rc}; messages '{text2.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else if (!ours.Keys.Order().SequenceEqual(expected.Keys.Order()))
                    failures.Add($"{id}: files [{string.Join(',', ours.Keys.Order())}] expected [{string.Join(',', expected.Keys.Order())}]");
                else if (expected.Keys.FirstOrDefault(k => ours[k] != expected[k]) is { } differing)
                    failures.Add($"{id}: file {differing} differs at {FirstDifference(ours[differing], expected[differing])}");
                else if (text2 != err)
                    failures.Add($"{id}: messages '{text2.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else if (expected.Values.Any(v => v.Length > 0))
                    withFiles[prog] = withFiles.GetValueOrDefault(prog) + 1;
            }
            catch (InvalidOperationException) when (refusal)
            {
                // the reference dies by a fault on this input
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count}:\n" + string.Join("\n", failures.Take(12)));
        foreach (string prog in new[] { "union", "seqretsplit" })
        {
            Assert.That(perProgram.GetValueOrDefault(prog), Is.GreaterThanOrEqualTo(200), $"runs of {prog}");
            Assert.That(withFiles.GetValueOrDefault(prog), Is.GreaterThanOrEqualTo(150), $"runs of {prog} that gave files");
        }

        Assert.That(overlaps, Is.GreaterThan(30), "union runs that found an overlap");
        Assert.That(renamed, Is.GreaterThan(30), "seqretsplit runs that renamed a duplicate");
    }
}
