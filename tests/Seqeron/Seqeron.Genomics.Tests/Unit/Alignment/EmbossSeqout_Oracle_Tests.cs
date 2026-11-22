using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 F252): the sequence WRITERS of EMBOSS 6.6.0 (<c>ajseqwrite.c</c>, <c>seqOutFormat</c>) and the output options of <c>seqret</c> through the installed program:
/// <c>seqret -auto -sequence in.dat [-sformat1 F] -outseq OUT -osformat2 FMT [-ossingle2 [-osextension2 E] [-osdirectory2 sub]] [-osname2 N] [-firstonly]</c>. TestData/Alignment/seqout_oracle.jsonl.gz,
/// scripts/alignment/generate_seqout_oracle.py: 200 inputs for each of 19 output formats (gcg, gcg8, nbrf, pir, codata, strider, staden, text, plain, raw, ig, gde, fitch, acedb, mase and the FASTA
/// names unknown / fasta / pearson); the inputs are the texts of the reader oracle (writer outputs of seqret in 24 formats, entries with feature tables, 15 % damaged), read with <c>-sformat1</c>, the output
/// format is named by <c>-osformat2</c> or by the <c>format::name</c> prefix of <c>-outseq</c> (also in capitals), with <c>-ossingle</c> (one file per sequence: <c>name.extension</c> in lower case, in
/// <c>-osdirectory</c>), <c>-osextension</c>, <c>-osname</c>, <c>-firstonly</c>, <c>-slower</c> / <c>-supper</c> / <c>-snucleotide</c> / <c>-sprotein</c>, an unknown output format. Every output file and the messages of
/// the library are compared byte for byte.
/// </summary>
[TestFixture]
public class EmbossSeqout_Oracle_Tests
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
    public void SequenceWriters_EqualSeqret_AllRows()
    {
        var assembly = System.Reflection.Assembly.GetExecutingAssembly();
        var failures = new List<string>();
        var perFormat = new Dictionary<string, int>();
        var written = new Dictionary<string, int>();
        var refusedBy = new Dictionary<string, int>();
        int rows = 0, singles = 0, refused = 0;
        var lines = new List<string>();
        foreach (string resource in assembly.GetManifestResourceNames().Where(n => n.Contains("TestData.Alignment.seqout_oracle", StringComparison.Ordinal)).Order(StringComparer.Ordinal))
        {
            using var stream = assembly.GetManifestResourceStream(resource)!;
            using var gz = new GZipStream(stream, CompressionMode.Decompress);
            using var reader = new StreamReader(gz);
            string? text0;
            while ((text0 = reader.ReadLine()) is not null)
                lines.Add(text0);
        }

        Assert.That(lines, Is.Not.Empty, "Embedded oracles 'seqout_oracle*.jsonl.gz' not found.");
        foreach (string line in lines)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string fmt = j.GetProperty("fmt").GetString()!;
            perFormat[fmt] = perFormat.GetValueOrDefault(fmt) + 1;
            string text = j.GetProperty("text").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            string? osformat = j.GetProperty("osformat").GetString();
            string outseq = j.GetProperty("outseq").GetString()!;
            bool single = j.GetProperty("single").GetBoolean();
            string? osext = j.GetProperty("osext").GetString();
            string? osdir = j.GetProperty("osdir").GetString();
            int rc = j.GetProperty("rc").GetInt32();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25;
            bool memerr = j.TryGetProperty("memerr", out var me) && me.GetBoolean();           // valgrind: the reference reads memory the library never wrote, its output has no defined value
            string? osdbname = j.TryGetProperty("osdbname", out var odn) && odn.ValueKind == JsonValueKind.String ? odn.GetString() : null;
            string? commandLine = j.TryGetProperty("cmd", out var cm) ? string.Join(' ', cm.EnumerateArray().Select(x => x.GetString())) : null;
            string err = j.GetProperty("err").GetString()!;
            string? sformat = null;
            bool firstOnly = false, lower = false, upper = false, protein = false, nucleotide = false;
            for (int i = 0; i < args.Count; i++)
            {
                switch (args[i])
                {
                    case "-sformat1":
                        sformat = args[i + 1];
                        break;
                    case "-firstonly":
                        firstOnly = true;
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

            string id = $"row {j.GetProperty("id")} [{fmt} {j.GetProperty("how").GetString()} {string.Join(' ', args)} -> {osformat} {outseq}{(single ? " single" : "")} {osext} {osdir}]";
            try
            {
                var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
                var o = EmbossSeqretPrograms.Seqret(new EmbossSeqretRequest("in.dat", new Dictionary<string, string> { ["in.dat"] = text },
                    new EmbossSequenceQualifiers { Format = sformat, Lower = lower, Upper = upper, Protein = protein, Nucleotide = nucleotide }, false, firstOnly, outseq, osformat, null, null, date)
                { OutSingle = single, OutDirectory = osdir, OutExtension = osext, CommandLine = commandLine, OutDbName = osdbname });
                if (memerr)
                    continue;
                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                var ours = new Dictionary<string, string>();
                var um = System.Text.RegularExpressions.Regex.Match(outseq, "^([A-Za-z0-9-]*)::?(.*)$");
                string mainName = um.Success ? um.Groups[2].Value : outseq;
                if (!o.Died)
                {
                    if (!single)
                        ours[mainName] = o.SequenceText;
                    foreach (var kv in o.SingleFiles)
                        ours[kv.Key] = kv.Value;
                }

                var expected = j.GetProperty("files").EnumerateObject().ToDictionary(p => p.Name, p => p.Value.GetString()!);
                string messages = o.Messages.Count > 0 ? string.Join("\n", o.Messages) + "\n" : "";
                if ((o.Died ? 1 : 0) != rc)
                    failures.Add($"{id}: exit code {(o.Died ? 1 : 0)} expected {rc}");
                else if (!ours.Keys.Order().SequenceEqual(expected.Keys.Order()))
                    failures.Add($"{id}: files [{string.Join(',', ours.Keys.Order())}] expected [{string.Join(',', expected.Keys.Order())}]");
                else if (expected.Keys.FirstOrDefault(k => ours[k] != expected[k]) is { } differing)
                    failures.Add($"{id}: file {differing} differs at {FirstDifference(ours[differing], expected[differing])}");
                else if (messages != err)
                    failures.Add($"{id}: messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else if (expected.Values.Any(v => v.Length > 0))
                {
                    written[fmt] = written.GetValueOrDefault(fmt) + 1;
                    if (single)
                        singles++;
                }
            }
            catch (InvalidOperationException) when (refusal || memerr)
            {
                refused++;
                refusedBy[fmt] = refusedBy.GetValueOrDefault(fmt) + 1;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(25)));
        string[] formats = ["unknown", "fasta", "pearson", "gcg", "gcg8", "nbrf", "pir", "codata", "strider", "staden", "text", "plain", "raw", "ig", "iguspto", "gde", "fitch", "acedb", "mase", "ncbi", "gifasta", "nib", "nibble", "dasdna", "das", "asn1", "fastq-sanger", "fastq-illumina", "fastq-solexa", "fastq", "sam", "gff2", "gff3", "gff", "msf", "clustal", "aln", "selex", "phylip", "phylipnon", "phylip3", "hennig86", "mega", "meganon", "nexus", "paup", "nexusnon", "paupnon", "jackknifer", "jackknifernon", "treecon", "scf"];
        foreach (string f in formats)
        {
            Assert.That(perFormat.GetValueOrDefault(f), Is.GreaterThanOrEqualTo(200), $"inputs of {f}");
            if (f == "iguspto")
                Assert.That(refusedBy.GetValueOrDefault(f), Is.GreaterThanOrEqualTo(100), "iguspto: the writer of the library dies by a segmentation fault (a null Fulldesc) on every sequence");
            else
                Assert.That(written.GetValueOrDefault(f), Is.GreaterThanOrEqualTo(100), $"inputs of {f} that gave files");
        }

        Assert.That(singles, Is.GreaterThan(200), "the -ossingle runs");
    }
}
