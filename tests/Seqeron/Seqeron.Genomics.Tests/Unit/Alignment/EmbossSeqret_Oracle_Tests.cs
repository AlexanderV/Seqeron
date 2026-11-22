using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3g (B09 F248 ...): the sequence READERS of EMBOSS 6.6.0 (<c>ajseqread.c</c>, <c>seqinFormatDef</c>) through the installed
/// <c>seqret -auto [-feature] -sequence in.dat [-sformat1 F] -outseq o.fa -osformat2 fasta -offormat2 gff3 -ofname2 o.gff</c> (the FASTA file of every sequence the reader returned, the GFF3 file of the feature tables, the messages
/// of the library on its error stream, byte for byte). TestData/Alignment/seqret_oracle.jsonl.gz, scripts/alignment/generate_seqret_oracle.py: 200 inputs of each of 49 input formats (the inputs are real writer outputs of the same
/// installation - header lines, line widths, numbering - and entries with feature tables converted by <c>seqret -feature</c>; 30 % damaged: a line deleted, duplicated, cut, replaced, CRLF line ends, blank lines), read with
/// <c>-sformat1 F</c>, <c>F::in.dat</c>, an upper-case name, no format at all (every format of the table is tried in turn), <c>-firstonly</c>, <c>-slower</c> / <c>-supper</c> / <c>-snucleotide</c> / <c>-sprotein</c>, with and without
/// <c>-feature</c>. The runs on which the reference dies by a fault (a GFF3 file with its features and none in it: <c>ajFeattableIsCircular</c> of a null table, valgrind: Invalid read of size 4) are refused with
/// <see cref="InvalidOperationException"/>.
/// </summary>
[TestFixture]
public class EmbossSeqret_Oracle_Tests
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
    public void SequenceReaders_EqualSeqret_AllRows()
    {
        var asm = System.Reflection.Assembly.GetExecutingAssembly();
        var names = asm.GetManifestResourceNames().Where(n => n.StartsWith("Seqeron.Genomics.Tests.TestData.Alignment.seqret_oracle", StringComparison.Ordinal) && n.EndsWith(".jsonl.gz", StringComparison.Ordinal)).Order().ToList();
        Assert.That(names, Is.Not.Empty, "the embedded seqret oracles");
        var lines = names.SelectMany(name =>
        {
            using var stream = asm.GetManifestResourceStream(name)!;
            using var gz = new GZipStream(stream, CompressionMode.Decompress);
            using var reader = new StreamReader(gz);
            return reader.ReadToEnd().Split('\n', StringSplitOptions.RemoveEmptyEntries);
        }).ToList();
        var failures = new List<string>();
        var perFormat = new Dictionary<string, int>();
        var read = new Dictionary<string, int>();
        var withFeatures = new Dictionary<string, int>();
        int rows = 0, refused = 0;
        foreach (string line in lines)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string fmt = j.GetProperty("fmt").GetString()!;
            perFormat[fmt] = perFormat.GetValueOrDefault(fmt) + 1;
            string text = j.GetProperty("text").GetString()!;
            string spec = j.GetProperty("spec").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            bool feature = !j.TryGetProperty("feature", out var fe) || fe.GetBoolean();
            int rc = j.GetProperty("rc").GetInt32();
            bool undefinedRun = j.TryGetProperty("memerr", out var me) && me.GetBoolean();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25;
            string fa = j.GetProperty("fa").GetString()!, gff = j.GetProperty("gff").GetString()!, err = j.GetProperty("err").GetString()!;
            string id = $"row {j.GetProperty("id")} [{fmt} {j.GetProperty("how").GetString()} {spec} {string.Join(' ', args)}{(feature ? "" : " (no -feature)")}]";
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

            try
            {
                var date = DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime;
                var o = EmbossSeqretPrograms.Seqret(new EmbossSeqretRequest(spec, new Dictionary<string, string> { ["in.dat"] = text },
                    new EmbossSequenceQualifiers { Format = sformat, Lower = lower, Upper = upper, Protein = protein, Nucleotide = nucleotide }, feature, firstOnly, "o.fa", "fasta", "o.gff", "gff3", date));
                if (undefinedRun)
                    continue;                               // the reference reads or writes memory it does not own (valgrind): the output has no defined value
                if (refusal)
                {
                    failures.Add($"{id}: the reference never stops or dies by a fault, the port returned");
                    continue;
                }

                string messages = o.Messages.Count > 0 ? string.Join("\n", o.Messages) + "\n" : "";
                if ((o.Died ? 1 : 0) != rc)
                    failures.Add($"{id}: exit code {(o.Died ? 1 : 0)} expected {rc}; messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else if (o.SequenceText != fa)
                    failures.Add($"{id}: FASTA differs at {FirstDifference(o.SequenceText, fa)}");
                else if (o.FeatureText != gff)
                    failures.Add($"{id}: GFF3 differs at {FirstDifference(o.FeatureText, gff)}");
                else if (messages != err)
                    failures.Add($"{id}: messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else if (fa.Length > 0)
                {
                    read[fmt] = read.GetValueOrDefault(fmt) + 1;
                    if (gff.Length > 0)
                        withFeatures[fmt] = withFeatures.GetValueOrDefault(fmt) + 1;
                }
            }
            catch (InvalidOperationException) when (refusal || undefinedRun)
            {
                refused++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        string[] formats = ["fasta", "ncbi", "pearson", "gifasta", "dbid", "gcg", "gcg8", "embl", "em", "swiss", "sw", "swissprot", "uniprot", "nbrf", "pir", "fastq", "fastq-sanger", "fastq-illumina", "fastq-solexa", "sam", "genbank", "gb",
            "ddbj", "refseq", "refseqp", "genpept", "codata", "strider", "ig", "igstrict", "iguspto", "staden", "textonly", "plain", "asis", "unknown", "gff2", "gff3", "gff", "raw", "gde", "clustal", "aln", "msf", "fitch", "mase", "biomart", "treecon", "acedb",
            "phylip", "phylipnon", "jackknifer", "hennig86", "nexus", "paup", "mega", "selex", "stockholm", "pfam", "experiment", "ace", "consed", "pdb", "pdbseq", "pdbnuc", "pdbnucseq"];
        // the readers whose reference runs die on every input that reaches their loop (ace / consed: a null table of names; hennig86: a loop over memory) give sequences only for the texts they refuse
        var fewReads = new Dictionary<string, int> { ["ace"] = 0, ["consed"] = 0, ["hennig86"] = 0 };
        foreach (string f in formats)
        {
            Assert.That(perFormat.GetValueOrDefault(f), Is.GreaterThanOrEqualTo(200), $"inputs of {f}");
            Assert.That(read.GetValueOrDefault(f), Is.GreaterThanOrEqualTo(fewReads.GetValueOrDefault(f, 60)), $"inputs of {f} that gave sequences (about half of the GFF inputs read with -feature are the runs on which the reference dies)");
        }

        foreach (string f in new[] { "embl", "genbank", "swiss", "nbrf", "gff3", "refseqp" })
            Assert.That(withFeatures.GetValueOrDefault(f), Is.GreaterThanOrEqualTo(20), $"inputs of {f} that gave features");
        Assert.That(refused, Is.GreaterThan(5), "the runs on which the reference dies by a fault");
    }
}
