using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3g (B09 F249): <c>maskfeat</c> of EMBOSS 6.6.0 (<c>emboss/maskfeat.c</c>: the regions of the features of the wanted types masked with a character or made lower case) equals the installed
/// <c>maskfeat -auto -sequence in.dat [-sformat1 F] -type TYPES [-tolower | -maskchar C] [-supper1 | -slower1] -outseq o.fa -osformat2 fasta</c> byte for byte on 1 200 runs
/// (TestData/Alignment/maskfeat_oracle.jsonl.gz, scripts/alignment/generate_maskfeat_oracle.py): entries with feature tables of embl, genbank, swiss, nbrf, gff3, refseqp and gff2 files and FASTA files, types
/// (single, wildcards, lists separated by blanks / commas / semicolons / bars, the default), mask characters (also a blank and none: lower case), damaged files, and the runs on which the reference dies by a fault.
/// </summary>
[TestFixture]
public class EmbossMaskfeat_Oracle_Tests
{
    [Test]
    public void Maskfeat_EqualsTheInstalledProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.maskfeat_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, masked = 0, refused = 0, died = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string text = j.GetProperty("text").GetString()!;
            string spec = j.GetProperty("spec").GetString()!;
            var args = j.GetProperty("args").EnumerateArray().Select(x => x.GetString()!).ToList();
            int rc = j.GetProperty("rc").GetInt32();
            bool refusal = j.GetProperty("hang").GetBoolean() || rc is -11 or -9 or -25;
            string fa = j.GetProperty("fa").GetString()!, err = j.GetProperty("err").GetString()!;
            string? sformat = null, type = "repeat*", maskchar = null;
            bool toLower = false, lower = false, upper = false;
            for (int i = 0; i < args.Count; i++)
            {
                switch (args[i])
                {
                    case "-sformat1":
                        sformat = args[i + 1];
                        break;
                    case "-type":
                        type = args[i + 1];
                        break;
                    case "-tolower":
                        toLower = true;
                        break;
                    case "-maskchar":
                        maskchar = args[i + 1];
                        break;
                    case "-slower1":
                        lower = true;
                        break;
                    case "-supper1":
                        upper = true;
                        break;
                }
            }

            string id = $"row {j.GetProperty("id")} [{j.GetProperty("fmt").GetString()} {string.Join(' ', args)}]";
            try
            {
                var o = EmbossMaskfeatPrograms.Maskfeat(new EmbossMaskfeatRequest(spec, new Dictionary<string, string> { ["in.dat"] = text }, new EmbossSequenceQualifiers { Format = sformat, Lower = lower, Upper = upper }, type, toLower, maskchar));
                if (refusal)
                {
                    failures.Add($"{id}: the reference dies by a fault, the port returned");
                    continue;
                }

                string messages = o.Messages.Count > 0 ? string.Join("\n", o.Messages) + "\n" : "";
                if ((o.Died ? 1 : 0) != rc)
                    failures.Add($"{id}: exit code {(o.Died ? 1 : 0)} expected {rc}; messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else if (o.SequenceText != fa)
                    failures.Add($"{id}: FASTA differs");
                else if (messages != err)
                    failures.Add($"{id}: messages '{messages.Replace("\n", "|")}' expected '{err.Replace("\n", "|")}'");
                else if (fa.Length > 0 && (fa.Contains('N') || fa.Contains('X') || fa.Any(char.IsAsciiLetterLower)))
                    masked++;
                if (rc != 0)
                    died++;
            }
            catch (InvalidOperationException) when (refusal)
            {
                refused++;
            }
            catch (Exception e)
            {
                failures.Add($"{id}: {e.GetType().Name}: {e.Message}");
            }
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(12)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(1200));
        Assert.That(masked, Is.GreaterThan(800));
        Assert.That(refused, Is.GreaterThan(10), "the runs on which the reference dies by a fault");
        Assert.That(died, Is.GreaterThan(5), "the runs the program refuses");
    }
}
