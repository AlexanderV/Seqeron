using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A2-28d (B09 F198): the files of EMBOSS 6.6.0 <c>wordmatch</c> (TestData/Alignment/wordmatch_oracle.jsonl.gz, scripts/alignment/generate_wordmatch_oracle.py): the alignments (-outfile in several -aformat), the statistics log
/// and the GFF3 feature files -aoutfeat / -boutfeat of nucleotide and protein pattern sets, byte for byte.
/// </summary>
[TestFixture]
public class EmbossWordmatch_Oracle_Tests
{
    [Test]
    public void Wordmatch_FilesEqualTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.wordmatch_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, nuc = 0, protein = 0, matchFormat = 0, otherFormat = 0, noDump = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            int id = j.GetProperty("id").GetInt32();
            bool isNuc = j.GetProperty("nuc").GetBoolean();
            bool dumpAlign = j.GetProperty("dumpalign").GetBoolean(), dumpFeat = j.GetProperty("dumpfeat").GetBoolean();
            string aformat = j.GetProperty("aformat").GetString()!;
            var cmd = j.GetProperty("cmd").EnumerateArray().Select(x => x.GetString()!).ToList();
            var request = new EmbossWordmatchRequest(j.GetProperty("fasta1").GetString()!, j.GetProperty("fasta2").GetString()!, j.GetProperty("wordsize").GetInt32(), dumpAlign, dumpFeat, cmd,
                DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime, "a.fa", "b.fa", "o.match", aformat);
            EmbossWordmatchOutput got;
            try
            {
                got = EmbossWordmatchPrograms.Run(request);
            }
            catch (Exception e)
            {
                failures.Add($"row {id}: {e.GetType().Name}: {e.Message} @ {e.StackTrace?.Split('\n').FirstOrDefault()?.Trim()}");
                continue;
            }

            if (isNuc) nuc++; else protein++;
            if (aformat == "match") matchFormat++; else otherFormat++;
            if (!dumpAlign || !dumpFeat) noDump++;
            string exAlign = j.GetProperty("align").GetString()!, exLog = j.GetProperty("log").GetString()!, exF1 = j.GetProperty("feat1").GetString()!, exF2 = j.GetProperty("feat2").GetString()!;
            if (got.Log != exLog)
                failures.Add($"row {id}: log differs at {Diff(got.Log, exLog)}");
            else if (got.Alignments != exAlign)
                failures.Add($"row {id} ({aformat}): alignments differ at {Diff(got.Alignments, exAlign)}");
            else if (got.FeatureFile1 != exF1 || got.FeatureFile2 != exF2)
                failures.Add($"row {id}: feature files differ: {Diff(got.FeatureFile1, exF1)} / {Diff(got.FeatureFile2, exF2)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(nuc, Is.GreaterThan(150));
        Assert.That(protein, Is.GreaterThan(40));
        Assert.That(matchFormat, Is.GreaterThan(100));
        Assert.That(otherFormat, Is.GreaterThan(100));
        Assert.That(noDump, Is.GreaterThan(20));
    }

    private static string Diff(string a, string b)
    {
        int i = 0;
        while (i < Math.Min(a.Length, b.Length) && a[i] == b[i])
            i++;
        string Ctx(string s) => i >= s.Length ? "<end>" : s.Substring(Math.Max(0, i - 30), Math.Min(s.Length - Math.Max(0, i - 30), 80)).Replace("\n", "\\n", StringComparison.Ordinal);
        return $"{i}: ours '{Ctx(a)}' expected '{Ctx(b)}'";
    }
}
