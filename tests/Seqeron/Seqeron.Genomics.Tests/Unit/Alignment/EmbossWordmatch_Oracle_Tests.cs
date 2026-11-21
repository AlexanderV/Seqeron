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

    /// <summary>
    /// A5-8c3d (B09 F225): the feature output formats of both feature files (-offormat4 / -offormat5 of every output feature format, the <c>format::file</c> form, -ofopenfile, other case, an unknown name) written with the same writers as
    /// diffseq (TestData/Alignment/wordmatch_featout_oracle.jsonl.gz, scripts/alignment/generate_wordmatch_featout_oracle.py): the alignments, the log and every file the program wrote, byte for byte.
    /// </summary>
    [Test]
    public void Wordmatch_FeatureOutputFormats_EqualTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.wordmatch_featout_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, nuc = 0, protein = 0, withFiles = 0, formats = 0;
        var formatNames = new HashSet<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            int id = j.GetProperty("id").GetInt32();
            bool isNuc = j.GetProperty("nuc").GetBoolean();
            bool dumpAlign = j.GetProperty("dumpalign").GetBoolean(), dumpFeat = j.GetProperty("dumpfeat").GetBoolean();
            var names = j.GetProperty("names").EnumerateArray().Select(x => x.GetString()!).ToList();
            string qual = j.GetProperty("qual").GetString()!;
            var cmd = new List<string> { "[-asequence] a.fa", "[-bsequence] b.fa", $"-wordsize {j.GetProperty("wordsize").GetInt32()}" };
            if (!dumpAlign) cmd.Add("-nodumpalign");
            if (!dumpFeat) cmd.Add("-nodumpfeat");
            var tokens = qual.Split(' ', StringSplitOptions.RemoveEmptyEntries);
            foreach (var tok in tokens)
            {
                if (tok.StartsWith('-') || cmd.Count == 0)
                    cmd.Add(tok);
                else
                    cmd[^1] += " " + tok;
            }

            if (dumpAlign) cmd.Add("-outfile o.match");
            cmd.AddRange([$"-aoutfeat {names[0]}", $"-boutfeat {names[1]}", "-logfile log.txt", "-auto"]);
            var request = new EmbossWordmatchRequest(j.GetProperty("fasta1").GetString()!, j.GetProperty("fasta2").GetString()!, j.GetProperty("wordsize").GetInt32(), dumpAlign, dumpFeat, cmd,
                DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime, "a.fa", "b.fa", "o.match", "match", default, qual, names[0], names[1]);
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
            foreach (System.Text.RegularExpressions.Match m in System.Text.RegularExpressions.Regex.Matches(qual, "-offormat[45] (\\S+)"))
                formatNames.Add(m.Groups[1].Value.ToLowerInvariant());
            var expected = new SortedDictionary<string, string>(StringComparer.Ordinal);
            foreach (var f in j.GetProperty("files").EnumerateObject())
                expected[f.Name] = f.Value.GetString() ?? "";
            var actual = new SortedDictionary<string, string>(StringComparer.Ordinal);
            if (got.FeatureFileName1 is { } n1) actual[n1] = got.FeatureFile1;
            if (got.FeatureFileName2 is { } n2) actual[n2] = got.FeatureFile2;
            if (expected.Count > 0) withFiles++;
            if (got.Log != j.GetProperty("log").GetString()! || (dumpAlign && got.Alignments != j.GetProperty("align").GetString()!))
                failures.Add($"row {id}: log {Diff(got.Log, j.GetProperty("log").GetString()!)} / alignments {Diff(got.Alignments, j.GetProperty("align").GetString()!)}");
            else if (!expected.Keys.SequenceEqual(actual.Keys))
                failures.Add($"row {id} ({qual}): files {string.Join(',', actual.Keys)} expected {string.Join(',', expected.Keys)}");
            else
            {
                foreach (var kv in expected)
                {
                    if (actual[kv.Key] != kv.Value)
                    {
                        failures.Add($"row {id} ({qual}) file {kv.Key}: {Diff(actual[kv.Key], kv.Value)}");
                        break;
                    }
                }
            }
        }

        formats = formatNames.Count;
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(500));
        Assert.That(nuc, Is.GreaterThan(250));
        Assert.That(protein, Is.GreaterThan(100));
        Assert.That(withFiles, Is.GreaterThan(400));
        Assert.That(formats, Is.GreaterThan(20));
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
