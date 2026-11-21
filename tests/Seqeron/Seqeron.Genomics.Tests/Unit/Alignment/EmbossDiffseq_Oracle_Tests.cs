using System.IO.Compression;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3 (B09 F196): the files of EMBOSS 6.6.0 <c>diffseq</c> (TestData/Alignment/diffseq_oracle.jsonl.gz, scripts/alignment/generate_diffseq_oracle.py): the 7-tag report in every <c>-rformat</c> with the report options and the
/// GFF3 feature files <c>-aoutfeat</c> / <c>-boutfeat</c> of nucleotide and protein pairs, byte for byte.
/// </summary>
[TestFixture]
public class EmbossDiffseq_Oracle_Tests
{
    private sealed record Row(int Id, bool Nuc, string Fasta1, string Fasta2, int WordSize, bool Glob, string Rq, List<string> Cmd, long Epoch, string Report, string Feat1, string Feat2, string? Error);

    private static List<Row> Load(string resource = "diffseq_oracle")
    {
        string name = $"Seqeron.Genomics.Tests.TestData.Alignment.{resource}.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows.Add(new Row(j.GetProperty("id").GetInt32(), j.GetProperty("nuc").GetBoolean(), j.GetProperty("fasta1").GetString()!, j.GetProperty("fasta2").GetString()!, j.GetProperty("wordsize").GetInt32(),
                j.GetProperty("glob").GetBoolean(), j.GetProperty("rq").GetString()!, [.. j.GetProperty("cmd").EnumerateArray().Select(x => x.GetString()!)], j.GetProperty("epoch").GetInt64(),
                j.GetProperty("report").GetString()!, j.GetProperty("feat1").GetString()!, j.GetProperty("feat2").GetString()!, j.GetProperty("error").ValueKind == JsonValueKind.Null ? null : j.GetProperty("error").GetString()));
        }

        return rows;
    }

    [Test]
    public void Diffseq_FilesEqualTheProgram_AllRows()
    {
        var rows = Load();
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(300));
        var failures = new List<string>();
        int nuc = 0, protein = 0, between = 0, glob = 0, refused = 0, formats = 0;
        var seenFormats = new HashSet<string>();
        foreach (var r in rows)
        {
            var request = new EmbossDiffseqRequest(r.Fasta1, r.Fasta2, r.WordSize, r.Glob, r.Rq, r.Cmd, DateTimeOffset.FromUnixTimeSeconds(r.Epoch).UtcDateTime, "a.fa", "b.fa", "o.rep", "a.gff", "b.gff");
            EmbossDiffseqOutput got;
            try
            {
                got = EmbossDiffseqPrograms.Run(request);
            }
            catch (Exception e) when (e is not EmbossReportException)
            {
                failures.Add($"row {r.Id} ({r.Rq}, w{r.WordSize}{(r.Glob ? " global" : "")}): {e.GetType().Name}: {e.Message}");
                continue;
            }
            catch (EmbossReportException e)
            {
                if (r.Error is null || !r.Error.Contains(e.Message, StringComparison.Ordinal))
                    failures.Add($"row {r.Id} ({r.Rq}): unexpected refusal '{e.Message}' (program: {r.Error})");
                else
                    refused++;
                continue;
            }

            if (r.Error is not null && r.Error.StartsWith("Error:", StringComparison.Ordinal))
            {
                failures.Add($"row {r.Id} ({r.Rq}): the program refuses: {r.Error}");
                continue;
            }

            seenFormats.Add(r.Rq.Split(' ')[1]);
            if (r.Nuc) nuc++; else protein++;
            if (r.Glob) glob++;
            if (r.Feat1.Contains("Insertion of", StringComparison.Ordinal) || r.Feat2.Contains("Insertion of", StringComparison.Ordinal) || r.Feat1.Contains("MISSING", StringComparison.Ordinal)) between++;      // rows with a gap (a between-feature in the report)
            if (got.Report != r.Report)
                failures.Add($"row {r.Id} ({r.Rq}, w{r.WordSize}{(r.Glob ? " global" : "")}): report differs at {FirstDifference(got.Report, r.Report)}");
            else if (got.FeatureFile1 != r.Feat1 || got.FeatureFile2 != r.Feat2)
                failures.Add($"row {r.Id} ({r.Rq}): feature file differs: {FirstDifference(got.FeatureFile1, r.Feat1)} / {FirstDifference(got.FeatureFile2, r.Feat2)}");
        }

        formats = seenFormats.Count;
        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(nuc, Is.GreaterThan(150));
        Assert.That(protein, Is.GreaterThan(60));
        Assert.That(between, Is.GreaterThan(50));
        Assert.That(glob, Is.GreaterThan(60));
        Assert.That(refused, Is.GreaterThan(15));
        Assert.That(formats, Is.GreaterThanOrEqualTo(22));
    }

    /// <summary>
    /// A5-8c3b (F216): the sequence qualifiers of both inputs (<c>-sbegin1/2 -send1/2 -sreverse1/2 -slower -supper -snucleotide -sprotein -sformat</c>) with the report formats and <c>-rusashow</c>: 260 runs
    /// (TestData/Alignment/diffseq_qualifiers_oracle.jsonl.gz, scripts/alignment/generate_diffseq_qualifiers_oracle.py), the report and both GFF3 files byte for byte, the refusals as the program refuses.
    /// </summary>
    [Test]
    public void Diffseq_SequenceQualifiers_EqualTheProgram_AllRows()
    {
        var rows = Load("diffseq_qualifiers_oracle");
        Assert.That(rows, Has.Count.GreaterThanOrEqualTo(250));
        var failures = new List<string>();
        int ranged = 0, reversed = 0, withDifferences = 0, refused = 0;
        foreach (var r in rows)
        {
            var request = new EmbossDiffseqRequest(r.Fasta1, r.Fasta2, r.WordSize, r.Glob, r.Rq, r.Cmd, DateTimeOffset.FromUnixTimeSeconds(r.Epoch).UtcDateTime, "a.fa", "b.fa", "o.rep", "a.gff", "b.gff");
            EmbossDiffseqOutput got;
            try
            {
                got = EmbossDiffseqPrograms.Run(request);
            }
            catch (EmbossReportException e)
            {
                if (r.Error is null || !r.Error.Contains(e.Message, StringComparison.Ordinal))
                    failures.Add($"row {r.Id} ({r.Rq}): unexpected refusal '{e.Message}' (program: {r.Error})");
                else
                    refused++;
                continue;
            }
            catch (Exception e)
            {
                failures.Add($"row {r.Id} ({r.Rq}): {e.GetType().Name}: {e.Message}");
                continue;
            }

            if (r.Rq.Contains("-sbegin", StringComparison.Ordinal) || r.Rq.Contains("-send", StringComparison.Ordinal)) ranged++;
            if (r.Rq.Contains("-sreverse", StringComparison.Ordinal)) reversed++;
            if (got.Result is { Differences.Count: > 0 }) withDifferences++;
            if (got.Report != r.Report)
                failures.Add($"row {r.Id} ({r.Rq}): report differs at {FirstDifference(got.Report, r.Report)}");
            else if (got.FeatureFile1 != r.Feat1 || got.FeatureFile2 != r.Feat2)
                failures.Add($"row {r.Id} ({r.Rq}): feature file differs: {FirstDifference(got.FeatureFile1, r.Feat1)} / {FirstDifference(got.FeatureFile2, r.Feat2)}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(ranged, Is.GreaterThan(150));
        Assert.That(reversed, Is.GreaterThan(30));
        Assert.That(withDifferences, Is.GreaterThan(80));
        Assert.That(refused, Is.GreaterThan(5));
    }

    private static string FirstDifference(string a, string b)
    {
        int i = 0;
        while (i < Math.Min(a.Length, b.Length) && a[i] == b[i])
            i++;
        string Ctx(string s) => i >= s.Length ? "<end>" : s.Substring(Math.Max(0, i - 30), Math.Min(s.Length - Math.Max(0, i - 30), 80)).Replace("\n", "\\n", StringComparison.Ordinal);
        return $"{i}: ours '{Ctx(a)}' expected '{Ctx(b)}'";
    }
}
