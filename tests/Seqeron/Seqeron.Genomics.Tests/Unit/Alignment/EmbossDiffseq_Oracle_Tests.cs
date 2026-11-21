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

    /// <summary>
    /// A5-8c3b (F218): GenBank and EMBL entries WITH FEATURES as the inputs of diffseq: the <c>first_feature</c> / <c>second_feature</c> notes (type, location, qualifiers) of the 7-tag report and of the other report formats, with the sequence qualifiers
    /// (the feature table is reversed and trimmed with the sequence), mixed with FASTA inputs, <c>-sformat1</c> (TestData/Alignment/diffseq_features_oracle.jsonl.gz, scripts/alignment/generate_diffseq_features_oracle.py, 420 runs; the GenBank entries
    /// also reproduce the program's loss of the features when the SOURCE block is followed directly by FEATURES): the report and both GFF3 files byte for byte.
    /// </summary>
    /// <remarks>
    /// A5-8c3c (F219): the second file (TestData/Alignment/diffseq_special_tags_oracle.jsonl.gz, the same generator with <c>special</c>, 450 runs) gives every feature 1 - 4 qualifiers of the tags with a special validator
    /// (anticodon, bio_material, citation, codon, collection_date, compare, cons_splice, estimated_length, inference, lat_lon, mobile_element, ncRNA_class, PCR_primers, replace, rpt_unit, rpt_unit_range, rpt_unit_seq, tag_peptide,
    /// transl_except, translation, protein_id, db_xref; valid, correctable, bad and other-case values, quoted and unquoted): the bad ones are notes (<c>*name: value</c>), the corrected ones are printed corrected, and
    /// <c>rpt_unit</c> / <c>rpt_unit_seq</c> share ONE compiled pattern of the program (the first of the two that is checked decides for both entries of the run).
    /// </remarks>
    [TestCase("diffseq_features_oracle.jsonl.gz", false)]
    [TestCase("diffseq_special_tags_oracle.jsonl.gz", true)]
    public void Diffseq_FeatureTableInputs_EqualTheProgram_AllRows(string resource, bool special)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + resource;
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        int rows = 0, withFeatures = 0, reversed = 0, trimmed = 0, genbank = 0, embl = 0, complement = 0, join = 0, badValues = 0, specialPrinted = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string rq = j.GetProperty("rq").GetString()!;
            var cmd = j.GetProperty("cmd").EnumerateArray().Select(x => x.GetString()!).ToList();
            string file1 = j.GetProperty("file1").GetString()!, file2 = j.GetProperty("file2").GetString()!;
            var request = new EmbossDiffseqRequest(file1, file2, j.GetProperty("wordsize").GetInt32(), j.GetProperty("glob").GetBoolean(), rq, cmd, DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime,
                j.GetProperty("name1").GetString()!, j.GetProperty("name2").GetString()!, "o.rep", "a.gff", "b.gff");
            string report = j.GetProperty("report").GetString()!;
            EmbossDiffseqOutput got;
            try
            {
                got = EmbossDiffseqPrograms.Run(request);
            }
            catch (EmbossReportException e) when (j.GetProperty("error").ValueKind != JsonValueKind.Null && j.GetProperty("error").GetString()!.Contains(e.Message, StringComparison.Ordinal))
            {
                continue;
            }
            catch (Exception e)
            {
                failures.Add($"row {j.GetProperty("id")} ({rq}): {e.GetType().Name}: {e.Message}");
                continue;
            }

            if (report.Contains("Feature: ", StringComparison.Ordinal) || report.Contains("first_feature", StringComparison.Ordinal) || report.Contains("second_feature", StringComparison.Ordinal))
                withFeatures++;
            if (report.Contains("note='*", StringComparison.Ordinal)) badValues++;
            if (report.Contains("lat_lon='", StringComparison.Ordinal) || report.Contains("rpt_unit_seq='", StringComparison.Ordinal) || report.Contains("inference='", StringComparison.Ordinal) || report.Contains("collection_date='", StringComparison.Ordinal))
                specialPrinted++;
            if (rq.Contains("-sreverse", StringComparison.Ordinal)) reversed++;
            if (rq.Contains("-sbegin", StringComparison.Ordinal) || rq.Contains("-send", StringComparison.Ordinal)) trimmed++;
            if (file1.StartsWith("LOCUS", StringComparison.Ordinal)) genbank++; else embl++;
            if (file1.Contains("complement(", StringComparison.Ordinal) || file2.Contains("complement(", StringComparison.Ordinal)) complement++;
            if (file1.Contains("join(", StringComparison.Ordinal) || file2.Contains("join(", StringComparison.Ordinal)) join++;
            if (got.Report != report)
                failures.Add($"row {j.GetProperty("id")} ({rq}): report differs at {FirstDifference(got.Report, report)}");
            else if (got.FeatureFile1 != j.GetProperty("feat1").GetString() || got.FeatureFile2 != j.GetProperty("feat2").GetString())
                failures.Add($"row {j.GetProperty("id")} ({rq}): feature file differs");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(400));
        Assert.That(withFeatures, Is.GreaterThan(150), "reports with feature notes");
        Assert.That(reversed, Is.GreaterThan(40));
        Assert.That(trimmed, Is.GreaterThan(80));
        Assert.That(genbank, Is.GreaterThan(100));
        Assert.That(embl, Is.GreaterThan(60));
        Assert.That(complement, Is.GreaterThan(80));
        Assert.That(join, Is.GreaterThan(100));
        if (special)
        {
            Assert.That(badValues, Is.GreaterThan(100), "reports with a bad special value (a note)");
            Assert.That(specialPrinted, Is.GreaterThan(60), "reports with a valid special value");
        }
    }

    /// <summary>
    /// A5-8c3c (F220): the feature output formats of both feature files (<c>-offormat4</c> / <c>-offormat5</c>: gff, gff2, gff3, bed, embl, em, genbank, gb, ddbj, refseq, refseqp, pir, nbrf, swiss, sw, swissprot, uniprot, dasgff, draw, debug, unknown,
    /// other case, a name that is no format; the <c>format::file</c> form of <c>-aoutfeat</c>; <c>-ofopenfile</c>; <c>-ofname</c> / <c>-ofextension</c> / <c>-ofsingle</c>) on nucleotide and protein pairs with the sequence qualifiers and the report formats
    /// (TestData/Alignment/diffseq_featout_oracle.jsonl.gz, scripts/alignment/generate_diffseq_featout_oracle.py, 500 runs): the report (its header names the files) and every file the program wrote, byte for byte, including the empty file
    /// of a format the program refuses for the kind of sequence and the text one writer leaves in the static string of the library for the next writer (gff2, pir print it in front of their first tag; the report is written after the feature files).
    /// </summary>
    [Test]
    public void Diffseq_FeatureOutputFormats_EqualTheProgram_AllRows()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.diffseq_featout_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var failures = new List<string>();
        var formats = new HashSet<string>(StringComparer.Ordinal);
        int rows = 0, protein = 0, nucleotide = 0, openfile = 0, prefixed = 0, refused = 0;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            rows++;
            string rq = j.GetProperty("rq").GetString()!;
            var cmd = j.GetProperty("cmd").EnumerateArray().Select(x => x.GetString()!).ToList();
            var names = j.GetProperty("names").EnumerateArray().Select(x => x.GetString()!).ToArray();
            var request = new EmbossDiffseqRequest(j.GetProperty("fasta1").GetString()!, j.GetProperty("fasta2").GetString()!, j.GetProperty("wordsize").GetInt32(), j.GetProperty("glob").GetBoolean(), rq, cmd,
                DateTimeOffset.FromUnixTimeSeconds(j.GetProperty("epoch").GetInt64()).UtcDateTime, "a.fa", "b.fa", "o.rep", names[0], names[1]);
            string? expectedError = j.GetProperty("error").ValueKind == JsonValueKind.Null ? null : j.GetProperty("error").GetString();
            EmbossDiffseqOutput got;
            try
            {
                got = EmbossDiffseqPrograms.Run(request);
            }
            catch (EmbossReportException e) when (expectedError is not null && expectedError.Contains(e.Message, StringComparison.Ordinal))
            {
                continue;
            }
            catch (Exception e)
            {
                failures.Add($"row {j.GetProperty("id")} ({rq}): {e.GetType().Name}: {e.Message}");
                continue;
            }

            foreach (System.Text.RegularExpressions.Match m in System.Text.RegularExpressions.Regex.Matches(rq, "-offormat[45] (\\S+)"))
                formats.Add(m.Groups[1].Value.ToLowerInvariant());
            if (j.GetProperty("nuc").GetBoolean()) nucleotide++; else protein++;
            if (rq.Contains("-ofopenfile", StringComparison.Ordinal)) openfile++;
            if (names.Any(n => n.Contains("::", StringComparison.Ordinal))) prefixed++;
            var expected = new SortedDictionary<string, string>(StringComparer.Ordinal);
            foreach (var p in j.GetProperty("files").EnumerateObject())
                expected[p.Name] = p.Value.ValueKind == JsonValueKind.Null ? "" : p.Value.GetString()!;
            var actual = new SortedDictionary<string, string>(StringComparer.Ordinal);
            if (got.FeatureFileName1 is not null) actual[got.FeatureFileName1] = got.FeatureFile1;
            if (got.FeatureFileName2 is not null) actual[got.FeatureFileName2] = got.FeatureFile2;
            refused += expected.Values.Count(v => v.Length == 0);
            string report = j.GetProperty("report").GetString()!;
            if (got.Report != report)
                failures.Add($"row {j.GetProperty("id")} ({rq}): report differs at {FirstDifference(got.Report, report)}");
            else if (!expected.Keys.SequenceEqual(actual.Keys))
                failures.Add($"row {j.GetProperty("id")} ({rq}): files {string.Join(',', actual.Keys)} instead of {string.Join(',', expected.Keys)}");
            else
                foreach (var key in expected.Keys.Where(k => expected[k] != actual[k]))
                    failures.Add($"row {j.GetProperty("id")} ({rq}): file {key} differs at {FirstDifference(actual[key], expected[key])}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(rows, Is.GreaterThanOrEqualTo(500));
        Assert.That(formats.Count, Is.GreaterThanOrEqualTo(20), "output formats and names tried");
        Assert.That(nucleotide, Is.GreaterThan(200));
        Assert.That(protein, Is.GreaterThan(100));
        Assert.That(openfile, Is.GreaterThan(10));
        Assert.That(prefixed, Is.GreaterThan(30));
        Assert.That(refused, Is.GreaterThan(100), "empty files");
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
