using System.IO.Compression;
using System.Text.Json;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8b (B09 F184): the report FILES of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem in every <c>-rformat</c> (simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif diffseq draw restrict excel feattable motif nametable
/// regions seqtable srs summary table tagseq), with the report options (<c>-raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq</c>), several sequences per file, the sequence qualifiers <c>-sbegin -send -sreverse</c> and the
/// <c>-origfile</c> of the tandem programs, written through <see cref="EmbossReportPrograms"/> and the ONE <see cref="EmbossReportWriter"/>, compared as TEXT (byte for byte, run date faked) with the installed EMBOSS 6.6.0 on 100 random runs per
/// program and format (TestData/Alignment/emboss_report_formats_oracle.jsonl.gz, scripts/alignment/generate_emboss_report_formats_oracle.py). A format the program refuses (embl for proteins, pir / swiss for nucleotides, diffseq: 7 tags)
/// must be refused with the program's message and an empty file.
/// </summary>
[TestFixture]
public class EmbossReportFormats_Oracle_Tests
{
    private sealed record Row(int Id, string Program, string Format, string Params, string Sq, string Rq, IReadOnlyList<string> Cmd, long Epoch, string Fasta, string Report, string? Orig, string? Error);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_report_formats_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            string? Opt(string n) => j.GetProperty(n).ValueKind == JsonValueKind.Null ? null : j.GetProperty(n).GetString();
            rows.Add(new Row(j.GetProperty("id").GetInt32(), j.GetProperty("program").GetString()!, j.GetProperty("format").GetString()!, j.GetProperty("params").GetString()!, j.GetProperty("sq").GetString()!, j.GetProperty("rq").GetString()!,
                [.. j.GetProperty("cmd").EnumerateArray().Select(x => x.GetString()!)], j.GetProperty("epoch").GetInt64(), j.GetProperty("fasta").GetString()!, j.GetProperty("report").GetString()!, Opt("orig"), Opt("error")));
        }

        return rows;
    });

    // the program's output (report, original file) or the refusal
    private static (string Report, string? Orig, string? Refusal) Produce(Row r)
    {
        var request = new EmbossReportRequest(r.Fasta, r.Sq, "-rformat " + r.Format + " " + r.Rq, r.Cmd, DateTimeOffset.FromUnixTimeSeconds(r.Epoch).UtcDateTime, "s.fa", "o.rep");
        try
        {
            EmbossProgramReport p;
            switch (r.Program)
            {
                case "fuzznuc":
                    {
                        var a = r.Params.Split('|');
                        p = EmbossReportPrograms.Fuzznuc(request, a[0], int.Parse(a[1]), a[2] == "1");
                        break;
                    }
                case "fuzzpro":
                    p = EmbossReportPrograms.Fuzzpro(request, r.Params);
                    break;
                case "fuzztran":
                    {
                        var a = r.Params.Split('|');
                        p = EmbossReportPrograms.Fuzztran(request, a[0], a[1], 0);
                        break;
                    }
                case "etandem":
                    {
                        var a = r.Params.Split(',').Select(int.Parse).ToArray();
                        p = EmbossReportPrograms.Etandem(request, a[0], a[1], a[2], a[3] == 1, a[4] == 1);
                        break;
                    }
                default:
                    {
                        var a = r.Params.Split(',').Select(int.Parse).ToArray();
                        p = EmbossReportPrograms.Equicktandem(request, a[0], a[1]);
                        break;
                    }
            }

            return (p.Report, p.OriginalFile, null);
        }
        catch (EmbossReportException e)
        {
            return ("", e.OriginalFile, "Error: " + e.Message);
        }
    }

    private static void Check(string program, string format)
    {
        var rows = Oracle.Value.Where(r => r.Program == program && r.Format == format).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(100), $"{program} {format}");
        var failures = new List<string>();
        int withContent = 0;
        foreach (var r in rows)
        {
            var (report, orig, refusal) = Produce(r);
            string? expectedRefusal = r.Error?.Split('\n')[0];
            bool refused = refusal is not null;
            bool expectedRefused = expectedRefusal is not null && expectedRefusal.StartsWith("Error:", StringComparison.Ordinal) && r.Report.Length == 0;
            if (refused != expectedRefused || (refused && refusal != expectedRefusal) || report != r.Report || (r.Orig is not null && orig != r.Orig))
                failures.Add($"id {r.Id} {r.Program} {r.Format} [{r.Params}] sq=[{r.Sq}] rq=[{r.Rq}]\n--- fasta\n{r.Fasta}--- expected\n{r.Report}--- actual\n{report}{(refusal is null ? "" : "REFUSED " + refusal)}");
            if (r.Report.Length > 0)
                withContent++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(2)));
        if (rows.Any(r => r.Report.Length > 0))
            Assert.That(withContent, Is.GreaterThanOrEqualTo(rows.Count / 5), "reports with content");
    }

    public static IEnumerable<TestCaseData> Cases()
    {
        foreach (string p in new[] { "fuzznuc", "fuzzpro", "fuzztran", "etandem", "equicktandem" })
            foreach (string f in EmbossReportFormats.Names.Append("gff3").Distinct())
                yield return new TestCaseData(p, f).SetName($"{p}_{f}_EqualsEmboss");
    }

    [TestCaseSource(nameof(Cases))]
    public void ReportFile_EqualsEmboss(string program, string format) => Check(program, format);

    [Test]
    public void FormatList_IsTheProgramsList()
    {
        Assert.That(EmbossReportFormats.Names, Is.EqualTo(new[] { "simple", "embl", "genbank", "gff", "gff3", "pir", "swiss", "dasgff", "debug", "listfile", "dbmotif", "diffseq", "draw", "restrict", "excel", "feattable", "motif", "nametable", "regions", "seqtable", "srs", "summary", "table", "tagseq" }));
        Assert.That(EmbossReportFormats.TryParse("TagSeq", out var f) && f == EmbossReportFormat.TagSeq, Is.True, "case-insensitive");
        Assert.That(EmbossReportFormats.TryParse("fasta", out _), Is.False, "fasta is a sequence format, not a report format");
    }

    [Test]
    public void UnknownFormat_IsRefusedAsTheProgramDoes()
    {
        var ex = Assert.Throws<EmbossReportException>(() => EmbossReportPrograms.Fuzznuc(new EmbossReportRequest(">a\nACGTACGT\n", ReportQualifiers: "-rformat nosuch"), "ACG"));
        Assert.That(ex!.Message, Is.EqualTo("Unknown report format 'nosuch'"));
    }

    [Test]
    public void NegativeCoordinates_OfTheReverseFrames_AreNotPaddedWithTwoToThe32Blanks()
    {
        // fuzztran -frame -1 -pattern x(8,12)-L -rformat motif: a hit from 1 to -28 (the program's own coordinates, read as unsigned) makes the program write ~2^32 blanks per feature
        // and run for minutes (the oracle generator drops such runs after 8 s): the port refuses instead of exhausting memory; the formats that print the numbers only are written
        var r = new EmbossReportRequest(">Q02352 kinase\nACGTAGCTGAATTCTGAAGTCGGTCGAGGAACAGC\n", "-send 47", "-rformat motif");
        Assert.Throws<NotSupportedException>(() => EmbossReportPrograms.Fuzztran(r, "x(8,12)-L", "-1"));
        Assert.That(EmbossReportPrograms.Fuzztran(r with { ReportQualifiers = "-rformat table" }, "x(8,12)-L", "-1").Report, Does.Contain("\t") .Or.Contain("Start"));
    }

    [Test]
    public void Fuzznuc_SimpleReport_ExampleOfTheProgram()
    {
        // fuzznuc a.fa -pattern ACGT -rformat simple (installed 6.6.0): the header, one block per hit, the totals
        var r = new EmbossReportRequest(">seq1 a test desc\nACGTACGTAAACGT\n", "", "-rformat simple", ["[-sequence] s.fa", "-pattern ACGT", "-rformat simple"], new DateTime(2026, 11, 21, 1, 52, 57, DateTimeKind.Utc), "s.fa", "o.rep");
        string text = EmbossReportPrograms.Fuzznuc(r, "ACGT").Report;
        Assert.That(text, Does.StartWith("########################################\n# Program: fuzznuc\n# Rundate: Sat 21 Nov 2026 01:52:57\n# Commandline: fuzznuc\n#    [-sequence] s.fa\n"));
        Assert.That(text, Does.Contain("# Sequence: seq1     from: 1   to: 14\n# HitCount: 3\n"));
        Assert.That(text, Does.Contain("Feature: 1\nName: seq1\nStart: 1\nEnd: 4\nLength: 4\nScore: 4\nStrand: +\nPattern: pattern:ACGT\nMismatch: .\n"));
        Assert.That(text, Does.EndWith("# Total_sequences: 1\n# Total_length: 14\n# Reported_sequences: 1\n# Reported_hitcount: 3\n#---------------------------------------\n"));
    }
}
