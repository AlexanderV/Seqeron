using System.IO.Compression;
using System.Text.Json;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c (B09 F192): the report files of fuzznuc, fuzzpro and fuzztran for FASTA headers with NCBI-style identifiers (<c>ajSeqParseNcbi</c>: gb| emb| dbj| sp| ref| tpg| gi| gnl| gnl|BL_ORD_ID| lcl| bbs| pir|| prf|| pdb| pat|
/// and unknown prefixes with 2 to 5 barred tokens; names with ':' ',' '/' '\' through <c>seqSetName</c>) and for pattern FILES (<c>-pattern @file</c>: one pattern per line or FASTA style with <c>&lt;mismatch=n&gt;</c>,
/// <c>-pname</c>, <c>-pformat</c>, <c>-pmismatch</c>), in every report format, compared as text with the installed EMBOSS 6.6.0
/// (TestData/Alignment/emboss_ncbi_pattern_oracle.jsonl.gz, scripts/alignment/generate_emboss_ncbi_pattern_oracle.py).
/// </summary>
[TestFixture]
public class EmbossNcbiPattern_Oracle_Tests
{
    private sealed record Row(int Id, string Kind, string Program, string Format, string Params, string? PatFile, string PArgs, string Rq, IReadOnlyList<string> Cmd, long Epoch, string Fasta, string Report, string? Error);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_ncbi_pattern_oracle.jsonl.gz";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var gz = new GZipStream(stream, CompressionMode.Decompress);
        using var reader = new StreamReader(gz);
        var rows = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            var j = JsonDocument.Parse(line).RootElement;
            string? Opt(string n) => j.GetProperty(n).ValueKind == JsonValueKind.Null ? null : j.GetProperty(n).GetString();
            rows.Add(new Row(j.GetProperty("id").GetInt32(), j.GetProperty("kind").GetString()!, j.GetProperty("program").GetString()!, j.GetProperty("format").GetString()!, j.GetProperty("params").GetString()!,
                Opt("patfile"), j.GetProperty("pargs").GetString()!, j.GetProperty("rq").GetString()!, [.. j.GetProperty("cmd").EnumerateArray().Select(x => x.GetString()!)], j.GetProperty("epoch").GetInt64(),
                j.GetProperty("fasta").GetString()!, j.GetProperty("report").GetString()!, Opt("error")));
        }

        return rows;
    });

    private static string Option(string tokens, string name)
    {
        var t = tokens.Split(' ', StringSplitOptions.RemoveEmptyEntries);
        int i = Array.IndexOf(t, name);
        return i >= 0 && i + 1 < t.Length ? t[i + 1] : "";
    }

    private static (string Report, string? Refusal) Produce(Row r)
    {
        var request = new EmbossReportRequest(r.Fasta, "", "-rformat " + r.Format + " " + r.Rq, r.Cmd, DateTimeOffset.FromUnixTimeSeconds(r.Epoch).UtcDateTime, "s.fa", "o.rep");
        string spec = r.Params.Split(' ')[0];
        string pname = Option(r.PArgs, "-pname"), pformat = Option(r.PArgs, "-pformat"), pmm = Option(r.PArgs, "-pmismatch");
        var patterns = EmbossPatternList.Read(spec, pname.Length > 0 ? pname : null, pformat.Length > 0 ? pformat : null, r.Program == "fuzznuc" && pmm.Length > 0 ? int.Parse(pmm) : 0, _ => r.PatFile ?? "");
        try
        {
            EmbossProgramReport p = r.Program switch
            {
                "fuzznuc" => EmbossReportPrograms.Fuzznuc(request, patterns, Option(r.Params, "-complement") == "Y"),
                "fuzzpro" => EmbossReportPrograms.Fuzzpro(request, patterns),
                _ => EmbossReportPrograms.Fuzztran(request, patterns, Option(r.Params, "-frame"), 0),
            };
            return (p.Report, null);
        }
        catch (EmbossReportException e)
        {
            return ("", "Error: " + e.Message);
        }
    }

    private static void Check(string kind, string program)
    {
        var rows = Oracle.Value.Where(r => r.Kind == kind && r.Program == program).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(kind == "slow" ? 3 : 150), $"{kind} {program}");
        var failures = new List<string>();
        int withContent = 0;
        foreach (var r in rows)
        {
            var (report, refusal) = Produce(r);
            string? expectedRefusal = r.Error?.Split('\n')[0];
            bool expectedRefused = expectedRefusal is not null && expectedRefusal.StartsWith("Error:", StringComparison.Ordinal) && r.Report.Length == 0;
            if ((refusal is not null) != expectedRefused || (refusal is not null && refusal != expectedRefusal) || report != r.Report)
                failures.Add($"id {r.Id} {r.Program} {r.Format} [{r.Params}] [{r.PArgs}] rq=[{r.Rq}]\n--- fasta\n{r.Fasta}--- patfile\n{r.PatFile}--- expected\n{r.Report}--- actual\n{report}{(refusal is null ? "" : "REFUSED " + refusal)}");
            if (r.Report.Length > 0)
                withContent++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(2)));
        Assert.That(withContent, Is.GreaterThanOrEqualTo(rows.Count / 5), "reports with content");
    }

    [TestCase("ncbi", "fuzznuc")]
    [TestCase("ncbi", "fuzzpro")]
    [TestCase("ncbi", "fuzztran")]
    [TestCase("patfile", "fuzznuc")]
    [TestCase("patfile", "fuzzpro")]
    [TestCase("patfile", "fuzztran")]
    public void ReportFile_EqualsEmboss(string kind, string program) => Check(kind, program);

    [TestCase("fuzzpro")]
    [TestCase("fuzznuc")]
    public void TwoBarredTokens_LeaveTheSequenceUnnamed(string program) => Check("slow", program);

    [Test]
    public void NcbiHeader_WorkedValues()
    {
        // measured with fuzzpro -rformat seqtable -raccshow -rdesshow
        Assert.That(EmbossFastaReader.ParseHeader(">sp|P12345|NAME_HUMAN Some protein OS=Homo"), Is.EqualTo(("NAME_HUMAN", "P12345", "Some protein OS=Homo")));
        Assert.That(EmbossFastaReader.ParseHeader(">gb|AB123456.1|locus1 desc here"), Is.EqualTo(("locus1", "AB123456", "desc here")));
        Assert.That(EmbossFastaReader.ParseHeader(">zz|AB123456.1|nm desc2"), Is.EqualTo(("nm", "AB123456", "desc2")));
        Assert.That(EmbossFastaReader.ParseHeader(">a|b|c|d|e desc"), Is.EqualTo(("e", "", "desc")));
        Assert.That(EmbossFastaReader.ParseHeader(">gi|12345 desc").Name, Is.EqualTo("12345"));
        Assert.That(EmbossFastaReader.ParseHeader(">lcl|chr1:100-200 x").Name, Is.EqualTo("100-200"), "seqSetName keeps the last ':' part");
        Assert.That(EmbossFastaReader.Read(">p1 d\nAC\n>xx|yy some desc\nAC\n>aa|bb\nAC\n>p4\nAC\n").Select(x => x.Name), Is.EqualTo(new[] { "p1", "EMBOSS_001", "EMBOSS_002", "p4" }));
    }

    [Test]
    public void PatternList_WorkedValues()
    {
        // measured with fuzznuc -pattern @file: blank lines count in the numbering and an empty pattern is dropped; a FASTA name is the whole header; <mismatch=n> sets the mismatches of its pattern
        Assert.That(EmbossPatternList.Read("@f", null, null, 0, _ => "ACGT\n\nTTGA\n  GGC  \n").Select(p => p.Name), Is.EqualTo(new[] { "pattern1", "pattern3", "pattern4" }));
        var fasta = EmbossPatternList.Read("@f", null, null, 2, _ => ">pa<mismatch=1>\nCCCGGG\n>pb desc x\nACGT\n>pc\nTT\nT\n");
        Assert.That(fasta, Is.EqualTo(new[] { new EmbossProgramPattern("pa", "CCCGGG", 1), new EmbossProgramPattern("pb desc x", "ACGT", 2), new EmbossProgramPattern("pc", "TTT", 2) }));
        Assert.That(EmbossPatternList.Read("ACGT", "x", null, 1), Is.EqualTo(new[] { new EmbossProgramPattern("x", "ACGT", 1) }));
        Assert.That(EmbossPatternList.Read("@f", "nm", "simple", 0, _ => ">a\nACGT\n").Select(p => p.Name + ":" + p.Pattern), Is.EqualTo(new[] { "nm1:>a", "nm2:ACGT" }));
    }
}
