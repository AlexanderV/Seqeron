using System.Text;
using System.Text.Json.Nodes;
using System.Text.RegularExpressions;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit R14-9, B09 F302): the extended CIGAR, the MD tag and the cs string of an <see cref="AlignmentResult"/> against the tools that write them:
/// parasail 1.3.4 (667 traceback alignments: <c>nw</c>, <c>sg</c>, <c>sg_qx</c>, <c>sg_dx</c>, <c>sg_qb</c>, <c>sg_db</c>, <c>sg_qe</c>, <c>sg_de</c>, ..., DNA and BLOSUM62 protein; the cigar the library decodes from the same traceback) and
/// minimap2 2.26 (120 alignments; PAF <c>cg</c> with <c>--eqx</c>, <c>MD</c>, <c>cs=short</c> and <c>cs=long</c>). <c>scripts/alignment/generate_alignment_result_cigar_oracle.py</c>.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class AlignmentResultCigar_Oracle_Tests
{
    private static readonly Lazy<List<JsonObject>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.alignment_result_cigar_oracle.jsonl";
        using var stream = TestDataResource.Open(name) ?? throw new InvalidOperationException(name);
        using var reader = new StreamReader(stream);
        var list = new List<JsonObject>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length > 0)
                list.Add((JsonObject)JsonNode.Parse(line)!);
        }

        return list;
    });

    private static IEnumerable<TestCaseData> Of(string kind) => Oracle.Value.Where(o => o["kind"]!.GetValue<string>() == kind).Select(o => new TestCaseData(o.ToJsonString()).SetName($"{kind}_{o["id"]!.GetValue<string>()}"));

    private static IEnumerable<TestCaseData> ParasailCases() => Of("parasail");

    private static IEnumerable<TestCaseData> Minimap2Cases() => Of("minimap2");

    [TestCaseSource(nameof(ParasailCases))]
    public void ExtendedCigar_EqualsParasail(string json)
    {
        var o = (JsonObject)JsonNode.Parse(json)!;
        // parasail aligns query against ref: the reference row is the target (row 1), a residue of the query against a gap is an insertion
        var r = new AlignmentResult(o["row1"]!.GetValue<string>(), o["row2"]!.GetValue<string>(), o["score"]!.GetValue<int>(), AlignmentType.Global, 0, 0, 0, 0);
        Assert.That(r.ToCigar(extended: true), Is.EqualTo(o["cigar"]!.GetValue<string>()));
        // the standard CIGAR is the same with the M runs merged
        string merged = Regex.Replace(o["cigar"]!.GetValue<string>(), "([0-9]+)([=X])", "$1M");
        Assert.That(r.ToCigar(), Is.EqualTo(MergeM(merged)));
    }

    private static string MergeM(string cigar)
    {
        var sb = new StringBuilder();
        int run = 0;
        foreach (Match m in Regex.Matches(cigar, "([0-9]+)([MIDNSHP=X])"))
        {
            int n = int.Parse(m.Groups[1].Value);
            if (m.Groups[2].Value == "M")
            {
                run += n;
                continue;
            }

            if (run > 0)
                sb.Append(run).Append('M');
            run = 0;
            sb.Append(n).Append(m.Groups[2].Value);
        }

        if (run > 0)
            sb.Append(run).Append('M');
        return sb.ToString();
    }

    /// <summary>The rows of a PAF alignment: the CIGAR (<c>=</c>, <c>X</c> or <c>M</c>, <c>I</c>, <c>D</c>) applied to the target from <paramref name="tstart"/> and the query from <paramref name="qstart"/>.</summary>
    private static AlignmentResult FromCigar(string target, string query, int tstart, int qstart, string cigar)
    {
        var r1 = new StringBuilder();
        var r2 = new StringBuilder();
        int t = tstart, q = qstart;
        foreach (Match m in Regex.Matches(cigar, "([0-9]+)([MIDX=])"))
        {
            int n = int.Parse(m.Groups[1].Value);
            switch (m.Groups[2].Value[0])
            {
                case 'I':
                    r1.Append('-', n);
                    r2.Append(query, q, n);
                    q += n;
                    break;
                case 'D':
                    r1.Append(target, t, n);
                    r2.Append('-', n);
                    t += n;
                    break;
                default:
                    r1.Append(target, t, n);
                    r2.Append(query, q, n);
                    t += n;
                    q += n;
                    break;
            }
        }

        return new AlignmentResult(r1.ToString(), r2.ToString(), 0, AlignmentType.Local, tstart, qstart, t - 1, q - 1);
    }

    [TestCaseSource(nameof(Minimap2Cases))]
    public void CigarMdAndCs_EqualMinimap2(string json)
    {
        var o = (JsonObject)JsonNode.Parse(json)!;
        string target = o["ref"]!.GetValue<string>(), query = o["query"]!.GetValue<string>();
        string cg = o["cg"]!.GetValue<string>();
        var r = FromCigar(target, query, o["tstart"]!.GetValue<int>(), o["qstart"]!.GetValue<int>(), cg);
        Assert.That(r.ToCigar(extended: true), Is.EqualTo(cg));
        Assert.That(r.ToCigar(), Is.EqualTo(MergeM(Regex.Replace(cg, "([0-9]+)[=X]", "$1M"))));
        Assert.That(r.ToMd(), Is.EqualTo(o["md"]!.GetValue<string>()));
        Assert.That(r.ToCs(), Is.EqualTo(o["cs_short"]!.GetValue<string>()));
        Assert.That(r.ToCs(longForm: true), Is.EqualTo(o["cs_long"]!.GetValue<string>()));
        // the clips of the SAM line: the query residues outside the PAF interval
        string clipped = r.ToCigar(target.Length, query.Length, extended: true);
        int qs = o["qstart"]!.GetValue<int>(), qe = o["qend"]!.GetValue<int>();
        Assert.That(clipped, Is.EqualTo((qs > 0 ? qs + "S" : "") + cg + (qe < query.Length ? (query.Length - qe) + "S" : "")));
    }

    [Test]
    public void Oracle_HasTheAnnouncedInventory()
    {
        Assert.That(Oracle.Value.Count(o => o["kind"]!.GetValue<string>() == "parasail"), Is.EqualTo(667));
        Assert.That(Oracle.Value.Count(o => o["kind"]!.GetValue<string>() == "minimap2"), Is.EqualTo(120));
        Assert.That(Oracle.Value.Where(o => o["kind"]!.GetValue<string>() == "parasail").Select(o => o["fn"]!.GetValue<string>()).Distinct().Count(), Is.EqualTo(10));
    }

    [Test]
    public void Cs_OfOwnAlignment_UsesTheMinimap2Notation()
    {
        var r = new AlignmentResult("ACGTAC-GT", "ACCTACAGT", 0, AlignmentType.Global, 0, 0, 8, 8);
        Assert.That(r.ToCs(), Is.EqualTo(":2*gc:3+a:2"));
        Assert.That(r.ToCs(longForm: true), Is.EqualTo("=AC*gc=TAC+a=GT"));
        Assert.That(r.ToMd(), Is.EqualTo("2G5"));
        Assert.That(r.ToCigar(extended: true), Is.EqualTo("2=1X3=1I2="));
    }
}
