// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.Motifs;
using Seqeron.Genomics.Emboss.Repeats;
using Seqeron.Genomics.Emboss.FeatureTools;
using System.Globalization;
using System.Text.RegularExpressions;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8 (B09 F176): the feature-report texts <c>-rformat gff</c> and <c>-rformat excel</c> of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem, written by the one canonical <see cref="EmbossReportWriter"/> through
/// <see cref="EmbossReports"/>, against the installed EMBOSS 6.6.0 on 120 random runs per program and format (TestData/Alignment/emboss_report_oracle.tsv, scripts/alignment/generate_emboss_report_oracle.py): the WHOLE
/// text is compared (the <c>#!Date</c> of the gff header is the date of the run and is masked).
/// </summary>
[TestFixture]
public class EmbossReport_Oracle_Tests
{
    private sealed record Row(int Line, string Program, string Format, string Params, string Sequence, string Expected);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.emboss_report_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        var rows = new List<Row>();
        reader.ReadLine();
        int line = 1;
        string? text;
        while ((text = reader.ReadLine()) is not null)
        {
            line++;
            if (text.Length == 0)
                continue;
            var f = text.Split('\t');
            rows.Add(new Row(line, f[1], f[2], f[3], f[4], f[5].Replace("\\n", "\n", StringComparison.Ordinal).Replace("\\t", "\t", StringComparison.Ordinal).Replace("\\\\", "\\", StringComparison.Ordinal)));
        }

        return rows;
    });

    private static string MaskDate(string text) => Regex.Replace(text, "#!Date [0-9-]+", "#!Date X");

    private static string Produce(Row r)
    {
        EmbossReport report;
        string seq = r.Sequence;
        switch (r.Program)
        {
            case "fuzznuc":
                {
                    var p = r.Params.Split('|');
                    int mismatches = int.Parse(p[1], CultureInfo.InvariantCulture);
                    report = EmbossMotifReports.Fuzznuc("s", seq, p[0], EmbossMotifPrograms.FindPatternHits(seq, p[0], mismatches, p[2] == "1", default, sequenceAmbiguity: true), "pattern", mismatches);
                    break;
                }
            case "fuzzpro":
                report = EmbossMotifReports.Fuzzpro("s", seq, r.Params, EmbossMotifPrograms.FindProteinPatternHits(seq, r.Params));
                break;
            case "fuzztran":
                {
                    var p = r.Params.Split('|');
                    report = EmbossMotifReports.Fuzztran("s", p[0], EmbossMotifPrograms.FindTranslatedPatternHits(seq, p[0], p[1]));
                    break;
                }
            case "etandem":
                {
                    var p = r.Params.Split(',').Select(x => int.Parse(x, CultureInfo.InvariantCulture)).ToArray();
                    report = EmbossRepeatReports.Etandem("s", seq, EmbossRepeatPrograms.FindTandemRepeatsEtandem(seq, p[0], p[1], p[2], p[3] == 1, p[4] == 1));
                    break;
                }
            default:
                {
                    var p = r.Params.Split(',').Select(x => int.Parse(x, CultureInfo.InvariantCulture)).ToArray();
                    report = EmbossRepeatReports.Equicktandem("s", seq, EmbossRepeatPrograms.FindTandemRepeatsEquicktandem(seq, p[0], p[1]));
                    break;
                }
        }

        return EmbossReportWriter.Write(report, r.Format == "gff" ? EmbossReportFormat.Gff : EmbossReportFormat.Excel, new DateTime(2026, 11, 20, 0, 0, 0, DateTimeKind.Utc));
    }

    private static void Check(string program, string format)
    {
        var rows = Oracle.Value.Where(r => r.Program == program && r.Format == format).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(100), $"{program} {format}");
        var failures = new List<string>();
        int withHits = 0;
        foreach (var r in rows)
        {
            string expected = MaskDate(r.Expected), got = MaskDate(Produce(r));
            if (got != expected)
                failures.Add($"line {r.Line} {r.Program} {r.Format} [{r.Params}] {r.Sequence}\n--- expected\n{expected}--- actual\n{got}");
            if (r.Expected.Length > 0)
                withHits++;
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows.Count}:\n" + string.Join("\n", failures.Take(2)));
        Assert.That(withHits, Is.GreaterThanOrEqualTo(rows.Count * 8 / 10), "reports with content");
    }

    [Test] public void Fuzznuc_Gff_EqualsEmboss() => Check("fuzznuc", "gff");
    [Test] public void Fuzznuc_Excel_EqualsEmboss() => Check("fuzznuc", "excel");
    [Test] public void Fuzzpro_Gff_EqualsEmboss() => Check("fuzzpro", "gff");
    [Test] public void Fuzzpro_Excel_EqualsEmboss() => Check("fuzzpro", "excel");
    [Test] public void Fuzztran_Gff_EqualsEmboss() => Check("fuzztran", "gff");
    [Test] public void Fuzztran_Excel_EqualsEmboss() => Check("fuzztran", "excel");
    [Test] public void Etandem_Gff_EqualsEmboss() => Check("etandem", "gff");
    [Test] public void Etandem_Excel_EqualsEmboss() => Check("etandem", "excel");
    [Test] public void Equicktandem_Gff_EqualsEmboss() => Check("equicktandem", "gff");
    [Test] public void Equicktandem_Excel_EqualsEmboss() => Check("equicktandem", "excel");

    [Test]
    public void FuzznucGff_ExampleOfTheAudit()
    {
        // fuzznuc a -pattern ACGTA -rformat2 gff (audit A5-8): a<TAB>fuzznuc<TAB>nucleotide_motif<TAB>1<TAB>5<TAB>5<TAB>+<TAB>.<TAB>ID=a.1;note=*pat pattern:ACGTA below the ##gff-version 3 header
        var report = EmbossMotifReports.Fuzznuc("a", "ACGTAGGG", "ACGTA", EmbossMotifPrograms.FindPatternHits("ACGTAGGG", "ACGTA"));
        string text = EmbossReportWriter.Write(report, EmbossReportFormat.Gff, new DateTime(2026, 11, 20, 0, 0, 0, DateTimeKind.Utc));
        Assert.That(text, Is.EqualTo("##gff-version 3\n##sequence-region a 1 8\n#!Date 2026-11-20\n#!Type DNA\n#!Source-version EMBOSS 6.6.0.0\na\tfuzznuc\tnucleotide_motif\t1\t5\t5\t+\t.\tID=a.1;note=*pat pattern:ACGTA\n"));
        Assert.That(EmbossReportWriter.Write(report, EmbossReportFormat.Excel, default), Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\tPattern\tMismatch\na\t1\t5\t5\t+\tpattern:ACGTA\t.\n"));
        Assert.That(EmbossReportWriter.Write(EmbossMotifReports.Fuzznuc("a", "AAAA", "CCC", []), EmbossReportFormat.Excel, default), Is.Empty);        // fuzznuc writes nothing for no hit
        Assert.That(EmbossReportWriter.Write(EmbossRepeatReports.Etandem("a", "ACGT", []), EmbossReportFormat.Excel, default), Is.EqualTo("SeqName\tStart\tEnd\tScore\tStrand\tSize\tCount\tIdentity\tConsensus\n"));    // etandem writes the header
    }
}
