// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Repeats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.Repeats;

using Seqeron.Genomics.Analysis;

/// <summary>
/// The report files of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem as the programs write them (<c>fuzznuc.c</c> ... <c>equicktandem.c</c>, EMBOSS-6.6.0): the search of the ported program, the feature tables of
/// <see cref="EmbossReports"/> and the one writer <see cref="EmbossReportWriter"/>, with every <c>-rformat</c>, the report options (<c>-raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq</c>), several sequences per file
/// and the extra <c>-origfile</c> of the tandem programs. Compared with the installed programs on TestData/Alignment/emboss_report_formats_oracle.jsonl.gz (audit A5-8b, B09 F184).
/// </summary>
public static class EmbossRepeatReportPrograms
{
    private static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static string Bool(bool b) => b ? "Yes" : "No";

    private static EmbossFastaRecord FirstRecord(EmbossReportRequest r)
    {
        var records = EmbossFastaReader.Read(r.Fasta);
        if (records.Count == 0)
            throw new ArgumentException("The input has no FASTA record.", nameof(r));
        return records[0];
    }

    /// <summary>The report file of <c>etandem</c> (one sequence; also the <c>-origfile</c> text); default format <c>table</c>.</summary>
    public static EmbossProgramReport Etandem(EmbossReportRequest r, int minRepeat, int maxRepeat, int threshold = 20, bool mismatch = false, bool uniform = false, IReadOnlyList<EmbossTandemRepeat>? hits = null)
    {
        var s = EmbossReportPrograms.Prepare("etandem", r, "table", $"Threshold: {Int(threshold)}\nMinrepeat: {Int(minRepeat)}\nMaxrepeat: {Int(maxRepeat)}\nMismatch: {Bool(mismatch)}\nUniform: {Bool(uniform)}\n", false, true);
        var rec = FirstRecord(r);
        var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Dna);
        var found = hits ?? EmbossRepeatPrograms.FindTandemRepeatsEtandem(rec.Text, s.Inputs[0], minRepeat, maxRepeat, threshold, mismatch, uniform, r.CancellationToken);
        var table = EmbossRepeatReports.Etandem(EmbossReportPrograms.ReportSequence(r, rec, view, s.Inputs[0], view.IsNucleotide), found);
        var orig = new StringBuilder();
        foreach (var h in table.Features)
        {
            // etandem_basicReport: "%6d %10d %10d %2d %3d %5.1f " and the consensus
            string size = h.Tags[1].Value[("*rpt_size ".Length)..], count = h.Tags[2].Value[("*rpt_count ".Length)..], ident = h.Tags[3].Value[("*identity ".Length)..], cons = h.Tags[4].Value[("*consensus ".Length)..];
            orig.Append(Int((long)h.Score).PadLeft(6)).Append(' ').Append(Int(h.Start).PadLeft(10)).Append(' ').Append(Int(h.End).PadLeft(10)).Append(' ').Append(size.PadLeft(2)).Append(' ').Append(count.PadLeft(3)).Append(' ').Append(ident.PadLeft(5)).Append(' ').Append(cons).Append('\n');
        }

        return Written(s, table, orig.ToString());
    }

    // the report of one table; the original file of the tandem programs is written before the report is refused
    private static EmbossProgramReport Written(EmbossReportPrograms.Setup s, EmbossReport table, string orig)
    {
        try
        {
            return new EmbossProgramReport(EmbossReportWriter.WriteFile(s.Options, [table]), orig);
        }
        catch (EmbossReportException e)
        {
            throw new EmbossReportException(e.Message) { OriginalFile = orig };
        }
    }

    /// <summary>The report file of <c>equicktandem</c> (one sequence; also the <c>-origfile</c> text); default format <c>table</c>.</summary>
    public static EmbossProgramReport Equicktandem(EmbossReportRequest r, int maxRepeat, int threshold = 20, IReadOnlyList<EmbossTandemRepeat>? hits = null)
    {
        var s = EmbossReportPrograms.Prepare("equicktandem", r, "table", $"Threshold: {Int(threshold)}\nMaxrepeat: {Int(maxRepeat)}\n", false, true);
        var rec = FirstRecord(r);
        var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Dna);
        var found = hits ?? EmbossRepeatPrograms.FindTandemRepeatsEquicktandem(rec.Text, s.Inputs[0], maxRepeat, threshold, r.CancellationToken);
        var table = EmbossRepeatReports.Equicktandem(EmbossReportPrograms.ReportSequence(r, rec, view, s.Inputs[0], view.IsNucleotide), found);
        var orig = new StringBuilder();
        foreach (var h in table.Features)
            orig.Append(Int((long)h.Score).PadLeft(6)).Append(' ').Append(Int(h.Start).PadLeft(10)).Append(' ').Append(Int(h.End).PadLeft(10)).Append(' ').Append(h.Tags[1].Value[("*rpt_size ".Length)..].PadLeft(2)).Append(' ').Append(h.Tags[2].Value[("*rpt_count ".Length)..].PadLeft(3)).Append('\n');
        return Written(s, table, orig.ToString());
    }
}
