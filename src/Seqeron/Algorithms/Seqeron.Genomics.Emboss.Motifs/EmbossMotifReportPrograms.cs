// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Motifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.Motifs;

using Seqeron.Genomics.Analysis;

/// <summary>
/// The report files of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem as the programs write them (<c>fuzznuc.c</c> ... <c>equicktandem.c</c>, EMBOSS-6.6.0): the search of the ported program, the feature tables of
/// <see cref="EmbossReports"/> and the one writer <see cref="EmbossReportWriter"/>, with every <c>-rformat</c>, the report options (<c>-raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq</c>), several sequences per file
/// and the extra <c>-origfile</c> of the tandem programs. Compared with the installed programs on TestData/Alignment/emboss_report_formats_oracle.jsonl.gz (audit A5-8b, B09 F184).
/// </summary>
public static class EmbossMotifReportPrograms
{
    private static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static string Bool(bool b) => b ? "Yes" : "No";

    // ajPatlistSeqDoc
    private static string PatternDoc(IReadOnlyList<EmbossProgramPattern> patterns)
    {
        var sb = new StringBuilder(string.Create(CultureInfo.InvariantCulture, $"{"Pattern_name",-12} {"Mismatch",8} {"Pattern"}\n"));
        foreach (var p in patterns)
            sb.Append(string.Create(CultureInfo.InvariantCulture, $"{p.Name,-12} {p.Mismatch,8} {p.Pattern}\n"));
        return sb.ToString();
    }

    /// <summary>The report file of <c>fuzznuc</c> (<c>fuzznuc.c</c>): header text = pattern documentation and <c>Complement</c>; default format <c>seqtable</c>.</summary>
    /// <param name="r">The request.</param>
    /// <param name="pattern">The pattern (<c>-pattern</c>).</param>
    /// <param name="mismatch">Mismatches (<c>-pmismatch</c>).</param>
    /// <param name="complement">Also search the complement strand (<c>-complement</c>).</param>
    public static EmbossProgramReport Fuzznuc(EmbossReportRequest r, string pattern, int mismatch = 0, bool complement = false, IReadOnlyList<PatternHit>? hits = null)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        return Fuzznuc(r, [new EmbossProgramPattern(r.PatternName, pattern, mismatch)], complement, hits is null ? null : [hits]);
    }

    /// <summary>
    /// The report file of <c>fuzznuc</c> for a pattern list (<c>-pattern @file</c>, <see cref="EmbossPatternList.Read"/>; A5-8c): the header table names every pattern with its mismatches, the features are the hits of
    /// every pattern (forward strand first, then the complement strand).
    /// </summary>
    public static EmbossProgramReport Fuzznuc(EmbossReportRequest r, IReadOnlyList<EmbossProgramPattern> patterns, bool complement = false, IReadOnlyList<IReadOnlyList<PatternHit>>? hits = null)
    {
        ArgumentNullException.ThrowIfNull(patterns);
        var s = EmbossReportPrograms.Prepare("fuzznuc", r, "seqtable", PatternDoc(patterns) + "\nComplement: " + Bool(complement) + "\n", true, false);
        var records = EmbossFastaReader.Read(r.Fasta);
        var tables = new List<EmbossReport>();
        var index = new List<int>();
        var lengths = new List<long>();
        for (int i = 0; i < records.Count; i++)
        {
            var rec = records[i];
            var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Dna);
            lengths.Add(view.Residues.Length);
            var groups = patterns.Select((p, k) => new EmbossPatternHits<PatternHit>(p.Name, p.Pattern, p.Mismatch,
                i == 0 && hits is not null ? hits[k] : EmbossMotifPrograms.FindPatternHits(rec.Text, s.Inputs[0], p.Pattern, p.Mismatch, complement, r.CancellationToken, programOrder: true))).ToList();
            var table = EmbossMotifReports.Fuzznuc(EmbossReportPrograms.ReportSequence(r, rec, view, s.Inputs[0], view.IsNucleotide), groups);
            if (table.Features.Count > 0)
            {
                tables.Add(table);
                index.Add(i);
            }
        }

        return new EmbossProgramReport(EmbossReportPrograms.Finish(s, tables, index, lengths, records.Count), null);
    }

    /// <summary>The report file of <c>fuzzpro</c>; default format <c>seqtable</c>.</summary>
    public static EmbossProgramReport Fuzzpro(EmbossReportRequest r, string pattern, IReadOnlyList<ProteinPatternHit>? hits = null, int mismatch = 0)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        return Fuzzpro(r, [new EmbossProgramPattern(r.PatternName, pattern, mismatch)], hits is null ? null : [hits]);
    }

    /// <summary>The report file of <c>fuzzpro</c> for a pattern list (<c>-pattern @file</c>; A5-8c).</summary>
    public static EmbossProgramReport Fuzzpro(EmbossReportRequest r, IReadOnlyList<EmbossProgramPattern> patterns, IReadOnlyList<IReadOnlyList<ProteinPatternHit>>? hits = null)
    {
        ArgumentNullException.ThrowIfNull(patterns);
        var s = EmbossReportPrograms.Prepare("fuzzpro", r, "seqtable", PatternDoc(patterns), true, false);
        var records = EmbossFastaReader.Read(r.Fasta);
        var tables = new List<EmbossReport>();
        var index = new List<int>();
        var lengths = new List<long>();
        for (int i = 0; i < records.Count; i++)
        {
            var rec = records[i];
            var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Protein);
            lengths.Add(view.Residues.Length);
            var groups = patterns.Select((p, k) => new EmbossPatternHits<ProteinPatternHit>(p.Name, p.Pattern, p.Mismatch,
                i == 0 && hits is not null ? hits[k] : EmbossMotifPrograms.FindProteinPatternHits(rec.Text, s.Inputs[0], p.Pattern, r.CancellationToken, p.Mismatch))).ToList();
            var table = EmbossMotifReports.Fuzzpro(EmbossReportPrograms.ReportSequence(r, rec, view, s.Inputs[0], false), groups, i == 0);
            if (table.Features.Count > 0)
            {
                tables.Add(table);
                index.Add(i);
            }
        }

        return new EmbossProgramReport(EmbossReportPrograms.Finish(s, tables, index, lengths, records.Count), null);
    }

    /// <summary>The report file of <c>fuzztran</c>; default format <c>table</c>; the header names the genetic code table and the frames as typed.</summary>
    /// <param name="r">The request.</param>
    /// <param name="pattern">The protein pattern.</param>
    /// <param name="frame">The frames (<c>-frame</c>: 1 2 3 F -1 -2 -3 R 6).</param>
    /// <param name="table">The genetic code (<c>-table</c>).</param>
    public static EmbossProgramReport Fuzztran(EmbossReportRequest r, string pattern, string frame = "1", int table = 0, IReadOnlyList<TranslatedPatternHit>? hits = null, int mismatch = 0)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        return Fuzztran(r, [new EmbossProgramPattern(r.PatternName, pattern, mismatch)], frame, table, hits is null ? null : [hits]);
    }

    /// <summary>The report file of <c>fuzztran</c> for a pattern list (<c>-pattern @file</c>; A5-8c).</summary>
    public static EmbossProgramReport Fuzztran(EmbossReportRequest r, IReadOnlyList<EmbossProgramPattern> patterns, string frame = "1", int table = 0, IReadOnlyList<IReadOnlyList<TranslatedPatternHit>>? hits = null)
    {
        ArgumentNullException.ThrowIfNull(patterns);
        ArgumentNullException.ThrowIfNull(frame);
        var s = EmbossReportPrograms.Prepare("fuzztran", r, "table", PatternDoc(patterns) + "TransTable: " + Int(table) + "\nFrames: " + frame + "\n", true, false);
        var records = EmbossFastaReader.Read(r.Fasta);
        var tables = new List<EmbossReport>();
        var index = new List<int>();
        var lengths = new List<long>();
        for (int i = 0; i < records.Count; i++)
        {
            var rec = records[i];
            var view = EmbossInputSequence.Read(rec.Text, s.Inputs[0], EmbossSequenceType.Dna);
            lengths.Add(view.Residues.Length);
            var groups = patterns.Select((p, k) => new EmbossPatternHits<TranslatedPatternHit>(p.Name, p.Pattern, p.Mismatch,
                i == 0 && hits is not null ? hits[k] : EmbossMotifPrograms.FindTranslatedPatternHits(rec.Text, s.Inputs[0], p.Pattern, frame, table, r.CancellationToken, p.Mismatch))).ToList();
            var t = EmbossMotifReports.Fuzztran(EmbossReportPrograms.ReportSequence(r, rec, view, s.Inputs[0], view.IsNucleotide), groups);
            if (t.Features.Count > 0)
            {
                tables.Add(t);
                index.Add(i);
            }
        }

        return new EmbossProgramReport(EmbossReportPrograms.Finish(s, tables, index, lengths, records.Count), null);
    }
}
