// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Repeats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.Repeats;

using Seqeron.Genomics.Analysis;

/// <summary>
/// The feature tables of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem (audit A5-8 / A5-8b, B09 F176 / F184): the program's hits as <see cref="EmbossReport"/> (what each program stores in its feature table: source, type, coordinates,
/// score, tags), written in every report format by the one canonical <see cref="EmbossReportWriter"/>. Verified against the installed programs
/// (TestData/Alignment/emboss_report_formats_oracle.jsonl.gz, scripts/alignment/generate_emboss_report_formats_oracle.py).
/// </summary>
public static class EmbossRepeatReports
{
    // (the pattern groups of the multi-pattern overloads are EmbossPatternHits<T>)
    private static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static EmbossReportSequence PlainSequence(string name, string text, bool nucleotide)
    {
        string clean = nucleotide ? EmbossSequenceReader.Clean(text) : text;
        return new EmbossReportSequence(name, name, "", "", clean, 1, clean.Length, 0, false, nucleotide);
    }

    private static EmbossFeatureTag Note(string text) => new("note", text);

    /// <summary>The table of <c>etandem</c> for a sequence given by name and text.</summary>
    public static EmbossReport Etandem(string sequenceName, string sequence, IReadOnlyList<EmbossTandemRepeat> hits)
    {
        ArgumentNullException.ThrowIfNull(sequenceName);
        ArgumentNullException.ThrowIfNull(sequence);
        return Etandem(PlainSequence(sequenceName, sequence, true), hits);
    }

    /// <summary>The table of <c>etandem</c> (the program's report order: highest score first): size, count, identity (one decimal) and the consensus.</summary>
    public static EmbossReport Etandem(EmbossReportSequence sequence, IReadOnlyList<EmbossTandemRepeat> hits)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(hits);
        var features = hits.Select(h => new EmbossFeature("SO:0000657", "etandem", h.Start, h.End, h.Score, '+', 0,
            [new("rpt_type", "TANDEM"), Note("*rpt_size " + Int(h.RepeatSize)), Note("*rpt_count " + Int(h.RepeatCount)), Note("*identity " + CFormat.Fixed(h.Identity, 1)), Note("*consensus " + h.Consensus)])).ToList();
        return new EmbossReport(sequence, features, false, "N", TableLength: sequence.Length, Tags: EmbossReportTagLists.TagList("etandem"), WriteWhenEmpty: true);
    }

    /// <summary>The table of <c>equicktandem</c> for a sequence given by name and text.</summary>
    public static EmbossReport Equicktandem(string sequenceName, string sequence, IReadOnlyList<EmbossTandemRepeat> hits)
    {
        ArgumentNullException.ThrowIfNull(sequenceName);
        ArgumentNullException.ThrowIfNull(sequence);
        return Equicktandem(PlainSequence(sequenceName, sequence, true), hits);
    }

    /// <summary>The table of <c>equicktandem</c>: size and count.</summary>
    public static EmbossReport Equicktandem(EmbossReportSequence sequence, IReadOnlyList<EmbossTandemRepeat> hits)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(hits);
        var features = hits.Select(h => new EmbossFeature("SO:0000657", "equicktandem", h.Start, h.End, h.Score, '+', 0,
            [new("rpt_type", "TANDEM"), Note("*rpt_size " + Int(h.RepeatSize)), Note("*rpt_count " + Int(h.RepeatCount))])).ToList();
        return new EmbossReport(sequence, features, false, "N", TableLength: sequence.Length, Tags: EmbossReportTagLists.TagList("equicktandem"), WriteWhenEmpty: true);
    }
}
