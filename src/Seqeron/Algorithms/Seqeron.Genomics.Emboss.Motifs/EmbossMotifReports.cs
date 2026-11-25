// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.Motifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.Motifs;

using Seqeron.Genomics.Analysis;

/// <summary>
/// The feature tables of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem (audit A5-8 / A5-8b, B09 F176 / F184): the program's hits as <see cref="EmbossReport"/> (what each program stores in its feature table: source, type, coordinates,
/// score, tags), written in every report format by the one canonical <see cref="EmbossReportWriter"/>. Verified against the installed programs
/// (TestData/Alignment/emboss_report_formats_oracle.jsonl.gz, scripts/alignment/generate_emboss_report_formats_oracle.py).
/// </summary>
public static class EmbossMotifReports
{
    // (the pattern groups of the multi-pattern overloads are EmbossPatternHits<T>)
    private static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    private static EmbossReportSequence PlainSequence(string name, string text, bool nucleotide)
    {
        string clean = nucleotide ? EmbossSequenceReader.Clean(text) : text;
        return new EmbossReportSequence(name, name, "", "", clean, 1, clean.Length, 0, false, nucleotide);
    }

    // the strand as the feature stores it: forward is NUL; ajFeatReverse of a reversed sequence turns NUL into '-' and '-' into '+'
    private static char Strand(bool outputReverse, bool reversedSequence)
    {
        if (reversedSequence)
            return outputReverse ? '-' : '+';
        return outputReverse ? '-' : '\0';
    }

    private static EmbossFeatureTag Note(string text) => new("note", text);

    // a pattern with a range element X(n,m): the program lists the hits of one frame / strand with the START DESCENDING (the hits of one start by ascending length), not in position order (measured on fuzznuc, fuzzpro and fuzztran 6.6.0)
    private static readonly System.Text.RegularExpressions.Regex RangeElement = new(@"\(\s*\d+\s*,\s*\d+\s*\)", System.Text.RegularExpressions.RegexOptions.Compiled | System.Text.RegularExpressions.RegexOptions.CultureInvariant);

    // the hits of the complement search are listed in the order of the reverse strand: END descending, then start ascending (the longer hit first); with mismatches by end ascending, then start (brute-force search).
    // `key` gives the start and length in the coordinates of the searched text (mirrored for a reversed sequence: the program searches the reverse complement)
    private static IEnumerable<PatternHit> ReverseOrder(IEnumerable<PatternHit> hits, string pattern, int maxMismatches, Func<PatternHit, (int Start, int Length)> key)
    {
        if (maxMismatches > 0)
            return RangeElement.IsMatch(pattern) ? hits : hits.OrderBy(h => key(h).Start + key(h).Length).ThenBy(h => key(h).Start);        // a range pattern: FindPatternHits lists the pass in the program's (reversed push) order
        return RangeElement.IsMatch(pattern) ? hits.OrderByDescending(h => key(h).Start + key(h).Length).ThenBy(h => key(h).Start) : hits;
    }

    private static IEnumerable<T> ReportOrder<T>(IEnumerable<T> hits, string pattern, Func<T, int> start, Func<T, int> length)
        => RangeElement.IsMatch(pattern) ? hits.OrderByDescending(start).ThenBy(length) : hits;

    /// <summary>The report of <c>fuzznuc</c> for a sequence given by name and text (<see cref="Fuzznuc(EmbossReportSequence, string, IReadOnlyList{PatternHit}, string, int)"/>).</summary>
    public static EmbossReport Fuzznuc(string sequenceName, string sequence, string pattern, IReadOnlyList<PatternHit> hits, string patternName = "pattern", int maxMismatches = 0)
    {
        ArgumentNullException.ThrowIfNull(sequenceName);
        ArgumentNullException.ThrowIfNull(sequence);
        return Fuzznuc(PlainSequence(sequenceName, sequence, true), pattern, hits, patternName, maxMismatches);
    }

    /// <summary>The table of <c>fuzznuc</c>: forward hits in position order (a range pattern searched without mismatches: start descending), then the hits of the complement strand; score = hit length minus mismatches; the pattern is noted as <c><paramref name="patternName"/>:<paramref name="pattern"/></c>.</summary>
    public static EmbossReport Fuzznuc(EmbossReportSequence sequence, string pattern, IReadOnlyList<PatternHit> hits, string patternName = "pattern", int maxMismatches = 0)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(hits);
        ArgumentNullException.ThrowIfNull(patternName);
        return Fuzznuc(sequence, [new EmbossPatternHits<PatternHit>(patternName, pattern, maxMismatches, hits)]);
    }

    /// <summary>
    /// The table of <c>fuzznuc</c> for a pattern list (<c>-pattern @file</c>, A5-8c): <c>embPatlistSeqSearchAll</c> searches the patterns one after the other, so the table holds the forward hits of every pattern (in the order of
    /// the list), then the hits of the complement strand of every pattern.
    /// </summary>
    public static EmbossReport Fuzznuc(EmbossReportSequence sequence, IReadOnlyList<EmbossPatternHits<PatternHit>> patterns)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(patterns);
        // a reversed sequence is searched as its reverse complement: the hits of the first pass come out on the minus strand and their positions are mirrored back
        bool rev = sequence.Reversed;
        int n = sequence.Length;
        (int Start, int Length) Key(PatternHit h) => rev ? (n - h.Position - h.Length, h.Length) : (h.Position, h.Length);
        var ordered = new List<(EmbossPatternHits<PatternHit> P, PatternHit H)>();
        foreach (var p in patterns)
        {
            string orderPattern = p.Mismatches > 0 ? "" : p.Pattern;          // a search with mismatches (patOUBrute) lists a range pattern's hits in position order
            ordered.AddRange(ReportOrder(p.Hits.Where(h => h.IsReverseStrand == rev), orderPattern, h => Key(h).Start, h => h.Length).Select(h => (p, h)));
        }

        foreach (var p in patterns)
            ordered.AddRange(ReverseOrder(p.Hits.Where(h => h.IsReverseStrand != rev), p.Pattern, p.Mismatches, Key).Select(h => (p, h)));
        var features = ordered.Select(x =>
        {
            var (p, h) = x;
            var notes = new List<EmbossFeatureTag> { Note("*pat " + p.Name + ":" + p.Pattern) };
            if (h.Mismatches > 0)
                notes.Add(Note("*mismatch " + Int(h.Mismatches)));
            return new EmbossFeature("SO:0000714", "fuzznuc", h.Position + 1, h.Position + h.Length, h.Length - h.Mismatches, Strand(h.IsReverseStrand, rev), 0, notes);
        }).ToList();
        return new EmbossReport(sequence, features, false, "N", TableLength: sequence.Length, Tags: EmbossReportTagLists.TagList("fuzznuc"));
    }

    /// <summary>The table of <c>fuzzpro</c> for a sequence given by name and text.</summary>
    public static EmbossReport Fuzzpro(string sequenceName, string sequence, string pattern, IReadOnlyList<ProteinPatternHit> hits, string patternName = "pattern")
    {
        ArgumentNullException.ThrowIfNull(sequenceName);
        ArgumentNullException.ThrowIfNull(sequence);
        return Fuzzpro(PlainSequence(sequenceName, sequence, false), pattern, hits, patternName);
    }

    /// <summary>The table of <c>fuzzpro</c>: the hits in position order, score = hit length; its length is not set (the GFF header names the largest hit end). <paramref name="first"/>: the first table of the run (later tables are reset tables, the GFF header then reads <c>DNA</c>).</summary>
    public static EmbossReport Fuzzpro(EmbossReportSequence sequence, string pattern, IReadOnlyList<ProteinPatternHit> hits, string patternName = "pattern", bool first = true)
    {
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(hits);
        ArgumentNullException.ThrowIfNull(patternName);
        return Fuzzpro(sequence, [new EmbossPatternHits<ProteinPatternHit>(patternName, pattern, 0, hits)], first);
    }

    /// <summary>The table of <c>fuzzpro</c> for a pattern list (<c>-pattern @file</c>, A5-8c): the hits of every pattern in the order of the list, each pattern's in position order.</summary>
    public static EmbossReport Fuzzpro(EmbossReportSequence sequence, IReadOnlyList<EmbossPatternHits<ProteinPatternHit>> patterns, bool first = true)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(patterns);
        // a search with mismatches (patOUBrute) lists a range pattern's hits in position order; the score is the length minus the mismatches and a hit with mismatches carries a *mismatch note (A5-8c2)
        var features = patterns.SelectMany(p => ReportOrder(p.Hits, p.Mismatches > 0 ? "" : p.Pattern, h => h.Position, h => h.Length).Select(h => new EmbossFeature("SO:0001067", "fuzzpro", h.Position + 1, h.Position + h.Length, h.Length - h.Mismatches, '\0', 0,
            h.Mismatches > 0 ? [Note("*pat " + p.Name + ":" + p.Pattern), Note("*mismatch " + Int(h.Mismatches))] : [Note("*pat " + p.Name + ":" + p.Pattern)]))).ToList();
        return new EmbossReport(sequence, features, true, first ? "P" : "p", TableLength: EmbossReportWriter.UnsignedRegionEnd(features), Tags: EmbossReportTagLists.TagList("fuzzpro"));
    }

    /// <summary>The table of <c>fuzztran</c> for a sequence given by name (the report sequence is the nucleotide text of the plain reader).</summary>
    public static EmbossReport Fuzztran(string sequenceName, string pattern, IReadOnlyList<TranslatedPatternHit> hits, string patternName = "pattern", string sequence = "")
    {
        ArgumentNullException.ThrowIfNull(sequenceName);
        ArgumentNullException.ThrowIfNull(sequence);
        return Fuzztran(PlainSequence(sequenceName, sequence, true), pattern, hits, patternName);
    }

    /// <summary>
    /// The table of <c>fuzztran</c>: the hits frame by frame (1, 2, 3, -1, -2, -3), the program's own coordinates (<see cref="TranslatedPatternHit.EmbossStart"/> / <see cref="TranslatedPatternHit.EmbossEnd"/>: negative for the reverse frames),
    /// score = residues of the hit, the frame number without its sign; the strand is stored as <c>+</c> / <c>-</c>; its length is not set.
    /// </summary>
    /// <exception cref="ArgumentException">A hit has no <see cref="TranslatedPatternHit.EmbossStart"/> / <see cref="TranslatedPatternHit.EmbossEnd"/> (a hit made by hand: the hits of <c>FindTranslatedPatternHits</c> carry the program's coordinates) or no translation (audit A6-4).</exception>
    public static EmbossReport Fuzztran(EmbossReportSequence sequence, string pattern, IReadOnlyList<TranslatedPatternHit> hits, string patternName = "pattern")
    {
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(hits);
        ArgumentNullException.ThrowIfNull(patternName);
        for (int i = 0; i < hits.Count; i++)
        {
            if (hits[i].EmbossStart is null || hits[i].EmbossEnd is null)
                throw new ArgumentException($"Hit {i} has no EmbossStart / EmbossEnd: the report prints the program's own coordinates (the hits of EmbossMotifPrograms.FindTranslatedPatternHits carry them).", nameof(hits));
            if (hits[i].Translation is null)
                throw new ArgumentException($"Hit {i} has no Translation.", nameof(hits));
        }

        return Fuzztran(sequence, [new EmbossPatternHits<TranslatedPatternHit>(patternName, pattern, 0, hits)]);
    }

    /// <summary>The table of <c>fuzztran</c> for a pattern list (<c>-pattern @file</c>, A5-8c): frame by frame (1, 2, 3, -1, -2, -3), in each frame the hits of every pattern in the order of the list.</summary>
    /// <exception cref="ArgumentException">A hit has no <see cref="TranslatedPatternHit.EmbossStart"/> / <see cref="TranslatedPatternHit.EmbossEnd"/> or no translation.</exception>
    public static EmbossReport Fuzztran(EmbossReportSequence sequence, IReadOnlyList<EmbossPatternHits<TranslatedPatternHit>> patterns)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(patterns);
        foreach (var p in patterns)
        {
            for (int i = 0; i < p.Hits.Count; i++)
            {
                if (p.Hits[i].EmbossStart is null || p.Hits[i].EmbossEnd is null)
                    throw new ArgumentException($"Hit {i} has no EmbossStart / EmbossEnd: the report prints the program's own coordinates (the hits of EmbossMotifPrograms.FindTranslatedPatternHits carry them).", nameof(patterns));
                if (p.Hits[i].Translation is null)
                    throw new ArgumentException($"Hit {i} has no Translation.", nameof(patterns));
            }
        }

        var frames = patterns.SelectMany(p => p.Hits.Select(h => h.Frame)).Distinct().OrderBy(f => f < 0 ? 3 - f : f).ToList();   // 1, 2, 3, -1, -2, -3
        var features = frames.SelectMany(frame => patterns.SelectMany(p => ReportOrder(p.Hits.Where(h => h.Frame == frame), p.Mismatches > 0 ? "" : p.Pattern, h => h.ProteinStart, h => h.ProteinEnd).Select(h =>
        {
            int number = Math.Abs(h.Frame);
            string pat = p.Name + ":" + p.Pattern;
            var notes = new List<EmbossFeatureTag> { Note("*pat " + pat) };
            if (h.Mismatches > 0)
                notes.Add(Note("*mismatch " + Int(h.Mismatches)));
            notes.AddRange([Note("*frame " + Int(number)), Note("*start " + Int(h.ProteinStart)), Note("*end " + Int(h.ProteinEnd)), Note("*translation " + h.Translation)]);
            return new EmbossFeature("SO:0000110", "fuzztran", h.EmbossStart!.Value, h.EmbossEnd!.Value, h.ProteinEnd - h.ProteinStart + 1 - h.Mismatches, h.Frame < 0 ? '-' : '+', number, notes);
        }))).ToList();
        return new EmbossReport(sequence, features, false, "N", TableLength: EmbossReportWriter.UnsignedRegionEnd(features), Tags: EmbossReportTagLists.TagList("fuzztran"));
    }
}
