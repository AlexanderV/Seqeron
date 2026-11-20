using System.Globalization;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The feature reports of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem (audit A5-8, B09 F176): the program's hits as <see cref="EmbossReport"/> (what each program stores in its feature table: source, type, coordinates, score,
/// notes), written by the one canonical <see cref="EmbossReportWriter"/> (<c>-rformat gff</c>, <c>excel</c>). Verified against the installed programs on 120 random runs per program and format
/// (TestData/Alignment/emboss_report_oracle.tsv, scripts/alignment/generate_emboss_report_oracle.py).
/// </summary>
public static class EmbossReports
{
    private static string Int(long v) => v.ToString(CultureInfo.InvariantCulture);

    // a pattern with a range element X(n,m): the program lists the hits of one frame / strand with the START DESCENDING (the hits of one start by ascending length), not in position order (measured on fuzznuc, fuzzpro and fuzztran 6.6.0)
    private static readonly System.Text.RegularExpressions.Regex RangeElement = new(@"\(\s*\d+\s*,\s*\d+\s*\)", System.Text.RegularExpressions.RegexOptions.Compiled | System.Text.RegularExpressions.RegexOptions.CultureInvariant);

    // the hits of the complement strand are listed in the order of the reverse strand: END descending, then start ascending (the longer hit first)
    private static IEnumerable<PatternHit> ReverseOrder(IEnumerable<PatternHit> hits, string pattern)
        => RangeElement.IsMatch(pattern) ? hits.OrderByDescending(h => h.Position + h.Length).ThenBy(h => h.Position) : hits;

    private static IEnumerable<T> ReportOrder<T>(IEnumerable<T> hits, string pattern, Func<T, int> start, Func<T, int> length)
        => RangeElement.IsMatch(pattern) ? hits.OrderByDescending(start).ThenBy(length) : hits;

    /// <summary>The report of <c>fuzznuc</c>: forward hits in position order (a range pattern searched without mismatches: start descending), then the hits of the complement strand; score = hit length minus mismatches; the pattern is printed as <c><paramref name="patternName"/>:<paramref name="pattern"/></c>.</summary>
    public static EmbossReport Fuzznuc(string sequenceName, string sequence, string pattern, IReadOnlyList<PatternHit> hits, string patternName = "pattern", int maxMismatches = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(hits);
        string orderPattern = maxMismatches > 0 ? "" : pattern;          // a search with mismatches (patOUBrute) lists a range pattern's hits in position order
        var features = ReportOrder(hits.Where(h => !h.IsReverseStrand), orderPattern, h => h.Position, h => h.Length).Concat(ReverseOrder(hits.Where(h => h.IsReverseStrand), orderPattern)).Select(h =>
        {
            string strand = h.IsReverseStrand ? "-" : "+";
            var notes = new List<string> { "note=*pat " + patternName + ":" + pattern };
            if (h.Mismatches > 0)
                notes.Add("note=*mismatch " + Int(h.Mismatches));
            return new EmbossReportFeature(h.Position + 1, h.Position + h.Length, h.Position + 1, h.Position + h.Length, Int(h.Length - h.Mismatches), strand, strand, ".",
                [patternName + ":" + pattern, h.Mismatches > 0 ? Int(h.Mismatches) : "."], notes);
        }).ToList();
        return new EmbossReport("fuzznuc", sequenceName, EmbossSequenceReader.Clean(sequence).Length, "nucleotide_motif", false, ["Pattern", "Mismatch"], features);
    }

    /// <summary>The report of <c>fuzzpro</c>: the hits in position order, score = hit length; the gff report has no strand and its <c>##sequence-region</c> end is the largest hit end (the feature table has no length).</summary>
    public static EmbossReport Fuzzpro(string sequenceName, string sequence, string pattern, IReadOnlyList<ProteinPatternHit> hits, string patternName = "pattern")
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(hits);
        var features = ReportOrder(hits, pattern, h => h.Position, h => h.Length).Select(h => new EmbossReportFeature(h.Position + 1, h.Position + h.Length, h.Position + 1, h.Position + h.Length, Int(h.Length), "+", ".", ".",
            [patternName + ":" + pattern, "."], ["note=*pat " + patternName + ":" + pattern])).ToList();
        return new EmbossReport("fuzzpro", sequenceName, EmbossReportWriter.UnsignedRegionEnd(features), "polypeptide_motif", true, ["Pattern", "Mismatch"], features);
    }

    /// <summary>
    /// The report of <c>fuzztran</c>: the hits frame by frame (1, 2, 3, -1, -2, -3), the program's own coordinates (<see cref="TranslatedPatternHit.EmbossStart"/> / <see cref="TranslatedPatternHit.EmbossEnd"/>: negative for the reverse frames),
    /// score = residues of the hit, the frame number without its sign (the strand column carries it). The gff report lists the smaller coordinate first unless the two have opposite signs, and its <c>##sequence-region</c> end is the largest
    /// coordinate read as an unsigned integer (the program's feature table has no length).
    /// </summary>
    public static EmbossReport Fuzztran(string sequenceName, string pattern, IReadOnlyList<TranslatedPatternHit> hits, string patternName = "pattern")
    {
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(hits);
        var features = hits.GroupBy(h => h.Frame).SelectMany(g => ReportOrder(g, pattern, h => h.ProteinStart, h => h.ProteinEnd)).Select(h =>
        {
            long s = h.EmbossStart!.Value, e = h.EmbossEnd!.Value;
            bool swap = s > e && !(e < 0 && s >= 0);
            int frame = Math.Abs(h.Frame);
            string pat = patternName + ":" + pattern;
            return new EmbossReportFeature(s, e, swap ? e : s, swap ? s : e, Int(h.ProteinEnd - h.ProteinStart + 1), h.Frame < 0 ? "-" : "+", h.Frame < 0 ? "-" : "+", Int(frame - 1),
                [pat, ".", Int(frame), Int(h.ProteinStart), Int(h.ProteinEnd), h.Translation],
                ["note=*pat " + pat, "note=*frame " + Int(frame), "note=*start " + Int(h.ProteinStart), "note=*end " + Int(h.ProteinEnd), "note=*translation " + h.Translation]);
        }).ToList();
        return new EmbossReport("fuzztran", sequenceName, EmbossReportWriter.UnsignedRegionEnd(features), "sequence_feature", false, ["Pattern", "Mismatch", "Frame", "PStart", "PEnd", "Translation"], features);
    }

    /// <summary>The report of <c>etandem</c> (the program's report order: highest score first): size, count, identity (one decimal) and the consensus.</summary>
    public static EmbossReport Etandem(string sequenceName, string sequence, IReadOnlyList<EmbossTandemRepeat> hits)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(hits);
        var features = hits.Select(h => new EmbossReportFeature(h.Start, h.End, h.Start, h.End, Int(h.Score), "+", "+", ".",
            [Int(h.RepeatSize), Int(h.RepeatCount), CFormat.Fixed(h.Identity, 1), h.Consensus ?? ""],
            ["rpt_type=TANDEM", "note=*rpt_size " + Int(h.RepeatSize), "note=*rpt_count " + Int(h.RepeatCount), "note=*identity " + CFormat.Fixed(h.Identity, 1), "note=*consensus " + h.Consensus])).ToList();
        return new EmbossReport("etandem", sequenceName, EmbossSequenceReader.Clean(sequence).Length, "repeat_region", false, ["Size", "Count", "Identity", "Consensus"], features, WriteWhenEmpty: true);
    }

    /// <summary>The report of <c>equicktandem</c>: size and count.</summary>
    public static EmbossReport Equicktandem(string sequenceName, string sequence, IReadOnlyList<EmbossTandemRepeat> hits)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(hits);
        var features = hits.Select(h => new EmbossReportFeature(h.Start, h.End, h.Start, h.End, Int(h.Score), "+", "+", ".",
            [Int(h.RepeatSize), Int(h.RepeatCount)], ["rpt_type=TANDEM", "note=*rpt_size " + Int(h.RepeatSize), "note=*rpt_count " + Int(h.RepeatCount)])).ToList();
        return new EmbossReport("equicktandem", sequenceName, EmbossSequenceReader.Clean(sequence).Length, "repeat_region", false, ["Size", "Count"], features, WriteWhenEmpty: true);
    }
}
