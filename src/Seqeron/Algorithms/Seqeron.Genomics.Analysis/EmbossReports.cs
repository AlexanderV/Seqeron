using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The feature tables of fuzznuc, fuzzpro, fuzztran, etandem and equicktandem (audit A5-8 / A5-8b, B09 F176 / F184): the program's hits as <see cref="EmbossReport"/> (what each program stores in its feature table: source, type, coordinates,
/// score, tags), written in every report format by the one canonical <see cref="EmbossReportWriter"/>. Verified against the installed programs
/// (TestData/Alignment/emboss_report_formats_oracle.jsonl.gz, scripts/alignment/generate_emboss_report_formats_oracle.py).
/// </summary>
public static class EmbossReports
{
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
            return hits.OrderBy(h => key(h).Start + key(h).Length).ThenBy(h => key(h).Start);
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
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(patternName);
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(hits);
        string orderPattern = maxMismatches > 0 ? "" : pattern;          // a search with mismatches (patOUBrute) lists a range pattern's hits in position order
        // a reversed sequence is searched as its reverse complement: the hits of the first pass come out on the minus strand and their positions are mirrored back
        bool rev = sequence.Reversed;
        int n = sequence.Length;
        (int Start, int Length) Key(PatternHit h) => rev ? (n - h.Position - h.Length, h.Length) : (h.Position, h.Length);
        var features = ReportOrder(hits.Where(h => h.IsReverseStrand == rev), orderPattern, h => Key(h).Start, h => h.Length).Concat(ReverseOrder(hits.Where(h => h.IsReverseStrand != rev), pattern, maxMismatches, Key)).Select(h =>
        {
            var notes = new List<EmbossFeatureTag> { Note("*pat " + patternName + ":" + pattern) };
            if (h.Mismatches > 0)
                notes.Add(Note("*mismatch " + Int(h.Mismatches)));
            return new EmbossFeature("SO:0000714", "fuzznuc", h.Position + 1, h.Position + h.Length, h.Length - h.Mismatches, Strand(h.IsReverseStrand, rev), 0, notes);
        }).ToList();
        return new EmbossReport(sequence, features, false, "N", TableLength: sequence.Length, Tags: TagList("fuzznuc"));
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
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(patternName);
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(hits);
        var features = ReportOrder(hits, pattern, h => h.Position, h => h.Length).Select(h => new EmbossFeature("SO:0001067", "fuzzpro", h.Position + 1, h.Position + h.Length, h.Length, '\0', 0,
            [Note("*pat " + patternName + ":" + pattern)])).ToList();
        return new EmbossReport(sequence, features, true, first ? "P" : "p", TableLength: EmbossReportWriter.UnsignedRegionEnd(features), Tags: TagList("fuzzpro"));
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
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(patternName);
        ArgumentNullException.ThrowIfNull(pattern);
        ArgumentNullException.ThrowIfNull(hits);
        for (int i = 0; i < hits.Count; i++)
        {
            if (hits[i].EmbossStart is null || hits[i].EmbossEnd is null)
                throw new ArgumentException($"Hit {i} has no EmbossStart / EmbossEnd: the report prints the program's own coordinates (the hits of GenomicAnalyzer.FindTranslatedPatternHits carry them).", nameof(hits));
            if (hits[i].Translation is null)
                throw new ArgumentException($"Hit {i} has no Translation.", nameof(hits));
        }

        var features = hits.GroupBy(h => h.Frame).SelectMany(g => ReportOrder(g, pattern, h => h.ProteinStart, h => h.ProteinEnd)).Select(h =>
        {
            int frame = Math.Abs(h.Frame);
            string pat = patternName + ":" + pattern;
            return new EmbossFeature("SO:0000110", "fuzztran", h.EmbossStart!.Value, h.EmbossEnd!.Value, h.ProteinEnd - h.ProteinStart + 1, h.Frame < 0 ? '-' : '+', frame,
                [Note("*pat " + pat), Note("*frame " + Int(frame)), Note("*start " + Int(h.ProteinStart)), Note("*end " + Int(h.ProteinEnd)), Note("*translation " + h.Translation)]);
        }).ToList();
        return new EmbossReport(sequence, features, false, "N", TableLength: EmbossReportWriter.UnsignedRegionEnd(features), Tags: TagList("fuzztran"));
    }

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
        return new EmbossReport(sequence, features, false, "N", TableLength: sequence.Length, Tags: TagList("etandem"), WriteWhenEmpty: true);
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
        return new EmbossReport(sequence, features, false, "N", TableLength: sequence.Length, Tags: TagList("equicktandem"), WriteWhenEmpty: true);
    }

    /// <summary>The report tags (the ACD <c>taglist</c>) of the programs.</summary>
    public static IReadOnlyList<EmbossReportTag> TagList(string program) => program switch
    {
        "fuzznuc" or "fuzzpro" => [new("str", "pat", "Pattern"), new("int", "mismatch", "Mismatch")],
        "fuzztran" => [new("str", "pat", "Pattern"), new("int", "mismatch", "Mismatch"), new("int", "frame", "Frame"), new("int", "start", "PStart"), new("int", "end", "PEnd"), new("str", "translation", "Translation")],
        "etandem" => [new("int", "rpt_size", "Size"), new("int", "rpt_count", "Count"), new("float", "identity", "Identity"), new("str", "consensus", "Consensus")],
        "equicktandem" => [new("int", "rpt_size", "Size"), new("int", "rpt_count", "Count")],
        _ => throw new ArgumentException($"Unknown program '{program}'.", nameof(program)),
    };
}
