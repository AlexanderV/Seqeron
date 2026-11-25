// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.FeatureTools, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S107, S3776, S1121, S3358, S3398, S1066, S1854

using System.Globalization;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.FeatureTools;

using Seqeron.Genomics.Analysis;

/// <summary>A run of <c>twofeat</c> (A5-8c3h, B09 F263): the options of its ACD file.</summary>
/// <param name="Sequence">The value of <c>-sequence</c> (the USA).</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers.</param>
/// <param name="Report">The report qualifiers (<c>-rformat2</c> ...); null = the defaults (format table).</param>
/// <param name="CommandLine">The entries of the command line the report header echoes.</param>
/// <param name="RunDate">The date of the run; null = now.</param>
/// <param name="ReportFile">The name of the report file (<c>Report_file:</c>).</param>
public sealed record EmbossTwofeatRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null, EmbossReportQualifiers? Report = null, IReadOnlyList<string>? CommandLine = null,
    DateTime? RunDate = null, string ReportFile = "stdout")
{
    /// <summary><c>-asource</c>.</summary>
    public string ASource { get; init; } = "*";

    /// <summary><c>-atype</c>.</summary>
    public string AType { get; init; } = "*";

    /// <summary><c>-asense</c>: 0, + or - (the minus sign is read as 0 by the program).</summary>
    public string ASense { get; init; } = "0";

    /// <summary><c>-aminscore</c> (the program takes the integer part).</summary>
    public float AMinScore { get; init; }

    /// <summary><c>-amaxscore</c>.</summary>
    public float AMaxScore { get; init; }

    /// <summary><c>-atag</c>.</summary>
    public string ATag { get; init; } = "*";

    /// <summary><c>-avalue</c>.</summary>
    public string AValue { get; init; } = "*";

    /// <summary><c>-bsource</c>.</summary>
    public string BSource { get; init; } = "*";

    /// <summary><c>-btype</c>.</summary>
    public string BType { get; init; } = "*";

    /// <summary><c>-bsense</c>.</summary>
    public string BSense { get; init; } = "0";

    /// <summary><c>-bminscore</c>.</summary>
    public float BMinScore { get; init; }

    /// <summary><c>-bmaxscore</c>.</summary>
    public float BMaxScore { get; init; }

    /// <summary><c>-btag</c>.</summary>
    public string BTag { get; init; } = "*";

    /// <summary><c>-bvalue</c>.</summary>
    public string BValue { get; init; } = "*";

    /// <summary><c>-overlap</c>: A any, O overlap, NO none, NW overlap but not within, AW A within B, BW B within A.</summary>
    public string Overlap { get; init; } = "A";

    /// <summary><c>-minrange</c>.</summary>
    public int MinRange { get; init; }

    /// <summary><c>-maxrange</c>.</summary>
    public int MaxRange { get; init; }

    /// <summary><c>-rangetype</c>: N nearest ends, L left ends, R right ends, F furthest ends.</summary>
    public string RangeType { get; init; } = "N";

    /// <summary><c>-sense</c>: A any, S same, O opposite.</summary>
    public string Sense { get; init; } = "A";

    /// <summary><c>-order</c>: accepted (the program does not use it).</summary>
    public string Order { get; init; } = "A";

    /// <summary><c>-twoout</c>: the two features of a pair are reported instead of one feature that covers them.</summary>
    public bool TwoOut { get; init; }

    /// <summary><c>-typeout</c>: the type of the new feature (default misc_feature).</summary>
    public string TypeOut { get; init; } = "misc_feature";

    /// <summary>The texts of the configuration files that define the databases of <c>db:id</c> USAs.</summary>
    public IReadOnlyList<string>? Definitions { get; init; }

    /// <summary>The report file is not named (no <c>-outfile</c>): it is <c>BASE.twofeat</c> with the base of <c>-rname</c> or of the name of the sequence and the extension of <c>-rextension</c> (<see cref="ReportFile"/> is not used).</summary>
    public bool DefaultReportFile { get; init; }

    /// <summary>The directories (without trailing slash) that exist where the program runs; a report file whose directory is not among them cannot be opened. Null: every directory exists.</summary>
    public IReadOnlyCollection<string>? ExistingDirectories { get; init; }
}

/// <summary>
/// <c>twofeat</c> of EMBOSS 6.6.0 (<c>emboss/twofeat.c</c>): the features of a sequence that match the criteria of feature A and the ones that match feature B are paired when they are in the required relation (overlap, distance, sense, order) and the pairs
/// are reported as features that cover them (with the tags startA, endA, startB, endB) or as the two features.
/// </summary>
public static class EmbossTwofeatProgram
{
    private const int OvAny = 0, OvOverlap = 1, OvNotOver = 2, OvNotIn = 3, OvAin = 4, OvBin = 5;
    private const int RaNearest = 0, RaLeft = 1, RaRight = 2, RaFurthest = 3;
    private const int SnAny = 0, SnSame = 1, SnOpposite = 2;

    private sealed record Hit(EmbossSeqFeature GfA, EmbossSeqFeature GfB, int Start, int End);

    /// <summary>Runs <c>twofeat</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    public static EmbossFeatureReportOutput Twofeat(EmbossTwofeatRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "twofeat" };
        var seqall = EmbossSeqall.Open("twofeat", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, true, context, "any", false, r.Definitions);
        if (seqall is null)
        {
            var died = EmbossSeqretPrograms.Died(context);
            return new EmbossFeatureReportOutput("", died.Messages) { Died = true };
        }

        var options = EmbossFeatureReports.Options("twofeat", "table", r.Report, r.CommandLine, r.RunDate, r.ReportFile, context, true);
        if (EmbossFeatureReports.ValidateTags(options, context))
            return new EmbossFeatureReportOutput("", [.. context.Messages, "Died: Report option -outfile: Validation failed"]) { Died = true };
        string? reportFile = EmbossOutPath.OpenReport(context, "twofeat", r.DefaultReportFile ? null : r.ReportFile, r.Report, EmbossOutPath.DefaultBaseName(seqall.Peek().Name, r.Sequence), r.ExistingDirectories);
        if (reportFile is null)
        {
            var died = EmbossSeqretPrograms.Died(context, "outfile");
            return new EmbossFeatureReportOutput("", died.Messages) { Died = true };
        }

        options = options with { ReportFile = reportFile };

        // the program converts the scores to integers and the senses and the relation codes to numbers
        float aminscore = (int)r.AMinScore, amaxscore = (int)r.AMaxScore, bminscore = (int)r.BMinScore, bmaxscore = (int)r.BMaxScore;
        int asensei = r.ASense == "+" ? 1 : 0;
        int bsensei = r.BSense == "+" ? 1 : 0;
        int overlapi = r.Overlap switch { "A" => OvAny, "O" => OvOverlap, "NO" => OvNotOver, "NW" => OvNotIn, "AW" => OvAin, "BW" => OvBin, _ => throw new InvalidOperationException($"Unknown Overlap code: {r.Overlap}") };
        int rangetypei = r.RangeType switch { "N" => RaNearest, "L" => RaLeft, "R" => RaRight, "F" => RaFurthest, _ => throw new InvalidOperationException($"Unknown rangetype code: {r.RangeType}") };
        int sensei = r.Sense switch { "A" => SnAny, "S" => SnSame, "O" => SnOpposite, _ => throw new InvalidOperationException($"Unknown sense code: {r.Sense}") };

        var tables = new List<EmbossReport>();
        long count = 0, total = 0;
        try
        {
            while (EmbossSeqretPrograms.TryNext(seqall, context, out var seq))
            {
                count++;
                total += EmbossSeqRange.GetEnd(seq) - EmbossSeqRange.GetBegin(seq) + 1;                // ajSeqGetLenTrimmed
                int begin = seqall.GetSeqBegin(seq), end = seqall.GetSeqEnd(seq);
                var tabA = NewTable(seq, context);
                var tabB = NewTable(seq, context);
                FindFeatures(seq, tabA, begin, end, r.ASource, r.AType, asensei, aminscore, amaxscore, r.ATag, r.AValue, context);
                FindFeatures(seq, tabB, begin, end, r.BSource, r.BType, bsensei, bminscore, bmaxscore, r.BTag, r.BValue, context);
                if (tabA.Features.Count > 0 && tabB.Features.Count > 0)
                {
                    var outtab = NewTable(seq, context);
                    Rippledown(tabA, tabB, overlapi, r.MinRange, r.MaxRange, rangetypei, sensei, outtab, r.TwoOut, r.TypeOut);
                    tables.Add(EmbossFeatureReports.Build(seq, outtab, options.Tags));
                }
            }
        }
        catch (EmbossFatalException e)
        {
            context.Messages.Add(e.Message);
            return new EmbossFeatureReportOutput("", [.. context.Messages]) { Died = true };
        }

        string text = EmbossFeatureReports.Write(options, tables, context, _ => (count, total), out bool failed);
        // ajReportClose: a report without a hit is a warning
        long reported = 0;
        foreach (var table in tables)
        {
            long n = table.Features.Count;
            if (options.MaxSeq != 0)
                n = Math.Min(n, options.MaxSeq);
            if (options.MaxAll != 0)
            {
                long remaining = options.MaxAll - reported;
                if (remaining <= 0)
                    break;
                n = Math.Min(n, remaining);
            }

            reported += n;
        }

        if (reported == 0 || failed)
            context.Warn($"No hits in output report file '{reportFile}'");
        return new EmbossFeatureReportOutput(text, [.. context.Messages]) { Died = context.Fatal, ReportFile = reportFile };
    }

    // ajFeattableNewSeq
    private static EmbossSeqFeatureTable NewTable(EmbossSeq seq, EmbossFeatureReadContext context)
        => new(context) { SeqId = seq.Name, Type = seq.Type == "P" ? "P" : "N", Length = seq.Seq.Length };

    // twofeat_find_features: the local features in the range of the stream that match the criteria are copied into the table
    private static void FindFeatures(EmbossSeq seq, EmbossSeqFeatureTable tab, int begin, int end, string source, string type, int sense, float minscore, float maxscore, string tag, string value, EmbossFeatureReadContext context)
    {
        _ = context;
        var seqtab = seq.Fttable;
        if (seqtab is null || seqtab.Features.Count == 0)
            return;
        bool tagsmatch = false;
        foreach (var gf in seqtab.Features)
        {
            if (!IsLocalRange(gf, begin, end))
                continue;
            if (MatchFeature(gf, source, type, sense, minscore, maxscore, tag, value, ref tagsmatch))
                tab.Add(tab.NewFeat(gf));
        }
    }

    // ajFeatIsLocalRange (unsigned)
    private static bool IsLocalRange(EmbossSeqFeature gf, int start, int end)
    {
        if (gf.Flags.HasFlag(EmbossLocationBits.RemoteId) || gf.Flags.HasFlag(EmbossLocationBits.Label))
            return false;
        return !(unchecked((uint)gf.End) < unchecked((uint)start) || unchecked((uint)gf.Start) > unchecked((uint)end));
    }

    // twofeat_MatchFeature
    private static bool MatchFeature(EmbossSeqFeature gf, string source, string type, int sense, float minscore, float maxscore, string tag, string value, ref bool tagsmatch)
    {
        bool scoreok = minscore < maxscore;
        if (!gf.Flags.HasFlag(EmbossLocationBits.Multiple))
            tagsmatch = MatchPatternTags(gf, tag, value);
        if (gf.Flags.HasFlag(EmbossLocationBits.RemoteId))
            return false;
        char strand = StrandOf(gf);
        if (!MatchPatternDelim(gf.Source, source, ",;|") || (strand == '+' && sense == -1) || (strand == '-' && sense == 1) || (scoreok && gf.Score < minscore) || (scoreok && gf.Score > maxscore) || !tagsmatch)
            return false;
        if (type.Length > 0)
        {
            bool val = false;
            foreach (string key in type.Split([' ', '\t', '\n', '\r', ',', ';', '|'], StringSplitOptions.RemoveEmptyEntries))
            {
                if (EmbossExtractfeatProgram.TypeMatchWild(gf, key))
                {
                    val = true;
                    break;
                }
            }

            if (!val)
                return false;
        }

        return true;
    }

    // twofeat_MatchPatternTags
    private static bool MatchPatternTags(EmbossSeqFeature feat, string tpattern, string vpattern)
    {
        if (tpattern == "*" && vpattern == "*")
            return true;
        foreach (var (tagnam, tagval) in feat.Tags)
        {
            bool tval = MatchPatternDelim(tagnam, tpattern, ",;|");
            bool vval = tagval.Length == 0 ? vpattern == "*" : MatchPatternDelim(tagval, vpattern, ",;|");
            if (tval && vval)
                return true;
        }

        return false;
    }

    // embMiscMatchPatternDelimC
    private static bool MatchPatternDelim(string text, string pattern, string delims)
    {
        char[] sep = [.. " \t\n\r", .. delims];
        foreach (string key in pattern.Split(sep, StringSplitOptions.RemoveEmptyEntries))
        {
            if (EmbossWild.MatchCase(text, key))
                return true;
        }

        return false;
    }

    // ajFeatGetStrand: '.' for a protein feature, else '-' or '+'
    private static char StrandOf(EmbossSeqFeature gf) => gf.Protein ? '.' : gf.Strand == '-' ? '-' : '+';

    private static string TypeOf(EmbossSeqFeature gf) => EmbossFeatureVocab.TypeExternal((gf.Protein ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna).Types, gf.Type);

    // twofeat_rippledown
    private static void Rippledown(EmbossSeqFeatureTable tabA, EmbossSeqFeatureTable tabB, int overlapi, int minrange, int maxrange, int rangetypei, int sensei, EmbossSeqFeatureTable outtab, bool twoout, string typeout)
    {
        var hits = new List<Hit>();
        foreach (var gfA in tabA.Features)
        {
            foreach (var gfB in tabB.Features)
            {
                if (CheckMatch(gfA, gfB, overlapi, minrange, maxrange, rangetypei, sensei) is { } hit)
                    hits.Insert(0, hit);                          // ajListPush adds at the head of the list: the hits are reported in the reverse order they were found
            }
        }

        ReportHits(hits, twoout, typeout, outtab);
    }

    // twofeat_report_hits
    private static void ReportHits(List<Hit> hits, bool twoout, string typeout, EmbossSeqFeatureTable outtab)
    {
        foreach (var detail in hits)
        {
            if (twoout)
            {
                outtab.Add(outtab.NewFeat(detail.GfA));
                outtab.Add(outtab.NewFeat(detail.GfB));
                continue;
            }

            // if both features are -ve then output this, else +ve
            char strand = detail.GfA.Strand == '-' && detail.GfB.Strand == '-' ? '-' : '+';
            var feature = outtab.Type == "P"
                ? outtab.NewProt("twofeat", typeout, detail.Start, detail.End, 0f, EmbossLocationBits.None, null)
                : outtab.NewNuc("twofeat", typeout, detail.Start, detail.End, 0f, strand, 0, 0, 0, 0, null, null, EmbossLocationBits.None, null);
            if (feature is null)
                continue;
            outtab.TagAdd(feature, "note", "*startA " + detail.GfA.GetStart().ToString(CultureInfo.InvariantCulture));
            outtab.TagAdd(feature, "note", "*endA " + detail.GfA.GetEnd().ToString(CultureInfo.InvariantCulture));
            outtab.TagAdd(feature, "note", "*startB " + detail.GfB.GetStart().ToString(CultureInfo.InvariantCulture));
            outtab.TagAdd(feature, "note", "*endB " + detail.GfB.GetEnd().ToString(CultureInfo.InvariantCulture));
        }
    }

    // twofeat_check_match
    private static Hit? CheckMatch(EmbossSeqFeature gfA, EmbossSeqFeature gfB, int overlapi, int minrange, int maxrange, int rangetypei, int sensei)
    {
        int sA = gfA.GetStart(), eA = gfA.GetEnd(), sB = gfB.GetStart(), eB = gfB.GetEnd();

        // distances from each end
        int ss = Math.Abs(sA - sB), ee = Math.Abs(eA - eB), se = Math.Abs(sA - eB), es = Math.Abs(eA - sB);
        int distance;
        int tmp1, tmp2;
        switch (rangetypei)
        {
            case RaNearest:
                tmp1 = ss < ee ? ss : ee;
                tmp2 = se < es ? se : es;
                distance = tmp1 < tmp2 ? tmp1 : tmp2;
                break;
            case RaLeft:
                distance = ss;
                break;
            case RaRight:
                distance = ee;
                break;
            default:
                tmp1 = ss > ee ? ss : ee;
                tmp2 = se > es ? se : es;
                distance = tmp1 > tmp2 ? tmp1 : tmp2;
                break;
        }

        // the distance criteria are ignored when the minimum is not less than the maximum
        if (minrange < maxrange)
        {
            if (distance < minrange)
                return null;
            if (distance > maxrange)
                return null;
        }

        switch (overlapi)
        {
            case OvAny:
                break;
            case OvOverlap:                                   // sA in B or sB in A
                if ((sA < sB || sA > eB) && (sB < sA || sB > eA))
                    return null;
                break;
            case OvNotOver:
                if ((sA >= sB && sA <= eB) || (sB >= sA && sB <= eA))
                    return null;
                break;
            case OvNotIn:
                if ((sA <= sB && eA >= eB) || (sB <= sA && eB >= eA))
                    return null;
                break;
            case OvAin:
                if (sA < sB || eA > eB)
                    return null;
                break;
            case OvBin:
                if (sB < sA || eB > eA)
                    return null;
                break;
        }

        // the sense criteria: the case of the opposite sense falls through to the default of the switch, which is a fatal error of the library
        switch (sensei)
        {
            case SnAny:
                break;
            case SnSame:
                if (StrandOf(gfA) != StrandOf(gfB))
                    return null;
                break;
            case SnOpposite:
                if (StrandOf(gfA) == StrandOf(gfB))
                    return null;
                throw new EmbossFatalException($"\n   EMBOSS An error in twofeat.c at line {FatalLineSense}:\nUnknown sense type: {sensei}");
        }

        // the order criteria are switched on the sense (the codes of the sense are the codes of the order): the same sense means "A then B", measured from the middle of each feature
        switch (sensei)
        {
            case SnAny:
                break;
            case 1:
                tmp1 = (eA - sA) / 2;
                tmp2 = (eB - sB) / 2;
                if (tmp1 > tmp2)
                    return null;
                break;
        }

        // a match of a feature to itself (same position, sense and type) is rejected
        if (ss == 0 && ee == 0 && StrandOf(gfA) == StrandOf(gfB) && TypeOf(gfA) == TypeOf(gfB))
            return null;
        return new Hit(gfA, gfB, sA < sB ? sA : sB, eA > eB ? eA : eB);
    }

    private const int FatalLineSense = 932;
}
