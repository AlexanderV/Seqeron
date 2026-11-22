#pragma warning disable S107, S3776, S1121, S3358, S3398, S1066, S2326

using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>A run of <c>showfeat</c> (A5-8c3h, B09 F260): the options of its ACD file.</summary>
/// <param name="Sequence">The value of <c>-sequence</c> (the USA).</param>
/// <param name="Files">The input files by name.</param>
/// <param name="Qualifiers">The sequence input qualifiers.</param>
public sealed record EmbossShowfeatRequest(string Sequence, IReadOnlyDictionary<string, string> Files, EmbossSequenceQualifiers? Qualifiers = null)
{
    /// <summary><c>-sourcematch</c>: sources to show (default <c>*</c>).</summary>
    public string SourceMatch { get; init; } = "*";

    /// <summary><c>-typematch</c>: types to show (default <c>*</c>).</summary>
    public string TypeMatch { get; init; } = "*";

    /// <summary><c>-tagmatch</c>: tags to show (default <c>*</c>).</summary>
    public string TagMatch { get; init; } = "*";

    /// <summary><c>-valuematch</c>: tag values to show (default <c>*</c>).</summary>
    public string ValueMatch { get; init; } = "*";

    /// <summary><c>-sort</c>: source, start, type or nosort.</summary>
    public string Sort { get; init; } = "start";

    /// <summary><c>-joinfeatures</c>: the parts of a multiple-location feature on one line.</summary>
    public bool JoinFeatures { get; init; }

    /// <summary><c>-html</c>.</summary>
    public bool Html { get; init; }

    /// <summary><c>-id</c> (default true).</summary>
    public bool Id { get; init; } = true;

    /// <summary><c>-description</c> (default true).</summary>
    public bool Description { get; init; } = true;

    /// <summary><c>-scale</c> (default true).</summary>
    public bool Scale { get; init; } = true;

    /// <summary><c>-width</c> (default 60).</summary>
    public int Width { get; init; } = 60;

    /// <summary><c>-collapse</c>.</summary>
    public bool Collapse { get; init; }

    /// <summary><c>-forward</c> (default true).</summary>
    public bool Forward { get; init; } = true;

    /// <summary><c>-reverse</c> (default true).</summary>
    public bool Reverse { get; init; } = true;

    /// <summary><c>-unknown</c> (default true).</summary>
    public bool Unknown { get; init; } = true;

    /// <summary><c>-strand</c>.</summary>
    public bool Strand { get; init; }

    /// <summary><c>-sourceshow</c>.</summary>
    public bool SourceShow { get; init; }

    /// <summary><c>-position</c>.</summary>
    public bool Position { get; init; }

    /// <summary><c>-typeshow</c> (default true).</summary>
    public bool TypeShow { get; init; } = true;

    /// <summary><c>-tagshow</c>.</summary>
    public bool TagShow { get; init; }

    /// <summary><c>-valueshow</c> (default true).</summary>
    public bool ValueShow { get; init; } = true;

    /// <summary><c>-stricttags</c>.</summary>
    public bool StrictTags { get; init; }

    /// <summary><c>-annotation</c>: ranges with a text to show as annotation (<c>1-10 text</c> ...).</summary>
    public string Annotation { get; init; } = "";

    /// <summary>The texts of the configuration files that define the databases of <c>db:id</c> USAs.</summary>
    public IReadOnlyList<string>? Definitions { get; init; }
}

/// <summary>The output of <c>showfeat</c>.</summary>
/// <param name="Text">The text of the output file.</param>
/// <param name="Sequences">The number of sequences shown.</param>
/// <param name="Messages">What the library printed on its error stream.</param>
public sealed record EmbossShowfeatOutput(string Text, int Sequences, IReadOnlyList<string> Messages)
{
    /// <summary>The program died ("Died: ...").</summary>
    public bool Died { get; init; }
}

/// <summary>
/// <c>showfeat</c> of EMBOSS 6.6.0 (<c>emboss/showfeat.c</c>): the features of the sequences drawn as lines of ASCII graphics under a scale, one line per feature, with the strand, source, positions, type and tags.
/// </summary>
public static class EmbossShowfeatProgram
{
    /// <summary>Runs <c>showfeat</c>.</summary>
    /// <exception cref="InvalidOperationException">The input is one on which the reference program never stops or dies by a fault.</exception>
    public static EmbossShowfeatOutput Showfeat(EmbossShowfeatRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.Sequence);
        ArgumentNullException.ThrowIfNull(r.Files);
        var context = new EmbossFeatureReadContext { Program = "showfeat" };
        var seqall = EmbossSeqall.Open("showfeat", r.Sequence, r.Files, r.Qualifiers ?? EmbossSequenceQualifiers.None, true, context, "any", false, r.Definitions);
        if (seqall is null)
        {
            var died = EmbossSeqretPrograms.Died(context);
            return new EmbossShowfeatOutput("", 0, died.Messages) { Died = true };
        }

        if (r.Width < 0)
        {
            // the ACD integer has the minimum 0
            context.Warn($"integer value out of range {r.Width} less than (reset to) 0");
            r = r with { Width = 0 };
        }

        var annotation = EmbossRangeSpec.Parse(r.Annotation);
        var sb = new StringBuilder();
        int count = 0;
        while (EmbossSeqretPrograms.TryNext(seqall, context, out var seq))
        {
            count++;
            int beg = EmbossSeqRange.GetBegin(seq) - 1;
            int end = EmbossSeqRange.GetEnd(seq) - 1;
            if (r.Id)
            {
                if (r.Html)
                    sb.Append("<H2>").Append(seq.Name).Append("</H2>\n");
                else
                    sb.Append(seq.Name).Append('\n');
            }

            if (r.Description)
            {
                if (r.Html)
                    sb.Append("<H3>").Append(seq.Desc).Append("</H3>\n");
                else
                    sb.Append(FmtWrap(seq.Desc, 80)).Append('\n');
            }

            if (r.Html)
                sb.Append("<PRE>");
            if (r.Scale && r.Width > 3)
            {
                sb.Append('|').Append('=', r.Width - 2).Append("| ");
                if (beg != 0 || end + 1 != seq.Seq.Length)
                    sb.Append(CultureInfo.InvariantCulture, $"{beg + 1}-{end + 1} of {seq.Seq.Length}\n");
                else
                    sb.Append(seq.Seq.Length).Append('\n');
            }

            ShowFeatSeq(sb, r, seq, beg, end, annotation);
            if (r.Html)
                sb.Append("</PRE>");
            sb.Append('\n');
        }

        return new EmbossShowfeatOutput(sb.ToString(), count, [.. context.Messages]) { Died = context.Fatal };
    }

    // ajStrFmtWrap
    private static string FmtWrap(string text, int width)
    {
        var s = new StringBuilder(text);
        if (width > s.Length)
            return text;
        int i = 0, imax = s.Length - width;
        while (i < imax)
        {
            int j = i + width + 1;
            if (j > s.Length)
                j = s.Length;
            int k = j;
            int cq = i;
            while (i < j)
            {
                char c = cq < s.Length ? s[cq] : '\0';
                if (c == ' ' || (c >= '\t' && c <= '\r'))
                {
                    k = i;
                    if (c == '\n')
                        break;
                }

                cq++;
                i++;
            }

            if ((cq < s.Length ? s[cq] : '\0') != '\n')
            {
                if (k == j)
                {
                    s.Insert(k, '\n');
                    imax++;
                }
                else
                {
                    s[k] = '\n';
                }
            }

            i = k + 1;
        }

        return s.ToString();
    }

    private sealed class Line
    {
        public StringBuilder Text = new();
        public char Strand;
        public string Source = "", Type = "";
        public readonly StringBuilder Tags = new();
        public readonly StringBuilder Pos = new();
    }

    // showfeat_ShowFeatSeq
    private static void ShowFeatSeq(StringBuilder outfile, EmbossShowfeatRequest r, EmbossSeq seq, int beg, int end, IReadOnlyList<EmbossRangeElement> annotation)
    {
        var feat = seq.Fttable;
        if (feat is null)
            return;
        var tagstmp = new StringBuilder();
        var cur = new Line { Text = new StringBuilder() };
        bool gotoutput = false;
        if (feat.Features.Count == 0)
            return;
        var fullfeats = feat.Features;
        var order = new int[fullfeats.Count];
        for (int i = 0; i < order.Length; i++)
            order[i] = i;
        Comparison<EmbossSeqFeature>? compare = r.Sort switch
        {
            "source" => CompareSource,
            "start" => ComparePos,
            "type" => CompareType,
            _ => null,
        };
        if (compare is not null)
            ToIndex(fullfeats, order, compare);
        IEnumerable<EmbossSeqFeature> sequence = r.Sort == "nosort" ? fullfeats : order.Select(i => fullfeats[i]);
        foreach (var gf in sequence)
        {
            tagstmp.Clear();
            if (!MatchPatternTags(gf, r.TagMatch, r.ValueMatch, r.StrictTags, tagstmp, r.ValueShow))
                continue;
            if (!MatchPatternDelim(gf.Source, r.SourceMatch, ",;|"))
                continue;
            if (r.TypeMatch.Length > 0)
            {
                bool val = false;
                foreach (string key in r.TypeMatch.Split([' ', '\t', '\n', '\r', ',', ';', '|'], StringSplitOptions.RemoveEmptyEntries))
                {
                    if (EmbossExtractfeatProgram.TypeMatchWild(gf, key))
                    {
                        val = true;
                        break;
                    }
                }

                if (!val)
                    continue;
            }

            // the feature is within the range to show
            if (beg + 1 > gf.GetEnd() || end + 1 < gf.GetStart())
                continue;

            if (gf.Subfeatures is { Count: > 0 } subs)
            {
                if (r.JoinFeatures)
                {
                    if (gotoutput)
                    {
                        FeatOut(outfile, cur, r);
                        gotoutput = false;
                    }

                    cur = New(r.Width, gf.Source, TypeOf(gf), StrandOf(gf));
                }

                foreach (var gfsub in subs)
                {
                    if (!r.Forward && StrandOf(gfsub) == '+')
                        continue;
                    if (!r.Reverse && StrandOf(gfsub) == '-')
                        continue;
                    if (!r.Unknown && StrandOf(gfsub) == '\0')
                        continue;
                    if (gfsub.Flags.HasFlag(EmbossLocationBits.RemoteId))
                        continue;
                    if (!r.JoinFeatures)
                    {
                        if (gotoutput)
                            FeatOut(outfile, cur, r);
                        cur = New(r.Width, gfsub.Source, TypeOf(gfsub), StrandOf(gfsub));
                    }

                    gotoutput = true;
                    cur.Tags.Append(tagstmp);
                    tagstmp.Clear();
                    AddPos(cur.Pos, gfsub.GetStart(), gfsub.GetEnd());
                    WriteFeat(cur, StrandCur(cur), gfsub.GetStart(), gfsub.GetEnd(), r.Width, beg, end);
                }

                if (gotoutput)
                {
                    FeatOut(outfile, cur, r);
                    gotoutput = false;
                }
            }
            else
            {
                if (!r.Forward && StrandOf(gf) == '+')
                    continue;
                if (!r.Reverse && StrandOf(gf) == '-')
                    continue;
                if (!r.Unknown && StrandOf(gf) == '\0')
                    continue;
                if (gf.Flags.HasFlag(EmbossLocationBits.RemoteId))
                    continue;
                if (!r.Collapse)
                {
                    if (gotoutput)
                        FeatOut(outfile, cur, r);
                    cur = New(r.Width, gf.Source, TypeOf(gf), StrandOf(gf));
                    gotoutput = true;
                }

                cur.Tags.Append(tagstmp);
                tagstmp.Clear();
                AddPos(cur.Pos, gf.GetStart(), gf.GetEnd());
                WriteFeat(cur, StrandCur(cur), gf.GetStart(), gf.GetEnd(), r.Width, beg, end);
            }
        }

        if (gotoutput)
            FeatOut(outfile, cur, r);

        if (annotation.Count > 0)
        {
            gotoutput = false;
            for (int count = 0; count < annotation.Count; count++)
            {
                var el = annotation[count];
                int rstart = el.Start, rend = el.End;
                if (beg + 1 > rend || end + 1 < rstart)
                    continue;
                string annText = el.Text;
                if (!r.Collapse || !string.Equals(annText, cur.Type, StringComparison.OrdinalIgnoreCase))
                {
                    if (gotoutput)
                        FeatOut(outfile, cur, r);
                    cur = New(r.Width, "Annotation", annText, '\0');
                    gotoutput = true;
                }

                AddPos(cur.Pos, rstart, rend);
                WriteFeat(cur, cur.Strand, rstart, rend, r.Width, beg, end);
            }

            if (gotoutput)
                FeatOut(outfile, cur, r);
        }
    }

    private static char StrandCur(Line l) => l.Strand;

    private static Line New(int width, string source, string type, char strand) => new() { Text = new StringBuilder(new string(' ', width)), Source = source, Type = type, Strand = strand };

    // ajFeatGetStrand: '.' for a protein feature, else '-' or '+'
    private static char StrandOf(EmbossSeqFeature gf) => gf.Protein ? '.' : gf.Strand == '-' ? '-' : '+';

    // ajFeatGetType: the type name of the vocabulary of the sequence kind
    private static string TypeOf(EmbossSeqFeature gf) => EmbossFeatureVocab.TypeExternal((gf.Protein ? EmbossFeatureVocab.Protein : EmbossFeatureVocab.Dna).Types, gf.Type);

    // showfeat_WriteFeat
    private static void WriteFeat(Line line, char strand, int fstart, int fend, int width, int beg, int end)
    {
        int len = end - beg + 1;
        float tf = ((float)(fstart - beg) / len) * width - 1f;
        int pos1 = (int)tf;
        tf = ((float)(fend - beg) / len) * width - 1f;
        int pos2 = (int)tf;
        var t = line.Text;
        if (t.Length == 0)
        {
            // the line of a collapsed feature that has no line yet is an empty string: ajStrPasteCountK on it writes to the end of the address space (the reference program dies by a fault)
            if ((pos1 >= 0 && pos1 < width) || (pos2 >= 0 && pos2 < width))
                throw new InvalidOperationException("showfeat -collapse draws on a line that does not exist: the reference program writes beyond an empty string and dies by a fault.");
            return;
        }

        for (int i = pos1; i < pos2; i++)
        {
            if (i >= 0 && i < width && t[i] == ' ')
                t[i] = '-';
        }

        if (pos1 >= 0 && pos1 < width)
            t[pos1] = strand == '-' ? '<' : '|';
        if (pos2 >= 0 && pos2 < width)
            t[pos2] = strand == '+' ? '>' : '|';
    }

    // showfeat_FeatOut
    private static void FeatOut(StringBuilder o, Line l, EmbossShowfeatRequest r)
    {
        bool strand = r.Strand, source = r.SourceShow, type = r.TypeShow, tags = r.TagShow, position = r.Position;
        if (r.Width > 3)
        {
            o.Append(l.Text);
            if (strand || source || type || tags || position)
                o.Append(' ');
        }

        if (strand)
        {
            o.Append(l.Strand == '+' ? '+' : l.Strand == '-' ? '-' : '0');
            if (source || type || tags || position)
                o.Append(' ');
        }

        if (source)
        {
            o.Append(l.Source);
            if (type || tags || position)
                o.Append(' ');
        }

        if (position)
        {
            o.Append(l.Pos);
            if (type || tags)
                o.Append(' ');
        }

        if (type)
        {
            o.Append(l.Type);
            if (tags)
                o.Append(' ');
        }

        if (tags)
            o.Append(l.Tags);
        o.Append('\n');
    }

    // showfeat_AddPos
    private static void AddPos(StringBuilder posout, int start, int end)
    {
        if (posout.Length > 0)
            posout.Append(',');
        posout.Append(CultureInfo.InvariantCulture, $"{start}-{end}");
    }

    // ajListToindex: a shell sort of the index with the compare function (not stable)
    private static void ToIndex(List<EmbossSeqFeature> list, int[] lind, Comparison<EmbossSeqFeature> compar)
    {
        int n = list.Count;
        var idx = new long[n];
        for (int i = 0; i < n; i++)
            idx[i] = i;
        for (int s = n / 2; s > 0; s /= 2)
        {
            for (int i = s; i < n; ++i)
            {
                for (int j = i - s; j >= 0 && compar(list[(int)idx[j]], list[(int)idx[j + s]]) > 0; j -= s)
                {
                    (lind[j], lind[j + s]) = (lind[j + s], lind[j]);
                    (idx[j], idx[j + s]) = (idx[j + s], idx[j]);
                }
            }
        }
    }

    private static int CmpCase(string a, string b) => string.CompareOrdinal(a.ToLowerInvariant(), b.ToLowerInvariant());

    private static int? StrandOrder(EmbossSeqFeature c, EmbossSeqFeature d)
    {
        char sc = StrandOf(c), sd = StrandOf(d);
        if (sc == sd)
            return null;
        if (sc == '+')
            return -1;
        if (sc == '\0' && sd == '-')
            return -1;
        return 1;
    }

    // showfeat_CompareFeatSource
    private static int CompareSource(EmbossSeqFeature c, EmbossSeqFeature d)
    {
        if (StrandOrder(c, d) is { } o)
            return o;
        int val = CmpCase(c.Source, d.Source);
        if (val != 0)
            return val;
        val = CmpCase(TypeOf(c), TypeOf(d));
        if (val != 0)
            return val;
        return c.GetStart() - d.GetStart();
    }

    // showfeat_CompareFeatType
    private static int CompareType(EmbossSeqFeature c, EmbossSeqFeature d)
    {
        if (StrandOrder(c, d) is { } o)
            return o;
        int val = CmpCase(TypeOf(c), TypeOf(d));
        if (val != 0)
            return val;
        val = CmpCase(c.Source, d.Source);
        if (val != 0)
            return val;
        return c.GetStart() - d.GetStart();
    }

    // showfeat_CompareFeatPos
    private static int ComparePos(EmbossSeqFeature c, EmbossSeqFeature d)
    {
        if (StrandOrder(c, d) is { } o)
            return o;
        int val = c.GetStart() - d.GetStart();
        if (val != 0)
            return val;
        val = CmpCase(TypeOf(c), TypeOf(d));
        if (val != 0)
            return val;
        return CmpCase(c.Source, d.Source);
    }

    // embMiscMatchPatternDelimC: the text matches (any case, wildcards) one of the patterns of the list (separated by white space and the delimiters)
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

    // showfeat_MatchPatternTags
    private static bool MatchPatternTags(EmbossSeqFeature feat, string tpattern, string vpattern, bool stricttags, StringBuilder tagstmp, bool values)
    {
        bool val = feat.Tags.Count == 0 && tpattern == "*" && vpattern == "*";
        foreach (var (tagnam, tagval) in feat.Tags)
        {
            bool tval = MatchPatternDelim(tagnam, tpattern, ",;|");
            bool vval = tagval.Length == 0 ? vpattern == "*" : tagval == vpattern || MatchPatternDelim(tagval, vpattern, ",;|");
            if (tval && vval)
            {
                val = true;
                if (tpattern != "*" || tagnam != "translation")
                {
                    if (values)
                        tagstmp.Append(' ').Append(tagnam).Append("=\"").Append(tagval).Append('"');
                    else
                        tagstmp.Append(' ').Append(tagnam);
                }
            }
            else if (!stricttags && tagnam != "translation")
            {
                if (values)
                    tagstmp.Append(' ').Append(tagnam).Append("=\"").Append(tagval).Append('"');
                else
                    tagstmp.Append(' ').Append(tagnam);
            }
        }

        if (!val)
            tagstmp.Clear();
        return val;
    }
}
