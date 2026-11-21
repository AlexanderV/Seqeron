using System.Text;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>The flags of a feature location (<c>AJFEATFLAG_*</c> of <c>ajfeat.h</c>; the values are private to this port).</summary>
[Flags]
internal enum EmbossLocationBits
{
    None = 0,
    StartBeforeSeq = 1,
    EndAfterSeq = 2,
    Point = 4,
    BetweenSeq = 8,
    RemoteId = 16,
    Label = 32,
    StartTwo = 64,
    EndTwo = 128,
    Multiple = 256,
    Order = 512,
    ComplementMain = 1024,
    StartUnsure = 2048,
    EndUnsure = 4096,
}

/// <summary>The state the feature library keeps between the features of ONE program run: <c>featRegSpecialRptunitSeq</c> is a static shared by the checks of <c>/rpt_unit</c> and <c>/rpt_unit_seq</c>, each of which compiles ITS pattern into it when it is still empty (the first of the two to run decides for both).</summary>
internal sealed class EmbossFeatureReadContext
{
    public System.Text.RegularExpressions.Regex? RptUnitSeq;
}

/// <summary>A feature of a feature table as the EMBOSS library holds it (<c>AjSFeature</c>): the (external) type, the location, the strand and the tags in the order they were added; a join / order location is one feature with its exons as <see cref="Subfeatures"/>.</summary>
internal sealed class EmbossSeqFeature
{
    public string Type = "";
    public int Start, End, Start2, End2, Frame, Exon;
    public EmbossLocationBits Flags;
    public char Strand;
    public List<(string Name, string Value)> Tags = [];
    public List<EmbossSeqFeature>? Subfeatures;

    /// <summary><c>ajFeatGetStart</c>: the second position of an uncertain start when it is the lower.</summary>
    public int GetStart() => Flags.HasFlag(EmbossLocationBits.StartTwo) && Start2 < Start ? Start2 : Start;

    /// <summary><c>ajFeatGetEnd</c>: the second position of an uncertain end when it is the lower.</summary>
    public int GetEnd() => Flags.HasFlag(EmbossLocationBits.EndTwo) && End2 < End ? End2 : End;
}

/// <summary>
/// A nucleotide feature table read from the features of a GenBank / EMBL entry (<c>ajfeatread.c</c> <c>featReadEmbl</c>, <c>featEmblFromLine</c>, <c>featEmblProcess</c> and the location / tag readers; <c>ajfeat.c</c> <c>ajFeatNewNucFlags</c>,
/// <c>ajFeatTagAddSS</c>, <c>ajFeattableTrimOff</c>, <c>ajFeattableReverse</c>, <c>ajFeatReverse</c>, <c>ajFeatTrimOffRange</c>; EMBOSS 6.6.0, kimrutherford/EMBOSS at raw.githubusercontent.com). A5-8c3b, B09 F218.
/// </summary>
internal sealed class EmbossSeqFeatureTable
{
    /// <summary>The features in table order (the parents of the join / order locations).</summary>
    public List<EmbossSeqFeature> Features = [];

    /// <summary>The length of the table (<c>ajFeattableSetLength</c>: the length of the sequence once it is read; while reading, the largest position seen).</summary>
    public int Length;

    /// <summary>The table was reversed (<c>thys-&gt;Reversed</c>).</summary>
    public bool Reversed;

    /// <summary>The entry name (<c>thys-&gt;Seqid</c>).</summary>
    public string SeqId = "";

    // ---------------------------------------------------------------------------------------------------------------- reading

    /// <summary>
    /// <c>featReadEmbl</c> over the lines of ONE entry: the <c>FT</c> lines of an EMBL entry, the feature lines (five blanks) after <c>FEATURES</c> of a GenBank entry. Reading stops at <c>//</c>.
    /// </summary>
    /// <param name="context">The state of the feature library shared by the entries of one program run (<see cref="EmbossFeatureReadContext"/>); null = a fresh one.</param>
    public static EmbossSeqFeatureTable Read(IEnumerable<string> entryLines, EmbossFeatureReadContext? context = null)
    {
        ArgumentNullException.ThrowIfNull(entryLines);
        var table = new EmbossSeqFeatureTable { _context = context ?? new EmbossFeatureReadContext() };
        string saveFeature = "", saveLoc = "", saveLine = "";
        bool isGenbank = false;
        foreach (string raw in entryLines)
        {
            string line = raw;
            if (line.StartsWith("//", StringComparison.Ordinal))
                break;
            if (line.StartsWith("ID   ", StringComparison.Ordinal))
            {
                isGenbank = false;
                string[] tokens = line.Split([' ', ';', '\t', '\n', '\r'], StringSplitOptions.RemoveEmptyEntries);
                table.SeqId = tokens.Length > 1 ? tokens[1] : "";
            }
            else if (line.StartsWith("LOCUS ", StringComparison.Ordinal) || line.StartsWith("FEATURES  ", StringComparison.Ordinal))
            {
                isGenbank = true;
            }
            else if (line.StartsWith("FT   ", StringComparison.Ordinal))
            {
                table.FromLine(line.TrimEnd(' ', '\t', '\n', '\r'), ref saveFeature, ref saveLoc, ref saveLine);
            }
            else if (isGenbank && line.StartsWith("     ", StringComparison.Ordinal))
            {
                table.FromLine(line.TrimEnd(' ', '\t', '\n', '\r'), ref saveFeature, ref saveLoc, ref saveLine);
            }
            else if (isGenbank)
            {
                isGenbank = false;
            }
        }

        table.FromLine(null, ref saveFeature, ref saveLoc, ref saveLine);
        return table;
    }

    // featEmblFromLine
    private void FromLine(string? origLine, ref string saveFeature, ref string saveLoc, ref string saveLine)
    {
        string line = "";
        bool newFeature = false, doFeature = false;
        if (origLine is not null)
        {
            line = origLine.Length >= 5 ? origLine[5..] : "";
            if (line.Length > 0 && line[0] != ' ')
            {
                newFeature = true;
                if (saveLoc.Length > 0)
                    doFeature = true;
            }
        }
        else if (saveLoc.Length > 0)
        {
            doFeature = true;
        }

        if (doFeature)
        {
            Process(saveFeature, ref saveLoc, ref saveLine);
            saveLoc = "";
            saveLine = "";
        }

        if (origLine is null)
            return;
        line = RemoveWhiteExcess(line);
        if (newFeature)
        {
            // ajStrTokenAssignC(" "), the key, then the location up to a blank or '/', then the rest of the string
            int p = 0;
            string key = NextToken(line, ref p, " ");
            saveFeature = key;
            string loc = NextToken(line, ref p, " /");
            if (loc.Length > 0 || p < line.Length)
            {
                saveLoc = loc;
                saveLine = p < line.Length ? line[p..] : "";
            }
            else
            {
                saveLoc = "";
                saveLine = "";                                              // location only
            }

            return;
        }

        if (saveLine.Length == 0 && (line.Length == 0 || line[0] != '/'))
        {
            // no tag-values yet: more of the location
            int p = 0;
            string temp = NextToken(line, ref p, " ");
            if (temp.Length > 0)
                saveLoc += temp;
            string rest = p < line.Length ? line[p..] : "";
            if (rest.Length > 0)
                saveLine += rest;
            return;
        }

        saveLine += " " + line;                                              // tag-values continued
    }

    // ajStrTokenNextParseC: skips the leading delimiters, takes the characters up to the next delimiter, and steps over ONE delimiter
    private static string NextToken(string s, ref int p, string delimiters)
    {
        while (p < s.Length && delimiters.Contains(s[p], StringComparison.Ordinal))
            p++;
        int start = p;
        while (p < s.Length && !delimiters.Contains(s[p], StringComparison.Ordinal))
            p++;
        string token = s[start..p];
        if (p < s.Length)
            p++;
        return token;
    }

    // ajStrRemoveWhiteExcess
    private static string RemoveWhiteExcess(string s)
    {
        var sb = new StringBuilder();
        bool space = false;
        foreach (char c in s)
        {
            char x = c is '\t' or '\n' or '\r' ? ' ' : c;
            if (x == ' ')
            {
                space = sb.Length > 0;
                continue;
            }

            if (space)
                sb.Append(' ');
            space = false;
            sb.Append(x);
        }

        return sb.ToString();
    }

    private static string RemoveWhite(string s) => string.Concat(s.Where(c => !char.IsWhiteSpace(c)));

    // featEmblProcess
    private void Process(string feature, ref string loc, ref string tags)
    {
        EmbossSeqFeature? ret = null, gf = null;
        bool fwd = true, join = false, order = false, parent = true;
        char strand;
        int frame = 0;
        int exon = 0;
        EmbossLocationBits noParentFlags = ~(EmbossLocationBits.Point | EmbossLocationBits.BetweenSeq | EmbossLocationBits.RemoteId | EmbossLocationBits.Label
            | EmbossLocationBits.StartTwo | EmbossLocationBits.EndTwo | EmbossLocationBits.StartBeforeSeq | EmbossLocationBits.EndAfterSeq);
        loc = RemoveWhite(loc);
        tags = RemoveWhiteExcess(tags);
        string opval = loc;
        if (opval.Contains(',', StringComparison.Ordinal))
        {
            while (opval.Length > 0 && OperOut(opval, out string opnam, out string inner))
            {
                if (!HasParentheses(inner))
                    break;
                if (opnam.Equals("complement", StringComparison.OrdinalIgnoreCase))
                    fwd = !fwd;
                else if (opnam.Equals("one_of", StringComparison.OrdinalIgnoreCase))
                    order = true;
                else if (opnam.Equals("join", StringComparison.OrdinalIgnoreCase))
                    join = true;
                else if (opnam.Equals("order", StringComparison.OrdinalIgnoreCase) || opnam.Equals("group", StringComparison.OrdinalIgnoreCase))
                    order = true;
                opval = inner;
            }
        }

        while (opval.Length > 0)
        {
            bool locFwd = fwd, begBound = false, endBound = false, simple = false, between = false, hasOper = false, remoteId = false, isLabel = false;
            uint begNum = 0, endNum = 0, beg2 = 0, end2 = 0;
            string featId = "", featLabel = "", locstr, rest;
            if (OperIn(opval, out string op2, out string inLoc, out string after))
            {
                if (op2.Equals("complement", StringComparison.OrdinalIgnoreCase))
                    locFwd = !locFwd;
                locstr = inLoc;
                opval = after;
                hasOper = true;
            }
            else
            {
                locstr = opval;
            }

            OperNone(locstr, out featId, out string tmp, out rest);
            if (featId.Length > 0)
                remoteId = true;
            string begstr, endstr;
            if (!Loc(tmp, out begstr, out between, out simple, out endstr))
            {
                begstr = tmp;
                endstr = begstr;
                simple = true;
            }

            tmp = rest;
            if (!hasOper)
                opval = tmp;
            if (LocNum(begstr, out begBound, out begNum))
            {
                beg2 = 0;
            }
            else if (LocRange(begstr, out begNum, out beg2))
            {
                begBound = false;
            }
            else
            {
                isLabel = true;
                simple = true;
                featLabel = begstr;
            }

            if (LocNum(endstr, out endBound, out endNum))
            {
                end2 = 0;
            }
            else if (LocRange(endstr, out end2, out endNum))
            {
                endBound = false;
            }
            else
            {
                isLabel = true;
                simple = true;
                featLabel = endstr;
            }

            var flags = EmbossLocationBits.None;
            strand = locFwd ? '+' : '-';
            if (simple) flags |= EmbossLocationBits.Point;
            if (between) flags |= EmbossLocationBits.BetweenSeq;
            if (end2 != 0) flags |= EmbossLocationBits.EndTwo;
            if (beg2 != 0) flags |= EmbossLocationBits.StartTwo;
            if (begBound) flags |= EmbossLocationBits.StartBeforeSeq;
            if (endBound) flags |= EmbossLocationBits.EndAfterSeq;
            if (remoteId) flags |= EmbossLocationBits.RemoteId;
            if (isLabel) flags |= EmbossLocationBits.Label;
            if (join || order) flags |= EmbossLocationBits.Multiple;
            if (order) flags |= EmbossLocationBits.Order;
            if (parent && !fwd) flags |= EmbossLocationBits.ComplementMain;
            string type = EmbossFeatureVocabulary.Type(feature);
            uint tmpBeg, tmpEnd;
            EmbossLocationBits tmpFlags;
            string tmpFeatId;
            if (parent)
            {
                if (opval.Length > 0 && remoteId)
                {
                    tmpBeg = 0;
                    tmpEnd = 0;
                    tmpFlags = flags & noParentFlags;
                    tmpFeatId = "";
                }
                else
                {
                    tmpBeg = begNum;
                    tmpEnd = endNum;
                    tmpFlags = flags;
                    tmpFeatId = featId;
                }

                gf = NewNuc(type, (int)tmpBeg, (int)tmpEnd, strand, frame, exon, (int)beg2, (int)end2, tmpFeatId, featLabel, tmpFlags);
            }

            if (!parent || opval.Length > 0)
            {
                exon++;
                var sub = new EmbossSeqFeature { Type = type, Strand = strand, Frame = frame, Start = (int)begNum, End = (int)endNum, Start2 = (int)beg2, End2 = (int)end2, Exon = exon, Flags = flags };
                gf!.Subfeatures ??= [];
                gf.Subfeatures.Add(sub);
                if (!remoteId)
                {
                    Length = Math.Max(Length, sub.Start);
                    Length = Math.Max(Length, sub.End);
                    if (fwd)
                    {
                        if (gf.Start == 0 || gf.Start > (int)begNum)
                            gf.Start = (int)begNum;
                        if (gf.End == 0 || gf.End < (int)endNum)
                            gf.End = (int)endNum;
                    }
                    else
                    {
                        if (gf.End == 0 || gf.End > (int)begNum)
                            gf.End = (int)begNum;
                        if (gf.Start == 0 || gf.Start < (int)endNum)
                            gf.Start = (int)endNum;
                    }
                }
            }

            if (parent)
            {
                ret = gf;
                parent = false;
                exon = 1;
            }
        }

        while (tags.Length > 0)
        {
            if (TvTagVal(ref tags, out string tag, out string val))
            {
                val = QuoteStrip(val);
                if (tag == "codon_start")
                {
                    frame = ToInt(val);
                    if (gf is not null)
                        gf.Frame = frame;
                }

                if (tag == "note" && val.Length > 0 && val[0] == '*')
                {
                    val = val[1..];
                    int sp = val.IndexOf(' ', StringComparison.Ordinal);
                    string tagName = sp < 0 ? val : val[..sp];
                    string tagRest = sp < 0 ? "" : val[(sp + 1)..].TrimStart();
                    if (tagName == "Type")
                    {
                        int sp2 = tagRest.IndexOf(' ', StringComparison.Ordinal);
                        if (gf is not null)
                            gf.Type = sp2 < 0 ? tagRest : tagRest[..sp2];
                    }
                    else if (ret is not null)
                    {
                        foreach (string piece in tagRest.Split(','))
                        {
                            if (!Gff3PredefinedTags.Contains(tagName))
                                ret.Tags.Add((tagName, piece));
                        }
                    }

                    continue;
                }

                if (ret is not null)
                    ret.Tags.Add(TagValue(tag, val, _context));
            }
            else if (TvRest(ref tags, out _))
            {
                // anything non-white up to '/' is bad: skipped
            }
            else
            {
                tags = "";
            }
        }
    }

    private static readonly HashSet<string> Gff3PredefinedTags = ["ID", "Name", "Alias", "Parent", "Target", "Gap", "Derives_from", "Note", "Dbxref", "Ontology_term"];

    private static int ToInt(string s)
    {
        // ajStrToInt: sscanf("%d")-like on the leading integer
        int i = 0;
        while (i < s.Length && char.IsWhiteSpace(s[i]))
            i++;
        int start = i;
        if (i < s.Length && s[i] is '+' or '-')
            i++;
        while (i < s.Length && char.IsAsciiDigit(s[i]))
            i++;
        return int.TryParse(s.AsSpan(start, i - start), out int v) ? v : 0;
    }

    // ajFeatNewNucFlags
    private EmbossSeqFeature NewNuc(string type, int start, int end, char strand, int frame, int exon, int start2, int end2, string entryId, string label, EmbossLocationBits flags)
    {
        var f = new EmbossSeqFeature { Type = type, Strand = strand is '+' or '-' ? strand : '\0', Frame = frame, Exon = exon, Start = start, End = end, Start2 = start2, End2 = end2, Flags = flags };
        if (entryId.Length == 0 && !flags.HasFlag(EmbossLocationBits.RemoteId) && !flags.HasFlag(EmbossLocationBits.Label))
        {
            Length = Math.Max(Length, f.Start);
            Length = Math.Max(Length, f.End);
        }

        Features.Add(f);
        return f;
    }

    // ajStrHasParentheses: the parentheses are balanced (true when there are none)
    private static bool HasParentheses(string s)
    {
        int left = 0;
        foreach (char c in s)
        {
            if (c == '(')
            {
                left++;
            }
            else if (c == ')')
            {
                if (left == 0)
                    return false;
                left--;
            }
        }

        return left == 0;
    }

    // featEmblOperOut: "oper(content)" covering the whole string; the operator starts with a lower-case letter
    private static bool OperOut(string loc, out string opnam, out string value)
    {
        opnam = value = "";
        int left = 0, right = 0, bracepos = 0;
        if (loc.Length == 0 || loc[0] == '(')
            return false;
        if (!char.IsAsciiLetter(loc[0]) || !char.IsAsciiLetterLower(loc[0]))
            return false;
        int ipos = 0;
        foreach (char c in loc)
        {
            ipos++;
            switch (c)
            {
                case ',':
                    if (left == 0)
                        return false;
                    break;
                case '(':
                    if (left == 0)
                        bracepos = ipos - 1;
                    left++;
                    break;
                case ')':
                    right++;
                    if (right == left)
                    {
                        if (ipos != loc.Length)
                            return false;
                    }
                    else if (right > left)
                    {
                        return false;
                    }

                    break;
            }
        }

        if (left == 0 || right != left)
            return false;
        opnam = loc[..bracepos];
        value = loc.Substring(bracepos + 1, ipos - 2 - bracepos);
        return true;
    }

    // featEmblOperIn: the first operator and its parenthesised content up to the first comma outside it; the rest after the comma
    private static bool OperIn(string loc, out string opnam, out string value, out string rest)
    {
        opnam = value = rest = "";
        int left = 0, right = 0, ipos = 0, bracepos = 0;
        if (loc.Length > 0 && loc[0] == '(')
            return false;
        foreach (char c in loc)
        {
            if (c == ',')
                break;
            ipos++;
            switch (c)
            {
                case '(':
                    if (left == 0)
                        bracepos = ipos - 1;
                    left++;
                    if (left > 3)
                        return false;
                    if (left - right > 2)
                        return false;
                    break;
                case ')':
                    right++;
                    if (right > 3 || right > left)
                        return false;
                    break;
            }
        }

        if (left == 0 || right != left)
            return false;
        opnam = loc[..bracepos];
        value = loc.Substring(bracepos + 1, ipos - 2 - bracepos);
        rest = ipos < loc.Length ? loc[(ipos + 1)..] : "";
        return true;
    }

    // featEmblOperNone: [entryid:]location up to the first comma; the rest after it
    private static void OperNone(string loc, out string entryId, out string value, out string rest)
    {
        int ipos = 0, idpos = 0;
        foreach (char c in loc)
        {
            if (c == ',')
                break;
            ipos++;
            if (c == ':')
                idpos = ipos;
        }

        if (idpos > 0)
        {
            entryId = loc[..(idpos - 1)];
            value = loc[idpos..ipos];
        }
        else
        {
            entryId = "";
            value = loc[..ipos];
        }

        rest = ipos < loc.Length ? loc[(ipos + 1)..] : "";
    }

    // featEmblLoc: "beg..end", "beg^end" (between), a single position
    private static bool Loc(string loc, out string begstr, out bool between, out bool simple, out string endstr)
    {
        begstr = endstr = "";
        between = simple = false;
        int ipos = 0, bracket = 0, iend = 0, ibeg = 0;
        bool numeric = true, dot = false, end = false;
        foreach (char c in loc)
        {
            ipos++;
            switch (c)
            {
                case '.':
                    if (bracket > 0)
                        break;
                    if (dot)
                    {
                        if (end)
                            return false;
                        end = true;
                        ibeg = ipos - 2;
                        iend = ipos - 1;
                    }
                    else
                    {
                        dot = true;
                    }

                    break;
                case '^':
                    ibeg = ipos - 1;
                    iend = ipos - 1;
                    end = true;
                    between = true;
                    break;
                case '(':
                    bracket++;
                    numeric = false;
                    break;
                case ')':
                    if (bracket == 0)
                        return false;
                    bracket--;
                    break;
                default:
                    if (!char.IsAsciiDigit(c))
                        numeric = false;
                    break;
            }
        }

        if (dot && !end)
            return false;
        if (end)
        {
            begstr = Sub(loc, 0, ibeg - 1);
            endstr = Sub(loc, iend + 1, ipos - 1);
        }
        else
        {
            simple = true;
            begstr = loc;
            endstr = numeric ? loc : "";
        }

        return true;
    }

    // ajStrAssignSubS(from, to): inclusive 0-based positions, the end clipped, an empty result when it is before the start
    private static string Sub(string s, int first, int last)
    {
        if (last < 0)
            last = s.Length + last;
        if (first < 0)
            first = s.Length + first;
        if (first < 0)
            first = 0;
        if (last >= s.Length)
            last = s.Length - 1;
        if (last < first)
            return "";
        return s.Substring(first, last - first + 1);
    }

    // featEmblLocNum: [<>]digits
    private static bool LocNum(string loc, out bool bound, out uint num)
    {
        bound = false;
        num = 0;
        int i = 0;
        if (loc.Length > 0 && loc[0] is '<' or '>')
        {
            bound = true;
            i++;
        }

        for (; i < loc.Length; i++)
        {
            if (!char.IsAsciiDigit(loc[i]))
                return false;
            num = unchecked((10 * num) + (uint)(loc[i] - '0'));
        }

        return true;
    }

    // featEmblLocRange: (num.num)
    private static bool LocRange(string loc, out uint num1, out uint num2)
    {
        num1 = num2 = 0;
        if (loc.Length == 0 || loc[0] != '(' || loc[^1] != ')')
            return false;
        bool dot = false;
        int ipos = 0;
        foreach (char c in loc)
        {
            ipos++;
            switch (c)
            {
                case '.':
                    if (dot)
                        return false;
                    dot = true;
                    break;
                case '(':
                    if (ipos != 1)
                        return false;
                    break;
                case ')':
                    if (ipos != loc.Length)
                        return false;
                    break;
                default:
                    if (!char.IsAsciiDigit(c))
                        return false;
                    if (dot)
                        num2 = unchecked((10 * num2) + (uint)(c - '0'));
                    else
                        num1 = unchecked((10 * num1) + (uint)(c - '0'));
                    break;
            }
        }

        return dot;
    }

    // ajStrQuoteStrip
    private static string QuoteStrip(string v)
    {
        if (v.Length > 0 && v[^1] == '"')
            v = v[..^1];
        if (v.Length > 0 && v[0] == '"')
            v = v[1..];
        v = v.Replace("\"\"", "\"", StringComparison.Ordinal);
        return v.Replace("\\\"", "\"", StringComparison.Ordinal);
    }

    private int _tvDepth;

    private EmbossFeatureReadContext _context = new();

    // featEmblTvTagVal: the next "/name=value" of the tag text
    private bool TvTagVal(ref string tags, out string name, out string value)
    {
        name = value = "";
        int cp = 0;
        string t = tags;
        while (cp < t.Length && t[cp] == ' ')
            cp++;
        if (cp >= t.Length || t[cp] != '/')
            return false;
        cp++;
        int cq = cp;
        while (cp < t.Length && (char.IsAsciiLetter(t[cp]) || t[cp] == '_'))
            cp++;
        name = t[cq..cp];
        while (cp < t.Length && t[cp] == ' ')
            cp++;
        if (name.Length == 0)
            return false;
        if (cp >= t.Length)
        {
            tags = "";
            return true;
        }

        if (t[cp] == '/')
        {
            tags = t[cp..];
            return true;
        }

        if (t[cp] != '=')
            return false;
        cp++;
        if (cp < t.Length && t[cp] == '"')
        {
            cq = cp;
            cp++;
            int i = 0;
            var sb = new StringBuilder();
            while (cp < t.Length)
            {
                i++;
                if (t[cp] == '"')
                {
                    sb.Append(t, cq, i);
                    cp++;
                    if (cp >= t.Length || t[cp] != '"')
                    {
                        while (cp < t.Length && t[cp] == ' ')
                            cp++;
                        tags = cp < t.Length ? t[cp..] : "";
                        value = sb.ToString();
                        return true;
                    }

                    // "" - is it really internal?
                    if (_tvDepth > 0)
                    {
                        value = sb.ToString();
                        return true;
                    }

                    _tvDepth++;
                    string testTags = t[cp..];
                    if (TvTagVal(ref testTags, out _, out _))
                    {
                        _tvDepth--;
                        cp++;
                        while (cp < t.Length && t[cp] == ' ')
                            cp++;
                        tags = cp < t.Length ? t[cp..] : "";
                        value = sb.ToString();
                        return true;
                    }

                    _tvDepth--;
                    cq = cp - 1;
                    i = 1;
                }

                cp++;
            }

            value = sb.ToString();
            tags = "";
            return false;
        }

        cq = cp;
        while (cp < t.Length && t[cp] != ' ' && t[cp] != '/')
            cp++;
        value = t[cq..cp];
        tags = t[cp..];
        return true;
    }

    // featEmblTvRest: skip to the next "/name" that reads as a tag
    private bool TvRest(ref string tags, out string skip)
    {
        var sb = new StringBuilder();
        for (int cp = 0; cp < tags.Length; cp++)
        {
            if (tags[cp] == '/')
            {
                string test = tags[cp..];
                if (TvTagVal(ref test, out _, out _))
                {
                    tags = tags[cp..];
                    skip = sb.ToString();
                    return true;
                }
            }

            sb.Append(tags[cp]);
        }

        tags = "";
        skip = sb.ToString();
        return true;
    }

    // ajFeatTagAddSS (DNA table): the canonical tag name and the value as the library stores them
    private static (string Name, string Value) TagValue(string tag, string value, EmbossFeatureReadContext context)
    {
        var (name, kind, known) = EmbossFeatureVocabulary.Tag(tag);
        string v = value;
        switch (kind)
        {
            case "QSPECIAL" or "SPECIAL":
                if (!SpecialValue(name, ref v, context))
                    return (EmbossFeatureVocabulary.DefaultTag, "*" + name + ": " + value);
                break;
            case "QTEXT":
                if (!known)
                    v = "*" + tag + ": " + value;
                break;
        }

        return (name, v);
    }

    private const System.Text.RegularExpressions.RegexOptions Rx = System.Text.RegularExpressions.RegexOptions.CultureInvariant;

    private static bool Matches(string pattern, string value) => System.Text.RegularExpressions.Regex.IsMatch(value, pattern, Rx);

    // ajFeattagSpecial: the checks (and the corrections) of the special tags; false = the value is bad and the tag becomes a note
    private static bool SpecialValue(string tag, ref string value, EmbossFeatureReadContext context)
    {
        switch (tag)
        {
            case "anticodon":
                return Matches("^[(]pos:([0-9]+)[.][.]([0-9]+),aa:([^)]+)[)]$", value);
            case "bio_material" or "culture_collection" or "specimen_voucher":
                return Matches("^([^:]+)(:([^:]+))?(:(.*))?$", value);
            case "citation":
                {
                    int cp = 0, i = 0;
                    bool ret = false;
                    if (value.Length > 0 && value[cp++] == '[')
                    {
                        while (cp < value.Length && char.IsAsciiDigit(value[cp]))
                        {
                            i++;
                            cp++;
                        }

                        if (cp < value.Length && value[cp++] == ']' && cp >= value.Length)
                            ret = true;
                        if (i == 0)
                            ret = false;
                    }

                    return ret;
                }

            case "codon":
                {
                    var m = System.Text.RegularExpressions.Regex.Match(value, "^[(]seq:\"([acgt][acgt][acgt])\",aa:([^)]+)[)]$", Rx);
                    if (m.Success)
                        return true;
                    var bad = System.Text.RegularExpressions.Regex.Match(value, "^[(]seq:([acgt][acgt][acgt]),aa:([^)]+)[)]$", Rx);
                    if (bad.Success)
                    {
                        value = "(seq:\"" + bad.Groups[1].Value + "\",aa:" + bad.Groups[2].Value + ")";
                        return true;
                    }

                    return false;
                }

            case "collection_date":
                return CollectionDate(value);
            case "compare":
                return Matches("^([A-Z]+[0-9]+([.][1-9][0-9]*)?)$", value);
            case "cons_splice":
                return ConsSplice(value);
            case "db_xref":
                {
                    int colon = value.IndexOf(':', StringComparison.Ordinal);
                    return colon >= 0 && colon + 1 < value.Length;
                }

            case "estimated_length":
                return Matches("^(unknown|[1-9][0-9]*)$", value);
            case "inference":
                return Inference(value);
            case "lat_lon":
                return Matches("^[0-9.]+ [NS] [0-9.]+ [EW]$", value);
            case "mobile_element":
                {
                    var m = System.Text.RegularExpressions.Regex.Match(value, "^([^:]+)(:(.*))?$", Rx);
                    return m.Success && Array.IndexOf(["transposon", "retrotransposon", "integron", "insertion sequence", "non-LTR retrotransposon", "SINE", "MITE", "LINE", "other"], m.Groups[1].Value) >= 0;
                }

            case "ncRNA_class":
                return Array.IndexOf(["antisense_RNA", "autocatalytically_spliced_intron", "ribozyme", "hammerhead_ribozyme", "RNase_P_RNA", "RNase_MRP_RNA", "telomerase_RNA", "guide_RNA", "rasiRNA",
                    "scRNA", "siRNA", "miRNA", "piRNA", "snoRNA", "snRNA", "SRP_RNA", "vault_RNA", "Y_RNA", "other"], value) >= 0;
            case "PCR_primers":
                return Matches("^((fwd|rev)_name: [^,]+, )?(fwd|rev)_seq: [acgtubdhkmnrsvwxy<>i]+(, ((fwd|rev)_name: [^,]+, )?(fwd|rev)_seq: [acgtubdhkmnrsvwxy<>i]+)*", value);
            case "protein_id":
                return ProteinId(value);
            case "replace":
                value = RemoveWhite(value);
                return Matches("^([abcdghkmnrstuvwxyABCDGHKMNRSTUVWXY/-]*)$", value);
            case "rpt_unit":
                {
                    context.RptUnitSeq ??= new System.Text.RegularExpressions.Regex("^([abcdghkmnrstuvwxyABCDGHKMNRSTUVWXY0-9/();-]+)$", Rx);
                    return Matches("^(([a-zA-Z0-9_]+:)?[0-9]+)[.][.]([0-9]+)$", value) || Matches("^(([a-zA-Z0-9_]+:)?[a-zA-Z0-9_]+)$", value)
                        || Matches("^(complement[(]([a-zA-Z0-9_]+:)?[0-9]+)[.][.]([0-9]+)[)]$", value) || context.RptUnitSeq.IsMatch(value);
                }

            case "rpt_unit_range" or "tag_peptide":
                return Matches("^(([a-zA-Z0-9_]+:)?[0-9]+)[.][.]([0-9]+)$", value) || Matches("^(([a-zA-Z0-9_]+:)?[a-zA-Z0-9_]+)$", value)
                    || Matches("^(complement[(]([a-zA-Z0-9_]+:)?[0-9]+)[.][.]([0-9]+)[)]$", value);
            case "rpt_unit_seq":
                return RptUnitSeq(ref value, context);
            case "transl_except":
                return Matches("^[(]pos:([0-9]+)[.][.]([0-9]+),aa:([^)]+)[)]$", value) || Matches("^[(]pos:complement\\(([0-9]+)[.][.]([0-9]+)\\),aa:([^)]+)[)]$", value)
                    || Matches("^[(]pos:([0-9]+),aa:([^)]+)[)]$", value) || Matches("^[(]pos:complement\\(([0-9]+)\\),aa:([^)]+)[)]$", value);
            case "translation":
                {
                    string w = RemoveWhite(value);
                    if (w.Any(c => !char.IsAsciiLetterUpper(c)))
                        return false;
                    value = w;
                    return true;
                }

            default:
                return false;                                           // "Unrecognised special EMBL feature tag"
        }
    }

    // featTagSpecialAllCollectiondate
    private static bool CollectionDate(string val)
    {
        var m = System.Text.RegularExpressions.Regex.Match(val, "^(([0-9][0-9])[-])?(([ADFJMNOS][A-Za-z][A-Za-z])[-])?([12][7890][0-9][0-9])$", Rx);
        if (!m.Success)
            return false;
        bool ret = true;
        string day = m.Groups[2].Value, mon = m.Groups[4].Value;
        if (day.Length > 0)
        {
            if (mon.Length == 0)
                ret = false;
            int d = int.Parse(day, System.Globalization.CultureInfo.InvariantCulture);
            if (d < 1 || d > 31)
                ret = false;
        }

        if (mon.Length > 0 && Array.FindIndex(["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], x => x.Equals(mon, StringComparison.OrdinalIgnoreCase)) < 0)
            ret = false;
        return ret;
    }

    // featTagSpecialAllConssplice
    private static bool ConsSplice(string val)
    {
        bool ret = false;
        if (!val.StartsWith("(5'site:", StringComparison.Ordinal) || val[^1] != ')')
            return false;
        int cp = 8;
        int cq = cp;
        while (cp < val.Length && val[cp] != ',')
            cp++;
        string beg = val[cq..cp];
        if (string.CompareOrdinal(val, cp, ",3'site:", 0, 8) != 0)
            return false;
        cp += 8;
        cq = cp;
        while (cp < val.Length && val[cp] != ')')
            cp++;
        string end = val[cq..cp];
        if (cp < val.Length && val[cp] == ')' && cp + 1 == val.Length && beg.Length > 0 && end.Length > 0)
            ret = true;
        return ret && SpliceWord(beg) && SpliceWord(end);
    }

    private static bool SpliceWord(string w)
        => w.Length > 0 && w[0] switch
        {
            'Y' => w.Equals("yes", StringComparison.OrdinalIgnoreCase),
            'N' => w.Equals("no", StringComparison.OrdinalIgnoreCase),
            'A' => w.Equals("absent", StringComparison.OrdinalIgnoreCase),
            _ => false,
        };

    private static readonly string[] InferenceTypes =
    [
        "non-experimental evidence, no additional details recorded", "similar to sequence", "similar to AA sequence", "similar to AA sequence sequence", "similar to DNA sequence", "similar to RNA sequence",
        "similar to RNA sequence, mRNA", "similar to RNA sequence,mRNA", "similar to RNA sequence, EST", "similar to RNA sequence, other RNA", "profile", "nucleotide motif", "protein motif", "ab initio prediction", "alignment",
    ];

    // featTagSpecialAllInference
    private static bool Inference(string val)
    {
        var m = System.Text.RegularExpressions.Regex.Match(val, "^([^:]+)(:([^:]+)?)*$", Rx);
        if (!m.Success)
            return false;
        string typ = m.Groups[1].Value;
        foreach (string t in InferenceTypes)
        {
            if (!typ.StartsWith(t, StringComparison.Ordinal))
                continue;
            if (typ == t)
                return true;
            if (SameSpeciesSuffix(typ) == typ.Length - t.Length)
                return true;
        }

        return false;
    }

    // the suffix of an inference type the checks of featTagSpecialAllInference accept (the first of the three endings that fits decides)
    private static int SameSpeciesSuffix(string typ)
    {
        if (typ.EndsWith(" (same species)", StringComparison.Ordinal))
            return 15;
        if (typ.EndsWith("(same species)", StringComparison.Ordinal))
            return 14;
        return typ.EndsWith(" (same species) ", StringComparison.Ordinal) ? 16 : -1;
    }

    private static readonly string[] RptUnitSeqSpecials =
    [
        "a-rich", "at-rich", "c-rich", "ct-rich", "CT-rich", "g-rich", "ga-rich", "GA-rich", "gc-rich", "t-rich", "CT_AT", "TA_CA", "AGAT_AGAC", "GAA-repeat", "GGA-repeat", "tgtc (imperfect)", "score", "direct", "dispersed", "inverted",
        "ca repeat", "cac/caa/cag compound repeat", "cata interrupted repeat", "CA dinucleotide microsatellite", "gaa/gac compound repeat", "dodecanucleotide", "heptamer", "gttt..gt", "telunit", "af(12)", "alu", "Alu", "AluJo", "AluSg",
        "AluY", "hif2DR", "B1element", "L1MC/D", "L1MC3", "L1MD3", "L1ME", "RMER17B", "RLTR11A", "poly A", "core-vntr1", "core-vntr2",
    ];

    // featTagSpecialAllRptunitseq
    private static bool RptUnitSeq(ref string value, EmbossFeatureReadContext context)
    {
        context.RptUnitSeq ??= new System.Text.RegularExpressions.Regex("^([ abcdghkmnrstuvwxyABCDGHKMNRSTUVWXY0-9_/\\[\\]\\(\\):;.,+-]+)$", Rx);
        if (context.RptUnitSeq.IsMatch(value))
        {
            value = RemoveWhite(value);
            return true;
        }

        if (Matches("^((complement[\\(])?[<]?[0-9]+[.][.]?[>]?[0-9]+[\\)]?)$", value))
            return true;
        return Array.IndexOf(RptUnitSeqSpecials, value) >= 0;
    }

    // featTagSpecialAllProteinid
    private static bool ProteinId(string val)
    {
        bool ret = false;
        if (val.StartsWith("ENSP", StringComparison.Ordinal))
        {
            int i = 0, cp = 4;
            while (cp < val.Length)
            {
                if (!char.IsAsciiDigit(val[cp++]))
                    break;
                i++;
            }

            ret = cp >= val.Length && i > 5 && val[4..].All(char.IsAsciiDigit);
            return ret;
        }

        int count = 0, idot = 0, p = 0;
        while (p < val.Length)
        {
            count++;
            char c = val[p];
            if (count <= 3)
            {
                if (c == '_')
                {
                    if (count != 3)
                        break;
                }
                else if (!char.IsAsciiLetterUpper(c))
                {
                    break;
                }
            }
            else if (c == '.')
            {
                if (idot > 0)
                    break;
                idot = count;
                if (idot < 9)
                    break;
            }
            else if (!char.IsAsciiDigit(c))
            {
                break;
            }

            p++;
        }

        return p >= val.Length && idot > 0 && count > idot;
    }

    // ---------------------------------------------------------------------------------------------------------------- trimming and reversing

    /// <summary><c>ajFeatReverse</c> for every feature that is not remote or a label (<c>ajFeattableReverse</c>) with the length <see cref="Length"/>.</summary>
    public void Reverse()
    {
        foreach (var f in Features)
        {
            if (f.Flags.HasFlag(EmbossLocationBits.RemoteId) || f.Flags.HasFlag(EmbossLocationBits.Label))
                continue;
            ReverseFeature(f, Length);
        }

        Reversed = !Reversed;
    }

    private static void ReverseFeature(EmbossSeqFeature t, int ilen)
    {
        var saved = t.Flags;
        if (t.Start != t.End)
            t.Strand = t.Strand == '-' ? '+' : '-';
        int itmp = t.Start;
        t.Start = t.End != 0 ? 1 + ilen - t.End : 0;
        t.End = itmp != 0 ? 1 + ilen - itmp : 0;
        itmp = t.Start2;
        t.Start2 = t.End2 != 0 ? 1 + ilen - t.End2 : 0;
        t.Flags = Swap(saved, t.Flags, EmbossLocationBits.StartBeforeSeq, EmbossLocationBits.EndAfterSeq);
        t.Flags = Swap(saved, t.Flags, EmbossLocationBits.StartTwo, EmbossLocationBits.EndTwo);
        t.Flags = Swap(saved, t.Flags, EmbossLocationBits.StartUnsure, EmbossLocationBits.EndUnsure);
        t.End2 = itmp != 0 ? 1 + ilen - itmp : 0;
        t.Frame = 0;
        if (t.Subfeatures is not null)
        {
            foreach (var s in t.Subfeatures)
                ReverseFeature(s, ilen);
        }
    }

    // "if(saveflags & A) flags |= B else flags &= ~B" and "if(saveflags & B) flags |= A else flags &= ~A"
    private static EmbossLocationBits Swap(EmbossLocationBits saved, EmbossLocationBits flags, EmbossLocationBits a, EmbossLocationBits b)
    {
        flags = saved.HasFlag(a) ? flags | b : flags & ~b;
        flags = saved.HasFlag(b) ? flags | a : flags & ~a;
        return flags;
    }

    /// <summary><c>ajFeattableTrimOff</c>: the features of a sequence that was cut to <paramref name="length"/> residues after <paramref name="offset"/> residues.</summary>
    public void TrimOff(int offset, int length)
    {
        int seqLength = length + offset;
        int begin = PositionII(Length, 1, 0);
        int end = Length;
        if (end > seqLength)
            end = seqLength;
        if (begin <= offset)
            begin = PositionII(Length, offset + 1, end);
        bool doBegin = begin > 1, doEnd = end < Length;
        Features.RemoveAll(f => !TrimOffRange(f, offset, begin, end, doBegin, doEnd));
    }

    // ajFeattablePosII
    private static int PositionII(int ilen, int imin, int ipos)
    {
        int jpos;
        if (ipos < 0)
            jpos = ilen + ipos + 1;
        else
            jpos = ipos != 0 ? ipos : 1;
        if (jpos > ilen)
            jpos = ilen;
        if (jpos < imin)
            jpos = imin;
        return jpos;
    }

    // ajFeatTrimOffRange: false deletes the feature
    private static bool TrimOffRange(EmbossSeqFeature ft, int offset, int begin, int end, bool doBegin, bool doEnd)
    {
        bool ok = true;
        int joffset = offset;
        if (doBegin && begin != 0)
            joffset += begin - 1;
        if (ft.Flags.HasFlag(EmbossLocationBits.RemoteId) || ft.Flags.HasFlag(EmbossLocationBits.Label))
            return true;
        if (doEnd)
        {
            if (ft.Start > end)
                return false;
            if (ft.Start2 > end)
                ft.Start2 = end;
            if (ft.End > end)
            {
                ft.End = end;
                ft.Flags |= EmbossLocationBits.EndAfterSeq;
            }

            if (ft.End2 > end)
                ft.End2 = end;
        }

        if (doBegin)
        {
            if (ft.End < begin)
                return false;
            if (ft.End2 != 0 && begin > offset + 1 && ft.End2 < begin)
                ft.End2 = begin;
            if (ft.Start != 0 && ft.Start < begin)
            {
                ft.Start = begin;
                ft.Flags |= EmbossLocationBits.StartBeforeSeq;
            }

            if (ft.Start2 != 0 && ft.Start2 < begin)
                ft.Start2 = begin;
        }

        if (joffset != 0)
        {
            if (ft.Start != 0)
                ft.Start -= joffset;
            if (ft.Start2 != 0)
                ft.Start2 -= joffset;
            if (ft.End != 0)
                ft.End -= joffset;
            if (ft.End2 != 0)
                ft.End2 -= joffset;
        }

        if (ft.Subfeatures is not null)
        {
            // every sub-feature is trimmed; one that lies outside the range makes the PARENT go (the loop only records the failure, it never removes the sub-feature)
            foreach (var sub in ft.Subfeatures)
            {
                if (!TrimOffRange(sub, offset, begin, end, doBegin, doEnd))
                    ok = false;
            }
        }

        return ok;
    }
}
