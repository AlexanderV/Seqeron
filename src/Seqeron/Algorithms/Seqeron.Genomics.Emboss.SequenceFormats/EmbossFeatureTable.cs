// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

/// <summary>The flags of a feature location (<c>AJFEATFLAG_*</c> of <c>ajfeatdata.h</c>, with the values of the library: the debug trace prints them in hexadecimal).</summary>
[Flags]
public enum EmbossLocationBits : uint
{
    None = 0,
    StartBeforeSeq = 0x0001,
    EndAfterSeq = 0x0002,
    Generated = 0x0004,
    BetweenSeq = 0x0008,
    StartTwo = 0x0010,
    EndTwo = 0x0020,
    Point = 0x0040,
    ComplementMain = 0x0080,
    Multiple = 0x0100,
    Group = 0x0200,
    Order = 0x0400,
    OneOf = 0x0800,
    RemoteId = 0x1000,
    Label = 0x2000,
    StartUnsure = 0x4000,
    EndUnsure = 0x8000,
}

/// <summary>The state the feature library keeps between the features and the tables of ONE program run: <c>featRegSpecialRptunitSeq</c> is a static shared by the checks of <c>/rpt_unit</c> and <c>/rpt_unit_seq</c>, each of which compiles ITS pattern into it when it is still empty (the first of the two to run decides for both); the statics of the GFF group readers and of the format readers (the feature the last call created) live here too.</summary>
/// <summary><c>ajFatal</c>: the program stops at once (the message is printed with the source line of the library).</summary>
public sealed class EmbossFatalException(string message) : Exception(message);

public sealed class EmbossFeatureReadContext
{
    public System.Text.RegularExpressions.Regex? RptUnitSeq;

    /// <summary>The messages the library printed (<c>ajWarn</c> / <c>ajErr</c>: <c>Warning: ...</c>, <c>Error: ...</c>) in the order they were printed.</summary>
    public List<string> Messages = [];

    /// <summary><c>ajFeatWarn</c> / <c>ajWarn</c>.</summary>
    public void Warn(string message) => Messages.Add("Warning: " + message);

    /// <summary><c>ajErr</c>.</summary>
    public void Error(string message) => Messages.Add("Error: " + message);

    /// <summary><c>featRegMore</c>: the pattern of the flag lists of the GFF readers, compiled by the first of <c>featFlagSet</c> / <c>featGff3FlagSet</c> that runs.</summary>
    public System.Text.RegularExpressions.Regex? RegMore;

    /// <summary>The count of <c>ajSeqSetNameMulti</c> (a static of <c>ajseq.c</c>): the sequences of the run that were given a name of the form EMBOSS_001.</summary>
    public int SeqNameCount;

    /// <summary><c>ajUtilGetProgram</c>: the source of a feature without one.</summary>
    public string Program = "featcopy";

    /// <summary><c>acdInTypeSeqName</c>: the ACD type (lower case) of the first sequence input of the program; the type of a sequence output that has none of its own (null: no input yet).</summary>
    public string? InputType;

    /// <summary>The library called <c>ajFatal</c>: the program ended at once with the message (no <c>Died:</c> line of the ACD follows).</summary>
    public bool Fatal;

    /// <summary>The input stores the text of the entries it reads (<c>Input-&gt;Text</c>: the ACD attribute <c>entry</c> of <c>feattext</c>).</summary>
    public bool StoreText;

    // featSaveGroup / featSaveGroupStr / featSaveExon of ajfeatread.c
    public uint SaveGroup;
    public string? SaveGroupStr;
    public uint SaveExon;

    // the "static AjPFeature gf" of featEmblFromLine / featSwissFromLine / featRefseqpFromLine / featPirFromLine
    public EmbossSeqFeature? EmblGf;
    public EmbossSeqFeature? SwissGf;
    public EmbossSeqFeature? RefseqpGf;
    public EmbossSeqFeature? PirGf;
}

/// <summary>A feature of a feature table as the EMBOSS library holds it (<c>AjSFeature</c>): the (internal) type, the source, the location, the strand and the tags in the order they were added; a join / order location is one feature with its exons as <see cref="Subfeatures"/>.</summary>
public sealed class EmbossSeqFeature
{
    public string Source = "";
    public string Type = "";
    public int Start, End, Start2, End2, Frame, Exon;
    public uint Group;
    public float Score;
    public bool Protein;
    public string? Remote, Label;
    public EmbossLocationBits Flags;
    public char Strand;
    public List<(string Name, string Value)> Tags = [];

    /// <summary>The GFF3 special tags (<c>ID</c>, <c>Name</c>, <c>Alias</c>, <c>Parent</c>, ...; <c>GffTags</c>).</summary>
    public List<(string Name, string Value)> GffTags = [];
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
public sealed partial class EmbossSeqFeatureTable
{
    /// <summary>The features in table order (the parents of the join / order locations).</summary>
    public List<EmbossSeqFeature> Features = [];

    /// <summary>The length of the table (<c>ajFeattableSetLength</c>: the length of the sequence once it is read; while reading, the largest position seen).</summary>
    public int Length;

    /// <summary>The table was reversed (<c>thys-&gt;Reversed</c>).</summary>
    public bool Reversed;

    /// <summary>The entry name (<c>thys-&gt;Seqid</c>).</summary>
    public string SeqId = "";

    /// <summary>The reader found at least one feature (<c>featReadEmbl</c> returned true): the library returns no table otherwise.</summary>
    public bool Found;

    /// <summary>The sequence type of the table (<c>thys-&gt;Type</c>: <c>N</c>, <c>P</c>, or empty while it is not known).</summary>
    public string Type = "";

    /// <summary>The first and last position used (<c>thys-&gt;Start</c>, <c>thys-&gt;End</c>) and the offset when trimmed.</summary>
    public uint Start, End, Offset;

    /// <summary>The largest group number (<c>thys-&gt;Groups</c>).</summary>
    public uint Groups;

    /// <summary>The name of the input format the table was read with (<c>thys-&gt;Formatstr</c>) and its index in the table of input formats (<c>thys-&gt;Format</c>).</summary>
    public string? FormatName;

    public int Format;

    /// <summary>The table is of a circular molecule (<c>thys-&gt;Circular</c>).</summary>
    public bool Circular;

    private EmbossFeatureReadContext _context = new();

    /// <summary>The state of the feature library the table is read with.</summary>
    public EmbossFeatureReadContext Context => _context;

    /// <summary>A new, empty table (<c>ajFeattableNew</c>).</summary>
    public EmbossSeqFeatureTable(EmbossFeatureReadContext? context = null)
    {
        _context = context ?? new EmbossFeatureReadContext();
    }

    // ---------------------------------------------------------------------------------------------------------------- reading

    /// <summary>
    /// <c>featReadEmbl</c> over the lines of ONE entry: the <c>FT</c> lines of an EMBL entry, the feature lines (five blanks) after <c>FEATURES</c> of a GenBank entry. Reading stops at <c>//</c>.
    /// </summary>
    /// <param name="entryLines">The lines of the entry.</param>
    /// <param name="context">The state of the feature library shared by the entries of one program run (<see cref="EmbossFeatureReadContext"/>); null = a fresh one.</param>
    public static EmbossSeqFeatureTable Read(IEnumerable<string> entryLines, EmbossFeatureReadContext? context = null)
    {
        ArgumentNullException.ThrowIfNull(entryLines);
        var table = new EmbossSeqFeatureTable(context);
        var lines = new EmbossFeatureLines([.. entryLines]);
        table.Found = table.ReadEmbl(lines);
        return table;
    }

    /// <summary>The text of the entry the reader has read (<c>ajFeattableGetEntry</c>: only the EMBL / GenBank and the GFF3 readers keep it, when the input has the ACD attribute <c>entry</c>); null = nothing stored.</summary>
    public string? TextPtr;

    // ajBuffreadLinePosStore: a line is appended to the text of the table when the input stores it
    private string? NextStored(EmbossFeatureLines file)
    {
        string? line = file.Next();
        if (line is not null && _context.StoreText)
            TextPtr = (TextPtr ?? "") + line + "\n";
        return line;
    }

    /// <summary><c>featReadEmbl</c>: false when the lines held no feature.</summary>
    internal bool ReadEmbl(EmbossFeatureLines file)
    {
        bool found = false, isGenbank = false;
        string saveFeature = "", saveLoc = "", saveLine = "";
        SetType("N");
        string? raw;
        while ((raw = NextStored(file)) is not null)
        {
            string line = raw;
            if (line.StartsWith("//", StringComparison.Ordinal))
                break;
            if (line.StartsWith("ID   ", StringComparison.Ordinal))
            {
                isGenbank = false;
                string[] tokens = line.Split([' ', ';', '\t', '\n', '\r'], StringSplitOptions.RemoveEmptyEntries);
                SeqId = tokens.Length > 1 ? tokens[1] : "";
            }
            else if (line.StartsWith("LOCUS ", StringComparison.Ordinal) || line.StartsWith("FEATURES  ", StringComparison.Ordinal))
            {
                isGenbank = true;
            }
            else if (line.StartsWith("FT   ", StringComparison.Ordinal) || (isGenbank && line.StartsWith("     ", StringComparison.Ordinal)))
            {
                if (FromLine(line.TrimEnd(' ', '\t', '\n', '\r'), ref saveFeature, ref saveLoc, ref saveLine, false))
                    found = true;
            }
            else if (isGenbank)
            {
                isGenbank = false;
            }
        }

        if (FromLine(null, ref saveFeature, ref saveLoc, ref saveLine, false))
            found = true;
        return found;
    }

    /// <summary><c>ajFeattableSetNuc</c> / <c>ajFeattableSetProt</c>: the type is set when it is still empty.</summary>
    public void SetType(string type)
    {
        if (Type.Length == 0)
        {
            Type = type;
            TypeNull = false;
        }
    }

    /// <summary>The type of the table is the NULL string of the library (a new table: only a reader that sets the type, or a failed attempt that did and was cleared, makes it the empty string) and prints as <c>&lt;null&gt;</c> where the writers print it with <c>%S</c>.</summary>
    public bool TypeNull = true;

    // featEmblFromLine / featRefseqpFromLine (the two differ only in the processing of a finished feature): the feature the call finished is returned (the "static" one of the library: the last one, even when this call finished none)
    internal bool FromLine(string? origLine, ref string saveFeature, ref string saveLoc, ref string saveLine, bool refseqp)
    {
        string line = "";
        bool newFeature = false, doFeature = false;
        if (origLine is not null)
        {
            line = origLine.Length >= 5 ? origLine[5..] : "";
            if (line.Length == 0 || line[0] != ' ')                 // ajStrGetCharFirst of an empty string is the end-of-string character: not a blank
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

        EmbossSeqFeature? gf = refseqp ? _context.RefseqpGf : _context.EmblGf;
        if (doFeature)
        {
            gf = refseqp ? ProcessProtein(saveFeature, ref saveLoc, ref saveLine) : Process(saveFeature, ref saveLoc, ref saveLine);
            if (refseqp)
                _context.RefseqpGf = gf;
            else
                _context.EmblGf = gf;
            saveLoc = "";
            saveLine = "";
        }

        if (origLine is null)
            return refseqp ? gf is not null : doFeature && gf is not null;
        line = RemoveWhiteExcess(line);
        if (newFeature)
        {
            int p = 0;
            saveFeature = NextToken(line, ref p, " ");
            string loc = NextToken(line, ref p, " /");
            if (loc.Length > 0 || p < line.Length)
            {
                saveLoc = loc;
                saveLine = p < line.Length ? line[p..] : "";
            }
            else
            {
                saveLoc = loc;
                saveLine = "";
            }

            return gf is not null;
        }

        if (saveLine.Length == 0 && (line.Length == 0 || line[0] != '/'))
        {
            int p = 0;
            string temp = NextToken(line, ref p, " ");
            if (temp.Length > 0)
                saveLoc += temp;
            string rest = p < line.Length ? line[p..] : "";
            if (rest.Length > 0)
                saveLine += rest;
            return gf is not null;
        }

        saveLine += " " + line;                                              // tag-values continued
        return gf is not null;
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

    private EmbossSeqFeature? Process(string feature, ref string loc, ref string tags) => ProcessCore(feature, ref loc, ref tags, false);

    private EmbossSeqFeature? ProcessProtein(string feature, ref string loc, ref string tags) => ProcessCore(feature, ref loc, ref tags, true);

    // featEmblProcess / featRefseqpProcess
    private EmbossSeqFeature? ProcessCore(string feature, ref string loc, ref string tags, bool protein)
    {
        EmbossSeqFeature? ret = null, gf = null;
        string source = protein ? "REFSEQP" : "EMBL";
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
                else
                    _context.Warn($"{SeqId}: unrecognised operator '{opnam}()' in '{opval}'");
                opval = inner;
            }
        }

        string previousOpval = "\u0001";
        while (opval.Length > 0)
        {
            if (opval == previousOpval)
                throw new InvalidOperationException($"The location '{loc}' is never consumed: the reference program loops forever on it (a location that begins with a comma).");
            previousOpval = opval;
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
            else if (LocRange(begstr, ref begNum, ref beg2))
            {
                begBound = false;
            }
            else
            {
                isLabel = true;
                simple = true;
                featLabel = begstr;
                _context.Warn($"{SeqId}: Simple feature location '{begstr}' in '{locstr}'");
            }

            if (LocNum(endstr, out endBound, out endNum))
            {
                end2 = 0;
            }
            else if (LocRange(endstr, ref end2, ref endNum))
            {
                endBound = false;
            }
            else
            {
                isLabel = true;
                simple = true;
                featLabel = endstr;
                _context.Error($"{SeqId}: Simple feature end '{begstr}' in '{locstr}'");
            }

            var flags = EmbossLocationBits.None;
            strand = locFwd ? '+' : '-';
            if (isLabel)
                _context.Warn($"{SeqId}: Feature location with label '{locstr}'");
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
            string type = protein ? EmbossFeatureVocab.TypeInternal(EmbossFeatureVocab.Get("refseqp").Types, feature) ?? "" : EmbossFeatureVocab.TypeInternal(EmbossFeatureVocab.Dna.Types, feature) ?? "";
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

                gf = protein ? NewProt(source, type, (int)tmpBeg, (int)tmpEnd, 0f, tmpFlags, null) : NewNuc(source, type, (int)tmpBeg, (int)tmpEnd, 0f, strand, frame, exon, (int)beg2, (int)end2, tmpFeatId, featLabel, tmpFlags, null);
            }

            if (!parent || opval.Length > 0)
            {
                exon++;
                if (protein)
                    NewProt(source, type, (int)begNum, (int)endNum, 0f, flags, gf);
                else
                    NewNuc(source, type, (int)begNum, (int)endNum, 0f, strand, frame, exon, (int)beg2, (int)end2, featId, featLabel, flags, gf);
                if (!remoteId)
                {
                    if (gf is null)
                        throw new InvalidOperationException("feature location without a feature");
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

        string? tagName = null, tagRest = null;
        while (tags.Length > 0)
        {
            if (TvTagVal(ref tags, out string tag, out string val))
            {
                val = QuoteStrip(val);
                if (!protein)
                {
                    if (tag == "codon_start")
                    {
                        frame = ToInt(val);
                        gf!.Frame = frame;
                    }

                    if (tag == "note" && val.Length > 0 && val[0] == '*')
                    {
                        val = val[1..];
                        if (ExtractFirst(val, out string word, out string restText))
                        {
                            tagName = word;
                            tagRest = restText;
                        }

                        if (tagName == "Type")
                        {
                            if (tagRest is not null && ExtractFirst(tagRest, out string typeWord, out string afterType))
                            {
                                gf!.Type = typeWord;
                                tagName = afterType;
                            }
                        }
                        else if (tagName is not null)
                        {
                            GfftagAdd(gf!, tagName, tagRest ?? "");
                        }

                        continue;
                    }

                    if (!TagAdd(ret!, tag, val))
                        _context.Warn($"{SeqId}: Bad value '{val}' for tag '/{tag}'");
                }
                else
                {
                    if (!TagAdd(ret!, tag, val))
                        _context.Warn($"{SeqId}: Bad value '{val}' for tag '/{tag}'");
                    if (tag == "codon_start")
                    {
                        frame = ToInt(val);
                        gf!.Frame = frame;
                    }
                }
            }
            else if (TvRest(ref tags, out string skipped))
            {
                // anything non-white up to '/' is bad: skipped
                _context.Warn($"Bad feature syntax {SeqId}: skipping '{skipped}'");
            }
            else
            {
                _context.Warn($"Bad feature syntax {SeqId}: giving up at '{tags}'");
                tags = "";
            }
        }

        return ret;
    }

    // ajStrExtractFirst: the first white-delimited word and the rest after the blanks that follow it; false (and nothing set) for an empty string or one that starts with white
    internal static bool ExtractFirst(string s, out string first, out string rest)
    {
        first = rest = "";
        if (s.Length == 0 || IsCSpace(s[0]))
            return false;
        int j = 1;
        while (j < s.Length && !IsCSpace(s[j]))
            j++;
        first = s[..j];
        while (j < s.Length && IsCSpace(s[j]))
            j++;
        rest = s[j..];
        return true;
    }

    internal static bool IsCSpace(char c) => c is ' ' or '\t' or '\n' or '\r' or '\f' or '\v';

    /// <summary><c>ajFeatTagAddSS</c>: the tag is stored under its canonical name of the vocabulary of the feature (protein or DNA); a value that fails its check, or a tag the vocabulary does not know, becomes a <c>note</c> <c>*tag: value</c>.</summary>
    /// <returns>False when the value was bad.</returns>
    public bool TagAdd(EmbossSeqFeature f, string tag, string value)
    {
        var tags = f.Protein ? EmbossFeatureVocab.Protein.Tags : EmbossFeatureVocab.Dna.Tags;
        string tmptag = EmbossFeatureVocab.TagNameOf(tags, tag, out bool known);
        string fmt = EmbossFeatureVocab.TagFormat(tags, tmptag);
        string valTmp = value, tagTmp = tmptag;
        bool ret = true;
        string kind = fmt.Length >= 2 ? fmt[..2] : fmt;
        switch (kind)
        {
            case "LI" or "QL" or "TE" or "VO" or "QW":
                break;
            case "QS" or "SP":
                if (!SpecialValue(tmptag, ref valTmp, _context))
                {
                    ret = false;
                    tagTmp = EmbossFeatureVocab.Dna.Tags.Fetch("") ?? "";
                    valTmp = "*" + tmptag + ": " + value;
                }

                break;
            case "QT":
                if (!known)
                    valTmp = "*" + tag + ": " + value;
                break;
        }

        string finalTag = EmbossFeatureVocab.TagNameOf(f.Protein ? EmbossFeatureVocab.Protein.Tags : EmbossFeatureVocab.Dna.Tags, tagTmp, out _);
        f.Tags.Add((finalTag, valTmp));
        return ret;
    }

    private static readonly HashSet<string> Gff3PredefinedTags = ["ID", "Name", "Alias", "Parent", "Target", "Gap", "Derives_from", "Note", "Dbxref", "Ontology_term", "Is_circular"];

    /// <summary><c>ajFeatGfftagAddSS</c>: the value is split at commas, one tag-value each; a predefined GFF3 tag goes to <see cref="EmbossSeqFeature.GffTags"/>, any other to the ordinary tags (as it is spelled).</summary>
    internal static int GfftagAdd(EmbossSeqFeature f, string tag, string value)
    {
        int n = 0;
        bool predefined = Gff3PredefinedTags.Contains(tag);
        foreach (string piece in SplitTokens(value, ','))
        {
            n++;
            (predefined ? f.GffTags : f.Tags).Add((tag, piece));
        }

        return n;
    }

    // ajStrTokenNewC(value, ","): ajStrTokenNextParse steps over the delimiters (an empty piece between two of them is skipped) and stops at the end
    private static IEnumerable<string> SplitTokens(string value, char delimiter)
    {
        int p = 0;
        while (true)
        {
            while (p < value.Length && value[p] == delimiter)
                p++;
            if (p >= value.Length)
                yield break;
            int q = p;
            while (p < value.Length && value[p] != delimiter)
                p++;
            yield return value[q..p];
        }
    }


    // ajStrToInt as the library uses it for a codon_start: a string that is not entirely an integer gives 0
    private static int ToInt(string s) => CStr.ToInt(s, out int v) ? v : 0;

    /// <summary><c>ajFeatNewNucFlags</c> (<paramref name="parent"/> null) or <c>ajFeatNewNucFlagsSub</c>: a nucleotide feature, appended to the table or to the sub-features of the parent.</summary>
    public EmbossSeqFeature? NewNuc(string source, string type, int start, int end, float score, char strand, int frame, int exon, int start2, int end2, string? entryId, string? label, EmbossLocationBits flags, EmbossSeqFeature? parent)
    {
        if (type.Length == 0)
            return null;
        Groups++;
        var f = new EmbossSeqFeature
        {
            Group = Groups,
            Exon = exon,
            Source = source.Length > 0 ? source : _context.Program,
            Type = Type == "P" ? ProteinInternalType(type) : DnaInternalType(type),
            Score = score,
            Flags = flags,
            Strand = strand is '+' or '-' ? strand : '\0',
            Frame = frame,
            Start = start,
            End = end,
            Start2 = start2,
            End2 = end2,
        };
        if (!string.IsNullOrEmpty(entryId))
        {
            f.Remote = entryId;
        }
        else if (!flags.HasFlag(EmbossLocationBits.RemoteId) && !flags.HasFlag(EmbossLocationBits.Label))
        {
            Length = Math.Max(Length, f.Start);
            Length = Math.Max(Length, f.End);
        }

        if (!string.IsNullOrEmpty(label))
        {
            f.Label = label;
            _context.Warn($"{SeqId}: Feature label '{label}' used");
        }
        Attach(f, parent);
        return f;
    }

    /// <summary><c>ajFeatNewProtFlags</c> / <c>ajFeatNewProtFlagsSub</c>: a protein feature.</summary>
    public EmbossSeqFeature? NewProt(string source, string type, int start, int end, float score, EmbossLocationBits flags, EmbossSeqFeature? parent)
    {
        if (type.Length == 0)
            return null;
        Groups++;
        var f = new EmbossSeqFeature
        {
            Group = Groups,
            Source = source.Length > 0 ? source : _context.Program,
            Type = ProteinInternalType(type),
            Score = score,
            Flags = flags,
            Strand = '\0',
            Start = start,
            End = end,
            Protein = true,
        };
        if (!flags.HasFlag(EmbossLocationBits.RemoteId) && !flags.HasFlag(EmbossLocationBits.Label))
        {
            Length = Math.Max(Length, f.Start);
            Length = Math.Max(Length, f.End);
        }

        Attach(f, parent);
        return f;
    }

    private void Attach(EmbossSeqFeature f, EmbossSeqFeature? parent)
    {
        if (parent is null)
        {
            Features.Add(f);                                             // ajFeattableAdd
            return;
        }

        parent.Subfeatures ??= [];
        parent.Subfeatures.Add(f);
    }

    // featTypeDna / featTypeProt: the internal type for a type name
    private static string DnaInternalType(string type) => EmbossFeatureVocab.TypeInternal(EmbossFeatureVocab.Dna.Types, type) ?? "";

    private static string ProteinInternalType(string type) => EmbossFeatureVocab.TypeInternal(EmbossFeatureVocab.Protein.Types, type) ?? "";

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

        // ajStrAssignSubS: a negative end counts from the end of the string (ipos = 0 gives the whole string)
        if (idpos > 0)
        {
            entryId = Sub(loc, 0, idpos - 2);
            value = Sub(loc, idpos, ipos - 1);
        }
        else
        {
            entryId = "";
            value = Sub(loc, 0, ipos - 1);
        }

        rest = ipos < loc.Length ? Sub(loc, ipos + 1, -1) : "";
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
    // ajStrAssignSubS with ajCvtSposToPos / ajCvtSposToPosStart: the positions are clamped to the string (a negative one counts from the end) and the result is NEVER shorter than one character
    internal static string Sub(string s, int first, int last)
    {
        int len = s.Length;
        if (len == 0)
            return "";
        int ibegin = first < 0 ? Math.Max(0, len + first) : first;
        if (ibegin >= len)
            ibegin = len - 1;
        int iend = last < 0 ? Math.Max(0, len + last) : last;
        if (iend >= len)
            iend = len - 1;
        if (iend < ibegin)
            iend = ibegin;
        return s.Substring(ibegin, iend - ibegin + 1);
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
    private static bool LocRange(string loc, ref uint num1, ref uint num2)
    {
        if (loc.Length == 0 || loc[0] != '(' || loc[^1] != ')')
            return false;
        num1 = num2 = 0;
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

    private const System.Text.RegularExpressions.RegexOptions Rx = System.Text.RegularExpressions.RegexOptions.CultureInvariant;

    private static bool Matches(string pattern, string value) => System.Text.RegularExpressions.Regex.IsMatch(value, pattern, Rx);

    // ajFeattagSpecial with its messages (ajFeatWarn of the validators)
    internal static bool SpecialValue(string tag, ref string value, EmbossFeatureReadContext context, string dialect = "EMBL")
    {
        string before = value;
        if (dialect == "GFF2" && tag == "ncRNA_class")
        {
            context.Warn($"Unrecognised special GFF feature tag '{tag}'");                       // ajFeattagSpecialGff2 knows no ncRNA_class
            return false;
        }

        bool ok = SpecialCheck(tag, ref value, context);
        switch (tag)
        {
            case "anticodon" or "citation" or "collection_date" or "compare" or "db_xref" or "lat_lon" or "mobile_element" or "ncRNA_class" or "PCR_primers" or "protein_id" or "rpt_unit" or "transl_except":
                if (!ok)
                    context.Warn($"bad /{tag} value '{before}'");
                break;
            case "bio_material" or "culture_collection" or "specimen_voucher":
                if (!ok)
                    context.Warn($"bad /bio_material value '{before}'");
                break;
            case "rpt_unit_range" or "tag_peptide":
                if (!ok)
                    context.Warn($"bad /rpt_unit_tange value '{before}'");
                break;
            case "codon":
                if (ok && value != before)
                    context.Warn($"unquoted /codon value corrected to '{value}'");
                if (!ok)
                    context.Warn($"bad /codon value '{before}'");
                break;
            case "cons_splice":
                if (!ok)
                    context.Warn($"bad /cons_splice value '{before}'");
                if (ConsSpliceIsLow(before, out string beg, out string end))
                {
                    value = "(5'site:" + beg + ",3'site:" + end + ")";
                    context.Warn($"bad /cons_splice value corrected to '{value}'");
                    if (!ok)
                        value = before;
                }

                break;
            case "estimated_length":
                if (!ok)
                    context.Warn($"bad /estimated_length value '{before}' corrected to 'unknown'");
                break;
            case "inference":
                if (!ok)
                {
                    var m = System.Text.RegularExpressions.Regex.Match(before, "^([^:]+)(:([^:]+)?)*$", Rx);
                    context.Warn($"bad /inference value '{before}' type: '{(m.Success ? m.Groups[1].Value : "<null>")}'");
                }

                break;
            case "replace":
                if (!ok)
                    context.Warn($"bad /replace value '{value}'");
                break;
            case "rpt_unit_seq":
                if (!ok)
                    context.Warn($"bad /rpt_unit_seq value '{value}'");
                break;
            case "translation":
                if (!ok)
                    context.Warn($"bad /translation value '{RemoveWhite(before)}'");
                break;
            default:
                context.Warn($"Unrecognised special {(dialect == "EMBL" ? "EMBL" : "GFF")} feature tag '{tag}'");
                break;
        }

        return ok;
    }

    // the lower-case test of featTagSpecialAllConssplice: true (with the two words) when a letter of either is lower case
    private static bool ConsSpliceIsLow(string val, out string beg, out string end)
    {
        beg = "";
        end = "<null>";                                                // the word that was never read prints as <null>
        if (!val.StartsWith("(5'site:", StringComparison.Ordinal) || val[^1] != ')')
            return false;
        int cp = 8;
        bool low = false;
        while (cp < val.Length && val[cp] != ',')
        {
            if (char.IsAsciiLetterLower(val[cp]))
                low = true;
            cp++;
        }

        beg = val[8..cp];
        if (string.CompareOrdinal(val, cp, ",3'site:", 0, 8) == 0)
        {
            cp += 8;
            int cq = cp;
            while (cp < val.Length && val[cp] != ')')
            {
                if (char.IsAsciiLetterLower(val[cp]))
                    low = true;
                cp++;
            }

            end = val[cq..cp];
        }

        return low;
    }

    // ajFeattagSpecial: the checks (and the corrections) of the special tags; false = the value is bad and the tag becomes a note
    internal static bool SpecialCheck(string tag, ref string value, EmbossFeatureReadContext context)
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
        int begin = PositionII(Length, 1, unchecked((int)Start));
        int end = End != 0 ? PositionII(Length, unchecked((int)End), seqLength) : Length;
        if (U(end) > U(seqLength))
            end = seqLength;
        if (U(begin) <= U(offset))
            begin = PositionII(Length, offset + 1, end);
        bool doBegin = begin > 1, doEnd = U(end) < U(Length);
        Features.RemoveAll(f => !TrimOffRange(f, offset, begin, end, doBegin, doEnd));
        Offset = unchecked((uint)offset);
        if (Rev)
            ReverseTable();
        Trimmed = true;
    }

    // ajFeattablePosII (unsigned: a negative result is a huge number, which is then cut to the length)
    private static int PositionII(int ilen, int imin, int ipos)
    {
        uint jpos;
        if (ipos < 0)
            jpos = unchecked((uint)(ilen + ipos + 1));
        else
            jpos = ipos != 0 ? (uint)ipos : 1u;
        if (jpos > unchecked((uint)ilen))
            jpos = unchecked((uint)ilen);
        if (jpos < unchecked((uint)imin))
            jpos = unchecked((uint)imin);
        return unchecked((int)jpos);
    }

    // ajFeatTrimOffRange: false deletes the feature
    // the positions of the library are unsigned: a position before the start of the sequence is a huge number
    private static uint U(int v) => unchecked((uint)v);

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
            if (U(ft.Start) > U(end))
                return false;
            if (U(ft.Start2) > U(end))
                ft.Start2 = end;
            if (U(ft.End) > U(end))
            {
                ft.End = end;
                ft.Flags |= EmbossLocationBits.EndAfterSeq;
            }

            if (U(ft.End2) > U(end))
                ft.End2 = end;
        }

        if (doBegin)
        {
            if (U(ft.End) < U(begin))
                return false;
            if (ft.End2 != 0 && U(begin) > U(offset + 1) && U(ft.End2) < U(begin))
                ft.End2 = begin;
            if (ft.Start != 0 && U(ft.Start) < U(begin))
            {
                ft.Start = begin;
                ft.Flags |= EmbossLocationBits.StartBeforeSeq;
            }

            if (ft.Start2 != 0 && U(ft.Start2) < U(begin))
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
