#pragma warning disable S1121, S3358, S2325 // the readers keep the control flow of ajfeatread.c (assignments in conditions, the static state of the library)
using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Infrastructure;

// The feature readers of ajfeatread.c other than the GenBank / EMBL one (featReadRefseq, featReadGenpept, featReadRefseqp, featReadPir, featReadSwiss, featReadGff2, featReadGff3, featReadGff3old, featReadDasgff, featReadChado)
// with their line parsers and the group / flag setters (EMBOSS 6.6.0, kimrutherford/EMBOSS at raw.githubusercontent.com). A5-8c3e, B09 F241.
internal sealed partial class EmbossSeqFeatureTable
{
    private const RegexOptions Ci = RegexOptions.CultureInvariant;

    private static readonly Regex SwNew = new("^FT   (([^ ]+) +([?<]?[0-9]+|[?]) +([?>]?[0-9]+|[?]) *)(.*)$", Ci);
    private static readonly Regex SwNext = new("^FT    +(.*)$", Ci);
    private static readonly Regex SwComment = new("^(.*)[(]([^)]+)[)]$", Ci);
    private static readonly Regex SwFtid = new("^(.*)/FTId=([^ .]+)$", Ci);
    private static readonly Regex PirAll = new("^F;([^/]+)/([^:]+):([^#]*)", Ci);
    private static readonly Regex PirCom = new("^#([^#]*)", Ci);
    private static readonly Regex PirLoc = new("^([^,]+),?", Ci);
    private static readonly Regex PirPos = new("^([^-]+)-?", Ci);
    private static readonly Regex GffBlank = new("^[ ]*$", Ci);
    private static readonly Regex GffVersion = new("^##gff-version[ ]+([0-9]+)", Ci);
    private static readonly Regex GffRegion = new("^##sequence-region[ ]+([0-9a-zA-Z]+)[ ]+([\\+-]?[0-9]+)[ ]+([\\+-]?[0-9]+)", Ci);
    private static readonly Regex GffComment = new("^#[ ]*(.*)", Ci);
    private static readonly Regex GffType = new("^##[Tt]ype +(\\S+)( +(\\S+))?", Ci);
    private static readonly Regex GffTvTagval = new(" *([^ =]+)[ =](\"[^\"]*\"|[^;]+)(;|$)", Ci);
    private static readonly Regex Gff3Directive = new("^##(.*)", Ci);
    private static readonly Regex Gff3Comment = new("^#(.*)", Ci);
    private static readonly Regex Gff3TvTagval = new(" *([^ =]+)[ =]([^;]+)*(;|$)", Ci);
    private static readonly Regex FlagRegex = new("[ \"]*(0x[0-9a-f]+)", Ci);
    private static readonly Regex More2 = new("[ \"]*([^:]+):([^: \"]+)", Ci);
    private static readonly Regex More3 = new("[,]*([^:]+):([^,]+)", Ci);
    private static readonly Regex Group2 = new("^\"(([^.]*)[.])?([0-9]+)", Ci);
    private static readonly Regex Group3 = new("^(([^.]*)[.])?([0-9]+)", Ci);

    // ajRegPost: the text after the match
    private static string Post(Match m, string s) => s[(m.Index + m.Length)..];

    private static string Trim4(string s) => s.Trim('\t', ' ', '\n', '\r');

    // ---------------------------------------------------------------------------------------------------------------- refseq / genpept / refseqp

    /// <summary><c>featReadRefseq</c> (nucleotide), <c>featReadGenpept</c> (protein) and <c>featReadRefseqp</c>: the feature lines (five blanks) of the file, to the end of it.</summary>
    internal bool ReadFiveBlanks(EmbossFeatureLines file, bool protein, bool refseqp)
    {
        bool found = false;
        string saveFeature = "", saveLoc = "", saveLine = "";
        SetType(protein ? "P" : "N");
        string? raw;
        while ((raw = file.Next()) is not null)
        {
            if (raw.StartsWith("     ", StringComparison.Ordinal) && FromLine(raw.TrimEnd(' ', '\t', '\n', '\r'), ref saveFeature, ref saveLoc, ref saveLine, refseqp))
                found = true;
        }

        if (FromLine(null, ref saveFeature, ref saveLoc, ref saveLine, refseqp))
            found = true;
        return found;
    }

    // ajFeattableSetDefname
    private void SetDefname(string name)
    {
        if (SeqId.Length > 0)
            return;
        if (name.Length > 0)
            SeqId = name;
    }

    // ---------------------------------------------------------------------------------------------------------------- swiss

    /// <summary><c>featReadSwiss</c>.</summary>
    internal bool ReadSwiss(EmbossFeatureLines file)
    {
        bool found = false, setname = false;
        string saveFeature = "", saveFrom = "", saveTo = "", saveLine = "";
        string? raw;
        while ((raw = file.Next()) is not null)
        {
            string line = Trim4(raw);
            if (line.StartsWith("//", StringComparison.Ordinal))
                break;
            if (line.StartsWith("ID   ", StringComparison.Ordinal))
            {
                var tok = new StrTok(line, " \n\r");
                tok.Next(out _);
                tok.Next(out string defname);
                if (!setname)
                {
                    SetDefname(defname);
                    setname = true;
                }
            }
            else if (line.StartsWith("AC   ", StringComparison.Ordinal))
            {
                var tok = new StrTok(line, " \n\r,;");
                tok.Next(out _);
                tok.Next(out string defname);
                if (!setname)
                {
                    SetDefname(defname);
                    setname = true;
                }
            }
            else if (line.StartsWith("FT   ", StringComparison.Ordinal) && SwissFromLine(line, ref saveFeature, ref saveFrom, ref saveTo, ref saveLine))
            {
                found = true;
            }
        }

        if (SwissFromLine(null, ref saveFeature, ref saveFrom, ref saveTo, ref saveLine))
            found = true;
        if (found)
            SetType("P");
        return found;
    }

    // featSwissFromLine
    private bool SwissFromLine(string? origline, ref string saveFeature, ref string saveFrom, ref string saveTo, ref string saveLine)
    {
        Match? mNew = origline is not null ? SwNew.Match(origline) : null;
        bool newft = mNew is { Success: true };
        EmbossSeqFeature? gf = _context.SwissGf;
        if (newft || origline is null)
        {
            gf = saveFrom.Length > 0 ? SwissProcess(saveFeature, saveFrom, saveTo, "SWISSPROT", saveLine) : null;
            _context.SwissGf = gf;
            saveFeature = saveFrom = saveTo = saveLine = "";
        }

        if (origline is null)
            return gf is not null;
        if (newft)
        {
            saveFeature = mNew!.Groups[2].Value;
            saveFrom = mNew.Groups[3].Value;
            saveTo = mNew.Groups[4].Value;
            saveLine = mNew.Groups[5].Value;
            saveFeature = Trim4(saveFeature);
            return gf is not null;
        }

        var mNext = SwNext.Match(origline);
        if (mNext.Success)
            saveLine = saveLine + " " + mNext.Groups[1].Value;
        else
            _context.Warn($"{SeqId}: Bad SwissProt feature line:\n{origline}");
        return gf is not null;
    }

    // featSwissProcess
    private EmbossSeqFeature? SwissProcess(string feature, string fromstr, string tostr, string source, string tags)
    {
        EmbossLocationBits flags = 0;
        int start, end;
        switch (fromstr.Length > 0 ? fromstr[0] : '\0')
        {
            case '?':
                flags |= EmbossLocationBits.StartUnsure;
                start = CStr.ToInt(fromstr[1..], out int a) ? a : 0;
                break;
            case '<' or '>':
                flags |= EmbossLocationBits.StartBeforeSeq;
                start = CStr.ToInt(fromstr[1..], out int b) ? b : 0;
                break;
            default:
                start = CStr.ToInt(fromstr, out int c) ? c : 0;
                break;
        }

        switch (tostr.Length > 0 ? tostr[0] : '\0')
        {
            case '?':
                flags |= EmbossLocationBits.EndUnsure;
                end = CStr.ToInt(tostr[1..], out int a) ? a : 0;
                break;
            case '<' or '>':
                flags |= EmbossLocationBits.EndAfterSeq;
                end = CStr.ToInt(tostr[1..], out int b) ? b : 0;
                break;
            default:
                end = CStr.ToInt(tostr, out int c) ? c : 0;
                break;
        }

        var ret = NewProt(source, feature, start, end, 0f, flags, null);
        if (ret is null)
            return null;
        string tagstr = TrimChars(tags, " .");
        string ftid = "";
        var mFtid = SwFtid.Match(tagstr);
        if (mFtid.Success)
        {
            ftid = mFtid.Groups[2].Value;
            tagstr = TrimChars(mFtid.Groups[1].Value, " .");
        }

        var mCom = SwComment.Match(tagstr);
        if (mCom.Success)
        {
            string note = TrimChars(mCom.Groups[1].Value, " .");
            string comment = mCom.Groups[2].Value;
            if (note.Length > 0)
                TagAdd(ret, "note", note);
            if (comment.Length > 0)
                TagAdd(ret, "comment", comment);
        }
        else if (tagstr.Length > 0)
        {
            TagAdd(ret, "note", tagstr);
        }

        if (ftid.Length > 0)
            TagAdd(ret, "ftid", ftid);
        return ret;
    }

    // ajStrTrimC: the characters of the set are removed from both ends
    private static string TrimChars(string s, string set) => s.Trim(set.ToCharArray());

    // ---------------------------------------------------------------------------------------------------------------- pir

    /// <summary><c>featReadPir</c>.</summary>
    internal bool ReadPir(EmbossFeatureLines file)
    {
        bool found = false;
        string? raw;
        while ((raw = file.Next()) is not null)
        {
            string line = Trim4(raw);
            if (line.StartsWith("F;", StringComparison.Ordinal) && PirFromLine(line))
                found = true;
        }

        return found;
    }

    // featPirFromLine
    private bool PirFromLine(string origline)
    {
        var m = PirAll.Match(origline);
        if (!m.Success)
            return false;
        string locstr = m.Groups[1].Value, typstr = m.Groups[2].Value, notestr = m.Groups[3].Value;
        typstr = typstr.Replace(' ', '_');
        string pirtype = EmbossFeatureVocab.TypeInternal(EmbossFeatureVocab.Get("pir").Types, typstr) ?? "";
        notestr = EmbossFeatureVocab.RemoveWhiteExcess(notestr);
        bool mother = true;
        EmbossSeqFeature? gf = _context.PirGf;
        string? posstr = null;
        uint start = 0, end = 0;
        Match ml;
        while ((ml = PirLoc.Match(locstr)).Success)
        {
            string exonstr = ml.Groups[1].Value;
            locstr = Post(ml, locstr);
            int i = 0;
            Match mp;
            while ((mp = PirPos.Match(exonstr)).Success)
            {
                posstr = mp.Groups[1].Value;
                if (i++ == 0 && !CStr.ToUint(posstr, out start))
                    start = 1;
                exonstr = Post(mp, exonstr);
            }

            if (!CStr.ToUint(posstr, out end))
                end = 1;
            if (mother)
                gf = NewProt("PIR", pirtype, unchecked((int)start), unchecked((int)end), 0f, 0, null);
            if (!mother || locstr.Length > 0)
            {
                NewProt("PIR", pirtype, unchecked((int)start), unchecked((int)end), 0f, 0, gf);
                if (gf!.Start == 0 || unchecked((uint)gf.Start) > start)
                    gf.Start = unchecked((int)start);
                if (gf.End == 0 || unchecked((uint)gf.End) < end)
                    gf.End = unchecked((int)end);
            }

            if (mother)
            {
                if (notestr.Length > 0)
                    TagAdd(gf!, "note", notestr);
                string temp = Post(m, origline);
                Match mc;
                while ((mc = PirCom.Match(temp)).Success)
                {
                    string comstr = EmbossFeatureVocab.RemoveWhiteExcess(mc.Groups[1].Value);
                    TagAdd(gf!, "comment", comstr);
                    temp = Post(mc, temp);
                }
            }

            mother = false;
        }

        if (mother)
            _context.Warn($"featPirFromLine: Did not understand location '{locstr}'");
        _context.PirGf = gf;
        return gf is not null;
    }

    // ---------------------------------------------------------------------------------------------------------------- gff

    // featGroupSet (gff2) / featGff3GroupSet
    private void GroupSet(EmbossSeqFeature gf, string? grouptag, bool gff3)
    {
        int len = grouptag?.Length ?? 0;
        if (len > 0 && _context.SaveGroupStr is not null && string.Equals(grouptag, _context.SaveGroupStr, StringComparison.OrdinalIgnoreCase))
        {
            gf.Group = _context.SaveGroup;
            gf.Exon = (int)++_context.SaveExon;
            return;
        }

        Match m = len > 0 ? (gff3 ? Group3 : Group2).Match(grouptag!) : Match.Empty;
        if (len > 0 && m.Success)
        {
            _context.SaveGroupStr = grouptag;
            string grpstr = m.Groups[3].Value;
            string namstr = m.Groups[2].Value;
            if (!gff3 && namstr.Length > 0 && !string.Equals(namstr, SeqId, StringComparison.OrdinalIgnoreCase))
                _context.Warn($"GFF group field '{grouptag}' for table '{SeqId}'");
            if (CStr.ToInt(grpstr, out int grpnum))
            {
                gf.Group = unchecked((uint)grpnum);
                _context.SaveGroup = unchecked((uint)grpnum);
            }
            else
            {
                gf.Group = ++Groups;
            }
        }
        else
        {
            gf.Group = ++Groups;
            _context.SaveGroup = gf.Group;
            gf.Exon = 0;
            _context.SaveExon = 0;
        }
    }

    // featFlagSet (gff2) / featGff3FlagSet
    private void FlagSet(EmbossSeqFeature gf, string flagstr, bool gff3)
    {
        // featRegMore is ONE static of the library: the first of the two functions to run compiles its pattern for both
        _context.RegMore ??= gff3 ? More3 : More2;
        string savstr = flagstr;
        var mf = FlagRegex.Match(savstr);
        if (mf.Success)
        {
            if (CStr.ToHex(mf.Groups[1].Value, out int flags))
            {
                var fl = (EmbossLocationBits)unchecked((uint)flags);
                if (fl.HasFlag(EmbossLocationBits.Group))
                    fl |= EmbossLocationBits.Order;
                if (fl.HasFlag(EmbossLocationBits.OneOf))
                    fl |= EmbossLocationBits.Order;
                if (gff3)
                    gf.Flags |= fl;
                else
                    gf.Flags = fl;
            }

            savstr = Post(mf, savstr);
        }

        Match mm;
        while ((mm = _context.RegMore.Match(savstr)).Success)
        {
            string typstr = mm.Groups[1].Value, valstr = mm.Groups[2].Value;
            bool Is(string name) => string.Equals(typstr, name, StringComparison.OrdinalIgnoreCase);
            if (gff3)
            {
                if (Is("start_before"))
                {
                    if (valstr == "true")
                        gf.Flags |= EmbossLocationBits.StartBeforeSeq;
                }
                else if (Is("end_after"))
                {
                    if (valstr == "true")
                        gf.Flags |= EmbossLocationBits.EndAfterSeq;
                }
                else if (Is("between"))
                {
                    if (valstr == "true")
                        gf.Flags |= EmbossLocationBits.BetweenSeq;
                }
                else if (Is("start2"))
                {
                    if (CStr.ToInt(valstr, out int n))
                    {
                        gf.Start2 = n;
                        gf.Flags |= EmbossLocationBits.StartTwo;
                    }
                }
                else if (Is("end2"))
                {
                    if (CStr.ToInt(valstr, out int n))
                    {
                        gf.End2 = n;
                        gf.Flags |= EmbossLocationBits.EndTwo;
                    }
                }
                else if (Is("label"))
                {
                    _context.Warn($"GFF3 label '{valstr}' used");
                    gf.Label = valstr;
                }
                else if (Is("start_unsure"))
                {
                    if (valstr == "true")
                        gf.Flags |= EmbossLocationBits.StartUnsure;
                }
                else if (Is("end_unsure"))
                {
                    if (valstr == "true")
                        gf.Flags |= EmbossLocationBits.EndUnsure;
                }
                else if (Is("type"))
                {
                    gf.Type = EmbossFeatureVocab.TypeInternal(EmbossFeatureVocab.Dna.Types, valstr) ?? "";
                }
                else
                {
                    _context.Warn($"Unknown GFF3 featflags type '{typstr}:{valstr}'");
                }
            }
            else if (Is("start2"))
            {
                if (CStr.ToInt(valstr, out int n))
                    gf.Start2 = n;
            }
            else if (Is("end2"))
            {
                if (CStr.ToInt(valstr, out int n))
                    gf.End2 = n;
            }
            else if (Is("remoteid"))
            {
                gf.Remote = valstr;
            }
            else if (Is("label"))
            {
                _context.Warn($"GFF label '{valstr}' used");
                gf.Label = valstr;
            }
            else
            {
                _context.Warn($"Unknown GFF FeatFlags type '{typstr}:{valstr}'");
            }

            savstr = Post(mm, savstr);
        }
    }

    // featGff2ProcessTagval
    private void Gff2ProcessTagval(EmbossSeqFeature gf, string groupfield, float version)
    {
        if (groupfield.Length == 0)
            return;
        if (Math.Abs(version - 1.0f) < 1e-5f)
        {
            GroupSet(gf, groupfield, false);
            return;
        }

        bool grpset = false;
        string tv = groupfield;
        if (tv.Contains("\\\"", StringComparison.Ordinal))
            tv = tv.Replace("\\\"", "\u0001", StringComparison.Ordinal);
        while (tv.Length > 0)
        {
            var m = GffTvTagval.Match(tv);
            if (!m.Success)
                throw new InvalidOperationException($"GFF tag parsing failed on '{tv}' in '{groupfield}': the reference program loops forever on such a group field.");
            string tmptag = m.Groups[1].Value;
            string val = Trim4(m.Groups[2].Value);
            tv = Post(m, tv);
            if (tmptag is "Sequence" or "ID")
            {
                GroupSet(gf, val, false);
                grpset = true;
            }
            else if (tmptag is "FeatFlags" or "featflags")
            {
                FlagSet(gf, val, false);
            }
            else
            {
                TagAdd(gf, tmptag, QuoteStripValue(val));
            }
        }

        if (!grpset)
            GroupSet(gf, null, false);
    }

    private static string QuoteStripValue(string v)
    {
        if (v.Length > 0 && v[^1] == '"')
            v = v[..^1];
        if (v.Length > 0 && v[0] == '"')
            v = v[1..];
        v = v.Replace("\"\"", "\"", StringComparison.Ordinal);
        return v.Replace("\\\"", "\"", StringComparison.Ordinal);
    }

    // ajStrFmtPercentDecode
    private static string PercentDecode(string s)
    {
        var sb = new StringBuilder();
        int p = 0;
        while (p < s.Length)
        {
            if (s[p] != '%')
            {
                sb.Append(s[p++]);
                continue;
            }

            int ia = 0;
            p++;
            if (p >= s.Length)
                break;
            ia += 16 * HexValue(s[p]);
            p++;
            if (p >= s.Length)
                break;
            ia += HexValue(s[p]);
            sb.Append((char)ia);
            p++;
        }

        return sb.ToString();
    }

    private static int HexValue(char c) => c switch { >= '0' and <= '9' => c - '0', >= 'A' and <= 'F' => 10 + c - 'A', >= 'a' and <= 'f' => 10 + c - 'a', _ => 0 };

    // featGff3ProcessTagval
    private string? Gff3ProcessTagval(EmbossSeqFeature gf, string groupfield, out bool parent)
    {
        parent = false;
        string? idstr = null;
        if (groupfield.Length == 0)
            return idstr;
        bool grpset = false;
        string tv = groupfield;
        while (tv.Length > 0)
        {
            var m = Gff3TvTagval.Match(tv);
            if (!m.Success)
            {
                _context.Warn($"GFF3 tag parsing failed on '{tv}' in\n'{groupfield}'");
                break;
            }

            string tmptag = m.Groups[1].Value;
            string val = Trim4(PercentDecode(m.Groups[2].Success ? m.Groups[2].Value : ""));
            tv = Post(m, tv);
            char cp = tmptag[0];
            if (cp is >= 'A' and <= 'Z')
            {
                if (tmptag == "ID")
                {
                    idstr = val;
                }
                else if (tmptag == "Parent")
                {
                    parent = true;
                    GroupSet(gf, val, true);
                    grpset = true;
                }
                else if (tmptag == "Is_circular")
                {
                    Circular = true;
                }

                GfftagAdd(gf, tmptag, val);
            }
            else if (tmptag == "featflags")
            {
                FlagSet(gf, val, true);
            }
            else
            {
                TagAdd(gf, tmptag, val);
            }
        }

        if (!grpset)
            GroupSet(gf, null, false);
        return idstr;
    }

    // featGff3oldProcessTagval
    private void Gff3oldProcessTagval(EmbossSeqFeature gf, string groupfield)
    {
        if (groupfield.Length == 0)
            return;
        bool grpset = false;
        string tv = groupfield;
        bool escapedquote = tv.Contains("\\\"", StringComparison.Ordinal);
        if (escapedquote)
            tv = tv.Replace("\\\"", "\u0001", StringComparison.Ordinal);
        while (tv.Length > 0)
        {
            var m = GffTvTagval.Match(tv);
            if (!m.Success)
                throw new InvalidOperationException($"GFF3 tag parsing failed on '{tv}' in '{groupfield}': the reference program loops forever on such a group field.");
            string tmptag = m.Groups[1].Value;
            string val = Trim4(m.Groups[2].Value);
            tv = Post(m, tv);
            if (tmptag is "Sequence" or "ID")
            {
                GroupSet(gf, val, false);
                grpset = true;
            }
            else if (tmptag is "FeatFlags" or "featflags")
            {
                FlagSet(gf, val, false);
            }
            else
            {
                if (escapedquote)
                    val = val.Replace("\u0001", "\"", StringComparison.Ordinal);
                TagAdd(gf, tmptag, QuoteStripValue(val));
            }
        }

        if (!grpset)
            GroupSet(gf, null, false);
    }

    private static int FrameOf(string token) => token switch { "0" => 1, "1" => 2, "2" => 3, _ => 0 };

    private static char StrandOf(string token) => token switch { "+" => '+', "-" => '-', _ => '\0' };

    // the columns of a GFF line up to the group field; null = a missing column (goto Error)
    private bool GffColumns(string line, out string seqid, out string source, out string feature, out int start, out int end, out float score, out char strand, out int frame, out StrTok? tok)
    {
        seqid = source = feature = "";
        start = end = frame = 0;
        score = 0f;
        strand = '\0';
        tok = new StrTok(line, "\t");
        if (!tok.Next(out seqid) || !tok.Next(out source) || !tok.Next(out feature))
            return false;
        if (!tok.Next(out string t))
            return false;
        if (!CStr.ToInt(t, out start))
            start = 0;
        if (!tok.Next(out t))
            return false;
        if (!CStr.ToInt(t, out end))
            end = 0;
        if (!tok.Next(out t))
            return false;
        if (!CStr.ToFloat(t, out score))
            score = 0f;
        if (!tok.Next(out t))
            return false;
        strand = StrandOf(t);
        if (!tok.Next(out t))
            return false;
        frame = FrameOf(t);
        return true;
    }

    // featGff2FromLine
    private EmbossSeqFeature? Gff2FromLine(string line, float version)
    {
        if (line.Length == 0)
            return null;
        if (!GffColumns(line, out _, out string source, out string feature, out int start, out int end, out float score, out char strand, out int frame, out var tok))
            return null;
        if (strand == '-' && start < end)
            (start, end) = (end, start);
        var gf = Type == "P" ? NewProt(source, feature, start, end, score, 0, null) : NewNuc(source, feature, start, end, score, strand, frame, 0, 0, 0, null, null, 0, null);
        if (gf is not null && tok!.Rest(out string group))
            Gff2ProcessTagval(gf, group, version);
        return gf;
    }

    // featGff3oldFromLine
    private EmbossSeqFeature? Gff3oldFromLine(string line)
    {
        if (line.Length == 0)
            return null;
        if (!GffColumns(line, out _, out string source, out string feature, out int start, out int end, out float score, out char strand, out int frame, out var tok))
            return null;
        var gf = Type == "P" ? NewProt(source, feature, start, end, score, 0, null) : NewNuc(source, feature, start, end, score, strand, frame, 0, 0, 0, null, null, 0, null);
        if (gf is not null && tok!.Rest(out string group))
            Gff3oldProcessTagval(gf, group);
        return gf;
    }

    // featGff3FromLine
    private EmbossSeqFeature? Gff3FromLine(string line, Dictionary<string, EmbossSeqFeature> idtable, List<EmbossSeqFeature> childlist)
    {
        if (line.Length == 0)
            return null;
        if (!GffColumns(line, out string seqid, out string source, out string feature, out int start, out int end, out float score, out char strand, out int frame, out var tok))
            return null;
        if (SeqId.Length == 0)
            SeqId = seqid;
        var gf = Type == "P" ? NewProt(source, feature, start, end, score, 0, null) : NewNuc(source, feature, start, end, score, strand, frame, 0, 0, 0, null, null, 0, null);
        if (gf is null)
            return null;
        if (seqid != SeqId)
        {
            gf.Remote = seqid;
            gf.Flags |= EmbossLocationBits.RemoteId;
        }

        if (tok!.Rest(out string group))
        {
            string? idstr = Gff3ProcessTagval(gf, group, out bool hasparent);
            if (idstr is not null)
            {
                if (idtable.TryGetValue(idstr, out var gfknown))
                {
                    EmbossSeqFeature gftop;
                    if (gfknown.Subfeatures is { Count: > 0 })
                    {
                        gftop = gfknown;
                    }
                    else
                    {
                        gftop = gfknown;
                        gfknown = NewFeat(gftop);
                        gftop.Flags |= EmbossLocationBits.Generated;
                        gftop.Flags |= EmbossLocationBits.Multiple;
                        gftop.Subfeatures ??= [];
                        GfftagAdd(gftop, "ID", idstr + ".emboss");
                        gftop.Subfeatures.Add(gfknown);
                        idtable[idstr] = gftop;
                    }

                    if (gf.Start < gftop.Start)
                        gftop.Start = gf.Start;
                    if (gf.End > gftop.End)
                        gftop.End = gf.End;
                    PopLast(gf);
                    gftop.Subfeatures!.Add(gf);
                    foreach (var (tag, value) in gf.GffTags)
                    {
                        if (tag is "ID" or "Parent")
                            continue;
                        GfftagAdd(gftop, tag, value);
                    }
                }
                else
                {
                    idtable[idstr] = gf;
                    if (hasparent)
                    {
                        PopLast(gf);
                        childlist.Add(gf);
                    }
                }
            }
            else if (hasparent)
            {
                PopLast(gf);
                childlist.Add(gf);
            }
        }

        return gf;
    }

    // ajListPopLast(thys->Features, &gf)
    private void PopLast(EmbossSeqFeature gf)
    {
        if (Features.Count > 0)
            Features.RemoveAt(Features.Count - 1);
        _ = gf;
    }

    /// <summary><c>ajFeatNewFeat</c>: a copy of a feature; the tags are added again, so they are checked again.</summary>
    internal EmbossSeqFeature NewFeat(EmbossSeqFeature orig)
    {
        var ret = new EmbossSeqFeature
        {
            Source = orig.Source,
            Type = orig.Type,
            Remote = orig.Remote,
            Label = orig.Label,
            Protein = orig.Protein,
            Start = orig.Start,
            End = orig.End,
            Start2 = orig.Start2,
            End2 = orig.End2,
            Score = orig.Score,
            Strand = orig.Strand,
            Frame = orig.Frame,
            Flags = orig.Flags,
            Group = orig.Group,
            Exon = orig.Exon,
        };
        if (orig.Subfeatures is not null)
        {
            ret.Subfeatures = [];
            foreach (var sub in orig.Subfeatures)
                ret.Subfeatures.Add(NewFeat(sub));
        }

        foreach (var (tag, value) in orig.Tags)
            TagAdd(ret, tag, value);
        foreach (var (tag, value) in orig.GffTags)
            GfftagAdd(ret, tag, value);
        return ret;
    }

    /// <summary><c>featReadGff2</c>.</summary>
    internal bool ReadGff2(EmbossFeatureLines file)
    {
        bool found = false;
        float version = 2.0f;
        string? raw;
        while ((raw = file.Next()) is not null)
        {
            string line = Trim4(raw);
            Match m;
            if (GffBlank.IsMatch(line))
            {
                version = 2.0f;
            }
            else if ((m = GffVersion.Match(line)).Success)
            {
                CStr.ToFloat(m.Groups[1].Value, out version);
            }
            else if ((m = GffRegion.Match(line)).Success)
            {
                SeqId = m.Groups[1].Value;
                CStr.ToUint(m.Groups[2].Value, out Start);
                CStr.ToUint(m.Groups[3].Value, out End);
            }
            else if ((m = GffType.Match(line)).Success)
            {
                if (string.Equals(m.Groups[1].Value, "Protein", StringComparison.OrdinalIgnoreCase))
                    SetType("P");
                else
                    SetType("N");
                SeqId = m.Groups[3].Value;
            }
            else if (GffComment.IsMatch(line))
            {
                version = 2.0f;
            }
            else if (Gff2FromLine(line, version) is not null)
            {
                found = true;
            }
        }

        return found;
    }

    /// <summary><c>featReadGff3old</c>.</summary>
    internal bool ReadGff3Old(EmbossFeatureLines file)
    {
        bool found = false;
        string? raw;
        while ((raw = file.Next()) is not null)
        {
            string line = Trim4(raw);
            Match m;
            if (GffBlank.IsMatch(line))
            {
                continue;
            }

            if ((m = GffVersion.Match(line)).Success)
            {
                CStr.ToFloat(m.Groups[1].Value, out float version);
                if (version < 3.0f)
                    return false;
            }
            else if ((m = GffRegion.Match(line)).Success)
            {
                SeqId = m.Groups[1].Value;
                CStr.ToUint(m.Groups[2].Value, out Start);
                CStr.ToUint(m.Groups[3].Value, out End);
            }
            else if (Gff3oldFromLine(line) is not null)
            {
                found = true;
            }
        }

        return found;
    }

    /// <summary><c>featReadGff3</c>: reading stops at a second <c>##gff-version</c> line, which is left for the next table.</summary>
    internal bool ReadGff3(EmbossFeatureLines file)
    {
        bool found = false, oldemboss = false, haveversion = false;
        var idtable = new Dictionary<string, EmbossSeqFeature>(StringComparer.Ordinal);
        var childlist = new List<EmbossSeqFeature>();
        string? raw;
        while ((raw = NextStored(file)) is not null)
        {
            string line = Trim4(raw);
            Match m;
            if (Gff3Directive.IsMatch(line))
            {
                if ((m = GffVersion.Match(line)).Success)
                {
                    if (haveversion)
                    {
                        file.Clear(1);
                        if (_context.StoreText && TextPtr is not null)
                            TextPtr = TextPtr[..Math.Max(0, TextPtr.Length - line.Length)];          // ajFilebuffClearStore cuts as many characters as the TRIMMED line has
                        Gff3Matchtable(childlist, idtable);
                        return found;
                    }

                    CStr.ToFloat(m.Groups[1].Value, out float version);
                    haveversion = true;
                    if (version < 3.0f)
                        return false;
                }
                else if ((m = GffRegion.Match(line)).Success)
                {
                    SeqId = m.Groups[1].Value;
                    CStr.ToUint(m.Groups[2].Value, out Start);
                    CStr.ToUint(m.Groups[3].Value, out End);
                }
            }
            else if (Gff3Comment.IsMatch(line))
            {
                if (line.StartsWith("#!Source-version EMBOSS 6.", StringComparison.Ordinal) && line.Length > 26 && line[26] is >= '0' and <= '3')
                    oldemboss = true;
                if (line.StartsWith("#!Type Protein", StringComparison.Ordinal))
                    SetType("P");
            }
            else if (oldemboss ? Gff3oldFromLine(line) is not null : Gff3FromLine(line, idtable, childlist) is not null)
            {
                found = true;
            }
        }

        Gff3Matchtable(childlist, idtable);
        return found;
    }

    // featGff3Matchtable: the children whose Parent is known become sub-features of it
    private void Gff3Matchtable(List<EmbossSeqFeature> childlist, Dictionary<string, EmbossSeqFeature> idtable)
    {
        foreach (var gf in childlist)
        {
            Features.Remove(gf);
            string? idparent = null;
            foreach (var (tag, value) in gf.GffTags)
            {
                if (tag == "Parent")
                {
                    idparent = value;
                    break;
                }
            }

            if (idparent is not null && idtable.TryGetValue(idparent, out var gfid))
            {
                if (gf.Type == gfid.Type)
                    gfid.Flags |= EmbossLocationBits.Multiple;
                gfid.Subfeatures ??= [];
                gfid.Subfeatures.Add(gf);
            }
            else
            {
                _context.Warn(idparent is not null ? $"featGff3Matchtable GFF3 Parent '{idparent}' no matching id" : $"Child feature with no parent {gf.Type} {(uint)gf.Start}..{(uint)gf.End}");
                Features.Add(gf);
            }
        }

        childlist.Clear();
        idtable.Clear();
    }
}
