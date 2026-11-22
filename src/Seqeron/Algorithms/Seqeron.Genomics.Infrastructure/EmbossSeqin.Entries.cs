#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589 // the readers keep the control flow of ajseqread.c

using System.Text;

namespace Seqeron.Genomics.Infrastructure;

// The database entry formats of ajseqread.c with features: embl, swiss, genbank (refseq), genpept, refseqp. The reference lists, cross references, keywords, comments, dates and taxonomy of an entry are not stored
// (no output of this repository shows them); the lines are read and skipped exactly as the library reads them, so the control flow, the name, accessions, version, description, sequence and features are those of the library.
internal sealed partial class EmbossSeqin
{
    private enum Sw { Unk, Ac, As, Av, Cc, Co, De, Dr, Dt, Ex, Fh, Ft, Gn, Id, Kw, Oc, Og, Oh, Os, Ox, Pe, Ra, Rc, Rg, Rl, Rn, Rp, Rt, Rx, Sq, Sv, Wp, Xx, End, More }

    // seqPrefixSwiss
    private static Sw PrefixSwiss(string s)
    {
        char a = s.Length > 0 ? s[0] : '\0';
        char b = s.Length > 1 ? s[1] : '\0';
        switch (a)
        {
            case 'A': return b switch { 'C' => Sw.Ac, 'H' => Sw.Fh, 'S' => Sw.As, 'V' => Sw.Av, _ => Sw.Unk };
            case 'C': return b switch { 'C' => Sw.Cc, 'O' => Sw.Co, _ => Sw.Unk };
            case 'D': return b switch { 'E' => Sw.De, 'R' => Sw.Dr, 'T' => Sw.Dt, _ => Sw.Unk };
            case 'E': return b is 'N' or 'X' ? Sw.Ex : Sw.Unk;
            case 'F': return b switch { 'H' => Sw.Fh, 'T' => Sw.Ft, _ => Sw.Unk };
            case 'G': return b == 'N' ? Sw.Gn : Sw.Unk;
            case 'I': return b switch { 'D' => Sw.Id, 'V' => Sw.Sv, _ => Sw.Unk };
            case 'K': return b == 'W' ? Sw.Kw : Sw.Unk;
            case 'O': return b switch { 'C' => Sw.Oc, 'G' => Sw.Og, 'H' => Sw.Oh, 'S' => Sw.Os, 'X' => Sw.Ox, _ => Sw.Unk };
            case 'P': return b switch { 'A' => Sw.Ac, 'E' => Sw.Pe, _ => Sw.Unk };
            case 'R': return b switch { 'A' => Sw.Ra, 'C' => Sw.Rc, 'G' => Sw.Rg, 'L' => Sw.Rl, 'N' => Sw.Rn, 'P' => Sw.Rp, 'T' => Sw.Rt, 'X' => Sw.Rx, _ => Sw.Unk };
            case 'S': return b switch { 'Q' => Sw.Sq, 'V' => Sw.Sv, _ => Sw.Unk };
            case 'T': return b == 'N' ? Sw.Ex : Sw.Unk;
            case 'W': return b == 'P' ? Sw.Wp : Sw.Unk;
            case 'X': return b == 'X' ? Sw.Xx : Sw.Unk;
            case '/': return b == '/' ? Sw.End : Sw.Unk;
            case ' ': return b == ' ' ? Sw.More : Sw.Unk;
            default: return Sw.Unk;
        }
    }

    private EmbossSeqFeatureTable? ReadFeatureLines(List<string> lines, string format, string seqname)
    {
        var table = new EmbossSeqFeatureTable(_context) { SeqId = seqname };
        table.Found = EmbossFeatureFormats.Read(table, EmbossFeatureFormats.Find(format), new EmbossFeatureLines(lines));
        table.FormatName = format;
        return table.Found ? table : null;
    }

    // a line as the feature readers take it (they read the lines without the line feed)
    private static string Chomp(string line) => line.EndsWith('\n') ? line[..^1] : line;

    private static bool PrefixC(string s, string p) => s.StartsWith(p, StringComparison.Ordinal);

    // ajStrPrefixCaseC
    private static bool PrefixCase(string s, string p) => s.StartsWith(p, StringComparison.OrdinalIgnoreCase);

    // ---------------------------------------------------------------------------------------------------------------- embl

    // seqReadEmbl (the lines of the references, cross-references, keywords, comments and dates are read and not stored; a CO (contig) entry is not supported: its pieces are read from databases)
    private bool ReadEmbl(EmbossSeq t)
    {
        int seqlen = 1024;
        if (!ReadLine())
            return false;
        Sw lineprefix = PrefixSwiss(_line);
        while (lineprefix == Sw.Wp)
        {
            if (!ReadLine())
                return false;
            lineprefix = PrefixSwiss(_line);
        }

        while (IsWhite(_line))
        {
            if (!ReadLine())
                return false;
            lineprefix = PrefixSwiss(_line);
        }

        if (lineprefix != Sw.Id)
        {
            StoreReset();
            return false;
        }

        AssignText();
        var h = new StrTok(_line, " ;\t\n\r");
        h.Next(out _);
        h.Next(out string token);
        SetName(t, token);
        h.Next(out token);
        if (token == "SV")
        {
            h.Next(out token);
            SvSave(t, t.Name + "." + token);
            h.Next(out token);
            if (token == "circular")
                t.Circular = true;
            h.Next(out _, ";\t\n\r");
            h.Next(out _);
            h.Next(out _);
            h.Next(out token);
            token = token.TrimEnd('B', 'P', '.').Trim(' ', '\t', '\n', '\r');
            if (CStr.ToUint(token, out uint n))
                seqlen = (int)n;
        }
        else
        {
            if (_line.Contains(" PRT; ", StringComparison.Ordinal) || _line.Contains(" Unreviewed; ", StringComparison.Ordinal) || _line.Contains(" Reviewed; ", StringComparison.Ordinal) || _line.Contains(" Preliminary; ", StringComparison.Ordinal))
            {
                StoreReset();
                return false;
            }
        }

        bool ok = ReadLine();
        lineprefix = PrefixSwiss(_line);
        bool dofeat = false;
        bool tryfeat = _opt.Features && !_opt.UfoGiven;
        List<string>? featureLines = null;
        var desc = new StringBuilder();
        var constr = new StringBuilder();
        while (ok && lineprefix != Sw.Sq && lineprefix != Sw.End)
        {
            if (lineprefix == Sw.Ex)
            {
                StoreReset();
                return false;
            }
            else if (lineprefix == Sw.Ac)
            {
                var a = new StrTok(_line, " ;\n\r");
                a.Next(out _);
                while (a.Next(out token))
                    AccSave(t, token);
            }
            else if (lineprefix == Sw.Sv)
            {
                var a = new StrTok(_line, " \n\r");
                a.Next(out _);
                a.Next(out token);
                SvSave(t, token);
            }
            else if (lineprefix == Sw.De)
            {
                var a = new StrTok(_line, " ");
                a.Next(out _);
                a.Next(out token, "\n\r");
                if (desc.Length > 0)
                    desc.Append(' ');
                desc.Append(token);
            }
            else if (tryfeat && lineprefix == Sw.Ft)
            {
                dofeat = true;
                (featureLines ??= []).Add(Chomp(_line));
            }
            else if (lineprefix == Sw.Co)
            {
                constr.Append(_line.TrimEnd(' ', '\t', '\n', '\r').Length > 5 ? _line.TrimEnd(' ', '\t', '\n', '\r')[5..] : "");
            }

            ok = ReadLine();
            lineprefix = PrefixSwiss(_line);
        }

        t.Desc = desc.ToString();
        if (lineprefix == Sw.End && constr.Length > 0)
            throw new InvalidOperationException("An EMBL entry with a CO (contig) line takes its sequence from other entries of a database; this port reads no database for it.");
        if (ok && lineprefix == Sw.Sq)
        {
            var a = new StrTok(_line, " ");
            a.Next(out _);
            a.Next(out _);
            a.Next(out token);
            if (CStr.ToUint(token, out uint tmplen) && (int)tmplen > seqlen)
                seqlen = (int)tmplen;
            if (dofeat)
            {
                t.Fttable = ReadFeatureLines(featureLines!, "embl", t.Name);
            }
            else if (tryfeat)
            {
                t.Fttable = new EmbossSeqFeatureTable(_context) { Found = true, SeqId = t.Name, Type = "P", TypeNull = false };      // ajFeattableNewSeq of a sequence with no residues yet: ajSeqIsProt is true
            }

            ok = ReadLine();
            lineprefix = PrefixSwiss(_line);
            while (ok && lineprefix != Sw.End)
            {
                Append(t, _line);
                ok = ReadLine();
                lineprefix = PrefixSwiss(_line);
            }
        }

        if (!IsNucleotideSeq(t))
        {
            StoreReset();
            return false;
        }

        t.SetNuc();
        if (t.Fttable is { } table)
            table.Length = t.Seq.Length;
        _buff.Clear(0);
        return true;
    }

    // ajSeqIsNuc: the type when the reader set one, else the letters (ajSeqTypeGapnucS)
    private static bool IsNucleotideSeq(EmbossSeq t)
    {
        if (t.Type == "N")
            return true;
        if (t.Type == "P")
            return false;
        return IsNucleotideText(t.SeqText);
    }

    // ---------------------------------------------------------------------------------------------------------------- swiss

    private sealed class SubDesc
    {
        public StringBuilder Name = new();
        public readonly List<StringBuilder> Short = [], Ec = [], Inn = [], Allergen = [], Biotech = [], Cdantigen = [];
    }

    private sealed class Desc
    {
        public StringBuilder Name = new();
        public readonly List<StringBuilder> Short = [], Ec = [];
        public readonly List<SubDesc> AltNames = [], SubNames = [];
        public readonly List<Desc> Includes = [], Contains = [];
        public bool Precursor;
        public int Fragments;
    }

    // seqDesSwiss: 0 unknown, 1 RecName, 2 AltName, 3 SubName, 4 Includes, 5 Contains, 6 Flags
    private static int DesSwiss(string s) => s switch { "AltName:" => 2, "Contains:" => 5, "Flags:" => 6, "Includes:" => 4, "RecName:" => 1, "SubName:" => 3, _ => 0 };

    // seqDessubSwiss: 0 unknown, 1 Full, 2 Short, 3 EC, 4 Allergen, 5 Biotech, 6 CD_antigen, 7 INN; the prefix is cut from the token
    private static int DessubSwiss(ref string s)
    {
        (string prefix, int code)[] subs = [("Allergen=", 4), ("Biotech=", 5), ("CD_antigen=", 6), ("EC=", 3), ("Full=", 1), ("INN=", 7), ("Short=", 2)];
        foreach (var (prefix, code) in subs)
        {
            if (s.StartsWith(prefix, StringComparison.Ordinal))
            {
                s = s[prefix.Length..];
                return code;
            }
        }

        return 0;
    }

    private static void Paren(StringBuilder desc, StringBuilder? item) => desc.Append(" (").Append(item).Append(')');

    // seqReadSwiss (the full description of the entry, DE lines of the RecName / AltName / SubName / Includes / Contains / Flags form; the gene, reference, comment, cross-reference and date lines are read and not stored)
    private bool ReadSwiss(EmbossSeq t)
    {
        if (!ReadLine())
            return false;
        Sw lineprefix = PrefixSwiss(_line);
        while (lineprefix == Sw.Wp)
        {
            if (!ReadLine())
                return false;
            lineprefix = PrefixSwiss(_line);
        }

        while (IsWhite(_line))
        {
            if (!ReadLine())
                return false;
            lineprefix = PrefixSwiss(_line);
        }

        if (lineprefix != Sw.Id)
        {
            StoreReset();
            return false;
        }

        var h = new StrTok(_line, " \n\r");
        h.Next(out _);
        h.Next(out string token);
        SetName(t, token);
        bool ok = ReadLine();
        lineprefix = PrefixSwiss(_line);
        var full = new Desc();
        t.FullDescMulti = [];
        var desc = new StringBuilder();
        bool dofeat = false;
        bool tryfeat = _opt.Features && !_opt.UfoGiven;
        List<string>? featureLines = null;
        Desc? desctop = null;
        Desc descmaster = full;
        SubDesc? subdesc = null;
        bool descistop = true;
        bool isdescflag = false;
        StringBuilder pdesc = desc;
        while (ok && lineprefix != Sw.Sq)
        {
            if (lineprefix == Sw.Ex)
            {
                StoreReset();
                return false;
            }
            else if (lineprefix == Sw.Ac)
            {
                var a = new StrTok(_line, " ;\n\r");
                a.Next(out _);
                while (a.Next(out token))
                    AccSave(t, token);
            }

            if (tryfeat && lineprefix == Sw.Ft)
            {
                dofeat = true;
                (featureLines ??= []).Add(Chomp(_line));
            }

            switch (lineprefix)
            {
                case Sw.De:
                    {
                        if (desctop is null)
                        {
                            desctop = full;
                            descmaster = full;
                            pdesc = desc;
                        }

                        var a = new StrTok(_line, " ");
                        a.Next(out _);
                        while (a.Next(out token, " ;\n\r"))
                        {
                            if (token.Length > 0 && token[^1] == ':')
                            {
                                isdescflag = false;
                                switch (DesSwiss(token))
                                {
                                    case 1:
                                        pdesc = descmaster.Name;
                                        descistop = true;
                                        break;
                                    case 2:
                                        subdesc = new SubDesc();
                                        descistop = false;
                                        pdesc = subdesc.Name;
                                        descmaster.AltNames.Add(subdesc);
                                        break;
                                    case 3:
                                        subdesc = new SubDesc();
                                        descistop = false;
                                        pdesc = subdesc.Name;
                                        descmaster.SubNames.Add(subdesc);
                                        break;
                                    case 4:
                                        descmaster = new Desc();
                                        descistop = true;
                                        full.Includes.Add(descmaster);
                                        pdesc = descmaster.Name;
                                        break;
                                    case 5:
                                        descmaster = new Desc();
                                        descistop = true;
                                        full.Contains.Add(descmaster);
                                        pdesc = descmaster.Name;
                                        break;
                                    case 6:
                                        isdescflag = true;
                                        break;
                                    default:
                                        if (pdesc.Length > 0)
                                            pdesc.Append(' ');
                                        pdesc.Append(token);
                                        break;
                                }
                            }
                            else if (token.Length > 1 && token.IndexOf('=', 1) >= 1)
                            {
                                switch (DessubSwiss(ref token))
                                {
                                    case 1:
                                        pdesc = descistop ? descmaster.Name : subdesc!.Name;
                                        pdesc.Clear().Append(token);
                                        break;
                                    case 2:
                                        pdesc = new StringBuilder(token);
                                        (descistop ? descmaster.Short : subdesc!.Short).Add(pdesc);
                                        break;
                                    case 3:
                                        pdesc = new StringBuilder(token);
                                        (descistop ? descmaster.Ec : subdesc!.Ec).Add(pdesc);
                                        break;
                                    case 4:
                                        pdesc = new StringBuilder(token);
                                        subdesc!.Allergen.Add(pdesc);
                                        break;
                                    case 5:
                                        pdesc = new StringBuilder(token);
                                        subdesc!.Biotech.Add(pdesc);
                                        break;
                                    case 6:
                                        pdesc = new StringBuilder(token);
                                        subdesc!.Cdantigen.Add(pdesc);
                                        break;
                                    case 7:
                                        pdesc = new StringBuilder(token.Length > 4 ? token[4..] : "");
                                        subdesc!.Inn.Add(pdesc);
                                        break;
                                    default:
                                        if (pdesc.Length > 0)
                                            pdesc.Append(' ');
                                        pdesc.Append(token);
                                        break;
                                }
                            }
                            else if (isdescflag)
                            {
                                if (token == "Precursor")
                                    full.Precursor = true;
                                else if (token == "Fragments")
                                    full.Fragments = 2;
                                else if (token == "Fragment")
                                    full.Fragments = 1;
                                else
                                {
                                    if (pdesc.Length > 0)
                                        pdesc.Append(' ');
                                    pdesc.Append("Flags: ").Append(token);
                                }
                            }
                            else
                            {
                                if (pdesc.Length > 0)
                                    pdesc.Append(' ');
                                pdesc.Append(token);
                            }
                        }

                        break;
                    }

                case Sw.Unk:
                case Sw.End:
                case Sw.More:
                case Sw.Xx:
                case Sw.Sv:
                    _context.Warn($"Unknown swissprot line type '{(_line.Length >= 2 ? _line[..2] : _line.PadLeft(2))}'");
                    break;
            }

            ok = ReadLine();
            lineprefix = PrefixSwiss(_line);
        }

        if (dofeat)
            t.Fttable = ReadFeatureLines(featureLines!, "swissprot", t.Name);
        ok = ReadLine();
        lineprefix = PrefixSwiss(_line);
        while (ok && lineprefix != Sw.End)
        {
            Append(t, _line);
            ok = ReadLine();
            lineprefix = PrefixSwiss(_line);
        }

        if (t.Fttable is { } table)
            table.Length = t.Seq.Length;
        if (desc.Length == 0)
        {
            desc.Append(full.Name);
            foreach (var s in full.Short)
                Paren(desc, s);
            foreach (var s in full.Ec)
                Paren(desc, s);
            foreach (var sd in full.AltNames)
            {
                if (sd.Name.Length > 0)
                    Paren(desc, sd.Name);
                foreach (var l in new[] { sd.Inn, sd.Short, sd.Ec, sd.Allergen, sd.Biotech, sd.Cdantigen })
                    foreach (var s in l)
                        Paren(desc, s);
            }

            foreach (var sd in full.SubNames)
            {
                Paren(desc, sd.Name);
                foreach (var l in new[] { sd.Short, sd.Ec })
                    foreach (var s in l)
                        Paren(desc, s);
            }

            foreach (var d in full.Includes)
                IncludedDesc(desc, d, false);
            foreach (var d in full.Contains)
                IncludedDesc(desc, d, true);
            if (full.Fragments != 0 || full.Precursor)
            {
                if (full.Fragments == 1)
                    desc.Append(" (Fragment)");
                if (full.Fragments == 2)
                    desc.Append(" (Fragments)");
                if (full.Precursor)
                    desc.Append(" (Precursor)");
            }

            if (desc.Length > 0 && desc[0] == ' ')
                desc.Remove(0, 1);
        }

        t.Desc = desc.ToString();
        t.SetProt();
        _buff.Clear(0);
        return true;
    }

    // the Includes / Contains blocks of the full description (the Contains block of the library prints the CD antigen names of its SubNames twice instead of their short and EC names)
    private static void IncludedDesc(StringBuilder desc, Desc d, bool contains)
    {
        Paren(desc, d.Name);
        foreach (var s in d.Short)
            Paren(desc, s);
        foreach (var s in d.Ec)
            Paren(desc, s);
        foreach (var sd in d.AltNames)
        {
            if (sd.Name.Length > 0)
                Paren(desc, sd.Name);
            foreach (var l in new[] { sd.Inn, sd.Short, sd.Ec, sd.Allergen, sd.Biotech, sd.Cdantigen })
                foreach (var s in l)
                    Paren(desc, s);
        }

        foreach (var sd in d.SubNames)
        {
            Paren(desc, sd.Name);
            if (contains)
            {
                foreach (var s in sd.Cdantigen)
                    Paren(desc, s);
                foreach (var s in sd.Cdantigen)
                    Paren(desc, s);
            }
            else
            {
                foreach (var s in sd.Short)
                    Paren(desc, s);
                foreach (var s in sd.Ec)
                    Paren(desc, s);
            }
        }
    }

    // ---------------------------------------------------------------------------------------------------------------- genbank, refseq

    private enum Gb { Unk, Ac, Base, Cc, Def, Feat, Key, Id, Ori, Ref, Src, Ver, Wp, End, More }

    // seqPrefixGenbank
    private static Gb PrefixGenbank(string s)
    {
        if (PrefixC(s, "ACCESSION")) return Gb.Ac;
        if (PrefixC(s, "BASE COUNT")) return Gb.Base;
        if (PrefixC(s, "COMMENT")) return Gb.Cc;
        if (PrefixC(s, "DEFINITION")) return Gb.Def;
        if (PrefixC(s, "FEATURES")) return Gb.Feat;
        if (PrefixC(s, "KEYWORDS")) return Gb.Key;
        if (PrefixC(s, "LOCUS")) return Gb.Id;
        if (PrefixC(s, "ORIGIN")) return Gb.Ori;
        if (PrefixC(s, "REFERENCE")) return Gb.Ref;
        if (PrefixC(s, "SEQVERSION")) return Gb.Ver;
        if (PrefixC(s, "SOURCE")) return Gb.Src;
        if (PrefixC(s, "WPCOMMENT")) return Gb.Wp;
        if (PrefixC(s, "//")) return Gb.End;
        if (PrefixC(s, "  ")) return Gb.More;
        return Gb.Unk;
    }

    // seqPrefixGenbankMore: only a line of nine blanks is GB_MORE_MORE (the library compares from the SECOND character, so the keywords AUTHORS, TITLE, JOURNAL, ORGANISM never match)
    private static bool IsMoreMore(string line) => PrefixC(line, "         ");

    // seqReadGenbank (refseq is the same reader)
    private bool ReadGenbank(EmbossSeq t)
    {
        bool done;
        int seqlen = 1024;
        if (!NextRaw())
            return false;
        Gb lineprefix = PrefixGenbank(_line);
        bool ok = true;
        while (IsWhite(_line))
        {
            if (!ReadLine())
                return false;
            lineprefix = PrefixGenbank(_line);
        }

        if (lineprefix == Gb.Wp)
        {
            ok = NextRaw();
            lineprefix = PrefixGenbank(_line);
            while (ok && lineprefix == Gb.More)
            {
                ok = NextRaw();
                lineprefix = PrefixGenbank(_line);
            }
        }

        if (_line.Contains("Genetic Sequence Data Bank", StringComparison.Ordinal))
        {
            while (ok && lineprefix != Gb.Id)
            {
                ok = NextRaw();
                lineprefix = PrefixGenbank(_line);
            }
        }

        if (!ok || lineprefix != Gb.Id)
        {
            StoreReset();
            return false;
        }

        var tokens = new StrTok(_line, " \n\r");
        int nfields = 0;
        while (tokens.Next(out _))
            nfields++;
        if (nfields == 9)
        {
            _buff.SetBuffered();                                    // ajFilebuffSetBuffered + ajTextinStoreReset: the reader of the other format starts from the first buffered line
            StoreReset();
            return ReadProteinEntry(t, false);
        }

        var h = new StrTok(_line, " \n\r");
        int i = 0;
        while (h.Next(out string token))
        {
            switch (++i)
            {
                case 2:
                    SetName(t, token);
                    break;
                case 3:
                    if (CStr.ToUint(token, out uint n))
                        seqlen = (int)n;
                    break;
                case 4:
                    if (token == "aa")
                    {
                        _buff.SetBuffered();
                        StoreReset();
                        return ReadProteinEntry(t, true);
                    }

                    if (token != "bp")
                        _context.Warn($"bad Genbank LOCUS line '{_line}'");
                    break;
                case 6:
                    if (token == "circular")
                        t.Circular = true;
                    break;
            }
        }

        AssignText();
        ok = ReadLine();
        lineprefix = PrefixGenbank(_line);
        bool dofeat = false;
        bool tryfeat = _opt.Features && !_opt.UfoGiven;
        List<string>? featureLines = null;
        var desc = new StringBuilder();
        while (ok && lineprefix != Gb.End && lineprefix != Gb.Ori && lineprefix != Gb.Base)
        {
            done = false;
            if (lineprefix == Gb.Def)
            {
                var a = new StrTok(_line, " ");
                a.Next(out _);
                a.Next(out string d, "\n\r");
                desc.Clear().Append(d);
                ok = ReadLine();
                lineprefix = PrefixGenbank(_line);
                done = true;
                while (ok && lineprefix == Gb.More)
                {
                    var c = new StrTok(_line, " ");
                    c.Next(out string more, "\n\r");
                    desc.Append(' ').Append(more);
                    ok = ReadLine();
                    lineprefix = PrefixGenbank(_line);
                }
            }
            else if (lineprefix == Gb.Ac)
            {
                var a = new StrTok(_line, " ;\n\r");
                a.Next(out _);
                while (a.Next(out string acc))
                    AccSave(t, acc);
            }
            else if (lineprefix == Gb.Ver)
            {
                var a = new StrTok(_line, " \n\r");
                a.Next(out _);
                a.Next(out string sv);
                SvSave(t, sv);
                if (a.Step(": \n\r"))
                {
                    a.Next(out string gi);
                    t.Gi = gi;
                }
            }
            else if (lineprefix == Gb.Src)
            {
                ok = ReadLine();
                lineprefix = PrefixGenbank(_line);
                done = true;
                while (ok && lineprefix == Gb.More)
                {
                    done = false;
                    ok = ReadLine();                                // the line after the lines of a SOURCE block is read again by the loop (a FEATURES line there is lost)
                    lineprefix = PrefixGenbank(_line);
                }
            }
            else if (tryfeat && lineprefix == Gb.Feat)
            {
                dofeat = true;
                (featureLines ??= []).Add(Chomp(_line));
                ok = ReadLine();
                lineprefix = PrefixGenbank(_line);
                done = true;
                while (ok && lineprefix == Gb.More)
                {
                    featureLines.Add(Chomp(_line));
                    ok = ReadLine();
                    lineprefix = PrefixGenbank(_line);
                }
            }
            else if (lineprefix == Gb.Ref || (ok && lineprefix == Gb.Cc))
            {
                ok = ReadLine();
                done = true;
                while (ok && IsMoreMore(_line))
                    ok = ReadLine();
                lineprefix = PrefixGenbank(_line);
            }
            else if (lineprefix == Gb.Key)
            {
                ok = ReadLine();
                lineprefix = PrefixGenbank(_line);
                done = true;
                while (ok && lineprefix == Gb.More)
                {
                    ok = ReadLine();
                    lineprefix = PrefixGenbank(_line);
                }
            }

            if (!done)
            {
                ok = ReadLine();
                lineprefix = PrefixGenbank(_line);
            }
        }

        t.Desc = desc.ToString();
        if (dofeat)
            t.Fttable = ReadFeatureLines(featureLines!, "genbank", t.Name);
        while (ok && lineprefix != Gb.End && lineprefix != Gb.Ori && lineprefix != Gb.Base)
        {
            ok = ReadLine();
            lineprefix = PrefixGenbank(_line);
            if (!ok)
                break;
        }

        if (ok && lineprefix != Gb.End)
        {
            ok = ReadLine();
            lineprefix = PrefixGenbank(_line);
        }

        while (ok && lineprefix != Gb.End)
        {
            if (lineprefix != Gb.Ori && lineprefix != Gb.Base)
                Append(t, _line);
            ok = ReadLine();
            lineprefix = PrefixGenbank(_line);
        }

        _ = seqlen;
        while (ok && lineprefix != Gb.End)
        {
            ok = ReadLine();
            lineprefix = PrefixGenbank(_line);
        }

        if (t.Fttable is { } table)
            table.Length = t.Seq.Length;
        _buff.Clear(0);
        return true;
    }

    // ajBuffreadLine: a line is read without counting it as a record
    private bool NextRaw()
    {
        string? l = _buff.Next();
        if (l is null)
        {
            _line = "";
            return false;
        }

        _line = l;
        return true;
    }

    // seqReadGenpept (refseqp = false) and seqReadRefseqp (true): the protein entries of GenBank format
    private bool ReadProteinEntry(EmbossSeq t, bool refseqp)
    {
        int seqlen = 1024;
        if (!NextRaw())
            return false;
        bool ok = true;
        while (IsWhite(_line))
        {
            if (!ReadLine())
                return false;
        }

        if (PrefixC(_line, "WPCOMMENT"))
        {
            ok = NextRaw();
            while (ok && PrefixC(_line, " "))
                ok = NextRaw();
        }

        if (!refseqp && _line.Contains("Genetic Sequence Data Bank", StringComparison.Ordinal))
        {
            while (ok && !PrefixC(_line, "LOCUS"))
                ok = NextRaw();
        }

        if (!ok || !PrefixC(_line, "LOCUS"))
        {
            StoreReset();
            return false;
        }

        if (refseqp)
            AssignText();
        if (!refseqp)
        {
            var count = new StrTok(_line, " \n\r");
            int nfields = 0;
            while (count.Next(out _))
                nfields++;
            if (nfields == 8)
            {
                _buff.SetBuffered();
                StoreReset();
                return ReadProteinEntry(t, true);
            }
        }

        var h = new StrTok(_line, " \n\r");
        int i = 0;
        while (h.Next(out string token))
        {
            switch (++i)
            {
                case 2:
                    SetName(t, token);
                    break;
                case 3:
                    if (CStr.ToUint(token, out uint n))
                        seqlen = (int)n;
                    break;
                case 4:
                    if (token != "aa")
                        _context.Warn($"bad {(refseqp ? "RefseqP" : "Genpept")} LOCUS line '{_line}'");
                    break;
            }
        }

        _ = seqlen;
        if (!refseqp)
            AssignText();
        ok = ReadLine();
        bool dofeat = false;
        bool tryfeat = refseqp && _opt.Features && !_opt.UfoGiven;
        List<string>? featureLines = null;
        var desc = new StringBuilder();
        while (ok && !PrefixC(_line, "ORIGIN") && !PrefixC(_line, "BASE COUNT"))
        {
            bool done = false;
            if (PrefixC(_line, "DEFINITION"))
            {
                var a = new StrTok(_line, " ");
                a.Next(out _);
                a.Next(out string d, "\n\r");
                desc.Clear().Append(d);
                ok = ReadLine();
                done = true;
                while (ok && PrefixC(_line, " "))
                {
                    var c = new StrTok(_line, " ");
                    c.Next(out string more, "\n\r");
                    desc.Append(' ').Append(more);
                    ok = ReadLine();
                }
            }
            else if (PrefixC(_line, "ACCESSION"))
            {
                var a = new StrTok(_line, " ;\n\r");
                a.Next(out _);
                while (a.Next(out string acc))
                    AccSave(t, acc);
            }
            else if (PrefixC(_line, "VERSION"))
            {
                var a = new StrTok(_line, " \n\r");
                a.Next(out _);
                a.Next(out string sv);
                SvSave(t, sv);
                if (a.Step(": \n\r"))
                {
                    a.Next(out string gi);
                    t.Gi = gi;
                }
            }
            else if (PrefixC(_line, "SOURCE"))
            {
                ok = ReadLine();
                done = true;
                while (ok && PrefixC(_line, " "))
                {
                    done = false;
                    if (PrefixC(_line, "  ORGANISM"))
                    {
                        ok = ReadLine();
                        done = true;
                        while (ok && PrefixC(_line, "    "))
                            ok = ReadLine();
                    }

                    if (!done)
                        ok = ReadLine();
                }
            }
            else if (tryfeat && PrefixC(_line, "FEATURES"))
            {
                dofeat = true;
                (featureLines ??= []).Add(Chomp(_line));
                ok = ReadLine();
                done = true;
                while (ok && PrefixC(_line, " "))
                {
                    featureLines.Add(Chomp(_line));
                    ok = ReadLine();
                }
            }
            else if (PrefixC(_line, "REFERENCE"))
            {
                ok = ReadLine();
                done = true;
                if (ok && PrefixC(_line, "  AUTHORS"))
                {
                    ok = ReadLine();
                    while (ok && PrefixC(_line, "          "))
                        ok = ReadLine();
                }

                if (ok && PrefixC(_line, "  TITLE"))
                {
                    ok = ReadLine();
                    while (ok && PrefixC(_line, "          "))
                        ok = ReadLine();
                }

                if (ok && PrefixC(_line, "  JOURNAL"))
                    ok = ReadLine();
                while (ok && PrefixC(_line, "  "))
                    ok = ReadLine();
            }
            else if (ok && PrefixC(_line, "COMMENT"))
            {
                ok = ReadLine();
                done = true;
                while (ok && PrefixC(_line, "          "))
                    ok = ReadLine();
            }
            else if (PrefixC(_line, "KEYWORDS"))
            {
                ok = ReadLine();
                done = true;
                while (ok && PrefixC(_line, " "))
                    ok = ReadLine();
            }
            else if (PrefixC(_line, "  ORGANISM"))
            {
                ok = ReadLine();
                done = true;
                while (ok && PrefixC(_line, "    "))
                    ok = ReadLine();
            }

            if (!done)
                ok = ReadLine();
        }

        t.Desc = desc.ToString();
        if (dofeat)
            t.Fttable = ReadFeatureLines(featureLines!, "refseqp", t.Name);
        while (!PrefixC(_line, "ORIGIN") && !PrefixC(_line, "BASE COUNT"))
        {
            if (!ReadLine())
                break;
        }

        ok = ReadLine();
        while (ok && !PrefixC(_line, "//"))
        {
            if (!PrefixC(_line, "ORIGIN") && !PrefixC(_line, "BASE COUNT"))
                Append(t, _line);
            ok = ReadLine();
        }

        while (ok && !PrefixC(_line, "//"))
            ok = ReadLine();
        if (t.Fttable is { } table)
            table.Length = t.Seq.Length;
        _buff.Clear(0);
        return true;
    }
}
