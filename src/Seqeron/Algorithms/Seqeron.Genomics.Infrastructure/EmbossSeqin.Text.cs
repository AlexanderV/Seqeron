#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589 // the readers keep the control flow of ajseqread.c

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Infrastructure;

// The sequence formats of one sequence per entry that need no tables: intelligenetics (ig, igstrict, iguspto), raw, fastq (fastq, fastq-sanger, fastq-illumina, fastq-solexa), sam, gff2 and gff3 with the sequence of the file.
internal sealed partial class EmbossSeqin
{
    // ajStrRemoveWhiteExcess
    internal static string RemoveWhiteExcess(string s)
    {
        if (s.Length == 0)
            return s;
        var chars = s.Replace('\t', ' ').Replace('\n', ' ').Replace('\r', ' ');
        var sb = new StringBuilder(chars.Length);
        bool lastSpace = true;
        foreach (char c in chars)
        {
            if (c == ' ')
            {
                if (!lastSpace)
                    sb.Append(' ');
                lastSpace = true;
            }
            else
            {
                sb.Append(c);
                lastSpace = false;
            }
        }

        if (sb.Length > 0 && sb[^1] == ' ')
            sb.Length--;
        return sb.ToString();
    }

    // ajStrExtractFirst: the first word (the string must not start with a blank) and the rest after the blanks that follow it
    private static bool ExtractFirst(string s, out string rest, out string word)
    {
        rest = word = "";
        if (s.Length == 0 || EmbossSeqFeatureTable.IsCSpace(s[0]))
            return false;
        int i = 1;
        while (i < s.Length && !EmbossSeqFeatureTable.IsCSpace(s[i]))
            i++;
        word = s[..i];
        int j = i;
        while (j < s.Length && EmbossSeqFeatureTable.IsCSpace(s[j]))
            j++;
        rest = j < s.Length ? s[j..] : "";
        return true;
    }

    // ---------------------------------------------------------------------------------------------------------------- ig

    // seqReadIg
    private bool ReadIg(EmbossSeq t)
    {
        bool ok;
        var desc = new StringBuilder();
        do
        {
            if (_records > 0)
            {
                string d = RemoveWhiteExcess(_line.Length > 0 ? _line[1..] : "");
                if (desc.Length > 0)
                    desc.Append(' ');
                desc.Append(d);
            }

            ok = ReadLine();
        }
        while (ok && PrefixC(_line, ";"));
        if (!ok)
        {
            StoreReset();
            return false;
        }

        _buff.SetUnbuffered();
        t.Desc = desc.ToString();
        SetName(t, _line);
        while (ReadLine() && !PrefixC(_line, ";"))
            Append(t, _line);
        if (PrefixC(_line, ";"))
            StoreClear(1);
        else
            _buff.Clear(0);
        return true;
    }

    // seqReadIgstrict
    private bool ReadIgstrict(EmbossSeq t)
    {
        bool ok;
        bool endnum = false;
        var desc = new StringBuilder();
        do
        {
            if (_records > 0)
            {
                string d = RemoveWhiteExcess(_line.Length > 0 ? _line[1..] : "");
                if (desc.Length > 0)
                    desc.Append(' ');
                desc.Append(d);
            }

            ok = ReadLine();
        }
        while (ok && PrefixC(_line, ";"));
        t.Desc = desc.ToString();
        if (!ok)
        {
            StoreReset();
            return false;
        }

        SetName(t, _line);
        while (ReadLine() && !PrefixC(_line, ";"))
        {
            _line = RemoveWhiteExcess(_line);
            endnum = _line.EndsWith('1') || _line.EndsWith('2');
            Append(t, _line);
        }

        if (!endnum)
        {
            StoreReset();
            return false;
        }

        if (PrefixC(_line, ";"))
            StoreClear(1);
        else
            _buff.Clear(0);
        return true;
    }

    // seqToken2 of the library (a static: the reason of the last bad first line of iguspto; the NULL string prints as <null>)
    private string _igToken2 = "<null>";

    // seqReadIguspto (the multi-line description of the US patent office sequences)
    private bool ReadIguspto(EmbossSeq t)
    {
        bool ok = true;
        bool seqok = false;
        bool isheader = true;
        bool firstline = true;
        var desc = new StringBuilder();
        while (ok && !seqok)
        {
            if (!PrefixC(_line, ";"))
            {
                StoreReset();
                return false;
            }

            if (firstline)
            {
                firstline = false;
                bool firstgood = true;
                if (!PrefixC(_line, "; Sequence "))
                    firstgood = false;
                int ipos = -1;
                if (firstgood)
                {
                    ipos = _line.IndexOf(", Application ", StringComparison.Ordinal);
                    if (ipos < 1)
                    {
                        _igToken2 = "', Application ' not found";
                        firstgood = false;
                    }
                }

                if (firstgood)
                {
                    string number = Sub(_line, 11, ipos - 1);
                    if (!IsInt(number))
                    {
                        _igToken2 = $"Sequence number '{number}' not an integer";
                        firstgood = false;
                    }
                }

                if (firstgood)
                {
                    string app = Sub(_line, ipos + 14, _line.Length - 2);
                    if (!IsWord(app))
                    {
                        _igToken2 = $"Application id '{app}' not a word";
                        firstgood = false;
                    }
                }

                if (!firstgood)
                {
                    if (_format != 0)
                        _context.Warn($"Iguspto: bad first line ({_igToken2}): {_line.TrimEnd(' ', '\t', '\n', '\r')}");
                    else
                        return false;
                }
            }

            do
            {
                if (_records > 0)
                {
                    string d = _line.EndsWith('\n') ? _line[..^1] : _line;
                    d = d.Length > 0 ? d[1..] : "";
                    if (d.Length > 0 && d[0] == ' ')
                        d = d[1..];
                    (t.FullDescMulti ??= []).Add(d);
                    if (PrefixC(d, "GENERAL INFORMATION"))
                        isheader = false;
                    if (isheader)
                    {
                        if (desc.Length > 0)
                            desc.Append("; ");
                        desc.Append(d);
                    }
                }

                ok = ReadLine();
            }
            while (ok && PrefixC(_line, ";"));
            if (!ok)
            {
                StoreReset();
                return false;
            }

            SetName(t, _line);
            bool endnum = false;
            ok = ReadLine();
            while (ok && !PrefixC(_line, ";") && !endnum)
            {
                _line = RemoveWhiteExcess(_line);
                endnum = _line.EndsWith('1') || _line.EndsWith('2');
                if (endnum)
                    _line = _line[..^1];
                AppendWarn(t, _line);
                ok = ReadLine();
            }

            if (endnum)
                seqok = true;
            else
            {
                if (_format != 0)
                    _context.Warn($"Sequence '{t.Name}' has bad iguspto sequence format");
                var fresh = new EmbossSeq();
                t.Name = fresh.Name;
                t.Seq.Clear();
                t.Desc = "";
                desc.Clear();
                t.Acc = t.Sv = "";
                t.Acclist.Clear();
            }
        }

        t.Desc = desc.ToString();
        if (ok)
        {
            _line = RemoveWhiteExcess(_line);
            while (ok && (_line.Length == 0 || _line == "\f"))
            {
                ok = NextRaw();
                if (ok)
                    _line = RemoveWhiteExcess(_line);
            }
        }

        if (PrefixC(_line, ";"))
            StoreClear(1);
        else
            _buff.Clear(0);
        return true;
    }

    // ajStrAssignSubS(&dest, str, begin, end): the characters from begin to end (both included, negative counts from the end); clipped
    private static string Sub(string s, int begin, int end)
    {
        if (end < 0)
            end += s.Length;
        if (begin < 0)
            begin += s.Length;
        if (begin < 0)
            begin = 0;
        if (end >= s.Length)
            end = s.Length - 1;
        if (begin > end || begin >= s.Length)
            return "";
        return s[begin..(end + 1)];
    }

    // ajStrIsInt
    private static bool IsInt(string s) => CStr.ToInt(s, out _);

    // ---------------------------------------------------------------------------------------------------------------- raw

    private bool _rawEnded;

    // seqReadRaw: the whole remaining file must consist of letters, digits, blanks and '*' '-'
    private bool ReadRaw(EmbossSeq t)
    {
        if (_rawEnded || _buff.IsEnded)
            return false;
        var sb = new StringBuilder();
        if (_opt.FullText is { } full)
        {
            // a regular file is read as a whole from its first byte (the offset of -ioffset only positions the line readers)
            _records = 0;
            foreach (char c in full)
            {
                if (!(char.IsAsciiLetterOrDigit(c) || c is ' ' or '\t' or '\n' or '\r' or '*' or '-'))
                {
                    StoreReset();
                    return false;
                }
            }

            sb.Append(full);
            _records = (full.Length + 2047) / 2048;
        }
        else
        {
            while (ReadLine())
            {
                foreach (char c in _line)
                {
                    if (!(char.IsAsciiLetterOrDigit(c) || c is ' ' or '\t' or '\n' or '\r' or '*' or '-'))
                    {
                        StoreReset();
                        return false;
                    }
                }

                sb.Append(_line);
            }
        }

        Append(t, sb.ToString());
        if (_opt.StoreText && (_opt.FullText ?? _fullText) is { } whole)
            _store.Clear().Append(whole);                           // the binary reader stores the chunks of the whole file, from its first byte
        _rawEnded = true;
        if (_records == 0)
        {
            StoreClear(-1);
            return false;
        }

        _buff.Clear(-1);
        return true;
    }

    // ---------------------------------------------------------------------------------------------------------------- fastq

    // ajSeqParseFastq: the id (the first word without the '@') and the description (the rest of the line)
    private static bool ParseFastq(string line, out string id, out string desc)
    {
        id = desc = "";
        if (!PrefixC(line, "@"))
            return false;
        _ = ExtractFirstWord(line, out string rest, out string word);
        desc = rest.Trim('\n');
        id = word.Length > 1 ? word[1..] : "";
        return true;
    }

    // ajStrExtractWord: the first word (leading blanks skipped) and the rest after the blanks that follow it
    private static bool ExtractFirstWord(string s, out string rest, out string word)
    {
        rest = word = "";
        if (s.Length == 0)
            return false;
        int i = 0;
        while (i < s.Length && EmbossSeqFeatureTable.IsCSpace(s[i]))
            i++;
        if (i >= s.Length)
            return false;
        int e = i;
        while (e < s.Length && !EmbossSeqFeatureTable.IsCSpace(s[e]))
            e++;
        word = s[i..e];
        int j = e;
        while (j < s.Length && EmbossSeqFeatureTable.IsCSpace(s[j]))
            j++;
        rest = j < s.Length ? s[j..] : "";
        return true;
    }

    // seqqualAppendWarn: the printable characters of the line
    private static void AppendQuality(StringBuilder qual, string line)
    {
        foreach (char c in line)
        {
            if (c is >= (char)33 and <= (char)126)
                qual.Append(c);
        }
    }

    // seqReadFastq (the quality values are read and not stored)
    private bool ReadFastq(EmbossSeq t)
    {
        bool ok = ReadLine();
        if (!ok)
            return false;
        if (_line.Length == 0 || _line[0] != '@')
        {
            StoreReset();
            return false;
        }

        if (!ParseFastq(_line, out string id, out string desc))
        {
            StoreReset();
            return false;
        }

        t.Name = NameNospace(id);
        t.Desc = desc;
        ok = ReadLine();
        while (ok && !PrefixC(_line, "+"))
        {
            AppendWarn(t, _line);
            ok = ReadLine();
        }

        if (!ok)
        {
            StoreReset();
            return false;
        }

        int seqlen = t.Seq.Length;
        ok = ReadLine();
        var qual = new StringBuilder();
        while (ok && (qual.Length < seqlen || _line.Length == 0 || _line[0] != '@'))
        {
            AppendQuality(qual, _line);
            ok = ReadLine();
        }

        if (qual.Length != seqlen)
        {
            StoreReset();
            return false;
        }

        if (ok)
            StoreClear(1);
        else
            _buff.Clear(0);
        return true;
    }

    // seqReadFastqSanger (sanger = true), seqReadFastqIllumina (qmin 64), seqReadFastqSolexa (qmin 59)
    private bool ReadFastqQuality(EmbossSeq t, string kind, int qmin)
    {
        const int qmax = 126;
        bool sanger = kind == "SANGER";
        bool ok = ReadLine();
        if (!ok)
            return false;
        string saveLine = _line;
        if (saveLine.Length == 0 || saveLine[0] != '@')
        {
            StoreReset();
            return false;
        }

        if (!ParseFastq(saveLine, out string id, out string desc))
        {
            StoreReset();
            return false;
        }

        t.Name = NameNospace(id);
        t.Desc = desc;
        int cntseq = 0, cntqual = 0, cntnewline = 0;
        if (sanger)
        {
            int i = saveLine.Length - 1;
            while (i >= 0 && (saveLine[i] == '\n' || saveLine[i] == '\r'))
            {
                cntnewline++;
                i--;
            }
        }

        ok = ReadLine();
        while (ok && (_line.Length == 0 || _line[0] != '+'))
        {
            if (sanger)
                cntseq += _line.Length - cntnewline;
            AppendWarn(t, _line);
            ok = ReadLine();
        }

        if (!ok)
        {
            StoreReset();
            return false;
        }

        if (sanger && _line.Length > cntnewline + 1)
        {
            string plus = "@" + _line[1..];
            if (plus != saveLine)
            {
                plus = "+" + plus[1..];
                _context.Warn($"Mismatch in file '{_opt.Filename}' + line does not match first line '{Prefix(saveLine, saveLine.Length - cntnewline)}' '{Prefix(plus, plus.Length - cntnewline)}'");
            }
        }

        int seqlen = t.Seq.Length;
        if (sanger && seqlen < cntseq)
            _context.Warn($"FASTQ format '{_opt.Filename}' sequence '{t.Name}' sequence skipped {cntseq - seqlen} character(s)");
        ok = ReadLine();
        var qual = new StringBuilder();
        string saveLine2 = "";
        while (ok && (qual.Length < seqlen || _line.Length == 0 || _line[0] != '@'))
        {
            if (sanger)
            {
                if (_line.Length > 0 && _line[0] == '@' && saveLine2.Length == 0)
                    saveLine2 = _line;
                cntqual += _line.Length - cntnewline;
            }

            AppendQuality(qual, _line);
            ok = ReadLine();
        }

        if (qual.Length != seqlen)
        {
            if (!sanger)
            {
                StoreReset();
                return false;
            }

            _context.Warn($"FASTQ quality length mismatch '{_opt.Filename}' '{t.Name}' expected: {seqlen} found: {qual.Length}");
            if (qual.Length > seqlen && saveLine2.Length > 0)
                _context.Warn($"(Possible short quality record before '{saveLine2.TrimEnd('\n', '\r')}')");
        }

        if (sanger && qual.Length < cntqual)
            _context.Warn($"FASTQ format '{_opt.Filename}' sequence '{t.Name}' quality skipped {cntqual - qual.Length} character(s)");
        if (ok)
            StoreClear(1);
        else
            _buff.Clear(0);
        if (seqlen > _qualsize)
            AccuracyResize(seqlen);
        _qualsize = seqlen;
        if (sanger && qual.Length > _qualsize)
            AccuracyResize(qual.Length);
        int ai = 0;
        foreach (char c in qual.ToString())
        {
            int iqual = c;
            if (c < qmin)
            {
                _context.Warn(sanger ? $"FASTQ-SANGER '{_opt.Filename}' sequence '{t.Name}' quality value '{c}' too low" : $"FASTQ-{kind} quality value too low '{_opt.Filename}' '{t.Name}' '{c}'");
                iqual = qmin;
            }

            if (c > qmax)
            {
                _context.Warn(sanger ? $"FASTQ-SANGER '{_opt.Filename}' sequence '{t.Name}' quality value '{c}' too high" : $"FASTQ-{kind} quality value too high '{_opt.Filename}' '{t.Name}' '{c}'");
                iqual = qmax;
            }

            _acc![ai++] = QualityValue(iqual, kind);
        }

        return true;
    }

    private static string Prefix(string s, int n) => n <= 0 ? "" : n >= s.Length ? s : s[..n];

    // ---------------------------------------------------------------------------------------------------------------- sam

    // ajStrTokenNextParseNoskip: the token up to the next delimiter (possibly empty); the position moves past the delimiter
    private sealed class NoskipTok(string s, string delim)
    {
        private int _pos;

        public void Next(out string token)
        {
            token = "";
            if (_pos >= s.Length)
                return;
            int e = _pos;
            while (e < s.Length && !delim.Contains(s[e], StringComparison.Ordinal))
                e++;
            token = s[_pos..e];
            _pos = e + 1;
        }
    }

    // seqReadSam
    private bool ReadSam(EmbossSeq t)
    {
        bool ok = ReadLine();
        if (_count == 1)
        {
            while (ok && _line.Length > 0 && _line[0] == '@')
            {
                var h = new StrTok(_line, "\t");
                h.Next(out string token);
                _seqToken = token;
                char kind = token.Length > 1 ? token[1] : '\0';
                bool bad = kind switch { 'H' => token != "@HD", 'S' => token != "@SQ", 'R' => token != "@RG", 'P' => token != "@PG", 'C' => token != "@CO", _ => true };
                if (bad)
                {
                    _context.Error($"bad sam format header record '{_line}'");
                    return false;
                }

                ok = ReadLine();
            }
        }

        if (!ok)
            return false;
        var count = new StrTok(_line, "\t");
        int fields = 0;
        while (count.Next(out _))
            fields++;
        if (fields < 11)
            return false;
        var tk = new NoskipTok(_line, "\t\n");
        tk.Next(out string tok);
        _seqToken = tok;
        t.Name = NameNospace(tok);
        tk.Next(out tok);
        _seqToken = tok;
        if (tok.Length > 0 && !CStr.ToUint(tok, out _))
        {
            _context.Error($"SAM {_opt.Filename} '{t.Name}' invalid FLAG value {tok}\n");
            return false;
        }

        tk.Next(out tok);
        _seqToken = tok;
        tk.Next(out tok);
        _seqToken = tok;
        if (tok.Length > 0 && !CStr.ToUint(tok, out _))
        {
            _context.Error($"SAM {_opt.Filename} '{t.Name}' invalid POS value {tok}\n");
            return false;
        }

        tk.Next(out tok);
        _seqToken = tok;
        tk.Next(out tok);
        _seqToken = tok;
        tk.Next(out tok);
        _seqToken = tok;
        tk.Next(out tok);
        _seqToken = tok;
        if (tok.Length > 0 && !CStr.ToUint(tok, out _))
        {
            _context.Error($"SAM {_opt.Filename} '{t.Name}' invalid MPOS value {tok}\n");
            return false;
        }

        tk.Next(out tok);
        _seqToken = tok;
        if (tok.Length > 0 && !CStr.ToInt(tok, out _))
        {
            _context.Error($"SAM {_opt.Filename} '{t.Name}' invalid ISIZE value {tok}\n");
            return false;
        }

        tk.Next(out tok);
        _seqToken = tok;
        Append(t, tok);
        int seqlen = tok.Length;
        tk.Next(out tok);
        _seqToken = tok;
        if (tok != "*")
        {
            if (tok.Length != seqlen)
            {
                _context.Error($"SAM quality length mismatch '{_opt.Filename}' '{t.Name}' expected: {seqlen} found: 0 '{t.SeqText}' '{tok}'");
                return false;
            }

            if (seqlen > _qualsize)
                AccuracyResize(seqlen);
            _qualsize = seqlen;
            int pos = 0;
            foreach (char c in tok)
            {
                int iqual = c;
                pos++;
                if (c < 33)
                {
                    _context.Warn($"SAM '{_opt.Filename}' sequence '{t.Name}' quality value {pos} '{c}' too low");
                    iqual = 33;
                }

                if (c > 126)
                {
                    _context.Warn($"SAM '{_opt.Filename}' sequence '{t.Name}' quality value '{c}' too high");
                    iqual = 126;
                }

                _acc![pos - 1] = QualityValue(iqual, "SANGER");
            }
        }

        _seqToken = "";
        return true;
    }

    // ---------------------------------------------------------------------------------------------------------------- gff

    private static readonly Regex RegGffTyp = new(@"^##([DR]NA|Protein) +([^ \t\r\n]+)", RegexOptions.CultureInvariant);
    private static readonly Regex RegGff3Typ = new(@"^#!Type (.*)", RegexOptions.CultureInvariant);

    // seqReadGff2
    private bool ReadGff2(EmbossSeq t)
    {
        bool ok = ReadLine();
        if (!ok)
            return false;
        if (!PrefixC(_line, "##gff-version "))
        {
            StoreReset();
            return false;
        }

        string verstr = _line;
        AssignText();
        string? outstr = null;
        string typstr = "";
        bool isseq = false;
        ok = ReadLine();
        while (ok && PrefixC(_line, "##"))
        {
            var m = RegGffTyp.Match(_line);
            if (m.Success)
            {
                isseq = true;
                typstr = m.Groups[1].Value;
                t.Name = m.Groups[2].Value;
                outstr = $"##Type {typstr} {t.Name}";
            }
            else if (PrefixC(_line, "##end-"))
                isseq = false;
            else if (isseq)
                Append(t, _line);
            ok = ReadLine();
        }

        if (t.Seq.Length == 0)
        {
            StoreReset();
            return false;
        }

        bool dofeat = false;
        List<string>? ftfile = null;
        if (ok && _opt.Features && !_opt.UfoGiven)
        {
            dofeat = true;
            ftfile = [Chomp(verstr)];
            if (outstr is not null)
                ftfile.Add(outstr);
            else
                ftfile.Add("");
            while (ok && !PrefixC(_line, "##"))
            {
                ftfile.Add(Chomp(_line));
                ok = ReadLine();
            }
        }

        if (dofeat)
        {
            var table = new EmbossSeqFeatureTable(_context) { SeqId = t.Name };
            table.Found = table.ReadGff3(new EmbossFeatureLines(ftfile!));          // the feature format of both GFF sequence readers is "gff" / "gff3": the GFF3 feature reader (the static seqFtFmtGff)
            table.FormatName = "gff";
            t.Fttable = table.Found ? table : null;
        }

        if (typstr == "Protein")
            t.SetProt();
        else if (IsNucleotideSeq(t))
            t.SetNuc();
        else
            t.SetProt();
        if (t.Fttable is { } tab)
            tab.Length = t.Seq.Length;
        _buff.Clear(0);
        return true;
    }

    // seqReadGff3
    private bool ReadGff3(EmbossSeq t)
    {
        bool ok = ReadLine();
        if (!ok)
            return false;
        _line = RemoveWhiteExcess(_line);
        if (_line != "##gff-version 3")
        {
            StoreReset();
            return false;
        }

        string verstr = _line;
        AssignText();
        string? outstr = null;
        string typstr = "";
        ok = ReadLine();
        while (ok && PrefixC(_line, "#"))
        {
            if (PrefixC(_line, "##sequence-region"))
            {
                var h = new StrTok(_line, " \t");
                h.Next(out _);
                h.Next(out string name);
                t.Name = name;
            }
            else if (PrefixC(_line, "##FASTA"))
                break;
            var m = RegGff3Typ.Match(_line);
            if (m.Success)
            {
                typstr = m.Groups[1].Value;
                outstr = $"#!Type {typstr}";
            }

            ok = ReadLine();
        }

        bool dofeat = false;
        List<string>? ftfile = null;
        if (ok && _opt.Features && !_opt.UfoGiven)
        {
            dofeat = true;
            ftfile = [verstr, outstr ?? ""];
        }

        while (ok)
        {
            if (PrefixC(_line, "##"))
            {
                if (_line.StartsWith("##FASTA", StringComparison.OrdinalIgnoreCase))
                    break;
                if (PrefixC(_line, "##gff-version "))
                    return false;
            }

            if (dofeat)
                ftfile!.Add(Chomp(_line));
            else if (t.Name.Length == 0)
            {
                if (ExtractFirst(_line, out _, out string word))
                    t.Name = word;
            }

            ok = ReadLine();
        }

        if (!_line.StartsWith("##FASTA", StringComparison.OrdinalIgnoreCase))
        {
            StoreReset();
            return false;
        }

        if (dofeat)
        {
            var table = new EmbossSeqFeatureTable(_context) { SeqId = t.Name };
            table.Found = table.ReadGff3(new EmbossFeatureLines(ftfile!));
            table.FormatName = "gff3";
            t.Fttable = table.Found ? table : null;
            if (t.Fttable is null)
                throw new InvalidOperationException("The reference program stops with a fault here: a GFF3 file read with its features and holding none leaves the feature table null, and the library asks it whether it is circular (valgrind: Invalid read of size 4 in ajFeattableIsCircular).");
            t.Fttable.Length = t.Seq.Length;
            if (t.Fttable.Circular)
                t.Circular = true;
        }

        ok = ReadLine();
        bool wantseq = false;
        while (ok && !PrefixC(_line, "##"))
        {
            while (ok && PrefixC(_line, ">"))
            {
                string header = _line[1..];
                if (wantseq)
                    wantseq = false;
                else
                {
                    _ = ExtractFirst(header, out string rest, out string token);
                    if (dofeat)
                    {
                        if (token == t.Fttable!.SeqId)
                        {
                            wantseq = true;
                            t.Name = token;
                        }
                    }
                    else if (token == t.Name)
                        wantseq = true;
                    if (wantseq)
                        t.Desc = RemoveWhiteExcess(rest);
                }

                ok = ReadLine();
            }

            if (wantseq)
                Append(t, _line);
            ok = ReadLine();
        }

        if (t.Seq.Length == 0)
        {
            StoreReset();
            return false;
        }

        if (typstr == "Protein")
            t.SetProt();
        else if (IsNucleotideSeq(t))
            t.SetNuc();
        else
            t.SetProt();
        _buff.Clear(0);
        return true;
    }

    private static readonly Regex RegPhylipTop = new(@"^ *([0-9]+) +([0-9]+)", RegexOptions.CultureInvariant);
    private static readonly Regex RegJackTop = new(@"^'(.*)'\s*$", RegexOptions.CultureInvariant);
    private static readonly Regex RegTreeconTop = new(@"^ *([0-9]+)", RegexOptions.CultureInvariant);
    private static readonly Regex RegFitchHead = new(@"^(\S+),\s+(\d+)\s+bases\n", RegexOptions.CultureInvariant);

    // The alignment-style and binary formats that are not read by this port: the entry of the format is recognised by the test the library makes on the first line(s) and the read stops there, any other text is refused by the reader
    // as it is by the library (so the formats that follow in the table are tried).
    private bool ReadUnported(EmbossSeq t, string format)
    {
        _ = t;
        if (_format != 0 && Table[_format].Name == format)
            throw new InvalidOperationException($"The sequence format '{format}' is not read by this port.");
        bool ok = ReadLine();
        string[] binary = ["nibble", "nib", "abi", "bam", "scf"];
        if (Array.IndexOf(binary, format) >= 0)
        {
            if (!ok)
                return false;
            if (format is "abi" && PrefixC(_line, "ABIF") || format is "scf" && PrefixC(_line, ".scf") || format is "bam" && _line.StartsWith('\u001f'))
                throw new InvalidOperationException($"The binary format '{format}' is not read by this port.");
            StoreReset();
            return false;
        }

        if (!ok)
            return false;
        bool match;
        switch (format)
        {
            case "clustal":
            case "aln":
                match = PrefixC(_line, "CLUSTAL");
                break;
            case "phylip":
            case "phylipnon":
                while (ok && IsWhite(_line))
                    ok = ReadLine();
                if (!ok)
                    return false;
                match = RegPhylipTop.IsMatch(_line);
                break;
            case "ace":
            case "consed":
                {
                    var h = new StrTok(_line, " \n\r");
                    h.Next(out string token, " \t");
                    match = string.Equals(token, "AS", StringComparison.OrdinalIgnoreCase);
                    break;
                }

            case "acedb":
                {
                    while (ok && (PrefixC(_line, "//") || PrefixC(_line, "\n")))
                        ok = ReadLine();
                    if (!ok)
                    {
                        StoreReset();
                        return false;
                    }

                    var h = new StrTok(_line, " \n\r");
                    h.Next(out string token, " \t");
                    match = string.Equals(token, "Peptide", StringComparison.OrdinalIgnoreCase) || string.Equals(token, "DNA", StringComparison.OrdinalIgnoreCase);
                    if (match)
                    {
                        h.Next(out token, " \t\"");
                        match = token == ":";
                        if (match)
                        {
                            h.Next(out token, "\"");
                            match = token.Length > 0;
                        }
                    }

                    break;
                }

            case "msf":
                match = RegGcgDot.IsMatch(_line) && RegGcgChk.IsMatch(_line) && RegGcgMsf.IsMatch(_line);
                if (!match)
                {
                    int n = 1;
                    while (n < 5000 && ReadLine())
                    {
                        n++;
                        if (!RegGcgDot.IsMatch(_line))
                            continue;
                        match = RegGcgChk.IsMatch(_line) && RegGcgMsf.IsMatch(_line);
                        break;
                    }
                }

                break;
            case "hennig86":
                match = PrefixC(_line, "xread");
                break;
            case "jackknifer":
                match = RegJackTop.IsMatch(_line);
                break;
            case "nexus":
            case "paup":
                match = _line.StartsWith("#NEXUS", StringComparison.OrdinalIgnoreCase);
                break;
            case "treecon":
                match = RegTreeconTop.IsMatch(_line);
                break;
            case "mega":
                match = _line.StartsWith("#MEGA", StringComparison.OrdinalIgnoreCase);
                break;
            case "stockholm":
            case "pfam":
                match = PrefixC(_line, "# STOCKHOLM 1.");
                break;
            case "fitch":
                match = RegFitchHead.IsMatch(_line);
                break;
            case "experiment":
                match = PrefixSwiss(_line) == Sw.Id;
                break;
            case "gde":
                match = _line.Length > 0 && _line[0] == '#';
                break;
            default:
                match = false;
                break;
        }

        if (match)
            throw new InvalidOperationException($"The sequence format '{format}' is not read by this port (only a file that is not of the format is refused by it).");
        StoreReset();
        return false;
    }

    // seqReadGde: "#name" starts an entry, the lines up to the next '#' line are its sequence
    private bool ReadGde(EmbossSeq t)
    {
        bool ok = ReadLine();
        if (!ok)
            return false;
        if (_line.Length == 0 || _line[0] != '#')
        {
            StoreReset();
            return false;
        }

        string token = _line[1..];                      // ajStrTokenAssign without delimiters: the token is the whole rest of the line
        _seqToken = token;
        t.Name = NameNospace(token);
        ok = ReadLine();
        while (ok && !PrefixC(_line, "#"))
        {
            AppendWarn(t, _line);
            ok = ReadLine();
        }

        if (ok)
            StoreClear(1);
        else
            _buff.Clear(0);
        return true;
    }
}
