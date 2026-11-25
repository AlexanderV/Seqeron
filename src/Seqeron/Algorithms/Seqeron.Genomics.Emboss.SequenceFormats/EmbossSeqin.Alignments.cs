// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589 // the readers keep the control flow of ajseqread.c

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

// The alignment-style formats of ajseqread.c that share the seqin->SeqData machinery (SeqPMsfData): the first call reads the whole alignment into a table of items, every call returns the next item as a sequence, the call after the last one
// ends the alignment: clustal (aln), msf.
public sealed partial class EmbossSeqin
{
    public sealed class AlnItem
    {
        public string Name = "";
        public readonly StringBuilder Seq = new();
        public float Weight = 1.0f;
        public uint Check;
        public uint Len;
        public string Desc = "";
    }

    public sealed class AlnData
    {
        public readonly Dictionary<string, AlnItem> Table = [];
        public readonly List<string> Names = [];
        public int Count;

        /// <summary>The names as the PHYLIP readers keep them (an array of the size of the header, entries not set are null) and the number of sequences of the header.</summary>
        public string?[]? NameArr;
        public int Nseq;
    }

    private AlnData? _seqData;

    // the call that returns an item of the alignment read by an earlier call (the end of the data deletes it and ends the alignment)
    private bool NextAlignmentItem(EmbossSeq t)
    {
        var data = _seqData!;
        if (data.Count >= data.Names.Count)
        {
            _buff.Clear(0);
            _seqData = null;
            return false;
        }

        string name = data.Names[data.Count];
        var item = data.Table[name];
        t.Name = name;
        t.Weight = item.Weight;
        t.Seq.Clear().Append(item.Seq.ToString());
        data.Count++;
        return true;
    }

    private static void AppendTo(StringBuilder sb, string line)
    {
        foreach (char c in line)
        {
            if (char.IsAsciiLetter(c) || c is '*' or '.' or '~' or '?' or '#' or '+' or '-')
                sb.Append(c);
        }
    }

    // seqReadClustal
    private bool ReadClustal(EmbossSeq t)
    {
        if (_seqData is null)
        {
            bool ok = ReadLine();
            if (!ok)
                return false;
            if (!PrefixC(_line, "CLUSTAL"))
            {
                StoreReset();
                return false;
            }

            while (ok)
            {
                ok = ReadLine();
                if (ok && !IsWhite(_line))
                    break;
            }

            if (!ok)
            {
                StoreReset();
                return false;
            }

            var data = new AlnData();
            _seqData = data;
            while (ok && ExtractFirst(_line, out string seqstr, out string name))
            {
                var item = new AlnItem { Name = name };
                AppendTo(item.Seq, seqstr);
                data.Table[name] = item;                       // ajTablePut: an item of the same name is replaced
                data.Names.Add(item.Name);
                ok = ReadLine();
            }

            while (ReadLine())
            {
                if (ExtractFirst(_line, out string rest, out string word) && data.Table.TryGetValue(word, out var found))
                    AppendTo(found.Seq, rest);
            }
        }

        return NextAlignmentItem(t);
    }

    private static readonly Regex RegGcgMsfnam = new(@"[Nn][Aa][Mm][Ee]:[ \t]*([^ \t]+)", RegexOptions.CultureInvariant);
    private static readonly Regex RegGcgMsflen = new(@"[Ll][Ee][Nn]:[ \t]*([0-9]+)", RegexOptions.CultureInvariant);
    private static readonly Regex RegGcgWgt = new(@"[Ww][Ee][Ii][Gg][Hh][Tt]:[ \t]*([0-9.]+)", RegexOptions.CultureInvariant);

    // seqGcgMsfDots
    private bool GcgMsfDots(int maxlines)
    {
        int nlines = 0;
        while (nlines < maxlines)
        {
            if (nlines++ > 0)
            {
                if (!ReadLine())
                    return false;
            }

            if (nlines > maxlines)
                return false;
            if (!RegGcgDot.IsMatch(_line))
                continue;
            if (!RegGcgChk.IsMatch(_line))
                return false;
            return RegGcgMsf.IsMatch(_line);
        }

        return false;
    }

    // seqGcgMsfHeader
    private static AlnItem? GcgMsfHeader(string line)
    {
        var nam = RegGcgMsfnam.Match(line);
        if (!nam.Success)
            return null;
        var chk = RegGcgChk.Match(line);
        if (!chk.Success)
            return null;
        var item = new AlnItem { Name = nam.Groups[1].Value };
        CStr.ToUint(chk.Groups[1].Value, out item.Check);
        var len = RegGcgMsflen.Match(line);
        if (len.Success)
            CStr.ToUint(len.Groups[1].Value, out item.Len);
        var wgt = RegGcgWgt.Match(line);
        if (wgt.Success)
            CStr.ToFloat(wgt.Groups[1].Value, out item.Weight);
        else
            item.Weight = 1.0f;
        return item;
    }

    // seqReadMsf
    private bool ReadMsf(EmbossSeq t)
    {
        if (_seqData is null)
        {
            bool ok = ReadLine();
            if (!ok)
                return false;
            if (PrefixC(_line, "!!"))
            {
                if (PrefixC(_line, "!!AA_MULTIPLE_ALIGNMENT"))
                    t.SetProt();
                if (PrefixC(_line, "!!NA_MULTIPLE_ALIGNMENT"))
                    t.SetNuc();
            }

            if (!GcgMsfDots(5000))
            {
                StoreReset();
                return false;
            }

            _buff.SetUnbuffered();
            var data = new AlnData();
            _seqData = data;
            ok = ReadLine();
            while (ok && !PrefixC(_line, "//"))
            {
                ok = ReadLine();
                var item = GcgMsfHeader(_line);
                if (item is not null)
                {
                    data.Table[item.Name] = item;
                    data.Names.Add(item.Name);
                }
            }

            while (ReadLine())
            {
                if (ExtractFirstWord(_line, out string rest, out string word) && data.Table.TryGetValue(word, out var found))
                    AppendTo(found.Seq, rest);
            }
        }

        if (NextAlignmentItem(t))
        {
            return true;
        }

        return false;
    }

    private static readonly Regex RegFitchHeadSeq = new(@"^(\S+),\s+(\d+)\s+bases\n", RegexOptions.CultureInvariant);
    private static readonly Regex RegMaseHead = new(@"^(;+)", RegexOptions.CultureInvariant);

    // seqReadFitch
    private bool ReadFitch(EmbossSeq t)
    {
        bool ok = ReadLine();
        var m = RegFitchHeadSeq.Match(_line);
        if (!m.Success)
        {
            StoreReset();
            return false;
        }

        SetName(t, m.Groups[1].Value);
        CStr.ToUint(m.Groups[2].Value, out uint ilen);
        _buff.SetUnbuffered();
        ok = ReadLine();
        while (ok && (uint)t.Seq.Length < ilen)
        {
            Append(t, _line);
            ok = ReadLine();
        }

        if (ok)
            _buff.Clear(1);
        else
            _buff.Clear(0);
        return true;
    }

    // seqReadMase
    private bool ReadMase(EmbossSeq t)
    {
        bool ok = ReadLine();
        if (!RegMaseHead.IsMatch(_line))
        {
            StoreReset();
            return false;
        }

        _buff.SetUnbuffered();
        StringBuilder? des = null;
        while (ok)
        {
            var m = RegMaseHead.Match(_line);
            if (!m.Success)
                break;
            if (m.Groups[1].Length == 1)
            {
                string post = _line[m.Length..];
                if (des is not null)
                    des.Append(' ');
                (des ??= new StringBuilder()).Append(post);
            }

            ok = ReadLine();
        }

        _line = RemoveWhiteExcess(_line);
        SetName(t, _line);
        t.Desc = des is null ? "" : RemoveWhiteExcess(des.ToString());
        ok = ReadLine();
        while (ok && !RegMaseHead.IsMatch(_line))
        {
            Append(t, _line);
            ok = ReadLine();
        }

        if (ok)
            StoreClear(1);
        else
            _buff.Clear(0);
        return true;
    }

    // seqReadBiomart
    private bool ReadBiomart(EmbossSeq t)
    {
        bool ok = ReadLine();
        if (!ok)
            return false;
        int ifields = _line.Count(c => c == '\t') + 1;
        if (ifields < 2)
            return false;
        var tk = new NoskipTok(_line, "\t\n");
        tk.Next(out string tok);
        _seqToken = tok;
        Append(t, tok);
        tk.Next(out tok);
        _seqToken = tok;
        SetName(t, tok);
        var desc = new StringBuilder();
        for (int i = 2; i < ifields; i++)
        {
            tk.Next(out tok);
            _seqToken = tok;
            if (tok.Length > 0)
            {
                if (i > 2)
                    desc.Append(' ');
                desc.Append(tok);
            }
        }

        t.Desc = desc.ToString();
        _buff.Clear(0);
        _seqToken = "";
        return true;
    }

    private static readonly Regex RegTreeconTopSeq = new(@"^ *([0-9]+)", RegexOptions.CultureInvariant);

    // seqitemSetName (the name of an alignment item: a word is split at ':' and the last non-empty part kept, any other text has its blanks turned to '_')
    private static string ItemName(string str)
    {
        if (IsWord(str))
        {
            string last = "";
            foreach (string part in str.Split(':'))
            {
                if (part.Length > 0)
                    last = part;
            }

            return last;
        }

        return RemoveWhiteExcess(str).Replace(' ', '_');
    }

    // seqReadTreecon
    private bool ReadTreecon(EmbossSeq t)
    {
        if (_seqData is null)
        {
            _multidone = false;
            bool ok = ReadLine();
            if (!ok)
                return false;
            var top = RegTreeconTopSeq.Match(_line);
            if (!top.Success)
            {
                StoreReset();
                return false;
            }

            CStr.ToInt(top.Groups[1].Value, out int len);
            var data = new AlnData();
            _seqData = data;
            AlnItem? item = null;
            ok = ReadLine();
            int ilen = -1;                                   // ilen = UINT_MAX of the library, held in a signed integer
            if (!ok)
            {
                StoreReset();
                _seqData = null;
                return false;
            }

            while (ok)
            {
                if (ilen < 0)
                {
                    _line = RemoveWhiteExcess(_line);
                    if (_line.Length == 0)
                    {
                        ok = ReadLine();
                        continue;
                    }

                    item = new AlnItem { Name = ItemName(_line) };
                    data.Table[item.Name] = item;
                    data.Names.Add(item.Name);
                    ilen = 0;
                }
                else
                {
                    _line = new string(_line.Where(c => !EmbossSeqFeatureTable.IsCSpace(c)).ToArray());
                    ilen += _line.Length;
                    AppendTo(item!.Seq, _line);
                    if (ilen > len)
                    {
                        StoreReset();
                        _seqData = null;
                        return false;
                    }

                    if (ilen == len)
                        ilen = -1;
                }

                ok = ReadLine();
            }

            if (ilen >= 0)
            {
                _seqData = null;
                return false;
            }

            if (data.Names.Count == 0)
                throw new InvalidOperationException("The reference program reads the name of a sequence that does not exist (a treecon file with the sequence count and no sequence): it reads past the end of an empty array and dies by a segmentation fault.");
        }

        if (NextAlignmentItem(t))
        {
            if (_seqData!.Count >= _seqData.Names.Count)
            {
                _multidone = true;
                _buff.Clear(0);
                _seqData = null;
            }

            return true;
        }

        return false;
    }

    // seqReadAcedb
    private bool ReadAcedb(EmbossSeq t)
    {
        bool ok;
        do
        {
            ok = ReadLine();
        }
        while (ok && (PrefixC(_line, "//") || PrefixC(_line, "\n")));
        if (!ok)
        {
            StoreReset();
            return false;
        }

        var h = new StrTok(_line, " \n\r");
        h.Next(out string token, " \t");
        _seqToken = token;
        if (string.Equals(token, "Peptide", StringComparison.OrdinalIgnoreCase))
            t.SetProt();
        else if (string.Equals(token, "DNA", StringComparison.OrdinalIgnoreCase))
            t.SetNuc();
        else
        {
            StoreReset();
            _seqToken = "";
            return false;
        }

        h.Next(out token, " \t\"");
        _seqToken = token;
        if (token != ":")
        {
            StoreReset();
            _seqToken = "";
            return false;
        }

        h.Next(out token, "\"");
        _seqToken = token;
        if (token.Length == 0)
        {
            StoreReset();
            _seqToken = "";
            return false;
        }

        _buff.SetUnbuffered();
        SetName(t, token);
        ok = ReadLine();
        while (ok && !PrefixC(_line, "\n"))
        {
            Append(t, _line);
            ok = ReadLine();
        }

        _buff.Clear(0);
        _seqToken = "";
        return true;
    }
}
