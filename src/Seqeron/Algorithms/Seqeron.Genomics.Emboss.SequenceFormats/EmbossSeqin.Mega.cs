// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the readers keep the control flow of ajseqread.c

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

// seqReadMega of ajseqread.c: the MEGA format (title, format commands, #name sequence blocks with [comments], "identical" and gap characters, gene / domain blocks that give one set of sequences each).
public sealed partial class EmbossSeqin
{
    private static readonly Regex RegMegaCommand = new(@"([^ =!]+)=([^ ;]+)", RegexOptions.CultureInvariant);
    private static readonly Regex RegMegaFeat = new("^(.*)\"[^\"]*\"(.*)$", RegexOptions.CultureInvariant);
    private static readonly Regex RegMegaSeq = new(@"^#([^ \t\n\r]+)(.*)$", RegexOptions.CultureInvariant);

    private sealed class MegaData
    {
        public Dictionary<string, AlnItem> Table = [];
        public string?[] Names = [];
        public int Nseq, Count;
        public bool Resume;
        public int CommentDepth;
        public char Identchar = '.', Indelchar = '-', Misschar = '?', Seqtype = ' ';
    }

    private MegaData? _mega;
    private bool _multidone;

    /// <summary><c>seqin-&gt;Multidone</c>: the reader finished a set of a file that holds several (phylip, phylipnon, treecon, jackknifer, nexus, mega).</summary>
    public bool Multidone => _multidone;

    // ajCvtSposToPosStart
    private static int PosStart(int len, int imin, int ipos)
    {
        int jpos = ipos < 0 ? len + ipos : ipos;
        if (jpos >= len)
            jpos = len - 1;
        if (jpos < imin)
            jpos = imin;
        return jpos;
    }

    // ajStrCutRange
    private static string CutRange(string s, int pos1, int pos2)
    {
        int ibegin = PosStart(s.Length, 0, pos1);
        int iend = PosStart(s.Length, ibegin, pos2) + 1;
        return s.Remove(ibegin, iend - ibegin);
    }

    private static string TrimStartSet(string s, string set) => s.TrimStart(set.ToCharArray());

    private static string TrimEndSet(string s, string set) => s.TrimEnd(set.ToCharArray());

    // seqReadMega
    private bool ReadMega(EmbossSeq t)
    {
        MegaData? phydata = _mega;
        bool resume = false;
        if (phydata is not null && _multidone)
            resume = phydata.Resume;
        if (phydata is null || (_multidone && resume))
        {
            int iseq = 0;
            _multidone = false;
            bool ok = true;
            var list = new List<string>();
            string tmpdesc = "";
            if (phydata is null)
            {
                ok = ReadLine();
                if (!ok)
                    return false;
                if (!PrefixCase(_line, "#MEGA"))
                {
                    StoreReset();
                    return false;
                }

                ok = ReadLine();
                if (!ok)
                    return false;
                if (!PrefixCase(_line, "TITLE") && !PrefixCase(_line, "!TITLE"))
                {
                    StoreReset();
                    return false;
                }

                int ib = PosStart(_line.Length, 0, 6);
                int ie = PosStart(_line.Length, ib, -1);
                tmpdesc = _line.Substring(ib, ie - ib + 1);
                tmpdesc = TrimEndSet(TrimStartSet(tmpdesc, ": \t"), "; \t\n\r");
                bool istitle = false, isformat = false, iscommand = false;
                char identchar = '.', indelchar = '-', misschar = '?', seqtype = ' ';
                if (_line[0] == '!')
                {
                    istitle = true;
                    if (_line.Contains(';'))
                        istitle = false;
                }

                while (ok && !PrefixC(_line, "#"))
                {
                    if (iscommand)
                    {
                        if (_line.Contains(';'))
                            iscommand = false;
                    }
                    else if (istitle)
                    {
                        string tmpstr = TrimEndSet(TrimStartSet(_line, ": \t"), "; \t\n\r");
                        tmpdesc = new StringBuilder(tmpdesc).Append(' ').Append(tmpstr).ToString();
                        if (_line.Contains(';'))
                            istitle = false;
                    }
                    else
                    {
                        if (PrefixCase(_line, "!FORMAT"))
                            isformat = true;
                        if (isformat)
                        {
                            string tmpstr = _line;
                            for (var m = RegMegaCommand.Match(tmpstr); m.Success; m = RegMegaCommand.Match(tmpstr))
                            {
                                string formatType = m.Groups[1].Value, formatValue = m.Groups[2].Value;
                                if (PrefixCase(formatType, "indel"))
                                    indelchar = formatValue[0];
                                if (PrefixCase(formatType, "ident"))
                                    identchar = formatValue[0];
                                if (PrefixCase(formatType, "match"))
                                    identchar = formatValue[0];
                                if (PrefixCase(formatType, "miss"))
                                    misschar = formatValue[0];
                                if (PrefixCase(formatType, "DataType"))
                                    seqtype = formatValue[0];
                                tmpstr = tmpstr[(m.Index + m.Length)..];
                            }

                            if (!_line.Contains(';'))
                                isformat = false;
                        }
                        else if (_line[0] == '!')
                        {
                            string tmpstr = _line;
                            for (var m = RegMegaCommand.Match(tmpstr); m.Success; m = RegMegaCommand.Match(tmpstr))
                            {
                                string formatType = m.Groups[1].Value, formatValue = m.Groups[2].Value;
                                // the gene and domain names of the header are kept by the library for the gene features of the sequence (not used by the readers' output)
                                tmpstr = tmpstr[(m.Index + m.Length)..];
                            }

                            if (!_line.Contains(';'))
                                iscommand = true;
                        }
                    }

                    ok = ReadLine();
                }

                if (isformat || istitle || iscommand)
                {
                    StoreReset();
                    return false;
                }

                phydata = new MegaData { Identchar = identchar, Indelchar = indelchar, Misschar = misschar, Seqtype = seqtype };
                _mega = phydata;
            }

            if (resume)
            {
                ok = ReadLine();
                resume = false;
                phydata.Resume = false;
                phydata.Table = [];
            }

            bool iscmd = false;
            string? tmpname = null;
            AlnItem? phyitem = null;
            while (ok)
            {
                string line = _line;
                int ipos = line.IndexOfAny(['[', ']']);
                int istart = 0;
                char ichar = ' ';
                while (ipos != -1 || (phydata.CommentDepth > 0 && istart < line.Length))
                {
                    int ilast = ipos;
                    if (ipos > -1)
                        ichar = line[ipos];
                    if (phydata.CommentDepth == 0)
                        istart = ipos;
                    if (ichar == '[')
                        phydata.CommentDepth++;
                    else if (ichar == ']' && phydata.CommentDepth > 0)
                        phydata.CommentDepth--;
                    line = CutRange(line, istart, ilast);
                    ipos = line.IndexOfAny(['[', ']']);
                    ichar = ' ';
                }

                _line = line;
                if (line.Length == 0)
                {
                    ok = ReadLine();
                    continue;
                }

                if (PrefixC(line, "!"))
                    iscmd = true;
                if (!iscmd)
                {
                    if (PrefixC(line, "#"))
                    {
                        var m = RegMegaSeq.Match(line);
                        if (!m.Success)
                        {
                            _mega = null;
                            return false;
                        }

                        SeqNameSetName(ref tmpname, m.Groups[1].Value);
                        if (tmpname is null)
                            throw new InvalidOperationException("The reference program dereferences a null name.");
                        if (!phydata.Table.TryGetValue(tmpname, out phyitem))
                        {
                            phyitem = new AlnItem { Name = tmpname, Desc = tmpdesc };
                            phydata.Table[tmpname] = phyitem;
                            list.Add(tmpname);
                            iseq++;
                        }

                        line = m.Groups[2].Value;
                    }

                    for (var m = RegMegaFeat.Match(line); m.Success; m = RegMegaFeat.Match(line))
                        line = m.Groups[1].Value + m.Groups[2].Value;
                    if (phyitem is null)
                        throw new InvalidOperationException("The reference program dereferences a null sequence (sequence text before the first name).");
                    AppendTo(phyitem.Seq, line);
                }
                else
                {
                    string tmpstr = line;
                    for (var m = RegMegaCommand.Match(tmpstr); m.Success; m = RegMegaCommand.Match(tmpstr))
                    {
                        string formatType = m.Groups[1].Value, formatValue = m.Groups[2].Value;
                        if (string.Equals(formatType, "gene", StringComparison.OrdinalIgnoreCase))
                        {
                            if (iseq > 0)
                                resume = true;
                        }

                        if (string.Equals(formatType, "domain", StringComparison.OrdinalIgnoreCase))
                        {
                            if (iseq > 0)
                                resume = true;
                        }

                        tmpstr = tmpstr[(m.Index + m.Length)..];
                    }

                    if (line.Contains(';'))
                        iscmd = false;
                }

                if (resume)
                    break;
                ok = ReadLine();
            }

            phydata.Names = [.. list];
            phydata.Nseq = iseq;
            phydata.Count = 0;
        }

        phydata = _mega!;
        if (phydata.Nseq == 0 || phydata.Names[0] is null)
            throw new InvalidOperationException("The reference program reads the name of a sequence that is not there (a MEGA file with no sequence).");
        var firstitem = phydata.Table[phydata.Names[0]!];
        int i = phydata.Count;
        string name = phydata.Names[i]!;
        var readitem = phydata.Table[name];
        t.Name = name;
        if (i > 0)
            phydata.Names[i] = null;
        if (resume)
        {
            phydata.Resume = true;
        }

        t.Desc = readitem.Desc;
        t.Seq.Clear().Append(readitem.Seq.ToString());
        if ("nNrRdD".Contains(phydata.Seqtype))
            t.SetNuc();
        else if ("pP".Contains(phydata.Seqtype))
            t.SetProt();
        string first = firstitem.Seq.ToString();
        for (int k = 0; k < t.Seq.Length; k++)
        {
            char c = t.Seq[k];
            if (c == phydata.Indelchar)
                t.Seq[k] = '-';
            else if (c == phydata.Identchar)
            {
                // the characters of the first sequence; the library reads the zero filled memory behind its end (the string grows by ajStrSetResRound, which zero fills) and the memory that follows it, which holds zeros in every case measured
                t.Seq[k] = k < first.Length ? first[k] : '\0';
            }
        }

        phydata.Count++;
        if (phydata.Count >= phydata.Nseq)
        {
            _multidone = true;
            phydata.Names[0] = null;
            if (!phydata.Resume)
            {
                _buff.Clear(0);
                _mega = null;
            }
        }

        return true;
    }
}
