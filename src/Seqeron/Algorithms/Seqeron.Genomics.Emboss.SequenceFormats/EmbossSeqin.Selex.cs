// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the readers keep the control flow of ajseqread.c

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

// seqReadSelex (with seqSelexHeader, seqSelexPos, seqSelexReadBlock, seqSelexAppend) and seqReadStockholm of ajseqread.c: the sequence sets of HMMER / Pfam alignments.
public sealed partial class EmbossSeqin
{
    private sealed class SelexData
    {
        public int N, Count;
        public string[] Name = [];
        public StringBuilder[] Str = [];

        /// <summary><c>sq[n]-&gt;wt</c>: the weight of the <c>#=SQ</c> line (0 when there is none).</summary>
        public float[] Wt = [];
    }

    private SelexData? _selex;
    private static readonly Regex RegStockholmSeq = new(@"^([^ \t\n]+)[ \t]+([^ \t\n]+)[ \t]+", RegexOptions.CultureInvariant);

    private static string[] SplitWhite(string s) => s.Split([' ', '\t', '\n', '\r', '\f', '\v'], StringSplitOptions.RemoveEmptyEntries);

    // ajStrTrimWhiteEnd
    private static string TrimWhiteEnd(string s) => s.TrimEnd(' ', '\t', '\n', '\r', '\f', '\v');

    // seqSelexHeader: false for a line that is no header line
    private static bool SelexHeader(SelexData selex, string line, ref bool named, ref uint sqcnt)
    {
        if (PrefixC(line, "#=ID") || PrefixC(line, "#=AC") || PrefixC(line, "#=DE") || PrefixC(line, "#=AU") || PrefixC(line, "#=GA") || PrefixC(line, "#=TC") || PrefixC(line, "#=NC"))
            return true;
        if (PrefixC(line, "#=SQ"))
        {
            if (sqcnt >= selex.N)
                throw new InvalidOperationException("The reference program writes past the end of the sequence array (more #=SQ lines than sequences).");
            var tokens = line.Split([' ', '\t', '\n'], StringSplitOptions.RemoveEmptyEntries);
            if (tokens.Length > 1)
                selex.Name[sqcnt] = tokens[1];
            if (tokens.Length > 2)
                CStr.ToFloat(tokens[2], out selex.Wt[sqcnt]);
            named = true;
            ++sqcnt;
            return true;
        }

        return false;
    }

    // seqSelexPos
    private static void SelexPos(string line, ref uint begin, ref uint end)
    {
        if (line.Length == 0)
            throw new InvalidOperationException("The reference program reads before a line it has not got (a Selex block that ends the file).");
        uint len = (uint)line.Length - 1;
        if (len == 0)
        {
            begin = 0;
            end = 0;
            return;
        }

        uint pos = len - 1;
        end = pos > end ? pos : end;
        int p = 0;
        while (p < line.Length && line[p] != ' ')
            ++p;
        while (p < line.Length && line[p] == ' ')
            ++p;
        pos = (uint)p;
        begin = pos < begin ? pos : begin;
    }

    // seqSelexAppend
    private static void SelexAppend(string src, StringBuilder dest, uint beg, uint end)
    {
        uint len = unchecked(end - beg + 1);
        if (beg >= src.Length)
        {
            for (uint i = 0; i < len; ++i)
                dest.Append('-');
            return;
        }

        uint pad = unchecked(end - (uint)src.Length + 2);
        for (int p = (int)beg; p < src.Length && src[p] != '\n'; ++p)
        {
            char c = src[p];
            if (c is '.' or '_' or ' ')
                c = '-';
            dest.Append(c);
        }

        if (pad > 1 << 24)
            throw new InvalidOperationException("The reference program pads a line with an unsigned count that wrapped.");
        for (uint i = 0; i < pad; ++i)
            dest.Append('-');
    }

    // seqSelexReadBlock
    private void SelexReadBlock(SelexData selex, ref bool named, uint n, ref string line)
    {
        uint begin = uint.MaxValue >> 1, end = 0;
        var seqs = new string[n];
        bool ok = true;
        uint cnt = 0;
        while (PrefixC(line, "\n"))
        {
            ok = ReadLine();
            line = _line;
        }

        while (ok)
        {
            SelexPos(line, ref begin, ref end);
            if (PrefixC(line, "#=SS") && cnt == 0)
                throw new InvalidOperationException("The reference program writes before the start of its array (a #=SS line before the first sequence).");
            if (!PrefixC(line, "#"))
            {
                if (cnt >= n)
                    throw new InvalidOperationException("The reference program writes past the end of the sequence array (more lines in a block than sequences).");
                var tok = SplitWhite(line);
                string first = tok.Length > 0 ? tok[0] : "";
                if (!named)
                    selex.Name[cnt] = first;
                else if (first.Length == 0 || !selex.Name[cnt].StartsWith(first, StringComparison.Ordinal))
                    _context.Warn($"Selex format sequence names do not match ['{selex.Name[cnt]}' '{first}']");
                seqs[cnt] = line;
                ++cnt;
            }

            ok = ReadLine();
            line = _line;
            if (PrefixC(line, "\n"))
                ok = false;
        }

        if (cnt != n)
            _context.Warn($"Selex format expected {n} sequences in block, found {cnt}");
        if (cnt > n)
            cnt = n;
        for (int i = 0; i < cnt; ++i)
            SelexAppend(seqs[i], selex.Str[i], begin, end);
        named = true;
    }

    // seqReadSelex
    private bool ReadSelex(EmbossSeq t)
    {
        SelexData? selex = _selex;
        string line = "";
        if (selex is null)
        {
            _buff.SetBuffered();
            bool ok = false, isseq = false, first = true;
            while (!isseq)
            {
                string? raw = _buff.Next();
                ok = raw is not null;
                if (raw is null)
                    break;
                line = raw;
                if (first)
                {
                    first = false;
                    if (!PrefixC(line, "#"))
                    {
                        StoreReset();
                        return false;
                    }
                }

                line = RemoveWhiteExcess(line);
                if (line.Length == 0 || line[0] == '#')
                    continue;
                isseq = true;
            }

            if (!ok && !isseq)
                return false;
            uint n = 1;
            ok = true;
            while (ok && _buff.Next() is { } l1)
            {
                line = RemoveWhiteExcess(l1);
                if (line.Length > 0 && line[0] == '#')
                    continue;
                if (line.Length == 0)
                    ok = false;
                else
                    ++n;
            }

            _buff.Rewind();
            StoreReset();
            selex = new SelexData { N = (int)n, Name = new string[n], Str = new StringBuilder[n], Wt = new float[n] };
            for (int i = 0; i < n; i++)
            {
                selex.Name[i] = "";
                selex.Str[i] = new StringBuilder();
            }

            bool named = false, head = true;
            uint sqcnt = 0;
            line = "";
            while (head && ReadLine())
            {
                line = _line;
                if (PrefixC(line, "#=RF") || PrefixC(line, "#=CS"))
                    break;
                if (PrefixC(line, "#="))
                {
                    head = SelexHeader(selex, line, ref named, ref sqcnt);
                    continue;
                }

                if (line.Length > 0 && line[0] >= '0')
                    head = false;
            }

            ok = true;
            while (ok && !PrefixC(line, "# ID"))
            {
                SelexReadBlock(selex, ref named, n, ref line);
                ok = ReadLine();
                line = _line;
            }

            if (ok)
                StoreClear(1);
            else
                _buff.Clear(0);
            _selex = selex;
        }

        if (selex.Count >= selex.N)
        {
            _selex = null;
            return false;
        }

        int idx = selex.Count;
        t.Seq.Clear().Append(selex.Str[idx].ToString());
        t.Name = selex.Name[idx];
        t.Weight = selex.Wt[idx];
        ++selex.Count;
        _buff.Clear(0);
        return true;
    }

    // ajFmtScanS("%*s%*s%d"): the integer after two words (an unsigned count; a text that does not scan leaves it alone)
    private static void ScanCountAfterTwoWords(string line, ref uint n)
    {
        var tokens = SplitWhite(line);
        if (tokens.Length < 3)
            return;
        var m = Regex.Match(tokens[2], @"^[+-]?[0-9]+", RegexOptions.CultureInvariant);
        if (m.Success && long.TryParse(m.Value, out long v))
            n = unchecked((uint)v);
    }

    private StockData? _stock;

    private sealed class StockData
    {
        public uint N, Count;
        public string[] Name = [];
        public StringBuilder[] Str = [];
    }

    // seqReadStockholm
    private bool ReadStockholm(EmbossSeq t)
    {
        StockData? stock = _stock;
        if (stock is null)
        {
            _buff.SetBuffered();
            bool ok = ReadLine();
            string line = TrimWhiteEnd(_line);
            if (!ok || !PrefixC(line, "# STOCKHOLM 1."))
            {
                StoreReset();
                return false;
            }

            uint n = 0;
            while (ok && (PrefixC(line, "#") || line.Length == 0))
            {
                if (PrefixC(line, "#=GF SQ"))
                    ScanCountAfterTwoWords(line, ref n);
                ok = ReadLine();
                line = TrimWhiteEnd(_line);
            }

            if (n == 0)
            {
                while (ok && line.Length > 0)
                {
                    if (!PrefixC(line, "#") && line != "//")
                        n++;
                    ok = ReadLine();
                    line = TrimWhiteEnd(_line);
                }
            }

            StoreReset();
            ok = ReadLine();
            line = TrimWhiteEnd(_line);
            if (n > 1 << 24)
                throw new InvalidOperationException("The reference program allocates an array of the size of an unsigned count that is out of range.");
            stock = new StockData { N = n, Name = new string[n], Str = new StringBuilder[n] };
            for (int i = 0; i < n; i++)
            {
                stock.Name[i] = "";
                stock.Str[i] = new StringBuilder();
            }

            uint scnt = int.MaxValue;
            while (ok && !PrefixC(line, "//"))
            {
                var m = RegStockholmSeq.Match(line);
                if (m.Success)
                {
                    // the #=GF / #=GS / #=GC lines are kept by the library for the data of the alignment (not used by the sequence)
                }
                else if (line.Length > 0)
                {
                    if (SplitWhite(line).Length > 1)
                    {
                        ++scnt;
                        if (scnt >= n)
                            scnt = 0;
                        if (n == 0)
                            throw new InvalidOperationException("The reference program writes into an array of no elements (a Stockholm alignment with no sequence count).");
                        var tok = SplitWhite(line);
                        string namstr = tok[0], seqstr = tok[1];
                        if (stock.Name[scnt].Length == 0)
                            stock.Name[scnt] = namstr;
                        else if (namstr != stock.Name[scnt])
                            _context.Warn($"Bad stockholm format found id {scnt} '{namstr}' expect '{stock.Name[scnt]}'");
                        stock.Str[scnt].Append(seqstr);
                    }
                    else
                    {
                        if (scnt >= n)
                            throw new InvalidOperationException("The reference program writes outside its array (a Stockholm line with one word before the first sequence).");
                        stock.Str[scnt].Append(line);
                    }
                }

                ok = ReadLine();
                line = TrimWhiteEnd(_line);
            }

            while (ok && !PrefixC(_line, "# STOCKHOLM 1."))
                ok = ReadLine();
            if (ok)
                StoreClear(1);
            else
                _buff.Clear(0);
            _stock = stock;
            _buff.Clear(0);
        }

        if (stock.Count >= stock.N)
        {
            _stock = null;
            return false;
        }

        int idx = (int)stock.Count;
        t.Seq.Clear().Append(stock.Str[idx].ToString());
        t.Name = stock.Name[idx];
        ++stock.Count;
        return true;
    }
}
