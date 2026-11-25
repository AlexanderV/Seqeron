// SPDX-License-Identifier: LGPL-2.1-or-later
// Port of the EMBOSS 6.6.0 AJAX / NUCLEUS libraries (library source files carry the GNU Lesser General Public License, version 2.1 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU Lesser General Public License, version 2.1 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the readers keep the control flow of ajseqread.c

using System.Text;
using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Emboss.SequenceFormats;

using Seqeron.Genomics.Infrastructure;

// The alignment formats of phylogeny programs (hennig86, jackknifer, nexus, mega) of ajseqread.c, on the seqin->SeqData machinery.
public sealed partial class EmbossSeqin
{
    private static readonly Regex RegJackSeq = new(@"^[(]([^)]+)(.*)$", RegexOptions.CultureInvariant);

    // seqnameSetName: the name of the last call is kept when the text gives none (a word of colons)
    private static void SeqNameSetName(ref string? name, string str)
    {
        if (IsWord(str))
        {
            foreach (string part in str.Split(':'))
            {
                if (part.Length > 0)
                    name = part;
            }

            return;
        }

        name = RemoveWhiteExcess(str).Replace(' ', '_');
    }

    // the item the readers of the sets return (Names[Count] of the table; the name is deleted for the readers that do so)
    private bool ReturnSetItem(EmbossSeq t, bool deleteName)
    {
        var data = _seqData!;
        int i = data.Count;
        string? name = data.NameArr is { } arr && i < arr.Length ? arr[i] : null;
        if (name is null || !data.Table.TryGetValue(name, out var item))
            throw new InvalidOperationException("The reference program dereferences a null sequence (the name of the alignment item is not set).");
        t.Name = name;
        if (deleteName)
            data.NameArr![i] = null;
        t.Seq.Clear().Append(item.Seq.ToString());
        data.Count++;
        if (data.Count >= data.Nseq)
        {
            _multidone = true;
            _buff.Clear(0);
            _seqData = null;
        }

        return true;
    }

    // seqReadJackknifer
    private bool ReadJackknifer(EmbossSeq t)
    {
        if (_seqData is null)
        {
            _multidone = false;
            bool ok = ReadLine();
            if (!ok)
                return false;
            if (!RegJackTop.IsMatch(_line))
            {
                StoreReset();
                return false;
            }

            ok = ReadLine();
            var data = new AlnData();
            _seqData = data;
            var list = new List<string>();
            string? tmpname = null;
            AlnItem? item = null;
            while (ok)
            {
                if (_line.Length == 0)
                {
                    ok = ReadLine();
                    continue;
                }

                if (PrefixC(_line, ";"))
                    break;
                string line = _line;
                if (PrefixC(line, "("))
                {
                    var m = RegJackSeq.Match(line);
                    if (!m.Success)
                    {
                        _seqData = null;
                        return false;
                    }

                    SeqNameSetName(ref tmpname, m.Groups[1].Value);
                    if (tmpname is null)
                        throw new InvalidOperationException("The reference program dereferences a null name.");
                    if (!data.Table.TryGetValue(tmpname, out item))
                    {
                        item = new AlnItem { Name = tmpname };
                        data.Table[tmpname] = item;
                        list.Add(tmpname);
                    }

                    line = m.Groups[2].Value;
                }

                if (item is null)
                    throw new InvalidOperationException("The reference program dereferences a null sequence (sequence text before the first name).");
                AppendTo(item.Seq, line);
                ok = ReadLine();
            }

            data.NameArr = [.. list];
            data.Nseq = list.Count;
            data.Count = 0;
        }

        return ReturnSetItem(t, true);
    }

    private static readonly Regex RegHennigHead = new(@"[^1-4? \t]", RegexOptions.CultureInvariant);
    private static readonly Regex RegHennigSeq = new(@"^([^ \t\n\r]+)", RegexOptions.CultureInvariant);

    // seqHennig86Readseq
    private static void Hennig86Readseq(string rdline, AlnData data)
    {
        var m = RegHennigSeq.Match(rdline);
        if (m.Success && data.Table.TryGetValue(m.Groups[1].Value, out var item))
            AppendTo(item.Seq, rdline[(m.Index + m.Length)..]);
    }

    // seqReadHennig86: the library converts the characters of a sequence line with a loop that has no end ("for(cp = ...; cp; cp++)" runs over the heap, writing, until it dies: valgrind Invalid write) and does not read the line
    // after it, so a line that starts with a token after a name is a fault of the reference (the run is refused); the rest of the reader is the code of the library
    private bool ReadHennig86(EmbossSeq t)
    {
        if (_seqData is null)
        {
            bool ok = ReadLine();
            if (!ok)
                return false;
            if (!PrefixC(_line, "xread"))
            {
                StoreReset();
                return false;
            }

            for (int i = 0; i < 2; i++)
            {
                ok = ReadLine();
                if (!ok)
                {
                    StoreReset();
                    return false;
                }
            }

            var top = RegPhylipTop.Match(_line);
            if (!top.Success)
                return false;
            long iseq = ToUInt(top.Groups[1].Value);
            var data = new AlnData();
            _seqData = data;
            var list = new List<string>();
            long jseq = 0;
            ok = ReadLine();
            while (ok && jseq < iseq)
            {
                if (!RegHennigHead.IsMatch(_line))
                    return false;
                var item = new AlnItem { Name = _line };
                ok = ReadLine();
                if (ok && RegHennigSeq.IsMatch(_line))
                    throw new InvalidOperationException("The reference program does not stop on a sequence line (seqReadHennig86 loops over memory without an end).");
                data.Table[item.Name] = item;
                list.Add(item.Name);
                jseq++;
                ok = ReadLine();
            }

            data.NameArr = new string?[Math.Min(iseq, 1 << 24)];
            for (int i = 0; i < data.NameArr.Length && i < list.Count; i++)
                data.NameArr[i] = list[i];
            while (ReadLine())
                Hennig86Readseq(_line, data);
            data.Nseq = (int)Math.Min(iseq, int.MaxValue);
            data.Count = 0;
        }

        var d = _seqData!;
        if (d.Count >= d.Nseq)
        {
            _buff.Clear(0);
            _seqData = null;
            return false;
        }

        int n = d.Count;
        string? name = n < d.NameArr!.Length ? d.NameArr[n] : null;
        if (name is null || !d.Table.TryGetValue(name, out var found))
            throw new InvalidOperationException("The reference program dereferences a null sequence (the name of the alignment item is not set).");
        t.Name = name;
        t.Seq.Clear().Append(found.Seq.ToString());
        d.Count++;
        return true;
    }
}
