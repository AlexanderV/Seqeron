#pragma warning disable S1121, S3358, S2325, S1066, S3776, S107, S1854, S1481, S2589, S3626 // the readers keep the control flow of ajseqread.c

using System.Text.RegularExpressions;

namespace Seqeron.Genomics.Infrastructure;

// seqReadPhylip (interleaved), seqReadPhylipnon (sequential) and seqPhylipReadseq of ajseqread.c: the PHYLIP formats share the alignment table of the msf reader (seqin->SeqData); the interleaved reader falls back to the sequential one
// whenever the layout does not fit (a repeated name, different lengths in the first block, an incomplete last block).
internal sealed partial class EmbossSeqin
{
    private static readonly Regex RegPhylipHead = new(@"^(.{10}) ?", RegexOptions.CultureInvariant);
    private static readonly Regex RegPhylipSeq2 = new(@"[^ \t\n\r]", RegexOptions.CultureInvariant);

    // the number a regular expression group gives ajStrToUint (an unsigned value: a text that does not scan is 0)
    private static long ToUInt(string s) => CStr.ToUint(s, out uint v) ? v : 0;

    // seqPhylipReadseq: the letters of the line are appended to the named sequence; false for a blank line, a name that is not in the table and a sequence that is longer than the header says
    private bool PhylipReadseq(string rdline, AlnData data, string? token, long len, ref long ilen, out bool done)
    {
        done = false;
        if (!RegPhylipSeq2.IsMatch(rdline))
            return false;
        if (token is null || !data.Table.TryGetValue(token, out var item))
            return false;
        AppendTo(item.Seq, rdline);
        ilen = item.Seq.Length;
        if (ilen == len)
            done = true;
        else if (ilen > len)
            return false;
        return true;
    }

    // the item of an alignment read by an earlier call (the tail of seqReadPhylip / seqReadPhylipnon)
    private bool PhylipReturn(EmbossSeq t, bool deleteName)
    {
        var data = _seqData!;
        int i = data.Count;
        string? name = i < data.NameArr!.Length ? data.NameArr[i] : null;
        if (name is null || !data.Table.TryGetValue(name, out var item))
            throw new InvalidOperationException("The reference program dereferences a null sequence (the name of the alignment item is not set).");
        t.Name = name;
        if (deleteName)
            data.NameArr[i] = null;
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

    // seqReadPhylipnon
    private bool ReadPhylipnon(EmbossSeq t)
    {
        if (_seqData is null)
        {
            _multidone = false;
            bool ok = ReadLine();
            if (!ok)
                return false;
            var top = RegPhylipTop.Match(_line);
            if (!top.Success)
            {
                StoreReset();
                return false;
            }

            long iseq = ToUInt(top.Groups[1].Value), len = ToUInt(top.Groups[2].Value);
            var data = new AlnData { NameArr = new string?[Math.Min(iseq, 1 << 24)] };
            _seqData = data;
            ok = ReadLine();
            long ilen = 0, jseq = 0;
            bool done = false;
            AlnItem? phyitem = null;
            while (ok && jseq < iseq)
            {
                var head = RegPhylipHead.Match(_line);
                if (!head.Success)
                {
                    StoreReset();
                    _seqData = null;
                    return false;
                }

                string tmp = head.Groups[1].Value;
                if (!IsWhite(tmp))
                {
                    if (jseq > 0 && ilen != len)
                    {
                        _seqData = null;
                        return false;
                    }

                    phyitem = new AlnItem { Name = ItemName(tmp) };
                    if (jseq < data.NameArr.Length)
                        data.NameArr[jseq] = phyitem.Name;
                    AppendTo(phyitem.Seq, _line[(head.Index + head.Length)..]);
                    ilen = phyitem.Seq.Length;
                    if (ilen == len)
                        done = true;
                    else if (ilen > len)
                    {
                        _seqData = null;
                        return false;
                    }

                    data.Table[phyitem.Name] = phyitem;
                }
                else
                {
                    if (phyitem is null)
                        throw new InvalidOperationException("The reference program dereferences a null sequence (a continuation line before the first name).");
                    if (PhylipReadseq(_line, data, phyitem.Name, len, ref ilen, out done) && done)
                        jseq++;
                }

                if (jseq < iseq)
                    ok = ReadLine();
            }

            if (ilen != len)
            {
                _seqData = null;
                return false;
            }

            data.Nseq = (int)Math.Min(iseq, int.MaxValue);
            data.Count = 0;
        }

        return PhylipReturn(t, true);
    }

    // seqReadPhylip
    private bool ReadPhylip(EmbossSeq t)
    {
        _buff.SetBuffered();
        if (_seqData is null)
        {
            _multidone = false;
            bool ok = ReadLine();
            while (ok && IsWhite(_line))
                ok = ReadLine();
            if (!ok)
                return false;
            var top = RegPhylipTop.Match(_line);
            if (!top.Success)
            {
                StoreReset();
                return false;
            }

            long iseq = ToUInt(top.Groups[1].Value), len = ToUInt(top.Groups[2].Value);
            var data = new AlnData();
            _seqData = data;
            var list = new List<string>();
            ok = ReadLine();
            long ilen = 0, jseq = 0, maxlen = 0;
            bool done = false;
            while (ok && jseq < iseq)
            {
                var head = RegPhylipHead.Match(_line);
                if (!head.Success)
                {
                    StoreReset();
                    _seqData = null;
                    return false;
                }

                var phyitem = new AlnItem { Name = ItemName(head.Groups[1].Value) };
                AppendTo(phyitem.Seq, _line[(head.Index + head.Length)..]);
                ilen = phyitem.Seq.Length;
                if (ilen == len)
                    done = true;
                else if (ilen > len)
                {
                    StoreReset();
                    _seqData = null;
                    return false;
                }

                if (IsWhite(phyitem.Name) || data.Table.ContainsKey(phyitem.Name))
                {
                    _buff.SetBuffered();
                    StoreReset();
                    _seqData = null;
                    return ReadPhylipnon(t);
                }

                data.Table[phyitem.Name] = phyitem;
                list.Add(phyitem.Name);
                if (jseq == 0)
                    maxlen = ilen;
                else if (ilen != maxlen)
                {
                    StoreReset();
                    _seqData = null;
                    if (ReadPhylipnon(t))
                        return true;
                    _context.Warn("phylip format length mismatch in header");
                    return false;
                }

                jseq++;
                if (jseq < iseq)
                    ok = ReadLine();
            }

            data.NameArr = new string?[Math.Min(iseq, 1 << 24)];
            for (int i = 0; i < data.NameArr.Length && i < list.Count; i++)
                data.NameArr[i] = list[i];

            if (ilen < len)
            {
                jseq = 0;
                while (ReadLine())
                {
                    string? nm = jseq < data.NameArr.Length ? data.NameArr[jseq] : null;
                    if (PhylipReadseq(_line, data, nm, len, ref ilen, out done))
                    {
                        if (jseq == 0)
                            maxlen = ilen;
                        else if (ilen != maxlen)
                        {
                            _buff.SetBuffered();
                            StoreReset();
                            _seqData = null;
                            return ReadPhylipnon(t);
                        }

                        jseq++;
                        if (jseq == iseq)
                            jseq = 0;
                        if (jseq == 0 && done)
                            break;
                        done = true;
                    }
                }

                if (!done)
                {
                    _buff.SetBuffered();
                    StoreReset();
                    _seqData = null;
                    return ReadPhylipnon(t);
                }

                if (jseq > 0)
                {
                    _buff.SetBuffered();
                    StoreReset();
                    _seqData = null;
                    return ReadPhylipnon(t);
                }
            }

            data.Nseq = (int)Math.Min(iseq, int.MaxValue);
            data.Count = 0;
        }

        return PhylipReturn(t, false);
    }
}
