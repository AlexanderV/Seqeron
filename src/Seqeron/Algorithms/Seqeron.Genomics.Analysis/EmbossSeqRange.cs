using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// The range and orientation of a sequence object of the AJAX library: <c>ajSeqSetRange</c>, <c>ajSeqSetRangeRev</c>, <c>ajSeqReverseDo</c>, <c>ajSeqGetBegin</c>, <c>ajSeqGetEnd</c> and the cut of <c>seqClone</c> that the writers apply.
/// A5-8c3h, B09 F255.
/// </summary>
internal static class EmbossSeqRange
{
    /// <summary><c>ajSeqSetRange</c>: the range of the stream is used where the sequence has none of its own (the one of its USA); a sequence flagged for reversal is reversed.</summary>
    public static void SetRange(EmbossSeq seq, int pos1, int pos2)
    {
        if (pos1 != 0 && seq.Begin == 0)
            seq.Begin = pos1;
        if (pos2 != 0 && seq.End == 0)
            seq.End = pos2;
        if (seq.Rev && !seq.Reversed)
            ReverseDo(seq);
    }

    /// <summary><c>ajSeqSetRangeRev</c>: as <see cref="SetRange"/>, the sequence is reversed.</summary>
    public static void SetRangeRev(EmbossSeq seq, int pos1, int pos2)
    {
        if (pos1 != 0 && seq.Begin == 0)
            seq.Begin = pos1;
        if (pos2 != 0 && seq.End == 0)
            seq.End = pos2;
        if (!seq.Rev)
        {
            seq.Rev = true;
            seq.Reversed = false;
        }

        ReverseDo(seq);
    }

    /// <summary><c>ajSeqReverseDo</c>: the sequence is reverse-complemented, the range is mirrored and its feature table reversed.</summary>
    public static void ReverseDo(EmbossSeq seq)
    {
        if (!seq.Rev)
            return;
        int ibegin = seq.Begin, iend = seq.End;
        seq.End = -ibegin;
        seq.Begin = -iend;
        (seq.Offset, seq.Offend) = (seq.Offend, seq.Offset);
        seq.Rev = false;
        seq.Reversed = !seq.Reversed;
        string reversed = ReverseComplement(seq.Seq.ToString());
        seq.Seq.Clear().Append(reversed);
        if (seq.Fttable is { } table)
        {
            table.SetReverse();
            table.ReverseTable();
        }
    }

    // ajSeqstrReverse: ajBaseAlphacharComp of every character ("ACGTURYWSMKBDHVNXacgturywsmkbdhvnx" -> "TGCAAYRWSKMVHDBNXtgcaayrwskmvhdbnx", every other character is kept) from the end
    private static string ReverseComplement(string s)
    {
        const string fwd = "ACGTURYWSMKBDHVNXacgturywsmkbdhvnx", rev = "TGCAAYRWSKMVHDBNXtgcaayrwskmvhdbnx";
        return string.Create(s.Length, s, static (dest, src) =>
        {
            for (int i = 0; i < src.Length; i++)
            {
                char c = src[src.Length - 1 - i];
                int k = fwd.IndexOf(c, StringComparison.Ordinal);
                dest[i] = k >= 0 ? rev[k] : c;
            }
        });
    }

    /// <summary>
    /// <c>ajSeqTrim</c>: the sequence is cut to its range (a program that works on the residues of the range: union ...); the offsets and the range of the feature table follow. Nothing is cut when the end is before the begin.
    /// </summary>
    public static void Trim(EmbossSeq seq)
    {
        int jbegin = seq.Begin > 0 ? seq.Begin - 1 : seq.Begin;
        int jend = seq.End > 0 ? seq.End - 1 : seq.End;
        int len = seq.Seq.Length;
        int begin = 1 + EmbossInputSequence.PositionStart(len, 0, jbegin);
        int end = 1 + EmbossInputSequence.PositionStart(len, begin - 1, jend);
        if (seq.End != 0)
        {
            if (end < begin)
                return;
            seq.Offend = len - end;
            seq.Seq.Length = end;
            seq.End = 0;
        }

        if (seq.Begin != 0)
        {
            seq.Seq.Remove(0, begin - 1);
            seq.Offset = begin - 1;
            seq.Begin = 0;
        }

        if (seq.Fttable is { } table)
            table.TrimOff(seq.Offset, seq.Seq.Length);
        if (seq.Rev)
            ReverseDo(seq);
    }

    /// <summary><c>ajSeqNewSeq</c> as <c>seqWriteListAppend</c> uses it: a copy of the sequence with its range (the feature table is not copied: the writers of the lists do not read it).</summary>
    public static EmbossSeq CopyForList(EmbossSeq seq)
    {
        var c = new EmbossSeq
        {
            Name = seq.Name, Desc = seq.Desc, Type = seq.Type, Begin = seq.Begin, End = seq.End, Offset = seq.Offset, Offend = seq.Offend, Rev = seq.Rev, Reversed = seq.Reversed, Weight = seq.Weight,
        };
        c.Seq.Append(seq.Seq);
        return c;
    }

    /// <summary><c>ajSeqGetBegin</c>: the 1-based start of the range in the sequence as it is now.</summary>
    public static int GetBegin(EmbossSeq seq)
    {
        if (seq.Begin == 0)
            return 1;
        int i = seq.Begin > 0 ? seq.Begin - 1 : seq.Begin;
        return 1 + EmbossInputSequence.PositionStart(seq.Seq.Length, 0, i);
    }

    /// <summary><c>ajSeqGetEnd</c>.</summary>
    public static int GetEnd(EmbossSeq seq)
    {
        if (seq.End == 0)
            return seq.Seq.Length;
        int i = seq.End > 0 ? seq.End - 1 : seq.End;
        return 1 + EmbossInputSequence.PositionStart(seq.Seq.Length, GetBegin(seq) - 1, i);
    }

    /// <summary>
    /// The cut of <c>seqClone</c>: the residues of the range of the sequence (the whole sequence when the range does not cut it) and the range as the writer sees it. The feature table of the sequence is given its range
    /// (<c>ajFeattableSetRange</c>) when the sequence has a begin or an end.
    /// </summary>
    public static string Cut(EmbossSeq seq, out int ibegin, out int iend, out bool ranged)
    {
        int len = seq.Seq.Length;
        ibegin = 1;
        iend = len;
        bool seqrange = false;
        if (seq.Begin != 0)
        {
            ibegin = GetBegin(seq);
            if (ibegin > 1)
                seqrange = true;
        }

        if (seq.End != 0)
        {
            iend = GetEnd(seq);
            if (iend < len)
                seqrange = true;
        }

        ranged = seq.Begin != 0 || seq.End != 0;
        string text = seq.Seq.ToString();
        string residues = "";
        if (iend >= ibegin)
            residues = seqrange ? text.Substring(ibegin - 1, iend - ibegin + 1) : text;
        if (seq.Fttable is { } table)
        {
            if (seq.Rev)
                table.SetReverse();
            if (seq.Circular)
                table.Circular = true;
            if (ranged)
                table.SetRange(ibegin, iend);
        }

        return residues;
    }
}
