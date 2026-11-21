using System.Globalization;

namespace Seqeron.Genomics.Alignment;

/// <summary>The part of an alignment object (<c>Alignment</c> of <c>postnuc.cc</c> / <c>postpro.cc</c>) that the extension routines use.</summary>
internal class MummerAlignment
{
    public long SA, SB, EA, EB;
    public List<long> Delta = [];
    public long DeltaApos;
    public long Errors, SimErrors, NonAlphas;
}

/// <summary>
/// <c>extendBackward</c> and <c>extendForward</c> of <c>postnuc.cc</c> and <c>postpro.cc</c> (the two files contain the same routines): extend an alignment backwards / forwards with the engine of <see cref="NucmerSwAlign"/>
/// and merge the delta values; shared by <see cref="NucmerPostnuc"/> (nucleotides) and <see cref="PromerPostpro"/> (amino acids, BLOSUM matrices).
/// </summary>
internal static class MummerExtension
{
    /// <summary>The engine settings of one run: <c>-t</c>, the break length, the banding and the matrix type (0 nucleotides, 1 - 3 BLOSUM 45 / 62 / 80).</summary>
    public readonly record struct Settings(bool ToSequenceEnd, int BreakLen, int Banding, int Matrix);

    public static bool ExtendBackward<T>(Settings st, List<T> alignments, int currAp, int targetAp, byte[] a, byte[] b)
        where T : MummerAlignment
    {
        bool targetReached;
        bool overflow = false, doubleFlag = false;
        var cur = alignments[currAp];
        uint mo = NucmerSwAlign.BackwardSearch;
        long targetA, targetB;
        bool noTarget = targetAp == alignments.Count;

        if (!noTarget)
        {
            targetA = alignments[targetAp].EA;
            targetB = alignments[targetAp].EB;
        }
        else
        {
            targetA = 1;
            targetB = 1;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (cur.SA - targetA + 1 > NucmerSwAlign.MaxAlignmentLength)
        {
            targetA = cur.SA - NucmerSwAlign.MaxAlignmentLength + 1;
            overflow = true;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (cur.SB - targetB + 1 > NucmerSwAlign.MaxAlignmentLength)
        {
            targetB = cur.SB - NucmerSwAlign.MaxAlignmentLength + 1;
            if (overflow)
                doubleFlag = true;
            else
                overflow = true;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (st.ToSequenceEnd && !doubleFlag)
            mo |= NucmerSwAlign.SeqEndBit;

        targetReached = NucmerSwAlign.AlignEngine(a, cur.SA, ref targetA, b, cur.SB, ref targetB, [], mo, st.BreakLen, st.Banding, st.Matrix);

        if (overflow || noTarget)
            targetReached = false;

        if (targetReached)
        {
            ExtendForward(st, alignments[targetAp], a, cur.SA, b, cur.SB, NucmerSwAlign.ForcedForwardAlign);
            alignments[targetAp].EA = cur.EA;
            alignments[targetAp].EB = cur.EB;
            alignments.RemoveAt(alignments.Count - 1);
        }
        else
        {
            long endA = cur.SA, endB = cur.SB;
            NucmerSwAlign.AlignEngine(a, targetA, ref endA, b, targetB, ref endB, cur.Delta, NucmerSwAlign.ForcedForwardAlign, st.BreakLen, st.Banding, st.Matrix);
            cur.SA = targetA;
            cur.SB = targetB;
            foreach (long d in cur.Delta)
                cur.DeltaApos += d > 0 ? d : Math.Abs(d) - 1;
        }

        return targetReached;
    }

    public static bool ExtendForward(Settings st, MummerAlignment cur, byte[] a, long targetA, byte[] b, long targetB, uint mo)
    {
        bool overflow = false, doubleFlag = false;
        int di = cur.Delta.Count;

        long valA = targetA - cur.EA + 1;
        long valB = targetB - cur.EB + 1;

        if (valA > NucmerSwAlign.MaxAlignmentLength)
        {
            targetA = cur.EA + NucmerSwAlign.MaxAlignmentLength - 1;
            overflow = true;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (valB > NucmerSwAlign.MaxAlignmentLength)
        {
            targetB = cur.EB + NucmerSwAlign.MaxAlignmentLength - 1;
            if (overflow)
                doubleFlag = true;
            else
                overflow = true;
            mo |= NucmerSwAlign.OptimalBit;
        }

        if (doubleFlag)
            mo &= ~NucmerSwAlign.SeqEndBit;

        bool targetReached = NucmerSwAlign.AlignEngine(a, cur.EA, ref targetA, b, cur.EB, ref targetB, cur.Delta, mo, st.BreakLen, st.Banding, st.Matrix);

        if (targetReached && overflow)
            targetReached = false;

        if (di < cur.Delta.Count)
        {
            valA = (cur.EA - cur.SA + 1) - cur.DeltaApos - 1;
            cur.Delta[di] += cur.Delta[di] > 0 ? valA : -valA;
            if (cur.Delta[di] == 0 || valA < 0)
                throw new FormatException("ERROR: failed to merge alignments at position " + cur.EA.ToString(CultureInfo.InvariantCulture) + "\n       Please file a bug report\n");

            for (int k = di; k < cur.Delta.Count; k++)
                cur.DeltaApos += cur.Delta[k] > 0 ? cur.Delta[k] : Math.Abs(cur.Delta[k]) - 1;
        }

        cur.EA = targetA;
        cur.EB = targetB;
        return targetReached;
    }
}
