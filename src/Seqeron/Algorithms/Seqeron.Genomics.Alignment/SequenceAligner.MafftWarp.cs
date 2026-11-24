#pragma warning disable S3626, S3871, S1643, S2325, IDE0052, IDE0051, IDE0059, S127, S3358, S1066, S4487, S1144, CA1823, CS0414, CA1870, S1125, S3776, S107, S138, S1541, S1854
namespace Seqeron.Genomics.Alignment;

public static partial class SequenceAligner
{
    #region MAFFT 7.505: --allowshift / --shiftpenalty (trywarp: the jump of the dynamic programs to the best cell of an earlier part of the matrix)

    /// <summary>
    /// The arrays that <c>A__align</c> / <c>G__align11</c> / <c>partA__align</c> keep while <c>trywarp</c> is on (<c>penalty_shift_factor &lt; 10</c>, <c>--allowshift</c> or <c>--shiftpenalty</c> below 10): for every column the best value of
    /// a match cell above it (<c>wmrecords</c>) with its position (<c>warpi</c>, <c>warpj</c>), the same for the row before (<c>prevwmrecords</c>, <c>prevwarpi</c>, <c>prevwarpj</c>), and the list of the jump targets that the trace-back reads
    /// (<c>warpis</c>, <c>warpjs</c>; a cell of the pointer matrix that holds <c>warpbase + n</c> jumps to the n-th target).
    /// </summary>
    internal sealed class WarpState
    {
        internal int Base { get; }

        internal List<int> Is { get; } = new();

        internal List<int> Js { get; } = new();

        private readonly double shift;
        private readonly double[] wmrecords, prevwmrecords;
        private readonly int[] warpi, warpj, prevwarpi, prevwarpj;

        internal WarpState(int lgth1, int lgth2, int penaltyShift)
        {
            Base = lgth1 + lgth2;
            shift = penaltyShift;
            wmrecords = new double[lgth2 + 1];
            prevwmrecords = new double[lgth2 + 1];
            warpi = new int[lgth2 + 1];
            warpj = new int[lgth2 + 1];
            prevwarpi = new int[lgth2 + 1];
            prevwarpj = new int[lgth2 + 1];
            Array.Fill(prevwarpi, -Base);
            Array.Fill(prevwarpj, -Base);
            Array.Fill(warpi, -Base);
            Array.Fill(warpj, -Base);
        }

        /// <summary>
        /// The part of the loop over the columns that follows the gap moves (cell <paramref name="i"/>, <paramref name="j"/>): the jump from the best cell of the row before is taken when it beats <paramref name="wm"/> (the pointer
        /// <c>warpbase + n</c> goes to <paramref name="ijpRow"/>[<paramref name="j"/>]), then the records of the column are updated with the match value <paramref name="cur"/> plus the new <c>wm</c>. Returns the new <c>wm</c>.
        /// </summary>
        internal double Step(int i, int j, double wm, double cur, double penaltyEx, int[] ijpRow)
        {
            double g = prevwmrecords[j - 1] + shift + penaltyEx * (i - prevwarpi[j - 1] + j - prevwarpj[j - 1]);
            if (g > wm)
            {
                int n = Is.Count;
                if (n > 0 && prevwarpi[j - 1] == Is[n - 1] && prevwarpj[j - 1] == Js[n - 1])
                {
                    ijpRow[j] = Base + n - 1;
                }
                else
                {
                    ijpRow[j] = Base + n;
                    Is.Add(prevwarpi[j - 1]);
                    Js.Add(prevwarpj[j - 1]);
                }
                wm = g;
            }
            double curm = cur + wm;
            if (wmrecords[j - 1] > wmrecords[j])
            {
                wmrecords[j] = wmrecords[j - 1];
                warpi[j] = warpi[j - 1];
                warpj[j] = warpj[j - 1];
            }
            if (curm > wmrecords[j])
            {
                wmrecords[j] = curm;
                warpi[j] = i;
                warpj[j] = j;
            }
            return wm;
        }

        /// <summary><see cref="Step"/> of <c>D__align</c> (<c>Dalignmm.c</c>): the jump costs the shift penalty only (no extension term for the distance of the best cell).</summary>
        internal double StepD(int i, int j, double wm, double cur, double penaltyShift, int[] ijpRow)
        {
            double g = prevwmrecords[j - 1] + penaltyShift;
            if (g > wm)
            {
                int n = Is.Count;
                if (n > 0 && prevwarpi[j - 1] == Is[n - 1] && prevwarpj[j - 1] == Js[n - 1])
                {
                    ijpRow[j] = Base + n - 1;
                }
                else
                {
                    ijpRow[j] = Base + n;
                    Is.Add(prevwarpi[j - 1]);
                    Js.Add(prevwarpj[j - 1]);
                }
                wm = g;
            }
            double curm = cur + wm;
            if (wmrecords[j - 1] > wmrecords[j])
            {
                wmrecords[j] = wmrecords[j - 1];
                warpi[j] = warpi[j - 1];
                warpj[j] = warpj[j - 1];
            }
            if (curm > wmrecords[j])
            {
                wmrecords[j] = curm;
                warpi[j] = i;
                warpj[j] = j;
            }
            return wm;
        }

        /// <summary>The end of a row: the records become the records of the row before (<c>fltncpy</c> / <c>intncpy</c> of <paramref name="count"/> entries).</summary>
        internal void EndRow(int count)
        {
            Array.Copy(wmrecords, prevwmrecords, count);
            Array.Copy(warpi, prevwarpi, count);
            Array.Copy(warpj, prevwarpj, count);
        }
    }

    #endregion
}
