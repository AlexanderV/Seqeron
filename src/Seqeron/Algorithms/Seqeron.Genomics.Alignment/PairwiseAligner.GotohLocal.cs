#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905, S2681, S4144
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public sealed partial class PairwiseAligner
{
    // ---------------------------------------------------------------------------------------------------------------- the local Gotoh recursion with gap scores the library's aligners refuse (R14-10b)
    // GOTOH_LOCAL_SCORE and GOTOH_LOCAL_ALIGN of _pairwisealigner.c (Biopython 1.88): the three states with the gap states floored at 0 in every step (so, with an open score above its extend score or a gap score above zero, a gap run
    // does not pay its opening: the model of the run-length recursion of the global mode does not hold), the gap states of the last row and column are 0, the end gap scores are not used (a local alignment has no end gaps).

    private sealed class GotohState(int nA, int nB)
    {
        public readonly double[,] M = new double[nA + 1, nB + 1], Ix = new double[nA + 1, nB + 1], Iy = new double[nA + 1, nB + 1];
        public readonly int[,] TraceM = new int[nA + 1, nB + 1], TraceIx = new int[nA + 1, nB + 1], TraceIy = new int[nA + 1, nB + 1];
    }

    private double GotohLocalScore(string a, string b)
    {
        int nA = a.Length, nB = b.Length;
        double oA = _open[0], eA = _extend[0], oB = _open[3], eB = _extend[3];
        var pair = PairScore();
        var M = new double[nA + 1, nB + 1];
        var Ix = new double[nA + 1, nB + 1];
        var Iy = new double[nA + 1, nB + 1];
        M[0, 0] = 0;
        Ix[0, 0] = Iy[0, 0] = NegativeInfinity;
        for (int j = 1; j <= nB; j++)
        {
            M[0, j] = NegativeInfinity;
            Ix[0, j] = NegativeInfinity;
            Iy[0, j] = 0;
        }

        double maximum = 0;
        for (int i = 1; i <= nA; i++)
        {
            M[i, 0] = NegativeInfinity;
            Ix[i, 0] = 0;
            Iy[i, 0] = NegativeInfinity;
            for (int j = 1; j <= nB; j++)
            {
                double score = Math.Max(M[i - 1, j - 1], Math.Max(Ix[i - 1, j - 1], Iy[i - 1, j - 1])) + pair(a[i - 1], b[j - 1]);
                if (score < 0)
                    score = 0;
                else if (score > maximum)
                    maximum = score;
                M[i, j] = score;
                if (i == nA || j == nB)
                {
                    Ix[i, j] = 0;
                    Iy[i, j] = 0;
                    continue;
                }

                score = Math.Max(M[i - 1, j] + oB, Math.Max(Ix[i - 1, j] + eB, Iy[i - 1, j] + oB));
                if (score < 0)
                    score = 0;
                else if (score > maximum)
                    maximum = score;
                Ix[i, j] = score;
                score = Math.Max(M[i, j - 1] + oA, Math.Max(Ix[i, j - 1] + oA, Iy[i, j - 1] + eA));
                if (score < 0)
                    score = 0;
                else if (score > maximum)
                    maximum = score;
                Iy[i, j] = score;
            }
        }

        return maximum;
    }

    private (double Score, List<AlignmentResult> Results) GotohLocalAlign(string a, string b)
    {
        int nA = a.Length, nB = b.Length;
        double eps = Epsilon;
        double oA = _open[0], eA = _extend[0], oB = _open[3], eB = _extend[3];
        var pair = PairScore();
        var st = new GotohState(nA, nB);
        var M = st.M;
        var Ix = st.Ix;
        var Iy = st.Iy;
        for (int i = 0; i <= nA; i++)
        {
            for (int j = 0; j <= nB; j++)
            {
                M[i, j] = 0;
                Ix[i, j] = Iy[i, j] = NegativeInfinity;
                st.TraceM[i, j] = TraceStart;
            }
        }

        double maximum = 0;
        int im = nA, jm = nB;
        for (int i = 1; i <= nA; i++)
        {
            for (int j = 1; j <= nB; j++)
            {
                // SELECT_TRACE_GOTOH_LOCAL_ALIGN
                int trace = TraceM;
                double score = M[i - 1, j - 1];
                double t = Ix[i - 1, j - 1];
                if (t > score + eps) { score = t; trace = TraceIx; }
                else if (t > score - eps) trace |= TraceIx;
                t = Iy[i - 1, j - 1];
                if (t > score + eps) { score = t; trace = TraceIy; }
                else if (t > score - eps) trace |= TraceIy;
                score += pair(a[i - 1], b[j - 1]);
                if (score < eps)
                {
                    score = 0;
                    trace = TraceStart;
                }
                else if (score > maximum - eps)
                {
                    if (score > maximum + eps)
                    {
                        maximum = score;
                        for (; im < i; im++, jm = 0)
                            for (; jm <= nB; jm++)
                                st.TraceM[im, jm] &= ~TraceEnd;
                        for (; jm < j; jm++)
                            st.TraceM[im, jm] &= ~TraceEnd;
                        im = i;
                        jm = j;
                    }

                    trace |= TraceEnd;
                }

                M[i, j] = score;
                st.TraceM[i, j] = trace;
                if (i == nA || j == nB)
                {
                    Ix[i, j] = 0;
                    Iy[i, j] = 0;
                    st.TraceIx[i, j] = 0;
                    st.TraceIy[i, j] = 0;
                    continue;
                }

                (Ix[i, j], st.TraceIx[i, j]) = GotohGap(M[i - 1, j] + oB, Ix[i - 1, j] + eB, Iy[i - 1, j] + oB, eps);
                (Iy[i, j], st.TraceIy[i, j]) = GotohGap(M[i, j - 1] + oA, Ix[i, j - 1] + oA, Iy[i, j - 1] + eA, eps);
            }
        }

        // no zero-score extensions: the reachability pass (the end of GOTOH_LOCAL_ALIGN)
        var path = new int[nA + 1, nB + 1];
        for (int j = 0; j <= nB; j++)
            path[0, j] = TraceM;
        for (int i = 1; i <= nA; i++)
        {
            path[i, 0] = TraceM;
            for (int j = 1; j <= nB; j++)
            {
                int trace = st.TraceM[i, j];
                if ((path[i - 1, j - 1] & TraceM) == 0) trace &= ~TraceM;
                if ((path[i - 1, j - 1] & TraceIx) == 0) trace &= ~TraceIx;
                if ((path[i - 1, j - 1] & TraceIy) == 0) trace &= ~TraceIy;
                if ((trace & (TraceStart | TraceM | TraceIx | TraceIy)) != 0)
                {
                    if ((trace & TraceEnd) != 0)
                        path[i, j] = 0;
                    else
                        path[i, j] |= TraceM;
                }
                else
                {
                    path[i, j] &= ~TraceM;
                    trace = 0;
                }

                st.TraceM[i, j] = trace;
                trace = st.TraceIx[i, j];
                if ((path[i - 1, j] & TraceM) == 0) trace &= ~TraceM;
                if ((path[i - 1, j] & TraceIx) == 0) trace &= ~TraceIx;
                if ((path[i - 1, j] & TraceIy) == 0) trace &= ~TraceIy;
                if ((trace & (TraceM | TraceIx | TraceIy)) != 0)
                    path[i, j] |= TraceIx;
                else
                {
                    path[i, j] &= ~TraceIx;
                    trace = 0;
                }

                st.TraceIx[i, j] = trace;
                trace = st.TraceIy[i, j];
                if ((path[i, j - 1] & TraceM) == 0) trace &= ~TraceM;
                if ((path[i, j - 1] & TraceIx) == 0) trace &= ~TraceIx;
                if ((path[i, j - 1] & TraceIy) == 0) trace &= ~TraceIy;
                if ((trace & (TraceM | TraceIx | TraceIy)) != 0)
                    path[i, j] |= TraceIy;
                else
                {
                    path[i, j] &= ~TraceIy;
                    trace = 0;
                }

                st.TraceIy[i, j] = trace;
            }
        }

        var results = new List<AlignmentResult>();
        var steps = new List<(int State, int I, int J)>();           // from the end of the alignment to its start
        void Finish(int startI, int startJ, int endI, int endJ)
        {
            var ra = new System.Text.StringBuilder();
            var rb = new System.Text.StringBuilder();
            for (int n = steps.Count - 1; n >= 0; n--)
            {
                var (state, i, j) = steps[n];
                switch (state)
                {
                    case TraceM:
                        ra.Append(a[i - 1]);
                        rb.Append(b[j - 1]);
                        break;
                    case TraceIx:
                        ra.Append(a[i - 1]);
                        rb.Append('-');
                        break;
                    default:
                        ra.Append('-');
                        rb.Append(b[j - 1]);
                        break;
                }
            }

            results.Add(new AlignmentResult(ra.ToString(), rb.ToString(), (int)Math.Round(maximum, MidpointRounding.AwayFromZero), AlignmentType.Local, startI, startJ, endI - 1, endJ - 1) { ExactScore = maximum });
        }

        bool Walk(int state, int i, int j, int endI, int endJ)
        {
            if (results.Count >= MaxAlignments)
                return false;
            int[,] table = state == TraceM ? st.TraceM : state == TraceIx ? st.TraceIx : st.TraceIy;
            int trace = table[i, j];
            if (state == TraceM && trace == TraceStart)
            {
                Finish(i, j, endI, endJ);
                return true;
            }

            if (i == 0 || j == 0)
                return true;
            int pi = state == TraceIy ? i : i - 1, pj = state == TraceIx ? j : j - 1;
            steps.Add((state, i, j));
            bool ok = true;
            foreach (int prev in new[] { TraceM, TraceIx, TraceIy })
            {
                if ((trace & prev) == 0)
                    continue;
                if (!Walk(prev, pi, pj, endI, endJ))
                {
                    ok = false;
                    break;
                }
            }

            steps.RemoveAt(steps.Count - 1);
            return ok;
        }

        if (maximum > 0)
        {
            for (int i = 0; i <= nA && results.Count < MaxAlignments; i++)
                for (int j = 0; j <= nB; j++)
                    if ((st.TraceM[i, j] & TraceEnd) != 0)
                        Walk(TraceM, i, j, i, j);
        }

        return (maximum, results);
    }

    // the Smith-Waterman recursion (SMITHWATERMAN_SCORE / SMITHWATERMAN_ALIGN: every gap score equal for opening and extending): one matrix, the flags of a cell are the moves that reach it
    private const int MoveDiagonal = 1, MoveHorizontal = 2, MoveVertical = 4;

    private bool UsesSmithWaterman => UsesGotohLocal && !IsAffine;

    private double SmithWatermanScore(string a, string b)
    {
        int nA = a.Length, nB = b.Length;
        double gA = _extend[0], gB = _extend[3];
        var pair = PairScore();
        var row = new double[nB + 1];
        double maximum = 0;
        for (int i = 1; i <= nA; i++)
        {
            double temp = 0;
            for (int j = 1; j <= nB; j++)
            {
                double score = temp + pair(a[i - 1], b[j - 1]);
                if (i < nA && j < nB)
                {
                    score = Math.Max(score, Math.Max(row[j - 1] + gA, row[j] + gB));
                }

                if (score < 0)
                    score = 0;
                else if (score > maximum)
                    maximum = score;
                temp = row[j];
                row[j] = score;
            }
        }

        return maximum;
    }

    private (double Score, List<AlignmentResult> Results) SmithWatermanAlign(string a, string b)
    {
        int nA = a.Length, nB = b.Length;
        double eps = Epsilon;
        double gA = _extend[0], gB = _extend[3];
        var pair = PairScore();
        var trace = new int[nA + 1, nB + 1];
        var row = new double[nB + 1];
        for (int j = 0; j <= nB; j++)
            trace[0, j] = TraceStart;
        double maximum = 0;
        int im = nA, jm = nB;
        for (int i = 1; i <= nA; i++)
        {
            double temp = 0;
            trace[i, 0] = TraceStart;
            for (int j = 1; j <= nB; j++)
            {
                int t = MoveDiagonal;
                double score = temp + pair(a[i - 1], b[j - 1]);
                if (i < nA && j < nB)
                {
                    double h = row[j - 1] + gA;
                    if (h > score + eps) { score = h; t = MoveHorizontal; }
                    else if (h > score - eps) t |= MoveHorizontal;
                    double v = row[j] + gB;
                    if (v > score + eps) { score = v; t = MoveVertical; }
                    else if (v > score - eps) t |= MoveVertical;
                }

                if (score < eps)
                {
                    score = 0;
                    t = TraceStart;
                }
                else if ((t & MoveDiagonal) != 0 && score > maximum - eps)
                {
                    if (score > maximum + eps)
                    {
                        for (; im < i; im++, jm = 0)
                            for (; jm <= nB; jm++)
                                trace[im, jm] &= ~TraceEnd;
                        for (; jm < j; jm++)
                            trace[im, jm] &= ~TraceEnd;
                        im = i;
                        jm = j;
                    }

                    t |= TraceEnd;
                }

                trace[i, j] = t;
                if (score > maximum)
                    maximum = score;
                temp = row[j];
                row[j] = score;
            }
        }

        // no zero-score extensions (the reachability pass)
        var reach = new bool[nA + 1, nB + 1];
        for (int j = 0; j <= nB; j++)
            reach[0, j] = true;
        for (int i = 1; i <= nA; i++)
        {
            reach[i, 0] = true;
            for (int j = 1; j <= nB; j++)
            {
                int t = trace[i, j];
                if (!reach[i - 1, j - 1]) t &= ~MoveDiagonal;
                if (!reach[i, j - 1]) t &= ~MoveHorizontal;
                if (!reach[i - 1, j]) t &= ~MoveVertical;
                if ((t & (TraceStart | MoveHorizontal | MoveVertical | MoveDiagonal)) != 0)
                    reach[i, j] = (t & TraceEnd) == 0;
                else
                {
                    reach[i, j] = false;
                    t = 0;
                }

                trace[i, j] = t;
            }
        }

        var results = new List<AlignmentResult>();
        var steps = new List<(int Move, int I, int J)>();
        void Finish(int startI, int startJ, int endI, int endJ)
        {
            var ra = new System.Text.StringBuilder();
            var rb = new System.Text.StringBuilder();
            for (int n = steps.Count - 1; n >= 0; n--)
            {
                var (move, i, j) = steps[n];
                if (move == MoveDiagonal)
                {
                    ra.Append(a[i - 1]);
                    rb.Append(b[j - 1]);
                }
                else if (move == MoveVertical)
                {
                    ra.Append(a[i - 1]);
                    rb.Append('-');
                }
                else
                {
                    ra.Append('-');
                    rb.Append(b[j - 1]);
                }
            }

            results.Add(new AlignmentResult(ra.ToString(), rb.ToString(), (int)Math.Round(maximum, MidpointRounding.AwayFromZero), AlignmentType.Local, startI, startJ, endI - 1, endJ - 1) { ExactScore = maximum });
        }

        bool Walk(int i, int j, int endI, int endJ, bool first = false)
        {
            if (results.Count >= MaxAlignments)
                return false;
            int t = trace[i, j];
            if (first)
                t &= ~(MoveHorizontal | MoveVertical);                       // an alignment ends with a pair: only the diagonal move leaves an ENDPOINT
            if (t == TraceStart)
            {
                Finish(i, j, endI, endJ);
                return true;
            }

            if ((t & (MoveDiagonal | MoveHorizontal | MoveVertical)) == 0)
                throw new InvalidOperationException("Unexpected trace in PathGenerator_next_smithwaterman");      // the reference raises RuntimeError at a cell without a move
            bool ok = true;
            foreach (int move in new[] { MoveDiagonal, MoveHorizontal, MoveVertical })
            {
                if ((t & move) == 0)
                    continue;
                steps.Add((move, i, j));
                ok = Walk(move == MoveHorizontal ? i : i - 1, move == MoveVertical ? j : j - 1, endI, endJ);
                steps.RemoveAt(steps.Count - 1);
                if (!ok)
                    break;
            }

            return ok;
        }

        if (maximum > 0)
        {
            for (int i = 0; i <= nA && results.Count < MaxAlignments; i++)
                for (int j = 0; j <= nB; j++)
                    if ((trace[i, j] & TraceEnd) != 0)
                        Walk(i, j, i, j, first: true);
        }

        return (maximum, results);
    }

    // SELECT_TRACE_GOTOH_LOCAL_GAP
    private static (double Score, int Trace) GotohGap(double fromM, double fromIx, double fromIy, double eps)
    {
        int trace = TraceM;
        double score = fromM;
        if (fromIx > score + eps) { score = fromIx; trace = TraceIx; }
        else if (fromIx > score - eps) trace |= TraceIx;
        if (fromIy > score + eps) { score = fromIy; trace = TraceIy; }
        else if (fromIy > score - eps) trace |= TraceIy;
        if (score < eps)
            return (NegativeInfinity, 0);
        return (score, trace);
    }
}
