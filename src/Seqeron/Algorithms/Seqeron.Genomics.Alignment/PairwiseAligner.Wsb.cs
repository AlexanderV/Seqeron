#pragma warning disable S3776, S107, S1854, S1066, S3358, S1125, S127, S1905, S2681, S4144
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

public sealed partial class PairwiseAligner
{
    // ---------------------------------------------------------------------------------------------------------------- Waterman-Smith-Beyer (R14-10b)
    // A port of the three-state recursion of _pairwisealigner.c of Biopython 1.88 (WATERMANSMITHBEYER_GLOBAL_SCORE / LOCAL_SCORE, the trace macros SELECT_TRACE_WATERMAN_SMITH_BEYER_*, PathGenerator_next_waterman_smith_beyer_*):
    //   M(i,j)  = max(M, Ix, Iy)(i-1,j-1) + s(a_i, b_j)
    //   Ix(i,j) = max over k = 1..i of max(M, Iy)(i-k, j) + deletion(j, k)        (k residues of the target against a gap)
    //   Iy(i,j) = max over k = 1..j of max(M, Ix)(i, j-k) + insertion(i, k)        (k residues of the query against a gap)
    // every gap of length k is one run whose score is a function of (position, length) - the end gaps too. Used for gap-score functions (Biopython: Waterman-Smith-Beyer) and for the gap scores above zero (which the Gotoh /
    // Needleman-Wunsch / Smith-Waterman code of the library refuses: the recursion with the run scores open + (k - 1) * extend is the same model).

    private const double NegativeInfinity = -double.MaxValue;
    private const int TraceM = 1, TraceIx = 2, TraceIy = 4, TraceStart = 8, TraceEnd = 16;

    private bool HasFunctions => InsertionScoreFunction is not null || DeletionScoreFunction is not null;

    // the settings the library's Needleman-Wunsch / Gotoh / Smith-Waterman code refuses: a gap score above zero, or an open score above its extend score
    private bool NeedsGeneralGaps
    {
        get
        {
            for (int i = 0; i < 6; i++)
                if (_open[i] > 0 || _extend[i] > 0 || _open[i] > _extend[i])
                    return true;
            return false;
        }
    }

    // the global recursion with run scores: functions, and gap scores the Gotoh code of the library refuses; the local mode of such Gotoh settings has its own recursion (floors in the gap states)
    private bool UsesWaterman => HasFunctions || (NeedsGeneralGaps && Mode == PairwiseAlignerMode.Global);

    private bool UsesGotohLocal => !HasFunctions && NeedsGeneralGaps && Mode == PairwiseAlignerMode.Local;

    // the score of a gap of length k at position i of the target (insertion) / position p of the query (deletion), nA / nB the lengths
    private Func<int, int, double> InsertionScoreFor(int nA, char strand)
    {
        if (InsertionScoreFunction is { } f)
            return (i, k) => f(i, k);
        // the Waterman-Smith-Beyer macros do not swap the ends of the target for strand '-', the Gotoh / Needleman-Wunsch macros do
        int left = strand == '+' || HasFunctions ? 1 : 2, right = 3 - left;
        return (i, k) =>
        {
            int slot = i == 0 ? left : i == nA ? right : 0;
            return _open[slot] + (k - 1) * _extend[slot];
        };
    }

    private Func<int, int, double> DeletionScoreFor(int nB, char strand)
    {
        if (DeletionScoreFunction is { } f)
            return (p, k) => f(p, k);
        int left = strand == '+' || HasFunctions ? 4 : 5, right = 9 - left;
        return (p, k) =>
        {
            int slot = p == 0 ? left : p == nB ? right : 3;
            return _open[slot] + (k - 1) * _extend[slot];
        };
    }

    private Func<char, char, double> PairScore()
    {
        if (SubstitutionMatrix is { } matrix)
            return (x, y) => (double)matrix.GetValue(x, y);
        char? wildcard = Wildcard is { } w ? char.ToUpperInvariant(w) : null;
        double match = MatchScore, mismatch = MismatchScore;
        return (x, y) => wildcard is { } c && (x == c || y == c) ? 0 : x == y ? match : mismatch;
    }

    private sealed class WsbState(int nA, int nB)
    {
        public readonly double[,] M = new double[nA + 1, nB + 1], Ix = new double[nA + 1, nB + 1], Iy = new double[nA + 1, nB + 1];
        public readonly int[,] Trace = new int[nA + 1, nB + 1];
        // the gap lengths (k) a gap state may come from: from M and from the other gap state (ended by the caller's reading: lists, not 0-terminated)
        public readonly List<int>?[,] MIx = new List<int>?[nA + 1, nB + 1], IyIx = new List<int>?[nA + 1, nB + 1], MIy = new List<int>?[nA + 1, nB + 1], IxIy = new List<int>?[nA + 1, nB + 1];
        public double Maximum;
        // Waterman-Smith-Beyer local alignment raises the maximum with its gap cells; the Gotoh / Smith-Waterman code of Biopython does not (while its score() does)
        public bool GapsRaiseMaximum = true;
    }

    // the recursion with the traces: local = the local variant (the gap states of the last row and column are 0, every state is floored at 0)
    private WsbState RunWsb(string a, string b, char strand, bool local, bool traces)
    {
        int nA = a.Length, nB = b.Length;
        double eps = Epsilon;
        var st = new WsbState(nA, nB) { GapsRaiseMaximum = !traces || InsertionScoreFunction is not null || DeletionScoreFunction is not null };
        var ins = InsertionScoreFor(nA, strand);
        var del = DeletionScoreFor(nB, strand);
        var pair = PairScore();
        int Pos(int j) => strand == '+' || !HasFunctions ? j : nB - j;                // query_gap_start (the macros of the Gotoh code swap the slots instead)
        var M = st.M;
        var Ix = st.Ix;
        var Iy = st.Iy;
        M[0, 0] = 0;
        Ix[0, 0] = Iy[0, 0] = NegativeInfinity;
        for (int i = 1; i <= nA; i++)
        {
            if (local && !traces)
            {
                // WATERMANSMITHBEYER_LOCAL_SCORE: M is -inf on the first column, the gap state Ix is 0
                M[i, 0] = NegativeInfinity;
                Ix[i, 0] = 0;
                Iy[i, 0] = NegativeInfinity;
            }
            else if (local)
            {
                M[i, 0] = 0;
                Ix[i, 0] = Iy[i, 0] = NegativeInfinity;
            }
            else
            {
                M[i, 0] = NegativeInfinity;
                Iy[i, 0] = NegativeInfinity;
                Ix[i, 0] = del(Pos(0), i);
            }
        }

        for (int j = 1; j <= nB; j++)
        {
            if (local && !traces)
            {
                M[0, j] = NegativeInfinity;
                Ix[0, j] = NegativeInfinity;
                Iy[0, j] = 0;
            }
            else if (local)
            {
                M[0, j] = 0;
                Ix[0, j] = Iy[0, j] = NegativeInfinity;
            }
            else
            {
                M[0, j] = NegativeInfinity;
                Ix[0, j] = NegativeInfinity;
                Iy[0, j] = ins(0, j);
            }
        }

        if (traces)
        {
            // the boundary cells of the traces (PathGenerator_create_WSB)
            for (int i = 0; i <= nA; i++)
                for (int j = 0; j <= nB; j++)
                    st.Trace[i, j] = local ? TraceStart : 0;
            if (!local)
                st.Trace[0, 0] = 0;
        }

        double maximum = 0;
        int im = nA, jm = nB;
        for (int i = 1; i <= nA; i++)
        {
            char ca = a[i - 1];
            for (int j = 1; j <= nB; j++)
            {
                double s = pair(ca, b[j - 1]);
                if (!local)
                {
                    int trace = TraceM;
                    double score = M[i - 1, j - 1];
                    double t = Ix[i - 1, j - 1];
                    if (t > score + eps) { score = t; trace = TraceIx; }
                    else if (t > score - eps) trace |= TraceIx;
                    t = Iy[i - 1, j - 1];
                    if (t > score + eps) { score = t; trace = TraceIy; }
                    else if (t > score - eps) trace |= TraceIy;
                    M[i, j] = score + s;
                    st.Trace[i, j] = trace;
                    GapRuns(st, i, j, nA, nB, traces, ins, del, Pos(j), local: false, eps, ref maximum);
                }
                else if (!traces)
                {
                    // WATERMANSMITHBEYER_LOCAL_SCORE with SELECT_SCORE_GOTOH_LOCAL_ALIGN
                    double score = Math.Max(M[i - 1, j - 1], Math.Max(Ix[i - 1, j - 1], Iy[i - 1, j - 1])) + s;
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

                    GapRuns(st, i, j, nA, nB, false, ins, del, Pos(j), local: true, eps, ref maximum);
                }
                else
                {
                    int trace = TraceM;
                    double score = M[i - 1, j - 1];
                    double t = Ix[i - 1, j - 1];
                    if (t > score + eps) { score = t; trace = TraceIx; }
                    else if (t > score - eps) trace |= TraceIx;
                    t = Iy[i - 1, j - 1];
                    if (t > score + eps) { score = t; trace = TraceIy; }
                    else if (t > score - eps) trace |= TraceIy;
                    score += s;
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
                            if (traces)
                            {
                                for (; im < i; im++, jm = 0)
                                    for (; jm <= nB; jm++)
                                        st.Trace[im, jm] &= ~TraceEnd;
                                for (; jm < j; jm++)
                                    st.Trace[im, jm] &= ~TraceEnd;
                                im = i;
                                jm = j;
                            }
                        }

                        trace |= TraceEnd;
                    }

                    M[i, j] = score;
                    st.Trace[i, j] = trace;
                    if (i == nA || j == nB)
                    {
                        Ix[i, j] = score;                                           // the library's macro stores the cell score here (the traces keep no gap lists)
                        Iy[i, j] = score;
                        continue;
                    }

                    GapRuns(st, i, j, nA, nB, traces, ins, del, Pos(j), local: true, eps, ref maximum);
                }
            }
        }

        st.Maximum = maximum;
        if (local && traces)
            PruneLocalTraces(st, nA, nB);
        return st;
    }

    // the end of WATERMANSMITHBEYER_LOCAL_ALIGN: "As we don't allow zero-score extensions to alignments, we need to remove all traces towards an ENDPOINT": a state reachable only through an ENDPOINT cell loses its traces, a cell that no
    // STARTPOINT reaches is no ENDPOINT
    private static void PruneLocalTraces(WsbState st, int nA, int nB)
    {
        var path = new int[nA + 1, nB + 1];
        for (int j = 0; j <= nB; j++)
            path[0, j] = TraceM;
        for (int i = 1; i <= nA; i++)
        {
            path[i, 0] = TraceM;
            for (int j = 1; j <= nB; j++)
            {
                int trace = st.Trace[i, j];
                if ((path[i - 1, j - 1] & TraceM) == 0) trace &= ~TraceM;
                if ((path[i - 1, j - 1] & TraceIx) == 0) trace &= ~TraceIx;
                if ((path[i - 1, j - 1] & TraceIy) == 0) trace &= ~TraceIy;
                if ((trace & (TraceStart | TraceM | TraceIx | TraceIy)) != 0)
                {
                    if ((trace & TraceEnd) != 0)
                        path[i, j] = 0;                                              // no extensions after an ENDPOINT
                    else
                        path[i, j] |= TraceM;
                }
                else
                {
                    path[i, j] &= ~TraceM;
                    trace = 0;
                }

                st.Trace[i, j] = trace;
                if (i == nA || j == nB)
                    continue;
                var mix = st.MIx[i, j]!;
                var iyix = st.IyIx[i, j]!;
                mix.RemoveAll(k => (path[i - k, j] & TraceM) == 0);
                iyix.RemoveAll(k => (path[i - k, j] & TraceIy) == 0);
                if (mix.Count == 0 && iyix.Count == 0) path[i, j] &= ~TraceIx; else path[i, j] |= TraceIx;
                var miy = st.MIy[i, j]!;
                var ixiy = st.IxIy[i, j]!;
                miy.RemoveAll(k => (path[i, j - k] & TraceM) == 0);
                ixiy.RemoveAll(k => (path[i, j - k] & TraceIx) == 0);
                if (miy.Count == 0 && ixiy.Count == 0) path[i, j] &= ~TraceIy; else path[i, j] |= TraceIy;
            }
        }
    }

    // the gap states of one cell: the best score over the run lengths and, with traces, the lists of the lengths that tie
    private static void GapRuns(WsbState st, int i, int j, int nA, int nB, bool traces, Func<int, int, double> ins, Func<int, int, double> del, int queryPos, bool local, double eps, ref double maximum)
    {
        _ = nA;
        _ = nB;
        var M = st.M;
        var Ix = st.Ix;
        var Iy = st.Iy;
        // Ix: runs of k residues of the target (deletion)
        double score = local ? 0.0 : NegativeInfinity;
        List<int>? gapM = null, gapXY = null;
        if (traces)
        {
            gapM = [];
            gapXY = [];
            score = NegativeInfinity;
        }

        if (!traces)
        {
            for (int k = 1; k <= i; k++)
            {
                double g = del(queryPos, k);
                double t = M[i - k, j] + g;
                if (t > score) score = t;
                t = Iy[i - k, j] + g;
                if (t > score) score = t;
            }

            if (local && score > maximum) maximum = score;
            Ix[i, j] = score;
            score = local ? 0.0 : NegativeInfinity;
            for (int k = 1; k <= j; k++)
            {
                double g = ins(i, k);
                double t = M[i, j - k] + g;
                if (t > score) score = t;
                t = Ix[i, j - k] + g;
                if (t > score) score = t;
            }

            if (local && score > maximum) maximum = score;
            Iy[i, j] = score;
            return;
        }

        for (int k = 1; k <= i; k++)
        {
            double g = del(queryPos, k);
            TraceGap(M[i - k, j] + g, Iy[i - k, j] + g, k, ref score, gapM!, gapXY!, eps);
        }

        if (local)
        {
            if (score < eps)
            {
                score = NegativeInfinity;
                gapM!.Clear();
                gapXY!.Clear();
            }
            else if (score > maximum && st.GapsRaiseMaximum)
            {
                maximum = score;
            }
        }

        Ix[i, j] = score;
        st.MIx[i, j] = gapM;
        st.IyIx[i, j] = gapXY;
        gapM = [];
        gapXY = [];
        score = NegativeInfinity;
        for (int k = 1; k <= j; k++)
        {
            double g = ins(i, k);
            TraceGap(M[i, j - k] + g, Ix[i, j - k] + g, k, ref score, gapM, gapXY, eps);
        }

        if (local)
        {
            if (score < eps)
            {
                score = NegativeInfinity;
                gapM.Clear();
                gapXY.Clear();
            }
            else if (score > maximum && st.GapsRaiseMaximum)
            {
                maximum = score;
            }
        }

        Iy[i, j] = score;
        st.MIy[i, j] = gapM;
        st.IxIy[i, j] = gapXY;
    }

    // SELECT_TRACE_WATERMAN_SMITH_BEYER_GAP
    private static void TraceGap(double fromM, double fromOther, int k, ref double score, List<int> gapM, List<int> gapXY, double eps)
    {
        if (fromM > score - eps)
        {
            if (fromM > score + eps)
            {
                score = fromM;
                gapM.Clear();
                gapXY.Clear();
            }

            gapM.Add(k);
        }

        if (fromOther > score - eps)
        {
            if (fromOther > score + eps)
            {
                score = fromOther;
                gapM.Clear();
                gapXY.Clear();
            }

            gapXY.Add(k);
        }
    }

    private double WsbScore(string a, string b, char strand)
    {
        bool local = Mode == PairwiseAlignerMode.Local;
        var st = RunWsb(a, b, strand, local, traces: false);
        if (local)
            return st.Maximum;
        int nA = a.Length, nB = b.Length;
        return Math.Max(st.M[nA, nB], Math.Max(st.Ix[nA, nB], st.Iy[nA, nB]));
    }

    private (double Score, List<AlignmentResult> Results) WsbAlign(string a, string b, char strand)
    {
        bool local = Mode == PairwiseAlignerMode.Local;
        int nA = a.Length, nB = b.Length;
        double eps = Epsilon;
        var st = RunWsb(a, b, strand, local, traces: true);
        var results = new List<AlignmentResult>();
        double score;
        // the explicit path walk: the steps (state, end cell, start cell) from the end of the alignment to its start
        var segments = new List<(int State, int Ei, int Ej, int Si, int Sj)>();
        void Finish(int startI, int startJ, int endI, int endJ)
        {
            var ra = new System.Text.StringBuilder();
            var rb = new System.Text.StringBuilder();
            for (int n = segments.Count - 1; n >= 0; n--)
            {
                var (state, ei, ej, si, sj) = segments[n];
                switch (state)
                {
                    case TraceM:
                        ra.Append(a[si]);
                        rb.Append(b[sj]);
                        break;
                    case TraceIx:
                        ra.Append(a, si, ei - si);
                        rb.Append('-', ei - si);
                        break;
                    default:
                        ra.Append('-', ej - sj);
                        rb.Append(b, sj, ej - sj);
                        break;
                }
            }

            results.Add(new AlignmentResult(ra.ToString(), rb.ToString(), 0, local ? AlignmentType.Local : AlignmentType.Global, startI, startJ, endI - 1, endJ - 1));
        }

        bool Walk(int state, int i, int j, int endI, int endJ)
        {
            // follows the traces from the cell (i, j) in the state back to the origin; false when the limit of alignments is reached
            if (results.Count >= MaxAlignments)
                return false;
            switch (state)
            {
                case TraceM:
                    {
                        int trace = st.Trace[i, j];
                        if (local && (trace & TraceStart) != 0 && (trace & (TraceM | TraceIx | TraceIy)) == 0)
                        {
                            // the alignment starts at the cell (i, j): no pair of this cell
                            Finish(i, j, endI, endJ);
                            return true;
                        }

                        if (!local && i == 0 && j == 0)
                        {
                            Finish(0, 0, endI, endJ);
                            return true;
                        }

                        if (i == 0 || j == 0)
                            return true;
                        segments.Add((TraceM, i, j, i - 1, j - 1));
                        bool ok = true;
                        foreach (int prev in new[] { TraceM, TraceIx, TraceIy })
                        {
                            if ((trace & prev) == 0)
                                continue;
                            if (!Walk(prev, i - 1, j - 1, endI, endJ))
                            {
                                ok = false;
                                break;
                            }
                        }

                        segments.RemoveAt(segments.Count - 1);
                        return ok;
                    }

                case TraceIx:
                    {
                        if (!local && j == 0)
                        {
                            // the first column: one run from the origin
                            segments.Add((TraceIx, i, 0, 0, 0));
                            bool done = Walk(TraceM, 0, 0, endI, endJ);
                            segments.RemoveAt(segments.Count - 1);
                            return done;
                        }

                        var fromM = st.MIx[i, j];
                        var fromIy = st.IyIx[i, j];
                        foreach (var (list, next) in new[] { (fromM, TraceM), (fromIy, TraceIy) })
                        {
                            if (list is null)
                                continue;
                            foreach (int k in list)
                            {
                                segments.Add((TraceIx, i, j, i - k, j));
                                bool ok = Walk(next, i - k, j, endI, endJ);
                                segments.RemoveAt(segments.Count - 1);
                                if (!ok)
                                    return false;
                            }
                        }

                        return true;
                    }

                default:
                    {
                        if (!local && i == 0)
                        {
                            segments.Add((TraceIy, 0, j, 0, 0));
                            bool done = Walk(TraceM, 0, 0, endI, endJ);
                            segments.RemoveAt(segments.Count - 1);
                            return done;
                        }

                        var fromM = st.MIy[i, j];
                        var fromIx = st.IxIy[i, j];
                        foreach (var (list, next) in new[] { (fromM, TraceM), (fromIx, TraceIx) })
                        {
                            if (list is null)
                                continue;
                            foreach (int k in list)
                            {
                                segments.Add((TraceIy, i, j, i, j - k));
                                bool ok = Walk(next, i, j - k, endI, endJ);
                                segments.RemoveAt(segments.Count - 1);
                                if (!ok)
                                    return false;
                            }
                        }

                        return true;
                    }
            }
        }

        if (local)
        {
            score = st.Maximum;
            for (int i = 0; i <= nA && results.Count < MaxAlignments; i++)
            {
                for (int j = 0; j <= nB; j++)
                {
                    if ((st.Trace[i, j] & TraceEnd) == 0)
                        continue;
                    Walk(TraceM, i, j, i, j);
                }
            }
        }
        else
        {
            double m = st.M[nA, nB], x = st.Ix[nA, nB], y = st.Iy[nA, nB];
            score = Math.Max(m, Math.Max(x, y));
            if (m >= score - eps)
                Walk(TraceM, nA, nB, nA, nB);
            if (x >= score - eps)
                Walk(TraceIx, nA, nB, nA, nB);
            if (y >= score - eps)
                Walk(TraceIy, nA, nB, nA, nB);
        }

        return (score, [.. results.Select(r => r with { Score = (int)Math.Round(score, MidpointRounding.AwayFromZero), ExactScore = score })]);
    }
}
