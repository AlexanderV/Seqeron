using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Alignment;

// ksw2-style two-piece (dual) affine gap alignment (audit A2-29, B09 F123): Li's ksw2 (minimap2) scores a gap of length k as the better of two affine functions, min(q + k e, q2 + k e2) as a cost, i.e. a short-gap and a
// long-gap regime (the long gap opens dearer but extends cheaper); the dynamic program is Gotoh's with two gap states per direction (H, E1, F1, E2, F2) in O(m n), optionally restricted to a diagonal band.

/// <summary>
/// Scoring of the dual affine gap model: <paramref name="Match"/> / <paramref name="Mismatch"/> for an identical / different pair and two affine gap functions. A gap of length k scores
/// <c>max(GapOpen1 + k * GapExtend1, GapOpen2 + k * GapExtend2)</c> (all four gap numbers are scores, at most 0: the Seqeron convention of <see cref="ScoringMatrix"/>; ksw2's <c>q + k e</c> cost is <c>GapOpen = -q</c>, <c>GapExtend = -e</c>).
/// </summary>
public sealed record DualAffineScoring(int Match, int Mismatch, int GapOpen1, int GapExtend1, int GapOpen2, int GapExtend2)
{
    /// <summary>The score of a gap of length <paramref name="k"/> (k at least 1).</summary>
    public long GapScore(int k) => Math.Max(GapOpen1 + ((long)k * GapExtend1), GapOpen2 + ((long)k * GapExtend2));

    /// <summary>minimap2's default (<c>-A 2 -B 4 -O 4,24 -E 2,1</c>): match 2, mismatch -4, gaps min(4 + 2k, 24 + k).</summary>
    public static DualAffineScoring Minimap2 { get; } = new(2, -4, -4, -2, -24, -1);
}

public static partial class SequenceAligner
{
    private const int DualNegative = int.MinValue / 4;

    /// <summary>
    /// Global alignment (end to end) with the dual affine gap model of ksw2 (<c>ksw2_extd2</c> in global mode), O(m n) time and memory; with <paramref name="band"/> at least 0 only the cells
    /// with |i - j| at most <paramref name="band"/> are computed (ksw2 <c>w</c>; the end cell must lie inside the band: |m - n| at most <paramref name="band"/>).
    /// </summary>
    /// <param name="sequence1">The first sequence (rows of the matrix, "target").</param>
    /// <param name="sequence2">The second sequence ("query").</param>
    /// <param name="scoring">Match, mismatch and the two gap functions.</param>
    /// <param name="band">Band half-width, or -1 for the full matrix.</param>
    /// <returns>The alignment; <see cref="AlignmentResult.EndPosition1"/> / <see cref="AlignmentResult.EndPosition2"/> are 0-based INCLUSIVE like those of every other aligner (length − 1 for a global result).</returns>
    public static AlignmentResult GlobalAlignDualAffine(string sequence1, string sequence2, DualAffineScoring scoring, int band = -1) =>
        DualAffineCore(sequence1, sequence2, scoring, band, local: false);

    /// <summary>Local alignment (Smith-Waterman, zero floor) with the dual affine gap model; see <see cref="GlobalAlignDualAffine"/> for the band.</summary>
    public static AlignmentResult LocalAlignDualAffine(string sequence1, string sequence2, DualAffineScoring scoring, int band = -1) =>
        DualAffineCore(sequence1, sequence2, scoring, band, local: true);

    private static AlignmentResult DualAffineCore(string a, string b, DualAffineScoring sc, int band, bool local)
    {
        ArgumentNullException.ThrowIfNull(a);
        ArgumentNullException.ThrowIfNull(b);
        ArgumentNullException.ThrowIfNull(sc);
        if (sc.GapOpen1 > 0 || sc.GapExtend1 > 0 || sc.GapOpen2 > 0 || sc.GapExtend2 > 0)
            throw new ArgumentException("The gap scores must be at most 0 (penalties).", nameof(sc));
        if (band < -1)
            throw new ArgumentOutOfRangeException(nameof(band), "band must be -1 (full) or at least 0.");
        a = a.ToUpperInvariant();
        b = b.ToUpperInvariant();
        int m = a.Length;
        int n = b.Length;
        if ((long)(m + 1) * (n + 1) > 25_000_000)
            throw new ArgumentException("The sequences are too long for the O(m n) dual affine program (limit 25 million cells).");
        if (!local && band >= 0 && Math.Abs(m - n) > band)
            throw new ArgumentException("The end cell lies outside the band (|m - n| > band): the global alignment does not exist.", nameof(band));
        long limit = (Math.Max(Math.Abs((long)sc.Match), Math.Abs((long)sc.Mismatch)) + Math.Abs((long)sc.GapOpen1) + Math.Abs((long)sc.GapExtend1) + Math.Abs((long)sc.GapOpen2) + Math.Abs((long)sc.GapExtend2) + 1) * ((long)m + n + 2);
        if (limit >= 1 << 28)
            throw new ArgumentException("Score range too large for the integer dynamic program.", nameof(sc));

        int w = n + 1;
        var h = new int[(m + 1) * w];
        var e1 = new int[(m + 1) * w];      // gap in sequence1 (consumes sequence2): horizontal move
        var e2 = new int[(m + 1) * w];
        var f1 = new int[(m + 1) * w];      // gap in sequence2 (consumes sequence1): vertical move
        var f2 = new int[(m + 1) * w];
        Array.Fill(h, DualNegative);
        Array.Fill(e1, DualNegative);
        Array.Fill(e2, DualNegative);
        Array.Fill(f1, DualNegative);
        Array.Fill(f2, DualNegative);
        bool Inside(int i, int j) => band < 0 || Math.Abs(i - j) <= band;
        h[0] = 0;
        int best = 0, bi = 0, bj = 0;
        for (int i = 0; i <= m; i++)
        {
            AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per row
            for (int j = 0; j <= n; j++)
            {
                if ((i == 0 && j == 0) || !Inside(i, j))
                    continue;
                int c = (i * w) + j;
                int ee1 = DualNegative, ee2 = DualNegative, ff1 = DualNegative, ff2 = DualNegative, hh = DualNegative;
                if (j > 0)
                {
                    int left = c - 1;
                    ee1 = Math.Max(Add(h[left], sc.GapOpen1 + sc.GapExtend1), Add(e1[left], sc.GapExtend1));
                    ee2 = Math.Max(Add(h[left], sc.GapOpen2 + sc.GapExtend2), Add(e2[left], sc.GapExtend2));
                }

                if (i > 0)
                {
                    int up = c - w;
                    ff1 = Math.Max(Add(h[up], sc.GapOpen1 + sc.GapExtend1), Add(f1[up], sc.GapExtend1));
                    ff2 = Math.Max(Add(h[up], sc.GapOpen2 + sc.GapExtend2), Add(f2[up], sc.GapExtend2));
                }

                if (i > 0 && j > 0)
                    hh = Add(h[c - w - 1], a[i - 1] == b[j - 1] ? sc.Match : sc.Mismatch);
                hh = Math.Max(hh, Math.Max(Math.Max(ee1, ee2), Math.Max(ff1, ff2)));
                if (local && hh < 0)
                    hh = 0;
                h[c] = hh;
                e1[c] = ee1;
                e2[c] = ee2;
                f1[c] = ff1;
                f2[c] = ff2;
                if (local && hh > best)
                {
                    best = hh;
                    bi = i;
                    bj = j;
                }
            }
        }

        int score;
        int ei, ej;
        if (local)
        {
            score = best;
            ei = bi;
            ej = bj;
        }
        else
        {
            score = h[(m * w) + n];
            ei = m;
            ej = n;
        }

        // traceback over the five states
        var r1 = new StringBuilder();
        var r2 = new StringBuilder();
        int ti = ei, tj = ej;
        int state = 0;      // 0 H, 1 E1, 2 E2, 3 F1, 4 F2
        while (ti > 0 || tj > 0)
        {
            int c = (ti * w) + tj;
            if (state == 0)
            {
                int cur = h[c];
                if (local && cur == 0)
                    break;
                if (ti > 0 && tj > 0 && cur == Add(h[c - w - 1], a[ti - 1] == b[tj - 1] ? sc.Match : sc.Mismatch))
                {
                    r1.Append(a[ti - 1]);
                    r2.Append(b[tj - 1]);
                    ti--;
                    tj--;
                }
                else if (cur == e1[c]) state = 1;
                else if (cur == e2[c]) state = 2;
                else if (cur == f1[c]) state = 3;
                else if (cur == f2[c]) state = 4;
                else throw new InvalidOperationException("Dual affine traceback lost the path.");
            }
            else
            {
                bool horizontal = state is 1 or 2;
                int[] gap = state switch { 1 => e1, 2 => e2, 3 => f1, _ => f2 };
                int open = state is 1 or 3 ? sc.GapOpen1 + sc.GapExtend1 : sc.GapOpen2 + sc.GapExtend2;
                int ext = state is 1 or 3 ? sc.GapExtend1 : sc.GapExtend2;
                int prev = horizontal ? c - 1 : c - w;
                if (horizontal)
                {
                    r1.Append('-');
                    r2.Append(b[tj - 1]);
                    tj--;
                }
                else
                {
                    r1.Append(a[ti - 1]);
                    r2.Append('-');
                    ti--;
                }

                if (gap[c] == Add(h[prev], open))
                    state = 0;
                else if (gap[c] != Add(gap[prev], ext))
                    throw new InvalidOperationException("Dual affine traceback lost the gap.");
            }
        }

        var x = Reverse(r1);
        var y = Reverse(r2);
        if (local && x.Length == 0)
            return new AlignmentResult("", "", 0, AlignmentType.Local, -1, -1, -1, -1);     // no positive-scoring region: the convention of LocalAlign
        return new AlignmentResult(x, y, score, local ? AlignmentType.Local : AlignmentType.Global, ti, tj, ei - 1, ej - 1);       // 0-based INCLUSIVE ends (R14-1)
    }

    private static int Add(int value, int delta) => value <= DualNegative / 2 ? DualNegative : value + delta;

    private static string Reverse(StringBuilder sb)
    {
        var a = new char[sb.Length];
        for (int i = 0; i < a.Length; i++)
            a[i] = sb[a.Length - 1 - i];
        return new string(a);
    }
}
