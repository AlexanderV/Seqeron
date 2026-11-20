namespace Seqeron.Genomics.Alignment;

/// <summary>
/// Which of the co-optimal local alignments a top-N run reports when several have the same score (audit A2-28d). Every
/// order reports a legal Waterman-Eggert sequence (each score is the true optimum of the matrix with the aligned pairs of the
/// earlier reports forbidden); only the choice among equal scores differs.
/// </summary>
public enum TopNTieOrder
{
    /// <summary>The fixed order of <c>LocalAlign</c> (first maximal cell in row-major order, traceback diagonal &gt; up &gt; left, stop at the first zero cell), applied round after round.</summary>
    WatermanEggert = 0,

    /// <summary>
    /// The order of EMBOSS <c>matcher</c> = Huang-Miller SIM: the K best alignments are collected during one forward pass by start point
    /// (<c>matcher_BigPass</c>, the K-vertex LIST with its replacement rule), reported in LIST order (<c>matcher_Findmax</c>), recovered by the
    /// linear-space divide and conquer with the program's own tie rules (<c>matcher_Diff</c>) and re-located after each report
    /// (<c>matcher_Locate</c> / <c>matcher_SmallPass</c>). Needs integer scores; the list depends on the requested number of alignments (the program's LIST size).
    /// </summary>
    Sim = 1,
}

public static partial class SequenceAligner
{
    #region EMBOSS matcher (Huang-Miller SIM) tie order

    /// <summary>
    /// <see cref="TopNTieOrder.Sim"/> run: a faithful port of EMBOSS 6.6.0 <c>matcher.c</c> (<c>matcher_Sim</c>, <c>matcher_BigPass</c>, <c>matcher_Locate</c>,
    /// <c>matcher_SmallPass</c>, <c>matcher_Addnode</c>, <c>matcher_NoCross</c>, <c>matcher_Findmax</c>, <c>matcher_Diff</c>, <c>matcher_Calcons</c>).
    /// The program aborts ("The number of alignments computed is too large") when fewer than K local alignments exist; this port then returns the alignments found so far.
    /// </summary>
    internal static IReadOnlyList<AlignmentResult> LocalTopNSimCore(string a, string b, PairScoring scoring, bool affine, int maxAlignments, long minScaled)
    {
        if (maxAlignments < 1)
            throw new ArgumentOutOfRangeException(nameof(maxAlignments), maxAlignments, "maxAlignments must be at least 1.");
        var results = new List<AlignmentResult>();
        if (a.Length == 0 || b.Length == 0)
            return results;
        if (scoring.Gaps is not null)
            throw new ArgumentException("TopNTieOrder.Sim has no position-specific gap scores.", nameof(scoring));
        scoring.EnsureRange(a.Length, b.Length);
        var sim = new SimMatcher(a, b, scoring, affine);
        sim.Run(maxAlignments, minScaled, results);
        return results;
    }

    private sealed class SimMatcher
    {
        private sealed class Vertex
        {
            public int Score, StarI, StarJ, EndI, EndJ, Top, Bot, Left, Right;
        }

        private readonly string _s1, _s2;            // original (upper-case) sequences
        private readonly char[] _a, _b;               // 1-based with a leading blank
        private readonly PairScoring _sc;
        private readonly int _m, _n;
        private readonly int _q, _r, _qr;

        private Vertex[] _list = Array.Empty<Vertex>();
        private Vertex? _low, _most;
        private int _lmin;
        private int _flag;
        private int _m1, _mm, _n1, _nn, _rl, _cl;

        private readonly int[] _cc, _dd, _rr, _ss, _ee, _ff;     // row vectors (length n + 1)
        private readonly int[] _hh, _ww, _ii, _jj, _xx, _yy;     // column vectors (length m + 1)
        private readonly HashSet<int>?[] _row;                    // used aligned pairs: row -> columns

        private int[] _script = Array.Empty<int>();
        private int _sp, _last, _curI, _curJ;

        public SimMatcher(string a, string b, PairScoring sc, bool affine)
        {
            _s1 = a;
            _s2 = b;
            _sc = sc;
            _m = a.Length;
            _n = b.Length;
            _a = (" " + a).ToCharArray();
            _b = (" " + b).ToCharArray();
            // matcher: -gapopen G -gapextend E, q = G - E, r = E; a gap of length k scores -(q + r k). Linear model: q = 0.
            _q = affine ? -sc.GapOpen : 0;
            _r = -sc.GapExtend;
            _qr = _q + _r;
            _cc = new int[_n + 1]; _dd = new int[_n + 1]; _rr = new int[_n + 1]; _ss = new int[_n + 1]; _ee = new int[_n + 1]; _ff = new int[_n + 1];
            _hh = new int[_m + 1]; _ww = new int[_m + 1]; _ii = new int[_m + 1]; _jj = new int[_m + 1]; _xx = new int[_m + 1]; _yy = new int[_m + 1];
            _row = new HashSet<int>?[_m + 1];
        }

        private int Sub(int i, int j) => _sc.Sub(_a[i], _b[j]);

        private bool Used(int i, int j) => _row[i] is { } s && s.Contains(j);

        private static void Order(ref int s1, ref int x1, ref int y1, int s2, int x2, int y2)
        {
            if (s1 < s2) { s1 = s2; x1 = x2; y1 = y2; }
            else if (s1 == s2)
            {
                if (x1 < x2) { x1 = x2; y1 = y2; }
                else if (x1 == x2 && y1 < y2) y1 = y2;
            }
        }

        private int Gap(int k) => k <= 0 ? 0 : _q + _r * k;

        public void Run(int k, long minScaled, List<AlignmentResult> results)
        {
            _script = new int[_m + _n + 2];
            _list = new Vertex[k];
            for (int i = 0; i < k; i++) _list[i] = new Vertex();
            int numnode = 0;
            _lmin = 0;
            BigPass(k, ref numnode);

            for (int count = k - 1; count >= 0; count--)
            {
                AlignmentCancellation.ThrowIfRequested();
                if (numnode == 0)
                    return; // the program aborts here ("The number of alignments computed is too large")
                var cur = Findmax(ref numnode);
                int score = cur.Score;
                int stari = ++cur.StarI;
                int starj = ++cur.StarJ;
                int endi = cur.EndI;
                int endj = cur.EndJ;
                _m1 = cur.Top;
                _mm = cur.Bot;
                _n1 = cur.Left;
                _nn = cur.Right;
                _rl = endi - stari + 1;
                _cl = endj - starj + 1;
                _curI = stari - 1;
                _curJ = starj - 1;
                _sp = 0;
                _last = 0;
                Diff(stari - 1, starj - 1, _rl, _cl, _q, _q);

                int min0 = stari - 1, min1 = starj - 1, max0 = stari + _rl - 1, max1 = starj + _cl - 1;
                if (score < minScaled)
                    return;
                results.Add(Calcons(min0, min1, max0, max1, score));

                if (count > 0)
                {
                    _flag = 0;
                    Locate(numnode);
                    if (_flag != 0)
                        SmallPass(count, ref numnode);
                }
            }
        }

        private AlignmentResult Calcons(int min0, int min1, int max0, int max1, int score)
        {
            var c0 = new System.Text.StringBuilder();
            var c1 = new System.Text.StringBuilder();
            int rp = 0, op = 0, i0 = min0, i1 = min1;
            while (i0 < max0 || i1 < max1)
            {
                if (op == 0 && _script[rp] == 0)
                {
                    op = _script[rp++];
                    c0.Append(_s1[i0++]);
                    c1.Append(_s2[i1++]);
                }
                else
                {
                    if (op == 0)
                        op = _script[rp++];
                    if (op > 0)
                    {
                        c0.Append('-');
                        c1.Append(_s2[i1++]);
                        op--;
                    }
                    else
                    {
                        c0.Append(_s1[i0++]);
                        c1.Append('-');
                        op++;
                    }
                }
            }
            return new AlignmentResult(c0.ToString(), c1.ToString(), score, AlignmentType.Local, min0, min1, max0 - 1, max1 - 1);
        }

        // ------------------------------------------------------------------------------------------------ LIST handling

        private int Addnode(int c, int ci, int cj, int i, int j, int K, int cost, ref int numnode)
        {
            bool found = false;
            if (_most is not null && _most.StarI == ci && _most.StarJ == cj)
                found = true;
            else
                for (int d = 0; d < numnode; d++)
                {
                    _most = _list[d];
                    if (_most.StarI == ci && _most.StarJ == cj)
                    {
                        found = true;
                        break;
                    }
                }

            if (found)
            {
                if (_most!.Score < c)
                {
                    _most.Score = c;
                    _most.EndI = i;
                    _most.EndJ = j;
                }
                if (_most.Top > i) _most.Top = i;
                if (_most.Bot < i) _most.Bot = i;
                if (_most.Left > j) _most.Left = j;
                if (_most.Right < j) _most.Right = j;
            }
            else
            {
                _most = numnode == K ? _low : _list[numnode++];
                _most!.Score = c;
                _most.StarI = ci;
                _most.StarJ = cj;
                _most.EndI = i;
                _most.EndJ = j;
                _most.Top = _most.Bot = i;
                _most.Left = _most.Right = j;
            }

            if (numnode == K)
            {
                if (_low == _most || _low is null)
                {
                    _low = _list[0];
                    for (int d = 1; d < numnode; d++)
                        if (_list[d].Score < _low.Score)
                            _low = _list[d];
                }
                return _low.Score;
            }
            return cost;
        }

        private Vertex Findmax(ref int numnode)
        {
            int j = 0;
            for (int i = 1; i < numnode; i++)
                if (_list[i].Score > _list[j].Score)
                    j = i;
            var cur = _list[j];
            if (j != --numnode)
            {
                _list[j] = _list[numnode];
                _list[numnode] = cur;
            }
            _most = _list[0];
            if (_low == cur)
                _low = _list[0];
            return cur;
        }

        private bool NoCross(int numnode)
        {
            int i;
            for (i = 0; i < numnode; i++)
            {
                var cur = _list[i];
                if (cur.StarI <= _mm && cur.StarJ <= _nn && cur.Bot >= _m1 - 1 && cur.Right >= _n1 - 1 && (cur.StarI < _rl || cur.StarJ < _cl))
                {
                    if (cur.StarI < _rl) _rl = cur.StarI;
                    if (cur.StarJ < _cl) _cl = cur.StarJ;
                    _flag = 1;
                    break;
                }
            }
            return i == numnode;
        }

        // ------------------------------------------------------------------------------------------------ passes (nseq == 2)

        private void BigPass(int K, ref int numnode)
        {
            int M = _m, N = _n;
            for (int j = 1; j <= N; j++)
            {
                _cc[j] = 0; _rr[j] = 0; _ee[j] = j;
                _dd[j] = -_q; _ss[j] = 0; _ff[j] = j;
            }
            for (int i = 1; i <= M; i++)
            {
                AlignmentCancellation.ThrowIfRequested();      // A4-1: one poll per row of the big pass
                int c = 0;
                int f = -_q;
                int ci = i, fi = i;
                int p = 0;
                int pi = i - 1;
                int cj = 0, fj = 0, pj = 0;
                for (int j = 1; j <= N; j++)
                {
                    f -= _r;
                    c -= _qr;
                    Order(ref f, ref fi, ref fj, c, ci, cj);
                    c = _cc[j] - _qr;
                    ci = _rr[j];
                    cj = _ee[j];
                    int d = _dd[j] - _r;
                    int di = _ss[j];
                    int dj = _ff[j];
                    Order(ref d, ref di, ref dj, c, ci, cj);
                    c = 0;
                    if (!Used(i, j)) c = p + Sub(i, j);
                    if (c <= 0) { c = 0; ci = i; cj = j; }
                    else { ci = pi; cj = pj; }
                    Order(ref c, ref ci, ref cj, d, di, dj);
                    Order(ref c, ref ci, ref cj, f, fi, fj);
                    p = _cc[j];
                    _cc[j] = c;
                    pi = _rr[j];
                    pj = _ee[j];
                    _rr[j] = ci;
                    _ee[j] = cj;
                    _dd[j] = d;
                    _ss[j] = di;
                    _ff[j] = dj;
                    if (c > _lmin)
                        _lmin = Addnode(c, ci, cj, i, j, K, _lmin, ref numnode);
                }
            }
        }

        private void Locate(int numnode)
        {
            int i, j, c, p, f, d, ci = 0, cj = 0, di = 0, dj = 0, fi = 0, fj = 0, pi, pj;
            bool cflag, rflag;
            // Reverse pass
            for (j = _nn; j >= _n1; j--)
            {
                _cc[j] = 0;
                _ee[j] = j;
                _dd[j] = -_q;
                _ff[j] = j;
                _rr[j] = _ss[j] = _mm + 1;
            }

            for (i = _mm; i >= _m1; i--)
            {
                c = p = 0;
                f = -_q;
                ci = fi = i;
                pi = i + 1;
                cj = fj = pj = _nn + 1;
                int limit = _n1;
                for (j = _nn; j >= limit; j--)
                {
                    f -= _r;
                    c -= _qr;
                    Order(ref f, ref fi, ref fj, c, ci, cj);
                    c = _cc[j] - _qr;
                    ci = _rr[j];
                    cj = _ee[j];
                    d = _dd[j] - _r;
                    di = _ss[j];
                    dj = _ff[j];
                    Order(ref d, ref di, ref dj, c, ci, cj);
                    c = 0;
                    if (!Used(i, j)) c = p + Sub(i, j);
                    if (c <= 0) { c = 0; ci = i; cj = j; }
                    else { ci = pi; cj = pj; }
                    Order(ref c, ref ci, ref cj, d, di, dj);
                    Order(ref c, ref ci, ref cj, f, fi, fj);
                    p = _cc[j];
                    _cc[j] = c;
                    pi = _rr[j];
                    pj = _ee[j];
                    _rr[j] = ci;
                    _ee[j] = cj;
                    _dd[j] = d;
                    _ss[j] = di;
                    _ff[j] = dj;
                    if (c > _lmin) _flag = 1;
                }
                _hh[i] = _cc[_n1];
                _ii[i] = _rr[_n1];
                _jj[i] = _ee[_n1];
                _ww[i] = _dd[_n1];
                _xx[i] = _ss[_n1];
                _yy[i] = _ff[_n1];
            }

            _rl = _m1;
            _cl = _n1;
            while (true)
            {
                rflag = cflag = true;
                while ((rflag && _m1 > 1) || (cflag && _n1 > 1))
                {
                    if (rflag && _m1 > 1) // compute one row
                    {
                        rflag = false;
                        _m1--;
                        c = p = 0;
                        f = -_q;
                        ci = fi = _m1;
                        pi = _m1 + 1;
                        cj = fj = pj = _nn + 1;
                        for (j = _nn; j >= _n1; j--)
                        {
                            f -= _r;
                            c -= _qr;
                            Order(ref f, ref fi, ref fj, c, ci, cj);
                            c = _cc[j] - _qr;
                            ci = _rr[j];
                            cj = _ee[j];
                            d = _dd[j] - _r;
                            di = _ss[j];
                            dj = _ff[j];
                            Order(ref d, ref di, ref dj, c, ci, cj);
                            c = 0;
                            if (!Used(_m1, j)) c = p + Sub(_m1, j);
                            if (c <= 0) { c = 0; ci = _m1; cj = j; }
                            else { ci = pi; cj = pj; }
                            Order(ref c, ref ci, ref cj, d, di, dj);
                            Order(ref c, ref ci, ref cj, f, fi, fj);
                            p = _cc[j];
                            _cc[j] = c;
                            pi = _rr[j];
                            pj = _ee[j];
                            _rr[j] = ci;
                            _ee[j] = cj;
                            _dd[j] = d;
                            _ss[j] = di;
                            _ff[j] = dj;
                            if (c > _lmin) _flag = 1;
                            if (!rflag && ((ci > _rl && cj > _cl) || (di > _rl && dj > _cl) || (fi > _rl && fj > _cl)))
                                rflag = true;
                        }
                        _hh[_m1] = _cc[_n1];
                        _ii[_m1] = _rr[_n1];
                        _jj[_m1] = _ee[_n1];
                        _ww[_m1] = _dd[_n1];
                        _xx[_m1] = _ss[_n1];
                        _yy[_m1] = _ff[_n1];
                        if (!cflag && ((ci > _rl && cj > _cl) || (di > _rl && dj > _cl) || (fi > _rl && fj > _cl)))
                            cflag = true;
                    }

                    if (cflag && _n1 > 1) // compute one column
                    {
                        cflag = false;
                        _n1--;
                        c = 0;
                        f = -_q;
                        cj = fj = _n1;
                        p = 0;
                        ci = fi = pi = _mm + 1;
                        pj = _n1 + 1;
                        int limit = _mm;
                        for (i = limit; i >= _m1; i--)
                        {
                            f -= _r;
                            c -= _qr;
                            Order(ref f, ref fi, ref fj, c, ci, cj);
                            c = _hh[i] - _qr;
                            ci = _ii[i];
                            cj = _jj[i];
                            d = _ww[i] - _r;
                            di = _xx[i];
                            dj = _yy[i];
                            Order(ref d, ref di, ref dj, c, ci, cj);
                            c = 0;
                            if (!Used(i, _n1)) c = p + Sub(i, _n1);
                            if (c <= 0) { c = 0; ci = i; cj = _n1; }
                            else { ci = pi; cj = pj; }
                            Order(ref c, ref ci, ref cj, d, di, dj);
                            Order(ref c, ref ci, ref cj, f, fi, fj);
                            p = _hh[i];
                            _hh[i] = c;
                            pi = _ii[i];
                            pj = _jj[i];
                            _ii[i] = ci;
                            _jj[i] = cj;
                            _ww[i] = d;
                            _xx[i] = di;
                            _yy[i] = dj;
                            if (c > _lmin) _flag = 1;
                            if (!cflag && ((ci > _rl && cj > _cl) || (di > _rl && dj > _cl) || (fi > _rl && fj > _cl)))
                                cflag = true;
                        }
                        _cc[_n1] = _hh[_m1];
                        _rr[_n1] = _ii[_m1];
                        _ee[_n1] = _jj[_m1];
                        _dd[_n1] = _ww[_m1];
                        _ss[_n1] = _xx[_m1];
                        _ff[_n1] = _yy[_m1];
                        if (!rflag && ((ci > _rl && cj > _cl) || (di > _rl && dj > _cl) || (fi > _rl && fj > _cl)))
                            rflag = true;
                    }
                }
                if ((_m1 == 1 && _n1 == 1) || NoCross(numnode))
                    break;
            }
            _m1--;
            _n1--;
        }

        private void SmallPass(int count, ref int numnode)
        {
            for (int j = _n1 + 1; j <= _nn; j++)
            {
                _cc[j] = 0;
                _rr[j] = _m1;
                _ee[j] = j;
                _dd[j] = -_q;
                _ss[j] = _m1;
                _ff[j] = j;
            }
            for (int i = _m1 + 1; i <= _mm; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                int c = 0;
                int f = -_q;
                int ci = i, fi = i;
                int p = 0;
                int pi = i - 1;
                int cj = _n1, fj = _n1, pj = _n1;
                int limit = _n1 + 1;
                for (int j = limit; j <= _nn; j++)
                {
                    f -= _r;
                    c -= _qr;
                    Order(ref f, ref fi, ref fj, c, ci, cj);
                    c = _cc[j] - _qr;
                    ci = _rr[j];
                    cj = _ee[j];
                    int d = _dd[j] - _r;
                    int di = _ss[j];
                    int dj = _ff[j];
                    Order(ref d, ref di, ref dj, c, ci, cj);
                    c = 0;
                    if (!Used(i, j)) c = p + Sub(i, j);
                    if (c <= 0) { c = 0; ci = i; cj = j; }
                    else { ci = pi; cj = pj; }
                    Order(ref c, ref ci, ref cj, d, di, dj);
                    Order(ref c, ref ci, ref cj, f, fi, fj);
                    p = _cc[j];
                    _cc[j] = c;
                    pi = _rr[j];
                    pj = _ee[j];
                    _rr[j] = ci;
                    _ee[j] = cj;
                    _dd[j] = d;
                    _ss[j] = di;
                    _ff[j] = dj;
                    if (c > _lmin)
                        _lmin = Addnode(c, ci, cj, i, j, count, _lmin, ref numnode);
                }
            }
        }

        // ------------------------------------------------------------------------------------------------ edit script (linear space)

        private void Del(int k)
        {
            _curI += k;
            if (_last < 0) _last = _script[_sp - 1] -= k;
            else _last = _script[_sp++] = -k;
        }

        private void Ins(int k)
        {
            _curJ += k;
            if (_last < 0) { _script[_sp - 1] = k; _script[_sp++] = _last; }
            else _last = _script[_sp++] = k;
        }

        private void Rep() => _last = _script[_sp++] = 0;

        /// <summary><c>matcher_Diff</c>: A and B are read at offsets (A[i] = _a[aOff + i]).</summary>
        private void Diff(int aOff, int bOff, int M, int N, int tb, int te)
        {
            if (N <= 0)
            {
                if (M > 0) Del(M);
                return;
            }

            if (M <= 1)
            {
                if (M <= 0)
                {
                    Ins(N);
                    return;
                }
                if (tb > te) tb = te;
                int midc1 = -(tb + _r + Gap(N));
                int midj1 = 0;
                for (int j = 1; j <= N; j++)
                {
                    if (!Used(_curI + 1, j + _curJ))
                    {
                        int c1 = Sub1(aOff + 1, bOff + j) - (Gap(j - 1) + Gap(N - j));
                        if (c1 > midc1)
                        {
                            midc1 = c1;
                            midj1 = j;
                        }
                    }
                }

                if (midj1 == 0)
                {
                    Ins(N);
                    Del(1);
                }
                else
                {
                    if (midj1 > 1) Ins(midj1 - 1);
                    Rep();
                    _curI++; _curJ++;
                    (_row[_curI] ??= new HashSet<int>()).Add(_curJ);
                    if (midj1 < N) Ins(N - midj1);
                }
                return;
            }

            int midi = M / 2;
            _cc[0] = 0;
            int t = -_q;
            for (int j = 1; j <= N; j++)
            {
                _cc[j] = t = t - _r;
                _dd[j] = t - _q;
            }
            t = -tb;
            int c, d, e, s;
            for (int i = 1; i <= midi; i++)
            {
                AlignmentCancellation.ThrowIfRequested();
                s = _cc[0];
                _cc[0] = c = t = t - _r;
                e = t - _q;
                for (int j = 1; j <= N; j++)
                {
                    if ((c = c - _qr) > (e = e - _r)) e = c;
                    if ((c = _cc[j] - _qr) > (d = _dd[j] - _r)) d = c;
                    if (!Used(i + _curI, j + _curJ)) c = s + Sub1(aOff + i, bOff + j);
                    if (c < d) c = d;
                    if (c < e) c = e;
                    s = _cc[j];
                    _cc[j] = c;
                    _dd[j] = d;
                }
            }
            _dd[0] = _cc[0];

            _rr[N] = 0;
            t = -_q;
            for (int j = N - 1; j >= 0; j--)
            {
                _rr[j] = t = t - _r;
                _ss[j] = t - _q;
            }
            t = -te;
            for (int i = M - 1; i >= midi; i--)
            {
                AlignmentCancellation.ThrowIfRequested();
                s = _rr[N];
                _rr[N] = c = t = t - _r;
                e = t - _q;
                for (int j = N - 1; j >= 0; j--)
                {
                    if ((c = c - _qr) > (e = e - _r)) e = c;
                    if ((c = _rr[j] - _qr) > (d = _ss[j] - _r)) d = c;
                    if (!Used(i + 1 + _curI, j + 1 + _curJ)) c = s + Sub1(aOff + i + 1, bOff + j + 1);
                    if (c < d) c = d;
                    if (c < e) c = e;
                    s = _rr[j];
                    _rr[j] = c;
                    _ss[j] = d;
                }
            }
            _ss[N] = _rr[N];

            int midc = _cc[0] + _rr[0];
            int midj = 0;
            int type = 1;
            for (int j = 0; j <= N; j++)
                if ((c = _cc[j] + _rr[j]) >= midc && (c > midc || (_cc[j] != _dd[j] && _rr[j] == _ss[j])))
                {
                    midc = c;
                    midj = j;
                }
            for (int j = N; j >= 0; j--)
                if ((c = _dd[j] + _ss[j] + _q) > midc)
                {
                    midc = c;
                    midj = j;
                    type = 2;
                }

            if (type == 1)
            {
                Diff(aOff, bOff, midi, midj, tb, _q);
                Diff(aOff + midi, bOff + midj, M - midi, N - midj, _q, te);
            }
            else
            {
                Diff(aOff, bOff, midi - 1, midj, tb, 0);
                Del(2);
                Diff(aOff + midi + 1, bOff + midj, M - midi - 1, N - midj, 0, te);
            }
        }

        private int Sub1(int ai, int bj) => _sc.Sub(_a[ai], _b[bj]);
    }

    #endregion
}
