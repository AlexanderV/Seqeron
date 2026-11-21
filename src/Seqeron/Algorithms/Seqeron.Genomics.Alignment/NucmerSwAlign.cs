namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The banded affine-gap alignment extender of MUMmer 3 (<c>src/tigr/sw_align.cc</c> / <c>sw_align.hh</c> / <c>sw_alignscore.hh</c>, Phillippy &amp; Delcher, TIGR): the engine <c>postnuc</c> uses to join the
/// matches of a cluster, to extend a cluster off its ends and to fuse neighbouring alignments. Port line by line: the anti-diagonal edit matrix whose nodes keep three scores (delete / insert / match, gap open
/// -10, gap extend -7, match +3 for equal A/C/G/T and -7 otherwise, anything else including N scores as a mismatch), the trimming of nodes that fall more than <c>3 x break length</c> below the best score,
/// the optional hard diagonal band (<c>-B</c>), the forward / backward <i>search</i> (finds where the score stops improving without keeping the matrix) and the targeted forward alignment that records the
/// delta vector (distance to the next insertion / deletion, negative for a deletion in the reference).
/// </summary>
/// <remarks>
/// Sequences are 1-based byte arrays with a zero at index 0 and after the last base (the layout of the C program: <c>A [1 .. N]</c>). Coordinates are 1-based and inclusive.
/// 64-bit arithmetic wraps (the C program's <c>long</c>, no overflow check), as the minimum score <c>-LONG_MAX</c> is added to in the program.
/// </remarks>
public static class NucmerSwAlign
{
    /// <summary>Modus operandi bit: forward (1) or backward (0) direction.</summary>
    public const uint DirectionBit = 0x1;
    /// <summary>Modus operandi bit: search only (no matrix kept, no delta).</summary>
    public const uint SearchBit = 0x2;
    /// <summary>Modus operandi bit: forced alignment (ignore the score, reach the target).</summary>
    public const uint ForcedBit = 0x4;
    /// <summary>Modus operandi bit: maximise the score by shrinking the coverage.</summary>
    public const uint OptimalBit = 0x8;
    /// <summary>Modus operandi bit: extend to the end of the shorter sequence when possible (<c>postnuc -t</c>).</summary>
    public const uint SeqEndBit = 0x10;

    /// <summary>Align forward until the target is reached or the score worsens (creates the delta).</summary>
    public const uint ForwardAlign = 0x1;
    /// <summary>Forward alignment with the score optimised.</summary>
    public const uint OptimalForwardAlign = 0x9;
    /// <summary>Forward alignment forced to reach the target.</summary>
    public const uint ForcedForwardAlign = 0x5;
    /// <summary>Backward search (no delta).</summary>
    public const uint BackwardSearch = 0x2;

    /// <summary>Longest alignment (bases in either sequence) one call may span.</summary>
    public const long MaxAlignmentLength = 10000;

    private const int Delete = 0, Insert = 1, Match = 2, Start = 3, None = 4;
    private const int GoodScore = 3, OpenGapScore = -10, ContGapScore = -7;
    private const char StopChar = 'O';

    private sealed class Diagonal
    {
        public long LBound, RBound;
        public long[]? Value;      // 3 per node
        public byte[]? Used;       // 3 per node
        public byte[]? Max;        // per node: index of the maximum score
    }

    private static bool IsAlpha(int c) => (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');

    private static int ToUpper(int c) => c >= 'a' && c <= 'z' ? c - 32 : c;

    /// <summary><c>MATCH_SCORE [NUCLEOTIDE] [toupper (a) - 'A'] [toupper (b) - 'A']</c> after the non-letters are replaced by the stop character.</summary>
    public static int MatchScore(int a, int b)
    {
        if (!IsAlpha(a)) a = StopChar;
        if (!IsAlpha(b)) b = StopChar;
        a = ToUpper(a);
        b = ToUpper(b);
        return a == b && (a == 'A' || a == 'C' || a == 'G' || a == 'T') ? 3 : -7;
    }

    /// <summary>
    /// <c>alignSearch</c> / <c>alignTarget</c> (<c>_alignEngine</c>): aligns <c>A [Astart ..]</c> with <c>B [Bstart ..]</c> until the target (<paramref name="aEnd"/>, <paramref name="bEnd"/>) is reached or the score has not
    /// improved for <paramref name="breakLength"/> diagonals; <paramref name="aEnd"/> / <paramref name="bEnd"/> are replaced by where the alignment ended. The delta is appended to <paramref name="delta"/> (not for a search).
    /// </summary>
    /// <returns>True when the target was reached.</returns>
    public static bool AlignEngine(byte[] a0, long aStart, ref long aEnd, byte[] b0, long bStart, ref long bEnd, List<long> delta, uint mo, int breakLength = 200, int banding = 0)
    {
        long minScore = -long.MaxValue;
        long highScore = minScore;
        long xhighScore = minScore;
        long maxDiff = (long)GoodScore * breakLength;

        long lbound = 0, rbound = 0;
        long finishCt = 0, finishCDi = 0, xFinishCt = 0, xFinishCDi = 0;
        long dl = 2;
        double dmid;
        double dband = banding / 2.0;

        long aBase, bBase, n, m;
        if ((mo & DirectionBit) != 0)
        {
            aBase = aStart - 1;
            bBase = bStart - 1;
            n = aEnd - aStart + 1;
            m = bEnd - bStart + 1;
        }
        else
        {
            aBase = aStart + 1;
            bBase = bStart + 1;
            n = aStart - aEnd + 1;
            m = bStart - bEnd + 1;
        }

        var diag = new List<Diagonal>(100);
        var d0 = new Diagonal { LBound = lbound, RBound = rbound++ };
        d0.Value = new long[3];
        d0.Used = new byte[3];
        d0.Max = new byte[1];
        d0.Value[Delete] = minScore;
        d0.Value[Insert] = minScore;
        d0.Value[Match] = 0;
        d0.Max[0] = Match;
        d0.Used[Delete] = None;
        d0.Used[Insert] = None;
        d0.Used[Match] = Start;
        diag.Add(d0);

        long l = n < m ? n : m;
        long dct;
        for (dct = 1; dct <= n + m && (dct - finishCt) <= breakLength && lbound <= rbound; dct++)
        {
            var cur = new Diagonal { LBound = lbound, RBound = rbound };
            long ds = rbound - lbound + 1;
            cur.Value = new long[ds * 3];
            cur.Used = new byte[ds * 3];
            cur.Max = new byte[ds];
            diag.Add(cur);

            int iadj, madj;
            if (dct <= n)
            {
                iadj = 0;
                madj = -1;
            }
            else
            {
                iadj = 1;
                madj = dct == n + 1 ? 0 : 1;
            }

            int dadj = iadj - 1;

            long pdct = dct - 1;
            var pd = diag[(int)pdct];
            long pds = pd.RBound - pd.LBound + 1;
            long pdi = lbound + dadj;
            pdi -= pd.LBound;

            long ppdct = dct - 2;
            long ppds, ppdi;
            Diagonal? ppd = null;
            if (ppdct >= 0)
            {
                ppd = diag[(int)ppdct];
                ppds = ppd.RBound - ppd.LBound + 1;
                ppdi = lbound + madj;
                ppdi -= ppd.LBound;
            }
            else
            {
                ppdi = ppds = 0;
            }

            if ((mo & ForcedBit) != 0)
                highScore = minScore;

            for (long cdi = lbound; cdi <= rbound; cdi++)
            {
                long di = cdi - cur.LBound;
                long o = di * 3;

                // DELETE
                if (pdi >= 0 && pdi < pds)
                {
                    long po = pdi * 3;
                    ScoreEdit(cur, o + Delete,
                        pd.Used![po + Delete] == None ? pd.Value![po + Delete] : pd.Value![po + Delete] + ContGapScore,
                        pd.Used[po + Insert] == None ? pd.Value[po + Insert] : pd.Value[po + Insert] + OpenGapScore,
                        pd.Used[po + Match] == None ? pd.Value[po + Match] : pd.Value[po + Match] + OpenGapScore);
                }
                else
                {
                    cur.Value[o + Delete] = minScore;
                    cur.Used[o + Delete] = None;
                }

                pdi++;

                // INSERT
                if (pdi >= 0 && pdi < pds)
                {
                    long po = pdi * 3;
                    ScoreEdit(cur, o + Insert,
                        pd.Used![po + Delete] == None ? pd.Value![po + Delete] : pd.Value![po + Delete] + OpenGapScore,
                        pd.Used[po + Insert] == None ? pd.Value[po + Insert] : pd.Value[po + Insert] + ContGapScore,
                        pd.Used[po + Match] == None ? pd.Value[po + Match] : pd.Value[po + Match] + OpenGapScore);
                }
                else
                {
                    cur.Value[o + Insert] = minScore;
                    cur.Used[o + Insert] = None;
                }

                // MATCH / MIS-MATCH
                if (ppdi >= 0 && ppdi < ppds)
                {
                    long po = ppdi * 3;
                    ScoreEdit(cur, o + Match, ppd!.Value![po + Delete], ppd.Value[po + Insert], ppd.Value[po + Match]);
                    cur.Value[o + Match] = unchecked(cur.Value[o + Match] + ScoreMatch(dct, cdi, a0, aBase, b0, bBase, n, mo));
                }
                else
                {
                    cur.Value[o + Match] = minScore;
                    cur.Used[o + Match] = None;
                }

                ppdi++;

                cur.Max[di] = MaxScore(cur.Value, o);
                if (cur.Value[o + cur.Max[di]] >= highScore)
                {
                    highScore = cur.Value[o + cur.Max[di]];
                    finishCt = dct;
                    finishCDi = cdi;
                }
            }

            // max non-optimal score
            if ((mo & SeqEndBit) != 0 && dct >= l)
            {
                if (l == n)
                {
                    if (lbound == 0)
                    {
                        long v = cur.Value[cur.Max[0]];
                        if (v >= xhighScore)
                        {
                            xhighScore = v;
                            xFinishCt = dct;
                            xFinishCDi = 0;
                        }
                    }
                }
                else
                {
                    if (rbound == m)
                    {
                        long idx = m - cur.LBound;
                        long v = cur.Value[idx * 3 + cur.Max[idx]];
                        if (v >= xhighScore)
                        {
                            xhighScore = v;
                            xFinishCt = dct;
                            xFinishCDi = m;
                        }
                    }
                }
            }

            // in the extender modus operandi the soon-to-be great-grandparent diagonal is dropped
            if ((mo & SearchBit) != 0 && dct > 1)
            {
                var dropped = diag[(int)ppdct];
                dropped.Value = null;
                dropped.Used = null;
                dropped.Max = null;
            }

            // trim hopeless diagonal nodes
            for (long di = 0; di < ds; di++)
            {
                if (unchecked(highScore - cur.Value[di * 3 + cur.Max[di]]) > maxDiff)
                    lbound++;
                else
                    break;
            }

            for (long di = ds - 1; di >= 0; di--)
            {
                if (unchecked(highScore - cur.Value[di * 3 + cur.Max[di]]) > maxDiff)
                    rbound--;
                else
                    break;
            }

            // grow the new diagonal and reset the boundaries
            if (dct < n && dct < m)
            {
                dl++;
                rbound++;
                dmid = (dct + 1) / 2.0;
            }
            else if (dct >= n && dct >= m)
            {
                dl--;
                lbound--;
                dmid = n - (dct + 1) / 2.0;
            }
            else if (dct >= n)
            {
                lbound--;
                dmid = n - (dct + 1) / 2.0;
            }
            else
            {
                rbound++;
                dmid = (dct + 1) / 2.0;
            }

            // hard band
            if (dband > 0)
            {
                long tlb = (long)Math.Ceiling(dmid - dband);
                if (lbound < tlb)
                    lbound = tlb;
                long trb = (long)Math.Floor(dmid + dband);
                if (rbound > trb)
                    rbound = trb;
            }

            if (lbound < 0)
                lbound = 0;
            if (rbound >= dl)
                rbound = dl - 1;
        }

        dct--;

        bool targetReached = false;
        if (dct == n + m)
        {
            if ((~mo & OptimalBit) != 0 || (mo & SeqEndBit) != 0)
            {
                targetReached = true;
                finishCt = n + m;
                finishCDi = 0;
            }
            else if (finishCt == dct)
            {
                targetReached = true;
            }
        }
        else if ((mo & SeqEndBit) != 0 && xFinishCt != 0)
        {
            finishCt = xFinishCt;
            finishCDi = xFinishCDi;
        }

        long aadj = finishCt <= n ? finishCt - finishCDi - 1 : n - finishCDi - 1;
        long badj = finishCt <= n ? finishCDi - 1 : finishCt - n + finishCDi - 1;
        if ((~mo & DirectionBit) != 0)
        {
            aadj *= -1;
            badj *= -1;
        }

        aEnd = aStart + aadj;
        bEnd = bStart + badj;

        if ((~mo & SearchBit) != 0)
            GenerateDelta(diag, finishCt, finishCDi, n, delta);

        return targetReached;
    }

    private static void GenerateDelta(List<Diagonal> diag, long finishCt, long finishCDi, long n, List<long> delta)
    {
        long dct = finishCt;
        long cdi = finishCDi;
        var reversePath = new List<byte>(100);

        long di = cdi - diag[(int)dct].LBound;
        int edit = diag[(int)dct].Max![di];

        while (dct >= 0)
        {
            var dg = diag[(int)dct];
            di = cdi - dg.LBound;
            byte currUsed = dg.Used![di * 3 + edit];

            reversePath.Add((byte)edit);
            switch (edit)
            {
                case Delete:
                    cdi = dct-- <= n ? cdi - 1 : cdi;
                    break;
                case Insert:
                    cdi = dct-- <= n ? cdi : cdi + 1;
                    break;
                case Match:
                    if (dct <= n)
                        cdi--;
                    else if (dct != n + 1)
                        cdi++;
                    dct -= 2;
                    break;
                case Start:
                    dct = -1;
                    break;
                default:
                    throw new InvalidOperationException("ERROR: Invalid edit matrix entry, please file a bug report");
            }

            edit = currUsed;
        }

        long count = 1;
        for (int pi = reversePath.Count - 2; pi >= 0; pi--)
        {
            switch (reversePath[pi])
            {
                case Delete:
                    delta.Add(-count);
                    count = 1;
                    break;
                case Insert:
                    delta.Add(count);
                    count = 1;
                    break;
                case Match:
                    count++;
                    break;
                case Start:
                    break;
                default:
                    throw new InvalidOperationException("ERROR: Invalid path matrix entry, please file a bug report");
            }
        }
    }

    private static byte MaxScore(long[] v, long o)
    {
        if (v[o + Delete] > v[o + Insert])
            return v[o + Delete] > v[o + Match] ? (byte)Delete : (byte)Match;
        return v[o + Insert] > v[o + Match] ? (byte)Insert : (byte)Match;
    }

    private static void ScoreEdit(Diagonal cur, long index, long del, long ins, long mat)
    {
        if (del > ins)
        {
            if (del > mat)
            {
                cur.Value![index] = del;
                cur.Used![index] = Delete;
            }
            else
            {
                cur.Value![index] = mat;
                cur.Used![index] = Match;
            }
        }
        else if (ins > mat)
        {
            cur.Value![index] = ins;
            cur.Used![index] = Insert;
        }
        else
        {
            cur.Value![index] = mat;
            cur.Used![index] = Match;
        }
    }

    private static long ScoreMatch(long dct, long cdi, byte[] a, long aBase, byte[] b, long bBase, long n, uint mo)
    {
        long dir = (mo & DirectionBit) != 0 ? 1 : -1;
        int ac, bc;
        if (dct <= n)
        {
            ac = Get(a, aBase + (dct - cdi) * dir);
            bc = Get(b, bBase + cdi * dir);
        }
        else
        {
            ac = Get(a, aBase + (n - cdi) * dir);
            bc = Get(b, bBase + (dct - n + cdi) * dir);
        }

        return MatchScore(ac, bc);
    }

    private static int Get(byte[] s, long index) => index >= 0 && index < s.Length ? s[index] : 0;
}
