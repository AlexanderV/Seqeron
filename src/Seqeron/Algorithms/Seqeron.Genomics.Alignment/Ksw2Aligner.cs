using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>Flags of the ksw2 extension functions (<c>KSW_EZ_*</c> of <c>ksw2.h</c>).</summary>
[Flags]
public enum Ksw2Option
{
    /// <summary>No flag.</summary>
    None = 0,
    /// <summary><c>KSW_EZ_SCORE_ONLY</c>: do not record the alignment path / CIGAR.</summary>
    ScoreOnly = 0x01,
    /// <summary><c>KSW_EZ_RIGHT</c>: right-align gaps.</summary>
    Right = 0x02,
    /// <summary><c>KSW_EZ_GENERIC_SC</c>: use the full matrix; without it only match / mismatch (<c>mat[0]</c> / <c>mat[1]</c>) and the last symbol as a wildcard (<c>mat[m*m-1]</c>) are used by the SIMD functions.</summary>
    GenericScore = 0x04,
    /// <summary><c>KSW_EZ_APPROX_MAX</c>: approximate maximum (faster with SSE).</summary>
    ApproxMax = 0x08,
    /// <summary><c>KSW_EZ_APPROX_DROP</c>: approximate Z-drop (with <see cref="ApproxMax"/>).</summary>
    ApproxDrop = 0x10,
    /// <summary><c>KSW_EZ_EXTZ_ONLY</c>: extension alignment only (the path ends at the best cell, or at the query end with the end bonus).</summary>
    ExtensionOnly = 0x40,
    /// <summary><c>KSW_EZ_REV_CIGAR</c>: the CIGAR is reported reversed (last operation first).</summary>
    ReverseCigar = 0x80,
    /// <summary><c>KSW_EZ_SPLICE_FOR</c>: look for the GT-AG signal on the forward strand (<see cref="Ksw2Aligner.ExtendS2"/>).</summary>
    SpliceForward = 0x100,
    /// <summary><c>KSW_EZ_SPLICE_REV</c>: look for the CT-AC signal (the reverse complement of GT-AG).</summary>
    SpliceReverse = 0x200,
    /// <summary><c>KSW_EZ_SPLICE_FLANK</c>: the flanking-sequence bonus of the old splice model.</summary>
    SpliceFlank = 0x400,
    /// <summary><c>KSW_EZ_SPLICE_CMPLX</c>: the miniprot splice model (the default of <c>minimap2 -J 1</c>).</summary>
    SpliceComplex = 0x800,
}

/// <summary>
/// Undefined behaviour of the REFERENCE C program on an input (the port itself is well defined and gives the value the C code produces when the undefined access happens to be harmless); the real-program oracle (B09 F310, ASan + valgrind of
/// <c>ksw_extz</c> / <c>ksw_extd</c> of lh3/ksw2) refuses such inputs because the printed result of an aborted or memory-corrupting run is not a reference value.
/// </summary>
[Flags]
public enum Ksw2ReferenceUndefined
{
    /// <summary>The C function is well defined on this input.</summary>
    None = 0,
    /// <summary>
    /// <c>ksw_extz</c> / <c>ksw_extd</c> line <c>eh[j].h = h1, eh[j].e = KSW_NEG_INF;</c> after a row whose band lies completely to the right of the query (<c>st &gt; qlen</c>, the loop ran no cell, <c>j = st &gt; qlen</c>) writes
    /// beyond the <c>qlen + 1</c> entries of <c>eh</c> (heap buffer overflow; AddressSanitizer: "0 bytes after a block of size 8 * (qlen + 1)", valgrind: invalid write). Reached when <c>tlen &gt; qlen + w + 1</c> and the Z-drop did not stop the loop before that row.
    /// </summary>
    WritesBeyondRowBuffer = 1,
    /// <summary>
    /// The global backtrack starts at the corner <c>(tlen - 1, qlen - 1)</c> although the band of the last row ends left of it (<c>tlen - 1 + w &lt; qlen - 1</c>): <c>ksw_backtrack</c> then reads cells of the <c>kmalloc</c>-ed (uninitialised)
    /// traceback matrix that the DP never wrote (and beyond its end), so the CIGAR is undefined.
    /// </summary>
    ReadsUnwrittenTraceback = 2,
}

/// <summary>The result of a ksw2 extension (<c>ksw_extz_t</c>).</summary>
/// <param name="Max">Best score of the extension (<c>max</c>), 0 when nothing scored.</param>
/// <param name="ZDropped">The extension stopped by the Z-drop rule (or because the band became empty).</param>
/// <param name="MaxQ">Query coordinate (0-based, inclusive) of the best cell, -1 when none.</param>
/// <param name="MaxT">Target coordinate of the best cell, -1 when none.</param>
/// <param name="Mqe">Best score when reaching the end of the query (<c>mqe</c>); <see cref="Ksw2Aligner.NegInf"/> when never reached.</param>
/// <param name="MqeT">Target coordinate of <paramref name="Mqe"/>.</param>
/// <param name="Mte">Best score when reaching the end of the target (<c>mte</c>).</param>
/// <param name="MteQ">Query coordinate of <paramref name="Mte"/>.</param>
/// <param name="Score">Score reaching both ends (global score); <see cref="Ksw2Aligner.NegInf"/> when the corner was not reached.</param>
/// <param name="ReachEnd">The path was taken to the end of the query because <c>mqe + end_bonus &gt; max</c> (SIMD functions, extension mode).</param>
/// <param name="Cigar">The BAM-style CIGAR operations (<c>M</c>, <c>I</c> = query-only, <c>D</c> = target-only) as (operation, length).</param>
public sealed record Ksw2Result(int Max, bool ZDropped, int MaxQ, int MaxT, int Mqe, int MqeT, int Mte, int MteQ, int Score, bool ReachEnd, IReadOnlyList<(char Operation, int Length)> Cigar)
{
    /// <summary>Undefined behaviour that the C reference program has on this input (<see cref="Ksw2ReferenceUndefined.None"/> for every SIMD function and for all inputs the scalar functions handle without an out-of-range access).</summary>
    public Ksw2ReferenceUndefined ReferenceUndefined { get; init; }

    /// <summary>The CIGAR as text (<c>*</c> when empty).</summary>
    public string CigarString => Cigar.Count == 0 ? "*" : string.Concat(Cigar.Select(c => c.Length.ToString(System.Globalization.CultureInfo.InvariantCulture) + c.Operation));
}

/// <summary>
/// The ksw2 extension-alignment functions of Heng Li's <c>ksw2</c> / <c>minimap2</c> (<c>ksw2_extz.c</c>, <c>ksw2_extz2_sse.c</c>, <c>ksw2_extd.c</c>, <c>ksw2_extd2_sse.c</c>): global and extension alignment with
/// affine (<c>q + k e</c>) or dual-affine (<c>min(q + k e, q2 + k e2)</c>) gap cost, a diagonal band (<c>w</c>), the Z-drop rule, the end bonus, ambiguous-base (wildcard) scoring and CIGAR output. The SIMD functions
/// (<c>ksw_extz2_sse</c>, <c>ksw_extd2_sse</c>) are ported lane by lane: the 8-bit difference recurrence of Suzuki &amp; Kasahara (2018) over anti-diagonals, including the 16-lane blocks that run past the band, the
/// byte wrap-around and the shared <c>mem</c> buffer layout, so that the results are those of the SSE4.1 code, not just those of the mathematical DP. Costs are positive penalties (<c>q</c>, <c>e</c>); the matrix is
/// <c>m * m</c> signed bytes indexed <c>mat[target * m + query]</c>.
/// </summary>
public static partial class Ksw2Aligner
{
    /// <summary><c>KSW_NEG_INF</c>.</summary>
    public const int NegInf = -0x40000000;

    private const int CigarMatch = 0;
    private const int CigarIns = 1;
    private const int CigarDel = 2;
    private const int CigarNSkip = 3;

    /// <summary>The <c>ksw_gen_simple_mat</c> matrix of minimap2: <c>a</c> on the diagonal, <c>-b</c> off it, <c>-ambiguous</c> in the last row and column (the wildcard).</summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="m"/> is outside 1 .. 256 (the residue codes are bytes) or a score does not fit the <c>int8_t</c> of the C function.</exception>
    public static sbyte[] SimpleMatrix(int m, int match, int mismatch, int ambiguous)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(m, 1);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(m, 256);
        if (match is < sbyte.MinValue or > sbyte.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(match), match, "A score of the matrix is an int8_t (ksw_gen_simple_mat).");
        if (mismatch is < sbyte.MinValue or > sbyte.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(mismatch), mismatch, "A score of the matrix is an int8_t (ksw_gen_simple_mat).");
        if (ambiguous is < -sbyte.MaxValue or > -sbyte.MinValue)
            throw new ArgumentOutOfRangeException(nameof(ambiguous), ambiguous, "The wildcard score is stored negated in an int8_t (ksw_gen_simple_mat).");
        var mat = new sbyte[m * m];
        for (int i = 0; i < m - 1; i++)
        {
            for (int j = 0; j < m - 1; j++)
                mat[(i * m) + j] = (sbyte)(i == j ? match : mismatch);
            mat[(i * m) + m - 1] = (sbyte)(-ambiguous);
        }

        for (int j = 0; j < m; j++)
            mat[((m - 1) * m) + j] = (sbyte)(-ambiguous);
        return mat;
    }

    /// <summary>The minimap2 nucleotide code of a base: A = 0, C = 1, G = 2, T = 3 (either case), anything else 4 (the wildcard N).</summary>
    public static byte[] EncodeNucleotides(string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var r = new byte[sequence.Length];
        for (int i = 0; i < r.Length; i++)
        {
            r[i] = char.ToUpperInvariant(sequence[i]) switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, _ => 4 };
        }

        return r;
    }

    private sealed class Ez
    {
        internal int Max;
        internal bool ZDropped;
        internal int MaxQ, MaxT, Mqe, MqeT, Mte, MteQ, Score;
        internal bool ReachEnd;
        internal Ksw2ReferenceUndefined Undefined;
        internal List<uint> Cigar = [];

        internal Ez()
        {
            MaxQ = MaxT = MqeT = MteQ = -1;
            Max = 0;
            Score = Mqe = Mte = NegInf;
        }

        internal Ksw2Result ToResult() => new(Max, ZDropped, MaxQ, MaxT, Mqe, MqeT, Mte, MteQ, Score, ReachEnd,
            [.. Cigar.Select(c => ("MIDN"[(int)(c & 0xf)], (int)(c >> 4)))]) { ReferenceUndefined = Undefined };
    }

    private static void PushCigar(List<uint> cigar, uint op, int len)
    {
        if (cigar.Count == 0 || op != (cigar[^1] & 0xf))
            cigar.Add(((uint)len << 4) | op);
        else
            cigar[^1] += (uint)len << 4;
    }

    // the C code reads the backtrack matrix unchecked: a cell outside it (the scalar functions backtracking from a corner that the band does not contain) is undefined there; it reads as 0 here
    private static uint CellOrZero(byte[] p, long index) => index >= 0 && index < p.Length ? p[index] : 0u;

    // ksw_backtrack
    private static List<uint> Backtrack(bool isRot, bool isRev, byte[] p, int[] off, int[]? offEnd, int nCol, int i0, int j0, int minIntronLen = 0)
    {
        var cigar = new List<uint>();
        int i = i0;
        int j = j0;
        int state = 0;
        while (i >= 0 && j >= 0)
        {
            int forceState = -1;
            uint tmp;
            if (isRot)
            {
                int r = i + j;
                if (i < off[r])
                    forceState = 2;
                if (offEnd is not null && i > offEnd[r])
                    forceState = 1;
                tmp = forceState < 0 ? CellOrZero(p, ((long)r * nCol) + i - off[r]) : 0u;
            }
            else
            {
                if (j < off[i])
                    forceState = 2;
                if (offEnd is not null && j > offEnd[i])
                    forceState = 1;
                tmp = forceState < 0 ? CellOrZero(p, ((long)i * nCol) + j - off[i]) : 0u;
            }

            if (state == 0)
                state = (int)(tmp & 7);
            else if (((tmp >> (state + 2)) & 1) == 0)
                state = 0;
            if (state == 0)
                state = (int)(tmp & 7);
            if (forceState >= 0)
                state = forceState;
            if (state == 0)
            {
                PushCigar(cigar, CigarMatch, 1);
                --i;
                --j;
            }
            else if (state == 1 || (state == 3 && minIntronLen <= 0))
            {
                PushCigar(cigar, CigarDel, 1);
                --i;
            }
            else if (state == 3)
            {
                PushCigar(cigar, CigarNSkip, 1);
                --i;
            }
            else
            {
                PushCigar(cigar, CigarIns, 1);
                --j;
            }
        }

        if (i >= 0)
            PushCigar(cigar, (uint)(minIntronLen > 0 && i >= minIntronLen ? CigarNSkip : CigarDel), i + 1);
        if (j >= 0)
            PushCigar(cigar, CigarIns, j + 1);
        if (!isRev)
            cigar.Reverse();
        return cigar;
    }

    // ksw_apply_zdrop
    private static bool ApplyZdrop(Ez ez, bool isRot, int h, int a, int b, int zdrop, int e)
    {
        int r;
        int t;
        if (isRot)
        {
            r = a;
            t = b;
        }
        else
        {
            r = a + b;
            t = a;
        }

        if (h > ez.Max)
        {
            ez.Max = h;
            ez.MaxT = t;
            ez.MaxQ = r - t;
        }
        else if (t >= ez.MaxT && r - t >= ez.MaxQ)
        {
            int tl = t - ez.MaxT;
            int ql = (r - t) - ez.MaxQ;
            int l = tl > ql ? tl - ql : ql - tl;
            if (zdrop >= 0 && ez.Max - h > zdrop + (l * e))
            {
                ez.ZDropped = true;
                return true;
            }
        }

        return false;
    }

    private static void CheckInputs(byte[] query, byte[] target, int m, sbyte[] mat)
    {
        ArgumentNullException.ThrowIfNull(query);
        ArgumentNullException.ThrowIfNull(target);
        ArgumentNullException.ThrowIfNull(mat);
        if (m < 0 || mat.Length < (long)m * m)       // 64-bit: m = int.MaxValue overflowed the product
            throw new ArgumentException("The matrix must have m * m entries.", nameof(mat));
        foreach (byte c in query)
            if (c >= m)
                throw new ArgumentException("Query symbols must be below m.", nameof(query));
        foreach (byte c in target)
            if (c >= m)
                throw new ArgumentException("Target symbols must be below m.", nameof(target));
    }

    /// <summary>
    /// The gap penalties of <c>ksw2.h</c> (<c>gapo</c> / <c>gape</c>: "a gap of length l costs -(gapo + l * gape)", prototype <c>int8_t q, int8_t e</c>): penalties, so <c>q &gt;= 0</c> and <c>e &gt;= 1</c>
    /// (the programs' differences of 8 bits and the oracle need a positive extension penalty), and within <c>int8_t</c> (the C call would truncate silently, here an exception).
    /// </summary>
    internal static void CheckGapPenalties(int q, int e)
    {
        if (q < 0 || q > sbyte.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(q), q, "ksw2: the gap opening penalty q is a penalty between 0 and 127 (int8_t in ksw2.h; a gap of length l costs q + l * e).");
        if (e < 1 || e > sbyte.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(e), e, "ksw2: the gap extension penalty e is a penalty between 1 and 127 (int8_t in ksw2.h; a gap of length l costs q + l * e).");
    }

    // the second piece of the dual affine cost: minimap2 accepts -E k,0 (a long gap at a fixed cost), so e2 may be 0
    internal static void CheckSecondGapPenalties(int q2, int e2)
    {
        if (q2 < 0 || q2 > sbyte.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(q2), q2, "ksw2: the second gap opening penalty q2 is a penalty between 0 and 127 (int8_t in ksw2.h).");
        if (e2 < 0 || e2 > sbyte.MaxValue)
            throw new ArgumentOutOfRangeException(nameof(e2), e2, "ksw2: the second gap extension penalty e2 is a penalty between 0 and 127 (int8_t in ksw2.h; minimap2 accepts 0).");
    }

    /// <summary>The memory of the traceback matrix (<c>ez->cigar</c> source): one byte per cell in the scalar functions, 16 bytes per anti-diagonal column in the SIMD ones; ksw2 allocates it with <c>kmalloc</c> and has no limit.</summary>
    internal static void RequireTraceback(string function, long bytes)
    {
        if (bytes > MaxTracebackBytes)
            throw SequenceAligner.ResourceLimit($"ksw2 {function}", "keeps a traceback matrix for the CIGAR", string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{bytes:N0} bytes"), MaxTracebackBytes, nameof(Ksw2Aligner) + "." + nameof(MaxTracebackBytes),
                "ksw2 itself has no limit and aborts with an out-of-memory error", "use the ScoreOnly flag, a band (w) or the linear-space aligners");
    }

    /// <summary>Largest traceback matrix of <c>ExtendZ / ExtendD / ExtendZ2 / ExtendD2</c> when a CIGAR is requested: 2^30 bytes (1 GiB); <see cref="Ksw2Option.ScoreOnly"/> needs none.</summary>
    public const long MaxTracebackBytes = 1L << 30;

    // ------------------------------------------------------------------ scalar ----

    /// <summary>
    /// <c>ksw_extz</c>: scalar extension / global alignment with affine gap cost <c>q + k e</c>, row-wise (target in the outer loop), band <paramref name="w"/> (&lt; 0: none), Z-drop (&lt; 0: none).
    /// </summary>
    public static Ksw2Result ExtendZ(byte[] query, byte[] target, int m, sbyte[] mat, int q, int e, int w, int zdrop, Ksw2Option flag = Ksw2Option.None)
    {
        CheckInputs(query, target, m, mat);
        CheckGapPenalties(q, e);
        int qlen = query.Length;
        int tlen = target.Length;
        var ez = new Ez();
        bool withCigar = (flag & Ksw2Option.ScoreOnly) == 0;
        int gapo = (sbyte)q;
        int gape = (sbyte)e;
        int gapoe = gapo + gape;
        if (w < 0 || w > Math.Max(qlen, tlen))
            w = tlen > qlen ? tlen : qlen;      // a band wider than the matrix is no band (and 2 * w + 1 must not overflow)
        int nCol = qlen < (2 * w) + 1 ? qlen : (2 * w) + 1;
        var eh = new int[qlen + 1];
        var ee = new int[qlen + 1];
        byte[]? z = null;
        int[]? off = null;
        if (withCigar)
        {
            RequireTraceback(nameof(ExtendZ), (long)nCol * tlen);
            z = new byte[(long)nCol * tlen];
            off = new int[tlen];
        }

        var qp = new sbyte[qlen * m];
        int ip = 0;
        for (int k = 0; k < m; ++k)
            for (int j = 0; j < qlen; ++j)
                qp[ip++] = mat[(k * m) + query[j]];

        eh[0] = 0;
        ee[0] = -gapoe - gapoe;
        int jj;
        for (jj = 1; jj <= qlen && jj <= w; ++jj)
        {
            eh[jj] = -(gapoe + (gape * (jj - 1)));
            ee[jj] = -(gapoe + gapoe + (gape * jj));
        }

        for (; jj <= qlen; ++jj)
            eh[jj] = ee[jj] = NegInf;

        int maxJ = 0;
        for (int i = 0; i < tlen; ++i)
        {
            AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per target row
            int max = NegInf;
            int qo = target[i] * qlen;
            int st = i > w ? i - w : 0;
            int en = i + w < qlen - 1 ? i + w : qlen - 1;
            int h1 = st > 0 ? NegInf : -(gapoe + (gape * i));
            int f = st > 0 ? NegInf : -(gapoe + gapoe + (gape * i));
            bool right = (flag & Ksw2Option.Right) != 0;
            if (withCigar)
                off![i] = st;
            int j = st;
            for (; j <= en; ++j)
            {
                int h = eh[j];
                int e1 = ee[j];
                int d = 0;
                eh[j] = h1;
                h += qp[qo + j];
                if (!withCigar)
                {
                    h = h >= e1 ? h : e1;
                    h = h >= f ? h : f;
                    h1 = h;
                    maxJ = max > h ? maxJ : j;
                    max = max > h ? max : h;
                    h -= gapoe;
                    e1 -= gape;
                    e1 = e1 > h ? e1 : h;
                    ee[j] = e1;
                    f -= gape;
                    f = f > h ? f : h;
                    continue;
                }

                if (!right)
                {
                    d = h >= e1 ? 0 : 1;
                    h = h >= e1 ? h : e1;
                    d = h >= f ? d : 2;
                    h = h >= f ? h : f;
                    h1 = h;
                    maxJ = max > h ? maxJ : j;
                    max = max > h ? max : h;
                    h -= gapoe;
                    e1 -= gape;
                    d |= e1 > h ? 0x08 : 0;
                    e1 = e1 > h ? e1 : h;
                    ee[j] = e1;
                    f -= gape;
                    d |= f > h ? 0x10 : 0;
                    f = f > h ? f : h;
                }
                else
                {
                    d = h > e1 ? 0 : 1;
                    h = h > e1 ? h : e1;
                    d = h > f ? d : 2;
                    h = h > f ? h : f;
                    h1 = h;
                    maxJ = max >= h ? maxJ : j;
                    max = max >= h ? max : h;
                    h -= gapoe;
                    e1 -= gape;
                    d |= e1 >= h ? 0x08 : 0;
                    e1 = e1 >= h ? e1 : h;
                    ee[j] = e1;
                    f -= gape;
                    d |= f >= h ? 0x10 : 0;
                    f = f >= h ? f : h;
                }

                z![((long)i * nCol) + j - st] = (byte)d;
            }

            if (j <= qlen)
            {
                eh[j] = h1;
                ee[j] = NegInf;
            }
            else
                ez.Undefined |= Ksw2ReferenceUndefined.WritesBeyondRowBuffer;     // the C code writes eh[j] here

            if (en == qlen - 1 && eh[qlen] > ez.Mqe)
            {
                ez.Mqe = eh[qlen];
                ez.MqeT = i;
            }

            if (i == tlen - 1)
            {
                ez.Mte = max;
                ez.MteQ = maxJ;
            }

            if (ApplyZdrop(ez, false, max, i, maxJ, zdrop, gape))
                break;
            if (i == tlen - 1 && en == qlen - 1)
                ez.Score = eh[qlen];
        }

        if (withCigar)
        {
            bool rev = (flag & Ksw2Option.ReverseCigar) != 0;
            if (!ez.ZDropped && (flag & Ksw2Option.ExtensionOnly) == 0)
            {
                if (tlen > 0 && qlen > 0 && Math.Min((long)tlen - 1 + w, qlen - 1) < qlen - 1)
                    ez.Undefined |= Ksw2ReferenceUndefined.ReadsUnwrittenTraceback;     // the corner is right of the band of the last row
                ez.Cigar = Backtrack(false, rev, z!, off!, null, nCol, tlen - 1, qlen - 1);
            }
            else if (ez.MaxT >= 0 && ez.MaxQ >= 0)
                ez.Cigar = Backtrack(false, rev, z!, off!, null, nCol, ez.MaxT, ez.MaxQ);
        }

        return ez.ToResult();
    }

    /// <summary>
    /// <c>ksw_extd</c>: scalar extension / global alignment with the dual-affine gap cost <c>min(q + k e, q2 + k e2)</c>.
    /// </summary>
    public static Ksw2Result ExtendD(byte[] query, byte[] target, int m, sbyte[] mat, int q, int e, int q2, int e2, int w, int zdrop, Ksw2Option flag = Ksw2Option.None)
    {
        CheckInputs(query, target, m, mat);
        CheckGapPenalties(q, e);
        CheckSecondGapPenalties(q2, e2);
        int qlen = query.Length;
        int tlen = target.Length;
        var ez = new Ez();
        bool withCigar = (flag & Ksw2Option.ScoreOnly) == 0;
        int gapo = (sbyte)q;
        int gape = (sbyte)e;
        int gapo2 = (sbyte)q2;
        int gape2 = (sbyte)e2;
        int gapoe = gapo + gape;
        int gapoe2 = gapo2 + gape2;
        if (w < 0 || w > Math.Max(qlen, tlen))
            w = tlen > qlen ? tlen : qlen;      // a band wider than the matrix is no band (and 2 * w + 1 must not overflow)
        int nCol = qlen < (2 * w) + 1 ? qlen : (2 * w) + 1;
        var eh = new int[qlen + 1];
        var ee = new int[qlen + 1];
        var ee2 = new int[qlen + 1];
        byte[]? z = null;
        int[]? off = null;
        if (withCigar)
        {
            RequireTraceback(nameof(ExtendD), (long)nCol * tlen);
            z = new byte[(long)nCol * tlen];
            off = new int[tlen];
        }

        var qp = new sbyte[qlen * m];
        int ip = 0;
        for (int k = 0; k < m; ++k)
            for (int j = 0; j < qlen; ++j)
                qp[ip++] = mat[(k * m) + query[j]];

        eh[0] = 0;
        ee[0] = -gapoe - gapoe;
        ee2[0] = -gapoe2 - gapoe2;
        int jj;
        for (jj = 1; jj <= qlen && jj <= w; ++jj)
        {
            eh[jj] = -(gapo + (gape * jj)) > -(gapo2 + (gape2 * jj)) ? -(gapo + (gape * jj)) : -(gapo2 + (gape2 * jj));
            int tmp0 = -(gapoe + (gape * jj)) > -(gapoe2 + (gape2 * jj)) ? -(gapoe + (gape * jj)) : -(gapoe2 + (gape2 * jj));
            ee[jj] = tmp0 - gapoe;
            ee2[jj] = tmp0 - gapoe2;
        }

        for (; jj <= qlen; ++jj)
            eh[jj] = ee[jj] = ee2[jj] = NegInf;

        int maxJ = 0;
        bool right = (flag & Ksw2Option.Right) != 0;
        for (int i = 0; i < tlen; ++i)
        {
            AlignmentCancellation.ThrowIfRequested();       // A3-5: one poll per target row
            int max = NegInf;
            int qo = target[i] * qlen;
            int st = i > w ? i - w : 0;
            int en = i + w < qlen - 1 ? i + w : qlen - 1;
            int tmp = -(gapoe + (gape * i)) > -(gapoe2 + (gape2 * i)) ? -(gapoe + (gape * i)) : -(gapoe2 + (gape2 * i));
            int h1 = st > 0 ? NegInf : tmp;
            int f = st > 0 ? NegInf : tmp - gapoe;
            int f2 = st > 0 ? NegInf : tmp - gapoe2;
            if (withCigar)
                off![i] = st;
            int j = st;
            for (; j <= en; ++j)
            {
                int h = eh[j];
                int e1 = ee[j];
                int e21 = ee2[j];
                int d = 0;
                eh[j] = h1;
                h += qp[qo + j];
                if (!withCigar)
                {
                    h = h >= e1 ? h : e1;
                    h = h >= f ? h : f;
                    h = h >= e21 ? h : e21;
                    h = h >= f2 ? h : f2;
                    h1 = h;
                    maxJ = max > h ? maxJ : j;
                    max = max > h ? max : h;
                    h -= gapoe;
                    e1 -= gape;
                    e1 = e1 > h ? e1 : h;
                    ee[j] = e1;
                    f -= gape;
                    f = f > h ? f : h;
                    int h2a = h1 - gapoe2;
                    e21 -= gape2;
                    e21 = e21 > h2a ? e21 : h2a;
                    ee2[j] = e21;
                    f2 -= gape2;
                    f2 = f2 > h2a ? f2 : h2a;
                    continue;
                }

                int h2;
                if (!right)
                {
                    d = h >= e1 ? 0 : 1;
                    h = h >= e1 ? h : e1;
                    d = h >= f ? d : 2;
                    h = h >= f ? h : f;
                    d = h >= e21 ? d : 3;
                    h = h >= e21 ? h : e21;
                    d = h >= f2 ? d : 4;
                    h = h >= f2 ? h : f2;
                    h1 = h;
                    maxJ = max > h ? maxJ : j;
                    max = max > h ? max : h;
                    h -= gapoe;
                    e1 -= gape;
                    d |= e1 > h ? 1 << 3 : 0;
                    e1 = e1 > h ? e1 : h;
                    ee[j] = e1;
                    f -= gape;
                    d |= f > h ? 1 << 4 : 0;
                    f = f > h ? f : h;
                    h2 = h1 - gapoe2;
                    e21 -= gape2;
                    d |= e21 > h2 ? 1 << 5 : 0;
                    e21 = e21 > h2 ? e21 : h2;
                    ee2[j] = e21;
                    f2 -= gape2;
                    d |= f2 > h2 ? 1 << 6 : 0;
                    f2 = f2 > h2 ? f2 : h2;
                }
                else
                {
                    d = h > e1 ? 0 : 1;
                    h = h > e1 ? h : e1;
                    d = h > f ? d : 2;
                    h = h > f ? h : f;
                    d = h > e21 ? d : 3;
                    h = h > e21 ? h : e21;
                    d = h > f2 ? d : 4;
                    h = h > f2 ? h : f2;
                    h1 = h;
                    maxJ = max > h ? maxJ : j;          // ksw_extd.c keeps "max > h" in its KSW_EZ_RIGHT loop (ksw_extz.c has "max >= h"): the last maximum of a row, as without the flag (B09 F310)
                    max = max > h ? max : h;
                    h -= gapoe;
                    e1 -= gape;
                    d |= e1 >= h ? 1 << 3 : 0;
                    e1 = e1 >= h ? e1 : h;
                    ee[j] = e1;
                    f -= gape;
                    d |= f >= h ? 1 << 4 : 0;
                    f = f >= h ? f : h;
                    h2 = h1 - gapoe2;
                    e21 -= gape2;
                    d |= e21 >= h2 ? 1 << 5 : 0;
                    e21 = e21 >= h2 ? e21 : h2;
                    ee2[j] = e21;
                    f2 -= gape2;
                    d |= f2 >= h2 ? 1 << 6 : 0;
                    f2 = f2 >= h2 ? f2 : h2;
                }

                z![((long)i * nCol) + j - st] = (byte)d;
            }

            if (j <= qlen)
            {
                eh[j] = h1;
                ee[j] = NegInf;
            }
            else
                ez.Undefined |= Ksw2ReferenceUndefined.WritesBeyondRowBuffer;     // the C code writes eh[j] here

            if (en == qlen - 1 && eh[qlen] > ez.Mqe)
            {
                ez.Mqe = eh[qlen];
                ez.MqeT = i;
            }

            if (i == tlen - 1)
            {
                ez.Mte = max;
                ez.MteQ = maxJ;
            }

            if (ApplyZdrop(ez, false, max, i, maxJ, zdrop, gape2))
                break;
            if (i == tlen - 1 && en == qlen - 1)
                ez.Score = eh[qlen];
        }

        if (withCigar)
        {
            bool rev = (flag & Ksw2Option.ReverseCigar) != 0;
            if (!ez.ZDropped && (flag & Ksw2Option.ExtensionOnly) == 0)
            {
                if (tlen > 0 && qlen > 0 && Math.Min((long)tlen - 1 + w, qlen - 1) < qlen - 1)
                    ez.Undefined |= Ksw2ReferenceUndefined.ReadsUnwrittenTraceback;     // the corner is right of the band of the last row
                ez.Cigar = Backtrack(false, rev, z!, off!, null, nCol, tlen - 1, qlen - 1);
            }
            else if (ez.MaxT >= 0 && ez.MaxQ >= 0)
                ez.Cigar = Backtrack(false, rev, z!, off!, null, nCol, ez.MaxT, ez.MaxQ);
        }

        return ez.ToResult();
    }
}
