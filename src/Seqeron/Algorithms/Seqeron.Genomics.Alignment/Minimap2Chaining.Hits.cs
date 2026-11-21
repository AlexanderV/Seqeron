namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The mapping options of minimap2 2.26 up to the end of the chaining stage (A4-7): seeding, chaining, long-join re-chaining and the primary / secondary selection. The defaults are those of <c>mm_mapopt_init</c> for
/// DNA-vs-DNA mapping (<c>map-ont</c>); <c>-c</c> (base-level alignment) is not part of this record. The option name of the program is given for every member.
/// </summary>
public sealed record MinimapMapOptions
{
    /// <summary>The seeding options (<c>-k -w -H -T -f -U -X</c>, <c>--q-occ-frac</c>, <c>--for-only</c> / <c>--rev-only</c>).</summary>
    public MinimapSeedOptions Seed { get; init; } = new();

    /// <summary><c>-g</c>: stop chain elongation if there are no minimizers in this many bases (5000); also the base of the strand-retention score <c>0.8 * max_gap</c>.</summary>
    public int MaxGap { get; init; } = 5000;

    /// <summary><c>-r</c> first value: chaining bandwidth (500).</summary>
    public int Bandwidth { get; init; } = 500;

    /// <summary><c>-r</c> second value: long-join bandwidth (20000); the long-join re-chaining (<c>mg_lchain_rmq</c>) runs only when it exceeds <see cref="Bandwidth"/>.</summary>
    public int BandwidthLong { get; init; } = 20000;

    /// <summary><c>--max-chain-skip</c> (25).</summary>
    public int MaxChainSkip { get; init; } = 25;

    /// <summary><c>--max-chain-iter</c> (5000).</summary>
    public int MaxChainIterations { get; init; } = 5000;

    /// <summary><c>-n</c>: fewest minimizers on a chain (3).</summary>
    public int MinCount { get; init; } = 3;

    /// <summary><c>-m</c>: smallest chaining score (40).</summary>
    public int MinChainScore { get; init; } = 40;

    /// <summary><c>--chain-gap-scale</c> (0.8).</summary>
    public float ChainGapScale { get; init; } = 0.8f;

    /// <summary><c>--chain-skip-scale</c> (0).</summary>
    public float ChainSkipScale { get; init; }

    /// <summary><c>--no-long-join</c>: no long-join re-chaining.</summary>
    public bool NoLongJoin { get; init; }

    /// <summary><c>rmq_inner_dist</c> (1000), <c>rmq_size_cap</c> (100000), <c>rmq_rescue_size</c> (1000) and <c>rmq_rescue_ratio</c> (0.1): the parameters of the RMQ chaining.</summary>
    public int RmqInnerDistance { get; init; } = 1000;

    /// <summary>See <see cref="RmqInnerDistance"/>.</summary>
    public int RmqSizeCap { get; init; } = 100000;

    /// <summary>See <see cref="RmqInnerDistance"/>.</summary>
    public int RmqRescueSize { get; init; } = 1000;

    /// <summary>See <see cref="RmqInnerDistance"/>.</summary>
    public float RmqRescueRatio { get; init; } = 0.1f;

    /// <summary><c>--rmq</c> (the <c>asm</c> presets): chain with the RMQ algorithm instead of the O(n * iterations) dynamic program.</summary>
    public bool UseRmq { get; init; }

    /// <summary><c>--all-chain</c> / <c>-P</c>: report every chain, no primary / secondary selection (the chains are all "primary").</summary>
    public bool AllChains { get; init; }

    /// <summary><c>--mask-level</c> (0.5): a hit whose overlap with a primary hit on the query exceeds this fraction (minus the uncovered fraction) is a secondary hit of it.</summary>
    public float MaskLevel { get; init; } = 0.5f;

    /// <summary><c>--mask-len</c> (<c>int.MaxValue</c>): a hit is secondary only when the length of its query part not covered by primary hits is at most this.</summary>
    public int MaskLength { get; init; } = int.MaxValue;

    /// <summary><c>--hard-mask-level</c>: the overlap test ignores the uncovered length.</summary>
    public bool HardMaskLevel { get; init; }

    /// <summary><c>-p</c> (0.8): secondary hits with a score below this fraction of their primary's are dropped (0 = keep every hit, no selection).</summary>
    public float PrimaryRatio { get; init; } = 0.8f;

    /// <summary><c>-N</c> (5): at most this many secondary hits are retained. (The program changes <c>-N 0</c> to the previous value with <c>--secondary=no</c>: use <see cref="OutputSecondary"/>.)</summary>
    public int BestN { get; init; } = 5;

    /// <summary><c>--secondary=no</c> sets this to false (default true): the secondary hits are selected as usual but not part of the result (the program does not print them).</summary>
    public bool OutputSecondary { get; init; } = true;

    /// <summary>The chaining options of <c>mg_lchain_dp</c> for this mapping.</summary>
    /// <summary><c>-G</c> (<c>max_gap_ref</c>): the longest gap on the reference in a chain, also used as <see cref="Bandwidth"/> and <see cref="BandwidthLong"/> by the command line (<c>-G</c> sets all three; <c>-x splice</c> uses 200 000). 0 = not set (the program's -1).</summary>
    public int MaxGapRef { get; init; }

    /// <summary><c>-F</c> (<c>max_frag_len</c>): the maximum fragment length, effective with <c>-x sr</c> (800): the reference gap of a chain is then <c>max(MaxGap, MaxFragmentLength - query length)</c> unless <see cref="MaxGapRef"/> is set. 0 = not set.</summary>
    public int MaxFragmentLength { get; init; }

    /// <summary><c>--sr</c> (<c>MM_F_SR</c>, set by <c>-x sr</c>): short-read mode: the query gap of a chain is at least the query length, no divergence estimate (<c>dv</c>) and no long-join re-chaining.</summary>
    public bool ShortRead { get; init; }

    /// <summary><c>--splice</c> (<c>MM_F_SPLICE</c>, set by <c>-x splice</c>): spliced alignment. <see cref="Minimap2Chaining.Map"/> chains a spliced query as <c>mg_lchain_dp</c> does with <c>is_cdna</c> (<see cref="ChainOptions"/> sets <c>IsCdna = Splice</c>; the long-join re-chaining is skipped for it, as in <c>mm_map_frag</c>): the <c>splice</c> / <c>splice:hq</c> presets equal <c>minimap2 -x splice</c> on 260 / 260 PAF texts each (A5-9, F173); the base-level stage (<c>Scoring</c>, the splice flags, <c>ksw_exts2_sse</c>) is read by <see cref="Cigar"/> (<c>-c</c>, A3-9c2, F187). See <see cref="MinimapPresets"/>.</summary>
    public bool Splice { get; init; }

    /// <summary>Scoring of the base-level alignment (<c>-A -B -O -E -z -s -C --score-N --junc-bonus</c>, <c>mm_mapopt_t</c> <c>a b q e q2 e2 sc_ambi zdrop zdrop_inv end_bonus min_dp_max noncan junc_bonus</c>): set by the presets, read by the base-level stage (<see cref="Cigar"/>).</summary>
    public MinimapAlignmentScoring Scoring { get; init; } = new();

    /// <summary><c>-c</c> / <c>-a</c> / <c>--cs</c> / <c>--MD</c> (<c>MM_F_CIGAR</c>): the base-level alignment of every chain (<c>mm_align_skeleton</c>: ksw2 between the anchors and the end extensions, <c>mm_fix_cigar</c>, <c>mm_update_extra</c>), the primary / secondary marking after it and, for paired reads, <c>mm_pair</c>.</summary>
    public bool Cigar { get; init; }

    /// <summary><c>--eqx</c> (<c>MM_F_EQX</c>): the CIGAR uses <c>=</c> and <c>X</c> instead of <c>M</c>.</summary>
    public bool Eqx { get; init; }

    /// <summary><c>-u f|r|b</c> (<c>MM_F_SPLICE_FOR</c>): match the GT-AG signal of the transcript strand (the splice presets set both directions).</summary>
    public bool SpliceForward { get; init; }

    /// <summary><c>-u r|b</c> (<c>MM_F_SPLICE_REV</c>): match CT-AC, the reverse complement of GT-AG.</summary>
    public bool SpliceReverse { get; init; }

    /// <summary><c>--splice-flank</c> (<c>MM_F_SPLICE_FLANK</c>, the splice presets set it): the flanking-sequence bonus of the old splice model.</summary>
    public bool SpliceFlank { get; init; }

    /// <summary><c>-J 0</c> (<c>MM_F_SPLICE_OLD</c>): the original minimap2 splice model instead of the miniprot model.</summary>
    public bool SpliceOld { get; init; }

    /// <summary><c>--min-dp-len</c> (<c>min_ksw_len</c>, 200): an anchor gap shorter than this is not aligned by ksw2 on its own.</summary>
    public int MinKswLength { get; init; } = 200;

    /// <summary><c>anchor_ext_len</c> (20): the extension of an end anchor tested by the splice end filter.</summary>
    public int AnchorExtensionLength { get; init; } = 20;

    /// <summary><c>--end-seed-pen</c> (<c>anchor_ext_shift</c>, 6).</summary>
    public int AnchorExtensionShift { get; init; } = 6;

    /// <summary><c>--max-clip-ratio</c> (1): an alignment clipped by more than this fraction of the query at BOTH ends is dropped.</summary>
    public float MaxClipRatio { get; init; } = 1.0f;

    /// <summary><c>rank_min_len</c> (500): queries at least this long get their DP maxima re-estimated (<c>mm_update_dp_max</c>) when the best alignments are similar.</summary>
    public int RankMinLength { get; init; } = 500;

    /// <summary><c>rank_frac</c> (0.9).</summary>
    public float RankFraction { get; init; } = 0.9f;

    /// <summary><c>pe_ori</c> of the paired-read mode (<c>-x sr</c> sets 1 = FR): bit 1 reverse-complements the first read, bit 0 the second; a negative value switches the pairing off.</summary>
    public int PairOrientation { get; init; }

    /// <summary><c>pe_bonus</c> (33): the bonus of a proper pair in <c>mm_pair</c>.</summary>
    public int PairBonus { get; init; } = 33;

    /// <summary><c>--no-end-flt</c> (<c>MM_F_NO_END_FLT</c>): do not drop the first / last anchors that look like noise before the base-level alignment.</summary>
    public bool NoEndFilter { get; init; }

    /// <summary><c>MM_F_NO_INV</c> (set by <c>--qstrand</c>): no inversion rescue (<c>mm_align1_inv</c>).</summary>
    public bool NoInversion { get; init; }

    /// <summary><c>--alt-drop</c> (0.15): the score of a hit on an ALT contig is multiplied by <c>1 - alt_drop</c> in the ranking.</summary>
    public float AltDrop { get; init; } = 0.15f;

    /// <summary><c>--split-prefix</c>: the mapping results of every index part are merged afterwards (<c>mm_split_merge</c>): the DP maxima are then re-estimated by the merge, not by the alignment stage.</summary>
    public bool SplitPrefix { get; init; }

    /// <summary><c>--max-qlen</c> (0 = no limit): a query longer than this is not mapped.</summary>
    public int MaxQueryLength { get; init; }

    /// <summary><c>--seed</c> (11): the seed of the hash that orders chains of equal score.</summary>
    public int HashSeed { get; init; } = 11;

    /// <summary><c>--no-hash-name</c>: the query name is not part of the hash.</summary>
    public bool NoHashName { get; init; }

    /// <summary><c>--frag=yes</c> (<c>MM_F_FRAG_MODE</c>, the <c>sr</c> preset): consecutive reads of a query file with the same name (after removing /1 /2) are the segments of one fragment.</summary>
    public bool FragmentMode { get; init; }

    /// <summary><c>--no-pairing</c> (<c>MM_F_INDEPEND_SEG</c>): the segments of a fragment are mapped independently.</summary>
    public bool NoPairing { get; init; }

    internal MinimapChainOptions ChainOptions(int queryLength = 0)
    {
        int queryGap = ShortRead ? Math.Max(queryLength, MaxGap) : MaxGap;
        int targetGap = MaxGap;
        if (MaxGapRef > 0)
            targetGap = MaxGapRef;
        else if (MaxFragmentLength > 0)
            targetGap = Math.Max(MaxGap, MaxFragmentLength - queryLength);
        return new MinimapChainOptions
        {
            KmerSize = Seed.KmerSize, MaxTargetGap = targetGap, MaxQueryGap = queryGap, Bandwidth = Bandwidth, MaxChainSkip = MaxChainSkip, MaxChainIterations = MaxChainIterations,
            MinCount = MinCount, MinScore = MinChainScore, ChainGapScale = ChainGapScale, ChainSkipScale = ChainSkipScale, IsCdna = Splice,
        };
    }
}

/// <summary>A hit of <see cref="Minimap2Chaining.Map"/>: a chain with the fields of a PAF line of <c>minimap2</c> without <c>-c</c>.</summary>
/// <param name="Id">Index of the hit in the result (the order of the PAF lines).</param>
/// <param name="Parent">Index of the primary hit it belongs to (equal to <paramref name="Id"/> for a primary hit; the PAF <c>tp:A:P</c> / <c>tp:A:S</c>).</param>
/// <param name="TargetId">Index of the target.</param>
/// <param name="Reverse">The hit is on the reverse strand of the query.</param>
/// <param name="QueryStart">Start on the query (0-based, PAF column 3).</param>
/// <param name="QueryEnd">End on the query (PAF column 4).</param>
/// <param name="TargetStart">Start on the target (PAF column 8).</param>
/// <param name="TargetEnd">End on the target (PAF column 9).</param>
/// <param name="MatchLength">The estimated number of matching bases (<c>mlen</c>, PAF column 10).</param>
/// <param name="BlockLength">The estimated alignment block length (<c>blen</c>, PAF column 11).</param>
/// <param name="MapQ">Mapping quality (PAF column 12).</param>
/// <param name="Score">Chain score (<c>s1</c>).</param>
/// <param name="SecondaryScore">Score of the best secondary hit of a primary hit (<c>s2</c>; 0 for a secondary hit, where the tag is not written).</param>
/// <param name="AnchorCount">Number of anchors (<c>cm</c>).</param>
/// <param name="Divergence">Divergence estimated from the minimizers (<c>dv</c>; negative: the tag is not written).</param>
/// <param name="Anchors">The anchors of the chain in increasing order along the target.</param>
/// <summary>The base-level alignment scoring of <c>mm_mapopt_t</c> (<c>-A</c> match, <c>-B</c> mismatch, <c>-O</c> gap open, <c>-E</c> gap extension (and the second pair of the dual affine gap cost), <c>-z</c> Z-drop, <c>-s</c> minimum peak DP score, <c>-u</c> / <c>-C</c> splice options): defaults of <c>mm_mapopt_init</c>.</summary>
public sealed record MinimapAlignmentScoring
{
    public int Match { get; init; } = 2;
    public int Mismatch { get; init; } = 4;
    public int GapOpen { get; init; } = 4;
    public int GapExtend { get; init; } = 2;
    public int GapOpen2 { get; init; } = 24;
    public int GapExtend2 { get; init; } = 1;
    public int ZDrop { get; init; } = 400;
    public int ZDropInverted { get; init; } = 200;
    public int EndBonus { get; init; } = -1;
    /// <summary><c>-s</c> (<c>min_dp_max</c>); null = 80, the value of <c>mm_mapopt_init</c> (<c>min_chain_score * a</c> of the DEFAULTS 40 * 2: the program does NOT recompute it after <c>-m</c> / <c>-A</c>; the presets set their own value or keep 80).</summary>
    public int? MinDpMax { get; init; }

    /// <summary><c>--score-N</c> (<c>sc_ambi</c>, 1): the penalty of an ambiguous base (N); the matrix holds <c>-|N|</c>.</summary>
    public int Ambiguous { get; init; } = 1;
    public int NonCanonical { get; init; }
    public int JunctionBonus { get; init; }
    public int MaxSmithWatermanMatrix { get; init; } = 100_000_000;
}

public sealed record MinimapHit(
    int Id, int Parent, int TargetId, bool Reverse, int QueryStart, int QueryEnd, int TargetStart, int TargetEnd, int MatchLength, int BlockLength, int MapQ,
    int Score, int SecondaryScore, int AnchorCount, float Divergence, IReadOnlyList<MinimapAnchor> Anchors)
{
    /// <summary>The base-level alignment (with <see cref="MinimapMapOptions.Cigar"/>), null for a hit without one (the chaining stage only).</summary>
    public MinimapAlignment? Alignment { get; init; }

    /// <summary>The hit is an inversion rescued between two parts of a chain (<c>tp:A:I</c> / <c>i</c>).</summary>
    public bool IsInversion { get; init; }

    /// <summary>The chain was split by a Z-drop (<c>zd:i:</c>): 1 = the first part, 2 = the second part.</summary>
    public int Split { get; init; }

    /// <summary>The hit is the primary hit of its read that SAM marks as non-supplementary (the first primary hit).</summary>
    public bool IsSamPrimary { get; init; }

    /// <summary>The hit is part of a proper pair chosen by <c>mm_pair</c> (SAM flag 0x2).</summary>
    public bool IsProperFragment { get; init; }

    /// <summary>Both reads of a pair cover the same reference interval (<c>pe_thru</c>).</summary>
    public bool IsPairedThrough { get; init; }

    /// <summary>The hit is a primary hit (<c>tp:A:P</c>).</summary>
    public bool IsPrimary => Id == Parent;
}

/// <summary>The result of <see cref="Minimap2Chaining.Map"/>.</summary>
/// <param name="Hits">The hits in the order of the PAF output.</param>
/// <param name="RepeatLength">Length of the query regions with frequent minimizers (the PAF <c>rl</c> tag).</param>
/// <param name="MidOcc">The occurrence limit that was applied.</param>
public sealed record MinimapMapResult(IReadOnlyList<MinimapHit> Hits, int RepeatLength, int MidOcc);

public static partial class Minimap2Chaining
{
    private const int ParentUnset = -1;

    // mm_reg1_t of minimap2 without a base-level alignment (the members that exist before -c)
    internal sealed class Reg
    {
        public int Id, Parent, Score, Score0, Cnt, As, Rid, Rs, Re, Qs, Qe, Mlen, Blen, Subsc, NSub, Mapq, Split;
        public uint Hash;
        public bool Rev, StrandRetained, Inv, SamPri, ProperFrag, PeThru, SegSplit, SplitInv, IsAlt;
        public float Div = -1f;
        public Extra? P;
        public Mm128[] A = [];

        public Reg Clone() => (Reg)MemberwiseClone();
    }

    // mm_extra_t: the base-level alignment of a hit
    internal sealed class Extra
    {
        public int DpScore, DpMax, DpMax2, NAmbi, TransStrand;
        public List<uint> Cigar = [];
    }

    private const int ParentTmpPri = -2;

    // __ac_X31_hash_string of khash.h: `h = (khint_t)*s; for (++s; *s; ++s) h = (h << 5) - h + (khint_t)*s` over a SIGNED `char` (a byte >= 0x80 is sign-extended, also the first one) and over a C string, which ends at the first NUL.
    private static uint X31Hash(string s)
    {
        var bytes = System.Text.Encoding.UTF8.GetBytes(s);
        int length = Array.IndexOf(bytes, (byte)0);
        if (length < 0)
            length = bytes.Length;
        if (length == 0)
            return 0;
        uint h = (uint)(sbyte)bytes[0];
        for (int i = 1; i < length; i++)
            h = (h << 5) - h + (uint)(sbyte)bytes[i];
        return h;
    }

    // __ac_Wang_hash of khash.h
    private static uint WangHash(uint key)
    {
        key += ~(key << 15);
        key ^= key >> 10;
        key += key << 3;
        key ^= key >> 6;
        key += ~(key << 11);
        key ^= key >> 16;
        return key;
    }

    // hit.c hash64 (no mask)
    private static ulong Hash64Full(ulong key)
    {
        key = ~key + (key << 21);
        key ^= key >> 24;
        key = key + (key << 3) + (key << 8);
        key ^= key >> 14;
        key = key + (key << 2) + (key << 4);
        key ^= key >> 28;
        key += key << 31;
        return key;
    }

    // mm_cal_fuzzy_len + mm_reg_set_coor
    private static void SetCoordinates(Reg r, int qlen, Mm128[] a, bool isQstrand = false)
    {
        int k = r.As, qSpan = (int)(a[k].Y >> 32 & 0xff);
        r.Rev = (a[k].X >> 63) != 0;
        r.Rid = (int)(a[k].X << 1 >> 33);
        r.Rs = (int)a[k].X + 1 > qSpan ? (int)a[k].X + 1 - qSpan : 0;
        r.Re = (int)a[k + r.Cnt - 1].X + 1;
        if (!r.Rev || isQstrand)
        {
            r.Qs = (int)a[k].Y + 1 - qSpan;
            r.Qe = (int)a[k + r.Cnt - 1].Y + 1;
        }
        else
        {
            r.Qs = qlen - ((int)a[k + r.Cnt - 1].Y + 1);
            r.Qe = qlen - ((int)a[k].Y + 1 - qSpan);
        }

        r.Mlen = r.Blen = 0;
        if (r.Cnt <= 0)
            return;
        r.Mlen = r.Blen = (int)(a[r.As].Y >> 32 & 0xff);
        for (int i = r.As + 1; i < r.As + r.Cnt; ++i)
        {
            int span = (int)(a[i].Y >> 32 & 0xff);
            int tl = (int)a[i].X - (int)a[i - 1].X;
            int ql = (int)a[i].Y - (int)a[i - 1].Y;
            r.Blen += tl > ql ? tl : ql;
            int step = Math.Min(tl, ql);
            r.Mlen += tl > span && ql > span ? span : step;
        }
    }

    /// <summary><c>mm_gen_regs</c>: the chains (<paramref name="u"/>: score &lt;&lt; 32 | anchor count, in the order of the compacted array <paramref name="a"/>) as hits sorted by score.</summary>
    private static List<Reg> GenRegs(uint hash, int qlen, List<ulong> u, Mm128[] a, bool isQstrand = false)
    {
        int nU = u.Count;
        if (nU == 0)
            return [];
        var z = new Mm128[nU];
        int k = 0;
        for (int i = 0; i < nU; ++i)
        {
            uint h = (uint)Hash64Full((Hash64Full(a[k].X) + Hash64Full(a[k].Y)) ^ hash);
            z[i] = new Mm128(u[i] ^ h, (ulong)k << 32 | (uint)(int)u[i]);
            k += (int)u[i];
        }

        RadixSort(z);
        Array.Reverse(z);
        var regs = new List<Reg>(nU);
        for (int i = 0; i < nU; ++i)
        {
            var r = new Reg { Id = i, Parent = ParentUnset, Score = (int)(z[i].X >> 32), Hash = (uint)z[i].X, Cnt = (int)z[i].Y, As = (int)(z[i].Y >> 32) };
            r.Score0 = r.Score;
            r.Div = -1f;
            SetCoordinates(r, qlen, a, isQstrand);
            regs.Add(r);
        }

        return regs;
    }

    private static int AltScore(int score, float altDiffFrac)
    {
        if (score < 0)
            return score;
        score = (int)((score * (1.0 - altDiffFrac)) + .499);
        return score > 0 ? score : 1;
    }

    /// <summary><c>mm_set_parent</c>: marks every hit as primary (parent = itself) or secondary of the first primary hit it overlaps on the query by more than the mask level; also the best secondary score (<c>subsc</c>), the number of
    /// suboptimal hits (<c>n_sub</c>) and, for aligned hits, <c>dp_max2</c>.</summary>
    private static void SetParent(float maskLevel, int maskLen, List<Reg> r, int subDiff, bool hardMaskLevel, float altDiffFrac)
    {
        int n = r.Count;
        if (n <= 0)
            return;
        for (int i = 0; i < n; ++i)
            r[i].Id = i;
        var cov = new ulong[n];
        var w = new int[n];
        w[0] = 0;
        r[0].Parent = 0;
        int k = 1;
        for (int i = 1; i < n; ++i)
        {
            var ri = r[i];
            int si = ri.Qs, ei = ri.Qe, nCov = 0, uncovLen = 0;
            bool isNewPrimary = false;
            if (!hardMaskLevel)
            {
                for (int j = 0; j < k; ++j)
                {
                    var rp = r[w[j]];
                    int sj = rp.Qs, ej = rp.Qe;
                    if (ej <= si || sj >= ei)
                        continue;
                    if (sj < si)
                        sj = si;
                    if (ej > ei)
                        ej = ei;
                    cov[nCov++] = (ulong)(uint)sj << 32 | (uint)ej;
                }

                if (nCov == 0)
                {
                    isNewPrimary = true; // no overlapping primary hit: i is a new primary hit (goto set_parent_test)
                }
                else
                {
                    int x = si;
                    Array.Sort(cov, 0, nCov);
                    for (int j = 0; j < nCov; ++j)
                    {
                        if ((int)(cov[j] >> 32) > x)
                            uncovLen += (int)(cov[j] >> 32) - x;
                        x = (int)(uint)cov[j] > x ? (int)(uint)cov[j] : x;
                    }

                    if (ei > x)
                        uncovLen += ei - x;
                }
            }

            if (!isNewPrimary)
            {
                int jj;
                for (jj = 0; jj < k; ++jj)
                {
                    var rp = r[w[jj]];
                    int sj = rp.Qs, ej = rp.Qe;
                    if (ej <= si || sj >= ei)
                        continue;
                    int min = ej - sj < ei - si ? ej - sj : ei - si;
                    int max = ej - sj > ei - si ? ej - sj : ei - si;
                    int ol = OverlapLength(si, ei, sj, ej);
                    if ((float)ol / min - (float)uncovLen / max > maskLevel && uncovLen <= maskLen)
                    {
                        int cntSub = 0, sci = ri.Score;
                        ri.Parent = rp.Parent;
                        if (!rp.IsAlt && ri.IsAlt)
                            sci = AltScore(sci, altDiffFrac);
                        rp.Subsc = rp.Subsc > sci ? rp.Subsc : sci;
                        if (ri.Cnt >= rp.Cnt)
                            cntSub = 1;
                        if (rp.P is not null && ri.P is not null && (rp.Rid != ri.Rid || rp.Rs != ri.Rs || rp.Re != ri.Re || ol != min))
                        {
                            sci = ri.P.DpMax;
                            if (!rp.IsAlt && ri.IsAlt)
                                sci = AltScore(sci, altDiffFrac);
                            rp.P.DpMax2 = rp.P.DpMax2 > sci ? rp.P.DpMax2 : sci;
                            if (rp.P.DpMax - ri.P.DpMax <= subDiff)
                                cntSub = 1;
                        }

                        if (cntSub != 0)
                            ++rp.NSub;
                        break;
                    }
                }

                isNewPrimary = jj == k;
            }

            if (isNewPrimary)
            {
                w[k++] = i;
                ri.Parent = i;
                ri.NSub = 0;
            }
        }
    }

    // the overlap of [si, ei) and [sj, ej) as hit.c writes it
    private static int OverlapLength(int si, int ei, int sj, int ej)
    {
        if (si < sj)
        {
            if (ei < sj)
                return 0;
            return Math.Min(ei, ej) - sj;
        }

        if (ej < si)
            return 0;
        return Math.Min(ej, ei) - si;
    }

    // mm_set_sam_pri
    private static void SetSamPri(List<Reg> r)
    {
        int nPri = 0;
        foreach (var x in r)
        {
            if (x.Id == x.Parent)
            {
                ++nPri;
                x.SamPri = nPri == 1;
            }
            else
            {
                x.SamPri = false;
            }
        }
    }

    // mm_sync_regs: keep id / parent in sync after hits were removed (and reset the ids)
    private static void SyncRegs(List<Reg> regs)
    {
        if (regs.Count == 0)
            return;
        int maxId = -1;
        foreach (var r in regs)
            maxId = Math.Max(maxId, r.Id);
        var tmp = new int[maxId + 1];
        Array.Fill(tmp, -1);
        for (int i = 0; i < regs.Count; ++i)
            if (regs[i].Id >= 0)
                tmp[regs[i].Id] = i;
        for (int i = 0; i < regs.Count; ++i)
        {
            var r = regs[i];
            r.Id = i;
            if (r.Parent == ParentTmpPri)
                r.Parent = i;
            else if (r.Parent >= 0 && tmp[r.Parent] >= 0)
                r.Parent = tmp[r.Parent];
            else
                r.Parent = ParentUnset;
        }

        SetSamPri(regs);
    }

    /// <summary><c>mm_select_sub</c>: keep the primary hits and at most <paramref name="bestN"/> secondary hits that score at least <paramref name="priRatio"/> of their primary (or within <paramref name="minDiff"/> of it); hits of the
    /// other strand with a score above <paramref name="minStrandScore"/> are retained too when <paramref name="checkStrand"/>.</summary>
    private static void SelectSub(float priRatio, int minDiff, int bestN, bool checkStrand, int minStrandScore, List<Reg> r)
    {
        if (!(priRatio > 0.0f) || r.Count <= 0)
            return;
        int n = r.Count, k = 0, n2nd = 0;
        for (int i = 0; i < n; ++i)
        {
            var ri = r[i];
            int p = ri.Parent;
            if (p == i || ri.Inv)
            {
                r[k++] = ri;
            }
            else if ((ri.Score >= r[p].Score * priRatio || ri.Score + minDiff >= r[p].Score) && n2nd < bestN)
            {
                var rp = r[p];
                if (!(ri.Qs == rp.Qs && ri.Qe == rp.Qe && ri.Rid == rp.Rid && ri.Rs == rp.Rs && ri.Re == rp.Re))
                {
                    r[k++] = ri;
                    ++n2nd;
                }
            }
            else if (checkStrand && n2nd < bestN && ri.Score > minStrandScore && ri.Rev != r[p].Rev)
            {
                ri.StrandRetained = true;
                r[k++] = ri;
                ++n2nd;
            }
        }

        if (k != n)
        {
            r.RemoveRange(k, n - k);
            SyncRegs(r);
        }
    }

    // esterr.c: the divergence of a hit estimated from the minimizers it shares with the query
    private static int ForwardQueryPosition(int qlen, Mm128 a)
    {
        int x = (int)a.Y, qSpan = (int)(a.Y >> 32 & 0xff);
        if ((a.X >> 63) != 0)
            x = qlen - 1 - (x + 1 - qSpan);
        return x;
    }

    private static int MiniIndex(int qlen, Mm128 a, int n, ulong[] miniPos)
    {
        int x = ForwardQueryPosition(qlen, a), lo = 0, hi = n - 1;
        while (lo <= hi)
        {
            int m = (int)(((ulong)lo + (ulong)hi) >> 1);
            int y = (int)miniPos[m];
            if (y < x)
                lo = m + 1;
            else if (y > x)
                hi = m - 1;
            else
                return m;
        }

        return -1;
    }

    private static void EstimateError(IReadOnlyList<int> targetLengths, int qlen, List<Reg> regs, Mm128[] a, ulong[] miniPos)
    {
        int n = miniPos.Length;
        if (n == 0)
            return;
        ulong sumK = 0;
        for (int i = 0; i < n; ++i)
            sumK += miniPos[i] >> 32 & 0xff;
        float avgK = (float)sumK / n;
        foreach (var r in regs)
        {
            r.Div = -1f;
            if (r.Cnt == 0)
                continue;
            int st = MiniIndex(qlen, r.Rev ? a[r.As + r.Cnt - 1] : a[r.As], n, miniPos);
            int en = st;
            if (st < 0)
                continue;
            int lRef = targetLengths[r.Rid];
            int nMatch = 1;
            for (int k = 1, j = st + 1; j < n && k < r.Cnt; ++j)
            {
                int x = ForwardQueryPosition(qlen, r.Rev ? a[r.As + r.Cnt - 1 - k] : a[r.As + k]);
                if (x == (int)miniPos[j])
                {
                    ++k;
                    en = j;
                    ++nMatch;
                }
            }

            int nTot = en - st + 1;
            if (r.Qs > avgK && r.Rs > avgK)
                ++nTot;
            if (qlen - r.Qs > avgK && lRef - r.Re > avgK)
                ++nTot;
            r.Div = nMatch >= nTot ? 0.0f : (float)(1.0 - Math.Pow((double)nMatch / nTot, 1.0 / avgK));
        }
    }

    // mm_set_inv_mapq
    private static void SetInvMapq(List<Reg> regs)
    {
        if (regs.Count < 3)
            return;
        if (!regs.Any(x => x.Inv))
            return;
        var aux = new List<Mm128>();
        for (int i = 0; i < regs.Count; ++i)
            if (regs[i].Parent == i || regs[i].Parent < 0)
                aux.Add(new Mm128((ulong)(uint)regs[i].Rid << 32 | (uint)regs[i].Rs, (ulong)i));
        var arr = aux.ToArray();
        RadixSort(arr);
        for (int i = 1; i < arr.Length - 1; ++i)
        {
            var inv = regs[(int)arr[i].Y];
            if (inv.Inv)
            {
                var l = regs[(int)arr[i - 1].Y];
                var r = regs[(int)arr[i + 1].Y];
                inv.Mapq = l.Mapq < r.Mapq ? l.Mapq : r.Mapq;
            }
        }
    }

    // hit.c mm_set_mapq
    private static void SetMapq(List<Reg> regs, int minChainScore, int matchSc, int repLen, bool isSr)
    {
        const float qCoef = 40.0f;
        if (regs.Count == 0)
            return;
        long sumSc = 0;
        foreach (var r in regs)
            if (r.Parent == r.Id)
                sumSc += r.Score;
        float uniqRatio = (float)sumSc / (sumSc + repLen);
        for (int i = 0; i < regs.Count; ++i)
        {
            var r = regs[i];
            if (r.Inv)
            {
                r.Mapq = 0;
            }
            else if (r.Parent == r.Id)
            {
                int mapq;
                float penS1 = (r.Score > 100 ? 1.0f : 0.01f * r.Score) * uniqRatio;
                float penCm = r.Cnt > 10 ? 1.0f : 0.1f * r.Cnt;
                penCm = penS1 < penCm ? penS1 : penCm;
                int subsc = r.Subsc > minChainScore ? r.Subsc : minChainScore;
                if (r.P is not null && r.P.DpMax2 > 0 && r.P.DpMax > 0)
                {
                    float identity = (float)r.Mlen / r.Blen;
                    float x = (float)r.P.DpMax2 * subsc / r.P.DpMax / r.Score0;
                    mapq = CInt(identity * penCm * qCoef * (1.0f - (x * x)) * MathF.Log((float)r.P.DpMax / matchSc));
                    if (!isSr)
                    {
                        int mapqAlt = CInt((6.02f * identity * identity * (r.P.DpMax - r.P.DpMax2) / matchSc) + .499f);
                        mapq = mapq < mapqAlt ? mapq : mapqAlt;
                    }
                }
                else
                {
                    float x = (float)subsc / r.Score0;
                    if (r.P is not null)
                    {
                        float identity = (float)r.Mlen / r.Blen;
                        mapq = CInt(identity * penCm * qCoef * (1.0f - x) * MathF.Log((float)r.P.DpMax / matchSc));
                    }
                    else
                    {
                        mapq = CInt(penCm * qCoef * (1.0f - x) * MathF.Log(r.Score));
                    }
                }

                mapq -= CInt((4.343f * MathF.Log(r.NSub + 1)) + .499f);
                mapq = mapq > 0 ? mapq : 0;
                r.Mapq = mapq < 60 ? mapq : 60;
                if (r.P is not null && r.P.DpMax > r.P.DpMax2 && r.Mapq == 0)
                    r.Mapq = 1;
            }
            else
            {
                r.Mapq = 0;
            }
        }

        SetInvMapq(regs);
    }

    // mm_hit_sort: sort by the DP maximum (or the chaining score without alignments) and the hash, best first; hits with cnt == 0 are removed
    private static void HitSort(List<Reg> r, float altDiffFrac)
    {
        int n = r.Count;
        if (n <= 1)
            return;
        var aux = new List<Mm128>(n);
        for (int i = 0; i < n; ++i)
        {
            if (r[i].Inv || r[i].Cnt > 0)
            {
                int score = r[i].P is not null ? r[i].P!.DpMax : r[i].Score;
                if (r[i].IsAlt)
                    score = AltScore(score, altDiffFrac);
                aux.Add(new Mm128((ulong)(uint)score << 32 | r[i].Hash, (ulong)i));
            }
        }

        var arr = aux.ToArray();
        RadixSort(arr);
        var t = new List<Reg>(arr.Length);
        for (int i = arr.Length - 1; i >= 0; --i)
            t.Add(r[(int)arr[i].Y]);
        r.Clear();
        r.AddRange(t);
    }

    /// <summary>The chains as <c>compact_a</c> leaves them: the anchors of the chains one after the other (ascending inside a chain) and <c>u[i] = score &lt;&lt; 32 | count</c>.</summary>
    private static (Mm128[] Anchors, List<ulong> U) Compact(Mm128[] a, IReadOnlyList<MinimapChain> chains)
    {
        var b = new List<Mm128>();
        var u = new List<ulong>(chains.Count);
        foreach (var c in chains)
        {
            foreach (int i in c.Anchors)
                b.Add(a[i]);
            u.Add((ulong)(uint)c.Score << 32 | (uint)c.Anchors.Count);
        }

        return (b.ToArray(), u);
    }

}
