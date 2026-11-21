namespace Seqeron.Genomics.Alignment;

/// <summary>The base-level alignment of a <see cref="MinimapHit"/> (<c>mm_extra_t</c> of minimap2 2.26): the CIGAR and the scores <c>minimap2 -c</c> prints in the tags <c>cg</c>, <c>AS</c>, <c>ms</c>, <c>NM</c> (with the hit), <c>nn</c> and <c>ts</c>.</summary>
/// <param name="Cigar">The CIGAR operations of the aligned part (<c>M I D N =</c> and <c>X</c>; no clipping).</param>
/// <param name="DpScore">The score of the whole alignment (<c>AS</c>).</param>
/// <param name="DpMax">The score of the best-scoring segment (<c>ms</c>).</param>
/// <param name="DpMax2">The best <c>ms</c> among the secondary hits of a primary hit.</param>
/// <param name="AmbiguousBases">The number of ambiguous bases in the alignment (<c>nn</c>).</param>
/// <param name="TranscriptStrand">0 unknown, 1 + , 2 - (<c>ts</c>, splice mode).</param>
public sealed record MinimapAlignment(IReadOnlyList<(char Operation, int Length)> Cigar, int DpScore, int DpMax, int DpMax2, int AmbiguousBases, int TranscriptStrand)
{
    /// <summary>The CIGAR as text (<c>cg:Z:</c>); <c>*</c> when empty.</summary>
    public string CigarString => Cigar.Count == 0 ? "*" : string.Concat(Cigar.Select(c => c.Length.ToString(System.Globalization.CultureInfo.InvariantCulture) + c.Operation));
}

/// <summary>The hits of every segment (read) of a fragment mapped by <see cref="Minimap2Chaining.MapFragment"/>, the reads in their original orientation.</summary>
/// <param name="Segments">The result of every segment, in the order of the input.</param>
/// <param name="FragmentGap">The longest gap on the reference of the fragment chain (<c>frag_gap</c> of the program).</param>
public sealed record MinimapFragmentResult(IReadOnlyList<MinimapMapResult> Segments, int FragmentGap);

public static partial class Minimap2Chaining
{
    private sealed class FragOut
    {
        internal required List<Reg>[] Regs { get; set; }
        internal required Mm128[][] Anchors { get; set; }
        internal int RepLen, FragGap, MidOcc;
    }

    private static int[] SegmentLengths(IReadOnlyList<string> segs) => [.. segs.Select(s => s.Length)];

    // align_regs: the base-level alignment (a no-op without -c), then primary / secondary marking
    private static List<Reg> AlignRegs(MinimapMapOptions o, MmIndex mi, int qlen, string seq, List<Reg> regs, Mm128[] a)
    {
        if (!o.Cigar)
            return regs;
        regs = AlignSkeleton(o, mi, qlen, seq, regs, a);
        if (!o.AllChains)
        {
            SetParent(o.MaskLevel, o.MaskLength, regs, (o.Scoring.Match * 2) + o.Scoring.Mismatch, o.HardMaskLevel, o.AltDrop);
            SelectSub(o.PrimaryRatio, o.Seed.KmerSize * 2, o.BestN, false, (int)(o.MaxGap * 0.8), regs);
            SetSamPri(regs);
        }

        return regs;
    }

    /// <summary><c>mm_map_frag</c> of map.c for 1 .. 255 segments: seeds, chaining (with the segment-aware step score), the long-join / re-chaining, the hits, the primary / secondary selection, the divergence estimate, the base-level alignment, the mapping quality and, for two reads under <c>-c</c>, the pairing.</summary>
    private static FragOut MapFragmentCore(MmIndex mi, SeedIndex ix, IReadOnlyList<string> segs, MinimapMapOptions o, IReadOnlyList<string>? targetNames, string? qname)
    {
        int nSegs = segs.Count;
        int[] qlens = SegmentLengths(segs);
        int qlenSum = qlens.Sum();
        var empty = new FragOut { Regs = [.. Enumerable.Range(0, nSegs).Select(_ => new List<Reg>())], Anchors = [.. Enumerable.Range(0, nSegs).Select(_ => Array.Empty<Mm128>())] };
        if (qlenSum == 0 || nSegs <= 0 || nSegs > 255)
            return empty;
        if (o.MaxQueryLength > 0 && qlenSum > o.MaxQueryLength)
            return empty;
        var chainOptions = o.ChainOptions(qlenSum) with { Segments = nSegs };
        ValidateOptions(chainOptions, nameof(o));
        int[] targetLengths = mi.Length;
        uint hash = qname is not null && !o.NoHashName ? X31Hash(qname) : 0;
        hash ^= WangHash((uint)qlenSum) + WangHash((uint)o.HashSeed);
        hash = WangHash(hash);

        var ctx = PrepareQuerySeeds(ix, segs, o.Seed);
        var seeds = CollectSeedArray(ctx, o.Seed, ctx.MidOcc, targetNames, targetLengths, qname);
        var a = seeds.A;
        int maxChainGapRef = chainOptions.MaxTargetGap;
        float penGap = (float)(o.ChainGapScale * 0.01 * o.Seed.KmerSize), penSkip = (float)(o.ChainSkipScale * 0.01 * o.Seed.KmerSize);
        var chains = o.UseRmq
            ? RmqChain(a, o.MaxGap, o.RmqInnerDistance, o.Bandwidth, o.MaxChainSkip, o.RmqSizeCap, o.MinCount, o.MinChainScore, penGap, penSkip)
            : ChainCore(a.Select(m => m.X).ToArray(), a.Select(m => m.Y).ToArray(), chainOptions);
        var compact = Compact(a, chains);

        int bwLong = Math.Max(o.BandwidthLong, o.Bandwidth); // mm_mapopt_update
        if (bwLong > o.Bandwidth && !o.NoLongJoin && !o.ShortRead && !o.Splice && nSegs == 1 && compact.U.Count > 1)
        {
            // re-chain / long-join for long sequences
            int st = (int)compact.Anchors[0].Y, en = (int)compact.Anchors[(int)(uint)compact.U[0] - 1].Y;
            if (qlenSum - (en - st) > o.RmqRescueSize || en - st > qlenSum * o.RmqRescueRatio)
            {
                int nA = 0;
                foreach (ulong u in compact.U)
                    nA += (int)(uint)u;
                var b = new Mm128[nA];
                Array.Copy(compact.Anchors, b, nA);
                RadixSort(b);
                compact = Compact(b, RmqChain(b, o.MaxGap, o.RmqInnerDistance, bwLong, o.MaxChainSkip, o.RmqSizeCap, o.MinCount, o.MinChainScore, penGap, penSkip));
            }
        }
        else if (o.Seed.MaxOcc > ctx.MidOcc && seeds.RepeatLength > 0 && !o.UseRmq)
        {
            // re-chain with a higher max_occ (mostly for short reads) unless the best chain has all the segments
            bool rechain = false;
            if (compact.U.Count > 0)
            {
                int nChainedSegs = 1, max = 0, maxI = -1, maxOff = -1, off = 0;
                for (int i = 0; i < compact.U.Count; ++i)
                {
                    if (max < (int)(compact.U[i] >> 32))
                    {
                        max = (int)(compact.U[i] >> 32);
                        maxI = i;
                        maxOff = off;
                    }

                    off += (int)(uint)compact.U[i];
                }

                for (int i = 1; i < (int)(uint)compact.U[maxI]; ++i)
                    if ((compact.Anchors[maxOff + i].Y & SeedSegMask) != (compact.Anchors[maxOff + i - 1].Y & SeedSegMask))
                        ++nChainedSegs;
                if (nChainedSegs < nSegs)
                    rechain = true;
            }
            else
            {
                rechain = true;
            }

            if (rechain)
            {
                var again = CollectSeedArray(ctx, o.Seed, o.Seed.MaxOcc, targetNames, targetLengths, qname);
                a = again.A;
                seeds = again;
                compact = Compact(a, ChainCore(a.Select(m => m.X).ToArray(), a.Select(m => m.Y).ToArray(), chainOptions));
            }
        }

        var regs0 = GenRegs(hash, qlenSum, compact.U, compact.Anchors, o.Seed.QueryStrand);
        if (mi.NAlt > 0)
        {
            for (int i = 0; i < regs0.Count; ++i)
                if (mi.IsAlt[regs0[i].Rid])
                    regs0[i].IsAlt = true;
            HitSort(regs0, o.AltDrop);
        }

        if (!o.AllChains)
        {
            SetParent(o.MaskLevel, o.MaskLength, regs0, (o.Scoring.Match * 2) + o.Scoring.Mismatch, o.HardMaskLevel, o.AltDrop);
            if (nSegs <= 1)
                SelectSub(o.PrimaryRatio, o.Seed.KmerSize * 2, o.BestN, true, (int)(o.MaxGap * 0.8), regs0);
            else
                SelectSubMulti(o.PrimaryRatio, 0.2f, 0.7f, maxChainGapRef, o.Seed.KmerSize * 2, o.BestN, nSegs, qlens, regs0);
        }

        if (!o.ShortRead && !o.Seed.QueryStrand)
        {
            EstimateError(targetLengths, qlenSum, regs0, compact.Anchors, seeds.MiniPos);
            int kept = 0;
            for (int i = 0; i < regs0.Count; ++i)
            {
                var r = regs0[i];
                int p = r.Parent;
                if (!r.StrandRetained || r.Div < regs0[p].Div * 5.0f || r.Div < 0.01f)
                    regs0[kept++] = r;
            }

            regs0.RemoveRange(kept, regs0.Count - kept);
        }

        var result = new FragOut { Regs = new List<Reg>[nSegs], Anchors = new Mm128[nSegs][], RepLen = seeds.RepeatLength, FragGap = maxChainGapRef, MidOcc = ctx.MidOcc };
        if (nSegs == 1)
        {
            var regs = AlignRegs(o, mi, qlens[0], segs[0], regs0, compact.Anchors);
            SetMapq(regs, o.MinChainScore, o.Scoring.Match, seeds.RepeatLength, o.ShortRead);
            result.Regs[0] = regs;
            result.Anchors[0] = compact.Anchors;
        }
        else
        {
            var (sregs, sa) = SegGen(hash, nSegs, qlens, regs0, compact.Anchors);
            for (int i = 0; i < nSegs; ++i)
            {
                SetParent(o.MaskLevel, o.MaskLength, sregs[i], (o.Scoring.Match * 2) + o.Scoring.Mismatch, o.HardMaskLevel, o.AltDrop);
                sregs[i] = AlignRegs(o, mi, qlens[i], segs[i], sregs[i], sa[i]);
                SetMapq(sregs[i], o.MinChainScore, o.Scoring.Match, seeds.RepeatLength, o.ShortRead);
            }

            result.Regs = sregs;
            result.Anchors = sa;
            if (nSegs == 2 && o.PairOrientation >= 0 && o.Cigar)
                Pair(maxChainGapRef, o.PairBonus, (o.Scoring.Match * 2) + o.Scoring.Mismatch, o.Scoring.Match, qlens, sregs);
        }

        for (int i = 0; i < nSegs; ++i)
            foreach (var r in result.Regs[i])
                r.A = result.Anchors[i];
        return result;
    }

    private static MinimapHit ToHit(Reg r)
    {
        var members = ToAnchors(r.A.AsSpan(r.As, r.Cnt).ToArray());
        MinimapAlignment? aln = null;
        if (r.P is not null)
            aln = new MinimapAlignment([.. r.P.Cigar.Select(c => ("MIDNSHP=XB"[(int)(c & 0xf)], (int)(c >> 4)))], r.P.DpScore, r.P.DpMax, r.P.DpMax2, r.P.NAmbi, r.P.TransStrand);
        return new MinimapHit(r.Id, r.Parent, r.Rid, r.Rev, r.Qs, r.Qe, r.Rs, r.Re, r.Mlen, r.Blen, r.Mapq, r.Score, r.Parent == r.Id ? r.Subsc : 0, r.Cnt, r.Div, members)
        {
            Alignment = aln, IsInversion = r.Inv, Split = r.Split, IsSamPrimary = r.SamPri, IsProperFragment = r.ProperFrag, IsPairedThrough = r.PeThru,
        };
    }

    private static MmIndex MakeIndex(IReadOnlyList<string> targets, IReadOnlyList<string>? targetNames, MinimapMapOptions o, string? junctionBed, string? altList)
    {
        var mi = MmIndex.Create(targets, targetNames is null ? null : [.. targetNames], o.Seed.KmerSize, o.Seed.HomopolymerCompressed);
        if (altList is not null)
            mi.ReadAlt(altList);
        if (junctionBed is not null)
            mi.ReadJunctionBed(junctionBed);
        return mi;
    }

    /// <summary>
    /// The mapping of a query by minimap2 (<c>mm_map_frag</c> of <c>map.c</c>): seeds (<see cref="CollectSeeds"/>), the chaining dynamic program, the long-join re-chaining for <c>-r b,B</c> with <c>B &gt; b</c>, <c>mm_gen_regs</c>, the primary /
    /// secondary marking (<c>mm_set_parent</c>, <c>mm_select_sub</c>), the divergence estimate (<c>mm_est_err</c>) and the mapping quality (<c>mm_set_mapq</c>): the fields of every PAF line of <c>minimap2 [options] ref.fa query.fa</c> (<c>tp</c>, <c>cm</c>, <c>s1</c>, <c>s2</c>,
    /// <c>dv</c>, <c>rl</c>). With <see cref="MinimapMapOptions.Cigar"/> (<c>-c</c>) also the base-level alignment of every hit (<c>mm_align_skeleton</c>: <see cref="MinimapHit.Alignment"/>, the <c>de</c> divergence, the re-ranking by the DP maximum) and
    /// the mapping quality of aligned hits. One query (one segment); see <see cref="MapFragment"/> for paired reads.
    /// </summary>
    /// <param name="targets">The target (reference) sequences.</param>
    /// <param name="query">The query sequence.</param>
    /// <param name="options">The options (null = the defaults of <c>minimap2 2.26</c>).</param>
    /// <param name="targetNames">Names of the targets (needed by <c>-X</c>), null = none.</param>
    /// <param name="queryName">Name of the query: its hash orders chains of equal score (the program hashes the name; null = hash 0).</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    /// <returns>The hits in the order of the PAF output (descending score, secondary hits filtered by <c>-p</c> / <c>-N</c>).</returns>
    public static MinimapMapResult Map(IReadOnlyList<string> targets, string query, MinimapMapOptions? options = null, IReadOnlyList<string>? targetNames = null, string? queryName = null, CancellationToken cancellationToken = default)
    {
        using var scope = AlignmentCancellation.Enter(cancellationToken);
        var o = options ?? new MinimapMapOptions();
        SequenceAligner.RequireNoNullElements(targets, nameof(targets), "Target");
        ArgumentNullException.ThrowIfNull(query);
        ValidateOptions(o.ChainOptions(query.Length), nameof(options));
        if (query.Length == 0)
            return new MinimapMapResult([], 0, 0);
        var mi = MakeIndex(targets, targetNames, o, null, null);
        var ix = BuildSeedIndex(targets, o.Seed);
        var frag = MapFragmentCore(mi, ix, [query], o, targetNames, queryName);
        return ToResult(frag, 0, o);
    }

    private static MinimapMapResult ToResult(FragOut frag, int segment, MinimapMapOptions o)
    {
        var hits = new List<MinimapHit>(frag.Regs[segment].Count);
        foreach (var r in frag.Regs[segment])
        {
            if (!o.OutputSecondary && r.Id != r.Parent)
                continue;
            hits.Add(ToHit(r));
        }

        return new MinimapMapResult(hits, frag.RepLen, frag.MidOcc);
    }

    /// <summary>
    /// The mapping of a fragment of one or more reads (<c>mm_map_frag</c> with <c>n_segs &gt; 1</c>, the multi-segment mode of <c>minimap2 -x sr ref.fa r1.fq r2.fq</c> or of an interleaved file): the reads are chained together (the segment ids of the anchors, the chaining rule of
    /// <c>lchain.c</c> for steps between reads, <c>mm_select_sub_multi</c>), the fragment chains are split into per-read chains (<c>mm_seg_gen</c>), every read gets its own primary / secondary marking and mapping quality and, with <c>-c</c>, the base-level alignment and
    /// the pairing (<c>mm_pair</c>: <see cref="MinimapMapOptions.PairBonus"/>, the proper-pair flag, the lifted mapping qualities). The reads are reverse-complemented before mapping and their hits flipped back as the program does for <see cref="MinimapMapOptions.PairOrientation"/> (the
    /// <c>--pe-ori</c> of the <c>sr</c> preset is FR, 1).
    /// </summary>
    /// <param name="targets">The target (reference) sequences.</param>
    /// <param name="segments">The reads of the fragment (one for a single read).</param>
    /// <param name="options">The options (null = the defaults of <c>minimap2 2.26</c>).</param>
    /// <param name="targetNames">Names of the targets.</param>
    /// <param name="queryName">Name of the fragment (the name of the first read; the program hashes it).</param>
    /// <param name="junctionBed"><c>--junc-bed</c>: the text of a BED file of annotated junctions (BED12 blocks or BED6 intervals), or null.</param>
    /// <param name="altList"><c>--alt</c>: the text of a list of ALT contig names (the first word of each line), or null.</param>
    /// <param name="cancellationToken">Cancellation token.</param>
    public static MinimapFragmentResult MapFragment(IReadOnlyList<string> targets, IReadOnlyList<string> segments, MinimapMapOptions? options = null, IReadOnlyList<string>? targetNames = null, string? queryName = null,
        string? junctionBed = null, string? altList = null, CancellationToken cancellationToken = default)
    {
        using var scope = AlignmentCancellation.Enter(cancellationToken);
        var o = options ?? new MinimapMapOptions();
        SequenceAligner.RequireNoNullElements(targets, nameof(targets), "Target");
        SequenceAligner.RequireNoNullElements(segments, nameof(segments), "Segment");
        if (segments.Count == 0)
            throw new ArgumentException("At least one read is required.", nameof(segments));
        var mi = MakeIndex(targets, targetNames, o, junctionBed, altList);
        var ix = BuildSeedIndex(targets, o.Seed);
        return MapFragmentIndexed(mi, ix, segments, o, targetNames, queryName);
    }

    // the reads of a pair that the orientation (pe_ori) reverse-complements before mapping
    internal static bool FlipsSegment(MinimapMapOptions o, int n, int j) => n == 2 && ((j == 0 && (o.PairOrientation >> 1 & 1) != 0) || (j == 1 && (o.PairOrientation & 1) != 0));

    /// <summary>The hits (internal) of a fragment with the strand and coordinates of the reverse-complemented reads flipped back, the repeat length of every read and the fragment gap (<c>worker_for</c> of map.c).</summary>
    internal sealed class FragRegs
    {
        internal required List<Reg>[] Regs { get; init; }
        internal required int[] RepLens { get; init; }
        internal int FragGap, MidOcc;
    }

    private static void FlipBack(List<Reg> regs, int ql)
    {
        foreach (var r in regs)
        {
            int t = r.Qs;
            r.Qs = ql - r.Qe;
            r.Qe = ql - t;
            r.Rev = !r.Rev;
        }
    }

    internal static FragRegs MapFragmentRegs(MmIndex mi, SeedIndex ix, IReadOnlyList<string> segments, MinimapMapOptions o, IReadOnlyList<string>? targetNames, string? queryName, IReadOnlyList<string>? segmentNames = null)
    {
        int n = segments.Count;
        var flipped = new string[n];
        for (int j = 0; j < n; j++)
            flipped[j] = FlipsSegment(o, n, j) ? ReverseComplement(segments[j]) : segments[j];
        var regs = new List<Reg>[n];
        var reps = new int[n];
        int gap = 0, mid = 0;
        if (o.NoPairing)
        {
            // --no-pairing (MM_F_INDEPEND_SEG): every read is mapped as a fragment of its own
            for (int j = 0; j < n; j++)
            {
                var one = MapFragmentCore(mi, ix, [flipped[j]], o, targetNames, segmentNames is not null ? segmentNames[j] : queryName);
                regs[j] = one.Regs[0];
                reps[j] = one.RepLen;
                gap = one.FragGap;
                mid = one.MidOcc;
            }
        }
        else
        {
            var frag = MapFragmentCore(mi, ix, flipped, o, targetNames, queryName);
            for (int j = 0; j < n; j++)
            {
                regs[j] = frag.Regs[j];
                reps[j] = frag.RepLen;
            }

            gap = frag.FragGap;
            mid = frag.MidOcc;
        }

        for (int j = 0; j < n; j++)
            if (FlipsSegment(o, n, j))
                FlipBack(regs[j], segments[j].Length);
        return new FragRegs { Regs = regs, RepLens = reps, FragGap = gap, MidOcc = mid };
    }

    internal static MinimapFragmentResult ToFragmentResult(FragRegs fr, MinimapMapOptions o)
    {
        var results = new MinimapMapResult[fr.Regs.Length];
        for (int j = 0; j < results.Length; j++)
        {
            var hits = new List<MinimapHit>(fr.Regs[j].Count);
            foreach (var r in fr.Regs[j])
                if (o.OutputSecondary || r.Id == r.Parent)
                    hits.Add(ToHit(r));
            results[j] = new MinimapMapResult(hits, fr.RepLens[j], fr.MidOcc);
        }

        return new MinimapFragmentResult(results, fr.FragGap);
    }

    internal static MinimapFragmentResult MapFragmentIndexed(MmIndex mi, SeedIndex ix, IReadOnlyList<string> segments, MinimapMapOptions o, IReadOnlyList<string>? targetNames, string? queryName, IReadOnlyList<string>? segmentNames = null) =>
        ToFragmentResult(MapFragmentRegs(mi, ix, segments, o, targetNames, queryName, segmentNames), o);

    private static string ReverseComplement(string s)
    {
        var c = new char[s.Length];
        for (int i = 0; i < c.Length; i++)
        {
            char ch = s[s.Length - 1 - i];
            c[i] = ch switch { 'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'U' => 'A', 'a' => 't', 'c' => 'g', 'g' => 'c', 't' => 'a', 'u' => 'a', _ => ch };
        }

        return new string(c);
    }
}

public static partial class Minimap2Chaining
{
    internal static MmIndex CreateIndex(IReadOnlyList<string> targets, IReadOnlyList<string?> names, MinimapMapOptions o, string? junctionBed, string? altList) =>
        MakeIndex(targets, names.Select(n => n ?? "").ToArray(), o, junctionBed, altList);

    internal static SeedIndex CreateSeedIndex(IReadOnlyList<string> targets, MinimapMapOptions o) => BuildSeedIndex(targets, o.Seed);

    internal static MinimapFragmentResult MapFragmentWithIndex(MmIndex mi, SeedIndex ix, IReadOnlyList<string> segments, MinimapMapOptions o, IReadOnlyList<string?> names, string queryName, IReadOnlyList<string> segmentNames) =>
        MapFragmentIndexed(mi, ix, segments, o, names.Select(n => n ?? "").ToArray(), queryName, segmentNames);
}

public static partial class Minimap2Chaining
{
    internal static int SeedIndexMidOcc(SeedIndex ix) => ix.MidOcc;

    /// <summary>
    /// <c>merge_hits</c> of map.c (<c>--split-prefix</c>): the hits of every index part of one fragment are concatenated (the target ids shifted by the number of targets of the earlier parts), the DP maxima re-estimated (<c>mm_update_dp_max</c>), the hits sorted, the
    /// primary / secondary marking and the mapping quality recomputed from the hits of all parts, and the pairing repeated with the fragment gap of the first part.
    /// </summary>
    internal static FragRegs MergeParts(IReadOnlyList<FragRegs> parts, IReadOnlyList<int> ridShift, int[] qlens, MinimapMapOptions o)
    {
        int n = qlens.Length;
        var regs = new List<Reg>[n];
        var reps = new int[n];
        for (int i = 0; i < n; i++)
        {
            regs[i] = [];
            for (int j = 0; j < parts.Count; j++)
            {
                foreach (var r in parts[j].Regs[i])
                {
                    r.Rid += ridShift[j];
                    regs[i].Add(r);
                }

                reps[i] = Math.Max(reps[i], parts[j].RepLens[i]);
            }

            if (!o.ShortRead && qlens[i] >= o.RankMinLength)
                UpdateDpMax(qlens[i], regs[i], o.RankFraction, o.Scoring.Match, o.Scoring.Mismatch);
            foreach (var r in regs[i])
            {
                if (r.P is not null)
                    r.P.DpMax2 = 0;
                r.Subsc = 0;
                r.NSub = 0;
            }

            HitSort(regs[i], o.AltDrop);
            SetParent(o.MaskLevel, o.MaskLength, regs[i], (o.Scoring.Match * 2) + o.Scoring.Mismatch, o.HardMaskLevel, o.AltDrop);
            if (!o.AllChains)
            {
                SelectSub(o.PrimaryRatio, o.Seed.KmerSize * 2, o.BestN, false, (int)(o.MaxGap * 0.8), regs[i]);
                SetSamPri(regs[i]);
            }

            SetMapq(regs[i], o.MinChainScore, o.Scoring.Match, reps[i], o.ShortRead);
        }

        if (n == 2 && o.PairOrientation >= 0 && o.Cigar)
            Pair(parts[0].FragGap, o.PairBonus, (o.Scoring.Match * 2) + o.Scoring.Mismatch, o.Scoring.Match, qlens, regs);
        return new FragRegs { Regs = regs, RepLens = reps, FragGap = parts[0].FragGap, MidOcc = parts[0].MidOcc };
    }

    internal static FragRegs MapFragmentRegsWithIndex(MmIndex mi, SeedIndex ix, IReadOnlyList<string> segments, MinimapMapOptions o, IReadOnlyList<string?> names, string queryName, IReadOnlyList<string> segmentNames) =>
        MapFragmentRegs(mi, ix, segments, o, names.Select(n => n ?? "").ToArray(), queryName, segmentNames);

    internal static MinimapFragmentResult ToResultWithIndex(FragRegs fr, MinimapMapOptions o) => ToFragmentResult(fr, o);
}
