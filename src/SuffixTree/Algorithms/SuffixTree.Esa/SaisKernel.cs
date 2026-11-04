namespace SuffixTree.Esa;

/// <summary>Which suffix-array implementation a harness runs: the safe permanent oracle or the pointer kernel (iteration 1o.1b).</summary>
internal enum SaisKernelKind
{
    /// <summary><see cref="SaisReference"/>: safe code, generic over the index type, never optimized.</summary>
    Reference,

    /// <summary><see cref="SaisKernel"/>: <c>int</c>-only kernel; all passes are pointer loops behind <see cref="Buf{T}"/>.</summary>
    Pointer,
}

/// <summary>
/// Iteration 1o.1b (REVIEW-S1a section 6): SA-IS with the LINEAR passes as pointer loops, the same algorithm, buffer layout and result as
/// <see cref="SaisReference{TIdx}"/> (bit for bit, equal recursion trace). Pointer passes: symbol counting, type classification fused with the LMS count, bucket
/// setup, placement of the LMS suffixes into their buckets (both rounds), compaction of the sorted LMS positions and of the names, the reduced-problem
/// bookkeeping (direct placement, LMS positions in text order, translation of the reduced ranks); iteration 1o.1c added the four induction passes and the naming
/// (<see cref="SaisInduction"/>, branch-free pointer loops). The symbols and types of a level are read through <see cref="ITypedText"/>: a byte text of fewer than
/// <see cref="PackedBytes.Limit"/> symbols is classified into ONE byte per position that holds symbol and type (<see cref="PackedBytes"/>, it replaces the type array,
/// so the memory is unchanged), every other text keeps a type byte per position next to the symbols (<see cref="SplitTyped{TText}"/>; <c>bool[]</c> in the reference).
/// Counts and buckets are native blocks, the suffix array is a caller-supplied <see cref="Buf{T}"/>; all blocks come from a <see cref="NativeArena"/> and are released when
/// their level returns.
/// <para>
/// Safety: every symbol is range-checked once when counting (<c>[0, sigma)</c>, <see cref="ArgumentException"/> like the reference), after which all bucket indices
/// are in range; suffix-array entries are produced by the kernel itself. Level layout invariants (<c>2m &lt;= n</c>, name slot <c>1 + m + (p &gt;&gt; 1) &lt;= n</c>, head and
/// tail disjoint) are those of the reference and are exercised by the guard bands and the checked build.
/// </para>
/// </summary>
internal static unsafe class SaisKernel
{
    private const int Empty = -1;

    /// <summary>Builds the suffix array with the selected implementation.</summary>
    internal static int[] Build(SaisKernelKind kind, PreparedText prepared, SaisTrace? trace = null)
    {
        ArgumentNullException.ThrowIfNull(prepared);
        if (kind == SaisKernelKind.Pointer) return Build(prepared, trace);
        return prepared.Width == SymbolWidth.Byte
            ? SaisReference.Build(prepared.ByteReader(), prepared.AlphabetSize, trace)
            : SaisReference.Build(prepared.UShortReader(), prepared.AlphabetSize, trace);
    }

    /// <summary>Builds the suffix array (n + 1 layout) of a prepared text with its own workspace.</summary>
    internal static int[] Build(PreparedText prepared, SaisTrace? trace = null)
    {
        using var arena = new NativeArena();
        return Build(prepared, arena, trace);
    }

    /// <summary>Builds the suffix array with the workspace taken from <paramref name="arena"/> (the tests pass a guarded arena and verify its canaries).</summary>
    internal static int[] Build(PreparedText prepared, NativeArena arena, SaisTrace? trace = null)
    {
        ArgumentNullException.ThrowIfNull(prepared);
        ArgumentNullException.ThrowIfNull(arena);
        int n = prepared.Length;
        EsaLimits.EnsureTextLength(n);
        int[] sa = GC.AllocateUninitializedArray<int>(n + 1);
        fixed (int* pSa = sa)
        {
            var saBuf = new Buf<int>(pSa, n + 1, "sa");
            if (prepared.Width == SymbolWidth.Byte)
            {
                fixed (byte* pText = prepared.Bytes) Build(new Buf<byte>(pText, n, "text"), n, prepared.AlphabetSize, saBuf, arena, default(NoProbe), trace);
            }
            else
            {
                fixed (ushort* pText = prepared.UShorts) Build(new Buf<ushort>(pText, n, "text"), n, prepared.AlphabetSize, saBuf, arena, default(NoProbe), trace);
            }
        }
        return sa;
    }

    /// <summary>The kernel on accessor-wrapped byte symbols over <c>[0, sigma)</c>: <paramref name="sa"/> must have n + 1 entries (any content: all are written).</summary>
    internal static void Build<TProbe>(Buf<byte> text, int n, int sigma, Buf<int> sa, NativeArena arena, TProbe probe, SaisTrace? trace = null)
        where TProbe : struct, ISaisProbe
    {
        Prepare(n, sigma, arena);
        Sais<BytePtrText, TProbe>(new BytePtrText(text, n), n, sigma, sa, arena, probe, trace, 0);
    }

    /// <summary>The kernel on accessor-wrapped 16-bit symbols over <c>[0, sigma)</c>.</summary>
    internal static void Build<TProbe>(Buf<ushort> text, int n, int sigma, Buf<int> sa, NativeArena arena, TProbe probe, SaisTrace? trace = null)
        where TProbe : struct, ISaisProbe
    {
        Prepare(n, sigma, arena);
        Sais<UShortPtrText, TProbe>(new UShortPtrText(text, n), n, sigma, sa, arena, probe, trace, 0);
    }

    private static void Prepare(int n, int sigma, NativeArena arena)
    {
        ArgumentNullException.ThrowIfNull(arena);
        EsaLimits.EnsureTextLength(n);
        ArgumentOutOfRangeException.ThrowIfNegative(sigma);
    }

    private static void Sais<TText, TProbe>(TText t, int n, int sigma, Buf<int> sa, NativeArena arena, TProbe probe, SaisTrace? trace, int level)
        where TText : struct, IPtrText
        where TProbe : struct, ISaisProbe
    {
        sa[0] = n;
        if (n == 0)
        {
            trace?.Levels.Add(new SaisLevel(0, sigma, 0));
            return;
        }

        int mark = arena.Mark();
        long tick = SaisTrace.Begin(trace);
        Buf<int> counts = arena.Allocate<int>(sigma, zero: true, "counts");
        CountSymbols(t, n, sigma, counts);
        trace?.End(level, SaisPhase.Count, tick);
        if (TProbe.Enabled) probe.Counts(level, counts.AsSpan(sigma));

        tick = SaisTrace.Begin(trace);
        if (typeof(TText) == typeof(BytePtrText) && sigma < PackedBytes.Limit)
        {
            Buf<byte> packed = arena.Allocate<byte>(n + 1, zero: false, "packed");
            int m = ClassifyTypes(t, n, new PackedSink(packed));
            trace?.End(level, SaisPhase.Classify, tick);
            Level(new PackedBytes(packed), n, sigma, m, counts, sa, arena, probe, trace, level);
        }
        else
        {
            Buf<byte> types = arena.Allocate<byte>(n + 1, zero: false, "types");
            int m = ClassifyTypes(t, n, new SplitSink(types));
            trace?.End(level, SaisPhase.Classify, tick);
            Level(new SplitTyped<TText>(t, types), n, sigma, m, counts, sa, arena, probe, trace, level);
        }
        arena.ReleaseTo(mark);
    }

    /// <summary>One level after classification: stage 1 sorts the LMS substrings, names them and solves the reduced problem, stage 2 induces the final order.</summary>
    private static void Level<TTyped, TProbe>(TTyped tt, int n, int sigma, int m, Buf<int> counts, Buf<int> sa, NativeArena arena, TProbe probe, SaisTrace? trace, int level)
        where TTyped : struct, ITypedText
        where TProbe : struct, ISaisProbe
    {
        trace?.Levels.Add(new SaisLevel(n, sigma, m));
        Buf<int> bkt = arena.Allocate<int>(sigma, zero: false, "bkt");

        // Stage 1: sort the LMS substrings by one induction round, name them, solve the reduced problem.
        long tick = SaisTrace.Begin(trace);
        sa.Slice(1).Fill(Empty, n);
        SetBucketTails(counts, bkt, sigma);
        if (TProbe.Enabled) probe.Tails(level, bkt.AsSpan(sigma));
        PlaceLmsSuffixes(tt, bkt, sa, n);
        trace?.End(level, SaisPhase.Round1Place, tick);
        if (TProbe.Enabled) probe.SaState(level, SaisPhase.Round1Place, sa.AsSpan(n + 1));
        tick = SaisTrace.Begin(trace);
        SetBucketHeads(counts, bkt, sigma);
        if (TProbe.Enabled) probe.Heads(level, bkt.AsSpan(sigma));
        SaisInduction.InduceL(tt, sa, bkt, n);
        trace?.End(level, SaisPhase.Round1InduceL, tick);
        if (TProbe.Enabled) probe.SaState(level, SaisPhase.Round1InduceL, sa.AsSpan(n + 1));
        tick = SaisTrace.Begin(trace);
        SetBucketTails(counts, bkt, sigma);
        SaisInduction.InduceS(tt, sa, bkt, n);
        trace?.End(level, SaisPhase.Round1InduceS, tick);
        if (TProbe.Enabled) probe.SaState(level, SaisPhase.Round1InduceS, sa.AsSpan(n + 1));

        tick = SaisTrace.Begin(trace);
        CompactSortedLms(tt, sa, n, m);
        trace?.End(level, SaisPhase.CompactLms, tick);
        if (TProbe.Enabled)
        {
            probe.LmsSorted(level, sa.Slice(1).AsSpan(m));
            probe.SaState(level, SaisPhase.CompactLms, sa.AsSpan(n + 1));
        }
        tick = SaisTrace.Begin(trace);
        int names = SaisInduction.NameLmsSubstrings(tt, sa, n, m);
        trace?.End(level, SaisPhase.Name, tick);
        if (TProbe.Enabled) probe.SaState(level, SaisPhase.Name, sa.AsSpan(n + 1));
        tick = SaisTrace.Begin(trace);
        CompactNames(sa, n, m);
        trace?.End(level, SaisPhase.CompactNames, tick);
        if (TProbe.Enabled) probe.Reduced(level, names, sa.Slice(n + 1 - m).AsSpan(m));
        SolveReduced(tt, sa, n, m, names, arena, probe, trace, level);

        // Stage 2: place the LMS suffixes in their final order and induce all others from them.
        tick = SaisTrace.Begin(trace);
        sa.Slice(m + 1).Fill(Empty, n - m);
        SetBucketTails(counts, bkt, sigma);
        PlaceSortedLms(tt, bkt, sa, m);
        trace?.End(level, SaisPhase.Round2Place, tick);
        if (TProbe.Enabled) probe.SaState(level, SaisPhase.Round2Place, sa.AsSpan(n + 1));
        tick = SaisTrace.Begin(trace);
        SetBucketHeads(counts, bkt, sigma);
        SaisInduction.InduceL(tt, sa, bkt, n);
        trace?.End(level, SaisPhase.Round2InduceL, tick);
        if (TProbe.Enabled) probe.SaState(level, SaisPhase.Round2InduceL, sa.AsSpan(n + 1));
        tick = SaisTrace.Begin(trace);
        SetBucketTails(counts, bkt, sigma);
        SaisInduction.InduceS(tt, sa, bkt, n);
        trace?.End(level, SaisPhase.Round2InduceS, tick);
        if (TProbe.Enabled) probe.SaState(level, SaisPhase.Round2InduceS, sa.AsSpan(n + 1));
    }

    /// <summary><c>counts[c]</c> = number of occurrences of symbol c; every symbol is range-checked, which makes all later bucket indices safe.</summary>
    private static void CountSymbols<TText>(TText t, int n, int sigma, Buf<int> counts)
        where TText : struct, IPtrText
    {
        for (int i = 0; i < n; i++)
        {
            int c = t[i];
            if ((uint)c >= (uint)sigma)
                throw new ArgumentException($"symbol {c} at index {i} is outside [0, {sigma})", nameof(t));
            counts[c]++;
        }
    }

    /// <summary>Type of every position (1 = S-type; position n, the virtual sentinel, is S-type, position n - 1 is L-type) written through <paramref name="sink"/>, and the number of LMS positions in <c>[1, n - 1]</c>, in one right-to-left pass.</summary>
    private static int ClassifyTypes<TText, TSink>(TText t, int n, TSink sink)
        where TText : struct, IPtrText
        where TSink : struct, ITypeSink
    {
        sink.SetSentinel(n);
        int next = t[n - 1];
        sink.Set(n - 1, next, 0);
        int m = 0;
        int nextIsS = 0;
        for (int i = n - 2; i >= 0; i--)
        {
            int c = t[i];
            int d = c - next; // symbols are below 2^30: no overflow
            int smaller = (d >> 31) & 1; // c < next
            int equal = ((d | -d) >> 31) + 1; // c == next
            int isS = smaller | (equal & nextIsS); // data-dependent branches here are coin flips on DNA: arithmetic instead
            sink.Set(i, c, isS);
            m += nextIsS & (isS ^ 1); // position i + 1 is S-type with an L-type predecessor: an LMS position
            nextIsS = isS;
            next = c;
        }
        return m;
    }

    private static void SetBucketHeads(Buf<int> counts, Buf<int> bkt, int sigma)
    {
        int sum = 1; // index 0 belongs to the sentinel suffix
        for (int c = 0; c < sigma; c++)
        {
            bkt[c] = sum;
            sum += counts[c];
        }
    }

    private static void SetBucketTails(Buf<int> counts, Buf<int> bkt, int sigma)
    {
        int sum = 1;
        for (int c = 0; c < sigma; c++)
        {
            sum += counts[c];
            bkt[c] = sum; // exclusive end
        }
    }

    /// <summary>
    /// Stage 1: LMS positions into the ends of their buckets (any order), scanning right to left. Branch-free: "is an LMS position" is a coin flip on DNA; a position that is
    /// not LMS stores to the dummy index 0 (restored at the end, <c>SA[0] = n</c>) and leaves the bucket cursor alone.
    /// </summary>
    private static void PlaceLmsSuffixes<TTyped>(TTyped tt, Buf<int> bkt, Buf<int> sa, int n)
        where TTyped : struct, ITypedText
    {
        tt.Get(n - 1, out int c, out int isS);
        for (int i = n - 1; i >= 1; i--)
        {
            tt.Get(i - 1, out int previousSymbol, out int previous);
            int lms = isS & (previous ^ 1); // i is S-type with an L-type predecessor
            int slot = bkt[c] - lms;
            sa[slot & -lms] = i;
            bkt[c] = slot;
            c = previousSymbol;
            isS = previous;
        }
        sa[0] = n;
    }

    /// <summary>Moves the sorted LMS positions to <c>SA[1..m]</c> (stable, in place: the write index never passes the read index). Branch-free; the entries behind the last LMS position are left undefined.</summary>
    private static void CompactSortedLms<TTyped>(TTyped tt, Buf<int> sa, int n, int m)
        where TTyped : struct, ITypedText
    {
        int k = 0;
        for (int i = 1; i <= n; i++)
        {
            int p = sa[i];
            int positive = (int)((uint)-p >> 31);
            int pc = p & -positive; // a valid position even for an empty slot
            int before = (pc - 1) & -positive; // its predecessor (0 for p <= 0)
            int lms = positive & tt.IsS(pc) & (tt.IsS(before) ^ 1);
            sa[1 + k] = p; // the write index is at most i: behind the read cursor; a non-LMS entry is overwritten by the next LMS position
            k += lms;
        }
        if (k != m)
            throw new InvalidOperationException($"internal error: {k} LMS positions after induction, expected {m}");
    }

    /// <summary>
    /// Gathers the names into the tail <c>SA[n + 1 - m .. n]</c> in text order (the reduced string s1). Branch-free: every entry is stored at the write cursor and the
    /// cursor moves down only for a name, so the entry at <c>SA[n - m]</c> is left undefined (it is not part of the tail, and nothing reads it).
    /// </summary>
    private static void CompactNames(Buf<int> sa, int n, int m)
    {
        int j = n;
        for (int i = n; i > m; i--)
        {
            int v = sa[i];
            sa[j] = v; // j >= i: the entry has been read already
            j -= (int)((uint)~v >> 31); // 1 for v >= 0
        }
    }

    /// <summary>
    /// Solves the reduced problem and replaces <c>SA[1..m]</c> by the LMS positions in sorted suffix order. The reduced text sits in the tail, its suffix
    /// array (n + 1 layout, so m + 1 entries) in the head; <c>SA[0]</c> is restored to n afterwards.
    /// </summary>
    private static void SolveReduced<TTyped, TProbe>(TTyped tt, Buf<int> sa, int n, int m, int names, NativeArena arena, TProbe probe, SaisTrace? trace, int level)
        where TTyped : struct, ITypedText
        where TProbe : struct, ISaisProbe
    {
        long tick = SaisTrace.Begin(trace);
        Buf<int> sa1 = sa.Slice(0, m + 1);
        Buf<int> s1 = sa.Slice(n + 1 - m, m);
        if (names < m)
        {
            trace?.End(level, SaisPhase.Reduced, tick); // the recursion is accounted by its own levels
            Sais<IntPtrText, TProbe>(new IntPtrText(s1, m), m, names, sa1, arena, probe, trace, level + 1);
            tick = SaisTrace.Begin(trace);
        }
        else
        {
            sa1[0] = m;
            for (int i = 0; i < m; i++) sa1[1 + s1[i]] = i; // all names distinct: the name is the rank
        }

        if (TProbe.Enabled) probe.SaState(level, SaisPhase.Reduced, sa1.AsSpan(m + 1)); // the suffix array of the reduced string (n + 1 layout)

        // reuse the tail for the LMS positions in text order, then translate ranks of s1 into text positions
        int tailBase = n + 1 - m;
        int idx = m;
        int isS = tt.IsS(n - 1);
        for (int i = n - 1; i >= 1; i--)
        {
            int previous = tt.IsS(i - 1);
            int lms = isS & (previous ^ 1);
            idx -= lms;
            sa[(tailBase + idx) & -lms] = i; // not an LMS position: the dummy index 0, restored below
            isS = previous;
        }
        if (TProbe.Enabled) probe.LmsPositions(level, sa.Slice(tailBase).AsSpan(m));
        for (int i = 1; i <= m; i++) sa[i] = sa[tailBase + sa[i]];
        sa[0] = n;
        trace?.End(level, SaisPhase.Reduced, tick);
    }

    /// <summary>Stage 2: the sorted LMS suffixes go to the ends of their buckets, largest first.</summary>
    private static void PlaceSortedLms<TTyped>(TTyped tt, Buf<int> bkt, Buf<int> sa, int m)
        where TTyped : struct, ITypedText
    {
        for (int r = m; r >= 1; r--)
        {
            int p = sa[r];
            sa[r] = Empty;
            int c = tt.Symbol(p);
            int slot = bkt[c] - 1;
            sa[slot] = p; // the target index is >= r, so no unread LMS entry is overwritten
            bkt[c] = slot;
        }
    }
}
