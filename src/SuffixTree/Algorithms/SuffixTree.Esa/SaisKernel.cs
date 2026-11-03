namespace SuffixTree.Esa;

/// <summary>Which suffix-array implementation a harness runs: the safe permanent oracle or the pointer kernel (iteration 1o.1b).</summary>
internal enum SaisKernelKind
{
    /// <summary><see cref="SaisReference"/>: safe code, generic over the index type, never optimized.</summary>
    Reference,

    /// <summary><see cref="SaisKernel"/>: <c>int</c>-only kernel; the linear passes are pointer loops behind <see cref="Buf{T}"/>.</summary>
    Pointer,
}

/// <summary>
/// Iteration 1o.1b (REVIEW-S1a section 6): SA-IS with the LINEAR passes as pointer loops, the same algorithm, buffer layout and result as
/// <see cref="SaisReference{TIdx}"/> (bit for bit, equal recursion trace). Pointer passes: symbol counting, type classification fused with the LMS count, bucket
/// setup, placement of the LMS suffixes into their buckets (both rounds), compaction of the sorted LMS positions and of the names, the reduced-problem
/// bookkeeping (direct placement, LMS positions in text order, translation of the reduced ranks). Types are one byte per position in a native block
/// (<c>bool[]</c> in the reference), counts and buckets are native blocks, the suffix array is a caller-supplied <see cref="Buf{T}"/>; all blocks come
/// from a <see cref="NativeArena"/> and are released when their level returns. The induction passes and the naming stay the reference algorithm
/// (<see cref="SaisSafePasses"/>, safe code) and are replaced in iteration 1o.1c.
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
        Sais<BytePtrText, ByteSymbols, TProbe>(new BytePtrText(text, n), new ByteSymbols(text.AsSpan(n)), n, sigma, sa, arena, probe, trace, 0);
    }

    /// <summary>The kernel on accessor-wrapped 16-bit symbols over <c>[0, sigma)</c>.</summary>
    internal static void Build<TProbe>(Buf<ushort> text, int n, int sigma, Buf<int> sa, NativeArena arena, TProbe probe, SaisTrace? trace = null)
        where TProbe : struct, ISaisProbe
    {
        Prepare(n, sigma, arena);
        Sais<UShortPtrText, UShortSymbols, TProbe>(new UShortPtrText(text, n), new UShortSymbols(text.AsSpan(n)), n, sigma, sa, arena, probe, trace, 0);
    }

    private static void Prepare(int n, int sigma, NativeArena arena)
    {
        ArgumentNullException.ThrowIfNull(arena);
        EsaLimits.EnsureTextLength(n);
        ArgumentOutOfRangeException.ThrowIfNegative(sigma);
    }

    private static void Sais<TText, TReader, TProbe>(TText t, TReader r, int n, int sigma, Buf<int> sa, NativeArena arena, TProbe probe, SaisTrace? trace, int level)
        where TText : struct, IPtrText
        where TReader : ISymbolReader, allows ref struct
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
        Buf<byte> types = arena.Allocate<byte>(n + 1, zero: false, "types");
        int m = ClassifyTypes(t, n, types);
        trace?.End(level, SaisPhase.Classify, tick);
        trace?.Levels.Add(new SaisLevel(n, sigma, m));
        Buf<int> bkt = arena.Allocate<int>(sigma, zero: false, "bkt");
        Span<int> bucketSpan = bkt.AsSpan(sigma);
        Span<int> saSpan = sa.AsSpan(n + 1);
        ReadOnlySpan<byte> typeSpan = types.AsSpan(n + 1);

        // Stage 1: sort the LMS substrings by one induction round, name them, solve the reduced problem.
        tick = SaisTrace.Begin(trace);
        sa.Slice(1).Fill(Empty, n);
        SetBucketTails(counts, bkt, sigma);
        if (TProbe.Enabled) probe.Tails(level, bucketSpan);
        PlaceLmsSuffixes(t, types, bkt, sa, n);
        trace?.End(level, SaisPhase.Round1Place, tick);
        tick = SaisTrace.Begin(trace);
        SetBucketHeads(counts, bkt, sigma);
        if (TProbe.Enabled) probe.Heads(level, bucketSpan);
        SaisSafePasses.InduceL(r, saSpan, typeSpan, bucketSpan);
        trace?.End(level, SaisPhase.Round1InduceL, tick);
        tick = SaisTrace.Begin(trace);
        SetBucketTails(counts, bkt, sigma);
        SaisSafePasses.InduceS(r, saSpan, typeSpan, bucketSpan);
        trace?.End(level, SaisPhase.Round1InduceS, tick);

        tick = SaisTrace.Begin(trace);
        CompactSortedLms(types, sa, n, m);
        trace?.End(level, SaisPhase.CompactLms, tick);
        if (TProbe.Enabled) probe.LmsSorted(level, sa.Slice(1).AsSpan(m));
        tick = SaisTrace.Begin(trace);
        int names = SaisSafePasses.NameLmsSubstrings(r, saSpan, typeSpan, m);
        trace?.End(level, SaisPhase.Name, tick);
        tick = SaisTrace.Begin(trace);
        CompactNames(sa, n, m);
        trace?.End(level, SaisPhase.CompactNames, tick);
        if (TProbe.Enabled) probe.Reduced(level, names, sa.Slice(n + 1 - m).AsSpan(m));
        SolveReduced(types, sa, n, m, names, arena, probe, trace, level);

        // Stage 2: place the LMS suffixes in their final order and induce all others from them.
        tick = SaisTrace.Begin(trace);
        sa.Slice(m + 1).Fill(Empty, n - m);
        SetBucketTails(counts, bkt, sigma);
        PlaceSortedLms(t, bkt, sa, m);
        trace?.End(level, SaisPhase.Round2Place, tick);
        tick = SaisTrace.Begin(trace);
        SetBucketHeads(counts, bkt, sigma);
        SaisSafePasses.InduceL(r, saSpan, typeSpan, bucketSpan);
        trace?.End(level, SaisPhase.Round2InduceL, tick);
        tick = SaisTrace.Begin(trace);
        SetBucketTails(counts, bkt, sigma);
        SaisSafePasses.InduceS(r, saSpan, typeSpan, bucketSpan);
        trace?.End(level, SaisPhase.Round2InduceS, tick);
        arena.ReleaseTo(mark);
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

    /// <summary>Type of every position (1 = S-type; position n, the virtual sentinel, is S-type, position n - 1 is L-type) and the number of LMS positions in <c>[1, n - 1]</c>, in one right-to-left pass.</summary>
    private static int ClassifyTypes<TText>(TText t, int n, Buf<byte> types)
        where TText : struct, IPtrText
    {
        types[n] = 1;
        types[n - 1] = 0;
        int m = 0;
        int next = t[n - 1];
        byte nextIsS = 0;
        for (int i = n - 2; i >= 0; i--)
        {
            int c = t[i];
            byte isS = c < next || (c == next && nextIsS != 0) ? (byte)1 : (byte)0;
            types[i] = isS;
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

    /// <summary>Stage 1: LMS positions into the ends of their buckets (any order), scanning right to left.</summary>
    private static void PlaceLmsSuffixes<TText>(TText t, Buf<byte> types, Buf<int> bkt, Buf<int> sa, int n)
        where TText : struct, IPtrText
    {
        byte isS = types[n - 1];
        for (int i = n - 1; i >= 1; i--)
        {
            byte previous = types[i - 1];
            if (isS != 0 && previous == 0)
            {
                int c = t[i];
                int slot = bkt[c] - 1;
                sa[slot] = i;
                bkt[c] = slot;
            }
            isS = previous;
        }
    }

    /// <summary>Moves the sorted LMS positions to <c>SA[1..m]</c> (stable, in place: the write index never passes the read index).</summary>
    private static void CompactSortedLms(Buf<byte> types, Buf<int> sa, int n, int m)
    {
        int k = 0;
        for (int i = 1; i <= n; i++)
        {
            int p = sa[i];
            if (p > 0 && types[p] != 0 && types[p - 1] == 0) sa[1 + k++] = p;
        }
        if (k != m)
            throw new InvalidOperationException($"internal error: {k} LMS positions after induction, expected {m}");
    }

    /// <summary>Gathers the names into the tail <c>SA[n + 1 - m .. n]</c> in text order (the reduced string s1).</summary>
    private static void CompactNames(Buf<int> sa, int n, int m)
    {
        int j = n;
        for (int i = n; i > m; i--)
        {
            int v = sa[i];
            if (v >= 0) sa[j--] = v;
        }
    }

    /// <summary>
    /// Solves the reduced problem and replaces <c>SA[1..m]</c> by the LMS positions in sorted suffix order. The reduced text sits in the tail, its suffix
    /// array (n + 1 layout, so m + 1 entries) in the head; <c>SA[0]</c> is restored to n afterwards.
    /// </summary>
    private static void SolveReduced<TProbe>(Buf<byte> types, Buf<int> sa, int n, int m, int names, NativeArena arena, TProbe probe, SaisTrace? trace, int level)
        where TProbe : struct, ISaisProbe
    {
        long tick = SaisTrace.Begin(trace);
        Buf<int> sa1 = sa.Slice(0, m + 1);
        Buf<int> s1 = sa.Slice(n + 1 - m, m);
        if (names < m)
        {
            trace?.End(level, SaisPhase.Reduced, tick); // the recursion is accounted by its own levels
            Sais<IntPtrText, IntSymbols, TProbe>(new IntPtrText(s1, m), new IntSymbols(s1.AsSpan(m)), m, names, sa1, arena, probe, trace, level + 1);
            tick = SaisTrace.Begin(trace);
        }
        else
        {
            sa1[0] = m;
            for (int i = 0; i < m; i++) sa1[1 + s1[i]] = i; // all names distinct: the name is the rank
        }

        // reuse the tail for the LMS positions in text order, then translate ranks of s1 into text positions
        int tailBase = n + 1 - m;
        int idx = m;
        byte isS = types[n - 1];
        for (int i = n - 1; i >= 1; i--)
        {
            byte previous = types[i - 1];
            if (isS != 0 && previous == 0) sa[tailBase + --idx] = i;
            isS = previous;
        }
        if (TProbe.Enabled) probe.LmsPositions(level, sa.Slice(tailBase).AsSpan(m));
        for (int i = 1; i <= m; i++) sa[i] = sa[tailBase + sa[i]];
        sa[0] = n;
        trace?.End(level, SaisPhase.Reduced, tick);
    }

    /// <summary>Stage 2: the sorted LMS suffixes go to the ends of their buckets, largest first.</summary>
    private static void PlaceSortedLms<TText>(TText t, Buf<int> bkt, Buf<int> sa, int m)
        where TText : struct, IPtrText
    {
        for (int r = m; r >= 1; r--)
        {
            int p = sa[r];
            sa[r] = Empty;
            int c = t[p];
            int slot = bkt[c] - 1;
            sa[slot] = p; // the target index is >= r, so no unread LMS entry is overwritten
            bkt[c] = slot;
        }
    }
}
