using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;

namespace SuffixTree.Esa;

/// <summary>
/// Iteration 1o.1c: the induction passes (L and S) and the naming of the LMS substrings of <see cref="SaisKernel"/> as pointer loops over <see cref="Buf{T}"/>, the same
/// algorithm and the same result as the safe passes of <see cref="SaisReference{TIdx}"/> (the in-situ test oracle <c>SaisSafePasses</c> in the test project is a verbatim copy;
/// every pass of every level is replayed on it, state for state). The passes are generic over <see cref="ITypedText"/>: for byte texts of up to 127 symbols
/// the symbol and the type of a position are ONE byte (<see cref="PackedBytes"/>), so each induction step makes one random memory access instead of two.
/// <para>
/// The induction passes have NO data-dependent branch. With the safe code, and with a straight translation of it, each iteration waits for the cache miss on
/// <c>types[p - 1]</c> to resolve the branch "is this predecessor of the right type" (a coin flip on DNA), and the mispredict throws away everything that the core had
/// started on the following iterations, so the misses are served one at a time. Here the decision is arithmetic: <c>take</c> is 0 or 1, the target index is
/// <c>slot &amp; -take</c> (the induced entry, or the harmless dummy index 0) and the bucket cursor advances by <c>take</c>. Nothing the next iteration needs depends on
/// the outcome except through the cursor, so the out-of-order core overlaps the misses of several iterations. The price is a dummy store to <c>SA[0]</c> per skipped entry,
/// which is undone when the pass ends (<c>SA[0] = n</c> is the one entry the L pass reads first and the S pass never reads).
/// </para>
/// <para>
/// Safety: all indices are in range for a valid state. Symbols are range-checked once at counting time, so every bucket index <c>c</c> is in <c>[0, sigma)</c>; entries
/// of the suffix array are produced by the kernel itself (positions in <c>[0, n]</c> or -1), the clamped predecessor index <c>jc</c> is in <c>[0, n - 1]</c>, the
/// bucket cursors stay inside their buckets by construction (counts equal the number of induced entries), and the dummy index is 0. The checked build and the guard bands
/// of the tests verify this on every corpus.
/// </para>
/// </summary>
internal static unsafe class SaisInduction
{
    /// <summary>Induces the L-type suffixes from the sorted entries, left to right; <paramref name="bkt"/> holds the bucket heads on entry.</summary>
    internal static void InduceL<TTyped>(TTyped tt, Buf<int> sa, Buf<int> bkt, int n)
        where TTyped : struct, ITypedText
    {
        for (int i = 0; i <= n; i++)
        {
            int ahead = sa[Math.Min(i + SaisPrefetch.Induce, n)]; // the entry that will be scanned SaisPrefetch.Induce iterations from now (a prefetch of a stale entry is wasted, never wrong)
            tt.Prefetch((ahead - 1) & -(int)((uint)-ahead >> 31));
            int p = sa[i];
            int positive = (int)((uint)-p >> 31); // 1 for p > 0; empty slots (-1) and suffix 0 induce nothing
            int j = p - 1;
            int jc = j & -positive; // a valid index even when nothing is induced
            tt.Get(jc, out int c, out int isS);
            int take = positive & (isS ^ 1); // the predecessor is L-type
            int slot = bkt[c];
            sa[slot & -take] = j;
            bkt[c] = slot + take;
        }
        sa[0] = n; // undo the dummy stores
    }

    /// <summary>Induces the S-type suffixes from the sorted entries, right to left; <paramref name="bkt"/> holds the exclusive bucket ends on entry.</summary>
    internal static void InduceS<TTyped>(TTyped tt, Buf<int> sa, Buf<int> bkt, int n)
        where TTyped : struct, ITypedText
    {
        for (int i = n; i >= 1; i--)
        {
            int ahead = sa[Math.Max(i - SaisPrefetch.Induce, 1)];
            tt.Prefetch((ahead - 1) & -(int)((uint)-ahead >> 31));
            int p = sa[i];
            int positive = (int)((uint)-p >> 31);
            int j = p - 1;
            int jc = j & -positive;
            tt.Get(jc, out int c, out int isS);
            int take = positive & isS; // the predecessor is S-type
            int slot = bkt[c] - take;
            sa[slot & -take] = j;
            bkt[c] = slot;
        }
        sa[0] = n;
    }

    /// <summary>
    /// Experiment 1o.4: the L pass with the (at most four) bucket cursors in one SIMD register instead of the <c>bkt</c> array, for packed byte texts with sigma &lt;= 4.
    /// </summary>
    internal static void InduceL<TTyped>(TTyped tt, Buf<int> sa, Buf<int> bkt, int n, int sigma)
        where TTyped : struct, ITypedText
    {
        if (Ind4 && typeof(TTyped) == typeof(PackedBytes) && sigma <= 4) InduceL4(Unsafe.As<TTyped, PackedBytes>(ref tt), sa, bkt, sigma, n);
        else InduceL(tt, sa, bkt, n);
    }

    /// <summary>Experiment 1o.4: the S pass, see the L overload.</summary>
    internal static void InduceS<TTyped>(TTyped tt, Buf<int> sa, Buf<int> bkt, int n, int sigma)
        where TTyped : struct, ITypedText
    {
        if (Ind4 && typeof(TTyped) == typeof(PackedBytes) && sigma <= 4) InduceS4(Unsafe.As<TTyped, PackedBytes>(ref tt), sa, bkt, sigma, n);
        else InduceS(tt, sa, bkt, n);
    }

    internal static readonly bool Ind4 = Environment.GetEnvironmentVariable("ESA_IND4") == "1"; // experiment switch

    private static void InduceL4(PackedBytes tt, Buf<int> sa, Buf<int> bkt, int sigma, int n)
    {
        Vector128<int> cursors = LoadCursors(bkt, sigma);
        Vector128<int> lanes = Vector128.Create(0, 1, 2, 3);
        for (int i = 0; i <= n; i++)
        {
            int ahead = sa[Math.Min(i + SaisPrefetch.Induce, n)];
            tt.Prefetch((ahead - 1) & -(int)((uint)-ahead >> 31));
            int p = sa[i];
            int positive = (int)((uint)-p >> 31);
            int j = p - 1;
            int jc = j & -positive;
            tt.Get(jc, out int c, out int isS);
            int take = positive & (isS ^ 1);
            Vector128<int> which = Vector128.Create(c);
            int slot = Vector128.Shuffle(cursors, which).ToScalar();
            sa[slot & -take] = j;
            cursors += Vector128.Equals(which, lanes) & Vector128.Create(take);
        }
        StoreCursors(bkt, sigma, cursors);
        sa[0] = n;
    }

    private static void InduceS4(PackedBytes tt, Buf<int> sa, Buf<int> bkt, int sigma, int n)
    {
        Vector128<int> cursors = LoadCursors(bkt, sigma);
        Vector128<int> lanes = Vector128.Create(0, 1, 2, 3);
        for (int i = n; i >= 1; i--)
        {
            int ahead = sa[Math.Max(i - SaisPrefetch.Induce, 1)];
            tt.Prefetch((ahead - 1) & -(int)((uint)-ahead >> 31));
            int p = sa[i];
            int positive = (int)((uint)-p >> 31);
            int j = p - 1;
            int jc = j & -positive;
            tt.Get(jc, out int c, out int isS);
            int take = positive & isS;
            Vector128<int> which = Vector128.Create(c);
            cursors -= Vector128.Equals(which, lanes) & Vector128.Create(take);
            int slot = Vector128.Shuffle(cursors, which).ToScalar();
            sa[slot & -take] = j;
        }
        StoreCursors(bkt, sigma, cursors);
        sa[0] = n;
    }

    private static Vector128<int> LoadCursors(Buf<int> bkt, int sigma)
    {
        int b1 = sigma > 1 ? bkt[1] : 0, b2 = sigma > 2 ? bkt[2] : 0, b3 = sigma > 3 ? bkt[3] : 0;
        return Vector128.Create(sigma > 0 ? bkt[0] : 0, b1, b2, b3);
    }

    private static void StoreCursors(Buf<int> bkt, int sigma, Vector128<int> cursors)
    {
        for (int c = 0; c < sigma; c++) bkt[c] = cursors.GetElement(c);
    }

    /// <summary>
    /// Names the sorted LMS substrings <c>SA[1..m]</c> (equal substrings share a name) and stores the name of the LMS position p at <c>SA[1 + m + (p &gt;&gt; 1)]</c>, all
    /// other entries of <c>SA[m + 1..n]</c> become -1. Returns the number of distinct names.
    /// </summary>
    internal static int NameLmsSubstrings<TTyped>(TTyped tt, Buf<int> sa, int n, int m)
        where TTyped : struct, ITypedText
    {
        sa.Slice(m + 1).Fill(-1, n - m);
        int name = -1;
        int previous = -1;
        for (int i = 1; i <= m; i++)
        {
            int ahead = sa[Math.Min(i + SaisPrefetch.Name, m)]; // the sorted LMS position that is compared SaisPrefetch.Name iterations from now, and the slot that gets its name
            tt.Prefetch(ahead);
            sa.Prefetch(1 + m + (ahead >> 1));
            int pos = sa[i];
            if (previous < 0 || !LmsSubstringsEqual(tt, n, previous, pos)) name++;
            sa[1 + m + (pos >> 1)] = name;
            previous = pos;
        }
        return name + 1;
    }

    /// <summary>
    /// Equality of the LMS substrings that start at the distinct LMS positions a and b: same length, same symbols, same types; the substring that reaches the sentinel is unique.
    /// A branch per symbol on purpose: the first mismatch comes after a few symbols, a word compare would make the next load depend on this one.
    /// </summary>
    internal static bool LmsSubstringsEqual<TTyped>(TTyped tt, int n, int a, int b)
        where TTyped : struct, ITypedText
    {
        int x = a, y = b;
        while (x < n && y < n)
        {
            if (!tt.Same(x, y)) return false;
            if (x != a && tt.IsS(x) != 0 && tt.IsS(x - 1) == 0) return true; // the end of the substring (an LMS position behind the first)
            x++;
            y++;
        }
        return false; // one of them reaches the sentinel
    }
}
