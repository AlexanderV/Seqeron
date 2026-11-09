using SuffixTree.Esa;
using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Build;

/// <summary>
/// Sorts the suffixes of one leaf bucket in RAM (PH3 of ADR-004): key = the first 64 bits of packed symbols, LSD radix sort of (key, position) pairs that skips digits on which
/// all keys agree, then every run of equal keys is ordered by the suffix comparator (the DCS comparator when ranks exist). The scratch lives in four arena blocks of K entries.
/// </summary>
internal sealed unsafe class LeafSorter
{
    private readonly Buf<ulong> _keys;
    private readonly Buf<ulong> _keys2;
    private readonly Buf<uint> _pos;
    private readonly Buf<uint> _pos2;
    private readonly int _capacity;
    private readonly TieSort _tieSort;

    internal LeafSorter(BudgetedArena arena, int capacity, TieSort tieSort)
    {
        _capacity = capacity;
        _tieSort = tieSort;
        _keys = arena.Allocate<ulong>(capacity, zero: false, "leaf.keys");
        _keys2 = arena.Allocate<ulong>(capacity, zero: false, "leaf.keys2");
        _pos = arena.Allocate<uint>(capacity, zero: false, "leaf.pos");
        _pos2 = arena.Allocate<uint>(capacity, zero: false, "leaf.pos2");
    }

    /// <summary>Sorts the bucket and writes its suffix array slice at <paramref name="saOffset"/>.</summary>
    internal void Sort(IPosSource bucket, long saOffset, SuffixOrder order, ISaSink sink)
    {
        int m = checked((int)bucket.Count);
        if (m > _capacity) throw new InvalidOperationException($"leaf of {m} positions exceeds the capacity {_capacity}");
        Span<uint> pos = _pos.AsSpan(m);
        bucket.Read(0, pos);
        Span<ulong> keys = _keys.AsSpan(m);
        var text = order.Text;
        for (int i = 0; i < m; i++) keys[i] = text.Key(pos[i]);

        if (m <= 24) InsertionSort(keys, pos);
        else RadixSort(m);

        keys = _keys.AsSpan(m);
        pos = _pos.AsSpan(m);
        WorkCounters c = order.Counters;
        c.LeafPositions += m;
        int s = 0;
        while (s < m)
        {
            int e = s + 1;
            while (e < m && keys[e] == keys[s]) e++;
            if (e - s >= 2)
            {
                c.TiePositions += e - s;
                if (e - s == 2)
                {
                    if (order.CompareTied(pos[s], pos[s + 1]) > 0) (pos[s], pos[s + 1]) = (pos[s + 1], pos[s]);
                }
                else if (_tieSort == TieSort.Merge)
                {
                    MergeSort(pos.Slice(s, e - s), _pos2.AsSpan(e - s), order);
                }
                else
                {
                    pos.Slice(s, e - s).Sort(new TiedComparer(order));
                }
            }
            s = e;
        }
        sink.Write(saOffset, pos);
    }

    /// <summary>
    /// Stable merge sort of the positions of a run of equal keys by the tied comparator. Chosen over the introsort of <see cref="MemoryExtensions"/> because the cost here is the number of
    /// comparator calls (each may compare hundreds of symbols or touch two sample ranks): a merge sort needs at most n - 1 per merge, only n / 2 when the run is in reverse order (the suffixes of a^n and
    /// of long tandem repeats arrive in text order, which is the reverse of their order), and a merge is skipped when the halves are already in order.
    /// </summary>
    internal static void MergeSort(Span<uint> a, Span<uint> tmp, SuffixOrder order)
    {
        int n = a.Length;
        if (n <= 6)
        {
            for (int i = 1; i < n; i++)
            {
                uint x = a[i];
                int lo = 0, hi = i; // binary insertion: first index whose element is greater than x
                while (lo < hi)
                {
                    int mid = (lo + hi) >> 1;
                    if (order.CompareTied(a[mid], x) <= 0) lo = mid + 1;
                    else hi = mid;
                }
                for (int j = i; j > lo; j--) a[j] = a[j - 1];
                a[lo] = x;
            }
            return;
        }
        int half = n >> 1;
        MergeSort(a[..half], tmp[..half], order);
        MergeSort(a[half..], tmp[half..], order);
        if (order.CompareTied(a[half - 1], a[half]) <= 0) return;
        a[..half].CopyTo(tmp);
        int p = 0, q = half, k = 0;
        while (p < half && q < n)
        {
            if (order.CompareTied(a[q], tmp[p]) < 0) a[k++] = a[q++];
            else a[k++] = tmp[p++];
        }
        while (p < half) a[k++] = tmp[p++];
    }

    private static void InsertionSort(Span<ulong> keys, Span<uint> pos)
    {
        for (int i = 1; i < keys.Length; i++)
        {
            ulong k = keys[i];
            uint p = pos[i];
            int j = i - 1;
            while (j >= 0 && keys[j] > k)
            {
                keys[j + 1] = keys[j];
                pos[j + 1] = pos[j];
                j--;
            }
            keys[j + 1] = k;
            pos[j + 1] = p;
        }
    }

    private void RadixSort(int m)
    {
        Span<ulong> k1 = _keys.AsSpan(m);
        Span<ulong> k2 = _keys2.AsSpan(m);
        Span<uint> p1 = _pos.AsSpan(m);
        Span<uint> p2 = _pos2.AsSpan(m);
        Span<int> hist = stackalloc int[8 * 256];
        hist.Clear();
        foreach (ulong k in k1)
        {
            hist[(int)(k & 0xFF)]++;
            hist[256 + (int)((k >> 8) & 0xFF)]++;
            hist[512 + (int)((k >> 16) & 0xFF)]++;
            hist[768 + (int)((k >> 24) & 0xFF)]++;
            hist[1024 + (int)((k >> 32) & 0xFF)]++;
            hist[1280 + (int)((k >> 40) & 0xFF)]++;
            hist[1536 + (int)((k >> 48) & 0xFF)]++;
            hist[1792 + (int)(k >> 56)]++;
        }
        ulong first = k1[0];
        bool swapped = false;
        for (int d = 0; d < 8; d++)
        {
            Span<int> h = hist.Slice(d * 256, 256);
            int shift = d * 8;
            if (h[(int)((first >> shift) & 0xFF)] == m) continue; // all keys have the same digit
            int sum = 0;
            for (int b = 0; b < 256; b++)
            {
                int c = h[b];
                h[b] = sum;
                sum += c;
            }
            Span<ulong> sk = swapped ? k2 : k1;
            Span<ulong> dk = swapped ? k1 : k2;
            Span<uint> sp = swapped ? p2 : p1;
            Span<uint> dp = swapped ? p1 : p2;
            for (int i = 0; i < m; i++)
            {
                ulong k = sk[i];
                int at = h[(int)((k >> shift) & 0xFF)]++;
                dk[at] = k;
                dp[at] = sp[i];
            }
            swapped = !swapped;
        }
        if (swapped)
        {
            k2.CopyTo(k1);
            p2.CopyTo(p1);
        }
    }
}
