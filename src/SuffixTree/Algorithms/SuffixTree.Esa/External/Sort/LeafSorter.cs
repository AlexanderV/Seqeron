namespace SuffixTree.Esa.External;

/// <summary>
/// Sorts the suffixes of one leaf bucket in RAM (PH3 of ADR-004): key = the first 64 bits of packed symbols, LSD radix sort of (key, position) pairs that skips digits on which all keys agree, then every run of equal keys is
/// ordered by the difference-cover comparator. The scratch lives in four arena blocks of K entries.
/// </summary>
internal sealed class LeafSorter
{
    private readonly Buf<ulong> _keys;
    private readonly Buf<ulong> _keys2;
    private readonly Buf<uint> _pos;
    private readonly Buf<uint> _pos2;
    private readonly int _capacity;

    internal LeafSorter(BudgetedArena arena, int capacity)
    {
        _capacity = capacity;
        _keys = arena.Allocate<ulong>(capacity, zero: false, "leaf.keys");
        _keys2 = arena.Allocate<ulong>(capacity, zero: false, "leaf.keys2");
        _pos = arena.Allocate<uint>(capacity, zero: false, "leaf.pos");
        _pos2 = arena.Allocate<uint>(capacity, zero: false, "leaf.pos2");
    }

    /// <summary>Sorts the bucket and writes its suffix array slice at <paramref name="saOffset"/>.</summary>
    internal void Sort(IPosSource bucket, long saOffset, SuffixOrder order, ISaSliceWriter writer)
    {
        int m = checked((int)bucket.Count);
        if (m > _capacity) throw new InvalidOperationException($"leaf of {m} positions exceeds the capacity {_capacity}");
        Span<uint> pos = _pos.AsSpan(m);
        bucket.Read(0, pos);
        Span<ulong> keys = _keys.AsSpan(m);
        MappedTextReader text = order.Text;
        for (int i = 0; i < m; i++) keys[i] = text.Key(pos[i]);

        if (m <= 24) KeySort.InsertionSort(keys, pos);
        else KeySort.RadixSort(_keys, _keys2, _pos, _pos2, m);

        keys = _keys.AsSpan(m);
        pos = _pos.AsSpan(m);
        SortWork c = order.Work;
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
                else
                {
                    MergeSort(pos.Slice(s, e - s), _pos2.AsSpan(e - s), order);
                }
            }
            s = e;
        }
        writer.Write(saOffset, pos);
    }

    /// <summary>
    /// Stable merge sort of the positions of a run of equal keys by the tied comparator. Chosen over the introsort of <see cref="MemoryExtensions"/> because the cost here is the number of comparator calls (each may compare
    /// hundreds of symbols or touch two sample ranks): a merge sort needs at most n - 1 per merge, only n / 2 when the run is in reverse order (the suffixes of a^n and of long tandem repeats arrive in text order, which is the
    /// reverse of their order), and a merge is skipped when the halves are already in order.
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
}
