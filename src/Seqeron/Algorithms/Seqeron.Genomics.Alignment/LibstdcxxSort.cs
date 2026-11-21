namespace Seqeron.Genomics.Alignment;

/// <summary>
/// <c>std::sort</c> of libstdc++ (GCC 4.x - 14, <c>bits/stl_algo.h</c>): introsort (median of three pivot, recursion depth limit 2 floor(log2 n), heap sort fallback) down to partitions of 16 elements and a final
/// insertion sort. The sort is not stable, so the order of elements that compare equal is whatever this exact algorithm produces; programs that print the result (MUMmer <c>postnuc</c> sorts clusters by their
/// first start) are reproduced only by running the same permutation.
/// </summary>
internal static class LibstdcxxSort
{
    private const int Threshold = 16;

    /// <summary>Sorts <paramref name="a"/> in place exactly as <c>std::sort (a.begin (), a.end (), comp)</c> does; <paramref name="less"/> is the strict-weak-order comparator.</summary>
    public static void Sort<T>(IList<T> a, Func<T, T, bool> less)
    {
        int n = a.Count;
        if (n == 0)
            return;
        IntroSortLoop(a, 0, n, 2 * Lg(n), less);
        FinalInsertionSort(a, 0, n, less);
    }

    private static int Lg(int n)
    {
        int k = 0;
        for (; n > 1; n >>= 1)
            k++;
        return k;
    }

    private static void Swap<T>(IList<T> a, int i, int j) => (a[i], a[j]) = (a[j], a[i]);

    private static void IntroSortLoop<T>(IList<T> a, int first, int last, int depthLimit, Func<T, T, bool> less)
    {
        while (last - first > Threshold)
        {
            if (depthLimit == 0)
            {
                PartialSort(a, first, last, last, less);
                return;
            }

            depthLimit--;
            int cut = UnguardedPartitionPivot(a, first, last, less);
            IntroSortLoop(a, cut, last, depthLimit, less);
            last = cut;
        }
    }

    private static int UnguardedPartitionPivot<T>(IList<T> a, int first, int last, Func<T, T, bool> less)
    {
        int mid = first + (last - first) / 2;
        MoveMedianToFirst(a, first, first + 1, mid, last - 1, less);
        return UnguardedPartition(a, first + 1, last, first, less);
    }

    private static void MoveMedianToFirst<T>(IList<T> a, int result, int x, int y, int z, Func<T, T, bool> less)
    {
        if (less(a[x], a[y]))
        {
            if (less(a[y], a[z]))
                Swap(a, result, y);
            else if (less(a[x], a[z]))
                Swap(a, result, z);
            else
                Swap(a, result, x);
        }
        else if (less(a[x], a[z]))
            Swap(a, result, x);
        else if (less(a[y], a[z]))
            Swap(a, result, z);
        else
            Swap(a, result, y);
    }

    private static int UnguardedPartition<T>(IList<T> a, int first, int last, int pivot, Func<T, T, bool> less)
    {
        while (true)
        {
            while (less(a[first], a[pivot]))
                first++;
            last--;
            while (less(a[pivot], a[last]))
                last--;
            if (first >= last)
                return first;
            Swap(a, first, last);
            first++;
        }
    }

    private static void FinalInsertionSort<T>(IList<T> a, int first, int last, Func<T, T, bool> less)
    {
        if (last - first > Threshold)
        {
            InsertionSort(a, first, first + Threshold, less);
            for (int i = first + Threshold; i != last; i++)
                UnguardedLinearInsert(a, i, less);
        }
        else
        {
            InsertionSort(a, first, last, less);
        }
    }

    private static void InsertionSort<T>(IList<T> a, int first, int last, Func<T, T, bool> less)
    {
        if (first == last)
            return;
        for (int i = first + 1; i != last; i++)
        {
            if (less(a[i], a[first]))
            {
                T val = a[i];
                for (int k = i; k > first; k--)
                    a[k] = a[k - 1];
                a[first] = val;
            }
            else
            {
                UnguardedLinearInsert(a, i, less);
            }
        }
    }

    private static void UnguardedLinearInsert<T>(IList<T> a, int last, Func<T, T, bool> less)
    {
        T val = a[last];
        int next = last - 1;
        while (less(val, a[next]))
        {
            a[last] = a[next];
            last = next;
            next--;
        }

        a[last] = val;
    }

    // std::partial_sort (first, middle, last) with middle == last: heap_select then sort_heap
    private static void PartialSort<T>(IList<T> a, int first, int middle, int last, Func<T, T, bool> less)
    {
        // __heap_select
        MakeHeap(a, first, middle, less);
        for (int i = middle; i < last; i++)
            if (less(a[i], a[first]))
                PopHeap(a, first, middle, i, less);
        // __sort_heap
        while (middle - first > 1)
        {
            middle--;
            PopHeap(a, first, middle, middle, less);
        }
    }

    private static void MakeHeap<T>(IList<T> a, int first, int last, Func<T, T, bool> less)
    {
        int len = last - first;
        if (len < 2)
            return;
        int parent = (len - 2) / 2;
        while (true)
        {
            T value = a[first + parent];
            AdjustHeap(a, first, parent, len, value, less);
            if (parent == 0)
                return;
            parent--;
        }
    }

    private static void PopHeap<T>(IList<T> a, int first, int last, int result, Func<T, T, bool> less)
    {
        T value = a[result];
        a[result] = a[first];
        AdjustHeap(a, first, 0, last - first, value, less);
    }

    private static void AdjustHeap<T>(IList<T> a, int first, int holeIndex, int len, T value, Func<T, T, bool> less)
    {
        int topIndex = holeIndex;
        int secondChild = holeIndex;
        while (secondChild < (len - 1) / 2)
        {
            secondChild = 2 * (secondChild + 1);
            if (less(a[first + secondChild], a[first + secondChild - 1]))
                secondChild--;
            a[first + holeIndex] = a[first + secondChild];
            holeIndex = secondChild;
        }

        if ((len & 1) == 0 && secondChild == (len - 2) / 2)
        {
            secondChild = 2 * (secondChild + 1);
            a[first + holeIndex] = a[first + secondChild - 1];
            holeIndex = secondChild - 1;
        }

        // __push_heap
        int parent = (holeIndex - 1) / 2;
        while (holeIndex > topIndex && less(a[first + parent], value))
        {
            a[first + holeIndex] = a[first + parent];
            holeIndex = parent;
            parent = (holeIndex - 1) / 2;
        }

        a[first + holeIndex] = value;
    }
}
