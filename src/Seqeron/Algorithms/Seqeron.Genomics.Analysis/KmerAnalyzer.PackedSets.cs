namespace Seqeron.Genomics.Analysis;

public static partial class KmerAnalyzer
{
    // B09 WP78 (R14-18, F360): packed k-mer sets for k <= 32. The set operations of the Jaccard / containment / Mash estimators only need the sizes of the two distinct k-mer
    // sets and of their intersection; a window of A/C/G/T is a 2-bit-per-base ulong (A = 0, C = 1, G = 2, T = 3: the numeric order is the ordinal order of the strings, so the
    // canonical representative min(w, RC(w)) is the same k-mer as the one of FoldToCanonical), the windows are collected in an array, sorted and counted, and the intersection is a
    // merge. The result is the number the string sets give (DistinctKmers); only the cost changes (no string per k-mer, no hash table).
    // Not applicable (the string route is used): k > 32, or a sequence with a symbol other than A/C/G/T in literal counting (such a window is a k-mer of its own there).
    private static bool TryPackedSetSizes(string a, string b, int k, KmerCountingOptions options, out (int Shared, int CountA, int CountB) sizes)
    {
        sizes = default;
        if (k > 32 || !TryPackedDistinct(a, k, options, out var setA, out int countA) || !TryPackedDistinct(b, k, options, out var setB, out int countB))
            return false;

        int shared = 0, i = 0, j = 0;
        while (i < countA && j < countB)
        {
            if (setA[i] < setB[j])
                i++;
            else if (setA[i] > setB[j])
                j++;
            else
            {
                shared++;
                i++;
                j++;
            }
        }

        sizes = (shared, countA, countB);
        return true;
    }

    private static bool TryPackedDistinct(string sequence, int k, KmerCountingOptions options, out ulong[] kmers, out int count)
    {
        kmers = Array.Empty<ulong>();
        count = 0;
        if (string.IsNullOrEmpty(sequence) || k > sequence.Length)
            return true;

        var windows = new ulong[sequence.Length - k + 1];
        int n = 0, valid = 0;
        ulong mask = k == 32 ? ulong.MaxValue : (1UL << (2 * k)) - 1;
        int rcShift = 2 * (k - 1);
        ulong forward = 0, reverse = 0;
        foreach (char c in sequence)
        {
            int code = c switch { 'A' or 'a' => 0, 'C' or 'c' => 1, 'G' or 'g' => 2, 'T' or 't' => 3, _ => -1 };
            if (code < 0)
            {
                if (!options.SkipsNonAcgt)
                    return false;
                valid = 0;
                continue;
            }

            forward = ((forward << 2) | (uint)code) & mask;
            reverse = (reverse >> 2) | ((ulong)(3 - code) << rcShift);
            if (valid < k)
                valid++;
            if (valid == k)
                windows[n++] = options.Canonical && reverse < forward ? reverse : forward;
        }

        Array.Sort(windows, 0, n);
        int distinct = 0;
        for (int i = 0; i < n; i++)
        {
            if (i == 0 || windows[i] != windows[i - 1])
                windows[distinct++] = windows[i];
        }

        kmers = windows;
        count = distinct;
        return true;
    }
}
