namespace SuffixTree.Esa.Bench.Links;

// Work package 1b.9f (the hand-off package, analysis/s1b-handoff). The three pieces of the query check of `sesa handoff queries`. They live in Links/ because that is the one folder of the Bench application that
// the test project compiles by a glob (see SuffixTree.Esa.Tests.csproj); they are not suffix-link prototypes. None of them uses the index under test: the oracle reads only the input text, the search reads
// the suffix array of the container and the input text.

/// <summary>One query of the set: its number, whether it was cut from the text, where (-1 for an absent one) and its bytes.</summary>
internal sealed record HandoffQuery(int Index, bool Present, long Origin, byte[] Bytes);

/// <summary>What the scan found: for every query the exact number of occurrences and the first <c>cap</c> start positions in text order.</summary>
internal sealed record ScanResult(long[] Counts, List<long>[] Positions);

/// <summary>The fixed deterministic query set: windows cut from the text itself ("present", at least one occurrence by construction) and random strings over the alphabet of the text ("absent", usually none).</summary>
internal static class HandoffQuerySet
{
    /// <summary>The distinct bytes of <paramref name="sample"/>, in ascending order.</summary>
    internal static byte[] AlphabetOf(ReadOnlySpan<byte> sample)
    {
        Span<bool> seen = stackalloc bool[256];
        foreach (byte b in sample) seen[b] = true;
        var letters = new List<byte>();
        for (int b = 0; b < 256; b++)
        {
            if (seen[b]) letters.Add((byte)b);
        }
        return [.. letters];
    }

    /// <summary>
    /// <paramref name="present"/> windows of every length cut at uniformly random positions (<paramref name="read"/> gives the text bytes at a position) and <paramref name="absent"/> random strings of every length over
    /// <paramref name="alphabet"/>, all from one SplitMix64 stream of <paramref name="seed"/>: the same text, parameters and seed give the same set on every machine. A length above the text gives no present query.
    /// </summary>
    internal static IReadOnlyList<HandoffQuery> Build(long n, Func<long, int, byte[]> read, ReadOnlySpan<byte> alphabet, IReadOnlyList<int> lengths, int present, int absent, ulong seed)
    {
        ArgumentNullException.ThrowIfNull(read);
        ArgumentNullException.ThrowIfNull(lengths);
        ulong state = seed;
        ulong Next()
        {
            unchecked
            {
                state += 0x9E3779B97F4A7C15UL;
                ulong z = state;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
                return z ^ (z >> 31);
            }
        }

        var list = new List<HandoffQuery>();
        foreach (int length in lengths)
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(length, 1);
            if (length <= n)
            {
                ulong windows = (ulong)(n - length + 1);
                for (int i = 0; i < present; i++)
                {
                    long origin = (long)(Next() % windows);
                    list.Add(new HandoffQuery(list.Count, true, origin, read(origin, length)));
                }
            }
            if (alphabet.IsEmpty) continue;
            for (int i = 0; i < absent; i++)
            {
                var bytes = new byte[length];
                for (int k = 0; k < length; k++) bytes[k] = alphabet[(int)(Next() % (ulong)alphabet.Length)];
                list.Add(new HandoffQuery(list.Count, false, -1, bytes));
            }
        }
        return list;
    }
}

/// <summary>
/// The oracle: one pass over the text, a rolling hash (a polynomial of the window modulo 2^64) for every distinct query length at once, and every window whose hash is in the table is compared with the
/// queries of that hash byte by byte, so that a hash collision can add nothing and the same window always has the same hash, so that nothing is missed. The result is the exact count of every query and the
/// first <c>cap</c> positions in text order. Cost: one multiplication and one table probe per symbol and length; memory: one block of the text.
/// </summary>
internal static class RollingScan
{
    private const ulong Base = 0x9E3779B97F4A7C15UL | 1;
    private const int FilterBits = 20;

    private sealed class Group
    {
        internal required int Length { get; init; }

        internal required ulong Power { get; init; }

        internal required Dictionary<ulong, List<int>> ByHash { get; init; }

        internal required ulong[] Filter { get; init; }

        internal ulong Hash;
    }

    /// <summary>Scans <paramref name="text"/> (read once, in blocks of at most <paramref name="blockBytes"/>) for the <paramref name="queries"/>.</summary>
    internal static ScanResult Run(Stream text, IReadOnlyList<HandoffQuery> queries, int cap, int blockBytes = 4 << 20)
    {
        ArgumentNullException.ThrowIfNull(text);
        ArgumentNullException.ThrowIfNull(queries);
        ArgumentOutOfRangeException.ThrowIfLessThan(blockBytes, 1);
        var counts = new long[queries.Count];
        var positions = new List<long>[queries.Count];
        for (int i = 0; i < positions.Length; i++) positions[i] = [];
        Group[] groups = [.. queries.Where(q => q.Bytes.Length > 0).GroupBy(q => q.Bytes.Length).Select(g => MakeGroup(g.Key, g))];
        if (groups.Length == 0) return new ScanResult(counts, positions);

        int carryMax = groups.Max(g => g.Length);
        var buffer = new byte[carryMax + blockBytes];
        long consumed = 0;
        int carry = 0;
        while (true)
        {
            int got = text.Read(buffer, carry, blockBytes);
            if (got <= 0) break;
            foreach (Group g in groups) ScanBlock(g, buffer, carry, got, consumed, queries, counts, positions, cap);
            consumed += got;
            int keep = (int)Math.Min(carryMax, consumed);
            Array.Copy(buffer, carry + got - keep, buffer, 0, keep);
            carry = keep;
        }
        return new ScanResult(counts, positions);
    }

    private static Group MakeGroup(int length, IEnumerable<HandoffQuery> members)
    {
        var byHash = new Dictionary<ulong, List<int>>();
        var filter = new ulong[1 << (FilterBits - 6)];
        ulong power = 1;
        unchecked
        {
            for (int i = 0; i < length; i++) power *= Base;
        }
        foreach (HandoffQuery q in members)
        {
            ulong h = HashOf(q.Bytes);
            if (!byHash.TryGetValue(h, out List<int>? list)) byHash[h] = list = [];
            list.Add(q.Index);
            int bit = Slot(h);
            filter[bit >> 6] |= 1UL << (bit & 63);
        }
        return new Group { Length = length, Power = power, ByHash = byHash, Filter = filter };
    }

    private static int Slot(ulong h) => (int)((h * 0xD6E8FEB86659FD93UL) >> (64 - FilterBits));

    private static ulong HashOf(ReadOnlySpan<byte> bytes)
    {
        ulong h = 0;
        unchecked
        {
            foreach (byte b in bytes) h = (h * Base) + b;
        }
        return h;
    }

    private static void ScanBlock(Group g, byte[] buffer, int carry, int got, long consumed, IReadOnlyList<HandoffQuery> queries, long[] counts, List<long>[] positions, int cap)
    {
        int length = g.Length;
        ulong h = g.Hash;
        ulong power = g.Power;
        ulong[] filter = g.Filter;
        for (int k = 0; k < got; k++)
        {
            int at = carry + k;
            long global = consumed + k;
            unchecked
            {
                h = (h * Base) + buffer[at];
                if (global >= length) h -= buffer[at - length] * power;
            }
            if (global < length - 1) continue;
            int bit = Slot(h);
            if ((filter[bit >> 6] & (1UL << (bit & 63))) == 0) continue;
            if (!g.ByHash.TryGetValue(h, out List<int>? candidates)) continue;
            ReadOnlySpan<byte> window = buffer.AsSpan(at - length + 1, length);
            foreach (int index in candidates)
            {
                if (!window.SequenceEqual(queries[index].Bytes)) continue;
                counts[index]++;
                if (positions[index].Count < cap) positions[index].Add(global - length + 1);
            }
        }
        g.Hash = h;
    }
}

/// <summary>The binary search of a pattern over a suffix array and the text it indexes, without any index structure: two reads per step.</summary>
internal static class SaSearch
{
    /// <summary>
    /// The ranks <c>[Lo, Hi)</c> of the suffixes that start with <paramref name="query"/>. <paramref name="sa"/> gives the position at a rank (the array has <c>n + 1</c> entries, the empty suffix first), <paramref name="readText"/> fills
    /// the buffer with the text from a position (at most to the end of the text) and returns the number of bytes. A suffix that is a proper prefix of the query sorts before it.
    /// </summary>
    internal static (long Lo, long Hi) Range(ReadOnlySpan<byte> query, long n, Func<long, long> sa, Func<long, byte[], int> readText)
    {
        ArgumentNullException.ThrowIfNull(sa);
        ArgumentNullException.ThrowIfNull(readText);
        var buffer = new byte[query.Length];

        long rows = n + 1;
        long lo = 0, hi = rows;
        while (lo < hi)
        {
            long mid = lo + ((hi - lo) / 2);
            if (Compare(query, buffer, sa(mid), readText) < 0) lo = mid + 1;
            else hi = mid;
        }
        long first = lo;
        hi = rows;
        while (lo < hi)
        {
            long mid = lo + ((hi - lo) / 2);
            if (Compare(query, buffer, sa(mid), readText) <= 0) lo = mid + 1;
            else hi = mid;
        }
        return (first, lo);
    }

    private static int Compare(ReadOnlySpan<byte> query, byte[] buffer, long position, Func<long, byte[], int> readText)
    {
        int got = readText(position, buffer);
        int c = buffer.AsSpan(0, got).SequenceCompareTo(query[..got]);
        if (c != 0) return c < 0 ? -1 : 1;
        return got < query.Length ? -1 : 0;
    }
}
