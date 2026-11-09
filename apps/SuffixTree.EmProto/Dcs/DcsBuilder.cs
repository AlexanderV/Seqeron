using System.Diagnostics;
using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Dcs;

/// <summary>What the construction of the sample ranks did.</summary>
internal sealed class DcsStats
{
    internal long SampleSize { get; set; }

    internal bool NamesUnique { get; set; }

    /// <summary>Recursion levels below this one (0 = the names were unique or the reduced string went to SA-IS directly).</summary>
    internal int Levels { get; set; }

    /// <summary>Runs of sample positions with equal 64-bit keys that needed the X-prefix comparator.</summary>
    internal long TieRuns { get; set; }
}

/// <summary>
/// PH1 of ADR-004, in RAM (E2): ranks of the difference-cover sample suffixes. The sample positions are sorted by their X-prefix and named; when all names are distinct the names are the ranks, otherwise the
/// names of every residue class (positions <c>d, d+X, d+2X, ...</c>) form a reduced string (classes concatenated, separated by 0) whose suffix array, by SA-IS up to <c>baseCase</c> symbols and by this builder
/// (recursively, with the same cover) above, orders the sample suffixes. The comparison of two reduced suffixes is decided before a separator is reached because the last name of a class is unique.
/// </summary>
internal static unsafe class DcsBuilder
{
    internal static DcsRankTable Build(PackedText text, DifferenceCover cover, WorkCounters counters, long baseCase, out DcsStats stats, int level = 0)
    {
        stats = new DcsStats();
        long n = text.Length;
        int x = cover.Period;
        int m = cover.Residues.Length;
        long slots = ((n / x) + 1) * m;
        long s = cover.SampleCount(n);
        stats.SampleSize = s;

        // sample positions in class-major order (the order of the reduced string), class index of each class start
        var pos = new uint[s];
        var classStart = new long[m + 1];
        long at = 0;
        for (int c = 0; c < m; c++)
        {
            classStart[c] = at;
            for (long p = cover.Residues[c]; p <= n; p += x) pos[at++] = (uint)p;
        }
        classStart[m] = at;

        int[] name = NameSample(text, x, pos, counters, stats, out int names);
        stats.NamesUnique = names == s;

        var ranks = new NativeBlock(Math.Max(slots, 1) * 4);
        uint* r = (uint*)ranks.Pointer;
        if (stats.NamesUnique)
        {
            for (int c = 0; c < m; c++)
            {
                for (long i = classStart[c]; i < classStart[c + 1]; i++) r[SlotOf(pos[i], x, m, c)] = (uint)name[i];
            }
        }
        else
        {
            // reduced string: names (1..), a 0 after every class
            var red = new int[s + m];
            for (int c = 0; c < m; c++)
            {
                for (long i = classStart[c]; i < classStart[c + 1]; i++) red[i + c] = name[i];
                red[classStart[c + 1] + c] = 0;
            }
            int[] saOfReduced = SuffixArrayOfNames(red, names + 1, cover, counters, baseCase, level, out int below);
            stats.Levels = below + 1;
            var rankOfReduced = new int[red.Length];
            for (int k = 1; k < saOfReduced.Length; k++) rankOfReduced[saOfReduced[k]] = k;
            for (int c = 0; c < m; c++)
            {
                for (long i = classStart[c]; i < classStart[c + 1]; i++) r[SlotOf(pos[i], x, m, c)] = (uint)rankOfReduced[i + c];
            }
        }
        return new DcsRankTable(cover, ranks, n);
    }

    /// <summary>Sorts the sample by (64-bit key, X-prefix) and gives equal X-prefixes the same name (names start at 1 and increase with the order).</summary>
    private static int[] NameSample(PackedText text, int x, uint[] pos, WorkCounters counters, DcsStats stats, out int names)
    {
        long s = pos.Length;
        var order = new uint[s];
        var keys = new ulong[s];
        for (long i = 0; i < s; i++)
        {
            order[i] = (uint)i;
            keys[i] = text.Key(pos[i]);
        }
        Array.Sort(keys, order);
        var prefix = new PrefixComparer(text, x, pos);
        long ties = 0;
        for (long i = 0; i < s;)
        {
            long j = i + 1;
            while (j < s && keys[j] == keys[i]) j++;
            if (j - i >= 2)
            {
                ties++;
                order.AsSpan((int)i, (int)(j - i)).Sort(prefix);
            }
            i = j;
        }
        stats.TieRuns = ties;
        var name = new int[s];
        names = 0;
        for (long i = 0; i < s; i++)
        {
            if (i == 0 || keys[i] != keys[i - 1] || prefix.Compare(order[i - 1], order[i]) != 0) names++;
            name[order[i]] = names;
        }
        counters.LceSteps += prefix.Steps;
        return name;
    }

    private static long SlotOf(uint p, int x, int m, int c) => ((p / x) * (long)m) + c;

    private static int[] SuffixArrayOfNames(int[] red, int sigma, DifferenceCover cover, WorkCounters counters, long baseCase, int level, out int below)
    {
        below = 0;
        if (red.Length <= baseCase) return SaisReference.Build(new IntSymbols(red), sigma);

        // recursion: the reduced string is a text over a 32-bit alphabet; sort it with this builder and the same cover
        int bits = PackedText.BitsFor(sigma);
        using PackedText t = PackedText.Pack(red, bits);
        var options = new BuildOptions { LeafCapacity = Math.Min(1 << 20, Math.Max(2, red.Length / 64)), FanOut = 64, Threads = 1, Seed = 7 + (ulong)level };
        using DcsRankTable sub = Build(t, cover, counters, baseCase, out DcsStats inner, level + 1);
        below = inner.Levels;
        using var sink = new MemorySaSink(red.Length + 1L);
        EmBuilder.Run(t, options, sink, sub);
        return [.. sink.ToArray().Select(v => (int)v)];
    }

    /// <summary>Order of the X-prefixes of sample positions (a position index into the sample array), with the end of the text smaller than every symbol; 0 means equal X-prefixes.</summary>
    private sealed class PrefixComparer(PackedText text, int x, uint[] pos) : IComparer<uint>
    {
        private long _steps;

        internal long Steps => _steps;

        public int Compare(uint a, uint b)
        {
            if (a == b) return 0;
            long pa = pos[a], pb = pos[b];
            long maxLen = text.Length - Math.Max(pa, pb);
            long lim = Math.Min(x, maxLen);
            long l = text.Lce(pa, pb, lim, ref _steps);
            if (l >= x) return 0;
            if (l == maxLen) return pa > pb ? -1 : 1;
            return text.Symbol(pa + l) < text.Symbol(pb + l) ? -1 : 1;
        }
    }
}
