using System.Diagnostics;
using System.Runtime.CompilerServices;

// The count and locate descent over a REAL container with the stored width swapped (ADR-005 section 5): a copy of MappedSearch.Run (src/SuffixTree/Algorithms/SuffixTree.Esa/MappedSearch.cs),
// Trusted access (no checks: they are the same for every width), generic over the reader of the suffix array and of the child cells. The suffix array and the child table are re-encoded from the
// container's 4-byte sections into anonymous memory in the layout of the variant; the LCP (U8 + overflow list) and the packed text are the container's.
//
//   variants (name = suffix array / child cells):
//     v0   w4 / w4 tagged, ints     the code of today (int arithmetic), the baseline of the cost of the long type
//     v1   w4 / w4 tagged, longs    (n <= 2^30 - 1 only: the tagged cell has 30 value bits)
//     v2   w4 / w4 packed, longs    (what the builder writes above 2^30: the kind of a cell is told by comparing LCP values)
//     v3   w5 / w5 tagged           40-bit entries, kind in bits 38-39
//     v4   w5 / w5 packed
//     v5   w8 / w8 tagged
//     v6   hyb / w5 tagged          hybrid SA (32-bit slots, 31 % of the entries escape: the threshold is scaled so that the share is that of n = 6.2e9), 40-bit child cells
//     v7   p34 / p34 packed         34 bits per entry, 15 per cache line
//   The packed child table is derived from the tagged one by the rule of the codec (a next link to k + 1 is not stored; the other fields keep their value), and checked equal by the results.

internal interface ILcpR
{
    long Get(long rank);
}

internal readonly unsafe struct Lcp8R(byte* values, uint* pairs, long pairCount) : ILcpR
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long Get(long rank)
    {
        int v = values[rank];
        return v != 255 ? v : Escaped(rank);
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private long Escaped(long rank)
    {
        long lo = 0, hi = pairCount - 1;
        while (lo <= hi)
        {
            long mid = (lo + hi) >> 1;
            uint r = pairs[2 * mid];
            if (r == rank) return pairs[(2 * mid) + 1];
            if (r < rank) lo = mid + 1; else hi = mid - 1;
        }
        throw new InvalidOperationException("no overflow entry for rank " + rank);
    }
}

internal interface ICellsR
{
    long First(long lb, long rb, bool lastChild);

    long Next(long k, long rb);
}

internal readonly struct TaggedCellsR<TR>(TR cells, int valueBits) : ICellsR where TR : struct, IReader
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long First(long lb, long rb, bool lastChild) => (long)(cells.Get(lastChild ? lb : rb) & ((1UL << valueBits) - 1));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long Next(long k, long rb)
    {
        ulong c = cells.Get(k);
        return (c >> valueBits) == 1 ? (long)(c & ((1UL << valueBits) - 1)) : -1;
    }
}

internal readonly struct PackedCellsR<TR, TL>(TR cells, TL lcp) : ICellsR where TR : struct, IReader where TL : struct, ILcpR
{
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long First(long lb, long rb, bool lastChild) => (long)cells.Get(lastChild ? lb : rb);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public long Next(long k, long rb)
    {
        if (k >= rb) return -1;
        long lk = lcp.Get(k);
        if (lcp.Get(k + 1) == lk) return k + 1;
        long v = (long)cells.Get(k);
        return lcp.Get(v) == lk ? v : -1;
    }
}

internal sealed unsafe class Patterns
{
    internal readonly int Length;
    internal readonly byte[] Symbols;   // count * Length dense symbols
    internal readonly int Count;

    internal Patterns(PackedText text, long n, int length, int count, ulong seed, bool present)
    {
        Length = length;
        Count = count;
        Symbols = new byte[(long)count * length];
        ulong x = seed;
        for (int q = 0; q < count; q++)
        {
            x = Rng.Mix(x + 1);
            if (present)
            {
                long start = (long)(((x & 0xFFFFFFFFUL) * (ulong)(n - length)) >> 32);
                for (int j = 0; j < length; j++) Symbols[((long)q * length) + j] = (byte)text.Symbol(start + j);
            }
            else
            {
                ulong r = Rng.Mix(x);
                for (int j = 0; j < length; j++) { Symbols[((long)q * length) + j] = (byte)(r & 3); r >>= 2; if (j % 30 == 29) r = Rng.Mix(r ^ x); }
            }
        }
    }
}

internal static unsafe class Descent
{
    private static (long lo, long hi) RunLong<TSa, TCells, TLcp>(TSa sa, TCells cells, TLcp lcp, PackedText text, long n, ReadOnlySpan<byte> pat, ulong[] words)
        where TSa : struct, IReader where TCells : struct, ICellsR where TLcp : struct, ILcpR
    {
        int bits = text.Bits;
        int m = pat.Length;
        long lb = 0, rb = n, ell = 0;
        int c = 0;
        bool isLast = true;
        bool packed = false;
        while (c < m)
        {
            int symbol = pat[c];
            long a = lb;
            long boundary = cells.First(lb, rb, isLast);
            long b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                long p = (long)sa.Get(a) + ell;
                int key = (ulong)p < (ulong)n ? text.Symbol(p) : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0) return (0, -1);
                a = boundary;
                boundary = cells.Next(a, rb);
            }

            bool childIsLast = boundary < 0;
            if (a == b)
            {
                long start = (long)sa.Get(a) + c;
                if ((ulong)start > (ulong)n || m - c > n - start) return (0, -1);
                if (!packed) { PackInto(pat, bits, words); packed = true; }
                if (!text.Equal(start, words, c, m - c)) return (0, -1);
                return (a, a);
            }

            long first = cells.First(a, b, childIsLast);
            long childEll = lcp.Get(first);
            long limit = Math.Min(childEll, m);
            long from = (long)sa.Get(a) + c;
            if (limit > c + 1)
            {
                if ((ulong)from >= (ulong)n || limit - c - 1 > n - (from + 1)) return (0, -1);
                if (!packed) { PackInto(pat, bits, words); packed = true; }
                if (!text.Equal(from + 1, words, c + 1, (int)(limit - c - 1))) return (0, -1);
            }
            c = (int)limit;
            lb = a;
            rb = b;
            ell = childEll;
            isLast = childIsLast;
        }
        return (lb, rb);
    }

    // the same loop with int arithmetic over the product's types (4-byte suffix array, tagged 4-byte cells, int ranks): the code of MappedSearch.Run
    private static (int lo, int hi) RunInt(int* sa, uint* cld, Lcp8R lcp8, PackedText text, int n, ReadOnlySpan<byte> pat, ulong[] words)
    {
        const uint ValueMask = (1u << 30) - 1;
        const uint KindMask = 3u << 30;
        const uint NextKind = 1u << 30;
        int bits = text.Bits;
        int m = pat.Length;
        int lb = 0, rb = n, ell = 0, c = 0;
        bool isLast = true;
        bool packed = false;
        while (c < m)
        {
            int symbol = pat[c];
            int a = lb;
            int boundary = (int)(cld[isLast ? lb : rb] & ValueMask);
            int b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                int p = sa[a] + ell;
                int key = (uint)p < (uint)n ? text.Symbol(p) : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0) return (0, -1);
                a = boundary;
                uint cell = cld[a];
                boundary = (cell & KindMask) == NextKind ? (int)(cell & ValueMask) : -1;
            }

            bool childIsLast = boundary < 0;
            if (a == b)
            {
                int start = sa[a] + c;
                if ((uint)start > (uint)n || m - c > n - start) return (0, -1);
                if (!packed) { PackInto(pat, bits, words); packed = true; }
                if (!text.Equal(start, words, c, m - c)) return (0, -1);
                return (a, a);
            }

            int first = (int)(cld[childIsLast ? a : b] & ValueMask);
            int childEll = (int)lcp8.Get(first);
            int limit = Math.Min(childEll, m);
            int from = sa[a] + c;
            if (limit > c + 1)
            {
                if ((uint)from >= (uint)n || limit - c - 1 > n - (from + 1)) return (0, -1);
                if (!packed) { PackInto(pat, bits, words); packed = true; }
                if (!text.Equal(from + 1, words, c + 1, limit - c - 1)) return (0, -1);
            }
            c = limit;
            lb = a;
            rb = b;
            ell = childEll;
            isLast = childIsLast;
        }
        return (lb, rb);
    }

    private static void PackInto(ReadOnlySpan<byte> pattern, int bits, ulong[] words)
    {
        Array.Clear(words);
        long bit = 0;
        for (int i = 0; i < pattern.Length; i++, bit += bits)
        {
            int w = (int)(bit >> 6);
            int room = 64 - (int)(bit & 63) - bits;
            ulong v = pattern[i];
            if (room >= 0) words[w] |= v << room;
            else
            {
                words[w] |= v >> -room;
                words[w + 1] |= v << (64 + room);
            }
        }
    }

    internal sealed class Ctx
    {
        internal Container C;
        internal PackedText Text;
        internal Lcp8R Lcp;
        internal long N;
        internal uint* SaSrc;      // the container's SA (4 bytes)
        internal uint* ChildSrc;   // the container's tagged child cells (4 bytes)
        internal bool Tagged;
        internal readonly Dictionary<string, Packed> Cache = new();
        internal readonly Dictionary<int, DeltaCells> DeltaCache = new();

        internal DeltaCells Delta(int width, Func<DeltaCells> make)
        {
            if (!DeltaCache.TryGetValue(width, out DeltaCells d)) DeltaCache[width] = d = make();
            return d;
        }

        /// <summary>An array of the variants, shared between the variants that use the same layout and freed with the context.</summary>
        internal Packed Arr(string key, Func<Packed> make)
        {
            if (!Cache.TryGetValue(key, out Packed p)) Cache[key] = p = make();
            return p;
        }
    }

    internal static int Run(Dictionary<string, string> o)
    {
        string path = Program.Get(o, "container");
        string[] variants = Program.Get(o, "variants", "v0,v1,v3,v5,v6").Split(',');
        int[] threads = Program.Get(o, "threads", "1,4").Split(',').Select(int.Parse).ToArray();
        int reps = int.Parse(Program.Get(o, "reps", "5"), Program.Inv);
        int q = int.Parse(Program.Get(o, "queries", "200000"), Program.Inv);
        int len = int.Parse(Program.Get(o, "length", "20"), Program.Inv);
        int locLen = int.Parse(Program.Get(o, "locate-length", "8"), Program.Inv);
        double virtualN = double.Parse(Program.Get(o, "virtual-n", "6.2e9"), Program.Inv);
        string outPath = Program.Get(o, "out");
        string tag = Program.Get(o, "tag", "");
        bool check = Program.Get(o, "check", "false") == "true";

        using var c = new Container(path);
        long n = c.N;
        Section child = c.Sections[6];
        var ctx = new Ctx
        {
            C = c,
            N = n,
            Text = new PackedText(c.At(1), (int)c.Sections[1].Param0),
            SaSrc = (uint*)c.At(3),
            ChildSrc = (uint*)c.At(6),
            Tagged = child.Codec == 2,
        };
        if (!ctx.Tagged) throw new InvalidOperationException("the container's child table is not tagged (codec 2); build a text of at most 2^30 - 1 symbols");
        ctx.Lcp = new Lcp8R(c.At(4), c.Has(5) ? (uint*)c.At(5) : null, c.Has(5) ? c.Sections[5].Count : 0);
        Console.WriteLine($"container {path}: n={n} bits={ctx.Text.Bits} sigma={c.Sigma} posWidth={c.PosWidth} childCodec={child.Codec} lcpKind={c.Sections[4].Codec} overflow={(c.Has(5) ? c.Sections[5].Count : 0)}");
        if (c.Sections[4].Codec != 1) throw new InvalidOperationException("only the U8 LCP is supported by this measurement");

        var present = new Patterns(ctx.Text, n, len, q, 11, true);
        var absent = new Patterns(ctx.Text, n, len, q, 12, false);
        var locate = new Patterns(ctx.Text, n, locLen, q / 20, 13, true);
        var rows = new List<Synthetic.Row>();
        long refSum = long.MinValue;

        var built = new List<Variant>();
        foreach (string v in variants)
        {
            var swBuild = Stopwatch.StartNew();
            Variant vr = Variant.Create(v, ctx, virtualN);
            Console.WriteLine($"{v}: built in {swBuild.Elapsed.TotalSeconds:F1} s, {vr.Bytes / 1e6:F0} MB for SA + child, {vr.Notes}");
            long sum = vr.Checksum(present, absent);
            Console.WriteLine($"{v}: checksum {sum}");
            if (refSum == long.MinValue) refSum = sum;
            else if (sum != refSum) { Console.WriteLine($"MISMATCH: checksum of {v} is {sum}, of the first variant {refSum}"); return 1; }
            if (check && built.Count == 0) CheckAgainstProduct(path, ctx, vr, present, absent);
            built.Add(vr);
        }
        string tagNotes = $"{tag} n={n} L={len} patterns={q}";
        foreach (int t in threads)
        {
            rows.AddRange(Interleave(built, "count20_present", t, reps, present, false, tagNotes, n));
            rows.AddRange(Interleave(built, "count20_absent", t, reps, absent, false, tagNotes, n));
            rows.AddRange(Interleave(built, "locate8", t, reps, locate, true, tagNotes, n));
        }
        foreach (Packed p in ctx.Cache.Values) p.Dispose();
        foreach (DeltaCells d in ctx.DeltaCache.Values) d.Dispose();
        Emit(outPath, rows);
        return 0;
    }

    private static void Emit(string path, List<Synthetic.Row> rows)
    {
        bool header = !File.Exists(path) || new FileInfo(path).Length == 0;
        using var w = new StreamWriter(path, true, System.Text.Encoding.ASCII);
        if (header) w.WriteLine("kind\tlayout\tn\tthreads\tmetric\tmedian\tmin\tmax\tunit\treps\tnotes");
        foreach (Synthetic.Row r in rows)
            w.WriteLine(string.Join('\t', r.Kind, r.Layout, r.N.ToString(Program.Inv), r.Threads.ToString(Program.Inv), r.Metric, r.Median.ToString("F1", Program.Inv), r.Min.ToString("F1", Program.Inv), r.Max.ToString("F1", Program.Inv), r.Unit, r.Reps.ToString(Program.Inv), r.Notes));
    }

    // One operation at one thread count for all variants, the repetitions interleaved (variant order rotates every repetition, so that a slow drift of the machine hits every variant alike);
    // two warm-up passes per variant first (tiered compilation with dynamic PGO, as the product runs). Result: ns per query per thread (wall time of a pass / (patterns x threads)).
    private static List<Synthetic.Row> Interleave(List<Variant> vs, string op, int threads, int reps, Patterns pats, bool locate, string notes, long n)
    {
        var samples = vs.Select(_ => new List<double>()).ToArray();
        for (int w = 0; w < 2; w++)
            foreach (Variant v in vs) Pass(v, threads, pats, locate);
        for (int r = 0; r < reps; r++)
            for (int k = 0; k < vs.Count; k++)
            {
                int vi = (k + r) % vs.Count;
                samples[vi].Add(Pass(vs[vi], threads, pats, locate));
            }
        var rows = new List<Synthetic.Row>();
        for (int vi = 0; vi < vs.Count; vi++)
        {
            samples[vi].Sort();
            List<double> sm = samples[vi];
            rows.Add(new Synthetic.Row("descent", vs[vi].Name, n, threads, op + "_ns_per_query_per_thread", sm[sm.Count / 2], sm[0], sm[^1], "ns/query", reps, $"{notes} bytes(SA+child)={vs[vi].Bytes} {vs[vi].Notes}"));
        }
        return rows;
    }

    private static double Pass(Variant vr, int threads, Patterns pats, bool locate)
    {
        var barrier = new Barrier(threads + 1);
        var ts = new Thread[threads];
        for (int t = 0; t < threads; t++)
        {
            int id = t;
            ts[t] = new Thread(() =>
            {
                ulong[] words = new ulong[16];
                barrier.SignalAndWait();
                int offset = (int)((long)pats.Count * id / threads);
                long acc = 0;
                for (int k = 0; k < pats.Count; k++)
                {
                    int qi = (offset + k) % pats.Count;
                    acc += vr.One(pats.Symbols.AsSpan(qi * pats.Length, pats.Length), words, locate);
                }
                GC.KeepAlive(acc);
            });
            ts[t].Start();
        }
        var sw = Stopwatch.StartNew();
        barrier.SignalAndWait();
        foreach (Thread t in ts) t.Join();
        sw.Stop();
        return sw.Elapsed.TotalMilliseconds * 1e6 / ((double)pats.Count * threads);
    }

    private static void CheckAgainstProduct(string path, Ctx ctx, Variant vr, Patterns present, Patterns absent)
    {
        global::SuffixTree.ISuffixTree tree = global::SuffixTree.Esa.SuffixIndexFactory.Open(path);
        char[] letters = ctx.C.Letters();
        int bad = 0;
        for (int i = 0; i < 3000; i++)
        {
            Patterns p = (i & 1) == 0 ? present : absent;
            int qi = (i * 7919) % p.Count;
            var span = p.Symbols.AsSpan(qi * p.Length, p.Length);
            var s = new string(span.ToArray().Select(x => letters[x]).ToArray());
            int want = tree.CountOccurrences(s);
            var words = new ulong[8];
            long got = vr.Count(span, words);
            if (want != got) { if (bad++ < 5) Console.WriteLine($"PRODUCT MISMATCH {s}: product {want} variant {got}"); }
        }
        Console.WriteLine(bad == 0 ? "check against the product: 3000 patterns equal" : $"check against the product: {bad} MISMATCHES");
        if (tree is IDisposable d) d.Dispose();
    }

    internal abstract class Variant : IDisposable
    {
        internal string Name;
        internal long Bytes;
        internal string Notes = "";
        internal abstract long One(ReadOnlySpan<byte> pat, ulong[] words, bool locate);
        internal abstract long Count(ReadOnlySpan<byte> pat, ulong[] words);
        public abstract void Dispose();

        internal long Checksum(Patterns present, Patterns absent)
        {
            long sum = 0;
            ulong[] words = new ulong[64];
            foreach (Patterns p in new[] { present, absent })
                for (int i = 0; i < Math.Min(p.Count, 50000); i++)
                    sum = (sum * 31) + Count(p.Symbols.AsSpan(i * p.Length, p.Length), words);
            return sum;
        }

        internal static Variant Create(string name, Ctx x, double virtualN)
        {
            long n = x.N;
            // the share of the suffix-array values that are "long" at the virtual text length: 1 - 2^32 / virtualN, here at a threshold in [0, n)
            ulong escapeAt = virtualN <= 4294967296.0 ? 0xFFFFFFFFUL : (ulong)Math.Min(4294967295.0, n * (4294967296.0 / virtualN));
            switch (name)
            {
                case "v0": return new V0(x);
                case "v1" or "v1b": return Make<R4, R4>(name, x, "w4", "w4", true, escapeAt, n, (p) => new R4(p.Main), (p) => new R4(p.Main), 30);
                case "v2": return Make<R4, R4>(name, x, "w4", "w4", false, escapeAt, n, (p) => new R4(p.Main), (p) => new R4(p.Main), 32);
                case "v3": return Make<R5, R5>(name, x, "w5", "w5", true, escapeAt, n, (p) => new R5(p.Main), (p) => new R5(p.Main), 38);
                case "v4": return Make<R5, R5>(name, x, "w5", "w5", false, escapeAt, n, (p) => new R5(p.Main), (p) => new R5(p.Main), 40);
                case "v5": return Make<R8, R8>(name, x, "w8", "w8", true, escapeAt, n, (p) => new R8(p.Main), (p) => new R8(p.Main), 62);
                case "v6": return MakeHyb(name, x, escapeAt, n);
                case "v8": return MakeDelta(name, x, 2, true, n, false);
                case "v9": return MakeDelta(name, x, 1, false, n, false);
                case "v8h": return MakeDelta(name, x, 2, true, n, true);
                case "v9h": return MakeDelta(name, x, 1, false, n, true);
                case "v7": return Make<R34, R34>(name, x, "p34", "p34", false, escapeAt, n, (p) => new R34(p.Main), (p) => new R34(p.Main), 34);
                default: throw new ArgumentException(name);
            }
        }

        private static Variant Make<TS, TC>(string name, Ctx x, string saLayout, string childLayout, bool tagged, ulong escapeAt, long n, Func<Packed, TS> mkSa, Func<Packed, TC> mkCell, int valueBits)
            where TS : struct, IReader where TC : struct, IReader
        {
            uint* sa4 = x.SaSrc;
            uint* ch4 = x.ChildSrc;
            Packed sa = x.Arr("sa:" + saLayout, () => Packed.CreateFrom(saLayout, n + 1, i => sa4[i], 0xFFFFFFFFUL));
            Func<long, ulong> cellValue = tagged ? i => Tag(ch4[i], valueBits, x.Tagged) : i => Untag(ch4, i);
            Packed ch = x.Arr($"cld:{childLayout}:{(tagged ? "t" : "p")}", () => Packed.CreateFrom(childLayout, n + 1, cellValue, 0xFFFFFFFFUL));
            return new VG<TS, TC>(name, x, sa, ch, mkSa(sa), mkCell(ch), tagged, valueBits, $"sa={saLayout} child={childLayout}{(tagged ? " tagged" : " packed")}");
        }

        private static Variant MakeDelta(string name, Ctx x, int width, bool tagged, long n, bool hashed)
        {
            Packed sa = x.Arr("sa:w4", () => Packed.CreateFrom("w4", n + 1, i => x.SaSrc[i], 0xFFFFFFFFUL));
            DeltaCells d = x.Delta(width + (hashed ? 10 : 0), () => new DeltaCells(width, x.ChildSrc, n, !tagged, hashed));
            string notes = $"sa=w4 child=d{8 * width}{(hashed ? "h" : "")} {(tagged ? "tagged" : "packed")} escapes={d.PairCount} ({100.0 * d.PairCount / (n + 1):F2} %) childBytes={d.Bytes}";
            return width == 2
                ? new VG<R4, RD16>(name, x, sa, null, new R4(sa.Main), new RD16(d), true, 40, notes) { Bytes = sa.Bytes + d.Bytes }
                : new VG<R4, RD8>(name, x, sa, null, new R4(sa.Main), new RD8(d), false, 0, notes) { Bytes = sa.Bytes + d.Bytes };
        }

        private static Variant MakeHyb(string name, Ctx x, ulong escapeAt, long n)
        {
            uint* sa4 = x.SaSrc;
            uint* ch4 = x.ChildSrc;
            Packed sa = x.Arr("sa:hyb", () => Packed.CreateFrom("hyb", n + 1, i => sa4[i], escapeAt));
            Packed ch = x.Arr("cld:w5:t", () => Packed.CreateFrom("w5", n + 1, i => Tag(ch4[i], 38, true), 0xFFFFFFFFUL));
            return new VG<RHyb, R5>(name, x, sa, ch, new RHyb(sa.Main, sa.Blocks, sa.Longs), new R5(ch.Main), true, 38, $"sa=hyb child=w5 tagged escapes={sa.Escapes} ({100.0 * sa.Escapes / (n + 1):F1} %)");
        }

        // a 4-byte tagged cell (kind in bits 30-31, value in bits 0-29) as a cell of valueBits value bits with the kind above them
        internal static ulong Tag(uint cell, int valueBits, bool _) => (ulong)(cell & ((1u << 30) - 1)) | ((ulong)(cell >> 30) << valueBits);

        // the packed form of cell i: the value of the cell; a next link to the following rank is not stored (the codec leaves it empty)
        private static ulong Untag(uint* ch4, long i)
        {
            uint cell = ch4[i];
            uint kind = cell >> 30;
            uint value = cell & ((1u << 30) - 1);
            if (kind == 1 && value == i + 1) return 0;
            return kind == 0 ? 0 : value;
        }
    }

    private sealed class V0 : Variant
    {
        private readonly Ctx _x;
        private readonly Packed _sa, _ch;
        private readonly int* _saP;
        private readonly uint* _cld;
        private readonly Lcp8R _lcp;
        private readonly int _n;

        internal V0(Ctx x)
        {
            _x = x;
            Name = "v0";
            _n = (int)x.N;
            uint* sa4 = x.SaSrc;
            uint* ch4 = x.ChildSrc;
            _sa = x.Arr("sa:w4", () => Packed.CreateFrom("w4", x.N + 1, i => sa4[i], 0xFFFFFFFFUL));
            _ch = x.Arr("cld:w4:t", () => Packed.CreateFrom("w4", x.N + 1, i => Variant.Tag(ch4[i], 30, true), 0xFFFFFFFFUL));
            _sa.Touch();
            _ch.Touch();
            _saP = (int*)_sa.Main;
            _cld = (uint*)_ch.Main;
            _lcp = x.Lcp;
            Bytes = _sa.Bytes + _ch.Bytes;
            Notes = "sa=w4 child=w4 tagged, int arithmetic (the code of MappedSearch.Run)";
        }

        internal override long One(ReadOnlySpan<byte> pat, ulong[] words, bool locate)
        {
            (int lo, int hi) = RunInt(_saP, _cld, _lcp, _x.Text, _n, pat, words);
            if (hi < lo) return 0;
            if (!locate) return hi - lo + 1L;
            long s = 0;
            for (int i = lo; i <= hi; i++) s += _saP[i];
            return s;
        }

        internal override long Count(ReadOnlySpan<byte> pat, ulong[] words)
        {
            (int lo, int hi) = RunInt(_saP, _cld, _lcp, _x.Text, _n, pat, words);
            return hi < lo ? 0 : hi - lo + 1L;
        }

        public override void Dispose() { }
    }

    private sealed class VG<TS, TC> : Variant where TS : struct, IReader where TC : struct, IReader
    {
        private readonly Ctx _x;
        private readonly Packed _sa, _ch;
        private readonly TS _saR;
        private readonly TC _chR;
        private readonly bool _tagged;
        private readonly int _vb;

        internal VG(string name, Ctx x, Packed sa, Packed ch, TS saR, TC chR, bool tagged, int vb, string notes)
        {
            Name = name;
            _x = x;
            _sa = sa;
            _ch = ch;
            _saR = saR;
            _chR = chR;
            _tagged = tagged;
            _vb = vb;
            Bytes = sa.Bytes + (ch?.Bytes ?? 0);
            Notes = notes;
            _sa.Touch();
            _ch?.Touch();
        }

        private (long lo, long hi) Find(ReadOnlySpan<byte> pat, ulong[] words)
            => _tagged
                ? RunLong(_saR, new TaggedCellsR<TC>(_chR, _vb), _x.Lcp, _x.Text, _x.N, pat, words)
                : RunLong(_saR, new PackedCellsR<TC, Lcp8R>(_chR, _x.Lcp), _x.Lcp, _x.Text, _x.N, pat, words);

        internal override long One(ReadOnlySpan<byte> pat, ulong[] words, bool locate)
        {
            (long lo, long hi) = Find(pat, words);
            if (hi < lo) return 0;
            if (!locate) return hi - lo + 1;
            return (long)_saR.SumRange(lo, hi - lo + 1);
        }

        internal override long Count(ReadOnlySpan<byte> pat, ulong[] words)
        {
            (long lo, long hi) = Find(pat, words);
            return hi < lo ? 0 : hi - lo + 1;
        }

        public override void Dispose() { }
    }
}
