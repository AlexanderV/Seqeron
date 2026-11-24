using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics;

// The read cost of a packed DNA text against the byte text (ADR-005 section 6): the same symbols (the dense codes 0..3 of an ACGT text) stored as one byte per symbol, as 2 bits and as 3 bits (what a text with N needs
// today: b = ceil(log2 sigma)); micro-benchmarks of
//   rand   one random symbol read, independent reads (the pattern of the induction and naming passes of SA-IS and of the leaf compare of a descent)
//   chase  one random read whose result decides the next index (a dependent chain: the latency)
//   scan   a sequential pass (symbols per second)
//   lce    the longest common extension of two positions (the compare of the PLCP loop of the LCP kernel and of the leaf compare), on the text itself (random: about 1 symbol) and on a text whose positions
//          are 100,000 apart in a repeated block (rep_10m: about 50 symbols)
// usage: text --text <ACGT file> --rep <ACGT file> --reps 7 --out <tsv> --tag <text> --ops 20000000

internal static unsafe class TextBench
{
    private interface ISym
    {
        int Get(long i);
    }

    private readonly struct SByte(byte* p) : ISym
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Get(long i) => p[i];
    }

    private readonly struct SBits(byte* p, int bits) : ISym
    {
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public int Get(long i) => new PackedText(p, bits).Symbol(i);
    }

    private sealed class Store : IDisposable
    {
        internal byte* Bytes;
        internal byte* B2;
        internal byte* B3;
        internal long N;

        internal static Store Load(string path)
        {
            byte[] raw = File.ReadAllBytes(path);
            var s = new Store { N = raw.Length };
            s.Bytes = Mem.Alloc(s.N);
            for (long i = 0; i < s.N; i++)
            {
                byte c = raw[i];
                s.Bytes[i] = c == (byte)'A' ? (byte)0 : c == (byte)'C' ? (byte)1 : c == (byte)'G' ? (byte)2 : (byte)3;
            }
            s.B2 = Pack(s.Bytes, s.N, 2);
            s.B3 = Pack(s.Bytes, s.N, 3);
            return s;
        }

        private static byte* Pack(byte* sym, long n, int bits)
        {
            long total = (((n * bits) + 7) / 8) + 64;
            byte* p = Mem.Alloc(total);
            Unsafe.InitBlock(p, 0, (uint)total);
            for (long i = 0; i < n; i++)
            {
                long bit = i * bits;
                for (int k = 0; k < bits; k++)
                {
                    if (((sym[i] >> (bits - 1 - k)) & 1) != 0) p[(bit + k) >> 3] |= (byte)(0x80 >> (int)((bit + k) & 7));
                }
            }
            return p;
        }

        public void Dispose()
        {
            Mem.Free(Bytes);
            Mem.Free(B2);
            Mem.Free(B3);
        }
    }

    internal static int Run(Dictionary<string, string> o)
    {
        string textPath = Program.Get(o, "text");
        string repPath = Program.Get(o, "rep");
        int reps = int.Parse(Program.Get(o, "reps", "7"), Program.Inv);
        long ops = long.Parse(Program.Get(o, "ops", "20000000"), Program.Inv);
        string outPath = Program.Get(o, "out");
        string tag = Program.Get(o, "tag", "");
        var rows = new List<Synthetic.Row>();
        using Store dna = Store.Load(textPath);
        using Store rep = Store.Load(repPath);
        // self test: the three stores agree
        for (long i = 0; i < dna.N; i += 9973)
            if (dna.Bytes[i] != new PackedText(dna.B2, 2).Symbol(i) || dna.Bytes[i] != new PackedText(dna.B3, 3).Symbol(i)) throw new InvalidOperationException("pack mismatch at " + i);
        foreach ((string kind, Store s) in new[] { ("dna", dna), ("rep", rep) })
        {
            string notes = $"{tag} {kind} n={s.N}";
            foreach (string layout in new[] { "byte", "b2", "b3" })
            {
                rows.Add(Time(layout, s, "rand", reps, ops, 0, notes));
                rows.Add(Time(layout, s, "chase", reps, ops / 4, 0, notes));
                rows.Add(Time(layout, s, "scan", reps, ops * 2, 0, notes));
            }
            foreach (string layout in new[] { "byte8", "byte256", "b2", "b3" })
            {
                rows.Add(Time(layout, s, "lce_self", reps, ops / 4, 0, notes));
                rows.Add(Time(layout, s, "lce_far", reps, ops / 4, 100_000, notes));
            }
        }
        bool header = !File.Exists(outPath) || new FileInfo(outPath).Length == 0;
        using var w = new StreamWriter(outPath, true, System.Text.Encoding.ASCII);
        if (header) w.WriteLine("kind\tlayout\tn\tthreads\tmetric\tmedian\tmin\tmax\tunit\treps\tnotes");
        foreach (Synthetic.Row r in rows)
            w.WriteLine(string.Join('\t', "text", r.Layout, r.N.ToString(Program.Inv), "1", r.Metric, r.Median.ToString("F3", Program.Inv), r.Min.ToString("F3", Program.Inv), r.Max.ToString("F3", Program.Inv), r.Unit, r.Reps.ToString(Program.Inv), r.Notes));
        return 0;
    }

    private static Synthetic.Row Time(string layout, Store s, string op, int reps, long ops, long distance, string notes)
    {
        var samples = new List<double>();
        for (int r = -1; r < reps; r++)
        {
            var sw = Stopwatch.StartNew();
            long sink = Dispatch(layout, s, op, ops, distance, (ulong)(r + 2));
            sw.Stop();
            GC.KeepAlive(sink);
            if (r >= 0) samples.Add(sw.Elapsed.TotalMilliseconds * 1e6 / ops);
        }
        samples.Sort();
        return new Synthetic.Row("text", layout, s.N, 1, op + "_ns_per_op", samples[samples.Count / 2], samples[0], samples[^1], "ns/op", reps, notes);
    }

    private static long Dispatch(string layout, Store s, string op, long ops, long distance, ulong seed)
    {
        switch (layout)
        {
            case "byte": return Loop(new SByte(s.Bytes), s.N, op, ops, seed);
            case "b2": return op.StartsWith("lce") ? Lce(new SBits(s.B2, 2), 2, s, s.B2, op, ops, distance, seed) : Loop(new SBits(s.B2, 2), s.N, op, ops, seed);
            case "b3": return op.StartsWith("lce") ? Lce(new SBits(s.B3, 3), 3, s, s.B3, op, ops, distance, seed) : Loop(new SBits(s.B3, 3), s.N, op, ops, seed);
            case "byte8": return LceByte(s, ops, distance, seed, false);
            case "byte256": return LceByte(s, ops, distance, seed, true);
            default: throw new ArgumentException(layout);
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static long Loop<TS>(TS t, long n, string op, long ops, ulong seed) where TS : struct, ISym
    {
        long s = 0;
        switch (op)
        {
            case "rand":
            {
                ulong x = seed | 1;
                for (long k = 0; k < ops; k++)
                {
                    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
                    s += t.Get((long)(((x & 0xFFFFFFFFUL) * (ulong)n) >> 32));
                }
                return s;
            }
            case "chase":
            {
                long idx = (long)(Rng.Mix(seed) % (ulong)n);
                for (long k = 0; k < ops; k++)
                {
                    int v = t.Get(idx);
                    s += v;
                    idx = (long)(((Rng.Mix((ulong)idx + (ulong)v) & 0xFFFFFFFFUL) * (ulong)n) >> 32);
                }
                return s;
            }
            case "scan":
            {
                long idx = 0;
                for (long k = 0; k < ops; k++)
                {
                    s += t.Get(idx);
                    if (++idx == n) idx = 0;
                }
                return s;
            }
            default: throw new ArgumentException(op);
        }
    }

    // the common extension of the suffixes at i and j = i + distance (distance 0: j is another random position), capped at 4,096 symbols
    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static long Lce<TS>(TS t, int bits, Store s, byte* packed, string op, long ops, long distance, ulong seed) where TS : struct, ISym
    {
        long total = 0;
        ulong x = seed | 1;
        long n = s.N;
        var pt = new PackedText(packed, bits);
        for (long k = 0; k < ops; k++)
        {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            long i = (long)(((x & 0xFFFFFFFFUL) * (ulong)(n - distance - 5000)) >> 32);
            long j = distance > 0 ? i + distance : (long)((((x >> 32) | 1) * (ulong)(n - 5000)) >> 32);
            int per = 64 / bits;
            long l = 0;
            while (l < 4096)
            {
                ulong d = pt.Window(i + l) ^ pt.Window(j + l);
                if (d != 0)
                {
                    l += BitOperations.LeadingZeroCount(d) / bits;
                    break;
                }
                l += per;
            }
            total += l;
        }
        return total;
    }

    [MethodImpl(MethodImplOptions.AggressiveOptimization)]
    private static long LceByte(Store s, long ops, long distance, ulong seed, bool vector)
    {
        long total = 0;
        ulong x = seed | 1;
        long n = s.N;
        byte* b = s.Bytes;
        for (long k = 0; k < ops; k++)
        {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            long i = (long)(((x & 0xFFFFFFFFUL) * (ulong)(n - distance - 5000)) >> 32);
            long j = distance > 0 ? i + distance : (long)((((x >> 32) | 1) * (ulong)(n - 5000)) >> 32);
            long l = 0;
            if (vector)
            {
                while (l < 4096)
                {
                    uint mask = Vector256.Equals(Vector256.Load(b + i + l), Vector256.Load(b + j + l)).ExtractMostSignificantBits();
                    if (mask != uint.MaxValue)
                    {
                        l += BitOperations.TrailingZeroCount(~mask);
                        break;
                    }
                    l += 32;
                }
            }
            else
            {
                while (l < 4096)
                {
                    ulong d = Unsafe.ReadUnaligned<ulong>(b + i + l) ^ Unsafe.ReadUnaligned<ulong>(b + j + l);
                    if (d != 0)
                    {
                        l += BitOperations.TrailingZeroCount(d) >> 3;
                        break;
                    }
                    l += 8;
                }
            }
            total += l;
        }
        return total;
    }
}
