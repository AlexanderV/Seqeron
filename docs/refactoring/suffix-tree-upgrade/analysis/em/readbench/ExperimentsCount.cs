#if EXPERIMENTS
using System.Globalization;
using SuffixTree.Esa;
using SuffixTree.Esa.External;

namespace ReadBench;

/// <summary>
/// Hypothesis 5 of the diagnosis of 1b.9a (access order): how many memory reads, distinct 64-byte cache lines and distinct 4 KB pages does one pattern descent touch, per array, in the in-memory layout (tagged child table, byte text)
/// and in the layouts of the container (packed child table, packed text)? The descent is replayed over the data of the mapped file with every read recorded by its virtual address (array number and byte offset); nothing is timed.
/// usage: ReadBench count &lt;container&gt; &lt;text&gt;
/// </summary>
internal static unsafe class ExperimentsCount
{
    private sealed class Recorder
    {
        internal readonly HashSet<long> Lines = [];
        internal readonly HashSet<long> Pages = [];
        internal readonly HashSet<long>[] LinesByArray = [[], [], [], []];
        internal readonly long[] ReadsByArray = new long[4];

        internal void Touch(int array, long offset, int length = 1)
        {
            long a = ((long)array << 44) + offset;
            long last = ((long)array << 44) + offset + length - 1;
            for (long line = a >> 6; line <= last >> 6; line++)
            {
                Lines.Add(line);
                LinesByArray[array].Add(line);
            }
            for (long page = a >> 12; page <= last >> 12; page++) Pages.Add(page);
            ReadsByArray[array]++;
        }

        internal void Clear()
        {
            Lines.Clear();
            Pages.Clear();
            foreach (var h in LinesByArray) h.Clear();
        }
    }

    internal static int Run(string file, string textPath)
    {
        string text = File.ReadAllText(textPath, System.Text.Encoding.Latin1);
        SesaFile f = SesaFile.Open(file, new EsaOpenOptions { Access = EsaAccessMode.Trusted });
        EsaIndex index = EsaIndex.FromFile(f);
        int n = index.Length;
        int* sa = index.Sa.Pointer;
        int* cldPacked = (int*)(f.Data + f.Parsed.Child.Offset);
        byte* lcp8 = index.Lcp.MappedBase;
        CompactLcp lcp = index.Lcp;
        int bits = f.TextReader.Bits;
        byte[] bytes = new byte[n + 64];
        f.TextReader.ReadSymbols(0, bytes.AsSpan(0, n));
        SetText(bytes);
        ChildTable t = ChildTable.Build(ChildLayoutKind.Tagged, n, lcp, out _);
        object boxed = t.Tagged;
        var tagged = (uint[])typeof(TaggedChild).GetField("_cld", System.Reflection.BindingFlags.NonPublic | System.Reflection.BindingFlags.Instance)!.GetValue(boxed)!;
        char[] alphabet = index.TextView.Alphabet;
        var map = new int[alphabet[^1] + 1];
        Array.Fill(map, -1);
        for (int s = 0; s < alphabet.Length; s++) map[alphabet[s]] = s;

        var rng = new SplitMix(1234);
        var sets = new (string Name, string[] Patterns)[]
        {
            ("count20 (present)", [.. Enumerable.Range(0, 20000).Select(_ => text.Substring(rng.NextInt(text.Length - 20), 20))]),
            ("locate8 (descent only)", [.. Enumerable.Range(0, 20000).Select(_ => text.Substring(rng.NextInt(text.Length - 8), 8))]),
            ("absent24", [.. Enumerable.Range(0, 20000).Select(_ => new string([.. Enumerable.Range(0, 24).Select(_ => alphabet[rng.NextInt(alphabet.Length)])]))]),
        };
        (string Name, bool TaggedChild, bool BytesText)[] configs =
        [
            ("in memory (tagged child, byte text)", true, true),
            ("container as written (packed child, packed text)", false, false),
            ("container + tagged child (packed text)", true, false),
            ("container + byte text (packed child)", false, true),
        ];
        Console.WriteLine("set\tlayout\treads/query\tlines/query\tpages/query\tlines sa\tlines child\tlines lcp\tlines text\treads sa\treads child\treads lcp\treads text");
        foreach ((string setName, string[] patterns) in sets)
        {
            foreach ((string name, bool tg, bool bt) in configs)
            {
                var rec = new Recorder();
                double lines = 0, pages = 0, reads = 0;
                double[] byArray = new double[4];
                double[] readsBy = new double[4];
                foreach (string p in patterns)
                {
                    rec.Clear();
                    Array.Clear(rec.ReadsByArray);
                    var pat = new byte[p.Length];
                    for (int i = 0; i < p.Length; i++) pat[i] = (byte)map[p[i]];
                    fixed (uint* tg32 = tagged)
                        Descend(pat, n, sa, cldPacked, tg32, lcp8, lcp, bits, tg, bt, rec);
                    lines += rec.Lines.Count;
                    pages += rec.Pages.Count;
                    reads += rec.ReadsByArray.Sum();
                    for (int k = 0; k < 4; k++)
                    {
                        byArray[k] += rec.LinesByArray[k].Count;
                        readsBy[k] += rec.ReadsByArray[k];
                    }
                }
                double q = patterns.Length;
                Console.WriteLine(string.Join('\t', setName, name, F(reads / q), F(lines / q), F(pages / q), F(byArray[0] / q), F(byArray[1] / q), F(byArray[2] / q), F(byArray[3] / q), F(readsBy[0] / q), F(readsBy[1] / q), F(readsBy[2] / q), F(readsBy[3] / q)));
            }
        }
        return 0;
    }

    private static string F(double v) => v.ToString("F2", CultureInfo.InvariantCulture);

    private static bool Descend(byte[] pat, int n, int* sa, int* cldPacked, uint* cldTagged, byte* lcp8, CompactLcp lcp, int bits, bool taggedChild, bool bytesText, Recorder r)
    {
        const uint ValueMask = (1u << 30) - 1;
        const uint NextTag = 1u << 30;
        int m = pat.Length;
        // the arrays: 0 sa, 1 child, 2 lcp, 3 text
        int Sa(int i) { r.Touch(0, 4L * i, 4); return sa[i]; }
        int Lcp(int i) { r.Touch(2, i); return lcp[i]; }
        int Cld(int i) { r.Touch(1, 4L * i, 4); return taggedChild ? (int)(cldTagged[i] & ValueMask) : cldPacked[i]; }
        int Sym(int p) { if (bytesText) r.Touch(3, p); else r.Touch(3, ((long)p * bits) >> 3, 2); return 0; }
        void Range(int start, int length) { if (bytesText) r.Touch(3, start, length); else r.Touch(3, ((long)start * bits) >> 3, (int)((((long)(start + length) * bits) + 7) >> 3) - (int)(((long)start * bits) >> 3) + 1); }
        int Next(int k, int rb, int ell)
        {
            if (taggedChild)
            {
                r.Touch(1, 4L * k, 4);
                uint c = cldTagged[k];
                return (c & ~ValueMask) == NextTag ? (int)(c & ValueMask) : -1;
            }
            if (k >= rb) return -1;
            int lk = Lcp(k);
            if (Lcp(k + 1) == lk) return k + 1;
            int v = Cld(k);
            return Lcp(v) == lk ? v : -1;
        }

        // the dense text for the comparison: read it from the byte copy kept by the caller through the same descent logic (symbols are needed for the control flow)
        int lb = 0, rb = n, ell = 0, c = 0;
        bool isLast = true;
        while (c < m)
        {
            int symbol = pat[c];
            int a = lb;
            int boundary = Cld(isLast ? lb : rb);
            int b;
            while (true)
            {
                b = boundary >= 0 ? boundary - 1 : rb;
                int p = Sa(a) + ell;
                int key = (uint)p < (uint)n ? SymbolAt(p, r, bytesText, bits) : -1;
                if (key == symbol) break;
                if (key > symbol || boundary < 0) return false;
                a = boundary;
                boundary = Next(boundary, rb, ell);
            }
            bool childIsLast = boundary < 0;
            if (a == b)
            {
                int start = Sa(a) + c;
                if ((uint)start > (uint)n || m - c > n - start) return false;
                Range(start, m - c);
                return CompareBytes(start, pat, c, m - c);
            }
            int first = Cld(childIsLast ? a : b);
            int childEll = Lcp(first);
            int limit = Math.Min(childEll, m);
            int from = Sa(a) + c;
            if (limit > c + 1)
            {
                Range(from + 1, limit - c - 1);
                if (!CompareBytes(from + 1, pat, c + 1, limit - c - 1)) return false;
            }
            c = limit;
            lb = a;
            rb = b;
            ell = childEll;
            isLast = childIsLast;
        }
        return true;
    }

    private static byte[]? s_text;

    private static int SymbolAt(int p, Recorder r, bool bytesText, int bits)
    {
        if (bytesText) r.Touch(3, p); else r.Touch(3, ((long)p * bits) >> 3, 2);
        return s_text![p];
    }

    private static bool CompareBytes(int start, byte[] pat, int from, int length) => s_text.AsSpan(start, length).SequenceEqual(pat.AsSpan(from, length));

    internal static void SetText(byte[] bytes) => s_text = bytes;

    private sealed class SplitMix(ulong seed)
    {
        private ulong _state = seed;

        internal int NextInt(int bound)
        {
            unchecked
            {
                _state += 0x9E3779B97F4A7C15UL;
                ulong z = _state;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
                return (int)((z ^ (z >> 31)) % (ulong)bound);
            }
        }
    }
}
#endif
