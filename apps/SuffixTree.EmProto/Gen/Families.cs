using SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Gen;

/// <summary>A text as dense order-preserving symbols (<c>0..Sigma-1</c>) with the ASCII letters of its alphabet (<c>Letters[symbol]</c>).</summary>
internal sealed record SymbolText(string Name, byte[] Symbols, string Letters)
{
    internal int Sigma => Letters.Length;

    internal long Length => Symbols.Length;
}

/// <summary>
/// The deterministic text families of the experiments (ADR-004 section 5): the adversarial families of the 1a.3 / 1a.7 generators (same definitions as
/// <c>StringFamilies</c> of the ESA test project, so that the committed libsais hashes apply), the 171-bp tandem repeat, a text with a long N run, random DNA and a genome-like family
/// (random background, interspersed repeats with 10-15 % divergence, 171-bp satellite arrays, N runs). Generators are pure functions of their arguments; randomness is <see cref="SplitMix64"/>.
/// </summary>
internal static class Families
{
    private const string Dna = "ACGT";

    internal static SymbolText AllSame(long n) => new("all-same", new byte[n], "a");

    /// <summary>Fibonacci word f1 = "a", f2 = "ab", f(k) = f(k-1) f(k-2), truncated to n.</summary>
    internal static SymbolText Fibonacci(long n)
    {
        var prev = new List<byte> { 0 };
        var cur = new List<byte> { 0, 1 };
        while (cur.Count < n)
        {
            var next = new List<byte>(cur.Count + prev.Count);
            next.AddRange(cur);
            next.AddRange(prev);
            prev = cur;
            cur = next;
        }
        return new SymbolText("fibonacci", cur.GetRange(0, (int)n).ToArray(), "ab");
    }

    internal static SymbolText ThueMorse(long n)
    {
        var s = new byte[n];
        for (long i = 0; i < n; i++) s[i] = (byte)(long.PopCount(i) & 1);
        return new SymbolText("thue-morse", s, "ab");
    }

    /// <summary>Uniform random over ACGT (<c>letters[rng.NextInt(4)]</c> as in <c>StringFamilies.RandomDna</c>).</summary>
    internal static SymbolText RandomDna(long n, ulong seed, string? name = null)
    {
        var rng = new SplitMix64(seed);
        var s = new byte[n];
        for (long i = 0; i < n; i++) s[i] = (byte)rng.NextInt(4);
        return new SymbolText(name ?? $"random-dna-s{seed}", s, Dna);
    }

    /// <summary>A random DNA unit of <paramref name="unit"/> symbols repeated <paramref name="copies"/> times (171 x 10^4 = 1.71 M symbols is the satellite-like adversary of E1).</summary>
    internal static SymbolText Tandem(int unit, int copies, ulong seed)
    {
        var rng = new SplitMix64(seed);
        var u = new byte[unit];
        for (int i = 0; i < unit; i++) u[i] = (byte)rng.NextInt(4);
        var s = new byte[(long)unit * copies];
        for (long i = 0; i < s.Length; i++) s[i] = u[i % unit];
        return new SymbolText($"tandem-{unit}x{copies}", s, Dna);
    }

    /// <summary>Random DNA with one run of N of <paramref name="runLength"/> symbols starting at <paramref name="runStart"/> (the alphabet becomes ACGNT: T moves to code 4).</summary>
    internal static SymbolText NRun(long n, long runStart, long runLength, ulong seed)
    {
        var rng = new SplitMix64(seed);
        var s = new byte[n];
        for (long i = 0; i < n; i++)
        {
            byte v = (byte)rng.NextInt(4);
            s[i] = v == 3 ? (byte)4 : v;
        }
        for (long i = runStart; i < Math.Min(n, runStart + runLength); i++) s[i] = 3;
        return new SymbolText($"nrun-{runLength}in{n}", s, "ACGNT");
    }

    /// <summary>
    /// Genome-like text over ACGNT (codes A0 C1 G2 N3 T4): segments are drawn until n symbols exist, in proportions of length of about 55 % random background, 40 % interspersed repeats
    /// (a copy of a random part of one of 24 family consensus sequences of 300-6000 bp, every position replaced with probability 10-15 %, 1 % indels) and 4 % satellite arrays (171-bp
    /// monomers organised in higher-order repeats of 12 monomers that differ by 25 % from each other; the array copies differ by 1 %); about 1 % of the length are six runs of N, the largest n/200.
    /// </summary>
    internal static SymbolText GenomeLike(long n, ulong seed)
    {
        var rng = new SplitMix64(seed);
        var s = new byte[n];
        var families = new byte[24][];
        for (int f = 0; f < families.Length; f++)
        {
            families[f] = new byte[300 + rng.NextInt(5700)];
            for (int i = 0; i < families[f].Length; i++) families[f][i] = (byte)rng.NextInt(4);
        }
        var monomer = new byte[171];
        for (int i = 0; i < monomer.Length; i++) monomer[i] = (byte)rng.NextInt(4);
        var nRuns = new List<(long Start, long Length)>
        {
            (n / 20, n / 200),
            (n / 3, n / 1000),
            (n / 2, n / 2000),
            (n * 3 / 5, n / 1000),
            (n * 4 / 5, n / 2000),
            (n - (n / 2000), n / 2000),
        };
        nRuns.Sort((a, b) => a.Start.CompareTo(b.Start));
        long satMean = Math.Max(2000, n / 400) + (Math.Max(1000, n / 800) / 2);
        double wBackground = 0.55 / 1000, wRepeat = 0.40 / 1500, wSatellite = 0.04 / satMean;
        double total = wBackground + wRepeat + wSatellite;
        long pos = 0;
        int nextRun = 0;
        while (pos < n)
        {
            if (nextRun < nRuns.Count && pos >= nRuns[nextRun].Start)
            {
                long end = Math.Min(n, pos + nRuns[nextRun].Length);
                for (long i = pos; i < end; i++) s[i] = 3;
                pos = end;
                nextRun++;
                continue;
            }
            long limit = nextRun < nRuns.Count ? Math.Min(n, nRuns[nextRun].Start) : n;
            double pick = rng.NextDouble() * total;
            if (pick < wBackground)
            {
                long len = Math.Min(limit - pos, 200 + rng.NextInt(1600));
                for (long i = 0; i < len; i++) s[pos + i] = Code((byte)rng.NextInt(4));
                pos += len;
            }
            else if (pick < wBackground + wRepeat)
            {
                long len = Math.Min(limit - pos, 500 + rng.NextInt(2000));
                RepeatCopy(s, pos, len, families, rng);
                pos += len;
            }
            else
            {
                long len = Math.Min(limit - pos, Math.Max(2000, n / 400) + rng.NextInt((int)Math.Max(1000, n / 800)));
                var hor = new byte[171 * 12];
                for (int m = 0; m < 12; m++)
                {
                    for (int i = 0; i < 171; i++) hor[(m * 171) + i] = rng.NextDouble() < 0.25 ? (byte)rng.NextInt(4) : monomer[i];
                }
                for (long k = 0; k < len; k++)
                {
                    byte v = hor[k % hor.Length];
                    if (rng.NextDouble() < 0.01) v = (byte)rng.NextInt(4);
                    s[pos + k] = Code(v);
                }
                pos += len;
            }
        }
        return new SymbolText($"genome-like-{n / 1_000_000}m-s{seed}", s, "ACGNT");
    }

    /// <summary>Final code of a base 0..3 (A C G T): T moves to 4 so that N (3) sorts between G and T.</summary>
    private static byte Code(byte b) => b == 3 ? (byte)4 : b;

    private static void RepeatCopy(byte[] s, long pos, long length, byte[][] families, SplitMix64 rng)
    {
        byte[] f = families[rng.NextInt(families.Length)];
        double divergence = 0.10 + (0.05 * rng.NextDouble());
        int at = rng.NextInt(Math.Max(1, f.Length - (int)Math.Min(length, f.Length) + 1));
        long written = 0;
        int src = at;
        while (written < length)
        {
            if (src >= f.Length) src = 0;
            double r = rng.NextDouble();
            if (r < 0.005)
            {
                src++; // deletion in the copy
                continue;
            }
            if (r < 0.01)
            {
                s[pos + written++] = Code((byte)rng.NextInt(4)); // insertion
                continue;
            }
            byte v = f[src++];
            if (rng.NextDouble() < divergence) v = (byte)((v + 1 + rng.NextInt(3)) & 3);
            s[pos + written++] = Code(v);
        }
    }

    /// <summary>The ASCII bytes of a text (<c>Letters[symbol]</c>), as in the benchmark inputs.</summary>
    internal static byte[] ToAscii(SymbolText t)
    {
        var map = new byte[t.Sigma];
        for (int i = 0; i < map.Length; i++) map[i] = (byte)t.Letters[i];
        var b = new byte[t.Symbols.Length];
        for (long i = 0; i < b.Length; i++) b[i] = map[t.Symbols[i]];
        return b;
    }

    /// <summary>Reads a benchmark input (one byte per symbol, no terminator) and remaps its distinct byte values to <c>0..sigma-1</c> in increasing byte order.</summary>
    internal static SymbolText FromFile(string path, string name)
    {
        byte[] raw = File.ReadAllBytes(path);
        var present = new bool[256];
        foreach (byte b in raw) present[b] = true;
        var map = new byte[256];
        var letters = new List<char>();
        for (int v = 0; v < 256; v++)
        {
            if (!present[v]) continue;
            map[v] = (byte)letters.Count;
            letters.Add((char)v);
        }
        if (letters.Count > 256) throw new InvalidOperationException("impossible");
        for (long i = 0; i < raw.Length; i++) raw[i] = map[raw[i]];
        return new SymbolText(name, raw, new string([.. letters]));
    }
}
