using System.Diagnostics;
using System.Globalization;
using System.Runtime.CompilerServices;
using SuffixTree.Esa.Tests.TestSupport;
using Prod = SuffixTree.Esa.External;
using Proto = SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>
/// Gate G3 of work package 1b.7b (cost parity of the production <c>MappedTextReader</c> with the prototype <c>PackedText</c>): random key extraction and common prefixes on a packed file, both variants in one process, the loops
/// compiled with full optimisation, every variant timed in turn for <c>G3_REPS</c> passes (the median is the figure), warm page cache (the file is read once before). Driven by
/// <c>docs/refactoring/suffix-tree-upgrade/analysis/em/run_1b7b_g3.sh</c>; not part of any run (<c>[Explicit]</c>).
/// </summary>
[TestFixture]
[Explicit("gate G3 of 1b.7b: run by analysis/em/run_1b7b_g3.sh, which sets the G3_* variables")]
[Category("G3Cost")]
public class TextCostTests
{
    private static string Env(string name) => Environment.GetEnvironmentVariable(name) ?? throw new InvalidOperationException($"{name} is not set");

    private static string EnvOr(string name, string fallback) => Environment.GetEnvironmentVariable(name) ?? fallback;

    [MethodImpl(MethodImplOptions.NoInlining | MethodImplOptions.AggressiveOptimization)]
    private static ulong ProtoKeys(Proto.PackedText t, uint[] pos)
    {
        ulong acc = 0;
        foreach (uint p in pos) acc ^= t.Key(p);
        return acc;
    }

    [MethodImpl(MethodImplOptions.NoInlining | MethodImplOptions.AggressiveOptimization)]
    private static ulong ProdKeys(Prod.MappedTextReader t, uint[] pos)
    {
        ulong acc = 0;
        foreach (uint p in pos) acc ^= t.Key(p);
        return acc;
    }

    [MethodImpl(MethodImplOptions.NoInlining | MethodImplOptions.AggressiveOptimization)]
    private static ulong ProtoSymbols(Proto.PackedText t, uint[] pos)
    {
        ulong acc = 0;
        foreach (uint p in pos) acc += t.Symbol(p);
        return acc;
    }

    [MethodImpl(MethodImplOptions.NoInlining | MethodImplOptions.AggressiveOptimization)]
    private static ulong ProdSymbols(Prod.MappedTextReader t, uint[] pos)
    {
        ulong acc = 0;
        foreach (uint p in pos) acc += t.Symbol(p);
        return acc;
    }

    [MethodImpl(MethodImplOptions.NoInlining | MethodImplOptions.AggressiveOptimization)]
    private static ulong ProtoLce(Proto.PackedText t, uint[] a, uint[] b, int[] cap)
    {
        ulong acc = 0;
        long steps = 0;
        for (int k = 0; k < a.Length; k++) acc += (ulong)t.Lce(a[k], b[k], cap[k], ref steps);
        return acc + (ulong)steps;
    }

    [MethodImpl(MethodImplOptions.NoInlining | MethodImplOptions.AggressiveOptimization)]
    private static ulong ProdLce(Prod.MappedTextReader t, uint[] a, uint[] b, int[] cap)
    {
        ulong acc = 0;
        long steps = 0;
        for (int k = 0; k < a.Length; k++) acc += (ulong)t.Lce(a[k], b[k], cap[k], ref steps);
        return acc + (ulong)steps;
    }

    private static void Measure(string label, string what, long ops, int reps, string rows, (string Name, Func<ulong> Body)[] variants)
    {
        var times = variants.ToDictionary(v => v.Name, _ => new List<double>());
        var sums = new Dictionary<string, ulong>();
        // warm-up by time: the runtime compiles a method at tier 0 first and promotes it after 30 calls once no new tier-0 method has appeared for 100 ms, so a few passes of a few milliseconds are not enough
        double warmSeconds = double.Parse(EnvOr("G3_WARM_S", "2"), CultureInfo.InvariantCulture);
        foreach ((string name, Func<ulong> body) in variants)
        {
            long w0 = Stopwatch.GetTimestamp();
            do
            {
                sums[name] = body();
            }
            while ((Stopwatch.GetTimestamp() - w0) / (double)Stopwatch.Frequency < warmSeconds);
        }
        for (int r = 0; r < reps; r++)
        {
            for (int k = 0; k < variants.Length; k++)
            {
                (string name, Func<ulong> body) = variants[(k + r) % variants.Length];
                long t0 = Stopwatch.GetTimestamp();
                ulong sum = body();
                double ns = (Stopwatch.GetTimestamp() - t0) * 1e9 / Stopwatch.Frequency / ops;
                times[name].Add(ns);
                Assert.That(sum, Is.EqualTo(sums[name]), $"{name}: the checksum changes between passes");
            }
        }
        foreach ((string name, _) in variants)
        {
            List<double> t = times[name];
            t.Sort();
            string row = string.Join('\t', label, what, name, ops.ToString(CultureInfo.InvariantCulture), t[t.Count / 2].ToString("F3", CultureInfo.InvariantCulture),
                t[0].ToString("F3", CultureInfo.InvariantCulture), t[^1].ToString("F3", CultureInfo.InvariantCulture), string.Join(',', t.Select(x => x.ToString("F3", CultureInfo.InvariantCulture))), sums[name].ToString(CultureInfo.InvariantCulture));
            File.AppendAllText(rows, row + "\n");
            TestContext.Out.WriteLine(row);
        }
        foreach ((string name, _) in variants)
        {
            if (name.StartsWith("proto", StringComparison.Ordinal)) Assert.That(sums[name.Replace("proto", "prod", StringComparison.Ordinal)], Is.EqualTo(sums[name]), $"{name}: the production and the prototype compute the same result ({what}), steps included");
        }
    }

    [Test]
    public unsafe void RandomKeysAndSymbols()
    {
        string pk = Env("G3_PK");
        long n = long.Parse(Env("G3_N"), CultureInfo.InvariantCulture);
        int bits = int.Parse(Env("G3_BITS"), CultureInfo.InvariantCulture);
        string rows = Env("G3_ROWS");
        string label = Env("G3_LABEL");
        int count = int.Parse(EnvOr("G3_KEYS", "8000000"), CultureInfo.InvariantCulture);
        int reps = int.Parse(EnvOr("G3_REPS", "5"), CultureInfo.InvariantCulture);
        var rng = new SplitMix64(1);
        var pos = new uint[count];
        for (int i = 0; i < count; i++) pos[i] = (uint)(rng.Next() % (ulong)n);
        File.ReadAllBytes(pk); // warm page cache
        using var mapped = new Proto.MappedFile(pk, writable: false);
        using var proto = new Proto.PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        using Prod.MappedTextReader prod = Prod.MappedTextReader.Open(pk, n, bits);
        Measure(label, "key", count, reps, rows, [("proto", () => ProtoKeys(proto, pos)), ("prod", () => ProdKeys(prod, pos))]);
        Measure(label, "symbol", count, reps, rows, [("proto", () => ProtoSymbols(proto, pos)), ("prod", () => ProdSymbols(prod, pos))]);
    }

    [Test]
    public unsafe void CommonPrefixes()
    {
        string pk = Env("G3_PK");
        long n = long.Parse(Env("G3_N"), CultureInfo.InvariantCulture);
        int bits = int.Parse(Env("G3_BITS"), CultureInfo.InvariantCulture);
        string rows = Env("G3_ROWS");
        string label = Env("G3_LABEL");
        string kind = Env("G3_PAIRS"); // random: independent positions; near: a short distance apart (long matches on repetitive texts)
        int count = int.Parse(EnvOr("G3_KEYS", "2000000"), CultureInfo.InvariantCulture);
        int reps = int.Parse(EnvOr("G3_REPS", "5"), CultureInfo.InvariantCulture);
        var rng = new SplitMix64(2);
        var a = new uint[count];
        var b = new uint[count];
        var cap = new int[count];
        for (int k = 0; k < count; k++)
        {
            long i = (long)(rng.Next() % (ulong)n);
            long j = kind == "near" ? Math.Min(n - 1, i + 1 + (long)(rng.Next() % 3000UL)) : (long)(rng.Next() % (ulong)n);
            a[k] = (uint)i;
            b[k] = (uint)j;
            cap[k] = (int)Math.Min(100_000, n - Math.Max(i, j));
        }
        File.ReadAllBytes(pk);
        using var mapped = new Proto.MappedFile(pk, writable: false);
        using var proto = new Proto.PackedText(mapped.Pointer, n, bits, mapped, ownsOwner: false);
        using Prod.MappedTextReader prod = Prod.MappedTextReader.Open(pk, n, bits);
        bool avx = System.Runtime.Intrinsics.X86.Avx2.IsSupported;
        Measure(label, "lce-" + kind, count, reps, rows,
        [
            ("proto-vector", () => { proto.UseVectors = avx; return ProtoLce(proto, a, b, cap); }),
            ("prod-vector", () => { prod.UseVectors = true; return ProdLce(prod, a, b, cap); }),
            ("proto-scalar", () => { proto.UseVectors = false; return ProtoLce(proto, a, b, cap); }),
            ("prod-scalar", () => { prod.UseVectors = false; return ProdLce(prod, a, b, cap); }),
        ]);
    }
}
