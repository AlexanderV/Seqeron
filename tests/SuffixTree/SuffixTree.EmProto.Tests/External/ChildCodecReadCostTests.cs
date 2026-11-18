using System.Diagnostics;
using System.Globalization;
using System.Reflection;
using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>
/// G5 of work package 1b.7c (reported, not a gate; the owner's open item of 2026-11-13): what reading the child table through the Packed codec costs against Tagged (valid below 2^30 symbols), on the files the
/// production stage wrote at 100 Mbp. The production <c>.cld</c> file is loaded into the existing in-RAM <c>PackedChild</c> / <c>ChildTable</c> / <c>EsaIndex</c> path (the private constructors, through reflection: the in-memory
/// structures stay unchanged), the Tagged table is built in memory from the same suffix array and LCP (<c>EsaIndex.FromArrays</c>), and the same query mix runs on both: a depth-first walk of the child table through
/// <c>FirstLIndex</c> / <c>NextLIndex</c> (the child steps), counts of 20-mers that occur and that do not, locates of 12-mers. Driven by <c>analysis/em/run_1b7c_g5.sh</c>; not part of any run (<c>[Explicit]</c>).
/// </summary>
[TestFixture]
[Explicit("G5 of 1b.7c: run by analysis/em/run_1b7c_g5.sh, which sets the G5_* variables")]
[Category("G5Cost")]
public class ChildCodecReadCostTests
{
    private static string Env(string name) => Environment.GetEnvironmentVariable(name) ?? throw new InvalidOperationException($"{name} is not set");

    private static EsaIndex WrapFile(EsaIndex donor, uint[] cld)
    {
        var packedCtor = typeof(PackedChild).GetConstructor(BindingFlags.NonPublic | BindingFlags.Instance, [typeof(int[]), typeof(CompactLcp)])!;
        object packed = packedCtor.Invoke([MemoryMarshal.Cast<uint, int>(cld).ToArray(), donor.Lcp]);
        var tableCtor = typeof(ChildTable).GetConstructor(BindingFlags.NonPublic | BindingFlags.Instance, [typeof(ChildLayoutKind), typeof(Reference3Child), typeof(PackedChild), typeof(TaggedChild)])!;
        var table = (ChildTable)tableCtor.Invoke([ChildLayoutKind.Packed, default(Reference3Child), packed, default(TaggedChild)]);
        var indexCtor = typeof(EsaIndex).GetConstructor(BindingFlags.NonPublic | BindingFlags.Instance, [typeof(string), typeof(PreparedText), typeof(int[]), typeof(CompactLcp), typeof(ChildTable), typeof(EsaStatistics)])!;
        return (EsaIndex)indexCtor.Invoke([donor.Text, donor.Prepared, donor.SuffixArray, donor.Lcp, table, donor.Statistics]);
    }

    /// <summary>Depth-first walk of the lcp-interval tree through the child table, at most <paramref name="budget"/> child steps; returns a checksum of the l-indices met.</summary>
    private static long Walk(ChildTable child, int n, long budget, out long steps)
    {
        var stack = new Stack<(int A, int B, bool Last)>();
        stack.Push((0, n, true));
        long sum = 0, count = 0;
        while (stack.Count > 0 && count < budget)
        {
            (int a, int b, bool last) = stack.Pop();
            int first = child.FirstLIndex(a, b, last);
            int start = a;
            int k = first;
            while (k >= 0 && count < budget)
            {
                count++;
                sum += k;
                if (k - 1 > start) stack.Push((start, k - 1, false));
                start = k;
                k = child.NextLIndex(k, b);
            }
            if (b > start) stack.Push((start, b, true));
        }
        steps = count;
        return sum;
    }

    [Test]
    public void RunOnce()
    {
        string textPath = Env("G5_TEXT");
        string saPath = Env("G5_SA");
        string lcpPrefix = Env("G5_LCP");
        string cldPath = Env("G5_CLD");
        string rows = Env("G5_ROWS");
        string label = Env("G5_LABEL");
        int passes = int.Parse(Environment.GetEnvironmentVariable("G5_PASSES") ?? "4", CultureInfo.InvariantCulture);
        long walkBudget = long.Parse(Environment.GetEnvironmentVariable("G5_WALK") ?? "40000000", CultureInfo.InvariantCulture);
        int queries = int.Parse(Environment.GetEnvironmentVariable("G5_QUERIES") ?? "300000", CultureInfo.InvariantCulture);

        byte[] raw = File.ReadAllBytes(textPath);
        int n = raw.Length;
        char[] letters = [.. raw.Distinct().OrderBy(b => b).Select(b => (char)b)];
        var rank = new byte[256];
        for (int i = 0; i < letters.Length; i++) rank[letters[i]] = (byte)i;
        var dense = new byte[n];
        for (int i = 0; i < n; i++) dense[i] = rank[raw[i]];
        string text = string.Create(n, raw, static (span, bytes) => { for (int i = 0; i < span.Length; i++) span[i] = (char)bytes[i]; });
        var prepared = new PreparedText(letters, dense, null);
        uint[] saU = MemoryMarshal.Cast<byte, uint>(File.ReadAllBytes(saPath)).ToArray();
        int[] sa = MemoryMarshal.Cast<uint, int>(saU).ToArray();
        byte[] lcpBytes = File.ReadAllBytes(lcpPrefix + ".lcp");
        uint[] ovf = MemoryMarshal.Cast<byte, uint>(File.ReadAllBytes(lcpPrefix + ".ovf")).ToArray();
        var ranks = new int[ovf.Length / 2];
        var values = new int[ovf.Length / 2];
        for (int i = 0; i < ranks.Length; i++)
        {
            ranks[i] = (int)ovf[2 * i];
            values[i] = (int)ovf[(2 * i) + 1];
        }
        Assert.That(lcpBytes, Has.Length.EqualTo(n + 1), "the test reads a U8 LCP (random DNA)");
        CompactLcp lcp = CompactLcp.FromParts(LcpKind.U8, n + 1, lcpBytes, null, null, ranks, values);
        uint[] cld = MemoryMarshal.Cast<byte, uint>(File.ReadAllBytes(cldPath)).ToArray();

        EsaIndex tagged = EsaIndex.FromArrays(text, prepared, sa, lcp, ChildLayoutKind.Tagged);
        EsaIndex packedMem = EsaIndex.FromArrays(text, prepared, sa, lcp, ChildLayoutKind.Packed);
        EsaIndex packedFile = WrapFile(tagged, cld);
        uint[] memArray = (uint[])(typeof(PackedChild).GetField("_cld", BindingFlags.NonPublic | BindingFlags.Instance)!.GetValue(packedMem.Child.Packed)!);
        Assert.That(memArray.AsSpan().SequenceEqual(cld), Is.True, "the production file is the in-RAM packed table, entry by entry");

        var rng = new SplitMix64(2026_11_18);
        string[] present = new string[queries];
        string[] absent = new string[queries];
        string[] twelve = new string[queries / 3];
        for (int i = 0; i < queries; i++)
        {
            present[i] = text.Substring(rng.NextInt(n - 20), 20);
            var a = new char[20];
            for (int k = 0; k < 20; k++) a[k] = letters[rng.NextInt(letters.Length)];
            absent[i] = new string(a);
        }
        for (int i = 0; i < twelve.Length; i++) twelve[i] = text.Substring(rng.NextInt(n - 12), 12);

        (string Name, EsaIndex Index)[] variants = [("packed", packedFile), ("tagged", tagged)];
        if (Environment.GetEnvironmentVariable("G5_FIRST") == "tagged") variants = [variants[1], variants[0]];
        long[] reference = new long[3];
        for (int pass = 0; pass < passes; pass++)
        {
            foreach ((string name, EsaIndex index) in pass % 2 == 0 ? variants : [variants[1], variants[0]])
            {
                long t0 = Stopwatch.GetTimestamp();
                long walk = Walk(index.Child, n, walkBudget, out long steps);
                long t1 = Stopwatch.GetTimestamp();
                long countPresent = 0;
                foreach (string p in present) countPresent += index.CountOccurrences(p);
                long t2 = Stopwatch.GetTimestamp();
                long countAbsent = 0;
                foreach (string p in absent) countAbsent += index.CountOccurrences(p);
                long t3 = Stopwatch.GetTimestamp();
                long located = 0;
                foreach (string p in twelve) located += index.FindAllOccurrences(p).Count;
                long t4 = Stopwatch.GetTimestamp();
                long[] got = [walk + countPresent, countAbsent, located];
                if (pass == 0 && reference[0] == 0) reference = got;
                Assert.That(got, Is.EqualTo(reference), $"{name} pass {pass}: the same answers on both codecs");
                double f = 1.0 / Stopwatch.Frequency;
                string row = string.Join('\t', label, name, pass == 0 ? "cold" : "warm", pass,
                    ((t1 - t0) * f).ToString("F4", CultureInfo.InvariantCulture), steps, ((t2 - t1) * f).ToString("F4", CultureInfo.InvariantCulture), ((t3 - t2) * f).ToString("F4", CultureInfo.InvariantCulture), ((t4 - t3) * f).ToString("F4", CultureInfo.InvariantCulture));
                File.AppendAllText(rows, row + "\n");
                TestContext.Out.WriteLine(row);
            }
        }
    }
}
