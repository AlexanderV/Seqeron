using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-10b3 (B09 F194): <see cref="RobinHoodMap"/> = the flat <c>robin_hood::unordered_map</c> of robin_hood.h 3.9.1 (the header bundled with Mash 2.3), whose iteration order is the order of the <c>counts32</c> of read-set sketch files.
/// Oracle: TestData/Alignment/robin_hood_order_oracle.tsv (scripts/alignment/generate_robin_hood_order_oracle.py compiles the real header with g++ and runs 150 deterministic operation streams): the size and the FNV-1a hash of the
/// (key, value) pairs in iteration order after a bottom-s MinHashHeap stream (32- and 64-bit keys) and after random inserts / erases.
/// </summary>
[TestFixture]
public class RobinHoodMap_Oracle_Tests
{
    private static ulong state;

    private static ulong Next64()
    {
        ulong z = state += 0x9E3779B97F4A7C15UL;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
        return z ^ (z >> 31);
    }

    private static (int Size, ulong Fnv) Run(int bits, int mode, ulong seed, long n, long param)
    {
        state = seed;
        var m = new RobinHoodMap();
        var ordered = new SortedSet<ulong>();
        ulong Mask(ulong x) => bits == 32 ? (uint)x : x;
        for (long i = 0; i < n; i++)
        {
            if (mode == 0)
            {
                ulong h = Mask(Next64() % (ulong)(param * 4));
                if (ordered.Count < param || h < ordered.Max)
                {
                    if (!m.Contains(h)) { m.Insert(h, 1); ordered.Add(h); } else m.Insert(h, 1);
                    if (m.Count > param) { ulong top = ordered.Max; m.Erase(top); ordered.Remove(top); }
                }
            }
            else
            {
                ulong r = Next64();
                ulong h = Mask((r >> 8) % (ulong)param * 0x9E3779B97F4A7C15UL);
                if (r % 3 == 0) m.Erase(h);
                else if (m.Contains(h)) m.Insert(h, 1);
                else m.Insert(h, (uint)(r & 0xFF));
            }
        }

        ulong fnv = 1469598103934665603UL;
        var keys = m.KeysInIterationOrder();
        var values = m.ValuesInIterationOrder();
        for (int i = 0; i < keys.Count; i++)
        {
            fnv = (fnv ^ keys[i]) * 1099511628211UL;
            fnv = (fnv ^ values[i]) * 1099511628211UL;
        }

        return (m.Count, fnv);
    }

    [Test]
    public void IterationOrder_EqualsRobinHoodHeader_150Streams()
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.robin_hood_order_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name) ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        int rows = 0, bits32 = 0, mode1 = 0;
        var failures = new List<string>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            int bits = int.Parse(f[0], CultureInfo.InvariantCulture), mode = int.Parse(f[1], CultureInfo.InvariantCulture);
            var (size, fnv) = Run(bits, mode, ulong.Parse(f[2], CultureInfo.InvariantCulture), long.Parse(f[3], CultureInfo.InvariantCulture), long.Parse(f[4], CultureInfo.InvariantCulture));
            rows++;
            if (bits == 32) bits32++;
            if (mode == 1) mode1++;
            if (size != int.Parse(f[5], CultureInfo.InvariantCulture) || fnv != ulong.Parse(f[6], CultureInfo.InvariantCulture))
                failures.Add($"row {rows}: {line} -> {size} {fnv}");
        }

        Assert.That(failures, Is.Empty, $"{failures.Count} of {rows}:\n" + string.Join("\n", failures.Take(5)));
        Assert.That(rows, Is.EqualTo(150));
        Assert.That(bits32, Is.GreaterThan(40));
        Assert.That(mode1, Is.GreaterThan(30));
    }

    [Test]
    public void HashInt_WorkedValue()
    {
        // robin_hood::hash_int(1): h1 = 0xA24BAED4963EE407, h2 = rotr(1, 32) * 0x9FB21C651E98DF25, result rotr(h1 + h2, 32)
        ulong h1 = 0xA24BAED4963EE407UL;
        ulong h2 = ulong.RotateRight(1UL, 32) * 0x9FB21C651E98DF25UL;
        Assert.That(RobinHoodMap.HashInt(1), Is.EqualTo(ulong.RotateRight(h1 + h2, 32)));
    }
}
