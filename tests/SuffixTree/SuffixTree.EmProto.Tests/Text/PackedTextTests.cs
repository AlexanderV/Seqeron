using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;

namespace SuffixTree.EmProto.Tests.Text;

[TestFixture]
public class PackedTextTests
{
    private static readonly int[] AllBits = [1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 16, 22, 31, 32];

    private static int[] RandomSymbols(int n, int bits, ulong seed, int repeatEvery = 0)
    {
        var rng = new SplitMix64(seed);
        var s = new int[n];
        ulong mask = bits == 32 ? uint.MaxValue : (1UL << bits) - 1;
        for (int i = 0; i < n; i++) s[i] = (int)(uint)(rng.Next() & mask);
        if (repeatEvery > 0)
        {
            // plant repeats so that the LCE has long matches
            for (int i = repeatEvery; i < n; i++)
            {
                if ((i / repeatEvery) % 3 != 0) s[i] = s[i - repeatEvery];
            }
        }
        return s;
    }

    private static int NaiveLce(int[] s, int i, int j, long max)
    {
        int k = 0;
        while (k < max && s[i + k] == s[j + k]) k++;
        return k;
    }

    [Test]
    public void SymbolAndWindow_RoundTrip_ForAllWidths()
    {
        foreach (int bits in AllBits)
        {
            int[] s = RandomSymbols(300, bits, (ulong)bits);
            using PackedText t = PackedText.Pack(s, bits);
            for (int i = 0; i < s.Length; i++) Assert.That((int)t.Symbol(i), Is.EqualTo(s[i]), $"b={bits} i={i}");
            for (int i = s.Length; i <= s.Length + 1; i++) Assert.That(t.Symbol(i), Is.EqualTo(0u), $"padding reads as symbol 0, b={bits}");
        }
    }

    [Test]
    public void Key_OrderIsTheLexicographicOrderOfTheFirstSymbols_PaddedWithZero()
    {
        foreach (int bits in new[] { 2, 3, 4, 8, 13 })
        {
            int[] s = RandomSymbols(200, bits, 77 + (ulong)bits, repeatEvery: 40);
            using PackedText t = PackedText.Pack(s, bits);
            int ks = t.KeySymbols;
            int Padded(int i, int k) => i + k < s.Length ? s[i + k] : 0;
            for (int i = 0; i < s.Length; i += 3)
            {
                for (int j = 0; j < s.Length; j += 5)
                {
                    int cmp = 0;
                    for (int k = 0; k < ks && cmp == 0; k++) cmp = Padded(i, k).CompareTo(Padded(j, k));
                    Assert.That(t.Key(i).CompareTo(t.Key(j)), Is.EqualTo(cmp), $"b={bits} i={i} j={j}");
                }
            }
        }
    }

    [Test]
    public void Lce_ScalarAndVector_EqualTheNaiveCount()
    {
        foreach (int bits in AllBits)
        {
            int n = 1500;
            int[] s = RandomSymbols(n, bits, 1000 + (ulong)bits, repeatEvery: 211);
            using PackedText t = PackedText.Pack(s, bits);
            var rng = new SplitMix64(5 + (ulong)bits);
            for (int round = 0; round < 400; round++)
            {
                int i = rng.NextInt(n);
                int j = rng.NextInt(n);
                if (round % 5 == 0) j = Math.Min(n - 1, i + 211 * (1 + rng.NextInt(3)));
                if (round % 7 == 0) i = n - 1 - rng.NextInt(40);
                long max = n - Math.Max(i, j);
                long cap = rng.NextInt(4) == 0 ? max : Math.Min(max, 1 + rng.NextInt(900));
                int expected = NaiveLce(s, i, j, cap);
                long steps = 0;
                t.UseVectors = false;
                Assert.That(t.Lce(i, j, cap, ref steps), Is.EqualTo(expected), $"scalar b={bits} i={i} j={j} cap={cap}");
                if (System.Runtime.Intrinsics.X86.Avx2.IsSupported)
                {
                    t.UseVectors = true;
                    Assert.That(t.Lce(i, j, cap, ref steps), Is.EqualTo(expected), $"vector b={bits} i={i} j={j} cap={cap}");
                }
            }
        }
    }

    [Test]
    public void Lce_OnAllSameText_RunsToTheEnd()
    {
        foreach (int bits in new[] { 1, 2, 3, 8 })
        {
            var s = new int[5000];
            using PackedText t = PackedText.Pack(s, bits);
            long steps = 0;
            Assert.That(t.Lce(0, 1, 4999, ref steps), Is.EqualTo(4999));
            Assert.That(t.Lce(10, 4000, 1000, ref steps), Is.EqualTo(1000));
        }
    }

    [Test]
    public void CompareSuffixes_EqualsTheNaiveSuffixOrder_IncludingTheEmptySuffix()
    {
        foreach (int bits in new[] { 1, 2, 3, 8 })
        {
            int n = 300;
            int[] s = RandomSymbols(n, bits, 31 + (ulong)bits, repeatEvery: 17);
            using PackedText t = PackedText.Pack(s, bits);
            for (int i = 0; i <= n; i += 7)
            {
                for (int j = 0; j <= n; j += 11)
                {
                    if (i == j) continue;
                    int expected = new ReadOnlySpan<int>(s, Math.Min(i, n), n - Math.Min(i, n)).SequenceCompareTo(new ReadOnlySpan<int>(s, Math.Min(j, n), n - Math.Min(j, n)));
                    long steps = 0;
                    int got = t.CompareSuffixes(i, j, 0, ref steps, out _);
                    Assert.That(Math.Sign(got), Is.EqualTo(Math.Sign(expected)), $"b={bits} i={i} j={j}");
                }
            }
        }
    }

    [Test]
    public void BitsFor_IsTheSmallestWidth()
    {
        Assert.That(PackedText.BitsFor(1), Is.EqualTo(1));
        Assert.That(PackedText.BitsFor(2), Is.EqualTo(1));
        Assert.That(PackedText.BitsFor(3), Is.EqualTo(2));
        Assert.That(PackedText.BitsFor(4), Is.EqualTo(2));
        Assert.That(PackedText.BitsFor(5), Is.EqualTo(3));
        Assert.That(PackedText.BitsFor(256), Is.EqualTo(8));
        Assert.That(PackedText.BitsFor(257), Is.EqualTo(9));
    }

    [Test]
    public void Pack_RejectsASymbolThatDoesNotFit()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => PackedText.Pack(new[] { 4 }, 2));
    }
}
