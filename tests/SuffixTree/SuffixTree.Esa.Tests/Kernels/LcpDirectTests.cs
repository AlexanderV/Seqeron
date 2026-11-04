using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>
/// 1o.3b: the suffix-array-order LCP kernel (<see cref="LcpKernel.BuildDirect{TSym}"/>) beyond the generic registry run: the vector-width boundaries, each reason
/// to give up (and the fallback result), the success path with escapes, the argument checks, and memory safety on corrupt input.
/// </summary>
[TestFixture]
public class LcpDirectTests
{
    private static LcpEngine Engine => KernelRegistry.DirectLcp;

    private static CompactLcp Reference(PreparedText p, int[] sa)
        => p.Width == SymbolWidth.Byte ? LcpReference.Build(p.ByteReader(), sa) : LcpReference.Build(p.UShortReader(), sa);

    /// <summary>Runs the direct attempt alone on a text with its true suffix array; null when it gave up.</summary>
    private static CompactLcp? Attempt(string text, out int[] sa, out PreparedText prepared)
    {
        prepared = TextPreparation.Prepare(text);
        sa = SaisReference.BuildSuffixArray(text);
        using var arena = new NativeArena(ArenaOptions.Guarded());
        int n = prepared.Length;
        Buf<int> saBuf = arena.Allocate<int>(n + 1, zero: false, "sa");
        sa.CopyTo(saBuf.AsSpan(n + 1));
        CompactLcp? result;
        if (prepared.Width == SymbolWidth.Byte)
        {
            Buf<byte> t = arena.Allocate<byte>(n, zero: false, "text");
            prepared.Bytes.CopyTo(t.AsSpan(n));
            result = LcpKernel.TryDirect(t, n, saBuf, arena);
        }
        else
        {
            Buf<ushort> t = arena.Allocate<ushort>(n, zero: false, "text");
            prepared.UShorts.CopyTo(t.AsSpan(n));
            result = LcpKernel.TryDirect(t, n, saBuf, arena);
        }
        arena.VerifyCanaries();
        return result;
    }

    private static void AssertAttempt(string text, bool succeeds)
    {
        CompactLcp? direct = Attempt(text, out int[] sa, out PreparedText p);
        Assert.That(direct is not null, Is.EqualTo(succeeds), $"direct attempt on {text.Length} symbols");
        if (direct is not null)
        {
            CompactLcp reference = Reference(p, sa);
            Assert.That(direct.Kind, Is.EqualTo(reference.Kind));
            Assert.That(direct.OverflowCount, Is.EqualTo(reference.OverflowCount));
            Assert.That(direct.ToArray(), Is.EqualTo(reference.ToArray()));
        }
        KernelHarness.AssertLcp(Engine, text, allPoisons: false); // with or without the attempt, the engine gives the reference result
    }

    [Test]
    public void TinyTexts_n0_To_n3_AndTheTwoEntryPoints()
    {
        foreach (string text in new[] { "", "a", "aa", "ab", "ba", "aaa", "aba", "abc", "cab" })
        {
            KernelHarness.AssertLcp(Engine, text);
            PreparedText p = TextPreparation.Prepare(text);
            int[] sa = SaisReference.BuildSuffixArray(text);
            Assert.That(LcpKernel.Build(LcpKernelKind.Direct, p, sa).ToArray(), Is.EqualTo(NaiveLcp.Build(text, sa)), $"selector on '{text}'");
        }
    }

    [Test]
    public void VectorWidthBoundaries_PeriodicTexts_AtEveryLengthAndPeriod()
    {
        // 16 symbols per vector (bytes) / 8 (16-bit): suffix pairs with 15/16/17 symbols left, matches that end in the first vector, exactly at its end, and in the next one
        for (int period = 1; period <= 20; period++)
            for (int n = 0; n <= 80; n++)
            {
                string s = StringFamilies.PeriodicAaB(n, period + 1);
                KernelHarness.AssertLcp(Engine, s, allPoisons: false);
                KernelHarness.AssertLcp(Engine, new string([.. s.Select(c => (char)(c + 0x200))]), allPoisons: false);
            }
        string distinct = new([.. Enumerable.Range(0, 300).Select(i => (char)(0x100 + i))]);
        Assert.That(TextPreparation.Prepare(distinct + "ab").Width, Is.EqualTo(SymbolWidth.UShort));
        for (int period = 1; period <= 12; period++)
            for (int n = 0; n <= 50; n++)
                KernelHarness.AssertLcp(Engine, distinct + StringFamilies.PeriodicAaB(n, period + 1), allPoisons: false);
    }

    [Test]
    public void RandomDna_SucceedsDirectly_AndLongRepeatsAreResolved()
    {
        AssertAttempt(StringFamilies.RandomDna(100_000, 5), succeeds: true);
        // a block of 1,000 twice inside random filler: the pair gives one LCP entry per offset, 746 escapes (16 * 746 <= 40,000: U8), 984 long ranks (<= n / 4), about 1,000^2 / 32 = 31,000 steps (<= 2 n = 80,000)
        string block = StringFamilies.Random(1000, 4, 7);
        string text = StringFamilies.Random(19_000, 4, 8) + block + StringFamilies.Random(19_000, 4, 9) + block;
        AssertAttempt(text, succeeds: true);
        Assert.That(Attempt(text, out _, out _)!.OverflowCount, Is.InRange(700, 800));
    }

    [Test]
    public void TooManyLongRanks_GiveUp_AndTheFallbackIsTheReference()
    {
        AssertAttempt(StringFamilies.AllSame(1000), succeeds: false);
        AssertAttempt(StringFamilies.AbRepeat(5000), succeeds: false);
        AssertAttempt(StringFamilies.Fibonacci(5000), succeeds: false);
    }

    [Test]
    public void TheLongRankCap_IsExactlyAQuarterOfTheRanks()
    {
        // a^m followed by a random tail: the ranks with an LCP of at least 16 are counted from the reference and compared with n / 4 on both sides of the limit
        int below = 0, above = 0;
        for (int m = 20; m <= 90; m += 7)
            for (int tail = 20; tail <= 400; tail += 13)
            {
                string text = StringFamilies.AllSame(m) + StringFamilies.Random(tail, 3, (ulong)(m * 1000 + tail)).Replace('a', 'b');
                PreparedText p = TextPreparation.Prepare(text);
                int[] sa = SaisReference.BuildSuffixArray(text);
                int[] lcp = Reference(p, sa).ToArray();
                int longs = lcp.Count(v => v >= 16);
                int n = text.Length;
                bool fitsCap = longs <= n / 4;
                CompactLcp? direct = Attempt(text, out _, out _);
                if (fitsCap) { Assert.That(direct, Is.Not.Null, $"m {m} tail {tail}: {longs} long of {n}"); below++; }
                else { Assert.That(direct, Is.Null, $"m {m} tail {tail}: {longs} long of {n}"); above++; }
            }
        Assert.That(below, Is.GreaterThan(5));
        Assert.That(above, Is.GreaterThan(5));
    }

    [Test]
    public void MoreEscapesThanTheU8Rule_GiveUp_AndTheFallbackPicksTheWiderKind()
    {
        // 1,500 symbols copied from the start of a random text of 10,000: 1,245 escapes against a limit of n / 16 = 718, 1,500 long ranks (<= n / 4)
        string x = StringFamilies.Random(10_000, 4, 21);
        string text = x + x[..1500];
        AssertAttempt(text, succeeds: false);
        PreparedText p = TextPreparation.Prepare(text);
        Assert.That(Reference(p, SaisReference.BuildSuffixArray(text)).Kind, Is.EqualTo(LcpKind.U16));
        // the limit itself passes (16 * escapes <= n), one escape more does not; a random text of f symbols followed by a copy of its first c symbols has c - 254 escapes:
        // n = 4,096 = 16 * 256 with (f, c) = (3,586, 510) is exactly on the limit, (3,585, 511) one escape over it (steps stay below 2 n: about c^2 / 32)
        foreach (var (filler, copy, u8) in new[] { (3586, 510, true), (3585, 511, false) })
        {
            string t = StringFamilies.Random(filler, 4, 77);
            t += t[..copy];
            Assert.That(t.Length, Is.EqualTo(4096));
            int escapes = Reference(TextPreparation.Prepare(t), SaisReference.BuildSuffixArray(t)).ToArray().Count(v => v >= 255);
            Assert.That(escapes, Is.EqualTo(copy - 254).Or.EqualTo(copy - 253), $"copy {copy}");
            Assert.That(Reference(TextPreparation.Prepare(t), SaisReference.BuildSuffixArray(t)).Kind, Is.EqualTo(u8 ? LcpKind.U8 : LcpKind.U16), $"copy {copy}: escapes {escapes}");
            AssertAttempt(t, succeeds: u8);
        }
    }

    [Test]
    public void TheStepBudget_GivesUp_AndTheFallbackIsTheReference()
    {
        // 2,000-symbol block twice in 40,000: 1,746 escapes (16 * 1,746 <= 40,000: U8 holds), 1,984 long ranks (<= n / 4), about 2,000^2 / 32 = 125,000 steps > 2 n = 80,000
        string block = StringFamilies.Random(2000, 4, 31);
        string text = StringFamilies.Random(18_000, 4, 32) + block + StringFamilies.Random(18_000, 4, 33) + block;
        Assert.That(text.Length, Is.EqualTo(40_000));
        PreparedText p = TextPreparation.Prepare(text);
        Assert.That(Reference(p, SaisReference.BuildSuffixArray(text)).Kind, Is.EqualTo(LcpKind.U8));
        AssertAttempt(text, succeeds: false);
    }

    [Test]
    public void ArgumentChecks_AreTheOnesOfTheReference()
    {
        PreparedText p = TextPreparation.Prepare("aba");
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 3, 2, 0 }));        // wrong length
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 0, 1, 2, 3 }));     // sa[0] != n
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 2, 0, 1, 2 }));     // sa[0] != n, every other entry in range
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 3, 0, 1, 3 }));     // entry n at rank >= 1
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 3, 0, -1, 2 }));    // negative entry
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, p, new[] { 3, 0, 1, 99 }));    // entry far outside the text
        // a bad entry far from the start, after many valid ones (the scan, not a pre-check, finds it) and in a text longer than the vector width
        string text = StringFamilies.RandomDna(500, 4);
        PreparedText q = TextPreparation.Prepare(text);
        int[] sa = SaisReference.BuildSuffixArray(text);
        sa[400] = 500;
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, q, sa));
        sa[400] = -1;
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(LcpKernelKind.Direct, q, sa));
    }

    [Test]
    public void CorruptSuffixArrays_NeverLeaveTheBuffers_AndEitherFailCleanlyOrReturnAnArray()
    {
        // arbitrary permutations and arrays with repeated entries: the result is unspecified, the accesses are not (guarded arena here, bounds-checked accessors in the checked build)
        var rng = new SplitMix64(0x9753UL);
        for (int i = 0; i < 4000; i++)
        {
            int n = 1 + rng.NextInt(i < 2000 ? 40 : 400);
            string text = StringFamilies.Random(n, 1 + rng.NextInt(3), rng.Next());
            PreparedText p = TextPreparation.Prepare(text);
            int[] sa = SaisReference.BuildSuffixArray(text);
            if (i % 2 == 0)
                for (int k = n; k > 1; k--) { int r = 1 + rng.NextInt(k); (sa[k], sa[r]) = (sa[r], sa[k]); }
            else
                for (int k = 0, copies = 1 + rng.NextInt(4); k < copies; k++)
                {
                    int from = 1 + rng.NextInt(n), to = 1 + rng.NextInt(n);
                    sa[to] = sa[from];
                }
            using var arena = new NativeArena(ArenaOptions.Guarded((byte)(i & 0xFF)));
            try
            {
                LcpRun run = Engine.Run(p, sa, arena);
                Assert.That(run.Lcp.Length, Is.EqualTo(n + 1));
            }
            catch (InvalidOperationException)
            {
                // the fallback path rejects arrays whose escape counts do not match, like the reference
            }
            arena.VerifyCanaries();
        }
    }
}
