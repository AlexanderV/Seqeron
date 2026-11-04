using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>
/// 1o.1a: the LCP pointer kernel (<see cref="LcpKernel"/>) beyond the generic registry run: the escape boundaries, word-compare boundaries, tiny texts,
/// the public entry points, the kernel selector, and memory safety on corrupt input (garbage in, garbage out, never an access outside the buffers).
/// </summary>
[TestFixture]
public class LcpKernelTests
{
    private static LcpEngine Engine => KernelRegistry.PointerLcp;

    [Test]
    public void TinyTexts_n0_n1_n2_n3_GoThroughTheEmptyAndNullPointerPaths()
    {
        foreach (string text in new[] { "", "a", "aa", "ab", "ba", "aaa", "aba", "abc", "cab" })
        {
            KernelHarness.AssertLcp(Engine, text);
            PreparedText p = TextPreparation.Prepare(text);
            int[] sa = SaisReference.BuildSuffixArray(text);
            Assert.That(LcpKernel.Build(p, sa).ToArray(), Is.EqualTo(NaiveLcp.Build(text, sa)), $"public entry on '{text}'");
        }
    }

    [Test]
    public void EscapeBoundaries_U8And_U16_AreHitExactly()
    {
        // two copies of a random block R: LCP values reach |R|; 254/255/256 and 65,534/65,535/65,536 straddle the two escape values
        foreach (int r in new[] { 254, 255, 256, 257, 65_534, 65_535, 65_536 })
        {
            string block = StringFamilies.Random(r, 4, (ulong)r);
            KernelHarness.AssertLcp(Engine, block + block, allPoisons: false);
        }
        KernelHarness.AssertLcp(Engine, StringFamilies.AllSame(70_000), allPoisons: false); // U16 with 4,465 escapes
        KernelHarness.AssertLcp(Engine, StringFamilies.Fibonacci(200_000), allPoisons: false); // U32
    }

    [Test]
    public void CommonPrefixesThatStraddleTheEightByteWord_AndTheTextEnd()
    {
        // periodic texts: the LCP of neighbouring suffixes is a multiple of the period near n, so word boundaries (8 bytes / 4 ushorts) and the end of the text are hit at every residue
        for (int period = 1; period <= 20; period++)
            for (int n = 0; n <= 70; n++)
            {
                KernelHarness.AssertLcp(Engine, StringFamilies.PeriodicAaB(n, period + 1), allPoisons: false);
                KernelHarness.AssertLcp(Engine, new string([.. StringFamilies.PeriodicAaB(n, period + 1).Select(c => (char)(c + 0x200))]), allPoisons: false); // ushort symbols
            }
    }

    [Test]
    public void SixteenBitSymbols_WithLongCommonPrefixes_PlaceTheMismatchInsideTheWord()
    {
        // 300 distinct leading symbols force the 16-bit reader (an alphabet of at most 256 symbols is remapped to bytes even for units above 255); the periodic tail
        // gives long common prefixes whose first mismatch falls at every position of a four-symbol word (mutant K13, mismatch offset `bits >> 3` instead of `>> 4`)
        string distinct = new([.. Enumerable.Range(0, 300).Select(i => (char)(0x100 + i))]);
        Assert.That(TextPreparation.Prepare(distinct + "ab").Width, Is.EqualTo(SymbolWidth.UShort));
        for (int period = 1; period <= 20; period++)
            for (int n = 0; n <= 70; n++)
                KernelHarness.AssertLcp(Engine, distinct + StringFamilies.PeriodicAaB(n, period + 1), allPoisons: false);
    }

    [Test]
    public void TheSelectorAndTheTwoEntryPoints_GiveTheReferenceResult()
    {
        string text = StringFamilies.RandomDna(5000, 3);
        PreparedText p = TextPreparation.Prepare(text);
        int[] sa = SaisReference.BuildSuffixArray(text);
        CompactLcp reference = LcpKernel.Build(LcpKernelKind.Reference, p, sa);
        CompactLcp pointer = LcpKernel.Build(LcpKernelKind.Pointer, p, sa);
        Assert.That(pointer.Kind, Is.EqualTo(reference.Kind));
        Assert.That(pointer.ToArray(), Is.EqualTo(reference.ToArray()));
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Assert.That(LcpKernel.Build(p, sa, arena).ToArray(), Is.EqualTo(reference.ToArray()));
        arena.VerifyCanaries();
        var trace = new LcpTrace();
        Assert.That(LcpKernel.Build(p, sa, trace).ToArray(), Is.EqualTo(reference.ToArray()));
        Assert.That(trace.PhiTicks + trace.PlcpTicks + trace.CountTicks + trace.OutputTicks, Is.GreaterThan(0));
    }

    [Test]
    public void ArgumentChecks_AreTheOnesOfTheReference()
    {
        PreparedText p = TextPreparation.Prepare("aba");
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(p, new[] { 3, 2, 0 }));        // wrong length
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(p, new[] { 0, 1, 2, 3 }));     // sa[0] != n
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(p, new[] { 2, 0, 1, 2 }));     // sa[0] != n, every other entry in range (only this check can reject it)
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(p, new[] { 3, 0, 1, 2, 2 }));  // wrong length (too long)
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(p, new[] { 3, 0, 1, 3 }));     // entry n at rank >= 1
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(p, new[] { 3, 0, -1, 2 }));    // negative entry
        Assert.Throws<ArgumentException>(() => LcpKernel.Build(p, new[] { 3, 0, 1, 99 }));    // entry far outside the text
    }

    [Test]
    public void AnArbitraryPermutationInsteadOfASuffixArray_GivesTheReferencePlcp_WithoutLeavingTheBuffers()
    {
        // review test for the reference (M25): garbage in, garbage out; the kernel has no bounds check in Release, so it is run in a guarded arena (and in the checked build)
        var rng = new SplitMix64(0x1357UL);
        for (int i = 0; i < 4000; i++)
        {
            int n = 1 + rng.NextInt(40);
            string text = StringFamilies.Random(n, 1 + rng.NextInt(3), rng.Next());
            PreparedText p = TextPreparation.Prepare(text);
            var sa = new int[n + 1];
            sa[0] = n;
            for (int k = 1; k <= n; k++) sa[k] = k - 1;
            for (int k = n; k > 1; k--)
            {
                int r = 1 + rng.NextInt(k);
                (sa[k], sa[r]) = (sa[r], sa[k]);
            }
            using var arena = new NativeArena(ArenaOptions.Guarded((byte)(i & 0xFF)));
            LcpRun run = Engine.Run(p, sa, arena);
            arena.VerifyCanaries();
            int[] expected = p.Width == SymbolWidth.Byte ? LcpReference.ComputePlcp(p.ByteReader(), sa) : LcpReference.ComputePlcp(p.UShortReader(), sa);
            Assert.That(run.Plcp, Is.EqualTo(expected), $"'{text}' sa [{string.Join(',', sa)}]");
            for (int j = 0; j < n; j++) Assert.That(run.Plcp![j], Is.InRange(0, n - j));
        }
    }

    [Test]
    public void ANonPermutationWithRepeatedEntries_NeverWritesOutsideTheEscapeLists()
    {
        // with a repeated entry the escape count announced from PLCP differs from the number of escapes pass 3 meets; the reference fails with
        // InvalidOperationException ("more escapes than announced" / "escapes written"), the kernel must fail the same way and never write past the lists
        var rng = new SplitMix64(0x2468UL);
        int differing = 0;
        for (int i = 0; i < 3000; i++)
        {
            // a random block of 256..285 symbols twice: 2..31 PLCP values reach 255 and the rule still picks U8, so a wrong copy count changes the number of escapes met in pass 3
            string block = StringFamilies.Random(256 + rng.NextInt(30), 4, rng.Next());
            string text = block + block;
            int n = text.Length;
            PreparedText p = TextPreparation.Prepare(text);
            int[] sa = SaisReference.BuildSuffixArray(text);
            for (int k = 0, copies = 1 + rng.NextInt(4); k < copies; k++)
            {
                int from = 1 + rng.NextInt(n), to = 1 + rng.NextInt(n);
                sa[to] = sa[from]; // a duplicate entry; all entries stay in [0, n)
            }
            using var arena = new NativeArena(ArenaOptions.Guarded());
            Exception? expected = null;
            CompactLcp? reference = null;
            try
            {
                reference = p.Width == SymbolWidth.Byte ? LcpReference.Build(p.ByteReader(), sa) : LcpReference.Build(p.UShortReader(), sa);
            }
            catch (InvalidOperationException ex)
            {
                expected = ex;
            }
            try
            {
                CompactLcp kernel = LcpKernel.Build(p, sa, arena);
                Assert.That(expected, Is.Null, "the reference rejected the array, the kernel did not");
                Assert.That(kernel.ToArray(), Is.EqualTo(reference!.ToArray()));
            }
            catch (InvalidOperationException) when (expected is not null)
            {
                differing++;
            }
            arena.VerifyCanaries();
        }
        TestContext.Out.WriteLine($"{differing} of 3000 corrupt arrays rejected by both");
    }
}
