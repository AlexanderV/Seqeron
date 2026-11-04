using SuffixTree.Esa.Tests.Kernels;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Review;

/// <summary>
/// Tests added by the independent review of the optimized kernels (docs/refactoring/suffix-tree-upgrade/adr/REVIEW-S1a-opt.md). The first group kills fresh mutants that
/// survived every earlier test (mutant table in the review), the second group is the boundary sweep of item 5 (sizes around block / unroll / vector boundaries for every
/// alphabet class), the third pins the contract "the output never depends on the initial content of the managed output arrays".
/// </summary>
[TestFixture]
public class ReviewS1aOptTests
{
    // ---- Buf: the checked accessor must reject a read that is only ONE element/byte past the buffer (mutants R82-R85) -------------

    [Test]
    public void CheckedBuf_ReadUnaligned_RejectsEveryByteOutsideTheLogicalLength_AcceptsExactFit()
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<byte> whole = arena.Allocate<byte>(64, zero: true, "whole");
        Buf<byte> b = whole.Slice(8, 32); // logical length 32 inside a larger physical block: a Release read stays inside the memory
        Assert.That(b.ReadUnaligned<ulong>(24), Is.EqualTo(0UL), "exact fit: bytes 24..31");
        Assert.That(b.ReadUnaligned<ushort>(30), Is.EqualTo((ushort)0), "exact fit: bytes 30..31");
        if (EsaBuild.BoundsChecked)
        {
            for (int index = 25; index <= 33; index++)
            {
                int i = index;
                Assert.Throws<EsaMemoryViolationException>(() => _ = b.ReadUnaligned<ulong>(i), $"ulong at {index} ends past the buffer");
            }

            Assert.Throws<EsaMemoryViolationException>(() => _ = b.ReadUnaligned<ushort>(31), "ushort at 31 ends one byte past");
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.ReadUnaligned<ulong>(-1), "ulong at -1 starts before the buffer");
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.ReadUnaligned<byte>(-1));
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.ReadUnaligned<byte>(32));
        }

        Assert.DoesNotThrow(arena.VerifyCanaries);
    }

    [Test]
    public void CheckedBuf_LoadVector128_RejectsAVectorThatEndsOneElementPastTheBuffer_AcceptsExactFit()
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<byte> bytes = arena.Allocate<byte>(96, zero: true, "bytes").Slice(16, 32);
        Assert.DoesNotThrow(() => _ = bytes.LoadVector128(16), "bytes 16..31 of 32: exact fit");
        Assert.DoesNotThrow(() => _ = bytes.LoadVector128(0));
        Buf<ushort> shorts = arena.Allocate<ushort>(48, zero: true, "shorts").Slice(8, 16);
        Assert.DoesNotThrow(() => _ = shorts.LoadVector128(8), "8 ushorts at 8 of 16: exact fit");
        if (EsaBuild.BoundsChecked)
        {
            Assert.Throws<EsaMemoryViolationException>(() => _ = bytes.LoadVector128(17), "one element past");
            Assert.Throws<EsaMemoryViolationException>(() => _ = bytes.LoadVector128(-1), "one element before");
            Assert.Throws<EsaMemoryViolationException>(() => _ = shorts.LoadVector128(9));
            Assert.Throws<EsaMemoryViolationException>(() => _ = shorts.LoadVector128(-1));
        }

        Assert.DoesNotThrow(arena.VerifyCanaries);
    }

    [Test]
    public void CheckedBuf_Slice_AndAsSpan_RejectExactlyOneElementTooMany()
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<int> b = arena.Allocate<int>(10, zero: true, "ints").Slice(0, 10);
        Assert.DoesNotThrow(() => _ = b.Slice(4, 6));
        Assert.DoesNotThrow(() => _ = b.Slice(10, 0));
        Assert.DoesNotThrow(() => _ = b.AsSpan(10));
        if (EsaBuild.BoundsChecked)
        {
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.Slice(4, 7));
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.Slice(5, 6));
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.Slice(11, 0));
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.Slice(11));
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.Slice(-1, 1));
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.AsSpan(11));
            Assert.Throws<EsaMemoryViolationException>(() => _ = b[10]);
            Assert.Throws<EsaMemoryViolationException>(() => _ = b[-1]);
        }
    }

    [Test]
    public void GuardBands_DetectWrittenZeroAndFfByte_AtEveryOffsetOfBothBands()
    {
        // a kernel that stores 0 or -1 (dummy stores, cleared cells) through a wrong index must not hide behind the pattern (mutant R90: constant canary)
        int detected = 0;
        foreach (byte value in new byte[] { 0x00, 0xFF })
            foreach (int offset in new[] { -256, -255, -129, -64, -17, -16, -1, 24, 25, 100, 255 + 24 })
            {
                using var arena = new NativeArena(ArenaOptions.Guarded());
                Buf<byte> b = arena.Allocate<byte>(24, zero: true, "block");
                Assert.DoesNotThrow(arena.VerifyCanaries);
                byte pattern;
                unsafe
                {
                    pattern = b.Pointer[offset];
                    b.Pointer[offset] = value;
                }

                if (pattern == value)
                {
                    Assert.DoesNotThrow(arena.VerifyCanaries, "writing the pattern byte itself changes nothing");
                    continue;
                }

                detected++;
                var ex = Assert.Throws<EsaMemoryViolationException>(arena.VerifyCanaries, $"offset {offset}, value 0x{value:X2}");
                Assert.That(ex!.Kind, Is.EqualTo(MemoryViolationKind.CanaryOverwritten));
            }

        Assert.That(detected, Is.GreaterThanOrEqualTo(20), "almost every written 0x00 / 0xFF must change a pattern byte");

        // the pattern is not a constant: few of the 256 bytes of a band are 0x00 or 0xFF
        using var probe = new NativeArena(ArenaOptions.Guarded());
        Buf<byte> p = probe.Allocate<byte>(1, zero: true, "probe");
        int zeros = 0, ones = 0;
        unsafe
        {
            for (int i = -256; i < 0; i++)
            {
                if (p.Pointer[i] == 0x00) zeros++;
                if (p.Pointer[i] == 0xFF) ones++;
            }
        }

        Assert.That(zeros, Is.LessThan(8), "the canary pattern must not be (mostly) 0x00");
        Assert.That(ones, Is.LessThan(8), "the canary pattern must not be (mostly) 0xFF");
    }

    [Test]
    public void PoisonCoversTheWholeBlock_IncludingTheLastBytes()
    {
        // mutant R91: the poison fill stops one byte early, so a kernel that reads the last byte of an uninitialised block would see allocator garbage
        foreach (int length in new[] { 1, 2, 3, 7, 8, 9, 15, 16, 17, 100 })
        {
            using var arena = new NativeArena(ArenaOptions.Guarded(0xCD));
            Buf<byte> b = arena.Allocate<byte>(length, zero: false, "poison");
            Assert.That(b.AsSpan(length).ToArray(), Is.All.EqualTo((byte)0xCD), $"length {length}");
            Buf<int> ints = arena.Allocate<int>(length, zero: false, "poison-int");
            Assert.That(ints.AsSpan(length).ToArray(), Is.All.EqualTo(unchecked((int)0xCDCDCDCD)), $"length {length}");
        }
    }

    // ---- the output arrays are managed, uninitialised arrays: nothing may depend on their initial content --------------------------

    [Test]
    public void ProductionEntryPoints_DoNotDependOnTheInitialContentOfTheirOutputArrays()
    {
        // GC.AllocateUninitializedArray hands out stale memory for arrays above ~2 KB (large-object heap above 85,000 bytes: never cleared). The registry harness poisons
        // the ARENA blocks only; the SA array and the LCP arrays of the production entry points (SaisKernel.Build(PreparedText), LcpKernel.Build) are not in the arena, so a
        // kernel that forgot to write one cell (for example dst[0]) is invisible on a fresh heap and shows up only on a dirty one (mutant R60).
        foreach (int n in new[] { 3_000, 25_000, 90_000, 200_000 })
        {
            foreach (string family in new[] { "random-dna", "all-same", "period1000-noise", "dup400-p15", "random-binary" })
            {
                string text = ReviewS1aOptWorstCaseTests.Family(family, n);
                PreparedText prepared = TextPreparation.Prepare(text);
                int[] referenceSa = SaisReference.BuildSuffixArray(text);
                CompactLcp referenceLcp = LcpKernel.Build(LcpKernelKind.Reference, prepared, referenceSa);
                int[] expected = referenceLcp.ToArray();
                for (int round = 0; round < 3; round++)
                {
                    Dirty(n);
                    int[] sa = SaisKernel.Build(SaisKernelKind.Pointer, prepared);
                    Assert.That(sa, Is.EqualTo(referenceSa), $"SA {family} n={n} round {round}");
                    foreach (LcpKernelKind kind in new[] { LcpKernelKind.Pointer, LcpKernelKind.Direct })
                    {
                        Dirty(n);
                        CompactLcp lcp = LcpKernel.Build(kind, prepared, sa);
                        Assert.That(lcp.Kind, Is.EqualTo(referenceLcp.Kind), $"{kind} {family} n={n}");
                        Assert.That(lcp.ToArray(), Is.EqualTo(expected), $"{kind} {family} n={n} round {round}");
                    }
                }
            }
        }
    }

    /// <summary>Allocates and drops arrays of the sizes the kernels allocate uninitialised, filled with 0xFF, and forces a compacting collection.</summary>
    private static void Dirty(int n)
    {
        for (int k = 0; k < 4; k++)
        {
            byte[] bytes = new byte[n + 1];
            Array.Fill(bytes, (byte)0xFF);
            ushort[] shorts = new ushort[n + 1];
            Array.Fill(shorts, ushort.MaxValue);
            int[] ints = new int[n + 1];
            Array.Fill(ints, 0x7F7F7F7F);
            GC.KeepAlive(bytes);
            GC.KeepAlive(shorts);
            GC.KeepAlive(ints);
        }

        GC.Collect(2, GCCollectionMode.Forced, blocking: true, compacting: true);
    }

    // ---- no shared mutable state: concurrent builds of different texts -----------------------------------------------------------------

    [Test]
    public void ConcurrentBuilds_OfDifferentTexts_GiveTheReferenceResults()
    {
        // the kernels keep no static state (NativeArena, Buf and the probes are per call): four threads building different texts repeatedly must all get the reference result
        string[] texts =
        [
            ReviewS1aOptWorstCaseTests.Family("random-dna", 60_000),
            ReviewS1aOptWorstCaseTests.Family("fibonacci", 60_000),
            ReviewS1aOptWorstCaseTests.Family("u16-random-1000", 60_000),
            ReviewS1aOptWorstCaseTests.Family("dup400-p15", 60_000),
        ];
        var failures = new System.Collections.Concurrent.ConcurrentBag<string>();
        Parallel.For(0, texts.Length * 3, new ParallelOptions { MaxDegreeOfParallelism = 4 }, k =>
        {
            string text = texts[k % texts.Length];
            PreparedText prepared = TextPreparation.Prepare(text);
            int[] expectedSa = SaisReference.BuildSuffixArray(text);
            for (int round = 0; round < 3; round++)
            {
                int[] sa = SaisKernel.Build(SaisKernelKind.Pointer, prepared);
                if (!sa.AsSpan().SequenceEqual(expectedSa)) failures.Add($"SA differs (thread {k}, round {round})");
                CompactLcp reference = LcpKernel.Build(LcpKernelKind.Reference, prepared, sa);
                foreach (LcpKernelKind kind in new[] { LcpKernelKind.Pointer, LcpKernelKind.Direct })
                {
                    CompactLcp lcp = LcpKernel.Build(kind, prepared, sa);
                    if (lcp.Kind != reference.Kind || !lcp.ToArray().AsSpan().SequenceEqual(reference.ToArray())) failures.Add($"LCP {kind} differs (thread {k}, round {round})");
                }
            }
        });
        Assert.That(failures, Is.Empty);
    }

    // ---- boundary sweep (item 5): sizes around block, unroll and vector boundaries, every alphabet class -------------------------

    private static readonly int[] SmallSizes =
        [0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 47, 48, 49, 63, 64, 65, 127, 128, 129, 255, 256, 257, 511, 512, 513, 1023, 1024, 1025];

    private static readonly int[] LargeSizes = [4095, 4096, 4097, 65_535, 65_536, 65_537, 131_071, 131_072, 131_073];

    private static readonly int[] Sigmas = [1, 2, 3, 4, 5, 16, 20, 126, 127, 128, 129, 200, 255, 256, 257, 1000, 4096, 65_535, 65_536];

    [Test]
    public void BoundarySweep_SmallSizes_EveryAlphabetClass_AllOptimizedEngines() => Sweep(SmallSizes, allPoisons: true);

    [Test]
    [Category("ReviewSweep")]
    [Explicit("Large sizes (up to 131,073 symbols, sigma up to 65,536): ~1 minute; run with --filter Category=ReviewSweep")]
    public void BoundarySweep_LargeSizes_EveryAlphabetClass_AllOptimizedEngines() => Sweep(LargeSizes, allPoisons: false);

    private static void Sweep(int[] sizes, bool allPoisons)
    {
        int cases = 0;
        foreach (int n in sizes)
            foreach (int sigma in Sigmas)
            {
                if (sigma > n && !(n == 0 && sigma == 1)) continue; // the realised alphabet is exactly sigma: every symbol occurs at least once
                foreach (string text in Texts(n, sigma))
                {
                    foreach (SaisEngine e in KernelRegistry.OptimizedSais) KernelHarness.AssertSais(e, text, allPoisons);
                    foreach (LcpEngine e in KernelRegistry.OptimizedLcp) KernelHarness.AssertLcp(e, text, allPoisons);
                    cases++;
                }
            }

        Assert.That(cases, Is.GreaterThan(sizes.Length));
    }

    /// <summary>Texts of length n whose alphabet is exactly sigma symbols (when n &gt;= sigma): shuffled alphabet then random, then repeat-heavy tails.</summary>
    private static IEnumerable<string> Texts(int n, int sigma)
    {
        if (n == 0) { yield return string.Empty; yield break; }
        var rng = new SplitMix64(((ulong)n << 20) ^ (ulong)sigma ^ 0xB0A7UL);
        char[] alphabet = new char[sigma];
        for (int k = 0; k < sigma; k++) alphabet[k] = (char)(sigma <= 256 ? 0x20 + k : 0x100 + k);
        char[] shuffled = [.. alphabet];
        for (int i = sigma - 1; i > 0; i--)
        {
            int j = rng.NextInt(i + 1);
            (shuffled[i], shuffled[j]) = (shuffled[j], shuffled[i]);
        }

        // 1: shuffled alphabet, then uniformly random symbols
        var a = new char[n];
        for (int i = 0; i < n; i++) a[i] = i < sigma ? shuffled[i] : alphabet[rng.NextInt(sigma)];
        yield return new string(a);

        // 2: shuffled alphabet, then a run of the smallest symbol, the last symbol different (a^k b tail: the longest possible LCPs at the end)
        var b = new char[n];
        for (int i = 0; i < n; i++) b[i] = i < sigma ? shuffled[i] : alphabet[0];
        if (n > sigma) b[n - 1] = alphabet[sigma - 1];
        yield return new string(b);

        // 3: shuffled alphabet, then a period-3 tail over random symbols
        var c = new char[n];
        char[] unit = [alphabet[rng.NextInt(sigma)], alphabet[rng.NextInt(sigma)], alphabet[rng.NextInt(sigma)]];
        for (int i = 0; i < n; i++) c[i] = i < sigma ? shuffled[i] : unit[(i - sigma) % 3];
        yield return new string(c);

        // 4: the alphabet in ascending order, repeated (period sigma): n / sigma equal blocks
        var d = new char[n];
        for (int i = 0; i < n; i++) d[i] = alphabet[i % sigma];
        yield return new string(d);

        // 5: descending then the single symbol run (the text is already in sorted-reverse order)
        var e = new char[n];
        for (int i = 0; i < n; i++) e[i] = i < sigma ? alphabet[sigma - 1 - i] : alphabet[sigma - 1];
        yield return new string(e);
    }
}
