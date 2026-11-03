namespace SuffixTree.Esa.Tests.Safety;

/// <summary>
/// The safety infrastructure itself (work package 1o.0, review MAJ-1): the accessor <see cref="Buf{T}"/> (a bare pointer in Release, bounds-checked with
/// <c>ESA_BOUNDS</c>) and the guarded <see cref="NativeArena"/> (canary bands, poison, zeroing). The same tests run in both build modes
/// (<c>dotnet test ... -p:EsaBounds=true</c> for the checked one); the expected outcome of an out-of-bounds write is stated per mode:
/// the checked build throws at the access, the Release build lets the write land in the guard band, where <see cref="NativeArena.VerifyCanaries"/> finds it.
/// </summary>
[TestFixture]
public class BufAndArenaTests
{
    private const int Length = 10;
    private const int GuardInts = 64; // ArenaOptions.Guarded: 256-byte bands

    [Test]
    public void BuildMode_IsReported_AndMatchesTheExpectedMode()
    {
        TestContext.Out.WriteLine($"ESA_BOUNDS (checked accessor): {EsaBuild.BoundsChecked}");
        // the checked run sets ESA_EXPECT_BOUNDS=true so that a run that silently fell back to the Release accessor cannot pass as "checked"
        string? expected = Environment.GetEnvironmentVariable("ESA_EXPECT_BOUNDS");
        if (expected is null) Assert.Pass("no expectation set");
        Assert.That(EsaBuild.BoundsChecked, Is.EqualTo(bool.Parse(expected!)));
    }

    [Test]
    public void InBoundsUse_ReadsWritesSlicesAndSpans_AndLeavesTheCanariesIntact()
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<int> b = arena.Allocate<int>(Length, zero: true, "b");
        for (int i = 0; i < Length; i++) b[i] = i * i;
        b[3]++;
        Buf<int> tail = b.Slice(6);
        Buf<int> mid = b.Slice(2, 4);
        Assert.That(b[3], Is.EqualTo(10));
        Assert.That(tail[0], Is.EqualTo(36));
        Assert.That(mid[1], Is.EqualTo(10));
        Assert.That(b.AsSpan(Length).ToArray(), Is.EqualTo(new[] { 0, 1, 4, 10, 16, 25, 36, 49, 64, 81 }));
        b.Fill(7, Length);
        Assert.That(b.AsSpan(Length).ToArray(), Is.All.EqualTo(7));
        b.Clear(Length);
        Assert.That(b.AsSpan(Length).ToArray(), Is.All.EqualTo(0));
        Buf<byte> bytes = arena.Allocate<byte>(16, zero: true, "bytes");
        bytes[1] = 0x34;
        bytes[2] = 0x12;
        Assert.That(bytes.ReadUnaligned<ushort>(1), Is.EqualTo((ushort)0x1234));
        Assert.That(bytes.ReadUnaligned<ulong>(8), Is.EqualTo(0UL)); // the last whole word
        Assert.DoesNotThrow(arena.VerifyCanaries);
    }

    [TestCase(-GuardInts)]
    [TestCase(-1)]
    [TestCase(Length)]
    [TestCase(Length + GuardInts - 1)]
    public void WriteOutsideTheLogicalLength_IsDetected_ByTheCheckedAccessorOrByTheCanary(int index)
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<int> b = arena.Allocate<int>(Length, zero: true, "victim");
        AssertDetected(arena, "victim", () => b[index] = 0x13579BDF, index < 0 ? "before" : "after");
    }

    [Test]
    public void SingleByteWrites_AtTheFirstAndLastGuardByte_AreDetected()
    {
        foreach (int index in new[] { -256, -1, 7, 7 + 255 })
        {
            using var arena = new NativeArena(ArenaOptions.Guarded());
            Buf<byte> b = arena.Allocate<byte>(7, zero: false, "bytes7");
            AssertDetected(arena, "bytes7", () => b[index] = 0, index < 0 ? "before" : "after");
        }
    }

    [Test]
    public unsafe void AWriteThatBypassesTheAccessor_IsStillFoundByTheCanary_InEveryBuildMode()
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<int> b = arena.Allocate<int>(Length, zero: true, "bypassed");
        b.Pointer[Length + 3] = 1; // raw pointer arithmetic: the checked accessor never sees it
        var ex = Assert.Throws<EsaMemoryViolationException>(arena.VerifyCanaries)!;
        Assert.That(ex.Kind, Is.EqualTo(MemoryViolationKind.CanaryOverwritten));
        Assert.That(ex.Message, Does.Contain("bypassed").And.Contain("after"));
    }

    [Test]
    public void TheViolationMessage_NamesTheBlock_TheSide_AndTheByteOffset()
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        _ = arena.Allocate<int>(4, zero: true, "first");
        Buf<int> second = arena.Allocate<int>(Length, zero: true, "second");
        string message = Capture(arena, () => second[-2] = 5);
        Assert.That(message, Does.Contain("second"));
        // checked build: "index -2, count 1, length 10"; Release: canary before the payload at byte offset -8 (two ints before it)
        Assert.That(message, EsaBuild.BoundsChecked ? Does.Contain("index -2") : Does.Contain("before").And.Contain("byte offset -8"));
    }

    [Test]
    public void ReadsPastTheEnd_AreOnlyDetectedByTheCheckedBuild()
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<int> b = arena.Allocate<int>(Length, zero: true, "reader");
        if (EsaBuild.BoundsChecked)
        {
            Assert.Throws<EsaMemoryViolationException>(() => _ = b[Length]);
            Assert.Throws<EsaMemoryViolationException>(() => _ = b[-1]);
            Assert.Throws<EsaMemoryViolationException>(() => _ = b.ReadUnaligned<ulong>(Length - 1));
            Assert.Throws<EsaMemoryViolationException>(() => b.Slice(5, 6));
            Assert.Throws<EsaMemoryViolationException>(() => b.Slice(Length + 1));
            Assert.Throws<EsaMemoryViolationException>(() => b.AsSpan(Length + 1));
            Assert.DoesNotThrow(() => _ = b.Slice(Length));
            Assert.DoesNotThrow(() => _ = b.ReadUnaligned<ulong>(Length - 2));
        }
        else
        {
            // a read lands in the guard band: memory safe, invisible to the canary (it checks writes), which is why the checked build exists
            Assert.DoesNotThrow(() => _ = b[Length]);
            Assert.DoesNotThrow(() => _ = b.ReadUnaligned<ulong>(Length - 1));
        }
        Assert.DoesNotThrow(arena.VerifyCanaries);
    }

    [Test]
    public void UninitialisedBlocks_CarryThePoison_ZeroedBlocksAreZero_InGuardedAndProductionArenas()
    {
        using var guarded = new NativeArena(ArenaOptions.Guarded(0xA7));
        Buf<int> p = guarded.Allocate<int>(Length, zero: false, "poisoned");
        Buf<int> z = guarded.Allocate<int>(Length, zero: true, "zeroed");
        Assert.That(p.AsSpan(Length).ToArray(), Is.All.EqualTo(unchecked((int)0xA7A7A7A7)));
        Assert.That(z.AsSpan(Length).ToArray(), Is.All.EqualTo(0));

        using var production = new NativeArena();
        Assert.That(production.Guarded, Is.False);
        Buf<long> big = production.Allocate<long>(1 << 20, zero: true, "big");
        Assert.That(big.AsSpan(1 << 20).ToArray().All(v => v == 0), Is.True);
        Assert.DoesNotThrow(production.VerifyCanaries); // no guard bands: nothing to verify
    }

    [Test]
    public void EmptyBlocks_AreValid_AndDisposeIsIdempotent_AndAFreedArenaRefusesNewBlocks()
    {
        var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<int> empty = arena.Allocate<int>(0, zero: true, "empty");
        Assert.That(empty.AsSpan(0).Length, Is.EqualTo(0));
        Assert.DoesNotThrow(() => empty.Fill(1, 0));
        Assert.DoesNotThrow(arena.VerifyCanaries);
        Assert.That(arena.BlockCount, Is.EqualTo(1));
        arena.Allocate<byte>(5, zero: false, "five");
        Assert.That(arena.AllocatedBytes, Is.EqualTo(5));
        arena.Dispose();
        Assert.DoesNotThrow(arena.Dispose);
        Assert.Throws<ObjectDisposedException>(() => arena.Allocate<int>(1, zero: true, "late"));
    }

    [Test]
    public unsafe void MarkAndReleaseTo_FreeTheNewestBlocks_VerifyTheirCanaries_AndKeepTheOlderOnes()
    {
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<int> keep = arena.Allocate<int>(4, zero: true, "keep");
        int mark = arena.Mark();
        _ = arena.Allocate<int>(Length, zero: true, "a");
        _ = arena.Allocate<byte>(33, zero: false, "b");
        Assert.That(arena.BlockCount, Is.EqualTo(3));
        arena.ReleaseTo(mark);
        Assert.That(arena.BlockCount, Is.EqualTo(1));
        keep[3] = 7;
        Assert.DoesNotThrow(arena.VerifyCanaries);

        // a damaged block is freed anyway, then the violation is reported (the damage is not lost by the release)
        int mark2 = arena.Mark();
        Buf<int> damaged = arena.Allocate<int>(Length, zero: true, "damaged");
        _ = arena.Allocate<int>(2, zero: true, "innocent");
        damaged.Pointer[Length + 1] = 1; // bypasses the accessor
        var ex = Assert.Throws<EsaMemoryViolationException>(() => arena.ReleaseTo(mark2))!;
        Assert.That(ex.Kind, Is.EqualTo(MemoryViolationKind.CanaryOverwritten));
        Assert.That(ex.Message, Does.Contain("damaged"));
        Assert.That(arena.BlockCount, Is.EqualTo(1));

        Assert.Throws<ArgumentOutOfRangeException>(() => arena.ReleaseTo(2));
        Assert.Throws<ArgumentOutOfRangeException>(() => arena.ReleaseTo(-1));
        Assert.DoesNotThrow(() => arena.ReleaseTo(arena.Mark())); // nothing newer than the mark
    }

    [Test]
    public void ArgumentChecks()
    {
        using var arena = new NativeArena();
        Assert.Throws<ArgumentOutOfRangeException>(() => arena.Allocate<int>(-1, zero: true, "x"));
        Assert.Throws<ArgumentException>(() => arena.Allocate<int>(1, zero: true, ""));
        Assert.Throws<ArgumentOutOfRangeException>(() => new NativeArena(new ArenaOptions(-1, null)));
    }

    private static void AssertDetected(NativeArena arena, string block, Action access, string side)
    {
        string message = Capture(arena, access);
        Assert.That(message, Does.Contain(block));
        if (!EsaBuild.BoundsChecked) Assert.That(message, Does.Contain(side));
    }

    [Test]
    public void Prefetch_NeverThrowsAndNeverTouchesMemory_ForAnyIndex_InBothBuildModes()
    {
        // 1o.2: a prefetch address is not checked (the hardware ignores an invalid one); the kernels prefetch positions that come from the data, also clamped ones
        using var arena = new NativeArena(ArenaOptions.Guarded());
        Buf<int> b = arena.Allocate<int>(Length, zero: true, "b");
        for (int i = 0; i < Length; i++) b[i] = i + 1;
        foreach (int index in new[] { 0, Length - 1, Length, Length + 1000, -1, -100_000, int.MaxValue / 8, int.MinValue / 8 })
            Assert.DoesNotThrow(() => b.Prefetch(index), $"index {index}");
        Assert.That(b.AsSpan(Length).ToArray(), Is.EqualTo(Enumerable.Range(1, Length)));
        arena.VerifyCanaries();
    }

    /// <summary>Runs the access; returns the message of the violation that the build mode reports for it (accessor exception, else canary), failing when none is reported.</summary>
    private static string Capture(NativeArena arena, Action access)
    {
        try
        {
            access();
        }
        catch (EsaMemoryViolationException ex)
        {
            Assert.That(EsaBuild.BoundsChecked, Is.True, "the Release accessor must not throw");
            Assert.That(ex.Kind, Is.EqualTo(MemoryViolationKind.OutOfBounds));
            return ex.Message;
        }
        Assert.That(EsaBuild.BoundsChecked, Is.False, "the checked accessor must reject the access");
        var canary = Assert.Throws<EsaMemoryViolationException>(arena.VerifyCanaries)!;
        Assert.That(canary.Kind, Is.EqualTo(MemoryViolationKind.CanaryOverwritten));
        return canary.Message;
    }
}
