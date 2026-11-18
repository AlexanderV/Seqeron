using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The in-process RAM cap (ADR-004 section 4.4): the sum of live blocks of every arena that shares the budget may not exceed the cap.</summary>
[TestFixture]
public class ArenaBudgetTests
{
    [Test]
    public void Reserve_UpToTheCapSucceeds_AndOneByteMoreThrowsWithoutReserving()
    {
        var b = new ArenaBudget(100);
        b.Reserve(60, "a");
        b.Reserve(40, "b"); // exactly the cap
        Assert.That(b.LiveBytes, Is.EqualTo(100));
        ArenaBudgetExceededException ex = Assert.Throws<ArenaBudgetExceededException>(() => b.Reserve(1, "one-more"))!;
        Assert.That(ex.Message, Does.Contain("one-more").And.Contain("100 bytes exceeded").And.Contain("1 bytes; 100 live"), "the message names the block, the cap, the request and what is live");
        Assert.That(b.LiveBytes, Is.EqualTo(100), "a refused reservation reserves nothing");
        Assert.That(b.PeakBytes, Is.EqualTo(100), "and does not raise the peak");
        b.Release(100);
        Assert.That(b.LiveBytes, Is.Zero);
        Assert.That(b.PeakBytes, Is.EqualTo(100), "the peak stays");
    }

    [Test]
    public void Reserve_ZeroBytes_IsAlwaysAllowed_EvenWithCapZero()
    {
        var b = new ArenaBudget(0);
        b.Reserve(0, "empty");
        Assert.That(b.LiveBytes, Is.Zero);
        Assert.Throws<ArenaBudgetExceededException>(() => b.Reserve(1, "x"));
        Assert.That(b.CapBytes, Is.Zero);
    }

    [Test]
    public void Unlimited_IsLongMaxValue()
    {
        var b = new ArenaBudget(long.MaxValue);
        b.Reserve(long.MaxValue, "all");
        Assert.That(b.LiveBytes, Is.EqualTo(long.MaxValue));
        Assert.Throws<ArenaBudgetExceededException>(() => b.Reserve(1, "overflow"), "no wrap-around past long.MaxValue");
        Assert.That(b.LiveBytes, Is.EqualTo(long.MaxValue));
    }

    [Test]
    public void Arguments_AreValidated()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new ArenaBudget(-1));
        var b = new ArenaBudget(10);
        Assert.Throws<ArgumentOutOfRangeException>(() => b.Reserve(-1, "neg"));
        Assert.Throws<ArgumentOutOfRangeException>(() => b.Release(-1));
        b.Reserve(5, "a");
        Assert.Throws<InvalidOperationException>(() => b.Release(6), "giving back more than is live is a bug of the caller, not a bonus of budget");
        Assert.That(b.LiveBytes, Is.EqualTo(5), "the refused release changes nothing");
    }

    [Test]
    public void RandomSequence_FollowsTheModel_LiveAndPeakAndRefusals()
    {
        var rng = new SplitMix64(77);
        for (int round = 0; round < 50; round++)
        {
            long cap = rng.NextInt(1000);
            var b = new ArenaBudget(cap);
            long live = 0;
            long peak = 0;
            for (int step = 0; step < 200; step++)
            {
                if (live > 0 && rng.NextInt(3) == 0)
                {
                    long give = 1 + rng.NextInt((int)live);
                    b.Release(give);
                    live -= give;
                }
                else
                {
                    long ask = rng.NextInt(300);
                    if (live + ask > cap)
                    {
                        Assert.Throws<ArenaBudgetExceededException>(() => b.Reserve(ask, "x"), $"cap {cap} live {live} ask {ask}");
                    }
                    else
                    {
                        b.Reserve(ask, "x");
                        live += ask;
                        peak = Math.Max(peak, live);
                    }
                }
                Assert.That(b.LiveBytes, Is.EqualTo(live));
                Assert.That(b.PeakBytes, Is.EqualTo(peak));
            }
        }
    }

    [Test]
    public void Threads_NeverBreakTheCap_AndEndAtZero()
    {
        const int cap = 1000;
        var b = new ArenaBudget(cap);
        long violations = 0;
        Parallel.For(0, 8, new ParallelOptions { MaxDegreeOfParallelism = 8 }, t =>
        {
            for (int i = 0; i < 20000; i++)
            {
                try
                {
                    b.Reserve(100 + (t % 3), "t");
                }
                catch (ArenaBudgetExceededException)
                {
                    continue;
                }
                if (b.LiveBytes > cap) Interlocked.Increment(ref violations);
                b.Release(100 + (t % 3));
            }
        });
        Assert.That(violations, Is.Zero);
        Assert.That(b.LiveBytes, Is.Zero);
        Assert.That(b.PeakBytes, Is.LessThanOrEqualTo(cap));
        Assert.That(b.PeakBytes, Is.GreaterThan(0));
    }
}

/// <summary>An arena that charges its blocks to the budget first (so a refused block is never allocated) and gives them back on release and dispose.</summary>
[TestFixture]
public class BudgetedArenaTests
{
    [Test]
    public void Allocate_ChargesThePayload_AndDisposeGivesItBack()
    {
        var budget = new ArenaBudget(10_000);
        using (var arena = new BudgetedArena(budget))
        {
            Buf<uint> a = arena.Allocate<uint>(100, zero: true, "a");
            Buf<ulong> b = arena.Allocate<ulong>(50, zero: false, "b");
            Assert.That(arena.LiveBytes, Is.EqualTo(400 + 400));
            Assert.That(budget.LiveBytes, Is.EqualTo(800));
            Assert.That(a.AsSpan(100).ToArray(), Is.All.Zero);
            b.AsSpan(50).Fill(7);
            Assert.That(arena.Budget, Is.SameAs(budget));
        }
        Assert.That(budget.LiveBytes, Is.Zero);
        Assert.That(budget.PeakBytes, Is.EqualTo(800));
    }

    [Test]
    public void Allocate_OverTheCap_ThrowsAndChargesNothing_AndTheArenaStaysUsable()
    {
        var budget = new ArenaBudget(100);
        using var arena = new BudgetedArena(budget);
        arena.Allocate<byte>(60, zero: true, "first");
        Assert.Throws<ArenaBudgetExceededException>(() => arena.Allocate<byte>(41, zero: true, "second"));
        Assert.That(budget.LiveBytes, Is.EqualTo(60));
        Assert.That(arena.LiveBytes, Is.EqualTo(60));
        arena.Allocate<byte>(40, zero: true, "fits"); // exactly the cap
        Assert.That(budget.LiveBytes, Is.EqualTo(100));
    }

    [Test]
    public void FailedAllocation_InsideTheNativeArena_GivesTheChargeBack()
    {
        var budget = new ArenaBudget(1000);
        using var arena = new BudgetedArena(budget);
        Assert.Throws<ArgumentException>(() => arena.Allocate<byte>(10, zero: true, ""), "an empty name is refused by the native arena after the charge");
        Assert.That(budget.LiveBytes, Is.Zero, "the charge was taken back");
        Assert.Throws<ArgumentOutOfRangeException>(() => arena.Allocate<byte>(-1, zero: true, "neg"));
        Assert.That(budget.LiveBytes, Is.Zero);
    }

    [Test]
    public void ZeroLength_ChargesNothing_AndGivesAUsableBuffer()
    {
        var budget = new ArenaBudget(0);
        using var arena = new BudgetedArena(budget);
        Buf<uint> empty = arena.Allocate<uint>(0, zero: false, "empty");
        Assert.That(empty.AsSpan(0).Length, Is.Zero);
        Assert.That(budget.LiveBytes, Is.Zero);
    }

    [Test]
    public void MarkAndReleaseTo_GiveBackExactlyTheNewerBlocks()
    {
        var budget = new ArenaBudget(10_000);
        using var arena = new BudgetedArena(budget);
        arena.Allocate<byte>(10, zero: true, "keep");
        int mark = arena.Mark();
        arena.Allocate<uint>(25, zero: true, "x");
        arena.Allocate<ulong>(5, zero: true, "y");
        Assert.That(budget.LiveBytes, Is.EqualTo(10 + 100 + 40));
        arena.ReleaseTo(mark);
        Assert.That(budget.LiveBytes, Is.EqualTo(10));
        Assert.That(arena.LiveBytes, Is.EqualTo(10));
        arena.ReleaseTo(mark); // nothing newer: a no-op
        Assert.That(budget.LiveBytes, Is.EqualTo(10));
        arena.Allocate<byte>(5, zero: true, "again");
        Assert.That(budget.LiveBytes, Is.EqualTo(15));
        Assert.That(budget.PeakBytes, Is.EqualTo(150));
    }

    [Test]
    public void TwoArenas_ShareOneBudget_AndDisposeIsIdempotent()
    {
        var budget = new ArenaBudget(100);
        var a = new BudgetedArena(budget);
        var b = new BudgetedArena(budget);
        a.Allocate<byte>(60, zero: true, "a");
        Assert.Throws<ArenaBudgetExceededException>(() => b.Allocate<byte>(50, zero: true, "b"), "the sum of the arenas is capped");
        b.Allocate<byte>(40, zero: true, "b");
        a.Dispose();
        Assert.That(budget.LiveBytes, Is.EqualTo(40));
        a.Dispose(); // a second dispose must not give the bytes back twice
        Assert.That(budget.LiveBytes, Is.EqualTo(40));
        b.Dispose();
        Assert.That(budget.LiveBytes, Is.Zero);
        Assert.Throws<ObjectDisposedException>(() => a.Allocate<byte>(1, zero: true, "late"));
        Assert.Throws<ObjectDisposedException>(() => a.Mark());
    }

    [Test]
    public unsafe void GuardedArena_DetectsAStrayWrite_ThroughVerifyCanaries()
    {
        var budget = new ArenaBudget(10_000);
        using var arena = new BudgetedArena(budget, ArenaOptions.Guarded(0xCD));
        Buf<uint> a = arena.Allocate<uint>(8, zero: true, "a");
        arena.VerifyCanaries(); // untouched guard bands
        a.AsSpan(8).Fill(1);
        arena.VerifyCanaries();
        Assert.That(budget.LiveBytes, Is.EqualTo(32), "guard bands are not charged");
        uint* p = a.Pointer;
        p[8] = 0xDEADBEEF; // one past the end
        Assert.Throws<EsaMemoryViolationException>(arena.VerifyCanaries);
    }

    [Test]
    public void Constructor_RejectsANullBudget()
    {
        Assert.Throws<ArgumentNullException>(() => _ = new BudgetedArena(null!));
    }
}
