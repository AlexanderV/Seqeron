using System.Runtime.InteropServices;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph5;
using SuffixTree.Esa;

namespace SuffixTree.EmProto.Tests.Ph5;

/// <summary>An entry as the interval scan keeps it: (l, lb, first, last), 16 bytes.</summary>
[StructLayout(LayoutKind.Sequential)]
internal struct TestEntry
{
    internal uint A;
    internal uint B;
    internal uint C;
    internal uint D;

    internal static TestEntry Of(uint x) => new() { A = x, B = x * 3 + 1, C = x ^ 0x5A5A5A5Au, D = ~x };
}

/// <summary>
/// The external stack of PH5 (ADR-004 section 3.2: the top two blocks in RAM, older blocks appended to a spill file, read back backwards when a pop leaves the window, at least one full block hot): a model-based test
/// against <see cref="Stack{T}"/> for every block size from 2 to 50, the spill arithmetic, the hysteresis at a block boundary, the dirty-block rule and the clean-up.
/// </summary>
[TestFixture]
public class SpillStackTests
{
    private string _dir = "";

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "spill-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private static BudgetedArena Arena(ArenaBudget? budget = null) => new(budget ?? new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));

    [Test]
    public void ModelFuzz_EveryBlockSize_TopAndHeightFollowTheModel()
    {
        var rng = new Random(12345);
        for (int block = 2; block <= 50; block++)
        {
            using BudgetedArena arena = Arena();
            using var stack = new SpillStack<TestEntry>(arena, Path.Combine(_dir, $"s{block}.spill"), block);
            var model = new Stack<TestEntry>();
            uint next = 1;
            long peak = 0;
            int phase = 0;
            for (int step = 0; step < 4000; step++)
            {
                // phases of mostly pushing and mostly popping so that the stack climbs through many blocks and comes back
                if (step % 400 == 0) phase = rng.Next(3);
                int pushShare = phase == 0 ? 70 : phase == 1 ? 30 : 50;
                bool push = model.Count == 0 || rng.Next(100) < pushShare;
                if (push)
                {
                    TestEntry e = TestEntry.Of(next++);
                    stack.Push(e);
                    model.Push(e);
                }
                else if (rng.Next(4) == 0)
                {
                    // update the top entry in place: the new value must survive a spill and a refill
                    ref TestEntry top = ref stack.TopForUpdate();
                    top.D = next++;
                    TestEntry m = model.Pop();
                    m.D = top.D;
                    model.Push(m);
                }
                else
                {
                    stack.Pop();
                    model.Pop();
                }
                peak = Math.Max(peak, model.Count);
                Assert.That(stack.Height, Is.EqualTo(model.Count), $"block {block}, step {step}: height");
                if (model.Count > 0)
                {
                    TestEntry t = model.Peek();
                    TestEntry g = stack.Top;
                    Assert.That((g.A, g.B, g.C, g.D), Is.EqualTo((t.A, t.B, t.C, t.D)), $"block {block}, step {step}: top");
                }
            }
            Assert.That(stack.PeakHeight, Is.EqualTo(peak), $"block {block}: peak height");
            while (model.Count > 0)
            {
                TestEntry t = model.Pop();
                Assert.That(stack.Top.A, Is.EqualTo(t.A), $"block {block}: drain");
                stack.Pop();
            }
            Assert.That(stack.Height, Is.Zero);
        }
    }

    [TestCase(2, 1)]
    [TestCase(2, 2)]
    [TestCase(2, 7)]
    [TestCase(3, 10)]
    [TestCase(4, 16)]
    [TestCase(10, 95)]
    [TestCase(10, 101)]
    [TestCase(50, 1000)]
    public void PushThenPopAll_SpillsThePredictedNumberOfBlocksAndReadsThemBack(int block, int height)
    {
        using BudgetedArena arena = Arena();
        using var stack = new SpillStack<TestEntry>(arena, Path.Combine(_dir, "s.spill"), block);
        for (uint i = 0; i < height; i++) stack.Push(TestEntry.Of(i));
        // blocks 0 .. t - 2 are spilled when the top block is t = (height - 1) / block
        long expected = Math.Max(0, ((height - 1) / block) - 1);
        Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(expected), "blocks written while climbing");
        Assert.That(stack.Stats.BytesWritten, Is.EqualTo(expected * block * 16), "bytes written");
        Assert.That(stack.Stats.PeakFileBytes, Is.EqualTo(expected * block * 16), "peak file size");
        for (int i = height - 1; i >= 0; i--)
        {
            Assert.That(stack.Top.A, Is.EqualTo((uint)i), $"entry {i} after the refills");
            stack.Pop();
        }
        Assert.That(stack.Stats.BlocksRead, Is.EqualTo(expected), "every spilled block is read back once");
        Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(expected), "nothing is written twice when nothing is modified");
    }

    [Test]
    public void Oscillation_AtEveryBlockBoundary_DoesNotThrash()
    {
        const int block = 8;
        foreach (int boundary in new[] { block, 2 * block, 3 * block, 4 * block })
        {
            using BudgetedArena arena = Arena();
            using var stack = new SpillStack<TestEntry>(arena, Path.Combine(_dir, $"o{boundary}.spill"), block);
            for (uint i = 0; i < boundary; i++) stack.Push(TestEntry.Of(i));
            SpillStats before = stack.Stats;
            for (int k = 0; k < 2000; k++)
            {
                stack.Push(TestEntry.Of(1000 + (uint)k));
                stack.Pop();
                stack.Pop();
                stack.Push(TestEntry.Of(uint.MaxValue - (uint)k));
            }
            SpillStats after = stack.Stats;
            Assert.That(after.BlocksWritten - before.BlocksWritten, Is.LessThanOrEqualTo(1), $"boundary {boundary}: writes");
            Assert.That(after.BlocksRead - before.BlocksRead, Is.LessThanOrEqualTo(1), $"boundary {boundary}: reads");
        }
    }

    [Test]
    public void CleanBlock_AfterARefill_IsNotWrittenAgain_AndAModifiedOneIs()
    {
        const int block = 4;
        foreach (bool modify in new[] { false, true })
        {
            using BudgetedArena arena = Arena();
            using var stack = new SpillStack<TestEntry>(arena, Path.Combine(_dir, $"c{modify}.spill"), block);
            for (uint i = 0; i < 6 * block; i++) stack.Push(TestEntry.Of(i)); // top block 5, blocks 0..3 spilled
            Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(4));
            for (int i = 0; i < 2 * block; i++) stack.Pop(); // block 3 is read back and becomes the top block: its file copy is current
            Assert.That(stack.Stats.BlocksRead, Is.EqualTo(1));
            Assert.That(stack.Top.A, Is.EqualTo(15u));
            if (modify) stack.TopForUpdate().D = 777;
            for (uint i = 0; i < 5; i++) stack.Push(TestEntry.Of(100 + i)); // the fifth push moves on from block 4 and evicts block 3
            if (modify)
            {
                Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(5), "a block that was modified after its refill is written");
                Assert.That(stack.Stats.CleanEvictions, Is.Zero);
            }
            else
            {
                Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(4), "an unmodified refilled block is evicted without a write");
                Assert.That(stack.Stats.CleanEvictions, Is.EqualTo(1));
            }
            for (int i = 0; i < 5; i++) stack.Pop();
            Assert.That(stack.Top.A, Is.EqualTo(15u), "block 3 is read back again");
            Assert.That(stack.Top.D, Is.EqualTo(modify ? 777u : ~15u), "the modification survived the eviction, or the original is intact");
            for (int i = 0; i < 2 * block; i++) stack.Pop(); // further down through the refills
            Assert.That(stack.Top.A, Is.EqualTo(7u), "the content below survived");
            Assert.That(stack.Top.D, Is.EqualTo(~7u));
        }
    }

    [Test]
    public void UpdateOfATopEntry_SurvivesSpillAndRefill()
    {
        const int block = 3;
        using BudgetedArena arena = Arena();
        using var stack = new SpillStack<TestEntry>(arena, Path.Combine(_dir, "u.spill"), block);
        for (uint i = 0; i < 4 * block; i++)
        {
            stack.Push(TestEntry.Of(i));
            stack.TopForUpdate().D = 7000 + i;
        }
        for (int i = 4 * block - 1; i >= 0; i--)
        {
            Assert.That(stack.Top.D, Is.EqualTo(7000u + (uint)i));
            stack.Pop();
        }
    }

    [Test]
    public void SpillFile_IsRemovedOnDispose_AndNoScratchRemainsTracked()
    {
        string path = Path.Combine(_dir, "d.spill");
        ScratchTracker.Reset();
        using (BudgetedArena arena = Arena())
        {
            var stack = new SpillStack<TestEntry>(arena, path, 4);
            for (uint i = 0; i < 100; i++) stack.Push(TestEntry.Of(i));
            Assert.That(File.Exists(path), Is.True);
            Assert.That(ScratchTracker.Live, Is.GreaterThan(0));
            stack.Dispose();
            Assert.That(File.Exists(path), Is.False);
            Assert.That(ScratchTracker.Live, Is.Zero);
        }
    }

    [Test]
    public void NoSpill_NoFile()
    {
        string path = Path.Combine(_dir, "n.spill");
        using BudgetedArena arena = Arena();
        using var stack = new SpillStack<TestEntry>(arena, path, 16);
        for (uint i = 0; i < 20; i++) stack.Push(TestEntry.Of(i));
        Assert.That(stack.Stats.BlocksWritten, Is.Zero);
        Assert.That(File.Exists(path), Is.False, "a stack that never leaves its two blocks creates no file");
    }

    [Test]
    public void TwoBlocksOfEntries_AreChargedToTheBudget()
    {
        var budget = new ArenaBudget(long.MaxValue);
        using (BudgetedArena arena = Arena(budget))
        using (var stack = new SpillStack<TestEntry>(arena, Path.Combine(_dir, "b.spill"), 1000))
        {
            Assert.That(budget.LiveBytes, Is.EqualTo(2L * 1000 * 16));
            stack.Push(TestEntry.Of(1));
        }
        Assert.That(budget.LiveBytes, Is.Zero);
    }

    [Test]
    public void SpillFlush_IsCountedAndBlockSizeBelowTwoIsRefused()
    {
        using BudgetedArena arena = Arena();
        using (var stack = new SpillStack<TestEntry>(arena, Path.Combine(_dir, "f.spill"), 4, flushBytes: 100))
        {
            for (uint i = 0; i < 80; i++) stack.Push(TestEntry.Of(i));
            Assert.That(stack.Stats.Flushes, Is.GreaterThan(0));
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => new SpillStack<TestEntry>(arena, Path.Combine(_dir, "x.spill"), 1));
    }

    [Test]
    public void EmptyStack_RefusesTopAndPop()
    {
        using BudgetedArena arena = Arena();
        using var stack = new SpillStack<TestEntry>(arena, Path.Combine(_dir, "e.spill"), 4);
        Assert.Throws<InvalidOperationException>(() => _ = stack.Top);
        Assert.Throws<InvalidOperationException>(stack.Pop);
        stack.Push(TestEntry.Of(1));
        stack.Pop();
        Assert.Throws<InvalidOperationException>(stack.Pop);
    }
}
