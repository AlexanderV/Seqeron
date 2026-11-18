using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>An entry as the interval scan keeps it: (l, lb, first, last), 16 bytes.</summary>
[StructLayout(LayoutKind.Sequential)]
internal struct Entry16
{
    internal uint A;
    internal uint B;
    internal uint C;
    internal uint D;

    internal static Entry16 Of(uint x) => new() { A = x, B = (x * 3) + 1, C = x ^ 0x5A5A5A5Au, D = ~x };
}

/// <summary>
/// The external stack of PH5 (ADR-004 section 3.2): the top two blocks in RAM, older blocks appended to a spill file, read back backwards when a pop leaves the window, at least one full block hot. Oracle: <see cref="Stack{T}"/>.
/// </summary>
[TestFixture]
public class SpillStackTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private static BudgetedArena Arena(ArenaBudget? budget = null) => new(budget ?? new ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));

    private string Spill(string name = "s.spill") => _dir.File(name);

    [Test]
    public void ModelFuzz_EveryBlockSize_TopHeightAndPeakFollowTheModel()
    {
        var rng = new SplitMix64(12345);
        for (int block = 2; block <= 50; block++)
        {
            using BudgetedArena arena = Arena();
            using var stack = new SpillStack<Entry16>(arena, Spill($"s{block}.spill"), block);
            var model = new Stack<Entry16>();
            uint next = 1;
            long peak = 0;
            int phase = 0;
            for (int step = 0; step < 3000; step++)
            {
                if (step % 300 == 0) phase = rng.NextInt(3);
                int pushShare = phase == 0 ? 70 : phase == 1 ? 30 : 50;
                bool push = model.Count == 0 || rng.NextInt(100) < pushShare;
                if (push)
                {
                    Entry16 e = Entry16.Of(next++);
                    stack.Push(e);
                    model.Push(e);
                }
                else if (rng.NextInt(4) == 0)
                {
                    ref Entry16 top = ref stack.TopForUpdate();
                    top.D = next++;
                    Entry16 m = model.Pop();
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
                    Entry16 t = model.Peek();
                    Entry16 g = stack.Top;
                    Assert.That((g.A, g.B, g.C, g.D), Is.EqualTo((t.A, t.B, t.C, t.D)), $"block {block}, step {step}: top");
                }
            }
            Assert.That(stack.PeakHeight, Is.EqualTo(peak), $"block {block}: peak height");
            while (model.Count > 0)
            {
                Entry16 t = model.Pop();
                Assert.That(stack.Top.A, Is.EqualTo(t.A), $"block {block}: drain");
                stack.Pop();
            }
            Assert.That(stack.Height, Is.Zero);
            arena.VerifyCanaries();
        }
    }

    [Test]
    public void ExhaustiveSmallSequences_OfPushAndPop_EveryBlockSizeToFour()
    {
        // every sequence of 12 operations (push or pop, 2^12 sequences, the invalid pops skipped) against the model: all the transitions around the block boundaries
        for (int block = 2; block <= 4; block++)
        {
            for (int bits = 0; bits < (1 << 12); bits++)
            {
                using BudgetedArena arena = Arena();
                using var stack = new SpillStack<ulong>(arena, Spill($"x{block}.spill"), block);
                var model = new Stack<ulong>();
                ulong next = 1;
                for (int k = 0; k < 12; k++)
                {
                    if ((bits >> k & 1) == 1 || model.Count == 0)
                    {
                        stack.Push(next);
                        model.Push(next++);
                    }
                    else
                    {
                        Assert.That(stack.Top, Is.EqualTo(model.Peek()));
                        stack.Pop();
                        model.Pop();
                    }
                    Assert.That(stack.Height, Is.EqualTo(model.Count));
                    if (model.Count > 0) Assert.That(stack.Top, Is.EqualTo(model.Peek()), $"block {block} bits {bits} step {k}");
                }
            }
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
        using var stack = new SpillStack<Entry16>(arena, Spill(), block);
        for (uint i = 0; i < height; i++) stack.Push(Entry16.Of(i));
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
        Assert.That(stack.Stats.BytesRead, Is.EqualTo(expected * block * 16));
        Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(expected), "nothing is written twice when nothing is modified");
    }

    [Test]
    public void TheSpillFile_HoldsBlockJAtOffsetJTimesBlockBytes_AsRawEntries()
    {
        // the spill format of the prototype: block j at j x B x sizeof(T), no header, entries as they are in memory
        const int block = 4;
        using BudgetedArena arena = Arena();
        string path = Spill("format.spill");
        using var stack = new SpillStack<ulong>(arena, path, block);
        for (ulong i = 0; i < 3 * block + 1; i++) stack.Push(0x1000 + i); // top block 3: blocks 0 and 1 spilled
        Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(2));
        byte[] file = File.ReadAllBytes(path);
        Assert.That(file, Has.Length.EqualTo(2 * block * 8));
        Assert.That(MemoryMarshal.Cast<byte, ulong>(file).ToArray(), Is.EqualTo(Enumerable.Range(0, 2 * block).Select(i => 0x1000UL + (ulong)i).ToArray()));
        // push on: block 2 is evicted next, appended after block 1
        for (ulong i = 3 * block + 1; i < 4 * block + 1; i++) stack.Push(0x1000 + i);
        Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(3));
        Assert.That(new FileInfo(path).Length, Is.EqualTo(3 * block * 8));
        Assert.That(MemoryMarshal.Cast<byte, ulong>(File.ReadAllBytes(path)).ToArray(), Is.EqualTo(Enumerable.Range(0, 3 * block).Select(i => 0x1000UL + (ulong)i).ToArray()));
    }

    [Test]
    public void ABlockOfOneByteEntries_UsesTheEntrySize_NotSixteen()
    {
        using BudgetedArena arena = Arena();
        string path = Spill("bytes.spill");
        using var stack = new SpillStack<byte>(arena, path, 8);
        for (int i = 0; i < 40; i++) stack.Push((byte)i);
        Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(3));
        Assert.That(stack.Stats.BytesWritten, Is.EqualTo(24));
        Assert.That(File.ReadAllBytes(path), Is.EqualTo(Enumerable.Range(0, 24).Select(i => (byte)i).ToArray()));
        for (int i = 39; i >= 0; i--)
        {
            Assert.That(stack.Top, Is.EqualTo((byte)i));
            stack.Pop();
        }
    }

    [Test]
    public void Oscillation_AtEveryBlockBoundary_DoesNotThrash()
    {
        const int block = 8;
        foreach (int boundary in new[] { block, 2 * block, 3 * block, 4 * block })
        {
            using BudgetedArena arena = Arena();
            using var stack = new SpillStack<Entry16>(arena, Spill($"o{boundary}.spill"), block);
            for (uint i = 0; i < boundary; i++) stack.Push(Entry16.Of(i));
            SpillStats before = stack.Stats;
            for (int k = 0; k < 2000; k++)
            {
                stack.Push(Entry16.Of(1000 + (uint)k));
                stack.Pop();
                stack.Pop();
                stack.Push(Entry16.Of(uint.MaxValue - (uint)k));
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
            using var stack = new SpillStack<Entry16>(arena, Spill($"c{modify}.spill"), block);
            for (uint i = 0; i < 6 * block; i++) stack.Push(Entry16.Of(i)); // top block 5, blocks 0..3 spilled
            Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(4));
            for (int i = 0; i < 2 * block; i++) stack.Pop(); // block 3 is read back and becomes the top block: its file copy is current
            Assert.That(stack.Stats.BlocksRead, Is.EqualTo(1));
            Assert.That(stack.Top.A, Is.EqualTo(15u));
            if (modify) stack.TopForUpdate().D = 777;
            for (uint i = 0; i < 5; i++) stack.Push(Entry16.Of(100 + i)); // the fifth push moves on from block 4 and evicts block 3
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
            for (int i = 0; i < 2 * block; i++) stack.Pop();
            Assert.That(stack.Top.A, Is.EqualTo(7u), "the content below survived");
            Assert.That(stack.Top.D, Is.EqualTo(~7u));
        }
    }

    [Test]
    public void UpdateOfATopEntry_SurvivesSpillAndRefill()
    {
        const int block = 3;
        using BudgetedArena arena = Arena();
        using var stack = new SpillStack<Entry16>(arena, Spill(), block);
        for (uint i = 0; i < 4 * block; i++)
        {
            stack.Push(Entry16.Of(i));
            stack.TopForUpdate().D = 7000 + i;
        }
        for (int i = (4 * block) - 1; i >= 0; i--)
        {
            Assert.That(stack.Top.D, Is.EqualTo(7000u + (uint)i));
            stack.Pop();
        }
    }

    [Test]
    public void AChain_TheWorstCase_SpillsSixteenBytesPerEntry_AndRestoresEveryEntry()
    {
        // the a^n scan pushes one entry per suffix and pops nothing until the end: height n, spill 16 B x (height - ~2 blocks)
        const int block = 1024;
        const int height = 300_000;
        using BudgetedArena arena = Arena();
        var counters = new IoCounters();
        using var stack = new SpillStack<Entry16>(arena, Spill("chain.spill"), block, counters: counters);
        for (uint i = 0; i < height; i++) stack.Push(Entry16.Of(i));
        long expectedBlocks = ((height - 1) / block) - 1;
        Assert.That(stack.Stats.BlocksWritten, Is.EqualTo(expectedBlocks));
        Assert.That(stack.Stats.PeakFileBytes, Is.EqualTo(expectedBlocks * block * 16));
        Assert.That(stack.Stats.PeakFileBytes, Is.GreaterThan(16L * (height - (3L * block))), "16 B x (max LCP + 1) up to three blocks");
        Assert.That(counters.WriteBytes, Is.EqualTo(stack.Stats.BytesWritten));
        Assert.That(counters.ScratchPeakBytes, Is.EqualTo(stack.Stats.PeakFileBytes));
        for (int i = height - 1; i >= 0; i--)
        {
            Entry16 e = stack.Top;
            if (e.A != (uint)i || e.D != ~(uint)i) Assert.Fail($"entry {i} came back as {e.A}");
            stack.Pop();
        }
        Assert.That(stack.Stats.BlocksRead, Is.EqualTo(expectedBlocks));
        Assert.That(counters.ReadBytes, Is.EqualTo(stack.Stats.BytesRead));
    }

    [Test]
    public void TheSpill_IsRemovedOnDispose_AndNoScratchStaysAccounted()
    {
        string path = Spill("d.spill");
        var counters = new IoCounters();
        using BudgetedArena arena = Arena();
        var stack = new SpillStack<Entry16>(arena, path, 4, counters: counters);
        for (uint i = 0; i < 100; i++) stack.Push(Entry16.Of(i));
        Assert.That(File.Exists(path), Is.True);
        Assert.That(counters.ScratchLiveBytes, Is.GreaterThan(0));
        stack.Dispose();
        Assert.That(File.Exists(path), Is.False);
        Assert.That(counters.ScratchLiveBytes, Is.Zero);
        stack.Dispose(); // idempotent
    }

    [Test]
    public void NoSpill_NoFile_AndNoSpillFileIsCreatedByTheConstructor()
    {
        string path = Spill("n.spill");
        using BudgetedArena arena = Arena();
        using var stack = new SpillStack<Entry16>(arena, path, 16);
        Assert.That(File.Exists(path), Is.False, "nothing is created before the first spill");
        for (uint i = 0; i < 20; i++) stack.Push(Entry16.Of(i));
        Assert.That(stack.Stats.BlocksWritten, Is.Zero);
        Assert.That(File.Exists(path), Is.False, "a stack that never leaves its two blocks creates no file");
    }

    [Test]
    public void TwoBlocksOfEntries_AreChargedToTheBudget_AndTheCapIsEnforced()
    {
        var budget = new ArenaBudget(long.MaxValue);
        using (BudgetedArena arena = Arena(budget))
        using (var stack = new SpillStack<Entry16>(arena, Spill("b.spill"), 1000))
        {
            Assert.That(budget.LiveBytes, Is.EqualTo(2L * 1000 * 16));
            stack.Push(Entry16.Of(1));
        }
        Assert.That(budget.LiveBytes, Is.Zero);
        var tight = new ArenaBudget((2L * 1000 * 16) - 1);
        using BudgetedArena small = Arena(tight);
        Assert.Throws<ArenaBudgetExceededException>(() => _ = new SpillStack<Entry16>(small, Spill("t.spill"), 1000));
        var exact = new ArenaBudget(2L * 1000 * 16);
        using BudgetedArena fits = Arena(exact);
        using var ok = new SpillStack<Entry16>(fits, Spill("t2.spill"), 1000);
    }

    [Test]
    public void FlushInterval_IsCounted_AndMeansFsyncEveryNBytes()
    {
        using BudgetedArena arena = Arena();
        const int block = 4; // 64 bytes per block
        using var stack = new SpillStack<Entry16>(arena, Spill("f.spill"), block, flushBytes: 128);
        for (uint i = 0; i < 40 * block; i++) stack.Push(Entry16.Of(i));
        long blocks = stack.Stats.BlocksWritten;
        Assert.That(blocks, Is.EqualTo(38));
        Assert.That(stack.Stats.Flushes, Is.EqualTo(blocks / 2), "every second block of 64 bytes reaches 128");
        using BudgetedArena arena2 = Arena();
        using var never = new SpillStack<Entry16>(arena2, Spill("g.spill"), block, flushBytes: 0);
        for (uint i = 0; i < 40 * block; i++) never.Push(Entry16.Of(i));
        Assert.That(never.Stats.Flushes, Is.Zero, "0 = never");
    }

    [Test]
    public void Arguments_AreValidated()
    {
        using BudgetedArena arena = Arena();
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new SpillStack<Entry16>(arena, Spill("x.spill"), 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new SpillStack<Entry16>(arena, Spill("x.spill"), 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => _ = new SpillStack<Entry16>(arena, Spill("x.spill"), -5));
        Assert.Throws<ArgumentNullException>(() => _ = new SpillStack<Entry16>(null!, Spill("x.spill"), 4));
        using var ok = new SpillStack<Entry16>(arena, Spill("x.spill"), 2);
        Assert.That(ok.BlockEntries, Is.EqualTo(2));
    }

    [Test]
    public void EmptyStack_RefusesTopAndPop_AndPeakFollows()
    {
        using BudgetedArena arena = Arena();
        using var stack = new SpillStack<Entry16>(arena, Spill("e.spill"), 4);
        Assert.Throws<InvalidOperationException>(() => _ = stack.Top);
        Assert.Throws<InvalidOperationException>(() => stack.TopForUpdate());
        Assert.Throws<InvalidOperationException>(stack.Pop);
        stack.Push(Entry16.Of(1));
        stack.Pop();
        Assert.Throws<InvalidOperationException>(stack.Pop);
        Assert.Throws<InvalidOperationException>(() => _ = stack.Top);
        Assert.That(stack.PeakHeight, Is.EqualTo(1));
    }

    // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // failure modes: surfaced, never a hang, never silent garbage

    [Test]
    public void ATruncatedSpillFile_IsSurfaced_WhenABlockIsReadBack_AndTheStackIsLeftConsistent()
    {
        const int block = 4;
        string path = Spill("trunc.spill");
        using BudgetedArena arena = Arena();
        using var stack = new SpillStack<Entry16>(arena, path, block);
        for (uint i = 0; i < 6 * block; i++) stack.Push(Entry16.Of(i)); // blocks 0..3 spilled
        using (var f = new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite)) f.SetLength((3 * block * 16) + 5); // block 3 is cut
        // pop through the two resident blocks, then the read of block 3 hits the truncated tail
        for (int i = 0; i < (2 * block) - 1; i++) stack.Pop();
        Assert.That(stack.Top.A, Is.EqualTo(16u));
        Assert.That(stack.Height, Is.EqualTo(17));
        Assert.Throws<EndOfStreamException>(stack.Pop, "the entry that would leave the window needs block 3");
        Assert.That(stack.Height, Is.EqualTo(17), "a failed pop does not pop");
        Assert.That(stack.Top.A, Is.EqualTo(16u), "and the top entry is still the one it was");
        Assert.Throws<EndOfStreamException>(stack.Pop, "retrying fails the same way: no hang, no garbage");
        // restore the file: the same pop now succeeds (the failure left nothing half done)
        using (var f = new FileStream(path, FileMode.Open, FileAccess.ReadWrite, FileShare.ReadWrite)) f.SetLength(4 * block * 16);
        // the truncated bytes were lost: only the structure of the stack is checked, not the data of block 3
        Assert.DoesNotThrow(stack.Pop);
        Assert.That(stack.Height, Is.EqualTo(16));
    }

    [Test]
    public void AMissingSpillDirectory_IsSurfacedAtTheFirstSpill_AndThePushDoesNotHappen()
    {
        string path = Path.Combine(_dir.Path, "gone", "x.spill");
        using BudgetedArena arena = Arena();
        using var stack = new SpillStack<Entry16>(arena, path, 4);
        for (uint i = 0; i < 8; i++) stack.Push(Entry16.Of(i)); // two blocks: still nothing to spill
        Assert.Throws<DirectoryNotFoundException>(() => stack.Push(Entry16.Of(8)), "the ninth entry starts block 2 and needs block 0 on disk");
        Assert.That(stack.Height, Is.EqualTo(8));
        Assert.That(stack.Top.A, Is.EqualTo(7u));
        Assert.Throws<DirectoryNotFoundException>(() => stack.Push(Entry16.Of(8)), "a retry fails the same way");
        Directory.CreateDirectory(Path.GetDirectoryName(path)!);
        stack.Push(Entry16.Of(8)); // the directory is there now: the stack carries on
        Assert.That(stack.Height, Is.EqualTo(9));
        for (int i = 8; i >= 0; i--)
        {
            Assert.That(stack.Top.A, Is.EqualTo((uint)i));
            stack.Pop();
        }
    }

    [Test]
    public void ADiskWriteError_IsSurfaced_AsAnIoException_AndThePushDoesNotHappen()
    {
        if (!File.Exists("/dev/full") || !OperatingSystem.IsLinux()) Assert.Ignore("no /dev/full on this platform");
        string path = Spill("full.spill");
        File.CreateSymbolicLink(path, "/dev/full"); // a link, so that disposing the stack deletes the link and never the device
        using BudgetedArena arena = Arena();
        var stack = new SpillStack<Entry16>(arena, path, 4);
        uint pushed = 0;
        IOException? error = null;
        try
        {
            for (; pushed < 100; pushed++) stack.Push(Entry16.Of(pushed));
        }
        catch (IOException ex)
        {
            error = ex;
        }
        Assert.That(error, Is.Not.Null, "ENOSPC reaches the caller");
        Assert.That(stack.Height, Is.EqualTo(pushed), "the push that needed the spill did not happen");
        Assert.That(stack.Top.A, Is.EqualTo(pushed - 1), "the top is the last entry that did");
        Assert.Throws<IOException>(() => stack.Push(Entry16.Of(999)), "a retry fails the same way");
        Assert.That(stack.Height, Is.EqualTo(pushed));
        stack.Dispose();
        Assert.That(File.Exists("/dev/full"), Is.True, "the device is still there");
        Assert.That(File.Exists(path), Is.False, "the link is gone");
    }
}
