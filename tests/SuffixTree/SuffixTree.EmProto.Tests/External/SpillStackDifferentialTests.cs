using System.Runtime.InteropServices;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;
using Prod = SuffixTree.Esa.External;
using Proto = SuffixTree.EmProto.Ph5;
using ProtoBuild = SuffixTree.EmProto.Build;
using ProtoHarness = SuffixTree.EmProto.Harness;

namespace SuffixTree.EmProto.Tests.External;

[StructLayout(LayoutKind.Sequential)]
internal struct DiffEntry
{
    internal uint A;
    internal uint B;
    internal uint C;
    internal uint D;

    internal static DiffEntry Of(uint x) => new() { A = x, B = (x * 3) + 1, C = x ^ 0x5A5A5A5Au, D = ~x };
}

/// <summary>
/// G2 of work package 1b.7a, part 2: the production <c>SpillStack</c> behaves exactly like the prototype <c>SpillStack</c> on the same push / pop / update sequences, including every spill and restore: the same top,
/// height, peak, statistics at every step and the same bytes in the spill file.
/// </summary>
[TestFixture]
public class SpillStackDifferentialTests
{
    private static void Step(Proto.SpillStack<DiffEntry> proto, Prod.SpillStack<DiffEntry> prod, string label)
    {
        Assert.That(prod.Height, Is.EqualTo(proto.Height), $"{label}: height");
        Assert.That(prod.PeakHeight, Is.EqualTo(proto.PeakHeight), $"{label}: peak");
        Proto.SpillStats a = proto.Stats;
        Prod.SpillStats b = prod.Stats;
        Assert.That((b.BlocksWritten, b.BlocksRead, b.BytesWritten, b.BytesRead, b.PeakFileBytes, b.Flushes, b.CleanEvictions),
            Is.EqualTo((a.BlocksWritten, a.BlocksRead, a.BytesWritten, a.BytesRead, a.PeakFileBytes, a.Flushes, a.CleanEvictions)), $"{label}: stats");
        if (proto.Height > 0)
        {
            DiffEntry x = proto.Top;
            DiffEntry y = prod.Top;
            Assert.That((y.A, y.B, y.C, y.D), Is.EqualTo((x.A, x.B, x.C, x.D)), $"{label}: top");
        }
    }

    [Test]
    public void RandomSequences_EveryBlockSizeAndFlushInterval_AreIdentical_StepByStep_AndInTheSpillFile()
    {
        var rng = new SplitMix64(2026_11_18);
        int sequences = 0;
        foreach (int block in new[] { 2, 3, 4, 5, 8, 16, 33 })
        {
            foreach (long flush in new long[] { 0, 100, 1000 })
            {
                using var dir = new DiffDir();
                string pathA = dir.File("proto.spill");
                string pathB = dir.File("prod.spill");
                using var arenaA = new ProtoHarness.BudgetedArena(new ProtoHarness.ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));
                using var arenaB = new Prod.BudgetedArena(new Prod.ArenaBudget(long.MaxValue), ArenaOptions.Guarded(0xCD));
                using var proto = new Proto.SpillStack<DiffEntry>(arenaA, pathA, block, flush);
                var counters = new Prod.IoCounters();
                using var prod = new Prod.SpillStack<DiffEntry>(arenaB, pathB, block, flush, counters);
                uint next = 1;
                int phase = 0;
                for (int step = 0; step < 2500; step++)
                {
                    if (step % 250 == 0) phase = rng.NextInt(3);
                    int pushShare = phase == 0 ? 75 : phase == 1 ? 25 : 50;
                    bool push = proto.Height == 0 || rng.NextInt(100) < pushShare;
                    if (push)
                    {
                        DiffEntry e = DiffEntry.Of(next++);
                        proto.Push(e);
                        prod.Push(e);
                    }
                    else if (rng.NextInt(4) == 0)
                    {
                        uint v = next++;
                        proto.TopForUpdate().D = v;
                        prod.TopForUpdate().D = v;
                    }
                    else
                    {
                        proto.Pop();
                        prod.Pop();
                    }
                    Step(proto, prod, $"block {block} flush {flush} step {step}");
                    if (step % 500 == 499 && File.Exists(pathA))
                    {
                        Assert.That(File.Exists(pathB), Is.True);
                        Assert.That(DiffInputs.SameBytes(pathA, pathB), Is.True, $"block {block} flush {flush} step {step}: spill file");
                    }
                }
                Assert.That(File.Exists(pathA), Is.EqualTo(File.Exists(pathB)));
                Assert.That(counters.WriteBytes, Is.EqualTo(prod.Stats.BytesWritten));
                Assert.That(counters.ReadBytes, Is.EqualTo(prod.Stats.BytesRead));
                while (proto.Height > 0)
                {
                    proto.Pop();
                    prod.Pop();
                    Step(proto, prod, $"block {block} flush {flush} drain");
                }
                sequences++;
            }
        }
        Assert.That(sequences, Is.EqualTo(21));
    }

    [Test]
    public void AChain_ThreeHundredThousandEntries_IsIdentical_IncludingTheSpillBytes()
    {
        using var dir = new DiffDir();
        string pathA = dir.File("proto.spill");
        string pathB = dir.File("prod.spill");
        using var arenaA = new ProtoHarness.BudgetedArena(new ProtoHarness.ArenaBudget(long.MaxValue));
        using var arenaB = new Prod.BudgetedArena(new Prod.ArenaBudget(long.MaxValue));
        using var proto = new Proto.SpillStack<DiffEntry>(arenaA, pathA, 1024);
        using var prod = new Prod.SpillStack<DiffEntry>(arenaB, pathB, 1024);
        for (uint i = 0; i < 300_000; i++)
        {
            proto.Push(DiffEntry.Of(i));
            prod.Push(DiffEntry.Of(i));
        }
        Step(proto, prod, "after the chain");
        Assert.That(DiffInputs.SameBytes(pathA, pathB), Is.True, "the spill file of the chain");
        for (int i = 299_999; i >= 0; i--)
        {
            proto.Pop();
            prod.Pop();
            if (i % 1024 == 0) Step(proto, prod, $"pop {i}");
        }
        Step(proto, prod, "drained");
    }

    [Test]
    public void TheArenaBudget_IsChargedTheSameWay()
    {
        var budgetA = new ProtoHarness.ArenaBudget(long.MaxValue);
        var budgetB = new Prod.ArenaBudget(long.MaxValue);
        using var dir = new DiffDir();
        using (var arenaA = new ProtoHarness.BudgetedArena(budgetA))
        using (var arenaB = new Prod.BudgetedArena(budgetB))
        using (var a = new Proto.SpillStack<DiffEntry>(arenaA, dir.File("a"), 4096))
        using (var b = new Prod.SpillStack<DiffEntry>(arenaB, dir.File("b"), 4096))
        {
            Assert.That(budgetB.LiveBytes, Is.EqualTo(budgetA.LiveBytes));
            Assert.That(budgetB.PeakBytes, Is.EqualTo(budgetA.PeakBytes));
            Assert.That(a.BlockEntries, Is.EqualTo(b.BlockEntries));
        }
        Assert.That(budgetB.LiveBytes, Is.Zero);
    }
}
