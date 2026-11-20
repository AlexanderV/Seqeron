using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The gap that the one-by-one re-run of the mutants that a watchdog or a crash had killed showed in <c>External/Dcs/NameSinks.cs</c> (N19 of <c>analysis/mutation/1b7d_dcs_namesinks_listed.json</c>): a buffer one name short is
/// an overrun of native memory that, in the production build, corrupts the neighbour at random (the earlier audit saw a hang or a crash) and in the quiet re-run did nothing at all, so no assertion noticed it.
/// </summary>
[TestFixture]
public class SaNamesGap1b9cTests
{
    [Test]
    public void TheReducedFileSink_HoldsOneBufferOfTheFullSizePerClass_AndWritesTheNamesClassMajor()
    {
        // N19 (the buffers are one name short): the arena of a test has guard bands whose canaries show a write past the end of a block, and it books exactly the bytes that were asked for; the last class is the one that runs past
        // the short block (its buffer ends the block), so a buffer of 4 names and 3 classes is filled to its end
        using var dir = new TempDir();
        DifferenceCover cover = DifferenceCover.Create(7);
        long n = 100;
        var g = new SampleGeometry(n, cover);
        const int ioBufferBytes = 16; // 4 names per class
        var budget = new ArenaBudget(long.MaxValue);
        using var arena = new BudgetedArena(budget, ArenaOptions.Guarded(0xCD));
        var counters = new IoCounters();
        var ledger = new ScratchLedger(counters);
        string path = dir.File("red.bin");
        uint[] names = [.. Enumerable.Range(0, (int)g.Samples).Select(i => (uint)(1000 + i))];
        using (var sink = new ReducedFileSink(path, g, ioBufferBytes, arena, ledger, counters))
        {
            Assert.That(budget.LiveBytes, Is.GreaterThanOrEqualTo((long)(ioBufferBytes * g.Classes)), "one buffer of four names for each of the three classes");
            sink.Accept(0, names);
            sink.Complete();
        }
        Assert.That(() => arena.VerifyCanaries(), Throws.Nothing, "nothing was written outside the block");

        // the file: class-major, big endian, a zero separator behind each class, the name of slot s at index s / classes of the class s % classes
        byte[] file = File.ReadAllBytes(path);
        for (int c = 0; c < g.Classes; c++)
        {
            for (long k = 0; k < g.ClassCount[c]; k++)
            {
                long slot = (k * g.Classes) + c;
                uint got = BinaryPrimitives.ReadUInt32BigEndian(file.AsSpan((int)((g.RedStart[c] + k) * 4), 4));
                Assert.That(got, Is.EqualTo(names[slot]), $"class {c}, member {k} (slot {slot})");
            }
            Assert.That(BinaryPrimitives.ReadUInt32BigEndian(file.AsSpan((int)((g.RedStart[c] + g.ClassCount[c]) * 4), 4)), Is.Zero, $"the separator of class {c}");
        }
    }
}
