using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The gaps that the audit of work package 1b.9c showed in the base case and in the statistics of PH1 (<c>External/Dcs/ExternalDcs.cs</c>; the mutants A05, A07 and L28 of <c>analysis/mutation/1b7d_dcs_externaldcs_listed.json</c>).
/// </summary>
[TestFixture]
[HangLimit(30)]
public class SaDcsGap1b9cTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    [Test]
    public void TheBaseCase_AllocatesNoMoreThanTheSlotsOfTheRankFile_WhenAChunkIsLarger()
    {
        // A05 (the chunk is not bounded by the slots): a chunk of 2^20 slots is 4 MiB, the rank file of this text has a few dozen slots; the arena of the base case must hold the slots and no more
        SaCase c = SaCase.Of(StringFamilies.AllSame(300));
        DifferenceCover cover = DifferenceCover.Create(13);
        var budget = new ArenaBudget(long.MaxValue);
        ExtDcsStats stats = ExternalDcsTests.RunAndCompare(_dir, c, cover, SaTiny.Dcs(_dir.Sub("scratch"), 2) with { BaseCase = 1 << 20, ChunkSlots = 1 << 20, Budget = budget }, "large chunk", out _);
        Assert.That(stats.Levels[0].BaseCase, Is.True, "the reduced string is sorted in RAM");
        Assert.That(budget.PeakBytes, Is.GreaterThan(0).And.LessThan(1 << 20), "no buffer of the base case is sized by the chunk");
    }

    [Test]
    public void TheBaseCaseOfLevelZero_ReportsWhatTheCreationOfTheRankFileDidAboutTheDiskSpace()
    {
        // A07 (the reservation is not recorded): the other tests of the reservation end level 0 in the sink of unique names; here the names are shared and the reduced string is small, so the rank file is created by the base case
        SaCase c = SaCase.Of(StringFamilies.AllSame(300));
        DifferenceCover cover = DifferenceCover.Create(13);
        var reservation = new FakeReservationSource { CanPreallocate = true };
        ExtDcsStats stats = ExternalDcsTests.RunAndCompare(_dir, c, cover, SaTiny.Dcs(_dir.Sub("scratch"), 2) with { BaseCase = 1 << 20, Files = new MappedFileOptions { Source = reservation } }, "base case", out _);
        Assert.That(stats.Levels[0].BaseCase, Is.True);
        Assert.That(reservation.PreallocateCalls, Is.EqualTo(1), "one file, the rank file, was created");
        Assert.That(stats.Reservation, Is.EqualTo(ReservationOutcome.Preallocated));
    }

    [Test]
    public void EveryLevel_CountsTheSecondsOfItsNamingAndOfItsMerge_AndEveryLevelBelowItsOwnWork()
    {
        // L28 (the naming statistics are not timed): the three phases of a level are measured one after the other
        SaCase c = SaCase.Of(StringFamilies.Fibonacci(6000));
        DifferenceCover cover = DifferenceCover.Create(13);
        ExtDcsStats stats = ExternalDcsTests.RunAndCompare(_dir, c, cover, SaTiny.Dcs(_dir.Sub("scratch"), 4) with { BaseCase = 8 }, "levels", out _);
        Assert.That(stats.Levels, Has.Count.GreaterThanOrEqualTo(2), "the recursion went below level 0");
        foreach (ExtDcsLevelStats level in stats.Levels)
        {
            Assert.That(level.NamingSeconds, Is.GreaterThan(0), $"level {level.Level}: naming");
            Assert.That(level.MergeSeconds, Is.GreaterThan(0), $"level {level.Level}: merge");
            if (!level.NamesUnique) Assert.That(level.BelowSeconds, Is.GreaterThan(0), $"level {level.Level}: the work below the level");
        }
    }

    [Test]
    public void TheRankTable_MapsEveryByteOfItsFile_EvenTheLastFourThatReachIntoAPageOfTheirOwn()
    {
        // R18 (the view is not read-only sized to the file): a read through the pointer of the view cannot tell a view that is four bytes short when those bytes share the last page; here the file is two pages and four bytes, so the
        // last slot is alone in the third page and the mapping of the file must be three pages long
        if (!OperatingSystem.IsLinux() || !File.Exists("/proc/self/maps")) Assert.Ignore("needs /proc/self/maps");
        DifferenceCover cover = DifferenceCover.Create(7); // 3 residues
        const long n = 4774; // 4774 / 7 + 1 = 683 blocks, 2049 slots, 8196 bytes
        int page = Environment.SystemPageSize;
        if (page != 4096) Assert.Ignore("the file is built for pages of 4,096 bytes");
        string path = _dir.File("ranks.bin");
        byte[] data = new byte[8196];
        data[8195] = 0x7F;
        File.WriteAllBytes(path, data);
        using RankTable table = RankTable.Open(path, cover, n);
        string full = Path.GetFullPath(path);
        long mapped = 0;
        foreach (string line in File.ReadLines("/proc/self/maps").Where(l => l.EndsWith(full, StringComparison.Ordinal)))
        {
            string[] range = line.Split(' ')[0].Split('-');
            mapped += Convert.ToInt64(range[1], 16) - Convert.ToInt64(range[0], 16);
        }
        Assert.That(mapped, Is.EqualTo(3L * page), "two pages and four bytes are three pages");
        Assert.That(table.Slots, Is.EqualTo(2049));
    }
}
