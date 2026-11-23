using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// PH1 of work package 1b.7d (the rank file of the difference-cover sample in external memory) against the definition: every rank is the 1-based rank of the sample suffix among the sample suffixes (plus |D| where the
/// recursion ran), on every string of the exhaustive tiers, on the families, for every thread count and bucket capacity, with the failure modes of the stage.
/// </summary>
[TestFixture]
[HangLimit(30)]
public class ExternalDcsTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private string Scratch => _dir.Sub("scratch");

    /// <summary>Runs PH1 on the case and compares the rank file with the oracle; returns the stats (and the bytes through <paramref name="rank"/>).</summary>
    internal static ExtDcsStats RunAndCompare(TempDir dir, SaCase c, DifferenceCover cover, ExtDcsOptions options, string label, out string rank, Action<MappedTextReader>? configure = null)
    {
        string textPath = SaRun.WriteText(dir, c, label.GetHashCode(StringComparison.Ordinal).ToString("X"));
        rank = dir.File("rank-" + Guid.NewGuid().ToString("N")[..8]);
        using (MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits))
        {
            configure?.Invoke(reader);
            ExtDcsStats stats = ExternalDcs.BuildRankFile(reader, cover, rank, options);
            uint[] expected = DcsOracle.ExpectedRanks(c.Symbols, cover, out bool unique, out long distinct, out long shared);
            Assert.That(new FileInfo(rank).Length, Is.EqualTo(4L * expected.Length), $"{label}: length of the rank file");
            uint[] got = FileBits.ReadUInts(rank);
            for (int k = 0; k < got.Length; k++)
            {
                if (got[k] != expected[k]) Assert.Fail($"{label}: rank[{k}] = {got[k]}, expected {expected[k]} (X = {cover.Period}, text '{(c.Text.Length > 50 ? c.Text[..50] + "..." : c.Text)}')");
            }
            Assert.That(stats.Levels[0].NamesUnique, Is.EqualTo(unique), $"{label}: names unique");
            Assert.That(stats.Levels[0].DistinctNames, Is.EqualTo(distinct), $"{label}: distinct names");
            Assert.That(stats.Levels[0].NonUniqueRecords, Is.EqualTo(shared), $"{label}: sample positions whose name is shared");
            Assert.That(stats.Levels[0].SampleSize, Is.EqualTo(cover.SampleCount(c.N)));
            Assert.That(stats.Levels[0].TextLength, Is.EqualTo(c.N));
            Assert.That(stats.Levels.Select(l => l.Level), Is.EqualTo(Enumerable.Range(0, stats.Levels.Count)), $"{label}: levels in order");
            Assert.That(options.Counters?.ScratchLiveBytes ?? 0, Is.Zero, $"{label}: no scratch byte left in the account");
            if (options.Budget is { } budget) Assert.That(budget.LiveBytes, Is.Zero, $"{label}: no arena byte left");
            return stats;
        }
    }

    private static void AssertScratchGone(string scratch, string label)
    {
        Assert.That(Directory.Exists(scratch), Is.False, $"{label}: the scratch directory is left behind");
        ProcFds.AssertClosed(label, [], [scratch]);
    }

    [Test]
    public void SmallTexts_AreTheOracle_WithRotatingTinyParameters_BothTheUniqueAndTheRecursiveCase()
    {
        int v = 0;
        foreach (string text in new[] { "a", "ab", "ba", "aa", "banana", "mississippi", "abracadabra", new string('a', 40), "abababababababababab", "abcabcabcabcabcabcabc", "aabaabaabaabaabaab", StringFamilies.Fibonacci(90), StringFamilies.ThueMorse(75), StringFamilies.RandomDna(120, 4), StringFamilies.RandomProtein(100, 5) })
        {
            SaCase c = SaCase.Of(text);
            foreach (int variant in new[] { v, v + 3, v + 8 })
            {
                using var dir = new TempDir();
                DifferenceCover cover = DifferenceCover.Create(SaTiny.Periods[variant % SaTiny.Periods.Length]);
                RunAndCompare(dir, c, cover, SaTiny.Dcs(dir.Sub("scratch"), variant, 1 + (variant % 3)), $"'{(text.Length > 24 ? text[..24] : text)}' variant {variant}", out _);
                AssertScratchGone(dir.Sub("scratch"), $"variant {variant}");
            }
            v++;
        }
    }

    [Test]
    public void TheEmptyText_AndTextsOfOneAndTwoSymbols_Work()
    {
        foreach (string text in new[] { string.Empty, "x", "xx", "xy", "yx" })
        {
            foreach (int x in new[] { 2, 3, 7, 993 })
            {
                SaCase c = SaCase.Of(text);
                using var dir = new TempDir();
                RunAndCompare(dir, c, DifferenceCover.Create(x), SaTiny.Dcs(dir.Sub("scratch"), x), $"'{text}' X = {x}", out string rank);
                Assert.That(File.Exists(rank));
            }
        }
    }

    [Test]
    [HangLimit(150)]
    public void EveryString_OverSigma1To3_OfTheDefaultTier_IsTheOracle()
    {
        int count = 0;
        foreach ((int sigma, int max) in new[] { (1, 14), (2, 9), (3, 5) })
        {
            foreach (string s in ExhaustiveStrings.All(sigma, max))
            {
                if (s.Length == 0) continue;
                SaCase c = SaCase.Of(s);
                using var dir = new TempDir();
                int v = count++;
                DifferenceCover cover = DifferenceCover.Create(SaTiny.Periods[v % SaTiny.Periods.Length]);
                RunAndCompare(dir, c, cover, SaTiny.Dcs(dir.Sub("scratch"), v, 1 + (v % 3)), $"'{s}' variant {v}", out _);
            }
        }
        Assert.That(count, Is.EqualTo(14 + 1022 + 363));
    }

    [Test]
    public void Families_AtEveryTinyVariant_AreTheOracle_AndTheRecursionReachesLevelTwo()
    {
        int v = 0;
        int deepest = 0;
        bool baseCaseSeen = false;
        foreach (int n in new[] { 700, 3000 })
        {
            foreach ((string name, string text) in SaRun.Families(n))
            {
                SaCase c = SaCase.Of(text);
                using var dir = new TempDir();
                DifferenceCover cover = DifferenceCover.Create(SaTiny.Periods[v % SaTiny.Periods.Length]);
                ExtDcsStats stats = RunAndCompare(dir, c, cover, SaTiny.Dcs(dir.Sub("scratch"), v, 1 + (v % 4)), $"{name} {n} variant {v}", out _);
                deepest = Math.Max(deepest, stats.MaxLevel);
                baseCaseSeen |= stats.Levels[^1].BaseCase;
                if (stats.MaxLevel > 0) Assert.That(stats.Levels.Take(stats.MaxLevel).All(l => !l.NamesUnique), Is.True, $"{name}: a level below exists only where the names are not unique");
                for (int i = 0; i < stats.Levels.Count - 1; i++) Assert.That(stats.Levels[i + 1].TextLength, Is.EqualTo(new SampleGeometry(stats.Levels[i].TextLength, cover).RedLength), "the text of a level is the reduced string above it");
                v++;
            }
        }
        Assert.That(deepest, Is.GreaterThanOrEqualTo(2), "a small N_base forces levels 1 and 2 on the repetitive families");
        Assert.That(baseCaseSeen, Is.True);
    }

    [Test]
    public void ARepetitiveText_WithASmallBaseCase_RecursesSeveralLevels_EachShorterThanTheOneAbove()
    {
        SaCase c = SaCase.Of(new string('a', 4000));
        DifferenceCover cover = DifferenceCover.Create(7);
        ExtDcsOptions o = SaTiny.Dcs(Scratch, 1) with { BaseCase = 12 };
        ExtDcsStats stats = RunAndCompare(_dir, c, cover, o, "a^n", out _);
        Assert.That(stats.MaxLevel, Is.GreaterThanOrEqualTo(4));
        for (int i = 1; i < stats.Levels.Count; i++) Assert.That(stats.Levels[i].TextLength, Is.LessThan(stats.Levels[i - 1].TextLength));
        Assert.That(stats.Levels[^1].BaseCase || stats.Levels[^1].NamesUnique, Is.True, "the recursion ends in the base case or where the names are unique");
        Assert.That(stats.Levels.Take(stats.Levels.Count - 1).All(l => !l.BaseCase), Is.True, "only the last level is the base case");
    }

    [Test]
    public void UniqueNames_NeedNoLevelBelow_AndNoBaseCase()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(5000, 9));
        ExtDcsStats stats = RunAndCompare(_dir, c, DifferenceCover.Create(31), SaTiny.Dcs(Scratch, 2), "random dna", out _);
        Assert.That(stats.MaxLevel, Is.Zero);
        Assert.That(stats.Levels[0].NamesUnique, Is.True);
        Assert.That(stats.Levels[0].BaseCase, Is.False);
        Assert.That(stats.Levels[0].NonUniqueRecords, Is.Zero);
        Assert.That(stats.Levels[0].DistinctNames, Is.EqualTo(stats.Levels[0].SampleSize));
    }

    [Test]
    public void ThePeriodTwo_MakesTheReducedStringLongerThanTheText_AndTheBaseCaseTakesItWithoutRecursion()
    {
        foreach (int x in new[] { 2 })
        {
            SaCase c = SaCase.Of(new string('a', 30));
            DifferenceCover cover = DifferenceCover.Create(x);
            ExtDcsStats stats = RunAndCompare(_dir, c, cover, SaTiny.Dcs(Scratch, 1) with { BaseCase = 4 }, $"X = {x}", out _);
            Assert.That(stats.MaxLevel, Is.Zero, $"X = {x}: the reduced string is not shorter than the text, so it is sorted in RAM at once");
            Assert.That(stats.Levels[0].BaseCase, Is.True);
            Assert.That(new SampleGeometry(30, cover).RedLength, Is.GreaterThanOrEqualTo(30));
        }
    }

    [Test]
    [HangLimit(290)]
    public void TheOutput_IsIdentical_ForEveryThreadCount_BucketCapacity_AndSeed()
    {
        foreach ((string name, string text) in SaRun.Families(12_000).Where(f => f.Name is "a^n" or "fibonacci" or "tandem" or "random-dna" or "genome-like"))
        {
            SaCase c = SaCase.Of(text);
            DifferenceCover cover = DifferenceCover.Create(31);
            string? reference = null;
            int v = 0;
            foreach (int threads in new[] { 1, 2, 4 })
            {
                foreach (long ks in new long[] { 7, 150, 4000, 1 << 20 })
                {
                    using var dir = new TempDir();
                    ExtDcsOptions o = SaTiny.Dcs(dir.Sub("scratch"), v, threads) with { BucketCapacity = ks, BaseCase = 25 + (v % 5 * 100), Seed = (ulong)(v + 100) };
                    if (reference is null)
                    {
                        RunAndCompare(dir, c, cover, o, $"{name} first", out string rank);
                        reference = SaRun.Sha(rank);
                    }
                    else
                    {
                        string textPath = SaRun.WriteText(dir, c);
                        string rank = dir.File("rank");
                        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
                        ExternalDcs.BuildRankFile(reader, cover, rank, o);
                        Assert.That(SaRun.Sha(rank), Is.EqualTo(reference), $"{name}: {threads} threads, K_s {ks}");
                    }
                    v++;
                }
            }
        }
    }

    [Test]
    public void TwoRuns_GiveTheSameBytes_AndTheSplitsAreMade()
    {
        SaCase c = SaCase.Of(StringFamilies.Fibonacci(6000));
        DifferenceCover cover = DifferenceCover.Create(13);
        ExtDcsStats first = RunAndCompare(_dir, c, cover, SaTiny.Dcs(Scratch, 4, 3), "first", out string a);
        using var dir2 = new TempDir();
        ExtDcsStats second = RunAndCompare(dir2, c, cover, SaTiny.Dcs(dir2.Sub("scratch"), 4, 3), "second", out string b);
        Assert.That(SaRun.Sha(a), Is.EqualTo(SaRun.Sha(b)));
        Assert.That(first.Levels[0].SplitScans, Is.GreaterThan(0), "a small K_s splits the root");
        Assert.That(first.Levels[0].Leaves, Is.EqualTo(second.Levels[0].Leaves), "the shape of the bucket tree is a function of the text and the seed");
        Assert.That(first.Levels[0].Scans, Is.GreaterThanOrEqualTo(first.Levels[0].Leaves + first.Levels[0].SplitScans));
    }

    [Test]
    [HangLimit(40)]
    public void ALeafThatOutgrowsItsCapacity_IsSplitAgain_AndTheRanksAreStillRight()
    {
        // the splitters come from a reservoir, so a bucket can hold more than K_s records: it overflows during its scan and is split once more
        long overflows = 0;
        int v = 0;
        foreach ((string name, string text) in SaRun.Families(1500))
        {
            SaCase c = SaCase.Of(text);
            foreach (int reservoir in new[] { 3, 4, 8 })
            {
                using var dir = new TempDir();
                ExtDcsOptions o = SaTiny.Dcs(dir.Sub("scratch"), v++, 1) with { ReservoirSize = reservoir, BucketCapacity = 5, FanOut = 2 };
                ExtDcsStats stats = RunAndCompare(dir, c, DifferenceCover.Create(9), o, $"{name} reservoir {reservoir}", out _);
                foreach (ExtDcsLevelStats l in stats.Levels) overflows += l.Scans - l.Leaves - l.SplitScans;
            }
        }
        Assert.That(overflows, Is.GreaterThan(0), "at least one bucket overflowed its capacity and was queued for another split");
    }

    [Test]
    public void RunsGoToScratchFiles_WhenTheSpillSizeIsZero_AndStayInMemoryWhenItIsLarge()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(4000, 3) + new string('a', 800));
        DifferenceCover cover = DifferenceCover.Create(13);
        var counters = new IoCounters();
        ExtDcsStats files = RunAndCompare(_dir, c, cover, SaTiny.Dcs(Scratch, 2) with { RunSpillBytes = 0, Counters = counters }, "files", out _);
        Assert.That(files.Levels[0].FileRuns, Is.GreaterThan(0).And.LessThanOrEqualTo(files.Levels[0].Leaves), "every non-empty bucket run is a file");
        Assert.That(counters.ScratchPeakBytes, Is.GreaterThan(0));
        Assert.That(counters.WriteBytes, Is.GreaterThan(0));
        Assert.That(counters.ReadBytes, Is.GreaterThan(0));
        using var dir2 = new TempDir();
        var counters2 = new IoCounters();
        ExtDcsStats memory = RunAndCompare(dir2, c, cover, SaTiny.Dcs(dir2.Sub("scratch"), 2) with { RunSpillBytes = 1 << 30, Counters = counters2 }, "memory", out _);
        Assert.That(memory.Levels[0].FileRuns, Is.Zero);
        Assert.That(counters2.WriteBytes, Is.LessThan(counters.WriteBytes));
    }

    [Test]
    public void TheRankFile_IsCreatedAtItsFinalLength_ThroughTheReservation_AndWrittenThroughWindows()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3000, 6));
        DifferenceCover cover = DifferenceCover.Create(13);
        var reservation = new FakeReservationSource { CanPreallocate = true };
        var counters = new IoCounters();
        ExtDcsStats stats = RunAndCompare(_dir, c, cover, SaTiny.Dcs(Scratch, 1) with { Files = new MappedFileOptions { Source = reservation }, Counters = counters, WindowBytes = 128 }, "reserve", out string rank);
        Assert.That(reservation.PreallocateCalls, Is.EqualTo(1), "one file, the rank file, was created");
        Assert.That(stats.Reservation, Is.EqualTo(ReservationOutcome.Preallocated));
        Assert.That(counters.MappedWriteBytes, Is.EqualTo(4 * stats.Levels[0].SampleSize), "every real slot of the rank file went through the windows, once (the unused slots of the last block stay zero)");
        Assert.That(counters.MappedMaxWindowBytes, Is.EqualTo(128));
    }

    [Test]
    public void TheCounters_OfARecursiveRun_CountTheScratchAndGiveItBack()
    {
        SaCase c = SaCase.Of(new string('a', 3000));
        var counters = new IoCounters();
        ExtDcsStats stats = RunAndCompare(_dir, c, DifferenceCover.Create(7), SaTiny.Dcs(Scratch, 1) with { Counters = counters, BaseCase = 20 }, "a^n", out _);
        Assert.That(stats.MaxLevel, Is.GreaterThanOrEqualTo(2));
        Assert.That(counters.ScratchLiveBytes, Is.Zero);
        Assert.That(counters.ScratchPeakBytes, Is.GreaterThan(4 * stats.Levels[1].TextLength), "the reduced string of level 1 is a scratch file of 4 bytes per symbol");
        Assert.That(counters.ReadBytes, Is.GreaterThan(0));
    }

    [Test]
    public void ThePrivateCounters_AreNotRequired()
    {
        SaCase c = SaCase.Of("mississippi");
        Assert.DoesNotThrow(() => RunAndCompare(_dir, c, DifferenceCover.Create(7), SaTiny.Dcs(Scratch, 1) with { Counters = null, Files = null, Budget = null, Trace = null }, "defaults", out _));
    }

    // ---- failure modes: no rank file, no scratch, no descriptor, no mapping

    private void AssertFails(SaCase c, DifferenceCover cover, ExtDcsOptions o, Type exception, string label, string? rank = null)
    {
        rank ??= _dir.File("fail.rank");
        string textPath = SaRun.WriteText(_dir, c, "f");
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        Exception? thrown = null;
        try
        {
            ExternalDcs.BuildRankFile(reader, cover, rank, o);
        }
        catch (Exception e)
        {
            thrown = e;
        }
        Assert.That(thrown, Is.Not.Null, $"{label}: the build must fail");
        Assert.That(thrown, Is.InstanceOf(exception), $"{label}: {thrown}");
        Assert.That(File.Exists(rank), Is.False, $"{label}: the rank file is left behind");
        AssertScratchGone(o.ScratchDir, label);
        ProcFds.AssertClosed(label, [rank], [o.ScratchDir]);
        Assert.That(ProcMaps.Contains(rank) ?? false, Is.False, $"{label}: the rank file is still mapped");
        if (o.Budget is { } b) Assert.That(b.LiveBytes, Is.Zero, $"{label}: arena bytes are still reserved");
        if (o.Counters is { } k) Assert.That(k.ScratchLiveBytes, Is.Zero, $"{label}: scratch bytes are still accounted");
    }

    [Test]
    public void AnArenaBudgetThatIsTooSmall_EndsWithTheBudgetException_AtEveryStage()
    {
        SaCase c = SaCase.Of(new string('a', 600));
        DifferenceCover cover = DifferenceCover.Create(7);
        foreach (long cap in new long[] { 0, 50, 95 })
        {
            using var dir = new TempDir();
            string textPath = SaRun.WriteText(dir, c);
            ExtDcsOptions o = SaTiny.Dcs(dir.Sub("scratch"), 3) with { Budget = new ArenaBudget(cap), BaseCase = 30 };
            string rank = dir.File("r");
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            Assert.Throws<ArenaBudgetExceededException>(() => ExternalDcs.BuildRankFile(reader, cover, rank, o), $"cap {cap}");
            Assert.That(File.Exists(rank), Is.False);
            Assert.That(Directory.Exists(o.ScratchDir), Is.False);
            Assert.That(o.Budget!.LiveBytes, Is.Zero);
        }
    }

    [Test]
    public void ASpillThatFails_IsAnIoException_NotWrapped_AndNothingIsLeft()
    {
        SaCase c = SaCase.Of(StringFamilies.Fibonacci(2500));
        DifferenceCover cover = DifferenceCover.Create(13);
        foreach (string stage in new[] { "naming", "merge" })
        {
            using var dir = new TempDir();
            string scratch = dir.Sub("scratch");
            bool done = false;
            ExtDcsOptions o = SaTiny.Dcs(scratch, 5, 2) with
            {
                RunSpillBytes = 0,
                Trace = step =>
                {
                    if (step == stage && !done)
                    {
                        done = true;
                        if (Directory.Exists(scratch)) Directory.Delete(scratch, true); // no file can be created in it any more
                    }
                },
            };
            string textPath = SaRun.WriteText(dir, c);
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            Assert.Throws(Is.InstanceOf<IOException>(), () => ExternalDcs.BuildRankFile(reader, cover, dir.File("r"), o), stage);
            Assert.That(File.Exists(dir.File("r")), Is.False, stage);
            Assert.That(Directory.Exists(scratch), Is.False, stage);
            ProcFds.AssertClosed(stage, [dir.File("r")], [scratch]);
        }
    }

    [Test]
    public void AHookThatThrowsInAWorker_ReachesTheCallerWithItsOwnType()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(1500, 2) + new string('a', 500));
        DifferenceCover cover = DifferenceCover.Create(13);
        foreach (int threads in new[] { 1, 2, 4 })
        {
            using var dir = new TempDir();
            int calls = 0;
            ExtDcsOptions o = SaTiny.Dcs(dir.Sub("scratch"), 2, threads) with { Trace = step => { if (step == "naming" && Interlocked.Increment(ref calls) == 3) throw new FormatException("hook"); } };
            AssertFailsIn(dir, c, cover, o, typeof(FormatException), $"{threads} threads");
        }
    }

    private static void AssertFailsIn(TempDir dir, SaCase c, DifferenceCover cover, ExtDcsOptions o, Type exception, string label)
    {
        string textPath = SaRun.WriteText(dir, c);
        string rank = dir.File("r");
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        Exception? e = Assert.Catch(() => ExternalDcs.BuildRankFile(reader, cover, rank, o), label);
        Assert.That(e, Is.InstanceOf(exception), $"{label}: {e}");
        Assert.That(e!.GetType(), Is.EqualTo(exception), $"{label}: the exception is not wrapped");
        Assert.That(File.Exists(rank), Is.False, label);
        Assert.That(Directory.Exists(o.ScratchDir), Is.False, label);
        ProcFds.AssertClosed(label, [rank], [o.ScratchDir]);
    }

    [Test]
    public void AFullDisk_AtTheCreationOfTheRankFile_OrOfAScratchFileOfALevel_IsAnIoException()
    {
        SaCase c = SaCase.Of(new string('a', 700));
        DifferenceCover cover = DifferenceCover.Create(7);
        var counting = new FailingAfter(int.MaxValue);
        using (var dir = new TempDir())
        {
            ExtDcsOptions o = SaTiny.Dcs(dir.Sub("scratch"), 1) with { BaseCase = 30, Files = new MappedFileOptions { Source = counting } };
            RunAndCompare(dir, c, cover, o, "counting", out _);
        }
        int creations = counting.Calls;
        Assert.That(creations, Is.GreaterThanOrEqualTo(3), "the rank file and the suffix array and the rank file of every level below");
        for (int succeed = 0; succeed < creations; succeed++)
        {
            using var dir = new TempDir();
            ExtDcsOptions o = SaTiny.Dcs(dir.Sub("scratch"), 1) with { BaseCase = 30, Files = new MappedFileOptions { Source = new FailingAfter(succeed) } };
            AssertFails(c, cover, o, typeof(IOException), $"creation {succeed + 1} of {creations} fails", dir.File("r"));
        }
    }

    [Test]
    public void ASimulatedPlatformWithoutFallocate_ChecksTheFreeSpace_AndAFullDiskIsStillAnIoException()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(800, 2));
        DifferenceCover cover = DifferenceCover.Create(13);
        foreach ((long free, ReservationOutcome expected) in new[] { (-1L, ReservationOutcome.Unchecked), (1L << 40, ReservationOutcome.FreeSpaceChecked) })
        {
            using var dir = new TempDir();
            var source = new FakeReservationSource { CanPreallocate = false, Free = free };
            ExtDcsStats stats = RunAndCompare(dir, c, cover, SaTiny.Dcs(dir.Sub("scratch"), 1) with { Files = new MappedFileOptions { Source = source } }, $"free {free}", out _);
            Assert.That(stats.Reservation, Is.EqualTo(expected));
            Assert.That(source.PreallocateCalls, Is.Zero);
            Assert.That(source.FreeCalls, Is.GreaterThanOrEqualTo(1));
        }
        using var dir2 = new TempDir();
        ExtDcsOptions full = SaTiny.Dcs(dir2.Sub("scratch"), 1) with { Files = new MappedFileOptions { Source = new FakeReservationSource { CanPreallocate = false, Free = 100 } } };
        AssertFails(c, cover, full, typeof(IOException), "full disk", dir2.File("r"));
    }

    [Test]
    public void TheOptions_AreChecked()
    {
        SaCase c = SaCase.Of("banana");
        DifferenceCover cover = DifferenceCover.Create(7);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        ExtDcsOptions ok = SaTiny.Dcs(Scratch, 1);
        foreach ((string what, ExtDcsOptions bad) in new (string, ExtDcsOptions)[]
        {
            ("bucket capacity 1", ok with { BucketCapacity = 1 }),
            ("bucket capacity too large", ok with { BucketCapacity = (int.MaxValue / 2) + 1L }),
            ("base case 0", ok with { BaseCase = 0 }),
            ("threads 0", ok with { Threads = 0 }),
            ("chunk slots 0", ok with { ChunkSlots = 0 }),
            ("fan-out 1", ok with { FanOut = 1 }),
            ("reservoir 2", ok with { ReservoirSize = 2 }),
            ("io buffer 8", ok with { IoBufferBytes = 8 }),
            ("window 0", ok with { WindowBytes = 0 }),
            ("merge budget 4", ok with { MergeBufferBudget = 4 }),
            ("sort capacity 1", ok with { SortLeafCapacity = 1 }),
            ("run spill -1", ok with { RunSpillBytes = -1 }),
        })
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => ExternalDcs.BuildRankFile(reader, cover, _dir.File("x"), bad), what);
            Assert.That(File.Exists(_dir.File("x")), Is.False, what);
            Assert.That(Directory.Exists(Scratch), Is.False, $"{what}: nothing was created");
        }
        Assert.Throws<ArgumentNullException>(() => ExternalDcs.BuildRankFile(null!, cover, "x", ok));
        Assert.Throws<ArgumentNullException>(() => ExternalDcs.BuildRankFile(reader, null!, "x", ok));
        Assert.Throws<ArgumentNullException>(() => ExternalDcs.BuildRankFile(reader, cover, null!, ok));
        Assert.Throws<ArgumentNullException>(() => ExternalDcs.BuildRankFile(reader, cover, "x", null!));
    }

    [Test]
    public void AScratchDirectoryThatExists_IsKept_WithItsForeignFile_WhenTheBuildEnds()
    {
        SaCase c = SaCase.Of(StringFamilies.Fibonacci(800));
        Directory.CreateDirectory(Scratch);
        File.WriteAllText(Path.Combine(Scratch, "foreign.txt"), "mine");
        RunAndCompare(_dir, c, DifferenceCover.Create(9), SaTiny.Dcs(Scratch, 2) with { RunSpillBytes = 0 }, "foreign", out _);
        Assert.That(Directory.GetFileSystemEntries(Scratch).Select(Path.GetFileName), Is.EqualTo(new[] { "foreign.txt" }), "only our files were removed");
    }

    [Test]
    public void ARankFileOfTheBuild_IsAValidTable_AndTheTableVerifiesIt()
    {
        foreach ((string name, string text) in SaRun.Families(2500))
        {
            SaCase c = SaCase.Of(text);
            DifferenceCover cover = DifferenceCover.Create(13);
            using var dir = new TempDir();
            RunAndCompare(dir, c, cover, SaTiny.Dcs(dir.Sub("scratch"), 6), name, out string rank);
            using RankTable t = RankTable.Open(rank, cover, c.N);
            Assert.DoesNotThrow(t.Verify, name);
        }
    }

    [Test]
    public void TheLeafBuffers_AreNoLargerThanTheSample_WhateverK_sIs()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(3000, 6));
        var budget = new ArenaBudget(long.MaxValue);
        ExtDcsStats stats = RunAndCompare(_dir, c, DifferenceCover.Create(13), new ExtDcsOptions { ScratchDir = Scratch, Budget = budget, ArenaGuards = null }, "default K_s", out _);
        Assert.That(stats.Levels[0].SampleSize, Is.LessThan(1500));
        Assert.That(budget.PeakBytes, Is.LessThan(24L * (stats.Levels[0].SampleSize + 1) + (4L * stats.Levels[0].SampleSize) + 100_000), "a text with fewer sample positions than K_s needs no larger buffers");
    }

    [Test]
    public void TheBytesOfTheRankFile_AreLittleEndianUInts_OneBySlot()
    {
        SaCase c = SaCase.Of("abracadabraabracadabra");
        DifferenceCover cover = DifferenceCover.Create(7);
        RunAndCompare(_dir, c, cover, SaTiny.Dcs(Scratch, 0), "bytes", out string rank);
        byte[] bytes = File.ReadAllBytes(rank);
        uint[] expected = DcsOracle.ExpectedRanks(c.Symbols, cover, out _);
        Assert.That(bytes.AsSpan().SequenceEqual(MemoryMarshal.AsBytes(expected.AsSpan())), Is.True);
    }
}
