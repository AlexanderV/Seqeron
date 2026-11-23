using System.Diagnostics;
using System.Globalization;
using System.Runtime.InteropServices;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The gaps that the orchestrator found in the audit of 1b.7d (ADR-004 section 5a, "Orchestrator verification (2026-11-19, 1b.7d"): the target of <see cref="Dist.Scatter"/> without a reservation option is created sparse (Y02),
/// the arena cap of PH3 is enforced on its own (Q52), the scratch root of the bucket sort is removed when it made it and only then (H08), the tie positions are counted (F06), and the accounting of the comparator
/// (the cost-only mutants of <c>SuffixOrder.CompareCore</c>). The pass-through of the arena guards to PH1 (Q76) changes nothing that a counter or an output shows on correct code and is argued in the audit, not tested.
/// </summary>
[TestFixture]
[HangLimit(30)]
public class SesaGapTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    /// <summary>Allocated bytes of a file (<c>stat %b</c> x 512), or null where there is no <c>stat</c>.</summary>
    private static long? AllocatedBytes(string path)
    {
        try
        {
            var info = new ProcessStartInfo("stat", ["-c", "%b", path]) { RedirectStandardOutput = true, RedirectStandardError = true, UseShellExecute = false };
            using Process p = Process.Start(info)!;
            string text = p.StandardOutput.ReadToEnd();
            p.WaitForExit();
            return p.ExitCode == 0 ? long.Parse(text.Trim(), CultureInfo.InvariantCulture) * 512 : null;
        }
        catch (Exception e) when (e is System.ComponentModel.Win32Exception or FormatException or InvalidOperationException)
        {
            return null;
        }
    }

    [Test]
    public void Y02_TheTargetOfScatterWithoutFiles_IsReserved_AsEveryOtherOutput_AndSparseOnlyOnRequest()
    {
        const int n = 4_000_000;
        string sa = _dir.File("sa.bin");
        uint[] reversed = new uint[n + 1];
        for (int i = 0; i <= n; i++) reversed[i] = (uint)(n - i);
        File.WriteAllBytes(sa, MemoryMarshal.AsBytes(reversed.AsSpan()).ToArray());
        string isa = _dir.File("isa.bin");
        // the budget lets pass 1 through and refuses the chunk array of pass 2 (16 MB), which is allocated right after the target file is created and before anything is written to it
        var options = new DistOptions { ScratchDir = _dir.Sub("scratch"), ChunkEntries = n + 1, IoBufferBytes = 4096, ReadBytes = 4096, Budget = new ArenaBudget(1 << 20) };
        Assert.Throws<ArenaBudgetExceededException>(() => Dist.BuildIsa(sa, n, isa, options));
        Assert.That(File.Exists(isa), Is.True, "the target is left (partial) after a failure: the documented contract");
        Assert.That(new FileInfo(isa).Length, Is.EqualTo(4L * (n + 1)));
        long? blocks = AllocatedBytes(isa);
        if (blocks is null) Assert.Ignore("no stat");
        if (OperatingSystem.IsLinux()) Assert.That(blocks, Is.GreaterThanOrEqualTo(4L * (n + 1)), "without Files the blocks are reserved at creation (fallocate), as for every other output of the build (work package 1b.9d: it was created sparse)");
        File.Delete(isa);
        var sparse = options with { Budget = new ArenaBudget(1 << 20), Files = new MappedFileOptions { Reservation = ReservationMode.Sparse } };
        Assert.Throws<ArenaBudgetExceededException>(() => Dist.BuildIsa(sa, n, isa, sparse));
        long? allocated = AllocatedBytes(isa);
        Assert.That(allocated, Is.LessThan(1 << 20), "created sparse on request: no block was reserved (a reservation would hold 16 MB)");
    }

    [Test]
    [HangLimit(130)]
    public void Q52_ThePh3Cap_IsEnforcedOnItsOwn_WhenPh1Fits()
    {
        string text = StringFamilies.RandomDna(30_000, 5);
        SaCase c = SaCase.Of(text);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        SaBuildOptions Options(string scratch, long cap) => SaTiny.Build(scratch, 1) with
        {
            ArenaCapBytes = cap,
            Tuning = SaTiny.Tuning(1) with { DcsBucketCapacity = 16, DcsBaseCase = 64, DcsSortCapacity = 16, BucketCapacity = 4000 },
        };
        SaBuildResult open = SaBuilder.Build(reader, _dir.File("a.sa"), Options(_dir.Sub("s0"), long.MaxValue));
        Assert.That(open.Ph23ArenaPeakBytes, Is.GreaterThan(open.Ph1ArenaPeakBytes), "the case needs PH3 to be the larger of the two");
        long cap = open.Ph23ArenaPeakBytes;
        Assert.DoesNotThrow(() => SaBuilder.Build(reader, _dir.File("b.sa"), Options(_dir.Sub("s1"), cap)), "exactly the peak fits");
        var e = Assert.Throws<ArenaBudgetExceededException>(() => SaBuilder.Build(reader, _dir.File("c.sa"), Options(_dir.Sub("s2"), cap - 1)), "one byte less is refused in PH3, which PH1 did not need");
        Assert.That(e!.Message, Is.Not.Empty);
        Assert.That(File.Exists(_dir.File("c.sa")), Is.False);
        Assert.That(open.Ph1ArenaPeakBytes, Is.LessThanOrEqualTo(cap - 1), "PH1 itself fitted under the cap that stopped PH3");
    }

    [Test]
    public void H08_TheBucketSortRemovesItsScratchRoot_OnlyWhenItMadeIt()
    {
        SaCase c = SaCase.Of(StringFamilies.RandomDna(400, 3));
        DifferenceCover cover = DifferenceCover.Create(13);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        string made = _dir.Sub("made-by-the-sort");
        BucketSorter.Run(reader, SaTiny.Sort(made, 2), new MemorySaSink(c.N + 1), ranks);
        Assert.That(Directory.Exists(made), Is.False, "a root that the sort made is removed when it is empty");
        string given = _dir.Sub("given");
        Directory.CreateDirectory(given);
        BucketSorter.Run(reader, SaTiny.Sort(given, 2), new MemorySaSink(c.N + 1), ranks);
        Assert.That(Directory.Exists(given), Is.True, "a root that existed before is not the sort's to remove");
        Assert.That(Directory.GetFileSystemEntries(given), Is.Empty);
        string foreign = _dir.Sub("foreign");
        Directory.CreateDirectory(foreign);
        File.WriteAllText(Path.Combine(foreign, "other.txt"), "not ours");
        BucketSorter.Run(reader, SaTiny.Sort(foreign, 2), new MemorySaSink(c.N + 1), ranks);
        Assert.That(File.ReadAllText(Path.Combine(foreign, "other.txt")), Is.EqualTo("not ours"));
    }

    [Test]
    public void F06_TiePositions_AreCounted_WhereKeysTie_AndNotWhereTheyDoNot()
    {
        foreach ((string label, string text, bool ties) in new[] { ("a^n", new string('a', 2000), true), ("random protein", StringFamilies.RandomProtein(2000, 9), false), ("tandem repeat", string.Concat(Enumerable.Repeat("ACGTTGCAAC", 300)), true) })
        {
            SaCase c = SaCase.Of(text);
            string textPath = SaRun.WriteText(_dir, c, label.Replace(' ', '-').Replace('^', '-'));
            using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
            SaBuildResult r = SaBuilder.Build(reader, _dir.File(label.Replace(' ', '-').Replace('^', '-') + ".sa"), SaTiny.Build(_dir.Sub("s" + label.Length), 3));
            Assert.That(r.Sort.Work.TiePositions > 0, Is.EqualTo(ties), $"{label}: tie positions {r.Sort.Work.TiePositions}");
            Assert.That(r.Sort.Work.TiePositions, Is.LessThanOrEqualTo(r.Sort.Work.LeafPositions), $"{label}: ties are positions of leaf buckets");
        }
    }

    [Test]
    public void TheComparator_AccountsTheLceSymbolsTheCallsAndTheRankLookups_OfEveryKindOfPair()
    {
        AssertComparatorAccounting(_dir, 13);
        // a cover whose witness for equal residues is not the witness for residues one apart (the first such period from 16 on): the equal residues must take their own
        int other = Enumerable.Range(16, 100).First(p => DifferenceCover.Create(p).Witness[0] != DifferenceCover.Create(p).Witness[1]);
        using var second = new TempDir();
        AssertComparatorAccounting(second, other);
    }

    private static void AssertComparatorAccounting(TempDir dir, int period)
    {
        // a text of 130 symbols (8 bits: the key holds 8 symbols) that starts with a long run, so that equal prefixes are as long as the comparator may look
        string text = new string('a', 300) + string.Concat(Enumerable.Range(0x80, 128).Select(k => (char)k)) + "b";
        SaCase c = SaCase.Of(text);
        Assert.That(c.Bits, Is.EqualTo(8));
        DifferenceCover cover = DifferenceCover.Create(period);
        string textPath = SaRun.WriteText(dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        Assert.That(reader.KeySymbols, Is.EqualTo(8));
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(dir, c, cover), cover, c.N);
        int x = cover.Period;
        int sawBeyondKey = 0, sawZero = 0, sawEqualResidues = 0;
        for (int i = 0; i < 40; i++)
        {
            for (int j = 0; j < 40; j++)
            {
                if (i == j) continue;
                int ri = i % x, rj = j % x;
                int delta = (rj - ri + x) % x;
                int l = (cover.Witness[delta] - ri + x) % x;
                var work = new SortWork();
                new SuffixOrder(reader, ranks, work).Compare((uint)i, (uint)j);
                string pair = $"({i}, {j}), l = {l}";
                Assert.That(work.LceSymbols, Is.EqualTo(l), $"{pair}: symbols compared");
                Assert.That(work.LceCalls, Is.EqualTo(l > 0 ? 1 : 0), $"{pair}: the LCE is called only when there is something to compare");
                Assert.That(work.RankLookups, Is.EqualTo(2), $"{pair}: two sample ranks decide the tie");
                var tied = new SortWork();
                new SuffixOrder(reader, ranks, tied).CompareTied((uint)i, (uint)j);
                int expectedTied = Math.Max(0, l - reader.KeySymbols);
                Assert.That(tied.LceSymbols, Is.EqualTo(expectedTied), $"{pair}: after a tie of the key the first {reader.KeySymbols} symbols are not compared again");
                Assert.That(tied.LceCalls, Is.EqualTo(expectedTied > 0 ? 1 : 0), $"{pair}: calls after a tie");
                Assert.That(tied.RankLookups, Is.EqualTo(2));
                if (l == 0) sawZero++;
                if (l > reader.KeySymbols) sawBeyondKey++;
                if (delta == 0) sawEqualResidues++;
            }
        }
        Assert.That(new[] { sawZero, sawBeyondKey, sawEqualResidues }, Is.All.GreaterThan(0), "the pairs cover a zero extension, an extension beyond the key and equal residues");

        // the pairs at the end of the text, where a suffix is as short as the key or shorter: the skip of the key applies only to a suffix that holds all of it
        int n = c.N;
        int sawShortAsKey = 0;
        for (int i = n - 14; i < n; i++)
        {
            for (int j = n - 14; j < n; j++)
            {
                if (i == j) continue;
                foreach (bool tied in new[] { false, true })
                {
                    int ri = i % x, rj = j % x;
                    int l = (cover.Witness[(rj - ri + x) % x] - ri + x) % x;
                    long maxLen = n - Math.Max(i, j);
                    long skip = tied && reader.KeySymbols <= maxLen ? reader.KeySymbols : 0;
                    long lim = Math.Min(l, maxLen);
                    long skipped = Math.Min(skip, lim);
                    int lcp = 0;
                    while (lcp < lim - skipped && c.Symbols[i + skipped + lcp] == c.Symbols[j + skipped + lcp]) lcp++;
                    long matched = skipped + (lim > skipped ? lcp : 0);
                    var work = new SortWork();
                    var order = new SuffixOrder(reader, ranks, work);
                    if (tied) order.CompareTied((uint)i, (uint)j);
                    else order.Compare((uint)i, (uint)j);
                    string pair = $"({i}, {j}), tied {tied}, l {l}, {maxLen} symbols in the shorter suffix";
                    Assert.That(work.LceCalls, Is.EqualTo(lim > skipped ? 1 : 0), $"{pair}: calls");
                    Assert.That(work.LceSymbols, Is.EqualTo(lim > skipped ? lcp : 0), $"{pair}: symbols");
                    Assert.That(work.RankLookups, Is.EqualTo(matched >= lim && l <= maxLen ? 2 : 0), $"{pair}: rank lookups");
                    if (maxLen == reader.KeySymbols) sawShortAsKey++;
                }
            }
        }
        Assert.That(sawShortAsKey, Is.GreaterThan(0), "a suffix with exactly as many symbols as the key was compared");
    }
}
