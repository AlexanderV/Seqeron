using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.External;

/// <summary>
/// The shape of the work of PH1 and of PH2 / PH3 pinned to numbers (work package 1b.7d, found by the mutation audit: a change of a seed, a fan-out, a sample, a buffer or a count of the partition leaves the suffix array
/// untouched, which is the design, but changes the tree of buckets, the scans and the counted I/O, and the differential tests against the prototype that would see it live in another project). The figures are those of
/// production with one thread and the tiny parameters of <see cref="SaTiny"/>; the differential tests of <c>EmProto.Tests</c> show that they are the prototype's (same seeds, same rule, same counted I/O). A line that changes
/// is a change of the algorithm, not of its output, and has to be understood before the figure is moved.
/// </summary>
[TestFixture]
[HangLimit(30)]
public class SaPinnedTests
{
    private TempDir _dir = null!;

    [SetUp]
    public void SetUp() => _dir = new TempDir();

    [TearDown]
    public void TearDown() => _dir.Dispose();

    private static string SortLine(SortStats s) =>
        $"buckets {s.Buckets} leaves {s.Leaves} partitions {s.Partitions} depth {s.MaxDepth} maxleaf {s.MaxLeaf} retries {s.Retries} compares {s.Work.Compares} tied {s.Work.TiedCompares} lce {s.Work.LceCalls}/{s.Work.LceSymbols} ranks {s.Work.RankLookups} tiepos {s.Work.TiePositions} leafpos {s.Work.LeafPositions} scanpos {s.Work.ScanPositions}";

    private static string IoLine(IoSnapshot io) => $"{io.ReadBytes}/{io.WriteBytes}/{io.ReadCalls}/{io.WriteCalls}/{io.MappedWriteBytes}/{io.MappedWindows}/{io.MappedMaxWindowBytes}";

    private static string BuildLine(SaBuildResult r) =>
        $"cover {r.CoverKind}/{r.CoverResidues}/{r.Period} rank {r.RankFileBytes} levels {r.RecursionLevels}\n"
        + string.Join('\n', r.Dcs.Levels.Select(l => $"L{l.Level} n {l.TextLength} s {l.SampleSize} d {l.DistinctNames} u {l.NonUniqueRecords} base {l.BaseCase} leaves {l.Leaves} split {l.SplitScans} scans {l.Scans} files {l.FileRuns}")) + "\n"
        + SortLine(r.Sort) + $"\nph1 {IoLine(r.Ph1Io)}\nph23 {IoLine(r.Ph23Io)}\nscratch {r.Ph1ScratchPeakBytes}/{r.Ph23ScratchPeakBytes} arena {r.Ph1ArenaPeakBytes}/{r.Ph23ArenaPeakBytes}";

    private static void Pin(string label, string actual, string expected)
    {
        TestContext.Out.WriteLine($"[{label}]\n{actual}");
        Assert.That(actual, Is.EqualTo(expected), $"{label}: the shape of the work changed:\n{actual}");
    }

    [TestCase(0)]
    [TestCase(1)]
    [TestCase(2)]
    [TestCase(3)]
    [TestCase(4)]
    [TestCase(5)]
    [HangLimit(50)]
    public void TheBuild_HasThePinnedShape_AndCounts(int which)
    {
        (string name, string text, int variant) = which switch
        {
            0 => ("genome-like 20000", SaRun.GenomeLike(20_000, 5), 3),
            1 => ("random dna 12000", StringFamilies.RandomDna(12_000, 9), 10),
            2 => ("a^n 2500", new string('a', 2500), 7),
            3 => ("fibonacci 5000", StringFamilies.Fibonacci(5000), 8),
            4 => ("genome-like 20000, N_base 3000", SaRun.GenomeLike(20_000, 5), 3),
            _ => ("random dna 9000, reservoir 4 fan-out 5", StringFamilies.RandomDna(9000, 4), 11), // a sample of four for five children: equal splitter indices
        };
        SaCase c = SaCase.Of(text);
        SaBuildOptions options = SaTiny.Build(Scratch, variant, 1) with { ArenaGuards = null, Counters = new IoCounters() };
        if (which == 4) options = options with { Tuning = options.Tuning with { DcsBaseCase = 3000 } }; // the recursion ends in the base case at level 2
        SaRun.Outcome o = SaRun.Run(_dir, c, options);
        SaRun.Compare(c, o, name);
        Pin(name, BuildLine(o.Result), Pins.Build[which]);
    }

    [TestCase(0)]
    [TestCase(1)]
    [TestCase(2)]
    public void TheSort_HasThePinnedShape_AndCounts(int which)
    {
        (string name, string text, int variant) = which switch
        {
            0 => ("genome-like 9000", SaRun.GenomeLike(9000, 2), 5),
            1 => ("a^n 1500", new string('a', 1500), 2),
            _ => ("tandem 4000", StringFamilies.Periodic(4000, 37, 5), 9),
        };
        SaCase c = SaCase.Of(text);
        DifferenceCover cover = DifferenceCover.Create(SaTiny.Periods[variant % SaTiny.Periods.Length]);
        string textPath = SaRun.WriteText(_dir, c);
        using MappedTextReader reader = MappedTextReader.Open(textPath, c.N, c.Bits);
        using RankTable ranks = RankTable.Open(SaInputs.WriteRankFile(_dir, c, cover), cover, c.N);
        var counters = new IoCounters();
        var sink = new MemorySaSink(c.N + 1);
        SortStats s = BucketSorter.Run(reader, SaTiny.Sort(Scratch, variant, 1) with { Counters = counters, ArenaGuards = null }, sink, ranks);
        Assert.That(sink.Data, Is.EqualTo(c.Sa), name);
        Pin(name, SortLine(s) + $"\nio {IoLine(counters.Snapshot())}", Pins.Sort[which]);
    }

    private string Scratch => _dir.Sub("scratch");
}

/// <summary>The pinned figures (see <see cref="SaPinnedTests"/>).</summary>
internal static class Pins
{
    internal static readonly string[] Build =
    [
        "cover Singer/4/13 rank 24624 levels 2\nL0 n 20000 s 6155 d 4749 u 1956 base False leaves 2697 split 815 scans 4326 files 0\nL1 n 6159 s 1895 d 1668 u 275 base False leaves 814 split 261 scans 1335 files 0\nL2 n 1899 s 585 d 585 u 0 base False leaves 271 split 88 scans 446 files 0\nbuckets 8987 leaves 4494 partitions 4493 depth 17 maxleaf 6 retries 0 compares 120291 tied 35180 lce 78491/357885 ranks 158308 tiepos 10325 leafpos 20000 scanpos 259416\nph1 161104/128832/10602/2976/66788/597/112\nph23 909380/850532/16551/14280/80004/715/112\nscratch 81560/145176 arena 4740/144",
        "cover Singer/8/57 rank 6752 levels 0\nL0 n 12000 s 1685 d 1685 u 0 base False leaves 234 split 84 scans 401 files 23\nbuckets 1453 leaves 969 partitions 484 depth 14 maxleaf 29 retries 0 compares 2420 tied 968 lce 1430/4994 ranks 206 tiepos 0 leafpos 12000 scanpos 97396\nph1 2184/2184/273/23/6740/40/171\nph23 374828/369032/2236/6180/48004/281/171\nscratch 2184/88088 arena 312/696",
        "cover Pruned/3/5 rank 6012 levels 4\nL0 n 2500 s 1500 d 4 u 1497 base False leaves 244 split 79 scans 401 files 0\nL1 n 1503 s 903 d 13 u 891 base False leaves 137 split 39 scans 214 files 0\nL2 n 906 s 544 d 31 u 515 base False leaves 93 split 29 scans 150 files 0\nL3 n 547 s 329 d 87 u 243 base False leaves 56 split 20 scans 95 files 0\nL4 n 332 s 200 d 200 u 0 base False leaves 26 split 7 scans 39 files 0\nbuckets 277 leaves 185 partitions 92 depth 5 maxleaf 26 retries 0 compares 29846 tied 26534 lce 3441/8696 ranks 59238 tiepos 2500 leafpos 2500 scanpos 11881\nph1 38224/38160/533/537/27096/457/60\nph23 44468/41192/954/710/10004/167/60\nscratch 20228/20228 arena 2564/624",
        "cover Pruned/7/25 rank 5628 levels 2\nL0 n 5000 s 1401 d 82 u 1384 base False leaves 110 split 109 scans 327 files 110\nL1 n 1408 s 396 d 331 u 126 base False leaves 33 split 32 scans 96 files 33\nL2 n 403 s 114 d 114 u 0 base False leaves 9 split 8 scans 24 files 9\nbuckets 316 leaves 231 partitions 85 depth 6 maxleaf 43 retries 0 compares 34607 tied 28987 lce 5405/58998 ranks 67366 tiepos 4931 leafpos 5000 scanpos 20807\nph1 36860/36796/2182/423/14932/157/97\nph23 84396/79420/1478/1368/20004/207/97\nscratch 16904/34288 arena 2240/1032",
        "cover Singer/4/13 rank 24624 levels 1\nL0 n 20000 s 6155 d 4749 u 1956 base False leaves 2697 split 815 scans 4326 files 0\nL1 n 6159 s 1895 d 1668 u 275 base True leaves 814 split 261 scans 1335 files 0\nbuckets 8987 leaves 4494 partitions 4493 depth 17 maxleaf 6 retries 0 compares 120291 tied 35180 lce 78491/357885 ranks 158308 tiepos 10325 leafpos 20000 scanpos 259416\nph1 123152/98496/8070/2260/56848/508/112\nph23 909380/850532/16551/14280/80004/715/112\nscratch 81560/145176 arena 15416/144",
        "cover Singer/3/7 rank 15432 levels 0\nL0 n 9000 s 3857 d 3857 u 0 base False leaves 448 split 169 scans 785 files 69\nbuckets 581 leaves 421 partitions 160 depth 8 maxleaf 46 retries 0 compares 3097 tied 420 lce 2284/5425 ranks 3266 tiepos 0 leafpos 9000 scanpos 44801\nph1 9752/9752/1219/69/15428/1929/8\nph23 176688/172080/1569/2901/36004/4501/8\nscratch 9752/57940 arena 708/1104",
    ];

    internal static readonly string[] Sort =
    [
        "buckets 728 leaves 530 partitions 198 depth 9 maxleaf 40 retries 0 compares 18892 tied 18062 lce 774/2925 ranks 35682 tiepos 5497 leafpos 9000 scanpos 47208\nio 182708/179816/1187/3071/0/0/0",
        "buckets 85 leaves 64 partitions 21 depth 3 maxleaf 37 retries 0 compares 14727 tied 13341 lce 1798/7511 ranks 29256 tiepos 1500 leafpos 1500 scanpos 4500\nio 18436/17476/315/311/0/0/0",
        "buckets 935 leaves 468 partitions 467 depth 12 maxleaf 12 retries 0 compares 43697 tied 24748 lce 18391/266385 ranks 83894 tiepos 3957 leafpos 4000 scanpos 36581\nio 114740/107740/1925/942/0/0/0",
    ];
}
