using System.Text;
using SuffixTree.EmProto.Tests.WholePipeline;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;
using Esa = global::SuffixTree.Esa;
using Prod = global::SuffixTree.Esa.External;
using ProtoDcs = SuffixTree.EmProto.Dcs;
using ProtoPipe = SuffixTree.EmProto.Pipeline;

namespace SuffixTree.EmProto.Tests.External;

/// <summary>Helpers of the differential tests of work package 1b.7d: the prototype <c>e6</c> makes the packed text and the reference suffix array, the production stage runs on the prototype's own packed file.</summary>
internal static class SaBuildDiff
{
    internal static ProtoPipe.EndToEndResult RunProto(DiffDir dir, string inputPath, string name, int threads = 4, ProtoPipe.EndToEndTuning? tuning = null, bool exact = false, int period = 1024, long scaleN = 0)
        => ProtoPipe.EndToEnd.Run(new ProtoPipe.EndToEndOptions
        {
            InputPath = inputPath,
            OutPrefix = dir.File(name),
            ScratchDir = dir.File(name + "-scratch"),
            Threads = threads,
            SyncOutputs = false,
            Tuning = tuning ?? new ProtoPipe.EndToEndTuning(),
            Period = period,
            ExactPeriod = exact,
            ScaleN = scaleN,
        });

    /// <summary>The same sizes for the production stage: the tuning of the prototype's pipeline, field by field.</summary>
    internal static Prod.SaBuildTuning ToProd(ProtoPipe.EndToEndTuning t) => new()
    {
        DcsBucketCapacity = t.DcsBucketCapacity,
        DcsBaseCase = t.DcsBaseCase,
        DcsChunkSlots = t.DcsChunkSlots,
        DcsSortCapacity = t.DcsSortCapacity,
        DcsRunSpillBytes = t.DcsRunSpillBytes,
        DcsFanOut = t.DcsFanOut,
        DcsReservoir = t.DcsReservoir,
        DcsIoBytes = t.DcsIoBytes,
        DcsMergeBudget = t.DcsMergeBudget,
        BucketCapacity = t.BucketCapacity,
        FanOut = t.FanOut,
        Oversample = t.Oversample,
        SpillBytes = t.SpillBytes,
        WindowBytes = t.WindowBytes,
    };

    /// <summary>The production stage on the packed text of a prototype run; the suffix array goes to <paramref name="outName"/> next to it.</summary>
    internal static (Prod.SaBuildResult Result, string SaPath) RunProd(DiffDir dir, ProtoPipe.EndToEndResult proto, string outName, Prod.SaBuildOptions options, Action<Prod.MappedTextReader>? configure = null)
    {
        string sa = dir.File(outName + ".sa");
        using Prod.MappedTextReader text = Prod.MappedTextReader.Open(proto.Files.Packed, proto.Prepare.N, proto.Prepare.Bits);
        configure?.Invoke(text);
        return (Prod.SaBuilder.Build(text, sa, options), sa);
    }

    internal static void AssertSameSa(ProtoPipe.EndToEndResult proto, string prodSa, string label)
        => Assert.That(DiffInputs.SameBytes(proto.Files.Sa, prodSa), Is.True, $"{label}: the suffix array differs from the prototype's");

    /// <summary>The counted I/O of a phase against the prototype's, within 2 % (the explicit reads and writes and the mapped writes, each).</summary>
    internal static void AssertSameIo(ProtoPipe.PhaseRecord phase, Prod.IoSnapshot io, string label, double tolerance = 0.02)
    {
        AssertClose(io.ReadBytes, phase.ExplicitReadBytes, tolerance, $"{label}: explicit reads");
        AssertClose(io.WriteBytes, phase.ExplicitWriteBytes, tolerance, $"{label}: explicit writes");
        AssertClose(io.MappedWriteBytes, phase.MappedWriteBytes, tolerance, $"{label}: mapped writes");
        TestContext.Out.WriteLine($"{label}: read {io.ReadBytes} / {phase.ExplicitReadBytes}, write {io.WriteBytes} / {phase.ExplicitWriteBytes}, mapped {io.MappedWriteBytes} / {phase.MappedWriteBytes} (production / prototype)");
    }

    internal static void AssertClose(long production, long prototype, double tolerance, string label)
    {
        double allowed = Math.Max(4096, prototype * tolerance);
        Assert.That(Math.Abs(production - prototype), Is.LessThanOrEqualTo(allowed), $"{label}: production {production}, prototype {prototype}");
    }
}

/// <summary>
/// G2 of work package 1b.7d: the production stage (PH1, PH2 and PH3) on the packed text of the prototype <c>e6</c> writes the suffix array byte for byte as <c>e6</c> did, its rank file is the prototype's byte for byte, the covers and the
/// bucket trees equal, and the counted I/O is the same within 2 %. The prototype code is unchanged.
/// </summary>
[TestFixture]
public class SaBuildDifferentialTests
{
    private static void Differential(DiffDir dir, string input, string label, long scaleN = 0, int threads = 4)
    {
        ProtoPipe.EndToEndResult proto = SaBuildDiff.RunProto(dir, input, "proto", threads, scaleN: scaleN);
        (Prod.SaBuildResult result, string sa) = SaBuildDiff.RunProd(dir, proto, "prod", new Prod.SaBuildOptions { ScratchDir = dir.File("prod-scratch"), Threads = threads, ScaleN = scaleN });
        SaBuildDiff.AssertSameSa(proto, sa, label);
        Assert.That(Directory.Exists(dir.File("prod-scratch")), Is.False, $"{label}: scratch");
        SaBuildDiff.AssertSameIo(proto.Phase(ProtoPipe.EndToEnd.PhaseDcs), result.Ph1Io, $"{label} PH1");
        SaBuildDiff.AssertSameIo(proto.Phase(ProtoPipe.EndToEnd.PhaseSa), result.Ph23Io, $"{label} PH2/PH3");
        Assert.That(result.Dcs.MaxLevel, Is.EqualTo(proto.Dcs.MaxLevel), $"{label}: levels of the recursion");
        Assert.That(result.Period, Is.EqualTo(proto.CoverPeriod), $"{label}: the period");
        Assert.That(result.CoverResidues, Is.EqualTo(proto.CoverResidues));
        Assert.That(result.Sort.Leaves, Is.EqualTo(proto.Build.Leaves), $"{label}: leaf buckets");
        Assert.That(result.Sort.Partitions, Is.EqualTo(proto.Build.Partitions), $"{label}: partitions");
        Assert.That(result.Sort.Work.LeafPositions, Is.EqualTo(proto.Build.Work.LeafPositions));
        Assert.That(result.Sort.Work.Compares, Is.EqualTo(proto.Build.Work.Compares), $"{label}: the comparator is called as often");
        Assert.That(result.Sort.Work.RankLookups, Is.EqualTo(proto.Build.Work.RankLookups));
        Assert.That(result.RankFileBytes, Is.EqualTo(4L * new ProtoDcs.SampleGeometry(proto.Prepare.N, ProtoDcs.DifferenceCover.Create(proto.CoverPeriod)).Slots));
    }

    [Test]
    public void TheBenchmarkFiles_AreByteIdentical_ToTheE6Outputs()
    {
        foreach (string name in new[] { "prot_1m", "dna_10m", "rep_10m" })
        {
            string input = LcpChildDiff.DataFile(name);
            using var dir = new DiffDir();
            Differential(dir, input, name);
        }
    }

    [Test]
    public void TheFamilies_AtTenToTheSixth_AreByteIdentical_ToTheE6Outputs()
    {
        foreach ((string name, string text) in DiffInputs.Families(1_000_000))
        {
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            Differential(dir, input, name);
        }
    }

    [Test]
    public void RandomDnaAndRepeatRichAndNRuns_AtTenToTheSixth_WithTheScaledGeometry_AreByteIdentical()
    {
        foreach ((string name, string text) in new[] { ("dna", StringFamilies.RandomDna(1_000_000, 1)), ("repeat-rich", RepeatRich(1_000_000)), ("n-runs", NRuns(1_000_000)), ("genome-like", GenomeLike(1_000_000)) })
        {
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            Differential(dir, input, name, scaleN: 3_000_000_000);
        }
    }

    internal static string RepeatRich(int n)
    {
        char[] t = StringFamilies.RandomDna(n, 19).ToCharArray();
        var rng = new SplitMix64(1019);
        for (int copy = 0; copy < 6; copy++)
        {
            int len = n / 12;
            int from = rng.NextInt(n - len), to = rng.NextInt(n - len);
            Array.Copy(t, from, t, to, len);
            for (int s = 0; s < len / 50; s++) t[to + rng.NextInt(len)] = "ACGT"[rng.NextInt(4)];
        }
        return new string(t);
    }

    internal static string NRuns(int n)
    {
        char[] t = StringFamilies.RandomDna(n, 20).ToCharArray();
        var rng = new SplitMix64(1020);
        for (int run = 0; run < 3; run++) Array.Fill(t, 'N', rng.NextInt(n - (n / 15)), n / 15);
        return new string(t);
    }

    internal static string GenomeLike(int n)
    {
        char[] t = StringFamilies.RandomDna(n, 21).ToCharArray();
        var rng = new SplitMix64(3021);
        for (int sat = 0; sat < 5; sat++)
        {
            int unit = 3 + rng.NextInt(12), len = n / 20;
            int at = rng.NextInt(n - len);
            for (int k = 0; k < len; k++) t[at + k] = t[at + (k % unit)];
        }
        for (int copy = 0; copy < 4; copy++)
        {
            int len = n / 10;
            int from = rng.NextInt(n - len), to = rng.NextInt(n - len);
            Array.Copy(t, from, t, to, len);
        }
        for (int run = 0; run < 2; run++) Array.Fill(t, 'N', rng.NextInt(n - (n / 30)), n / 30);
        return new string(t);
    }

    [Test]
    public void SmallTexts_WithRotatingTinyParameters_AreByteIdentical_ToTheE6Outputs()
    {
        var rng = new SplitMix64(1_818);
        string[] letters = ["a", "ab", "acgt", "abcdefgh"];
        for (int v = 0; v < 120; v++)
        {
            int n = 1 + rng.NextInt(v % 4 == 0 ? 600 : 90);
            string alphabet = letters[rng.NextInt(letters.Length)];
            var sb = new StringBuilder();
            for (int i = 0; i < n; i++) sb.Append(alphabet[rng.NextInt(alphabet.Length)]);
            string text = v % 5 == 0 ? new string('a', n) : v % 5 == 1 ? StringFamilies.Fibonacci(n) : sb.ToString();
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            ProtoPipe.EndToEndTuning tuning = PipelineChecks.Tiny(v);
            ProtoPipe.EndToEndResult proto = SaBuildDiff.RunProto(dir, input, "proto", 1 + (v % 3), tuning, exact: true, period: PipelineChecks.Period(v));
            (Prod.SaBuildResult result, string sa) = SaBuildDiff.RunProd(dir, proto, "prod", new Prod.SaBuildOptions
            {
                ScratchDir = dir.File("prod-scratch"),
                Threads = 1 + (v % 4),
                Period = PipelineChecks.Period(v),
                ExactPeriod = true,
                Tuning = SaBuildDiff.ToProd(tuning),
                ArenaGuards = Esa.ArenaOptions.Guarded(0xCD),
            });
            SaBuildDiff.AssertSameSa(proto, sa, $"v={v} n={n}");
            Assert.That(result.Dcs.MaxLevel, Is.EqualTo(proto.Dcs.MaxLevel), $"v={v} n={n}: levels");
            Assert.That(result.Dcs.Levels.Select(l => (l.TextLength, l.SampleSize, l.DistinctNames, l.NonUniqueRecords, l.BaseCase)), Is.EqualTo(proto.Dcs.Levels.Select(l => (l.TextLength, l.SampleSize, l.DistinctNames, l.NonUniqueRecords, l.BaseCase))), $"v={v} n={n}: the levels");
        }
    }

    [Test]
    public unsafe void TheRankFile_IsTheE3bFile_ByteForByte_AndTheBucketTreesOfBothLevelsAreTheSame()
    {
        var rng = new SplitMix64(77);
        for (int v = 0; v < 60; v++)
        {
            int n = 20 + rng.NextInt(v % 3 == 0 ? 5000 : 400);
            string text = v % 4 == 0 ? new string('a', n) : v % 4 == 1 ? StringFamilies.Fibonacci(n) : v % 4 == 2 ? StringFamilies.ThueMorse(n) : StringFamilies.RandomDna(n, (ulong)v);
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            ProtoPipe.EndToEndTuning tuning = PipelineChecks.Tiny(v);
            ProtoPipe.EndToEndResult proto = SaBuildDiff.RunProto(dir, input, "proto", 1, tuning, exact: true, period: PipelineChecks.Period(v));
            int x = PipelineChecks.Period(v);
            using var mapped = new MappedFile(proto.Files.Packed, writable: false);
            using var packed = new PackedText(mapped.Pointer, proto.Prepare.N, proto.Prepare.Bits, mapped, ownsOwner: false);
            ProtoDcs.ExtDcsStats protoStats = ProtoDcs.ExternalDcs.BuildRankFile(packed, ProtoDcs.DifferenceCover.Create(x), dir.File("proto.rank"), new ProtoDcs.ExtDcsOptions
            {
                ScratchDir = dir.File("ps"),
                BucketCapacity = tuning.DcsBucketCapacity!.Value,
                BaseCase = tuning.DcsBaseCase!.Value,
                ChunkSlots = tuning.DcsChunkSlots!.Value,
                SortLeafCapacity = tuning.DcsSortCapacity!.Value,
                RunSpillBytes = tuning.DcsRunSpillBytes!.Value,
                FanOut = tuning.DcsFanOut!.Value,
                ReservoirSize = tuning.DcsReservoir!.Value,
                IoBufferBytes = tuning.DcsIoBytes!.Value,
                MergeBufferBudget = tuning.DcsMergeBudget!.Value,
                WindowBytes = tuning.WindowBytes!.Value,
                Seed = 1,
            });
            using Prod.MappedTextReader text2 = Prod.MappedTextReader.Open(proto.Files.Packed, proto.Prepare.N, proto.Prepare.Bits);
            Prod.ExtDcsStats prodStats = Prod.ExternalDcs.BuildRankFile(text2, Prod.DifferenceCover.Create(x), dir.File("prod.rank"), new Prod.ExtDcsOptions
            {
                ScratchDir = dir.File("qs"),
                BucketCapacity = tuning.DcsBucketCapacity!.Value,
                BaseCase = tuning.DcsBaseCase!.Value,
                ChunkSlots = tuning.DcsChunkSlots!.Value,
                SortLeafCapacity = tuning.DcsSortCapacity!.Value,
                RunSpillBytes = tuning.DcsRunSpillBytes!.Value,
                FanOut = tuning.DcsFanOut!.Value,
                ReservoirSize = tuning.DcsReservoir!.Value,
                IoBufferBytes = tuning.DcsIoBytes!.Value,
                MergeBufferBudget = tuning.DcsMergeBudget!.Value,
                WindowBytes = tuning.WindowBytes!.Value,
                Seed = 1,
            });
            Assert.That(DiffInputs.SameBytes(dir.File("proto.rank"), dir.File("prod.rank")), Is.True, $"v={v} n={n} X={x}: the rank file differs from the prototype's");
            Assert.That(prodStats.Levels.Select(l => (l.TextLength, l.SampleSize, l.DistinctNames, l.NonUniqueRecords, l.BaseCase, l.Leaves, l.SplitScans, l.Scans, l.FileRuns)),
                Is.EqualTo(protoStats.Levels.Select(l => (l.TextLength, l.SampleSize, l.DistinctNames, l.NonUniqueRecords, l.BaseCase, l.Leaves, l.SplitScans, l.Scans, l.FileRuns))), $"v={v} n={n} X={x}: the bucket trees");
        }
    }

    [Test]
    public void TheDifferenceCovers_AreThePrototypes_ForEveryPeriodUpTo2300()
    {
        for (int x = 2; x <= 2300; x++)
        {
            Assert.That(Prod.DifferenceCover.Create(x).Residues, Is.EqualTo(ProtoDcs.DifferenceCover.Create(x).Residues), $"X = {x}");
        }
        foreach (int target in new[] { 7, 100, 512, 1000, 1024, 2048, 3000 })
            Assert.That(Prod.DifferenceCover.NearestSingerPeriod(target), Is.EqualTo(ProtoDcs.DifferenceCover.NearestSingerPeriod(target)));
    }

    [Test]
    public void TheCountedIo_AndTheScratchPeaks_AreThePrototypes_InEachPhase()
    {
        foreach ((string name, string text) in new[] { ("dna", StringFamilies.RandomDna(600_000, 3)), ("a^n", new string('a', 300_000)), ("genome-like", GenomeLike(600_000)) })
        {
            using var dir = new DiffDir();
            string input = LcpChildDiff.TextFile(dir, "in", text);
            ProtoPipe.EndToEndResult proto = SaBuildDiff.RunProto(dir, input, "proto", scaleN: 3_000_000_000);
            (Prod.SaBuildResult result, _) = SaBuildDiff.RunProd(dir, proto, "prod", new Prod.SaBuildOptions { ScratchDir = dir.File("prod-scratch"), ScaleN = 3_000_000_000, Counters = new Prod.IoCounters() });
            SaBuildDiff.AssertSameIo(proto.Phase(ProtoPipe.EndToEnd.PhaseDcs), result.Ph1Io, $"{name} PH1");
            SaBuildDiff.AssertSameIo(proto.Phase(ProtoPipe.EndToEnd.PhaseSa), result.Ph23Io, $"{name} PH2/PH3");
            ProtoPipe.PhaseRecord ph1 = proto.Phase(ProtoPipe.EndToEnd.PhaseDcs);
            ProtoPipe.PhaseRecord ph3 = proto.Phase(ProtoPipe.EndToEnd.PhaseSa);
            TestContext.Out.WriteLine($"{name}: scratch peak PH1 {result.Ph1ScratchPeakBytes} / {ph1.ScratchPeakBytes}, PH3 {result.Ph23ScratchPeakBytes} / {ph3.ScratchPeakBytes}");
            Assert.That(result.Ph1ScratchPeakBytes, Is.LessThanOrEqualTo((long)(ph1.ScratchPeakBytes * 1.02) + 4096), $"{name}: PH1 scratch peak");
            Assert.That(result.Ph23ScratchPeakBytes, Is.LessThanOrEqualTo((long)(ph3.ScratchPeakBytes * 1.02) + 4096), $"{name}: PH3 scratch peak");
        }
    }

    [Test]
    public void TheThreadCounts_AndTheScalarPath_GiveTheSameFileAsThePrototype()
    {
        using var dir = new DiffDir();
        string input = LcpChildDiff.TextFile(dir, "in", GenomeLike(400_000));
        ProtoPipe.EndToEndResult proto = SaBuildDiff.RunProto(dir, input, "proto", 4, scaleN: 3_000_000_000);
        foreach ((int threads, bool vectors) in new[] { (1, true), (2, true), (4, false) })
        {
            (Prod.SaBuildResult _, string sa) = SaBuildDiff.RunProd(dir, proto, $"prod{threads}{vectors}", new Prod.SaBuildOptions { ScratchDir = dir.File($"prod-scratch{threads}{vectors}"), Threads = threads, ScaleN = 3_000_000_000 }, r => r.UseVectors = vectors);
            SaBuildDiff.AssertSameSa(proto, sa, $"{threads} threads, vectors {vectors}");
        }
    }
}
