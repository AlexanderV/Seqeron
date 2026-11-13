using SuffixTree.EmProto.Experiments;
using SuffixTree.EmProto.Pipeline;

namespace SuffixTree.EmProto.Tests.WholePipeline;

/// <summary>
/// The scale-down rule of E6 (ADR-004 section 5): the values that E3b (K_s, N_base, R, K_sort, run spill), E3 (K = n / 3000) and E4 (R, B2 workers) derive from <c>--scale-n 3000000000</c>, applied to every stage of the pipeline.
/// The expected numbers are computed by hand from the reference values (<c>1 Mi, 16 Mi, 4 Mi, 1 Mi, 64 KiB, 10^6, 8,000,000</c> times n / 3x10^9), not by the code under test.
/// </summary>
[TestFixture]
public class ScaleRulesTests
{
    private static EndToEndOptions Opts(long scaleN, EndToEndTuning? tuning = null, int threads = 4, long cap = 160L << 20) => new()
    {
        InputPath = "in",
        OutPrefix = "out",
        ScratchDir = "scratch",
        ScaleN = scaleN,
        Threads = threads,
        ArenaCapBytes = cap,
        Tuning = tuning ?? new EndToEndTuning(),
    };

    [TestCase(100_000_000L, 34953L, 559241L, 139810L, 34953L, 2185, 33333L, 266667L)]
    [TestCase(10_000_000L, 3495L, 55924L, 13981L, 3495L, 218, 3333L, 26667L)]
    [TestCase(50_000_000L, 17476L, 279620L, 69905L, 17476L, 1092, 16666L, 133334L)]
    [TestCase(200_000_000L, 69905L, 1118481L, 279620L, 69905L, 4369, 66666L, 533334L)]
    [TestCase(3_000_000_000L, 1_048_576L, 16_777_216L, 4_194_304L, 1_048_576L, 65_536, 1_000_000L, 8_000_000L)]
    public void ScaledBy_ScaleN_3e9_GivesThe3e9ReferenceTimesNOver3e9(long n, long ks, long nbase, long r1, long ksort, int runSpill, long k, long r4)
    {
        EndToEndParameters p = ScaleRules.Resolve(n, Opts(3_000_000_000));
        Assert.That(p.Scaled, Is.True);
        Assert.That(p.Factor, Is.EqualTo(n / 3e9).Within(1e-15));
        Assert.That(p.DcsBucketCapacity, Is.EqualTo(ks), "PH1 K_s");
        Assert.That(p.DcsBaseCase, Is.EqualTo(nbase), "PH1 N_base");
        Assert.That(p.DcsChunkSlots, Is.EqualTo(r1), "PH1 R");
        Assert.That(p.DcsSortCapacity, Is.EqualTo(ksort), "PH1 K_sort");
        Assert.That(p.DcsRunSpillBytes, Is.EqualTo(runSpill), "PH1 run spill");
        Assert.That(p.BucketCapacity, Is.EqualTo(k), "PH3 K = n / 3000 (E3), the reference K is 10^6");
        Assert.That(p.LcpChunk, Is.EqualTo(r4), "PH4 R = ceil(8,000,000 f)");
        Assert.That(p.LcpWorkers, Is.EqualTo(2), "PH4 B2 tasks from the arena cap at the reference geometry (E4: 2 of 4 threads)");
        Assert.That((n + p.LcpChunk - 1) / p.LcpChunk, Is.LessThanOrEqualTo(375).And.GreaterThanOrEqualTo(375 - 1).Or.EqualTo(375), "P = n / R chunks, 375 at 3x10^9");
    }

    [Test]
    public void At100Mbp_TheChunkCountIs375_AndTheBucketCountIsThatOf3e9()
    {
        EndToEndParameters p = ScaleRules.Resolve(100_000_000, Opts(3_000_000_000));
        Assert.That((100_000_000 + p.LcpChunk - 1) / p.LcpChunk, Is.EqualTo(375));
        Assert.That(100_000_000 / p.BucketCapacity, Is.EqualTo(3000).Within(1), "n / K = 3000 buckets, as at 3x10^9 with K = 10^6");
    }

    [Test]
    public void Unscaled_UsesTheReferenceValues_WhateverN()
    {
        foreach (long n in new[] { 1000L, 100_000_000L })
        {
            EndToEndParameters p = ScaleRules.Resolve(n, Opts(0));
            Assert.That(p.Scaled, Is.False);
            Assert.That(p.Factor, Is.EqualTo(1.0));
            Assert.That((p.DcsBucketCapacity, p.DcsBaseCase, p.DcsChunkSlots, p.DcsSortCapacity, p.DcsRunSpillBytes), Is.EqualTo((1L << 20, 16L << 20, 4L << 20, 1L << 20, 64 * 1024)));
            Assert.That(p.BucketCapacity, Is.EqualTo(1L << 20));
            Assert.That(p.LcpChunk, Is.EqualTo(8_000_000));
            Assert.That(p.FanOut, Is.EqualTo(64));
            Assert.That(p.DcsFanOut, Is.EqualTo(1024));
            Assert.That(p.BlockEntries, Is.EqualTo(256 * 1024));
            Assert.That((p.ChildWindowBytes, p.ChildSlots, p.WindowBytes), Is.EqualTo((1 << 20, 8, 1 << 20)));
            Assert.That(p.SpillFlushBytes, Is.EqualTo(16L << 20));
        }
    }

    [Test]
    public void TinyN_HitsTheFloors()
    {
        EndToEndParameters q = ScaleRules.Resolve(100, Opts(3_000_000_000));
        Assert.That((q.DcsBucketCapacity, q.DcsBaseCase, q.DcsChunkSlots, q.DcsSortCapacity, q.DcsRunSpillBytes), Is.EqualTo((3L, 4L, 16L, 3L, 0)), "N_base has a floor of 4");
        EndToEndParameters p = ScaleRules.Resolve(1000, Opts(3_000_000_000));
        Assert.That((p.DcsBucketCapacity, p.DcsBaseCase, p.DcsChunkSlots, p.DcsSortCapacity, p.DcsRunSpillBytes), Is.EqualTo((3L, 6L, 16L, 3L, 0)));
        Assert.That(p.BucketCapacity, Is.EqualTo(2), "K has a floor of 2 (a bucket of one position cannot be split)");
        Assert.That(p.LcpChunk, Is.EqualTo(3), "R has a floor of 3");
    }

    [Test]
    public void AnExplicitValue_Wins_AndTheOthersStillScale()
    {
        EndToEndParameters p = ScaleRules.Resolve(100_000_000, Opts(3_000_000_000, new EndToEndTuning { BucketCapacity = 16_666, DcsBucketCapacity = 17_476, LcpChunk = 133_334, LcpWorkers = 3, BlockEntries = 5 }));
        Assert.That(p.BucketCapacity, Is.EqualTo(16_666));
        Assert.That(p.DcsBucketCapacity, Is.EqualTo(17_476));
        Assert.That(p.LcpChunk, Is.EqualTo(133_334));
        Assert.That(p.LcpWorkers, Is.EqualTo(3));
        Assert.That(p.BlockEntries, Is.EqualTo(5));
        Assert.That(p.DcsBaseCase, Is.EqualTo(559_241), "not overridden: still scaled");
        Assert.That(p.DcsChunkSlots, Is.EqualTo(139_810));
    }

    [Test]
    public void TheWorkersOfPh4_FollowTheArenaCapAndTheThreads_AtTheReferenceGeometry()
    {
        Assert.That(ScaleRules.Resolve(1_000_000, Opts(3_000_000_000, threads: 4)).LcpWorkers, Is.EqualTo(2));
        Assert.That(ScaleRules.Resolve(1_000_000, Opts(3_000_000_000, threads: 1)).LcpWorkers, Is.EqualTo(1), "never more than the threads");
        Assert.That(ScaleRules.Resolve(1_000_000, Opts(3_000_000_000, threads: 2)).LcpWorkers, Is.EqualTo(2));
        Assert.That(ScaleRules.Resolve(1_000_000, Opts(3_000_000_000, cap: 1L << 20)).LcpWorkers, Is.EqualTo(1), "never fewer than one");
        Assert.That(ScaleRules.Resolve(1_000_000, Opts(3_000_000_000, cap: 4L << 30, threads: 8)).LcpWorkers, Is.EqualTo(8));
        // the number does not depend on n (it is the reference geometry's)
        Assert.That(ScaleRules.Resolve(5_000, Opts(3_000_000_000)).LcpWorkers, Is.EqualTo(ScaleRules.Resolve(2_000_000_000, Opts(3_000_000_000)).LcpWorkers));
    }

    [Test]
    public void TheConstants_AreThoseOfTheStageCommands()
    {
        Assert.That((E3b.ReferenceN, E3b.KsRef, E3b.BaseRef, E3b.ChunkRef, E3b.SortKRef), Is.EqualTo((3_000_000_000L, 1L << 20, 16L << 20, 4L << 20, 1L << 20)));
        Assert.That((E4.ReferenceN, E4.ChunkRef), Is.EqualTo((3_000_000_000L, 8_000_000L)));
        Assert.That(ScaleRules.BucketReference, Is.EqualTo(3_000_000_000L / 3000));
    }

    [Test]
    public void ANegativeScale_IsRefused()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => ScaleRules.Resolve(10, Opts(-1)));
        Assert.Throws<ArgumentOutOfRangeException>(() => ScaleRules.Resolve(-1, Opts(0)));
    }

    [Test]
    public void ToString_NamesTheValuesThatTheResultsPrint()
    {
        string s = ScaleRules.Resolve(100_000_000, Opts(3_000_000_000)).ToString();
        Assert.That(s, Does.Contain("K_s=34953").And.Contain("N_base=559241").And.Contain("K=33333").And.Contain("R=266667").And.Contain("P=375").And.Contain("B2 workers=2"));
    }
}
