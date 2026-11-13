using System.Globalization;
using SuffixTree.EmProto.Pipeline;
using SuffixTree.EmProto.Tests.Ph5;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.WholePipeline;

/// <summary>
/// Differential fuzz of the whole pipeline (ADR-004 section 5, "Oracles per stage"): scaled-down K = 3..50, F = 2..4, X = 7..64, tiny blocks, windows and chunks on texts of 10^3 to 10^5 symbols from the adversarial families,
/// 1 to 4 threads, against the in-RAM structures. A short deterministic run is part of the default tests; the long one is <c>[Explicit]</c>, category <c>Fuzz</c>, environment <c>EMPROTO_FUZZ_SEED</c> (default 1) and
/// <c>EMPROTO_FUZZ_SECONDS</c> (default 120).
/// </summary>
[TestFixture]
public class EndToEndFuzzTests : TempDirFixture
{
    private static EndToEndOptions RandomOptions(string dir, string text, SplitMix64 rng)
    {
        int v = rng.NextInt(1000);
        EndToEndOptions o = PipelineChecks.Options(dir, text, v, new[] { 1, 2, 3, 4 }[rng.NextInt(4)]);
        EndToEndTuning t = o.Tuning with
        {
            BucketCapacity = 3 + rng.NextInt(48),
            FanOut = 2 + rng.NextInt(3),
            Oversample = 3 + rng.NextInt(14),
            SpillBytes = 32 + rng.NextInt(200),
            DcsBucketCapacity = 3 + rng.NextInt(48),
            DcsBaseCase = rng.NextInt(3) == 0 ? 16_000_000 : 4 + rng.NextInt(197),
            DcsChunkSlots = 16 + rng.NextInt(100),
            DcsSortCapacity = 3 + rng.NextInt(48),
            DcsRunSpillBytes = 16 + rng.NextInt(300),
            DcsFanOut = 2 + rng.NextInt(4),
            DcsReservoir = 8 + rng.NextInt(200),
            DcsIoBytes = 16 + rng.NextInt(400),
            DcsMergeBudget = 8 + rng.NextInt(3000),
            LcpChunk = 3 + rng.NextInt(60) + (rng.NextInt(3) == 0 ? rng.NextInt(5000) : 0),
            LcpWorkers = rng.NextInt(3) == 0 ? null : 1 + rng.NextInt(4),
            LcpReadBytes = 16 + rng.NextInt(900),
            LcpTupleBytes = 48 + rng.NextInt(300),
            BlockEntries = 2 + rng.NextInt(49),
            ChildWindowBytes = 16 << rng.NextInt(9),
            ChildSlots = 1 + rng.NextInt(4),
            SpillFlushBytes = rng.NextInt(2) == 0 ? 0 : 50 + rng.NextInt(2000),
            WindowBytes = 16 + rng.NextInt(4000),
            PrepareBlockSymbols = 8 * (1 + rng.NextInt(40)),
        };
        return o with { Period = 7 + rng.NextInt(58), ExactPeriod = true, Tuning = t };
    }

    private (long Cases, long Multi, long Spill) Fuzz(ulong seed, Func<long, bool> goOn, int[] sizes)
    {
        var rng = new SplitMix64(seed * 0x9E3779B97F4A7C15UL);
        long cases = 0, multi = 0, spill = 0;
        while (goOn(cases))
        {
            string text;
            int pick = rng.NextInt(30);
            if (pick == 0)
            {
                string x = StringFamilies.RandomDna(300 + rng.NextInt(800), rng.Next());
                text = x + x + x;
            }
            else if (pick == 1)
            {
                text = new string('a', 1000 + rng.NextInt(5000)) + (rng.NextInt(2) == 0 ? "b" : string.Empty);
            }
            else
            {
                int n = sizes[rng.NextInt(sizes.Length)];
                List<(string Name, string Text)> all = [.. StringFamilies.Adversarial(n, rng.Next()).Where(f => f.Text.Length > 0 && f.Text.All(c => c <= 255))];
                text = all[rng.NextInt(all.Count)].Text;
            }
            string d = Run("z");
            EndToEndOptions o = RandomOptions(d, text, rng);
            EndToEndResult res = PipelineChecks.RunAndCompare(ChildReference.Of(text), o, $"fuzz case {cases} (seed {seed}): n={text.Length} K={o.Tuning.BucketCapacity} X={o.Period} threads={o.Threads} R={o.Tuning.LcpChunk}");
            Directory.Delete(d, true);
            cases++;
            if (res.Dcs.MaxLevel >= 1 || res.Lcp.TextChunks > 4 || res.Build.Partitions > 0) multi++;
            if (res.Scan.Spill.BlocksWritten > 0) spill++;
        }
        return (cases, multi, spill);
    }

    [Test]
    public void ASeededSample_Texts1kTo20k_EqualTheInRamStructures()
    {
        (long cases, long multi, long spill) = Fuzz(5, c => c < 40, [1000, 3000, 10_000, 20_000]);
        Assert.That(cases, Is.EqualTo(40));
        Assert.That(multi, Is.GreaterThan(cases * 3 / 4), "most cases cross a partition, a recursion level or a chunk boundary");
        TestContext.Out.WriteLine($"{cases} cases, {multi} multi-bucket / multi-level / multi-chunk, {spill} with a stack spill");
    }

    [Test]
    [Category("Fuzz")]
    [Explicit("Differential fuzz of the whole pipeline; run with --filter Category=Fuzz (see PROGRESS-S1.md, Log E6)")]
    public void TinyParameters_Texts1kTo100k_EqualTheInRamStructures()
    {
        ulong seed = ulong.Parse(Environment.GetEnvironmentVariable("EMPROTO_FUZZ_SEED") ?? "1", CultureInfo.InvariantCulture);
        double seconds = double.Parse(Environment.GetEnvironmentVariable("EMPROTO_FUZZ_SECONDS") ?? "120", CultureInfo.InvariantCulture);
        var clock = System.Diagnostics.Stopwatch.StartNew();
        (long cases, long multi, long spill) = Fuzz(seed, _ => clock.Elapsed.TotalSeconds < seconds, [1000, 3000, 10_000, 30_000, 100_000]);
        TestContext.Out.WriteLine($"fuzz: {cases} cases ({multi} multi-bucket / multi-level / multi-chunk, {spill} with a spill), {clock.Elapsed.TotalSeconds:F0} s");
        Assert.That(cases, Is.GreaterThan(0));
    }
}
