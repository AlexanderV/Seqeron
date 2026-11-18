using System.Globalization;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// Seeded differential fuzz of the stage of work package 1b.7d (PH1 and PH2 / PH3 on a packed text): the text of a case comes from the generator families of the S1a harness (<see cref="FuzzTexts"/>, alphabets up to 256
/// symbols), and the period, K_s, N_base, K, the fan-out, the oversample, the windows, the thread count and the seed are drawn per case. Category <c>Fuzz</c>, <c>[Explicit]</c>; environment <c>ESA_FUZZ_SEED</c>,
/// <c>ESA_FUZZ_SCALE</c> as for the other fuzz classes (<c>TMPDIR=/dev/shm</c> makes the many small files cheap). A failure names the case id.
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Differential fuzz of PH1 and PH2 / PH3; run with --filter FullyQualifiedName~SaBuildFuzzTests (see PROGRESS-S1.md, Log 1b.7d)")]
public class SaBuildFuzzTests
{
    [Test]
    public void Texts_UpTo300_TinyParameters_1500Cases() => Run(0, 300, 1_500, 51, tiny: true);

    [Test]
    public void Texts_UpTo20000_RandomParameters_200Cases() => Run(301, 20_000, 200, 52, tiny: false);

    private static void Run(int minN, int maxN, long cases, ulong salt, bool tiny)
    {
        ulong seed = ulong.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SEED") ?? "1", CultureInfo.InvariantCulture);
        double scale = double.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SCALE") ?? "1", CultureInfo.InvariantCulture);
        long count = Math.Max(1, (long)(cases * scale));
        int skipped = 0;
        for (long i = 0; i < count; i++)
        {
            ulong h = seed ^ (salt * 0xA24BAED4963EE407UL) ^ ((ulong)i * 0x9FB21C651E98DF25UL);
            var rng = new SplitMix64(new SplitMix64(h).Next());
            int n = minN + rng.NextInt(maxN - minN + 1);
            FuzzCase c = FuzzTexts.Pick(rng, n);
            string text = c.Text();
            if (text.Distinct().Count() > 256)
            {
                skipped++;
                continue;
            }
            SaCase oracle = SaCase.Of(text);
            int threads = 1 + rng.NextInt(4);
            int variant = rng.NextInt(1000);
            SaBuildTuning tuning = tiny
                ? SaTiny.Tuning(variant)
                : new SaBuildTuning
                {
                    DcsBucketCapacity = 3 + rng.NextInt(2000),
                    DcsBaseCase = 4 + rng.NextInt(3000),
                    DcsChunkSlots = 16 + rng.NextInt(3000),
                    DcsSortCapacity = 3 + rng.NextInt(2000),
                    DcsRunSpillBytes = rng.NextInt(2000),
                    DcsFanOut = 2 + rng.NextInt(60),
                    DcsReservoir = 3 + rng.NextInt(500),
                    DcsIoBytes = 16 + rng.NextInt(3000),
                    DcsMergeBudget = 8 + rng.NextInt(4000),
                    BucketCapacity = 2 + rng.NextInt(2000),
                    FanOut = 2 + rng.NextInt(64),
                    Oversample = 1 + rng.NextInt(20),
                    SpillBytes = 16 + rng.NextInt(3000),
                    WindowBytes = 1 + rng.NextInt(5000),
                };
            using var dir = new TempDir();
            var options = new SaBuildOptions
            {
                ScratchDir = dir.Sub("scratch"),
                Threads = threads,
                Period = SaTiny.Periods[rng.NextInt(SaTiny.Periods.Length)],
                ExactPeriod = true,
                Tuning = tuning,
                Seed = rng.Next(),
                ArenaGuards = ArenaOptions.Guarded(0xCD),
            };
            string label = $"case {c.Id} X={options.Period} threads {threads} variant {variant}";
            SaRun.Outcome o = SaRun.Run(dir, oracle, options);
            SaRun.Compare(oracle, o, label);
        }
        Assert.That(skipped, Is.LessThan(count), "at least one case ran");
    }
}
