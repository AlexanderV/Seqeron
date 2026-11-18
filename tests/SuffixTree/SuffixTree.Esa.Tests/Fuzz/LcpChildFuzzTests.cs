using System.Globalization;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// Seeded differential fuzz of the stage of work package 1b.7c (PH4 + PH5 on a given suffix array): the text of a case comes from the generator families of the S1a harness (<see cref="FuzzTexts"/>, alphabets up to 256
/// symbols), the chunk size, thread count, stack block, window and codec are drawn per case. Category <c>Fuzz</c>, <c>[Explicit]</c>; environment <c>ESA_FUZZ_SEED</c>, <c>ESA_FUZZ_SCALE</c> as for the other fuzz classes
/// (<c>TMPDIR=/dev/shm</c> makes the many small files cheap). A failure names the case id.
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Differential fuzz of PH4 + PH5; run with --filter FullyQualifiedName~LcpChildFuzzTests (see PROGRESS-S1.md, Log 1b.7c)")]
public class LcpChildFuzzTests
{
    [Test]
    public void Texts_UpTo300_TinyParameters_1500Cases() => Run(0, 300, 1_500, 41, tiny: true);

    [Test]
    public void Texts_UpTo20000_RandomParameters_200Cases() => Run(301, 20_000, 200, 42, tiny: false);

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
            LcpChildCase oracle = LcpChildCase.Of(text);
            ChildCodec codec = rng.NextInt(3) == 0 ? ChildCodec.Tagged : ChildCodec.Packed;
            int threads = 1 + rng.NextInt(4);
            LcpChildTuning tuning = tiny
                ? LcpChildRun.TinyTuning(rng.NextInt(1000)) with { LcpChunk = 1 + rng.NextInt(Math.Max(1, n)) }
                : new LcpChildTuning { LcpChunk = 1 + rng.NextInt(n), LcpTupleBytes = 24 + rng.NextInt(4000), LcpReadBytes = 16 + rng.NextInt(5000), BlockEntries = 2 + rng.NextInt(5000), ChildWindowBytes = 16 << rng.NextInt(12), ChildSlots = 1 + rng.NextInt(8), LcpWorkers = 1 + rng.NextInt(4) };
            using var dir = new TempDir();
            string label = $"case {c.Id} R={tuning.LcpChunk} threads {threads} {codec}";
            LcpChildRun.Outcome o = LcpChildRun.Run(dir, oracle, new LcpChildOptions { ScratchDir = dir.Sub("scratch"), Threads = threads, Codec = codec, Tuning = tuning, Sigma = oracle.Sigma, ArenaGuards = ArenaOptions.Guarded(0xCD) });
            LcpChildRun.Compare(oracle, o, codec, label);
        }
        Assert.That(skipped, Is.LessThan(count), "at least one case ran");
    }
}
