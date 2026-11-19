using System.Globalization;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// Seeded differential fuzz of the container writer of work package 1b.7e: the text of a case comes from the generator families of the S1a harness (<see cref="FuzzTexts"/>, only those whose code units are all below 256),
/// the thread count, the sizes of every stage and the PSI option are drawn per case, and every section and the header are compared with the in-RAM index. Category <c>Fuzz</c>, <c>[Explicit]</c>; environment
/// <c>ESA_FUZZ_SEED</c>, <c>ESA_FUZZ_SCALE</c> as for the other fuzz classes. A failure names the case id.
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Differential fuzz of the container writer; run with --filter FullyQualifiedName~SesaBuildFuzzTests (see PROGRESS-S1.md, Log 1b.7e)")]
public class SesaBuildFuzzTests
{
    [Test]
    public void Texts_UpTo300_TinyParameters_1000Cases() => Run(0, 300, 1_000, 71);

    [Test]
    public void Texts_UpTo20000_TinyParameters_120Cases() => Run(301, 20_000, 120, 72);

    private static void Run(int minN, int maxN, long cases, ulong salt)
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
            if (text.Any(ch => ch > 255))
            {
                skipped++;
                continue;
            }
            int variant = rng.NextInt(1000);
            bool psi = rng.NextInt(4) != 0;
            using var dir = new TempDir();
            SesaRun.Outcome o = SesaRun.Build(dir, text, SesaRun.Tiny(variant, 1 + rng.NextInt(4), psi));
            SesaRun.Compare(o, $"case {c.Id} variant {variant} psi {psi}", psi);
        }
        Assert.That(skipped, Is.LessThan(count), "at least one case ran");
    }
}
