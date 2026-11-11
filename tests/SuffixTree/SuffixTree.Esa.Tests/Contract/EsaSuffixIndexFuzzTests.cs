using System.Globalization;
using SuffixTree.Esa.Tests.Fuzz;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>
/// Seeded differential fuzz of <see cref="EsaSuffixIndex"/> against the definitional oracle: the text of a case comes from the generator families of the S1a harness
/// (<see cref="FuzzTexts"/>), the child layout is drawn per case, short texts are checked at every locus, longer ones on a seeded sample. Category <c>Fuzz</c>, <c>[Explicit]</c>;
/// environment <c>ESA_FUZZ_SEED</c>, <c>ESA_FUZZ_SCALE</c> as for <see cref="DifferentialFuzzTests"/>. A failure names the case id and the layout.
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Differential fuzz of the index contract; run with --filter Category=Fuzz (see PROGRESS-S1.md, Log 1b.5)")]
public class EsaSuffixIndexFuzzTests
{
    [Test]
    public void Texts_UpTo45_EveryLocus_3000Cases() => Run(0, 45, 3_000, 11, CheckOptions.Thorough);

    [Test]
    public void Texts_UpTo300_SampledLoci_250Cases() => Run(46, 300, 250, 12, CheckOptions.Sampled);

    private static void Run(int minN, int maxN, long cases, ulong salt, CheckOptions options)
    {
        ulong seed = ulong.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SEED") ?? "1", CultureInfo.InvariantCulture);
        double scale = double.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SCALE") ?? "1", CultureInfo.InvariantCulture);
        long count = Math.Max(1, (long)(cases * scale));
        for (long i = 0; i < count; i++)
        {
            ulong h = seed ^ (salt * 0xA24BAED4963EE407UL) ^ ((ulong)i * 0x9FB21C651E98DF25UL);
            var rng = new SplitMix64(new SplitMix64(h).Next());
            int n = minN + rng.NextInt(maxN - minN + 1);
            FuzzCase c = FuzzTexts.Pick(rng, n);
            ChildLayoutKind layout = EsaSuffixIndexTests.Layouts[rng.NextInt(EsaSuffixIndexTests.Layouts.Length)];
            string text = c.Text();
            var (index, _) = EsaSuffixIndexTests.Make(text, layout);
            EsaContractChecker.CheckAll(index, new IndexContractOracle(text), $"case {c.Id} [{layout}]", options with { Seed = h });
        }
    }
}
