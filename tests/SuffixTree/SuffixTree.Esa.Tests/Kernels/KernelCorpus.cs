using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>The texts every registered engine is run on (default run: small and medium; the explicit categories add the large ones).</summary>
internal static class KernelCorpus
{
    /// <summary>All strings over sigma 2 up to length 10, sigma 3 up to 7, sigma 4 up to 5, 30 random texts n &lt;= 2000 per alphabet size, the adversarial families at n = 0..5, 64, 2000 and 20,000.</summary>
    internal static IEnumerable<(string Name, string Text)> Default()
    {
        foreach (var (sigma, maxLen) in new[] { (2, 10), (3, 7), (4, 5) })
            foreach (string s in ExhaustiveStrings.All(sigma, maxLen))
                yield return ($"exh{sigma}", s);

        foreach (int sigma in new[] { 1, 2, 3, 4, 20, 256, 257, 1000 })
        {
            var rng = new SplitMix64(0xE5A0UL + (ulong)sigma);
            for (int k = 0; k < 30; k++)
            {
                int n = rng.NextInt(2001);
                yield return ($"random-s{sigma}-{k}", StringFamilies.Random(n, sigma, rng.Next()));
            }
        }

        foreach (int n in new[] { 0, 1, 2, 3, 4, 5, 64, 2000, 20_000 })
            foreach (var (name, text) in StringFamilies.Adversarial(n, 1))
                yield return ($"{name}@{n}", text);
    }

    /// <summary>The families at 10^5 characters and sigma = n = 65,536 (explicit categories and the A/B acceptance).</summary>
    internal static IEnumerable<(string Name, string Text)> Large()
    {
        foreach (var (name, text) in StringFamilies.Adversarial(100_000, 1)) yield return ($"{name}@1e5", text);
        yield return ("all-distinct@65536", StringFamilies.AllDistinct(65_536, 1));
    }
}
