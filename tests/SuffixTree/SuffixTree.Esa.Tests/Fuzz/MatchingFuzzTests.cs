using System.Globalization;
using SuffixTree.Esa.Tests.Navigation;
using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// Differential fuzz of the matching operations of work package 1b.4 (MEM, MUM, anchors, the LCS family on <see cref="EsaNavigator"/>): the text of a case comes from the generator
/// families of the S1a fuzz harness (<see cref="FuzzTexts"/>), the query from a seeded recipe (independent text of the same family, mutated window, mutated copy, reversed window, text over a wider
/// alphabet). Short texts are checked against the definition-level implementations and the Ukkonen tree, long ones against the Ukkonen tree. Category <c>Fuzz</c>, <c>[Explicit]</c>;
/// environment: <c>ESA_FUZZ_SEED</c>, <c>ESA_FUZZ_SCALE</c> as for <see cref="DifferentialFuzzTests"/>. A failure names the case id and the recipe.
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Differential fuzz of the matching operations; run with --filter Category=Fuzz (see PROGRESS-S1.md, Log 1b.4)")]
public class MatchingFuzzTests
{
    [Test]
    public void Texts_UpTo200_AgainstDefinitionAndUkkonen_20000Cases() => Run(0, 200, 20_000, 1, true);

    [Test]
    public void Texts_UpTo4000_AgainstUkkonen_2000Cases() => Run(201, 4000, 2_000, 2, false);

    private static void Run(int minN, int maxN, long cases, ulong salt, bool definition)
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
            string text = c.Text();
            int recipe = rng.NextInt(5);
            string query = recipe switch
            {
                0 => FuzzTexts.Generate(new FuzzCase(c.Family, rng.Next(), 1 + rng.NextInt(Math.Max(1, n)), c.Sigma, c.Param)),
                1 when text.Length > 0 => MatchingChecks.Mutate(Window(text, rng), 0.03, rng),
                2 => MatchingChecks.Mutate(text, 0.1, rng, "\u0001"),
                3 when text.Length > 0 => new string(Window(text, rng).Reverse().ToArray()),
                _ => FuzzTexts.Generate(new FuzzCase("uniform", rng.Next(), 1 + rng.NextInt(60), Math.Min(c.Sigma + 1, 1000), 0)),
            };
            int[] minLengths = [1 + rng.NextInt(4), 4 + rng.NextInt(12)];
            try
            {
                EsaIndex esa = EsaIndex.Build(text);
                UkkonenTree ukkonen = UkkonenTree.Build(text);
                if (definition)
                {
                    MatchingChecks.CheckAll(esa, ukkonen, text, query, minLengths, $"case {i}");
                }
                else
                {
                    foreach (int minLength in minLengths)
                    {
                        var mems = MatchingChecks.Triples(esa.FindMaximalExactMatches(query, minLength));
                        if (!mems.SequenceEqual(MatchingChecks.Triples(ukkonen.FindMaximalExactMatches(query, minLength)))) Assert.Fail($"MEM L {minLength} differs from Ukkonen");
                        foreach (MumUniqueness mode in new[] { MumUniqueness.Both, MumUniqueness.Reference })
                            if (!MatchingChecks.Triples(esa.FindMaximalUniqueMatches(query, minLength, mode)).SequenceEqual(MatchingChecks.Triples(ukkonen.FindMaximalUniqueMatches(query, minLength, mode)))) Assert.Fail($"MUM {mode} L {minLength} differs from Ukkonen");
                        if (!esa.FindExactMatchAnchors(query, minLength).Select(a => (a.PositionInQuery, a.Length)).SequenceEqual(ukkonen.FindExactMatchAnchors(query, minLength).Select(a => (a.PositionInQuery, a.Length)))) Assert.Fail($"anchors L {minLength} differ from Ukkonen");
                    }
                    var lcs = esa.FindAllDistinctLongestCommonSubstrings(query);
                    var reference = ukkonen.FindAllDistinctLongestCommonSubstrings(query);
                    if (lcs.Count != reference.Count || lcs.Where((r, k) => r.Substring != reference[k].Substring || !r.PositionsInText.SequenceEqual(reference[k].PositionsInText) || !r.PositionsInOther.SequenceEqual(reference[k].PositionsInOther)).Any())
                        Assert.Fail("distinct longest common substrings differ from Ukkonen");
                }
            }
            catch (AssertionException ex)
            {
                throw new AssertionException($"fuzz case {i} (master seed {seed}, salt {salt}): text {c.Id}, recipe {recipe}, minLengths {string.Join(",", minLengths)}, query length {query.Length}: {ex.Message}", ex);
            }
        }
    }

    private static string Window(string text, SplitMix64 rng)
    {
        int start = rng.NextInt(text.Length);
        return text.Substring(start, 1 + rng.NextInt(text.Length - start));
    }
}
