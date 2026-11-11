using System.Globalization;
using SuffixTree.Esa.Tests.Fuzz;
using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// Seeded differential fuzz of the facade (every member of <see cref="ISuffixTree"/>): the text of a case comes from the generator families of the S1a harness (<see cref="FuzzTexts"/>),
/// the child layout is drawn per case, short texts are checked against the definitions and the Ukkonen and v6 trees, longer ones against the Ukkonen tree. Category <c>Fuzz</c>,
/// <c>[Explicit]</c>; environment <c>ESA_FUZZ_SEED</c>, <c>ESA_FUZZ_SCALE</c> as for <see cref="DifferentialFuzzTests"/>. A failure names the case id and the layout.
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Differential fuzz of the facade; run with --filter Category=Fuzz (see PROGRESS-S1.md, Log 1b.6)")]
public class FacadeFuzzTests
{
    [Test]
    public void Texts_UpTo120_EveryMember_AgainstDefinitionsAndBothOracles_2500Cases() => Run(0, 120, 2_500, 21, true);

    [Test]
    public void Texts_UpTo3000_AgainstUkkonen_300Cases() => Run(121, 3000, 300, 22, false);

    private static void Run(int minN, int maxN, long cases, ulong salt, bool definitions)
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
            SuffixTreeBackend backend = SuffixTreeBackends.EsaBackends[rng.NextInt(SuffixTreeBackends.EsaBackends.Count)];
            string text = c.Text();
            ISuffixTree tree = backend.FromString(text);
            try
            {
                UkkonenTree ukkonen = UkkonenTree.Build(text);
                string context = $"case {c.Id} [{backend.Name}]";
                if (definitions)
                {
                    using var v6 = PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(text));
                    FacadeChecks.CheckEverything(tree, ukkonen, v6, text, h, context, maxPatternLength: 4);
                }
                else
                {
                    FacadeChecks.CheckAgainstTree(tree, ukkonen, "ukkonen", text, context);
                    for (int k = 0; k < 30; k++)
                    {
                        int start = rng.NextInt(text.Length);
                        string pattern = text.Substring(start, Math.Min(1 + rng.NextInt(30), text.Length - start));
                        if (tree.CountOccurrences(pattern) != ukkonen.CountOccurrences(pattern) || !tree.FindAllOccurrences(pattern).Order().SequenceEqual(ukkonen.FindAllOccurrences(pattern).Order()))
                            Assert.Fail($"{context}: search of '{pattern}' differs from the Ukkonen tree's");
                    }
                    if (tree.CountDistinctSubstrings() != ukkonen.CountDistinctSubstrings()) Assert.Fail($"{context}: CountDistinctSubstrings");
                    if (!tree.CountDistinctSubstringsByLength(20).SequenceEqual(ukkonen.CountDistinctSubstringsByLength(20))) Assert.Fail($"{context}: CountDistinctSubstringsByLength");
                    if (tree.LongestRepeatedSubstring().Length != ukkonen.LongestRepeatedSubstring().Length) Assert.Fail($"{context}: LongestRepeatedSubstring length");
                    var a = tree.FindAllLongestRepeatedSubstrings();
                    var b = ukkonen.FindAllLongestRepeatedSubstrings();
                    if (a.Count != b.Count || a.Where((x, j) => x.Substring != b[j].Substring || !x.Positions.SequenceEqual(b[j].Positions)).Any()) Assert.Fail($"{context}: FindAllLongestRepeatedSubstrings");
                    if (!tree.GetAllSuffixes().SequenceEqual(ukkonen.GetAllSuffixes())) Assert.Fail($"{context}: GetAllSuffixes");
                    foreach (string query in FacadeChecks.Queries(text, rng).Take(4))
                    {
                        int minLength = 1 + rng.NextInt(10);
                        if (!tree.FindMaximalExactMatches(query, minLength).SequenceEqual(ukkonen.FindMaximalExactMatches(query, minLength))) Assert.Fail($"{context}: MEM L {minLength} query {query.Length}");
                        if (!tree.FindMaximalUniqueMatches(query, minLength).SequenceEqual(ukkonen.FindMaximalUniqueMatches(query, minLength))) Assert.Fail($"{context}: MUM L {minLength}");
                        if (!tree.FindExactMatchAnchors(query, minLength).Select(x => (x.PositionInQuery, x.Length)).SequenceEqual(ukkonen.FindExactMatchAnchors(query, minLength).Select(x => (x.PositionInQuery, x.Length)))) Assert.Fail($"{context}: anchors L {minLength}");
                        var lcsA = tree.FindAllDistinctLongestCommonSubstrings(query);
                        var lcsB = ukkonen.FindAllDistinctLongestCommonSubstrings(query);
                        if (lcsA.Count != lcsB.Count || lcsA.Where((x, j) => x.Substring != lcsB[j].Substring || !x.PositionsInText.SequenceEqual(lcsB[j].PositionsInText) || !x.PositionsInOther.SequenceEqual(lcsB[j].PositionsInOther)).Any()) Assert.Fail($"{context}: LCS family");
                    }
                }
            }
            catch (AssertionException ex)
            {
                throw new AssertionException($"fuzz case {i} (master seed {seed}, salt {salt}): {ex.Message}", ex);
            }
            finally
            {
                (tree as IDisposable)?.Dispose();
            }
        }
    }
}
