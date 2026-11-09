using SuffixTree.Esa.Tests.TestSupport;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Search;

/// <summary>
/// Heavy checks of work package 1b.2. Category <c>Exhaustive</c>, <c>[Explicit]</c>:
/// <c>dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter Category=Exhaustive</c>.
/// Every string over sigma 2 up to length 14, sigma 3 up to 9, sigma 4 up to 7: count and positions of every distinct substring, and of a few absent patterns, equal
/// the Ukkonen tree for all three child-table layouts (the persistent v6 tree up to length 11 / 7 / 5); random and adversarial texts up to 2000 symbols.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log 1b.2)")]
public class EsaSearchExhaustiveTests
{
    private static void CheckString(string s, bool withV6)
    {
        UkkonenTree ukkonen = UkkonenTree.Build(s);
        using PersistentSuffixTree? v6 = withV6 ? PersistentSuffixTreeFactory.CreatePersistent(new StringTextSource(s)) : null;
        foreach (ChildLayoutKind layout in Enum.GetValues<ChildLayoutKind>())
        {
            EsaIndex esa = EsaIndex.Build(s, layout);
            var seen = new HashSet<string>(StringComparer.Ordinal);
            var patterns = new List<string>();
            for (int i = 0; i < s.Length; i++)
                for (int l = 1; i + l <= s.Length; l++)
                {
                    string pattern = s.Substring(i, l);
                    if (seen.Add(pattern)) patterns.Add(pattern);
                }
            patterns.AddRange([s + "a", "z" + s, "az", string.Concat(s, s)]);
            foreach (string pattern in patterns.Where(p => p.Length > 0))   // the empty pattern is pinned in EsaSearchTests
            {
                int expected = ukkonen.CountOccurrences(pattern);
                if (esa.CountOccurrences(pattern) != expected) Assert.Fail($"count of '{pattern}' in '{s}' with {layout}");
                if (esa.Contains(pattern) != expected > 0) Assert.Fail($"Contains '{pattern}' in '{s}' with {layout}");
                if (layout == ChildLayoutKind.Tagged)
                {
                    if (!esa.FindAllOccurrences(pattern).OrderBy(x => x).SequenceEqual(ukkonen.FindAllOccurrences(pattern).OrderBy(x => x))) Assert.Fail($"positions of '{pattern}' in '{s}'");
                    if (v6 is not null && v6.CountOccurrences(pattern) != expected) Assert.Fail($"v6 count of '{pattern}' in '{s}' (oracle disagreement)");
                }
            }
        }
    }

    [Test]
    public void AllBinaryStrings_UpTo14() => CheckAll(2, 14, 11);

    [Test]
    public void AllTernaryStrings_UpTo9() => CheckAll(3, 9, 7);

    [Test]
    public void AllFourLetterStrings_UpTo7() => CheckAll(4, 7, 5);

    private static void CheckAll(int sigma, int maxLength, int v6MaxLength)
    {
        foreach (string s in ExhaustiveStrings.All(sigma, maxLength)) CheckString(s, withV6: s.Length <= v6MaxLength);
    }

    [Test]
    public void RandomTexts_UpTo2000_EverySampledSubstring()
    {
        var rng = new SplitMix64(0x1B31);
        for (int i = 0; i < 12; i++)
        {
            int n = 200 + rng.NextInt(1801);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 9 == 0))
            {
                UkkonenTree ukkonen = UkkonenTree.Build(text);
                EsaIndex esa = EsaIndex.Build(text);
                var seen = new HashSet<string>(StringComparer.Ordinal);
                int checkedPatterns = 0;
                for (int start = 0; start < text.Length && checkedPatterns < 6000; start += 1 + rng.NextInt(3))
                    for (int len = 1; len <= Math.Min(60, text.Length - start); len++)
                    {
                        string pattern = text.Substring(start, len);
                        if (!seen.Add(pattern)) continue;
                        checkedPatterns++;
                        if (esa.CountOccurrences(pattern) != ukkonen.CountOccurrences(pattern)) Assert.Fail($"{name} n {n}: count of '{pattern}'");
                        if (!esa.FindAllOccurrences(pattern).OrderBy(x => x).SequenceEqual(ukkonen.FindAllOccurrences(pattern).OrderBy(x => x))) Assert.Fail($"{name} n {n}: positions of '{pattern}'");
                    }
            }
        }
    }
}
