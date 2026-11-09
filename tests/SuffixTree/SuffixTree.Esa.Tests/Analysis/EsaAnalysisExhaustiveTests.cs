using SuffixTree.Esa.Tests.TestSupport;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Analysis;

/// <summary>
/// Heavy checks of work package 1b.3. Category <c>Exhaustive</c>, <c>[Explicit]</c>:
/// <c>dotnet test tests/SuffixTree/SuffixTree.Esa.Tests/SuffixTree.Esa.Tests.csproj -c Release --filter Category=Exhaustive</c>.
/// Every string over sigma 2 up to length 16, sigma 3 up to 10, sigma 4 up to 8: every native analysis operation against the definition, the Ukkonen tree and the
/// <c>Traverse</c>-default path; the maximal repeated pairs through the interface default; then random and adversarial texts up to 3000 symbols against Ukkonen and the
/// <c>Traverse</c>-default path.
/// </summary>
[TestFixture]
[Category("Exhaustive")]
[Explicit("Heavy enumeration; run with --filter Category=Exhaustive (see PROGRESS-S1.md, Log 1b.3)")]
public class EsaAnalysisExhaustiveTests
{
    private static void CheckString(string s)
    {
        EsaIndex esa = EsaIndex.Build(s);
        UkkonenTree ukkonen = UkkonenTree.Build(s);
        var shim = new EsaTreeShim(esa);
        string where = $"'{s}'";

        // distinct substrings
        long total = esa.CountDistinctSubstrings();
        if (total != AnalysisReference.DistinctTotal(s) || total != ukkonen.CountDistinctSubstrings() || total != SuffixTreeAlgorithms.CountDistinctSubstrings(shim)) Assert.Fail($"total of {where}");
        foreach (int m in new[] { 1, 2, s.Length, s.Length + 1 }.Where(m => m >= 1).Distinct())
        {
            long[] actual = esa.CountDistinctSubstringsByLength(m);
            if (!actual.SequenceEqual(AnalysisReference.DistinctByLength(s, m)) || !actual.SequenceEqual(ukkonen.CountDistinctSubstringsByLength(m)) ||
                !actual.SequenceEqual(SuffixTreeAlgorithms.CountDistinctSubstringsByLength(shim, m))) Assert.Fail($"ByLength({m}) of {where}");
        }

        // longest repeats
        var all = esa.FindAllLongestRepeatedSubstrings();
        var reference = AnalysisReference.AllLongestRepeated(s);
        if (all.Count != reference.Count) Assert.Fail($"FindAllLRS of {where}: {all.Count} entries, expected {reference.Count}");
        for (int i = 0; i < all.Count; i++)
            if (all[i].Substring != reference[i].Substring || !all[i].Positions.SequenceEqual(reference[i].Positions)) Assert.Fail($"FindAllLRS of {where}: entry {i}");
        var viaUkkonen = ukkonen.FindAllLongestRepeatedSubstrings();
        var viaTraverse = SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings(shim);
        for (int i = 0; i < all.Count; i++)
        {
            if (all[i].Substring != viaUkkonen[i].Substring || !all[i].Positions.SequenceEqual(viaUkkonen[i].Positions)) Assert.Fail($"FindAllLRS of {where}: Ukkonen entry {i}");
            if (all[i].Substring != viaTraverse[i].Substring || !all[i].Positions.SequenceEqual(viaTraverse[i].Positions)) Assert.Fail($"FindAllLRS of {where}: Traverse entry {i}");
        }
        if (viaUkkonen.Count != all.Count || viaTraverse.Count != all.Count) Assert.Fail($"FindAllLRS of {where}: entry counts of the oracles");
        string lrs = esa.LongestRepeatedSubstring();
        if (lrs.Length != ukkonen.LongestRepeatedSubstring().Length) Assert.Fail($"LRS length of {where}");
        if (all.Count == 0 ? lrs.Length != 0 : lrs != all.Select(e => e.Substring).Min(StringComparer.Ordinal)) Assert.Fail($"LRS of {where}");

        // suffixes
        string[] suffixes = AnalysisReference.SortedSuffixes(s);
        if (!esa.GetAllSuffixes().SequenceEqual(suffixes) || !esa.EnumerateSuffixes().SequenceEqual(suffixes) || !ukkonen.GetAllSuffixes().SequenceEqual(suffixes)) Assert.Fail($"suffixes of {where}");

        // maximal repeated pairs through the interface default
        ISuffixTree viaInterface = new EsaNativeTreeShim(esa);
        foreach (int minLength in new[] { 1, 2, 4 })
            foreach (Func<char, bool>? unique in new Func<char, bool>?[] { null, c => c == 'a', c => c == 'b' })
            {
                var actual = viaInterface.FindMaximalRepeatedPairs(minLength, unique).Select(p => (p.FirstPosition, p.SecondPosition, p.Length));
                if (!actual.SequenceEqual(AnalysisReference.MaximalRepeatedPairs(s, minLength, unique))) Assert.Fail($"pairs L {minLength} of {where}");
            }
    }

    [Test]
    public void AllBinaryStrings_UpTo16() => CheckAll(2, 16);

    [Test]
    public void AllTernaryStrings_UpTo10() => CheckAll(3, 10);

    [Test]
    public void AllFourLetterStrings_UpTo8() => CheckAll(4, 8);

    private static void CheckAll(int sigma, int maxLength)
    {
        long count = 0;
        foreach (string s in ExhaustiveStrings.All(sigma, maxLength))
        {
            CheckString(s);
            count++;
        }
        long expected = 0;
        for (int len = 0, size = 1; len <= maxLength; len++, size *= sigma) expected += size;
        Assert.That(count, Is.EqualTo(expected), "every string of every length up to the bound was checked");
    }

    [Test]
    public void RandomAndAdversarialTexts_UpTo3000_AgainstUkkonenAndTheTraverseDefaultPath()
    {
        var rng = new SplitMix64(0x1B3B01);
        for (int i = 0; i < 14; i++)
        {
            int n = 200 + rng.NextInt(2801);
            ulong seed = rng.Next();
            foreach (var (name, text) in StringFamilies.Adversarial(n, seed).Where((_, k) => (k + i) % 4 == 0))
            {
                EsaIndex esa = EsaIndex.Build(text);
                UkkonenTree ukkonen = UkkonenTree.Build(text);
                var shim = new EsaTreeShim(esa);
                string where = $"{name} n {n}";
                if (esa.CountDistinctSubstrings() != ukkonen.CountDistinctSubstrings() || esa.CountDistinctSubstrings() != SuffixTreeAlgorithms.CountDistinctSubstrings(shim)) Assert.Fail($"total of {where}");
                foreach (int m in new[] { 1, 5, 64, n, n + 3 })
                {
                    long[] actual = esa.CountDistinctSubstringsByLength(m);
                    if (!actual.SequenceEqual(ukkonen.CountDistinctSubstringsByLength(m)) || !actual.SequenceEqual(SuffixTreeAlgorithms.CountDistinctSubstringsByLength(shim, m))) Assert.Fail($"ByLength({m}) of {where}");
                }
                var all = esa.FindAllLongestRepeatedSubstrings();
                var viaUkkonen = ukkonen.FindAllLongestRepeatedSubstrings();
                var viaTraverse = SuffixTreeAlgorithms.FindAllLongestRepeatedSubstrings(shim);
                if (all.Count != viaUkkonen.Count || all.Count != viaTraverse.Count) Assert.Fail($"FindAllLRS count of {where}");
                for (int k = 0; k < all.Count; k++)
                    if (all[k].Substring != viaUkkonen[k].Substring || !all[k].Positions.SequenceEqual(viaUkkonen[k].Positions) || all[k].Substring != viaTraverse[k].Substring ||
                        !all[k].Positions.SequenceEqual(viaTraverse[k].Positions)) Assert.Fail($"FindAllLRS entry {k} of {where}");
                if (esa.LongestRepeatedSubstring().Length != ukkonen.LongestRepeatedSubstring().Length) Assert.Fail($"LRS length of {where}");
                if (n <= 1500 && !esa.GetAllSuffixes().SequenceEqual(ukkonen.GetAllSuffixes())) Assert.Fail($"suffixes of {where}");
            }
        }
    }
}
