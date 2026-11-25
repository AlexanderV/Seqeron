using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.Search;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Navigation;

/// <summary>
/// Work package S2.5: MEM, MUM, the anchors and the longest common substrings (the generic algorithms of <c>SuffixTree.Core</c>, 32-bit) over an index that the 32-bit navigator does not serve (5-byte entries, a wide LCP array; the
/// navigator is <see cref="EsaNavigatorWide"/>, built on the long navigator): the answers equal those of the index built in memory, for every layout of the matrix, both access modes, and every combination of the sections that the pack reads in place.
/// </summary>
[TestFixture]
public class WideMatchingTests
{
    private static readonly (bool Psi, bool Rmq)[] Sections = [(false, false), (true, true)];

    public static IEnumerable<string> Texts()
    {
        var rng = new SplitMix64(0x77A1C);
        yield return "mississippimississippi";
        yield return "abracadabraabracadabra";
        yield return "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
        yield return "abababababababababababababab";
        yield return RandomText(rng, 400, 4);
        yield return RandomText(rng, 300, 2);
        string block = RandomText(rng, 50, 4);
        yield return RandomText(rng, 100, 4) + block + RandomText(rng, 100, 4) + block + new string('A', 90);
    }

    private static string RandomText(SplitMix64 rng, int n, int sigma)
    {
        var chars = new char[n];
        for (int i = 0; i < n; i++) chars[i] = (char)('A' + rng.NextInt(sigma));
        return new string(chars);
    }

    private static IEnumerable<string> Queries(string text, SplitMix64 rng)
    {
        yield return text;
        yield return text[(text.Length / 3)..];
        for (int i = 0; i < 6; i++)
        {
            int len = Math.Min(text.Length, 20 + rng.NextInt(60));
            string piece = text.Substring(rng.NextInt(text.Length - len + 1), len);
            var chars = piece.ToCharArray();
            chars[rng.NextInt(chars.Length)] = 'Z'; // a mismatch in the middle: two maximal matches
            yield return new string(chars);
            yield return piece + "ZZ" + text[..Math.Min(text.Length, 15)];
        }
    }

    [Test]
    public void TheMatchingStatistics_OverEveryLayoutAndEverySectionCombination_EqualTheMemoryIndex([ValueSource(nameof(Texts))] string text)
    {
        EsaIndex memory = EsaIndex.Build(text);
        LcpChildCase c = LcpChildCase.Of(text);
        char[] alphabet = [.. text.Distinct().Order()];
        var rng = new SplitMix64(0xCAFE);
        string[] queries = [.. Queries(text, rng)];
        foreach ((int width, LcpKind kind, bool tagged) in WideImage.Layouts)
        {
            using var image = new WideImage(c, width, kind, tagged);
            foreach (EsaAccessMode mode in new[] { EsaAccessMode.Checked, EsaAccessMode.Trusted })
            {
                foreach ((bool psi, bool rmq) in Sections)
                {
                    EsaIndex wide = image.Index(mode == EsaAccessMode.Checked, alphabet, psi, rmq);
                    Assert.That(wide.IsNarrow, Is.False, "the views are of the wide forms: the long navigator serves them");
                    string label = $"width {width} {kind}{(tagged ? " tagged" : string.Empty)} {mode} psi={psi} rmq={rmq}";
                    foreach (string q in queries)
                    {
                        Assert.That(wide.FindMaximalExactMatches(q, 4), Is.EqualTo(memory.FindMaximalExactMatches(q, 4)), $"{label} MEM");
                        Assert.That(wide.FindMaximalUniqueMatches(q, 4), Is.EqualTo(memory.FindMaximalUniqueMatches(q, 4)), $"{label} MUM");
                        Assert.That(wide.FindMaximalUniqueMatches(q, 4, MumUniqueness.Reference), Is.EqualTo(memory.FindMaximalUniqueMatches(q, 4, MumUniqueness.Reference)), $"{label} MUM reference");
                        Assert.That(wide.FindExactMatchAnchors(q, 4), Is.EqualTo(memory.FindExactMatchAnchors(q, 4)), $"{label} anchors");
                        Assert.That(wide.LongestCommonSubstringInfo(q), Is.EqualTo(memory.LongestCommonSubstringInfo(q)), $"{label} LCS");
                    }
                    string other = queries[3];
                    var (s1, p1, o1) = wide.FindAllLongestCommonSubstrings(other);
                    var (s2, p2, o2) = memory.FindAllLongestCommonSubstrings(other);
                    Assert.That((s1, p1, o1), Is.EqualTo((s2, p2, o2)), $"{label} all LCS");
                    Assert.That(wide.FindAllDistinctLongestCommonSubstrings(other).Select(t => (t.Substring, t.PositionsInText, t.PositionsInOther)), Is.EqualTo(memory.FindAllDistinctLongestCommonSubstrings(other).Select(t => (t.Substring, t.PositionsInText, t.PositionsInOther))), $"{label} distinct LCS");
                }
            }
        }
    }

    [Test]
    public void ANavigatorOverAnIndexThatTheThirtyTwoBitNavigatorServes_IsTheNavigatorOfS1_AndTheOtherWayRound()
    {
        EsaIndex memory = EsaIndex.Build("mississippi");
        Assert.That(memory.IsNarrow, Is.True);
        _ = memory.CreateNavigator();
        LcpChildCase c = LcpChildCase.Of("mississippi");
        using var image = new WideImage(c, 5, LcpKind.U8, tagged: false);
        EsaIndex wide = image.Index(true, [.. "mississippi".Distinct().Order()], psi: true, rmq: true);
        Assert.That(wide.IsNarrow, Is.False);
        NotSupportedException e = Assert.Throws<NotSupportedException>(() => wide.CreateNavigator())!;
        Assert.That(e.Message, Does.Contain("EsaNavigatorWide"));
        _ = wide.CreateWideNavigator();
    }
}
