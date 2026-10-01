using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

/// <summary>
/// suffix_tree_find_mems — delegation to <c>ISuffixTreeAnalysis.FindMaximalExactMatches</c>.
/// Expected literals are the MUMmer 3.23 outputs locked in SuffixTree.Tests MaximalMatchTests
/// (mummer -maxmatch -l L, converted to 0-based (text, query, length)).
/// </summary>
[TestFixture]
[Category("McpCore")]
public class SuffixTreeFindMemsTests
{
    private static (int, int, int)[] Items(SuffixTreeMaximalMatchesResult r)
        => r.Matches.Select(m => (m.PositionInText, m.PositionInQuery, m.Length)).ToArray();

    [Test]
    public void SuffixTreeFindMems_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindMems("", "ACGT", 3));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindMems(null!, "ACGT", 3));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindMems("ACGT", "", 3));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindMems("ACGT", null!, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeCoreTools.SuffixTreeFindMems("ACGT", "ACGT", 0));
    }

    [Test]
    public void SuffixTreeFindMems_GattacaRepeat_MatchesMummer()
    {
        // mummer -maxmatch -l 3: 2 1 7 / 9 1 6 / 1 8 7 / 8 8 7
        var r = SuffixTreeCoreTools.SuffixTreeFindMems("GATTACAGATTACA", "ATTACAGGATTACAT", 3);
        Assert.That(Items(r), Is.EqualTo(new[] { (1, 0, 7), (8, 0, 6), (0, 7, 7), (7, 7, 7) }));
    }

    [Test]
    public void SuffixTreeFindMems_ThreeSetsDiffer_MatchesMummer_AndEqualsLibrary()
    {
        const string reference = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG";
        const string query = "TTTCCTCATGCAAATAGCCTCATGCATCCGTAATGTAG";
        var r = SuffixTreeCoreTools.SuffixTreeFindMems(reference, query, 4);
        Assert.That(Items(r), Is.EqualTo(new[]
            { (0, 0, 13), (21, 6, 4), (15, 10, 4), (3, 17, 9), (21, 20, 4), (25, 26, 12), (22, 32, 4) }));

        var library = global::SuffixTree.SuffixTree.Build(reference).FindMaximalExactMatches(query, 4)
            .Select(m => (m.PositionInText, m.PositionInQuery, m.Length)).ToArray();
        Assert.That(Items(r), Is.EqualTo(library));
    }

    [Test]
    public void SuffixTreeFindMems_DefaultMinLength_Is20()
    {
        // Short inputs have no MEM of length >= 20 (MUMmer default -l 20).
        Assert.That(SuffixTreeCoreTools.SuffixTreeFindMems("GATTACAGATTACA", "ATTACAGGATTACAT").Matches, Is.Empty);
    }
}
