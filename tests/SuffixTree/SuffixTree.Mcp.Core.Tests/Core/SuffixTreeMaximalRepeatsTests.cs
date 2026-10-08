using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

/// <summary>
/// suffix_tree_maximal_repeats — delegation to <c>ISuffixTree.FindMaximalRepeatedPairs</c>.
/// Expected literals are the MUMmer 3.23 <c>repeat-match -f -n L</c> outputs locked in SuffixTree.Tests
/// MaximalRepeatedPairsTests (converted to 0-based, ordered by first then second position).
/// </summary>
[TestFixture]
[Category("McpCore")]
public class SuffixTreeMaximalRepeatsTests
{
    private static (int, int, int)[] Items(SuffixTreeMaximalRepeatsResult r)
        => r.Pairs.Select(p => (p.FirstPosition, p.SecondPosition, p.Length)).ToArray();

    [Test]
    public void SuffixTreeMaximalRepeats_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeMaximalRepeats("", 3));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeMaximalRepeats(null!, 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeCoreTools.SuffixTreeMaximalRepeats("ACGT", 0));
    }

    [Test]
    public void SuffixTreeMaximalRepeats_GusfieldExample_MatchesRepeatMatch()
    {
        // repeat-match -f -n 3: 2 10 4 / 6 10 3 / 2 6 3
        var r = SuffixTreeCoreTools.SuffixTreeMaximalRepeats("xabcyabcwabcyz", 3);
        Assert.That(Items(r), Is.EqualTo(new[] { (1, 5, 3), (1, 9, 4), (5, 9, 3) }));
    }

    [Test]
    public void SuffixTreeMaximalRepeats_UniqueSymbols_AndNonAcgt()
    {
        // repeat-match -f -n 3: 1 6 9 / 1 11 4 (N matches N)
        Assert.That(Items(SuffixTreeCoreTools.SuffixTreeMaximalRepeats("ACGTNACGTNACGT", 3)),
            Is.EqualTo(new[] { (0, 5, 9), (0, 10, 4) }));
        var expected = new[] { (0, 5, 4), (0, 10, 4), (5, 10, 4) };
        Assert.That(Items(SuffixTreeCoreTools.SuffixTreeMaximalRepeats("ACGTNACGTNACGT", 3, uniqueSymbols: "N")), Is.EqualTo(expected));
        Assert.That(Items(SuffixTreeCoreTools.SuffixTreeMaximalRepeats("ACGTNACGTNACGT", 3, nonAcgtUnique: true)), Is.EqualTo(expected));

        var library = global::SuffixTree.SuffixTree.Build("ACGTNACGTNACGT").FindMaximalRepeatedPairs(3, c => c == 'N')
            .Select(p => (p.FirstPosition, p.SecondPosition, p.Length)).ToArray();
        Assert.That(library, Is.EqualTo(expected));
    }

    [Test]
    public void SuffixTreeMaximalRepeats_DefaultMinLength_Is20()
        => Assert.That(SuffixTreeCoreTools.SuffixTreeMaximalRepeats("xabcyabcwabcyz").Pairs, Is.Empty);
}
