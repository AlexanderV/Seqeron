using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

/// <summary>
/// suffix_tree_find_mums — delegation to <c>ISuffixTreeAnalysis.FindMaximalUniqueMatches</c>.
/// Expected literals are the MUMmer 3.23 outputs locked in SuffixTree.Tests MaximalMatchTests
/// (mummer -mum / -mumreference -l L).
/// </summary>
[TestFixture]
[Category("McpCore")]
public class SuffixTreeFindMumsTests
{
    private const string Reference = "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCG";
    private const string Query = "TTTCCTCATGCAAATAGCCTCATGCATCCGTAATGTAG";

    private static (int, int, int)[] Items(SuffixTreeMaximalMatchesResult r)
        => r.Matches.Select(m => (m.PositionInText, m.PositionInQuery, m.Length)).ToArray();

    [Test]
    public void SuffixTreeFindMums_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindMums("", "ACGT", 3));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindMums("ACGT", "", 3));
        Assert.Throws<ArgumentOutOfRangeException>(() => SuffixTreeCoreTools.SuffixTreeFindMums("ACGT", "ACGT", 0));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindMums("ACGT", "ACGT", 1, "query"));
        Assert.Throws<ArgumentException>(() => SuffixTreeCoreTools.SuffixTreeFindMums("ACGT", "ACGT", 1, null!));
    }

    [Test]
    public void SuffixTreeFindMums_Both_MatchesMummerMum()
    {
        // mummer -mum -l 4: 1 1 13 / 16 11 4 / 26 27 12
        var r = SuffixTreeCoreTools.SuffixTreeFindMums(Reference, Query, 4, "both");
        Assert.That(Items(r), Is.EqualTo(new[] { (0, 0, 13), (15, 10, 4), (25, 26, 12) }));
        Assert.That(Items(SuffixTreeCoreTools.SuffixTreeFindMums(Reference, Query, 4)), Is.EqualTo(Items(r)), "default = both");
    }

    [Test]
    public void SuffixTreeFindMums_Reference_MatchesMummerMumreference_AndEqualsLibrary()
    {
        // mummer -mumreference -l 4: 1 1 13 / 16 11 4 / 4 18 9 / 26 27 12
        var r = SuffixTreeCoreTools.SuffixTreeFindMums(Reference, Query, 4, "Reference");
        Assert.That(Items(r), Is.EqualTo(new[] { (0, 0, 13), (15, 10, 4), (3, 17, 9), (25, 26, 12) }));

        var library = global::SuffixTree.SuffixTree.Build(Reference)
            .FindMaximalUniqueMatches(Query, 4, MumUniqueness.Reference)
            .Select(m => (m.PositionInText, m.PositionInQuery, m.Length)).ToArray();
        Assert.That(Items(r), Is.EqualTo(library));
    }

    [Test]
    public void SuffixTreeFindMums_LengthOne_KeepsMumAtReferenceStart()
    {
        // Definition kept over MUMmer's -mum -l 1 sweep artefact (MaximalMatchTests).
        Assert.That(Items(SuffixTreeCoreTools.SuffixTreeFindMums("ACGT", "TTAG", 1, "both")),
            Is.EqualTo(new[] { (0, 2, 1), (2, 3, 1) }));
    }
}
