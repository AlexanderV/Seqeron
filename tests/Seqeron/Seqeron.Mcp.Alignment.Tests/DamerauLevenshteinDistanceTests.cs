using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// damerau_levenshtein_distance — delegation to <c>ApproximateMatcher.DamerauLevenshteinDistance</c>
/// (variant "unrestricted") and <c>OptimalStringAlignmentDistance</c> (variant "osa").
/// Expected values: rapidfuzz / jellyfish cases locked in ApproximateMatcher_EditAlignment_Tests.
/// </summary>
[TestFixture]
public class DamerauLevenshteinDistanceTests
{
    [Test]
    public void DamerauLevenshteinDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.DamerauLevenshteinDistance("CA", "ABC"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.DamerauLevenshteinDistance("", "ABC"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.DamerauLevenshteinDistance("CA", null!));
        Assert.Throws<ArgumentException>(() => AlignmentTools.DamerauLevenshteinDistance("CA", "ABC", "levenshtein"));
    }

    [TestCase("CA", "ABC", 3, 2)]
    [TestCase("ab", "ba", 1, 1)]
    [TestCase("a cat", "an act", 2, 2)]
    [TestCase("ACGT", "TGCA", 3, 3)]
    [TestCase("kitten", "sitting", 3, 3)]
    public void DamerauLevenshteinDistance_Variants_EqualRapidfuzz(string a, string b, int osa, int dl)
    {
        Assert.Multiple(() =>
        {
            var u = AlignmentTools.DamerauLevenshteinDistance(a, b);
            Assert.That(u.Distance, Is.EqualTo(dl));
            Assert.That(u.Variant, Is.EqualTo("unrestricted"));
            var o = AlignmentTools.DamerauLevenshteinDistance(a, b, "OSA");
            Assert.That(o.Distance, Is.EqualTo(osa));
            Assert.That(o.Variant, Is.EqualTo("osa"));
            Assert.That(o.Distance, Is.EqualTo(global::Seqeron.Genomics.Alignment.ApproximateMatcher.OptimalStringAlignmentDistance(a, b)));
            Assert.That(u.Distance, Is.EqualTo(global::Seqeron.Genomics.Alignment.ApproximateMatcher.DamerauLevenshteinDistance(a, b)));
        });
    }
}

/// <summary>
/// damerau_levenshtein_distance weighted costs and damerau_alignment — delegation to the weighted
/// Lowrance–Wagner / OSA methods and the transposition-aware tracebacks. Expected values: R stringdist
/// 0.9.12 (osa) and exhaustive Dijkstra (unrestricted), locked in ApproximateMatcher_WeightedDamerau_Tests.
/// </summary>
[TestFixture]
public class DamerauWeightedAndAlignmentToolTests
{
    [TestCase("CA", "ABC", 2, 3, 4, 3, 7, 5)]
    [TestCase("AXB", "BYA", 2, 1, 3, 2, 6, 5)]
    [TestCase("CAAAACA", "CA", 1, 7, 5, 5, 35, 35)]
    public void DamerauLevenshteinDistance_WeightedCosts(string a, string b, int ins, int del, int sub, int tr, int osa, int dl)
    {
        var u = AlignmentTools.DamerauLevenshteinDistance(a, b, "unrestricted", ins, del, sub, tr);
        var o = AlignmentTools.DamerauLevenshteinDistance(a, b, "osa", ins, del, sub, tr);
        Assert.Multiple(() =>
        {
            Assert.That(u.Distance, Is.EqualTo(dl));
            Assert.That(o.Distance, Is.EqualTo(osa));
            var costs = new global::Seqeron.Genomics.Alignment.DamerauCosts(ins, del, sub, tr);
            Assert.That(u.Distance, Is.EqualTo(global::Seqeron.Genomics.Alignment.ApproximateMatcher.DamerauLevenshteinDistance(a, b, costs)));
            Assert.That(o.Distance, Is.EqualTo(global::Seqeron.Genomics.Alignment.ApproximateMatcher.OptimalStringAlignmentDistance(a, b, costs)));
        });
    }

    [Test]
    public void DamerauLevenshteinDistance_InvalidCosts_Throw()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.DamerauLevenshteinDistance("CA", "ABC", transpositionCost: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.DamerauLevenshteinDistance("CA", "ABC", insertionCost: -1));
        Assert.Throws<ArgumentException>(() => AlignmentTools.DamerauLevenshteinDistance("ABC", "BCA", "unrestricted", 2, 2, 4, 1));
        Assert.That(AlignmentTools.DamerauLevenshteinDistance("ABC", "BCA", "osa", 2, 2, 4, 1).Distance, Is.EqualTo(4));
    }

    [Test]
    public void DamerauAlignment_CaToAbc()
    {
        var u = AlignmentTools.DamerauAlignment("CA", "ABC");
        var o = AlignmentTools.DamerauAlignment("CA", "ABC", "osa");
        Assert.Multiple(() =>
        {
            Assert.That(u.Distance, Is.EqualTo(2));
            Assert.That(u.Variant, Is.EqualTo("unrestricted"));
            Assert.That(u.Script, Is.EqualTo("Td"));
            Assert.That(u.TranspositionCount, Is.EqualTo(1));
            Assert.That(u.Operations, Is.EqualTo(new[] { new DamerauOperationDto("transposition", 0, 2, 0, 3, 2) }));
            Assert.That(o.Distance, Is.EqualTo(3));
            Assert.That(o.Variant, Is.EqualTo("osa"));
            Assert.That(o.Script, Is.EqualTo(global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetOptimalStringAlignment("CA", "ABC").Script));
        });
    }

    [Test]
    public void DamerauAlignment_Weighted_DelegatesAndValidates()
    {
        var a = AlignmentTools.DamerauAlignment("CA", "ABC", "unrestricted", 2, 3, 4, 3);
        var lib = global::Seqeron.Genomics.Alignment.ApproximateMatcher.GetDamerauLevenshteinAlignment(
            "CA", "ABC", new global::Seqeron.Genomics.Alignment.DamerauCosts(2, 3, 4, 3));
        Assert.That(a.Distance, Is.EqualTo(5));
        Assert.That(a.Script, Is.EqualTo(lib.Script));
        Assert.That(a.Operations.Sum(op => op.Cost), Is.EqualTo(5));
        Assert.Throws<ArgumentException>(() => AlignmentTools.DamerauAlignment("", "ABC"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.DamerauAlignment("CA", "ABC", "levenshtein"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.DamerauAlignment("ABC", "BCA", "unrestricted", 2, 2, 4, 1));
    }
}

/// <summary>Doc example locks for damerau_alignment (docs/mcp/tools/alignment/damerau_alignment.md).</summary>
[TestFixture]
public class DamerauAlignmentDocExampleTests
{
    [Test]
    public void DamerauAlignment_OsaAdjacentSwap_DocExample()
    {
        var o = AlignmentTools.DamerauAlignment("ACGT", "AGCT", "osa");
        Assert.Multiple(() =>
        {
            Assert.That(o.Distance, Is.EqualTo(1));
            Assert.That(o.Script, Is.EqualTo("=T="));
            Assert.That(o.Operations, Is.EqualTo(new[]
            {
                new DamerauOperationDto("match", 0, 1, 0, 1, 0),
                new DamerauOperationDto("transposition", 1, 2, 1, 2, 1),
                new DamerauOperationDto("match", 3, 1, 3, 1, 0),
            }));
            Assert.That(AlignmentTools.DamerauLevenshteinDistance("CA", "ABC", "osa", 2, 3, 4, 3).Distance, Is.EqualTo(7));
        });
    }
}
