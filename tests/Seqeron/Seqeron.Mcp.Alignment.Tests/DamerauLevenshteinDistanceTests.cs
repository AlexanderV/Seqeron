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
