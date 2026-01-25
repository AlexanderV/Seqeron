using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class CalculateSimilarityTests
{
    [Test]
    public void CalculateSimilarity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.CalculateSimilarity("ATGCATGC", "ATGCATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CalculateSimilarity("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CalculateSimilarity(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CalculateSimilarity("ATGC", ""));
    }

    [Test]
    public void CalculateSimilarity_Binding_InvokesSuccessfully()
    {
        // Identical sequences should have high similarity
        var identical = SuffixTreeTools.CalculateSimilarity("ATGCATGCATGC", "ATGCATGCATGC");
        Assert.That(identical, Is.Not.Null);
        Assert.That(identical.Similarity, Is.GreaterThan(0.9));

        // Different sequences should have similarity between 0 and 1
        var different = SuffixTreeTools.CalculateSimilarity("ATGCATGC", "AAAATTTT");
        Assert.That(different.Similarity, Is.GreaterThanOrEqualTo(0.0));
        Assert.That(different.Similarity, Is.LessThanOrEqualTo(1.0));
    }
}
