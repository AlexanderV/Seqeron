using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class CountApproximateOccurrencesTests
{
    [Test]
    public void CountApproximateOccurrences_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.CountApproximateOccurrences("ATGCATGC", "ATGC", 1));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CountApproximateOccurrences("", "ATGC", 1));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CountApproximateOccurrences(null!, "ATGC", 1));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CountApproximateOccurrences("ATGC", "", 1));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.CountApproximateOccurrences("ATGC", "AT", -1));
    }

    [Test]
    public void CountApproximateOccurrences_Binding_InvokesSuccessfully()
    {
        // Exact matches
        var exact = SuffixTreeTools.CountApproximateOccurrences("ATGATGATG", "ATG", 0);
        Assert.That(exact.Count, Is.EqualTo(3));

        // With 1 mismatch allowed
        var approx = SuffixTreeTools.CountApproximateOccurrences("ATGCTG", "ATG", 1);
        Assert.That(approx.Count, Is.GreaterThanOrEqualTo(1));

        // No matches
        var none = SuffixTreeTools.CountApproximateOccurrences("AAAA", "GGG", 0);
        Assert.That(none.Count, Is.EqualTo(0));
    }
}
