using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class EditDistanceTests
{
    [Test]
    public void EditDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SuffixTreeTools.EditDistance("ATGC", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.EditDistance("", "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.EditDistance(null!, "ATGC"));
        Assert.Throws<ArgumentException>(() => SuffixTreeTools.EditDistance("ATGC", ""));
    }

    [Test]
    public void EditDistance_Binding_InvokesSuccessfully()
    {
        var identical = SuffixTreeTools.EditDistance("ATGC", "ATGC");
        Assert.That(identical.Distance, Is.EqualTo(0));

        var oneEdit = SuffixTreeTools.EditDistance("ATGC", "ATGG");
        Assert.That(oneEdit.Distance, Is.EqualTo(1));

        var differentLengths = SuffixTreeTools.EditDistance("ATGC", "ATG");
        Assert.That(differentLengths.Distance, Is.EqualTo(1)); // One deletion

        var classic = SuffixTreeTools.EditDistance("kitten", "sitting");
        Assert.That(classic.Distance, Is.EqualTo(3)); // k->s, e->i, +g
    }
}
