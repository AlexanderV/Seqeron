using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class FrequentKmersWithMismatchesTests
{
    [Test]
    public void FrequentKmersWithMismatches_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.FrequentKmersWithMismatches("AAAAA", 2, 0));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FrequentKmersWithMismatches("", 2, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.FrequentKmersWithMismatches("AAAAA", 0, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.FrequentKmersWithMismatches("AAAAA", 2, -1));
    }

    [Test]
    public void FrequentKmersWithMismatches_Binding_InvokesSuccessfully()
    {
        // d=0 reduces to plain frequent words. "AAAAA" has four 2-mer windows, all "AA".
        var r = AlignmentTools.FrequentKmersWithMismatches("AAAAA", 2, 0);
        Assert.Multiple(() =>
        {
            Assert.That(r.Items, Has.Length.EqualTo(1));
            Assert.That(r.Items[0].Kmer, Is.EqualTo("AA"));
            Assert.That(r.Items[0].Count, Is.EqualTo(4));
        });
    }
}
