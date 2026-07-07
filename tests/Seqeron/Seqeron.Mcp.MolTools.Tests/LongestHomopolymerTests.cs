using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class LongestHomopolymerTests
{
    [Test]
    public void LongestHomopolymer_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.longest_homopolymer("AAATTGC"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.longest_homopolymer(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.longest_homopolymer(null!));
    }

    [Test]
    public void LongestHomopolymer_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Longest run is "AAA" = 3.
            Assert.That(MolToolsTools.longest_homopolymer("AAATTGC").Length, Is.EqualTo(3));
            // No repeats -> longest run is 1.
            Assert.That(MolToolsTools.longest_homopolymer("ACGT").Length, Is.EqualTo(1));
            // Case-insensitive: "gggg" = 4.
            Assert.That(MolToolsTools.longest_homopolymer("acgggg").Length, Is.EqualTo(4));
            // Whole sequence one base -> its length.
            Assert.That(MolToolsTools.longest_homopolymer("TTTTT").Length, Is.EqualTo(5));
        });
    }
}
