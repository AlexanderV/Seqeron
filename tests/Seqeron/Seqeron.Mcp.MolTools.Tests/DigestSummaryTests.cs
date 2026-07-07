using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class DigestSummaryTests
{
    [Test]
    public void DigestSummary_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.digest_summary("AAAGAATTCAAA", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.digest_summary("", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.digest_summary("AAAGAATTCAAA", Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.digest_summary("AAAGAATTCAAA", null!));
    }

    [Test]
    public void DigestSummary_Binding_InvokesSuccessfully()
    {
        // EcoRI on "AAAGAATTCAAA" -> fragments of length 4 and 8.
        var summary = MolToolsTools.digest_summary("AAAGAATTCAAA", new[] { "EcoRI" });

        Assert.Multiple(() =>
        {
            Assert.That(summary.TotalFragments, Is.EqualTo(2));
            // Sizes are sorted descending.
            Assert.That(summary.FragmentSizes, Is.EqualTo(new[] { 8, 4 }));
            Assert.That(summary.LargestFragment, Is.EqualTo(8));
            Assert.That(summary.SmallestFragment, Is.EqualTo(4));
            Assert.That(summary.AverageFragmentSize, Is.EqualTo(6.0).Within(1e-9));
            Assert.That(summary.EnzymesUsed, Is.EqualTo(new[] { "EcoRI" }));
        });
    }
}
