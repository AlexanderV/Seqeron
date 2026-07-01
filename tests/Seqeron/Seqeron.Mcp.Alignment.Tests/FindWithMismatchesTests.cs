using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class FindWithMismatchesTests
{
    [Test]
    public void FindWithMismatches_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.FindWithMismatches("ACGTACGT", "ACG", 0));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FindWithMismatches("", "ACG", 0));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FindWithMismatches("ACGT", null!, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.FindWithMismatches("ACGT", "ACG", -1));
    }

    [Test]
    public void FindWithMismatches_Binding_InvokesSuccessfully()
    {
        // Exact (Hamming distance 0) occurrences of "ACG" in "ACGTACGT": positions 0 and 4.
        var exact = AlignmentTools.FindWithMismatches("ACGTACGT", "ACG", 0);
        Assert.Multiple(() =>
        {
            Assert.That(exact.Items, Has.Length.EqualTo(2));
            Assert.That(exact.Items[0].Position, Is.EqualTo(0));
            Assert.That(exact.Items[1].Position, Is.EqualTo(4));
            Assert.That(exact.Items[0].Distance, Is.EqualTo(0));
        });

        // Allowing 1 mismatch, "AA" matches every window of "AAAA" (positions 0,1,2).
        var approx = AlignmentTools.FindWithMismatches("AAAA", "AA", 1);
        Assert.That(approx.Items, Has.Length.EqualTo(3));
    }
}
