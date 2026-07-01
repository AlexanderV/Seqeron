using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class FindBestMatchTests
{
    [Test]
    public void FindBestMatch_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.FindBestMatch("ACGTACGT", "ACG"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FindBestMatch("", "ACG"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FindBestMatch("ACGT", null!));
    }

    [Test]
    public void FindBestMatch_Binding_InvokesSuccessfully()
    {
        // Exact match short-circuits at leftmost window (Hamming distance 0).
        var exact = AlignmentTools.FindBestMatch("ACGTACGT", "ACG");
        Assert.Multiple(() =>
        {
            Assert.That(exact.Match, Is.Not.Null);
            Assert.That(exact.Match!.Position, Is.EqualTo(0));
            Assert.That(exact.Match.MatchedSequence, Is.EqualTo("ACG"));
            Assert.That(exact.Match.Distance, Is.EqualTo(0));
        });

        // No good window: best (leftmost minimum) window "AA" has Hamming distance 2 to "TT".
        var far = AlignmentTools.FindBestMatch("AAAA", "TT");
        Assert.Multiple(() =>
        {
            Assert.That(far.Match, Is.Not.Null);
            Assert.That(far.Match!.Position, Is.EqualTo(0));
            Assert.That(far.Match.Distance, Is.EqualTo(2));
            Assert.That(far.Match.MismatchPositions, Is.EqualTo(new[] { 0, 1 }));
        });
    }
}
