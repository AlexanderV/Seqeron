using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class ErrorCorrectReadsTests
{
    [Test]
    public void ErrorCorrectReads_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.ErrorCorrectReads(new[] { "ACGT" }, 4, 2));
        // k-mer size below 1 is rejected by the underlying spectrum method.
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.ErrorCorrectReads(new[] { "ACGT" }, 0, 2));
        // Empty input returns empty output.
        Assert.That(AlignmentTools.ErrorCorrectReads(Array.Empty<string>(), 4, 2).Corrected, Is.Empty);
    }

    [Test]
    public void ErrorCorrectReads_Binding_InvokesSuccessfully()
    {
        // Three copies of ACGTACGT make its 4-mers trusted (freq >= 2); the fourth read
        // "ACGTTCGT" has a single error at the middle base that the k-mer spectrum method
        // corrects back to "ACGTACGT" (Musket/Quake two-sided substitution).
        var r = AlignmentTools.ErrorCorrectReads(
            new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTTCGT" }, 4, 2);
        Assert.Multiple(() =>
        {
            Assert.That(r.Corrected, Has.Length.EqualTo(4));
            Assert.That(r.Corrected[3], Is.EqualTo("ACGTACGT"));
            Assert.That(r.Corrected[0], Is.EqualTo("ACGTACGT"));
        });
    }
}
