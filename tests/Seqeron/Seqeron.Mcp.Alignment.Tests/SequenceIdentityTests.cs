using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class SequenceIdentityTests
{
    [Test]
    public void SequenceIdentity_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.SequenceIdentity("ACGT", "ACGT"));
        // Unequal-length inputs return 0 by contract (no throw).
        Assert.DoesNotThrow(() => AlignmentTools.SequenceIdentity("ACGT", "ACG"));
        Assert.That(AlignmentTools.SequenceIdentity("ACGT", "ACG").Identity, Is.EqualTo(0.0));
    }

    [Test]
    public void SequenceIdentity_Binding_InvokesSuccessfully()
    {
        // 3 of 4 positions match -> 0.75.
        var partial = AlignmentTools.SequenceIdentity("ACGT", "ACGA");
        Assert.That(partial.Identity, Is.EqualTo(0.75).Within(1e-9));

        // Identical -> 1.0 (case-insensitive).
        var full = AlignmentTools.SequenceIdentity("ACGT", "acgt");
        Assert.That(full.Identity, Is.EqualTo(1.0).Within(1e-9));
    }
}
