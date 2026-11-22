using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class SemiGlobalAlignTests
{
    [Test]
    public void SemiGlobalAlign_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.SemiGlobalAlign("ACGT", "TTACGTTT"));
        // R14-2: any letters are accepted like global_align / local_align (the strict ACGT DnaSequence is gone); only null is refused.
        Assert.DoesNotThrow(() => AlignmentTools.SemiGlobalAlign("ACXT", "TTACGTTT"));
        Assert.DoesNotThrow(() => AlignmentTools.SemiGlobalAlign("ACGT", "TTZZTT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SemiGlobalAlign(null!, "TTACGTTT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SemiGlobalAlign("ACGT", null!));
    }

    [Test]
    public void SemiGlobalAlign_Binding_InvokesSuccessfully()
    {
        // Fitting the query "ACGT" into reference "TTACGTTT": free end gaps in the reference
        // let the query match its exact occurrence (score 4). End gaps carry the reference
        // flanks, so AlignedSequence2 reproduces the full reference.
        var r = AlignmentTools.SemiGlobalAlign("ACGT", "TTACGTTT");
        Assert.Multiple(() =>
        {
            Assert.That(r.Score, Is.EqualTo(4));
            Assert.That(r.AlignmentType, Is.EqualTo("SemiGlobal"));
            Assert.That(r.AlignedSequence1, Is.EqualTo("--ACGT--"));
            Assert.That(r.AlignedSequence2, Is.EqualTo("TTACGTTT"));
        });
    }
}
