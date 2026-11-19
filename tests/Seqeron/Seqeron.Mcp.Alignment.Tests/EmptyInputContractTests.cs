using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// Audit A2-39 (B09 F82): the pairwise tools DELEGATE to the library, whose contract (F39) accepts an empty sequence - the alignment of the
/// other sequence against all gaps (no top-N / empty co-optimal alignments for the local tools). Before, the tools re-implemented the check
/// ("Sequence cannot be null or empty", error 1001) and contradicted the documented library behaviour. A missing (null) value is still an error.
/// </summary>
[TestFixture]
public class EmptyInputContractTests
{
    [Test]
    public void GlobalFamily_EmptySequence_IsTheAllGapAlignment_LikeTheLibrary()
    {
        var lib = global::Seqeron.Genomics.Alignment.SequenceAligner.GlobalAlign("", "ACGT");
        var tool = AlignmentTools.GlobalAlign("", "ACGT");
        Assert.That(tool.AlignedSequence1, Is.EqualTo("----"));
        Assert.That(tool.AlignedSequence2, Is.EqualTo("ACGT"));
        Assert.That(tool.Score, Is.EqualTo(lib.Score));

        var swapped = AlignmentTools.GlobalAlign("ACGT", "");
        Assert.That(swapped.AlignedSequence1, Is.EqualTo("ACGT"));
        Assert.That(swapped.AlignedSequence2, Is.EqualTo("----"));

        Assert.That(AlignmentTools.GlobalAlignAffine("", "ACGT").AlignedSequence1, Is.EqualTo("----"));
        Assert.That(AlignmentTools.GlobalAlign("", "").AlignedSequence1, Is.EqualTo(""));
        Assert.That(AlignmentTools.GlobalAlign("", "MKV", matrix: "BLOSUM62").AlignedSequence1, Is.EqualTo("---"));
        Assert.That(AlignmentTools.GlobalAlignCoOptimal("", "ACGT").Alignments, Has.Length.EqualTo(1));
    }

    [Test]
    public void LocalFamily_EmptySequence_DelegatesToTheLibrary()
    {
        var local = AlignmentTools.LocalAlign("", "ACGT");
        Assert.That(local.Score, Is.EqualTo(0));
        Assert.That(AlignmentTools.LocalAlignAffine("ACGT", "").Score, Is.EqualTo(0));
        Assert.That(AlignmentTools.LocalAlignTopN("", "ACGT").Count, Is.EqualTo(0));
        Assert.That(AlignmentTools.LocalAlignCoOptimal("ACGT", "").Count, Is.EqualTo(0));
    }

    [Test]
    public void SemiGlobalAndGeneralGap_EmptySequence_DelegateToTheLibrary()
    {
        var semi = AlignmentTools.SemiGlobalAlign("", "ACGT");
        Assert.That(semi.AlignedSequence2.Replace("-", ""), Is.EqualTo("ACGT"));
        Assert.That(AlignmentTools.SemiGlobalAlignAffine("ACGT", "").AlignedSequence1.Replace("-", ""), Is.EqualTo("ACGT"));
        var gg = AlignmentTools.GeneralGapAlign("", "ACGT", new[] { -1.0, -2.0, -3.0, -4.0 });
        Assert.That(gg.AlignedSequence1, Is.EqualTo("----"));
    }

    [Test]
    public void MissingSequence_IsStillRejected()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign(null!, "ACGT"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlign("ACGT", null!));
        Assert.Throws<ArgumentException>(() => AlignmentTools.SemiGlobalAlign(null!, "ACGT", matrix: "EDNAFULL"));
    }
}
