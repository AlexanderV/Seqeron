using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class QualityTrimReadsTests
{
    [Test]
    public void QualityTrimReads_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.QualityTrimReads(
            new[] { new QualityReadInput("ACGT", "IIII") }, 20, 1));
        // Sequence/quality length mismatch is rejected.
        Assert.Throws<ArgumentException>(() => AlignmentTools.QualityTrimReads(
            new[] { new QualityReadInput("ACGT", "II") }, 20, 1));
    }

    [Test]
    public void QualityTrimReads_Binding_InvokesSuccessfully()
    {
        // Phred+33: 'I' = 40 (>= cutoff 20), so a uniformly high-quality read is kept whole.
        var kept = AlignmentTools.QualityTrimReads(
            new[] { new QualityReadInput("ACGTACGTAC", "IIIIIIIIII") }, 20, 5);
        Assert.Multiple(() =>
        {
            Assert.That(kept.Trimmed, Has.Length.EqualTo(1));
            Assert.That(kept.Trimmed[0], Is.EqualTo("ACGTACGTAC"));
        });

        // '!' = Phred 0 (< cutoff 20): the whole read is trimmed away, then dropped for
        // being below minLength.
        var dropped = AlignmentTools.QualityTrimReads(
            new[] { new QualityReadInput("ACGTACGTAC", "!!!!!!!!!!") }, 20, 5);
        Assert.That(dropped.Trimmed, Is.Empty);
    }
}
