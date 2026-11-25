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
    // cutadapt 5.2 quality_trim_index on the cutadapt worked example (Phred 42,40,26,27,8,7,11,4,2,3):
    // Phred+33 "KI;<)(,%#$" -> (0,4); the same qualities in Phred+64 "jhZ[GFJCAB" -> (0,4);
    // reversed "$#%,()<;IK" with cutoffs (10,0) -> (6,10), with 3'-only (0,10) -> (0,10).
    [Test]
    public void QualityTrimReads_CutadaptOptions_DelegateToCore()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.QualityTrimReads(
                new[] { new QualityReadInput("ACGTACGTAC", "KI;<)(,%#$") }, 10, 1).Trimmed,
                Is.EqualTo(new[] { "ACGT" }));
            Assert.That(AlignmentTools.QualityTrimReads(
                new[] { new QualityReadInput("ACGTACGTAC", "jhZ[GFJCAB") }, 10, 1, phredOffset: 64).Trimmed,
                Is.EqualTo(new[] { "ACGT" }));
            Assert.That(AlignmentTools.QualityTrimReads(
                new[] { new QualityReadInput("ACGTACGTAC", "$#%,()<;IK") }, 10, 1, cutoffFront: 10, cutoffBack: 0).Trimmed,
                Is.EqualTo(new[] { "GTAC" }), "5' only (cutadapt -q 10,0)");
            // A8-1 default flip: a single cutoff is 3'-only (cutadapt -q 10 = -q 0,10; BWA aln -q 10).
            Assert.That(AlignmentTools.QualityTrimReads(
                new[] { new QualityReadInput("ACGTACGTAC", "$#%,()<;IK") }, 10, 1).Trimmed,
                Is.EqualTo(new[] { "ACGTACGTAC" }), "default keeps the low-quality 5' prefix");
            Assert.That(AlignmentTools.QualityTrimReads(
                new[] { new QualityReadInput("ACGTAC", "IIII##") }, 10, 1).Trimmed,
                Is.EqualTo(new[] { "ACGT" }), "default trims the 3' tail (cutadapt (0,4))");
            Assert.That(AlignmentTools.QualityTrimReads(
                new[] { new QualityReadInput("ACGTACGTAC", "$#%,()<;IK") }, 10, 1, cutoffFront: 10).Trimmed,
                Is.EqualTo(new[] { "GTAC" }), "-q 10,10 trims both ends");
            Assert.That(AlignmentTools.QualityTrimReads(
                new[] { new QualityReadInput("ACGTACGTAC", "$#%,()<;IK") }, 10, 1, cutoffFront: 0).Trimmed,
                Is.EqualTo(new[] { "ACGTACGTAC" }));
            // A character outside the Phred+33 range is malformed FASTQ.
            Assert.Throws<ArgumentException>(() => AlignmentTools.QualityTrimReads(
                new[] { new QualityReadInput("ACGT", "II\u001fI") }, 10, 1));
        });
    }

    [Test]
    public void QualityTrimReads_BwaMinLength_KeepsAtLeast35Bases()
    {
        // A8-2: 30 x Q40 + 20 x Q2 (50 bases), cutoff 20: cutadapt cuts to 30, BWA bwa_trim_read (BWA_MIN_RDLEN 35) to 35.
        var read = new QualityReadInput(new string('A', 50), new string('I', 30) + new string('#', 20));
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.QualityTrimReads(new[] { read }, 20, 1).Trimmed.Single(), Has.Length.EqualTo(30));
            Assert.That(AlignmentTools.QualityTrimReads(new[] { read }, 20, 1, bwaMinLength: 35).Trimmed.Single(), Has.Length.EqualTo(35));
            Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.QualityTrimReads(new[] { read }, 20, 1, bwaMinLength: -1));
        });
    }
}
