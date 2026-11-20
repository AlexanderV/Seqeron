using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// WP42 step 0 (B09 for B04 / B05 / B06): the secondary public entry points that the reflection fuzz of A5-3 found throwing <see cref="NullReferenceException"/> /
/// <see cref="IndexOutOfRangeException"/> / <see cref="OutOfMemoryException"/> / <see cref="OverflowException"/> in the FINISHED batches' files now keep the library's contract (an <see cref="ArgumentException"/> family
/// exception for an invalid argument, the documented result otherwise). Findings in unfinished batches' files are registered in QUEUE.md, not fixed here.
/// </summary>
[TestFixture]
public class CrossBatchContract_B09_Tests
{
    [Test]
    public void ApproximateMatcher_NullDnaSequence_IsArgumentNull()
    {
        Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.FindWithMismatches((DnaSequence)null!, "ACGT", 1));
        Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.FindWithMismatches((DnaSequence)null!, "ACGT", 1, CancellationToken.None));
        Assert.Throws<ArgumentNullException>(() => ApproximateMatcher.FindWithEdits((DnaSequence)null!, "ACGT", 1));
    }

    [Test]
    public void ApproximateMatcher_KLargerThanTheSequence_HasNoFrequentKmersAndNoScratchBuffer()
    {
        long before = GC.GetAllocatedBytesForCurrentThread();
        Assert.That(ApproximateMatcher.FindFrequentKmersWithMismatches("ACGTACGT", int.MaxValue, 1), Is.Empty);
        Assert.That(ApproximateMatcher.FindFrequentKmersWithMismatchesAndReverseComplements("ACGTACGT", int.MaxValue, 1), Is.Empty);
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(1 << 20));
    }

    [Test]
    public void SequenceComplexity_HugeStepOrWindow_DoesNotOverflowOrAllocateTheWindow()
    {
        const string seq = "ACGTACGTTAGCCGATTAGGCATCGATCGGATCGAAACCCTTTGGG";
        // one window at 0, then the next start (0 + int.MaxValue) is far beyond the sequence: exactly one point
        Assert.That(SequenceComplexity.CalculateWindowedComplexity(seq, 8, int.MaxValue).Count(), Is.EqualTo(1));
        long before = GC.GetAllocatedBytesForCurrentThread();
        Assert.That(SequenceComplexity.FindLowComplexityIntervals(seq, int.MaxValue), Is.Empty);
        Assert.That(SequenceComplexity.MaskLowComplexity(seq, int.MaxValue), Is.EqualTo(seq));
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(8 << 20));
    }

    [Test]
    public void KmerAnalyzer_HugeArguments_AreArgumentExceptions()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => KmerAnalyzer.GenerateAllKmers(int.MaxValue).First());
        Assert.That(KmerAnalyzer.GenerateAllKmers(3).Count(), Is.EqualTo(64));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.MarkovOrderBics("ACGTACGTAC", int.MaxValue));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.SelectMarkovOrder("ACGTACGTAC", int.MaxValue));
        Assert.Throws<ArgumentException>(() => KmerAnalyzer.MarkovOrderBics("ACGTACGTAC", 10));
        Assert.That(KmerAnalyzer.MarkovOrderBics("ACGTACGTAC", 2), Has.Length.EqualTo(3));
    }
}
