using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class GcContentTests
{
    [Test]
    public void GcContent_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.GcContent("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.GcContent(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.GcContent(null!));
    }

    [Test]
    public void GcContent_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.GcContent("ATGCGATCGATCG");
        Assert.That(result.GcContent, Is.InRange(0, 100));
        Assert.That(result.GcCount, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.TotalCount, Is.EqualTo(13));

        // 50% GC content
        var fiftyPercent = SequenceTools.GcContent("ATGC");
        Assert.That(fiftyPercent.GcContent, Is.EqualTo(50));
        Assert.That(fiftyPercent.GcCount, Is.EqualTo(2));

        // 100% GC content
        var hundredPercent = SequenceTools.GcContent("GCGC");
        Assert.That(hundredPercent.GcContent, Is.EqualTo(100));
    }

    /// <summary>
    /// A1-11 / F24: the tool's counts and percentage are the Core canonical
    /// <c>CountGcAndValidNucleotides</c> / <c>CalculateGcContent</c> values (bit-identical), and
    /// totalCount is the valid-nucleotide denominator, not the sequence length.
    /// Reference: Biopython 1.88 <c>gc_fraction(s, "remove")</c> × 100.
    /// </summary>
    [TestCase("ATGC", 50.0, 2, 4)]
    [TestCase("ATGCGATCGATCG", 53.84615384615385, 7, 13)]   // gc_fraction = 0.5384615384615384
    [TestCase("acguACGU", 50.0, 4, 8)]                    // lower case + RNA U
    [TestCase("ACGTNNNN", 50.0, 2, 4)]                    // N removed from the denominator
    [TestCase("ACGT-RY", 50.0, 2, 4)]                     // gaps / R / Y removed
    [TestCase("AAUUGGCCNN", 50.0, 4, 8)]
    [TestCase("NNNN", 0.0, 0, 0)]                         // gc_fraction("NNNN", "remove") = 0
    public void GcContent_MatchesCanonicalAndBiopythonRemove(string seq, double expectedPercent, int gc, int valid)
    {
        var r = SequenceTools.GcContent(seq);
        var (cGc, cValid) = global::Seqeron.Genomics.Core.SequenceExtensions.CountGcAndValidNucleotides(seq.AsSpan());
        Assert.Multiple(() =>
        {
            Assert.That(r.GcContent, Is.EqualTo(expectedPercent).Within(1e-12));
            Assert.That(r.GcCount, Is.EqualTo(gc));
            Assert.That(r.TotalCount, Is.EqualTo(valid));
            Assert.That(r.GcCount, Is.EqualTo(cGc));
            Assert.That(r.TotalCount, Is.EqualTo(cValid));
            Assert.That(r.GcContent, Is.EqualTo(global::Seqeron.Genomics.Core.SequenceExtensions.CalculateGcContentFast(seq)));
        });
    }

    /// <summary>
    /// Documented difference from Biopython: S (G|C) / W (A|T) are excluded like every other
    /// ambiguity code. Biopython 1.88 gc_fraction("ACGTSS", "remove") = 0.6667; canonical = 2/4.
    /// </summary>
    [Test]
    public void GcContent_StrongWeakCodes_ExcludedLikeCanonical()
    {
        var r = SequenceTools.GcContent("ACGTSS");
        Assert.Multiple(() =>
        {
            Assert.That(r.GcContent, Is.EqualTo(50.0));
            Assert.That(r.GcCount, Is.EqualTo(2));
            Assert.That(r.TotalCount, Is.EqualTo(4));
        });
    }
}
