using System;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for SEQ-GC-001: GC Content Calculation.
/// 
/// Evidence Sources:
/// - Wikipedia: https://en.wikipedia.org/wiki/GC-content
/// - Biopython: https://biopython.org/docs/latest/api/Bio.SeqUtils.html
/// 
/// Formula: GC% = (G + C) / (A + T + G + C) × 100
/// </summary>
[TestFixture]
public class SequenceExtensions_CalculateGcContent_Tests
{
    #region MUST: Empty Sequence (Evidence: Biopython docs)

    [Test]
    public void CalculateGcContent_EmptySequence_ReturnsZero()
    {
        // Evidence: Biopython "Note that this will return zero for an empty sequence"
        ReadOnlySpan<char> empty = "";

        double result = empty.CalculateGcContent();

        Assert.That(result, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateGcFraction_EmptySequence_ReturnsZero()
    {
        ReadOnlySpan<char> empty = "";

        double result = empty.CalculateGcFraction();

        Assert.That(result, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateGcContentFast_EmptyString_ReturnsZero()
    {
        double result = "".CalculateGcContentFast();

        Assert.That(result, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateGcFractionFast_EmptyString_ReturnsZero()
    {
        double result = "".CalculateGcFractionFast();

        Assert.That(result, Is.EqualTo(0.0));
    }

    #endregion

    #region MUST: All GC Returns 100% (Evidence: Formula derivation)

    [Test]
    public void CalculateGcContent_AllG_Returns100()
    {
        // Formula: 4G / 4 total = 100%
        ReadOnlySpan<char> allG = "GGGG";

        double result = allG.CalculateGcContent();

        Assert.That(result, Is.EqualTo(100.0));
    }

    [Test]
    public void CalculateGcContent_AllC_Returns100()
    {
        ReadOnlySpan<char> allC = "CCCC";

        double result = allC.CalculateGcContent();

        Assert.That(result, Is.EqualTo(100.0));
    }

    [Test]
    public void CalculateGcContent_MixedGC_Returns100()
    {
        ReadOnlySpan<char> mixedGc = "GCGCGC";

        double result = mixedGc.CalculateGcContent();

        Assert.That(result, Is.EqualTo(100.0));
    }

    #endregion

    #region MUST: All AT Returns 0% (Evidence: Formula derivation)

    [Test]
    public void CalculateGcContent_AllA_ReturnsZero()
    {
        ReadOnlySpan<char> allA = "AAAA";

        double result = allA.CalculateGcContent();

        Assert.That(result, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateGcContent_AllT_ReturnsZero()
    {
        ReadOnlySpan<char> allT = "TTTT";

        double result = allT.CalculateGcContent();

        Assert.That(result, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateGcContent_MixedAT_ReturnsZero()
    {
        ReadOnlySpan<char> mixedAt = "ATATAT";

        double result = mixedAt.CalculateGcContent();

        Assert.That(result, Is.EqualTo(0.0));
    }

    #endregion

    #region MUST: Equal ACGT Returns 50% (Evidence: Formula)

    [Test]
    public void CalculateGcContent_EqualACGT_Returns50()
    {
        // Formula: 2 GC / 4 total = 50%
        ReadOnlySpan<char> acgt = "ACGT";

        double result = acgt.CalculateGcContent();

        Assert.That(result, Is.EqualTo(50.0));
    }

    [Test]
    public void CalculateGcFraction_EqualACGT_Returns05()
    {
        ReadOnlySpan<char> acgt = "ACGT";

        double result = acgt.CalculateGcFraction();

        Assert.That(result, Is.EqualTo(0.5));
    }

    #endregion

    #region MUST: Mixed Case Handling (Evidence: Biopython "Copes with mixed case")

    [Test]
    public void CalculateGcContent_LowercaseInput_MatchesUppercase()
    {
        ReadOnlySpan<char> upper = "ACGT";
        ReadOnlySpan<char> lower = "acgt";

        double upperResult = upper.CalculateGcContent();
        double lowerResult = lower.CalculateGcContent();

        Assert.That(lowerResult, Is.EqualTo(upperResult));
    }

    [Test]
    public void CalculateGcContent_MixedCaseInput_MatchesUppercase()
    {
        ReadOnlySpan<char> upper = "GCGCGC";
        ReadOnlySpan<char> mixed = "GcGcGc";

        double upperResult = upper.CalculateGcContent();
        double mixedResult = mixed.CalculateGcContent();

        Assert.That(mixedResult, Is.EqualTo(upperResult));
    }

    [Test]
    public void CalculateGcContentFast_MixedCase_MatchesUppercase()
    {
        double upper = "ACGT".CalculateGcContentFast();
        double lower = "acgt".CalculateGcContentFast();
        double mixed = "AcGt".CalculateGcContentFast();

        Assert.Multiple(() =>
        {
            Assert.That(lower, Is.EqualTo(upper));
            Assert.That(mixed, Is.EqualTo(upper));
        });
    }

    #endregion

    #region MUST: Invariant - Fraction = Percentage / 100 (Evidence: Formula)

    [Test]
    [TestCase("ACGT", 50.0)]
    [TestCase("GCGC", 100.0)]
    [TestCase("ATAT", 0.0)]
    [TestCase("ACGTACGT", 50.0)]
    [TestCase("GGGAAA", 50.0)]
    public void CalculateGcContent_FractionMatchesPercentage(string sequence, double expectedPercentage)
    {
        ReadOnlySpan<char> span = sequence;

        double percentage = span.CalculateGcContent();
        double fraction = span.CalculateGcFraction();

        Assert.Multiple(() =>
        {
            Assert.That(percentage, Is.EqualTo(expectedPercentage).Within(0.0001));
            Assert.That(fraction, Is.EqualTo(expectedPercentage / 100.0).Within(0.000001));
        });
    }

    #endregion

    #region MUST: Single Nucleotide Boundary (Evidence: Formula derivation)

    [Test]
    public void CalculateGcContent_SingleG_Returns100()
    {
        ReadOnlySpan<char> singleG = "G";

        double result = singleG.CalculateGcContent();

        Assert.That(result, Is.EqualTo(100.0));
    }

    [Test]
    public void CalculateGcContent_SingleC_Returns100()
    {
        ReadOnlySpan<char> singleC = "C";

        double result = singleC.CalculateGcContent();

        Assert.That(result, Is.EqualTo(100.0));
    }

    [Test]
    public void CalculateGcContent_SingleA_ReturnsZero()
    {
        ReadOnlySpan<char> singleA = "A";

        double result = singleA.CalculateGcContent();

        Assert.That(result, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateGcContent_SingleT_ReturnsZero()
    {
        ReadOnlySpan<char> singleT = "T";

        double result = singleT.CalculateGcContent();

        Assert.That(result, Is.EqualTo(0.0));
    }

    #endregion

    #region SHOULD: Delegate Methods Match Canonical

    [Test]
    public void CalculateGcContentFast_MatchesSpanVersion()
    {
        const string sequence = "ACGTACGTACGT";
        ReadOnlySpan<char> span = sequence;

        double spanResult = span.CalculateGcContent();
        double fastResult = sequence.CalculateGcContentFast();

        Assert.That(fastResult, Is.EqualTo(spanResult));
    }

    [Test]
    public void CalculateGcFractionFast_MatchesSpanVersion()
    {
        const string sequence = "GCGCATATATAT";
        ReadOnlySpan<char> span = sequence;

        double spanResult = span.CalculateGcFraction();
        double fastResult = sequence.CalculateGcFractionFast();

        Assert.That(fastResult, Is.EqualTo(spanResult));
    }

    #endregion

    #region SHOULD: Sequence Type Wrappers Match Canonical

    [Test]
    public void DnaSequence_GcContent_MatchesCanonical()
    {
        const string sequence = "ACGTACGT";
        var dna = new DnaSequence(sequence);

        double canonicalResult = sequence.CalculateGcContentFast();
        double dnaResult = dna.GcContent();

        Assert.That(dnaResult, Is.EqualTo(canonicalResult));
    }

    [Test]
    public void DnaSequence_GcContentFast_MatchesCanonical()
    {
        const string sequence = "GCGCGC";
        var dna = new DnaSequence(sequence);
        ReadOnlySpan<char> span = sequence;

        double spanResult = span.CalculateGcContent();
        double dnaFastResult = dna.GcContentFast();

        Assert.That(dnaFastResult, Is.EqualTo(spanResult));
    }

    [Test]
    public void RnaSequence_GcContent_MatchesCanonical()
    {
        // RNA uses U instead of T, but GC content calculation is the same
        const string sequence = "GCGCGC";
        var rna = new RnaSequence(sequence);

        double canonicalResult = sequence.CalculateGcContentFast();
        double rnaResult = rna.GcContent();

        Assert.That(rnaResult, Is.EqualTo(canonicalResult));
    }

    #endregion

    #region SHOULD: Accurate Calculation for Various Ratios

    [Test]
    [TestCase("G", 100.0)]
    [TestCase("GC", 100.0)]
    [TestCase("GA", 50.0)]
    [TestCase("GCA", 66.666666666666667)]
    [TestCase("GCAA", 50.0)]
    [TestCase("GCAAA", 40.0)]
    [TestCase("GCAAAA", 33.333333333333333)]
    [TestCase("GCAAAAA", 28.571428571428571)]
    public void CalculateGcContent_VariousRatios_ReturnsCorrectPercentage(string sequence, double expected)
    {
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcContent();

        Assert.That(result, Is.EqualTo(expected).Within(0.0000001));
    }

    #endregion

    #region SHOULD: Long Sequence Accuracy

    [Test]
    public void CalculateGcContent_LongSequence_AccurateResult()
    {
        // Create sequence with exactly 500 G/C and 500 A/T
        string sequence = new string('G', 250) + new string('C', 250) +
                          new string('A', 250) + new string('T', 250);
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcContent();

        Assert.That(result, Is.EqualTo(50.0));
    }

    [Test]
    public void CalculateGcContent_10KSequence_AccurateResult()
    {
        // 7500 GC, 2500 AT = 75%
        string sequence = new string('G', 3750) + new string('C', 3750) +
                          new string('A', 1250) + new string('T', 1250);
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcContent();

        Assert.That(result, Is.EqualTo(75.0));
    }

    #endregion

    #region Invariant: Result Bounds (Evidence: Mathematical constraint)

    [Test]
    [TestCase("A")]
    [TestCase("AAAA")]
    [TestCase("G")]
    [TestCase("GGGG")]
    [TestCase("ACGT")]
    [TestCase("ACGTACGTACGTACGT")]
    public void CalculateGcContent_AnyValidInput_ResultInRange0To100(string sequence)
    {
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcContent();

        Assert.That(result, Is.InRange(0.0, 100.0));
    }

    [Test]
    [TestCase("A")]
    [TestCase("AAAA")]
    [TestCase("G")]
    [TestCase("GGGG")]
    [TestCase("ACGT")]
    [TestCase("ACGTACGTACGTACGT")]
    public void CalculateGcFraction_AnyValidInput_ResultInRange0To1(string sequence)
    {
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcFraction();

        Assert.That(result, Is.InRange(0.0, 1.0));
    }

    #endregion

    #region Biological Reference Values (Evidence: Wikipedia)

    [Test]
    public void CalculateGcContent_SimulatedHumanLike_InExpectedRange()
    {
        // Human genome: 35-60% GC (mean ~41%) - Wikipedia ref 20
        // Simulate ~40% GC content
        string sequence = new string('G', 20) + new string('C', 20) +
                          new string('A', 30) + new string('T', 30);
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcContent();

        Assert.That(result, Is.EqualTo(40.0));
    }

    [Test]
    public void CalculateGcContent_SimulatedHighGC_InExpectedRange()
    {
        // Streptomyces coelicolor: 72% GC - Wikipedia ref 29
        string sequence = new string('G', 36) + new string('C', 36) +
                          new string('A', 14) + new string('T', 14);
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcContent();

        Assert.That(result, Is.EqualTo(72.0));
    }

    [Test]
    public void CalculateGcContent_SimulatedLowGC_InExpectedRange()
    {
        // Plasmodium falciparum: ~20% GC - Wikipedia ref 23
        string sequence = new string('G', 10) + new string('C', 10) +
                          new string('A', 40) + new string('T', 40);
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcContent();

        Assert.That(result, Is.EqualTo(20.0));
    }

    #endregion
}
