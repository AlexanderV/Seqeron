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

    #endregion

    #region MUST: All GC Returns 100%, All AT Returns 0% (Evidence: Formula derivation)

    [Test]
    [TestCase("GGGG", 100.0, Description = "All G → 100%")]
    [TestCase("CCCC", 100.0, Description = "All C → 100%")]
    [TestCase("GCGCGC", 100.0, Description = "Mixed G+C → 100%")]
    [TestCase("AAAA", 0.0, Description = "All A → 0%")]
    [TestCase("TTTT", 0.0, Description = "All T → 0%")]
    [TestCase("ATATAT", 0.0, Description = "Mixed A+T → 0%")]
    public void CalculateGcContent_HomogeneousSequences_ReturnsExtrema(string sequence, double expected)
    {
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcContent();

        Assert.That(result, Is.EqualTo(expected));
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
    public void CalculateGcContent_LowercaseInput_Returns50()
    {
        // Evidence: Biopython "Copes with mixed case sequences"
        ReadOnlySpan<char> lower = "acgt";

        double result = lower.CalculateGcContent();

        Assert.That(result, Is.EqualTo(50.0));
    }

    [Test]
    public void CalculateGcContent_MixedCaseInput_Returns100()
    {
        ReadOnlySpan<char> mixed = "GcGcGc";

        double result = mixed.CalculateGcContent();

        Assert.That(result, Is.EqualTo(100.0));
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

    // INV-1/INV-2 (result bounds 0-100 / 0-1) verified by property tests in GcContentProperties.cs

    #region Biological Reference Values (Evidence: Wikipedia)

    [Test]
    [TestCase(20, 20, 30, 30, 40.0, Description = "Human-like ~41% GC — Wikipedia ref 20")]
    [TestCase(36, 36, 14, 14, 72.0, Description = "Streptomyces coelicolor 72% GC — Wikipedia ref 29")]
    [TestCase(10, 10, 40, 40, 20.0, Description = "Plasmodium falciparum ~20% GC — Wikipedia ref 23")]
    public void CalculateGcContent_BiologicalReference_ReturnsExactPercentage(
        int gCount, int cCount, int aCount, int tCount, double expected)
    {
        string sequence = new string('G', gCount) + new string('C', cCount) +
                          new string('A', aCount) + new string('T', tCount);
        ReadOnlySpan<char> span = sequence;

        double result = span.CalculateGcContent();

        Assert.That(result, Is.EqualTo(expected));
    }

    #endregion

    #region Non-nucleotide Character Handling (Evidence: Wikipedia formula, Biopython "remove" mode)

    /// <summary>
    /// Wikipedia formula: GC% = (G+C) / (A+T+G+C) × 100
    /// Non-nucleotide characters (N, R, Y, etc.) are excluded from both numerator and denominator.
    /// Biopython gc_fraction default "remove" mode: "will only count GCS and will only include ACTGSWU
    /// when calculating the sequence length."
    /// 
    /// Cross-verified with Biopython:
    ///   gc_fraction("ACTGN", "remove") == 0.50  → 2/4 = 0.50 (N excluded)
    /// Our result: (G+C)/(A+T+G+C) = 2/4 = 50%
    /// </summary>
    [Test]
    public void CalculateGcContent_SequenceWithN_ExcludesNFromDenominator()
    {
        // "ACTGN": valid = A,C,T,G (4 nucleotides), GC = C,G (2) → 50%
        // Biopython: gc_fraction("ACTGN", "remove") = 0.50
        ReadOnlySpan<char> seq = "ACTGN";

        double result = seq.CalculateGcContent();

        Assert.That(result, Is.EqualTo(50.0));
    }

    [Test]
    public void CalculateGcFraction_SequenceWithN_ExcludesNFromDenominator()
    {
        // Biopython: gc_fraction("ACTGN", "remove") = 0.50
        ReadOnlySpan<char> seq = "ACTGN";

        double result = seq.CalculateGcFraction();

        Assert.That(result, Is.EqualTo(0.5));
    }

    [Test]
    public void CalculateGcContent_SequenceWithMultipleN_ExcludesAllN()
    {
        // "CCTGNN": valid = C,C,T,G (4), GC = C,C,G (3) → 75%
        // Biopython: gc_fraction("CCTGNN", "remove") = 0.75
        ReadOnlySpan<char> seq = "CCTGNN";

        double result = seq.CalculateGcContent();

        Assert.That(result, Is.EqualTo(75.0));
    }

    [Test]
    public void CalculateGcContent_OnlyNonNucleotides_ReturnsZero()
    {
        // No valid nucleotides → return 0 (same as empty sequence)
        ReadOnlySpan<char> seq = "NNNNN";

        double result = seq.CalculateGcContent();

        Assert.That(result, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateGcContent_AllGC_WithAmbiguousBases_ExcludesAmbiguous()
    {
        // "GCNN": valid = G,C (2), GC = G,C (2) → 100%
        // Biopython: gc_fraction("GDVV", "remove") = 1.00 (only G is valid, and it's GC)
        ReadOnlySpan<char> seq = "GCNN";

        double result = seq.CalculateGcContent();

        Assert.That(result, Is.EqualTo(100.0));
    }

    [Test]
    public void CalculateGcContent_WithRNAUracil_TreatsUAsValid()
    {
        // RNA uses U instead of T; U is a valid nucleotide in denominator
        // "GCAU": valid = G,C,A,U (4), GC = G,C (2) → 50%
        // Biopython: gc_fraction("GCAU", "remove") = 0.50
        ReadOnlySpan<char> seq = "GCAU";

        double result = seq.CalculateGcContent();

        Assert.That(result, Is.EqualTo(50.0));
    }

    [Test]
    public void CalculateGcFraction_WithRNAUracil_TreatsUAsValid()
    {
        ReadOnlySpan<char> seq = "GCAU";

        double result = seq.CalculateGcFraction();

        Assert.That(result, Is.EqualTo(0.5));
    }

    [Test]
    public void CalculateGcContent_RNASequence_Returns50()
    {
        // Biopython: gc_fraction("GGAUCUUCGGAUCU", "remove") = 0.50
        // Valid: G,G,A,U,C,U,U,C,G,G,A,U,C,U (14), GC: G,G,C,C,G,G,C (7) → 50%
        ReadOnlySpan<char> seq = "GGAUCUUCGGAUCU";

        double result = seq.CalculateGcContent();

        Assert.That(result, Is.EqualTo(50.0));
    }

    [Test]
    public void CalculateGcContent_BiopythonGDVV_Returns100()
    {
        // Biopython: gc_fraction("GDVV", "remove") == 1.00
        // Only G is a valid nucleotide (D, V are ambiguous IUPAC codes).
        // valid = G (1), GC = G (1) → 100%
        ReadOnlySpan<char> seq = "GDVV";

        double result = seq.CalculateGcContent();

        Assert.That(result, Is.EqualTo(100.0));
    }

    [Test]
    public void CalculateGcContent_SingleU_ReturnsZero()
    {
        // U is a valid nucleotide but not GC → 0%
        ReadOnlySpan<char> seq = "U";

        double result = seq.CalculateGcContent();

        Assert.That(result, Is.EqualTo(0.0));
    }

    [Test]
    public void CalculateGcFraction_SequenceWithMultipleN_ExcludesAllN()
    {
        // "CCTGNN": valid = C,C,T,G (4), GC = C,C,G (3) → 0.75
        // Biopython: gc_fraction("CCTGNN", "remove") = 0.75
        ReadOnlySpan<char> seq = "CCTGNN";

        double result = seq.CalculateGcFraction();

        Assert.That(result, Is.EqualTo(0.75));
    }

    #endregion
}
