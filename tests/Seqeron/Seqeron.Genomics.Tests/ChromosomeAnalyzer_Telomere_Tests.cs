using NUnit.Framework;
using Seqeron.Genomics;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for ChromosomeAnalyzer telomere analysis methods.
/// Test Unit: CHROM-TELO-001
/// 
/// Sources:
/// - Wikipedia (Telomere): https://en.wikipedia.org/wiki/Telomere
/// - Meyne et al. (1989): Conservation of human telomere sequence (TTAGGG)n among vertebrates
/// - Cawthon (2002): Telomere measurement by quantitative PCR (T/S ratio method)
/// </summary>
[TestFixture]
[Category("Chromosome")]
[Category("Telomere")]
public class ChromosomeAnalyzer_Telomere_Tests
{
    #region Constants

    private const string VertebrateTelomereRepeat = "TTAGGG";
    private const string VertebrateTelomereRC = "CCCTAA"; // Reverse complement at 5' end
    private const int RepeatLength = 6;

    #endregion

    #region AnalyzeTelomeres - 3' End Detection

    /// <summary>
    /// Validates that TTAGGG repeats at the 3' end are correctly detected.
    /// Source: Wikipedia - Telomere sequences table (TTAGGG for vertebrates)
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_With3PrimeTelomere_DetectsCorrectly()
    {
        // Arrange: 200 TTAGGG repeats = 1200 bp telomere at 3' end
        const int repeatCount = 200;
        string telomereRepeats = string.Concat(Enumerable.Repeat(VertebrateTelomereRepeat, repeatCount));
        string sequence = new string('A', 1000) + telomereRepeats;
        int expectedLength = repeatCount * RepeatLength;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.True, "Should detect 3' telomere");
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(expectedLength),
                "Length should exactly match repeat count × repeat length");
            Assert.That(result.RepeatPurity3Prime, Is.EqualTo(1.0),
                "Perfect repeats must have purity = 1.0");
            Assert.That(result.Chromosome, Is.EqualTo("chr1"),
                "Chromosome name should be preserved");
        });
    }

    /// <summary>
    /// Validates that telomere length matches the number of repeat units.
    /// Source: Algorithm definition - length = repeats × 6
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_3Prime_LengthMatchesRepeatCount()
    {
        // Arrange: Exact number of repeats for predictable length
        const int repeatCount = 100;
        string telomereRepeats = string.Concat(Enumerable.Repeat(VertebrateTelomereRepeat, repeatCount));
        string sequence = new string('A', 1000) + telomereRepeats;
        int expectedLength = repeatCount * RepeatLength;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 100);

        // Assert: Perfect repeats give exact length
        Assert.That(result.TelomereLength3Prime, Is.EqualTo(expectedLength),
            $"Length should be exactly {expectedLength}");
    }

    #endregion

    #region AnalyzeTelomeres - 5' End Detection

    /// <summary>
    /// Validates that CCCTAA (reverse complement) repeats at the 5' end are detected.
    /// Source: Wikipedia - Telomere structure (5' end has RC of telomere repeat)
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_With5PrimeTelomere_DetectsCorrectly()
    {
        // Arrange: CCCTAA repeats at 5' end
        const int repeatCount = 200;
        string telomereRepeats = string.Concat(Enumerable.Repeat(VertebrateTelomereRC, repeatCount));
        string sequence = telomereRepeats + new string('A', 1000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has5PrimeTelomere, Is.True, "Should detect 5' telomere");
            Assert.That(result.TelomereLength5Prime, Is.EqualTo(repeatCount * RepeatLength),
                "Length should exactly match repeat count × repeat length");
            Assert.That(result.RepeatPurity5Prime, Is.EqualTo(1.0),
                "Perfect repeats must have purity = 1.0");
        });
    }

    #endregion

    #region AnalyzeTelomeres - Both Ends

    /// <summary>
    /// Validates that telomeres at both ends are detected independently.
    /// Source: Wikipedia - Chromosome structure has telomeres at both ends
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_BothEnds_BothDetectedIndependently()
    {
        // Arrange: Telomeres at both ends
        const int repeatCount = 150;
        string telomere5 = string.Concat(Enumerable.Repeat(VertebrateTelomereRC, repeatCount));
        string telomere3 = string.Concat(Enumerable.Repeat(VertebrateTelomereRepeat, repeatCount));
        string sequence = telomere5 + new string('A', 2000) + telomere3;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has5PrimeTelomere, Is.True, "Should detect 5' telomere");
            Assert.That(result.Has3PrimeTelomere, Is.True, "Should detect 3' telomere");
            Assert.That(result.TelomereLength5Prime, Is.EqualTo(repeatCount * RepeatLength),
                "5' length should exactly match repeat count × repeat length");
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(repeatCount * RepeatLength),
                "3' length should exactly match repeat count × repeat length");
        });
    }

    #endregion

    #region AnalyzeTelomeres - Edge Cases

    /// <summary>
    /// Validates handling of empty sequence.
    /// Source: Edge case - empty input should not crash
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_EmptySequence_ReturnsNoTelomereCriticallyShort()
    {
        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chr1", "");

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has5PrimeTelomere, Is.False);
            Assert.That(result.Has3PrimeTelomere, Is.False);
            Assert.That(result.TelomereLength5Prime, Is.EqualTo(0));
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(0));
            Assert.That(result.IsCriticallyShort, Is.True,
                "Empty sequence should be marked critically short");
        });
    }

    /// <summary>
    /// Validates handling of sequence without telomeric repeats.
    /// Source: Edge case - random sequence
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_NoRepeats_ReturnsNoTelomere()
    {
        // Arrange: Random non-telomeric sequence
        string sequence = new string('A', 1000);

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 500);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has5PrimeTelomere, Is.False);
            Assert.That(result.Has3PrimeTelomere, Is.False);
            Assert.That(result.TelomereLength5Prime, Is.EqualTo(0),
                "Non-telomeric sequence must have 0 length at 5' end");
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(0),
                "Non-telomeric sequence must have 0 length at 3' end");
        });
    }

    /// <summary>
    /// Validates that telomere below minimum threshold is not reported as present.
    /// Source: API contract - minTelomereLength parameter
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_BelowMinThreshold_HasTelomereFalse()
    {
        // Arrange: 50 repeats = 300 bp, threshold = 500 bp
        const int repeatCount = 50;
        string telomereRepeats = string.Concat(Enumerable.Repeat(VertebrateTelomereRepeat, repeatCount));
        string sequence = new string('A', 1000) + telomereRepeats;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 500);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.False,
                "Telomere below threshold should report false");
            Assert.That(result.TelomereLength3Prime, Is.GreaterThan(0),
                "But length should still be measured");
        });
    }

    /// <summary>
    /// Validates that very short sequence is handled gracefully.
    /// Source: Edge case - robustness
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_VeryShortSequence_HandledGracefully()
    {
        // Arrange: Sequence shorter than repeat length
        string sequence = "ATG";

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chr1", sequence);

        // Assert: Should not throw, should return empty result
        Assert.Multiple(() =>
        {
            Assert.That(result.Has5PrimeTelomere, Is.False);
            Assert.That(result.Has3PrimeTelomere, Is.False);
        });
    }

    #endregion

    #region AnalyzeTelomeres - Critical Length

    /// <summary>
    /// Validates critically short telomere detection.
    /// Source: Wikipedia - Critical telomere length ~3000 bp triggers senescence
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_CriticallyShort_FlaggedCorrectly()
    {
        // Arrange: 100 repeats = 600 bp, critical threshold = 1000 bp
        const int repeatCount = 100;
        string telomereRepeats = string.Concat(Enumerable.Repeat(VertebrateTelomereRepeat, repeatCount));
        string sequence = new string('A', 1000) + telomereRepeats;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 100, criticalLength: 1000);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.True);
            Assert.That(result.IsCriticallyShort, Is.True,
                "Telomere below critical threshold should be flagged");
        });
    }

    /// <summary>
    /// Validates that long telomere is not marked as critically short.
    /// Source: Inverse of critical length test
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_LongTelomere_NotCriticallyShort()
    {
        // Arrange: 500 repeats = 3000 bp, critical threshold = 2000 bp
        const int repeatCount = 500;
        string telomereRepeats = string.Concat(Enumerable.Repeat(VertebrateTelomereRepeat, repeatCount));
        string sequence = new string('A', 1000) + telomereRepeats;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 100, criticalLength: 2000);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.True);
            Assert.That(result.IsCriticallyShort, Is.False,
                "Long telomere should not be critically short");
        });
    }

    #endregion

    #region AnalyzeTelomeres - Case Sensitivity

    /// <summary>
    /// Validates that lowercase sequence is handled correctly.
    /// Source: Robustness requirement
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_LowercaseSequence_DetectedCorrectly()
    {
        // Arrange: Lowercase telomere repeats
        const int repeatCount = 200;
        string telomereRepeats = string.Concat(Enumerable.Repeat("ttaggg", repeatCount));
        string sequence = new string('a', 1000) + telomereRepeats;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.True,
                "Should detect lowercase telomere repeats");
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(repeatCount * RepeatLength),
                "Length must match uppercase result exactly");
            Assert.That(result.RepeatPurity3Prime, Is.EqualTo(1.0),
                "Purity must be 1.0 for perfect repeats regardless of case");
        });
    }

    /// <summary>
    /// Validates that mixed case sequence works.
    /// Source: Robustness requirement
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_MixedCase_DetectedCorrectly()
    {
        // Arrange: Mixed case repeats
        const int repeatCount = 200;
        string telomereRepeats = string.Concat(Enumerable.Repeat("TtAgGg", repeatCount));
        string sequence = new string('A', 1000) + telomereRepeats;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.True,
                "Should detect mixed-case telomere repeats");
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(repeatCount * RepeatLength),
                "Length must match canonical case exactly");
            Assert.That(result.RepeatPurity3Prime, Is.EqualTo(1.0),
                "Purity must be 1.0 for perfect repeats regardless of case");
        });
    }

    #endregion

    #region AnalyzeTelomeres - Custom Repeat

    /// <summary>
    /// Validates detection with non-default repeat (Arabidopsis).
    /// Source: Wikipedia - Species variation table (TTTAGGG for Arabidopsis)
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_CustomRepeat_Detected()
    {
        // Arrange: Arabidopsis repeat TTTAGGG (7 bp)
        const string arabidopsisRepeat = "TTTAGGG";
        const int repeatCount = 150;
        string telomereRepeats = string.Concat(Enumerable.Repeat(arabidopsisRepeat, repeatCount));
        string sequence = new string('A', 1000) + telomereRepeats;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, telomereRepeat: arabidopsisRepeat, minTelomereLength: 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.True,
                "Should detect custom Arabidopsis repeat");
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(repeatCount * arabidopsisRepeat.Length),
                "Length should be repeatCount × 7 for Arabidopsis repeat");
            Assert.That(result.RepeatPurity3Prime, Is.EqualTo(1.0),
                "Perfect repeats must have purity = 1.0");
        });
    }

    #endregion

    #region AnalyzeTelomeres - Invariants

    /// <summary>
    /// Property-based test for telomere result invariants.
    /// Source: Algorithm definition
    /// </summary>
    [Test]
    [TestCase(0, Description = "Zero repeats")]
    [TestCase(10, Description = "Few repeats")]
    [TestCase(100, Description = "Normal repeats")]
    [TestCase(500, Description = "Many repeats")]
    public void AnalyzeTelomeres_ResultInvariants_AlwaysHold(int repeatCount)
    {
        // Arrange
        string telomereRepeats = repeatCount > 0
            ? string.Concat(Enumerable.Repeat(VertebrateTelomereRepeat, repeatCount))
            : "";
        string sequence = new string('A', 1000) + telomereRepeats;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres("chr1", sequence, minTelomereLength: 100);

        // Assert invariants
        Assert.Multiple(() =>
        {
            // Length non-negative
            Assert.That(result.TelomereLength5Prime, Is.GreaterThanOrEqualTo(0),
                "5' length must be non-negative");
            Assert.That(result.TelomereLength3Prime, Is.GreaterThanOrEqualTo(0),
                "3' length must be non-negative");

            // Purity in valid range
            Assert.That(result.RepeatPurity5Prime, Is.InRange(0.0, 1.0),
                "5' purity must be in [0, 1]");
            Assert.That(result.RepeatPurity3Prime, Is.InRange(0.0, 1.0),
                "3' purity must be in [0, 1]");

            // HasTelomere consistency with length (both ends)
            if (result.Has5PrimeTelomere)
            {
                Assert.That(result.TelomereLength5Prime, Is.GreaterThanOrEqualTo(100),
                    "Has5Prime implies length >= minTelomereLength");
            }
            if (result.Has3PrimeTelomere)
            {
                Assert.That(result.TelomereLength3Prime, Is.GreaterThanOrEqualTo(100),
                    "Has3Prime implies length >= minTelomereLength");
            }
        });
    }

    #endregion

    #region EstimateTelomereLengthFromTSRatio

    /// <summary>
    /// Validates T/S ratio = 1.0 returns reference length.
    /// Source: Cawthon (2002) - T/S ratio is proportional to telomere length
    /// </summary>
    [Test]
    public void EstimateTelomereLengthFromTSRatio_Ratio1_ReturnsReferenceLength()
    {
        // Act
        double length = ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(
            tsRatio: 1.0, referenceRatio: 1.0, referenceLength: 7000);

        // Assert
        Assert.That(length, Is.EqualTo(7000).Within(0.001));
    }

    /// <summary>
    /// Validates zero T/S ratio returns zero length.
    /// Source: Edge case - mathematical boundary
    /// </summary>
    [Test]
    public void EstimateTelomereLengthFromTSRatio_ZeroRatio_ReturnsZero()
    {
        // Act
        double length = ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(
            tsRatio: 0.0, referenceRatio: 1.0, referenceLength: 7000);

        // Assert
        Assert.That(length, Is.EqualTo(0).Within(0.001));
    }

    /// <summary>
    /// Validates T/S ratio formula: length = refLen × (tsRatio / refRatio).
    /// Source: Cawthon (2002)
    /// </summary>
    [Test]
    [TestCase(1.5, 1.0, 7000, 10500)]
    [TestCase(0.5, 1.0, 7000, 3500)]
    [TestCase(2.0, 1.0, 7000, 14000)]
    [TestCase(1.0, 2.0, 7000, 3500)]
    [TestCase(3.0, 1.5, 6000, 12000)]
    public void EstimateTelomereLengthFromTSRatio_LinearProportionality(
        double tsRatio, double refRatio, double refLength, double expected)
    {
        // Act
        double length = ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(
            tsRatio, refRatio, refLength);

        // Assert
        Assert.That(length, Is.EqualTo(expected).Within(0.001),
            $"Formula: {refLength} × ({tsRatio} / {refRatio}) = {expected}");
    }

    /// <summary>
    /// Property-based test for T/S ratio invariants.
    /// Source: Algorithm definition
    /// </summary>
    [Test]
    [TestCase(0.1)]
    [TestCase(0.5)]
    [TestCase(1.0)]
    [TestCase(2.0)]
    [TestCase(5.0)]
    public void EstimateTelomereLengthFromTSRatio_ResultNonNegative(double tsRatio)
    {
        // Act
        double length = ChromosomeAnalyzer.EstimateTelomereLengthFromTSRatio(
            tsRatio, referenceRatio: 1.0, referenceLength: 7000);

        // Assert
        Assert.That(length, Is.GreaterThanOrEqualTo(0),
            "Estimated length must be non-negative");
    }

    #endregion

    #region AnalyzeTelomeres - Divergent Repeats

    /// <summary>
    /// Validates that imperfect repeats reduce purity below 1.0.
    /// Source: Evidence doc - biological telomeres show some divergence;
    ///   70% per-window threshold → 1 mismatch per 6bp allowed → purity = 5/6.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_DivergentRepeats_LowerPurity()
    {
        // Arrange: TTAGGA differs from TTAGGG in last base → 5/6 = 83.3% similarity per window
        const int repeatCount = 200;
        string divergentRepeats = string.Concat(Enumerable.Repeat("TTAGGA", repeatCount));
        string sequence = new string('A', 1000) + divergentRepeats;

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, minTelomereLength: 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.True,
                "Divergent repeats above 70% threshold should still be detected");
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(repeatCount * RepeatLength),
                "All windows pass threshold, so full length should be measured");
            Assert.That(result.RepeatPurity3Prime, Is.EqualTo(5.0 / 6.0).Within(0.0001),
                "Purity = matchingBases/totalBases = (200×5)/(200×6) = 5/6");
        });
    }

    #endregion

    #region AnalyzeTelomeres - Long Telomere

    /// <summary>
    /// Validates that very long telomere (>10 kb) is fully measured.
    /// Source: Boundary test - human telomeres can be many kilobases.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_LongTelomere_FullyMeasured()
    {
        // Arrange: 2000 repeats = 12000 bp (>10 kb), require searchLength > 12000
        const int repeatCount = 2000;
        string telomereRepeats = string.Concat(Enumerable.Repeat(VertebrateTelomereRepeat, repeatCount));
        string sequence = new string('A', 1000) + telomereRepeats;
        int expectedLength = repeatCount * RepeatLength; // 12000

        // Act
        var result = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, searchLength: 15000, minTelomereLength: 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.True);
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(expectedLength),
                "Full 12000 bp telomere should be measured");
            Assert.That(result.RepeatPurity3Prime, Is.EqualTo(1.0),
                "Perfect repeats must have purity = 1.0");
        });
    }

    #endregion

    #region AnalyzeTelomeres - SearchLength Limits

    /// <summary>
    /// Validates that searchLength parameter limits the telomere search region.
    /// Source: API contract - only search within searchLength from chromosome ends.
    /// </summary>
    [Test]
    public void AnalyzeTelomeres_SearchLength_LimitsDetection()
    {
        // Arrange: 200 repeats = 1200 bp telomere at 3' end
        const int repeatCount = 200;
        string telomereRepeats = string.Concat(Enumerable.Repeat(VertebrateTelomereRepeat, repeatCount));
        string sequence = new string('A', 1000) + telomereRepeats;
        int fullLength = repeatCount * RepeatLength; // 1200

        // Act with small searchLength (600 bp, all within telomere region)
        var resultLimited = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, searchLength: 600, minTelomereLength: 100);

        // Act with large searchLength (covers entire telomere + body)
        var resultFull = ChromosomeAnalyzer.AnalyzeTelomeres(
            "chr1", sequence, searchLength: 2500, minTelomereLength: 100);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(resultLimited.TelomereLength3Prime, Is.EqualTo(600),
                "Limited search should measure only 600 bp (100 repeats within 600 bp window)");
            Assert.That(resultFull.TelomereLength3Prime, Is.EqualTo(fullLength),
                "Full search should measure entire 1200 bp telomere");
        });
    }

    #endregion

    #region Constants Validation

    /// <summary>
    /// Validates the HumanTelomereRepeat constant is correct.
    /// Source: Wikipedia, Meyne et al. (1989)
    /// </summary>
    [Test]
    public void HumanTelomereRepeat_IsCorrectVertebrate()
    {
        // Assert
        Assert.That(ChromosomeAnalyzer.HumanTelomereRepeat, Is.EqualTo("TTAGGG"),
            "Human telomere repeat should be TTAGGG per Meyne et al. (1989)");
    }

    #endregion
}
