using NUnit.Framework;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for CodonOptimizer.FindRareCodons method.
/// Test Unit: CODON-RARE-001
/// 
/// References:
/// - Shu et al. (2006) - Rare codon inhibition in E. coli
/// - Sharp &amp; Li (1987) - Codon Adaptation Index
/// - Kazusa Codon Usage Database
/// </summary>
[TestFixture]
public class CodonOptimizer_FindRareCodons_Tests
{
    #region MUST Tests - Critical functionality

    /// <summary>
    /// M01: Empty sequence should return empty enumerable.
    /// </summary>
    [Test]
    public void FindRareCodons_EmptySequence_ReturnsEmpty()
    {
        var result = CodonOptimizer.FindRareCodons("", CodonOptimizer.EColiK12).ToList();

        Assert.That(result, Is.Empty);
    }

    /// <summary>
    /// M02: Single rare codon (AGA, freq=0.07) is detected at threshold 0.10.
    /// AGA is a known rare arginine codon in E. coli (Shu et al. 2006).
    /// </summary>
    [Test]
    public void FindRareCodons_SingleRareCodon_Detected()
    {
        // AUGAGA = Met + Arg(AGA rare, freq 0.07)
        string sequence = "AUGAGA";

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1));
            Assert.That(result[0].Codon, Is.EqualTo("AGA"));
            Assert.That(result[0].Position, Is.EqualTo(3));
        });
    }

    /// <summary>
    /// M03: Multiple rare codons (AGA, AGG, CGA) all detected.
    /// All are rare arginine codons in E. coli (Kazusa database).
    /// </summary>
    [Test]
    public void FindRareCodons_MultipleRareCodons_AllDetected()
    {
        // AUGAGAAGGCGA = Met + Arg(AGA) + Arg(AGG) + Arg(CGA)
        // AGA=0.07, AGG=0.04, CGA=0.07 - all < 0.10
        string sequence = "AUGAGAAGGCGA";

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(3));
            Assert.That(result.Select(r => r.Codon), Is.EquivalentTo(new[] { "AGA", "AGG", "CGA" }));
        });
    }

    /// <summary>
    /// M04: Position is reported as nucleotide index (codon_index * 3).
    /// </summary>
    [Test]
    public void FindRareCodons_ReportsNucleotidePosition()
    {
        // CUGAGA = Leu(CUG common) + Arg(AGA rare at codon index 1)
        string sequence = "CUGAGA";

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1));
            Assert.That(result[0].Position, Is.EqualTo(3), "Position should be nucleotide index (1 * 3 = 3)");
            Assert.That(result[0].Codon, Is.EqualTo("AGA"));
        });
    }

    /// <summary>
    /// M05: Codon with frequency just below threshold is detected.
    /// UAG (stop) has freq 0.09, should be detected at threshold 0.10.
    /// </summary>
    [Test]
    public void FindRareCodons_ThresholdBoundary_BelowDetected()
    {
        // UAG has frequency 0.09 in E. coli K12
        string sequence = "AUGUAG"; // Met + Stop(UAG)

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.That(result.Any(r => r.Codon == "UAG"), Is.True,
            "UAG (freq 0.09) should be detected at threshold 0.10");
    }

    /// <summary>
    /// M06: Codon with frequency exactly at threshold is NOT detected.
    /// Threshold comparison is strictly less than.
    /// </summary>
    [Test]
    public void FindRareCodons_ThresholdBoundary_AtThresholdNotDetected()
    {
        // CUA has frequency 0.04 in E. coli K12
        // With threshold 0.04, it should NOT be detected (< not ≤)
        string sequence = "AUGCUA"; // Met + Leu(CUA)

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.04).ToList();

        Assert.That(result.Any(r => r.Codon == "CUA"), Is.False,
            "CUA (freq 0.04) should NOT be detected at threshold 0.04 (strict <)");
    }

    /// <summary>
    /// M07: Sequence with only common codons returns empty.
    /// </summary>
    [Test]
    public void FindRareCodons_OnlyCommonCodons_ReturnsEmpty()
    {
        // AUGCUGUGG = Met(1.0) + Leu(0.47) + Trp(1.0) - all high frequency
        string sequence = "AUGCUGUGG";

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.15).ToList();

        Assert.That(result, Is.Empty);
    }

    /// <summary>
    /// M08: Reported amino acid matches standard genetic code translation.
    /// </summary>
    [Test]
    public void FindRareCodons_AminoAcidTranslation_Correct()
    {
        // AGA codes for Arginine (R)
        // CUA codes for Leucine (L)
        string sequence = "AGACUA";

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.Multiple(() =>
        {
            var aga = result.First(r => r.Codon == "AGA");
            Assert.That(aga.AminoAcid, Is.EqualTo("R"), "AGA should translate to R (Arginine)");

            var cua = result.First(r => r.Codon == "CUA");
            Assert.That(cua.AminoAcid, Is.EqualTo("L"), "CUA should translate to L (Leucine)");
        });
    }

    /// <summary>
    /// M09: Reported frequency matches codon usage table lookup.
    /// </summary>
    [Test]
    public void FindRareCodons_FrequencyValue_MatchesTable()
    {
        string sequence = "AUGAGA"; // AGA freq = 0.07 in E. coli K12

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1));
            Assert.That(result[0].Frequency, Is.EqualTo(0.07).Within(0.001),
                "Frequency should match E. coli K12 table value for AGA");
        });
    }

    /// <summary>
    /// M10: DNA input with T is converted to RNA (U) before processing.
    /// </summary>
    [Test]
    public void FindRareCodons_DnaInput_ConvertedToRna()
    {
        // DNA: ATGAGA (T instead of U)
        string dnaSequence = "ATGAGA";

        var result = CodonOptimizer.FindRareCodons(dnaSequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1));
            Assert.That(result[0].Codon, Is.EqualTo("AGA"), "Should detect AGA after T→U conversion");
        });
    }

    #endregion

    #region SHOULD Tests - Important behavior

    /// <summary>
    /// S01: Default threshold when not specified is 0.15.
    /// </summary>
    [Test]
    public void FindRareCodons_DefaultThreshold_Is015()
    {
        // AUA has frequency 0.11 - should be detected with default 0.15 but not 0.10
        string sequence = "AUGAUA"; // Met + Ile(AUA)

        var resultDefault = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12).ToList();

        Assert.That(resultDefault.Any(r => r.Codon == "AUA"), Is.True,
            "AUA (freq 0.11) should be detected with default threshold 0.15");
    }

    /// <summary>
    /// S02: Incomplete codons (trailing nucleotides) are ignored.
    /// </summary>
    [Test]
    public void FindRareCodons_IncompleteCodon_Ignored()
    {
        // 10 nucleotides = 3 complete codons + 1 trailing nucleotide
        string sequence = "AUGAGAAGAG"; // 10 nt

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        // Should process only AUG, AGA, AGA (positions 0, 3, 6)
        // The trailing 'G' is ignored
        Assert.That(result.All(r => r.Position <= 6), Is.True,
            "Should not report position beyond complete codons");
    }

    /// <summary>
    /// S03: Input is case-insensitive.
    /// </summary>
    [Test]
    public void FindRareCodons_LowercaseInput_Handled()
    {
        string sequence = "augaga"; // lowercase

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.That(result, Has.Count.EqualTo(1), "Should detect rare codon regardless of case");
    }

    /// <summary>
    /// S04: Sequence where all codons are rare returns all positions.
    /// </summary>
    [Test]
    public void FindRareCodons_AllRareCodons_AllDetected()
    {
        // All rare arginine codons: AGA(0.07), AGG(0.04), CGA(0.07)
        string sequence = "AGAAGGCGA";

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(3));
            Assert.That(result[0].Position, Is.EqualTo(0));
            Assert.That(result[1].Position, Is.EqualTo(3));
            Assert.That(result[2].Position, Is.EqualTo(6));
        });
    }

    /// <summary>
    /// S05: Mixed common/rare leucine codons - only rare detected.
    /// CUG (0.47) common, CUA (0.04) rare - per Shu et al. 2006.
    /// </summary>
    [Test]
    public void FindRareCodons_MixedCommonRare_OnlyRareDetected()
    {
        // CUG=0.47 (common), CUA=0.04 (rare)
        string sequence = "CUGCUA";

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1));
            Assert.That(result[0].Codon, Is.EqualTo("CUA"));
            Assert.That(result[0].Position, Is.EqualTo(3));
        });
    }

    #endregion

    #region COULD Tests - Edge cases

    /// <summary>
    /// C01: Threshold of 0 should not report any codons (all have freq > 0 or = 0).
    /// </summary>
    [Test]
    public void FindRareCodons_ThresholdZero_NoneDetected()
    {
        string sequence = "AUGAGA";

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.0).ToList();

        // Threshold 0 means freq < 0, which is impossible for valid frequencies
        // AGA has freq 0.07 which is not < 0
        Assert.That(result.Count(r => r.Codon == "AGA" && r.Frequency > 0), Is.EqualTo(0),
            "No codons with frequency > 0 should be detected at threshold 0");
    }

    /// <summary>
    /// C02: Threshold of 1.0 should report all codons as rare.
    /// </summary>
    [Test]
    public void FindRareCodons_ThresholdOne_AllDetected()
    {
        // Even common codons should be detected at threshold 1.0
        // Exception: codons with freq = 1.0 exactly (AUG, UGG)
        string sequence = "CUGCUC"; // Both common leucine codons

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 1.0).ToList();

        Assert.That(result, Has.Count.EqualTo(2), "All codons with freq < 1.0 should be detected");
    }

    /// <summary>
    /// C03: Different organism tables have different rare codons.
    /// </summary>
    [Test]
    public void FindRareCodons_DifferentOrganism_DifferentRare()
    {
        // AGA in E. coli: 0.07 (rare)
        // AGA in Yeast: 0.48 (common!)
        string sequence = "AUGAGA";

        var ecoliResult = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();
        var yeastResult = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.Yeast, 0.10).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(ecoliResult.Any(r => r.Codon == "AGA"), Is.True,
                "AGA should be rare in E. coli");
            Assert.That(yeastResult.Any(r => r.Codon == "AGA"), Is.False,
                "AGA should NOT be rare in Yeast (freq 0.48)");
        });
    }

    #endregion

    #region Invariant Tests

    /// <summary>
    /// Invariant: All reported positions are multiples of 3.
    /// </summary>
    [Test]
    public void FindRareCodons_Invariant_PositionsAreMultiplesOf3()
    {
        string sequence = "AGAAGGCGAAUACUA"; // Multiple rare codons

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.15).ToList();

        Assert.That(result.All(r => r.Position % 3 == 0), Is.True,
            "All positions must be multiples of 3");
    }

    /// <summary>
    /// Invariant: All reported frequencies are below threshold.
    /// </summary>
    [Test]
    public void FindRareCodons_Invariant_AllFrequenciesBelowThreshold()
    {
        const double threshold = 0.12;
        string sequence = "CUGAGAAGGCUACGA"; // Mix of common and rare

        var result = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, threshold).ToList();

        Assert.That(result.All(r => r.Frequency < threshold), Is.True,
            $"All reported frequencies must be < {threshold}");
    }

    /// <summary>
    /// Invariant: Same input produces same output (deterministic).
    /// </summary>
    [Test]
    public void FindRareCodons_Invariant_Deterministic()
    {
        string sequence = "AUGAGAAGGCUA";

        var result1 = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();
        var result2 = CodonOptimizer.FindRareCodons(sequence, CodonOptimizer.EColiK12, 0.10).ToList();

        Assert.That(result1.SequenceEqual(result2), Is.True, "Results must be deterministic");
    }

    #endregion
}
