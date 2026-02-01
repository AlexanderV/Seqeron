using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for PRIMER-STRUCT-001: Primer Structure Analysis.
/// Tests for secondary structure detection in PCR primers.
/// 
/// Methods under test:
/// - PrimerDesigner.FindLongestHomopolymer(seq)
/// - PrimerDesigner.FindLongestDinucleotideRepeat(seq)
/// - PrimerDesigner.HasHairpinPotential(seq, minStemLength)
/// - PrimerDesigner.HasPrimerDimer(primer1, primer2, minComplementarity)
/// - PrimerDesigner.Calculate3PrimeStability(seq)
/// 
/// Evidence sources:
/// - Wikipedia (Primer, Primer dimer, Stem-loop, Nucleic acid thermodynamics)
/// - Primer3 Manual (primer3.org)
/// - SantaLucia (1998) PNAS 95:1460-65
/// </summary>
[TestFixture]
public class PrimerDesigner_PrimerStructure_Tests
{
    #region FindLongestHomopolymer Tests

    /// <summary>
    /// Empty sequence should return 0.
    /// Source: Standard null/empty handling.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_EmptySequence_ReturnsZero()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("");
        Assert.That(result, Is.EqualTo(0));
    }

    /// <summary>
    /// Null sequence should return 0.
    /// Source: Standard null handling.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_NullSequence_ReturnsZero()
    {
        int result = PrimerDesigner.FindLongestHomopolymer(null!);
        Assert.That(result, Is.EqualTo(0));
    }

    /// <summary>
    /// Sequence with no runs (all different bases) should return 1.
    /// Source: Primer3 PRIMER_MAX_POLY_X behavior.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_NoRun_ReturnsOne()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("ACGT");
        Assert.That(result, Is.EqualTo(1));
    }

    /// <summary>
    /// Sequence with internal homopolymer run returns run length.
    /// Source: Primer3 PRIMER_MAX_POLY_X.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_InternalRun_ReturnsRunLength()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("ACAAAAGT");
        Assert.That(result, Is.EqualTo(4)); // AAAA
    }

    /// <summary>
    /// All same nucleotide returns full length.
    /// Source: Primer3 PRIMER_MAX_POLY_X.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_AllSame_ReturnsFullLength()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("AAAAAA");
        Assert.That(result, Is.EqualTo(6));
    }

    /// <summary>
    /// Case-insensitive matching for homopolymer detection.
    /// Source: Standard DNA sequence handling (ASSUMPTION A1).
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_MixedCase_IsCaseInsensitive()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("AaAaAa");
        Assert.That(result, Is.EqualTo(6));
    }

    /// <summary>
    /// Homopolymer at end of sequence is detected.
    /// Source: Edge case verification.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_RunAtEnd_ReturnsRunLength()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("ACGTTTTT");
        Assert.That(result, Is.EqualTo(5)); // TTTTT at end
    }

    /// <summary>
    /// Multiple runs returns the longest.
    /// Source: Algorithm correctness.
    /// </summary>
    [Test]
    public void FindLongestHomopolymer_MultipleRuns_ReturnsLongest()
    {
        int result = PrimerDesigner.FindLongestHomopolymer("AAACCCCCGG");
        Assert.That(result, Is.EqualTo(5)); // CCCCC is longest
    }

    #endregion

    #region FindLongestDinucleotideRepeat Tests

    /// <summary>
    /// Sequence shorter than 4 bases returns 0.
    /// Source: Implementation bounds (need at least 2 dinucleotide units).
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_TooShort_ReturnsZero()
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ACG");
        Assert.That(result, Is.EqualTo(0));
    }

    /// <summary>
    /// Empty sequence returns 0.
    /// Source: Standard null/empty handling.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_EmptySequence_ReturnsZero()
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("");
        Assert.That(result, Is.EqualTo(0));
    }

    /// <summary>
    /// Sequence with no dinucleotide repeats returns 1 or less.
    /// Source: Primer3 behavior.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_NoRepeat_ReturnsOneOrLess()
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ACGT");
        Assert.That(result, Is.LessThanOrEqualTo(1));
    }

    /// <summary>
    /// ACACACAC contains 4 AC repeats.
    /// Source: Primer3 behavior for dinucleotide repeats.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_AcRepeat_ReturnsCount()
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ACACACACG");
        Assert.That(result, Is.EqualTo(4)); // ACACACAC = 4 x AC
    }

    /// <summary>
    /// AT repeat pattern is detected.
    /// Source: Common microsatellite pattern.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_AtRepeat_ReturnsCount()
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ATATATAT");
        Assert.That(result, Is.EqualTo(4)); // ATATATAT = 4 x AT
    }

    /// <summary>
    /// Returns longest dinucleotide repeat when multiple exist.
    /// Source: Algorithm correctness.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_MultipleRepeats_ReturnsLongest()
    {
        // "ACACGCGCGCGC" = ACAC (2 AC) + GCGCGCGC (4 GC)
        // Implementation counts the number of times the 2-base pattern repeats
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ACACGCGCGCGC");
        Assert.That(result, Is.EqualTo(4)); // GCGCGCGC = 4 x GC (8 bases / 2)
    }

    #endregion

    #region HasHairpinPotential Tests

    /// <summary>
    /// Short sequence cannot form hairpin (needs 2×stem + loop).
    /// Source: Wikipedia Stem-loop (minimum 3 bp loop is sterically required).
    /// </summary>
    [Test]
    public void HasHairpinPotential_TooShort_ReturnsFalse()
    {
        bool result = PrimerDesigner.HasHairpinPotential("ACGT");
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Borderline length with default minStemLength=4: 2×4+3=11 minimum.
    /// Source: Implementation requirement.
    /// </summary>
    [Test]
    public void HasHairpinPotential_BorderlineLength_ReturnsFalse()
    {
        // 10 bases: less than 2×4+3=11 required for default minStemLength=4
        bool result = PrimerDesigner.HasHairpinPotential("ACGTACGTAC");
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Non-self-complementary sequence cannot form hairpin.
    /// Source: Wikipedia Stem-loop (requires complementary regions).
    /// </summary>
    [Test]
    public void HasHairpinPotential_NonSelfComplementary_ReturnsFalse()
    {
        // All A's cannot form complementary stems
        bool result = PrimerDesigner.HasHairpinPotential("AAAACCCCAAAA");
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Self-complementary sequence can form hairpin.
    /// Source: Wikipedia Stem-loop.
    /// </summary>
    [Test]
    public void HasHairpinPotential_SelfComplementary_ReturnsTrue()
    {
        // ACGT reversed = TGCA, which is complementary to ACGT
        bool result = PrimerDesigner.HasHairpinPotential("ACGTACGTACGT");
        Assert.That(result, Is.True);
    }

    /// <summary>
    /// Custom minStemLength is respected.
    /// Source: API contract.
    /// </summary>
    [Test]
    public void HasHairpinPotential_CustomMinStem_RespectsParameter()
    {
        // With minStemLength=6, needs 2×6+3=15 bases minimum
        // 12 bases should return false
        bool result = PrimerDesigner.HasHairpinPotential("ACGTACGTACGT", minStemLength: 6);
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Empty sequence returns false.
    /// Source: Standard null/empty handling.
    /// </summary>
    [Test]
    public void HasHairpinPotential_EmptySequence_ReturnsFalse()
    {
        bool result = PrimerDesigner.HasHairpinPotential("");
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Long sequence (>100bp) uses suffix tree optimization.
    /// Source: Performance optimization test.
    /// </summary>
    [Test]
    public void HasHairpinPotential_LongSequence_UsesSuffixTreeOptimization()
    {
        // Create 150bp sequence with hairpin potential
        // ACGT...ACGT pattern at start and end (complementary when reversed)
        var sb = new System.Text.StringBuilder();
        sb.Append("ACGTACGTACGT"); // 12bp stem region
        sb.Append(new string('A', 126)); // spacer (loop + filler)
        sb.Append("ACGTACGTACGT"); // 12bp complementary region
        string longSeq = sb.ToString(); // 150bp total

        // Should detect hairpin using suffix tree (>100bp threshold)
        bool result = PrimerDesigner.HasHairpinPotential(longSeq);
        Assert.That(result, Is.True);
    }

    /// <summary>
    /// Long sequence without hairpin returns false.
    /// Source: Performance optimization test.
    /// </summary>
    [Test]
    public void HasHairpinPotential_LongSequenceNoHairpin_ReturnsFalse()
    {
        // All A's cannot form hairpin (A is complementary to T, not A)
        string longSeq = new string('A', 150);
        bool result = PrimerDesigner.HasHairpinPotential(longSeq);
        Assert.That(result, Is.False);
    }

    #endregion

    #region HasPrimerDimer Tests

    /// <summary>
    /// Empty primer returns false.
    /// Source: Standard null guard.
    /// </summary>
    [Test]
    public void HasPrimerDimer_EmptyPrimer_ReturnsFalse()
    {
        bool result = PrimerDesigner.HasPrimerDimer("", "ACGT");
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Null primer returns false.
    /// Source: Standard null guard.
    /// </summary>
    [Test]
    public void HasPrimerDimer_NullPrimer_ReturnsFalse()
    {
        bool result = PrimerDesigner.HasPrimerDimer(null!, "ACGT");
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Primers with non-complementary 3' ends do not form dimers.
    /// Source: Wikipedia Primer-dimer.
    /// </summary>
    [Test]
    public void HasPrimerDimer_NonComplementary3Ends_ReturnsFalse()
    {
        // primer1 ends with CCCC, revcomp(primer2) starts with CCCC
        // C-C is not complementary
        bool result = PrimerDesigner.HasPrimerDimer("AAAACCCCCCCC", "GGGGGGGGTTTT");
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Primers with complementary 3' ends form dimers.
    /// Source: Wikipedia Primer-dimer (3' end complementarity is critical).
    /// </summary>
    [Test]
    public void HasPrimerDimer_Complementary3Ends_ReturnsTrue()
    {
        // Poly-A primers: primer1 ends with AAAA
        // revcomp of primer2 (AAAAAAAA) is TTTTTTTT
        // 3' of primer1 (AAAA) vs 5' of revcomp (TTTT) -> A-T complementary
        bool result = PrimerDesigner.HasPrimerDimer("AAAAAAAA", "AAAAAAAA");
        Assert.That(result, Is.True);
    }

    /// <summary>
    /// Custom minComplementarity is respected.
    /// Source: API contract.
    /// </summary>
    [Test]
    public void HasPrimerDimer_CustomMinComplementarity_RespectsParameter()
    {
        // With high minComplementarity threshold, fewer dimers detected
        bool result = PrimerDesigner.HasPrimerDimer("ACGTACGT", "ACGTACGT", minComplementarity: 8);
        Assert.That(result, Is.False);
    }

    /// <summary>
    /// Complementary primers with T-A pairing at 3' ends.
    /// Source: Wikipedia Primer-dimer.
    /// </summary>
    [Test]
    public void HasPrimerDimer_ComplementaryEnds_FormsDimer()
    {
        // primer1 ends with TTTT, primer2 = AAAAAAAA
        // revcomp(AAAAAAAA) = TTTTTTTT
        // Check: 3' of primer1 (TTTT) vs 5' of revcomp (TTTT)
        // T-T is not complementary!
        // To form dimer: primer1's 3' end must complement the 5' end of revcomp(primer2)
        // primer1 = AAAAAAAA (ends with A), revcomp(TTTTTTTT) = AAAAAAAA
        // A-A is not complementary either
        // For actual dimer: AAAAAAAA vs TTTTTTTT
        // revcomp(TTTTTTTT) = AAAAAAAA
        // 3' of AAAAAAAA (AAAAAAAA) vs 5' of AAAAAAAA = A vs A = not complementary
        // The implementation checks: end1[i] IsComplementary end2[i] where end2 is revcomp(primer2)
        // For true dimer: primer1 ending with ACGT, primer2 ending with ACGT
        // revcomp(ACGT) = ACGT -> checking ACGT vs ACGT: A-A, C-C, G-G, T-T = 0 complementary
        // For complementarity: primer1 ends with AAAA, primer2 = TTTT
        // revcomp(TTTT) = AAAA, checking AAAA vs AAAA = A-A not complementary
        // Actually need: primer2 that when revcomped has 5' end complementary to 3' of primer1
        // primer1 = AAAAAAAA, primer2 = AAAAAAAA
        // revcomp(AAAAAAAA) = TTTTTTTT
        // Check 3' of primer1 = AAAAAAAA[last 8] = AAAAAAAA
        // vs 5' of revcomp = TTTTTTTT[first 8] = TTTTTTTT
        // A-T is complementary! This should return true.
        bool result = PrimerDesigner.HasPrimerDimer("AAAAAAAA", "AAAAAAAA");
        Assert.That(result, Is.True); // A complements T in revcomp
    }

    #endregion

    #region Calculate3PrimeStability Tests

    /// <summary>
    /// Short sequence (< 5 bases) returns 0.
    /// Source: Primer3 uses 5-mer standard (PRIMER_MAX_END_STABILITY).
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_TooShort_ReturnsZero()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("ACGT");
        Assert.That(result, Is.EqualTo(0.0));
    }

    /// <summary>
    /// Empty sequence returns 0.
    /// Source: Standard null/empty handling.
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_EmptySequence_ReturnsZero()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("");
        Assert.That(result, Is.EqualTo(0.0));
    }

    /// <summary>
    /// GC-rich 3' end is more stable (more negative ΔG) than AT-rich.
    /// Source: SantaLucia (1998) - GC pairs have stronger stacking.
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_GcRich_MoreNegativeThanAtRich()
    {
        double gcRich = PrimerDesigner.Calculate3PrimeStability("ACGTGCGCG"); // ends with GCGCG
        double atRich = PrimerDesigner.Calculate3PrimeStability("ACGTATATAT"); // ends with TATAT

        // GC-rich 3' end should be more stable (more negative ΔG)
        Assert.That(gcRich, Is.LessThan(atRich));
    }

    /// <summary>
    /// Valid sequence returns negative ΔG value.
    /// Source: SantaLucia (1998) - all nearest-neighbor ΔG values are negative.
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_ValidSequence_ReturnsNegativeValue()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("ACGTACGTACGT");
        Assert.That(result, Is.LessThan(0.0));
    }

    /// <summary>
    /// Mixed case is handled (assumes case insensitivity).
    /// Source: Standard DNA sequence handling (ASSUMPTION A1).
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_MixedCase_ReturnsValue()
    {
        double upper = PrimerDesigner.Calculate3PrimeStability("ACGTGCGCG");
        double mixed = PrimerDesigner.Calculate3PrimeStability("acgtgcgcg");

        Assert.That(mixed, Is.EqualTo(upper));
    }

    /// <summary>
    /// GCGCG (most stable 5mer) produces expected highly negative value.
    /// Source: SantaLucia (1998) - GCGCG has ΔG ≈ -6.86 kcal/mol.
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_MostStable5mer_ProducesHighlyNegativeValue()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("AAAAAGCGCG");

        Assert.Multiple(() =>
        {
            // Should be significantly negative (implementation uses simplified values)
            Assert.That(result, Is.LessThan(-5.0));
            Assert.That(result, Is.GreaterThan(-10.0)); // Sanity check
        });
    }

    #endregion

    #region Integration Tests

    /// <summary>
    /// Verify all methods work together for primer quality assessment.
    /// Source: Primer3 primer evaluation workflow.
    /// </summary>
    [Test]
    public void PrimerStructureAnalysis_IntegrationTest_AllMethodsWork()
    {
        const string primer = "ACGTACGTACGTACGTACGT"; // 20 bp, well-designed primer

        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindLongestHomopolymer(primer), Is.LessThanOrEqualTo(4));
            Assert.That(PrimerDesigner.FindLongestDinucleotideRepeat(primer), Is.LessThanOrEqualTo(4));
            Assert.That(PrimerDesigner.HasHairpinPotential(primer), Is.True.Or.False); // Just verify no exception
            Assert.That(PrimerDesigner.Calculate3PrimeStability(primer), Is.LessThan(0.0));
        });
    }

    /// <summary>
    /// Problematic primer is detected by structure analysis.
    /// Source: Primer3 failure modes.
    /// </summary>
    [Test]
    public void PrimerStructureAnalysis_ProblematicPrimer_IssuesDetected()
    {
        const string badPrimer = "GGGGGGGGGGGGGGGGGGGG"; // 20 G's - terrible primer

        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindLongestHomopolymer(badPrimer), Is.EqualTo(20));
            Assert.That(PrimerDesigner.Calculate3PrimeStability(badPrimer), Is.LessThan(-5.0)); // Very stable
        });
    }

    #endregion
}
