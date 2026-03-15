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
    /// Source: Universal DNA sequence handling convention.
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
    /// Null, empty, or short (&lt;4 bp) sequence returns 0.
    /// Source: Implementation bounds (need at least 2 dinucleotide units = 4 bp).
    /// </summary>
    [TestCase(null, 0)]
    [TestCase("", 0)]
    [TestCase("ACG", 0)]
    public void FindLongestDinucleotideRepeat_InvalidInput_ReturnsZero(string? sequence, int expected)
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat(sequence!);
        Assert.That(result, Is.EqualTo(expected));
    }

    /// <summary>
    /// Sequence with no dinucleotide repeats returns 1.
    /// Source: Primer3 behavior — single dinucleotide = 1 repeat unit.
    /// </summary>
    [Test]
    public void FindLongestDinucleotideRepeat_NoRepeat_ReturnsOne()
    {
        int result = PrimerDesigner.FindLongestDinucleotideRepeat("ACGT");
        Assert.That(result, Is.EqualTo(1));
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
    /// Null, empty, or too-short sequences cannot form hairpin.
    /// Source: Wikipedia Stem-loop (minimum structure = 2×stem + loop).
    /// </summary>
    [TestCase(null)]
    [TestCase("")]
    [TestCase("ACGT")]
    [TestCase("ACGTACGTAC")] // 10 bp < 2×4+3=11 for default minStemLength=4
    public void HasHairpinPotential_InvalidOrTooShort_ReturnsFalse(string? sequence)
    {
        bool result = PrimerDesigner.HasHairpinPotential(sequence!);
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
    /// Custom minLoopLength is respected.
    /// Source: Wikipedia Stem-loop — minimum 3 nt loop is steric constraint.
    /// </summary>
    [Test]
    public void HasHairpinPotential_CustomMinLoopLength_RespectsParameter()
    {
        // GCGCTTTTGCGC: stem=4 (GCGC), loop=4 (TTTT)
        // With minLoopLength=5, the 4-nt loop is insufficient (and 12 < 4*2+5=13)
        bool withLoop5 = PrimerDesigner.HasHairpinPotential("GCGCTTTTGCGC", minStemLength: 4, minLoopLength: 5);
        bool withLoop3 = PrimerDesigner.HasHairpinPotential("GCGCTTTTGCGC", minStemLength: 4, minLoopLength: 3);

        Assert.Multiple(() =>
        {
            Assert.That(withLoop5, Is.False, "Loop=4 < minLoopLength=5 → no hairpin");
            Assert.That(withLoop3, Is.True, "Loop=4 ≥ minLoopLength=3 → hairpin detected");
        });
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
    /// Null or empty primer returns false.
    /// Source: Standard null guard.
    /// </summary>
    [TestCase(null, "ACGT")]
    [TestCase("", "ACGT")]
    [TestCase("ACGT", null)]
    [TestCase("ACGT", "")]
    public void HasPrimerDimer_NullOrEmptyPrimer_ReturnsFalse(string? p1, string? p2)
    {
        bool result = PrimerDesigner.HasPrimerDimer(p1!, p2!);
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

    #endregion

    #region Calculate3PrimeStability Tests

    /// <summary>
    /// Null, empty, or short (&lt;5 bp) sequence returns 0.
    /// Source: Primer3 uses 5-mer standard (PRIMER_MAX_END_STABILITY).
    /// </summary>
    [TestCase(null, 0.0)]
    [TestCase("", 0.0)]
    [TestCase("ACGT", 0.0)]
    public void Calculate3PrimeStability_InvalidInput_ReturnsZero(string? sequence, double expected)
    {
        double result = PrimerDesigner.Calculate3PrimeStability(sequence!);
        Assert.That(result, Is.EqualTo(expected));
    }

    /// <summary>
    /// Exact 5-base input produces correct ΔG.
    /// Source: SantaLucia (1998) + Primer3 Manual.
    /// TACGT = TA(-0.58) + AC(-1.44) + CG(-2.17) + GT(-1.44) + Init(A·T)(+1.03) + Init(A·T)(+1.03) = -3.57
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_Exact5Bases_ProducesCorrectDeltaG()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("TACGT");
        Assert.That(result, Is.EqualTo(-3.57).Within(0.01));
    }

    /// <summary>
    /// GC-rich 3' end is more stable (more negative ΔG) than AT-rich — with exact values.
    /// Source: SantaLucia (1998) + Primer3 Manual.
    /// GCGCG = -6.86, TATAT = -0.86.
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_GcRich_MoreNegativeThanAtRich()
    {
        double gcRich = PrimerDesigner.Calculate3PrimeStability("ACGTGCGCG"); // ends with GCGCG
        double atRich = PrimerDesigner.Calculate3PrimeStability("ACGTATATAT"); // ends with TATAT

        Assert.Multiple(() =>
        {
            Assert.That(gcRich, Is.EqualTo(-6.86).Within(0.01));
            Assert.That(atRich, Is.EqualTo(-0.86).Within(0.01));
            Assert.That(gcRich, Is.LessThan(atRich));
        });
    }

    /// <summary>
    /// Case insensitive: upper and lower case return identical ΔG.
    /// Source: Universal DNA convention. Verified with exact GCGCG = -6.86.
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_MixedCase_ReturnsSameExactValue()
    {
        double upper = PrimerDesigner.Calculate3PrimeStability("ACGTGCGCG");
        double lower = PrimerDesigner.Calculate3PrimeStability("acgtgcgcg");

        Assert.Multiple(() =>
        {
            Assert.That(upper, Is.EqualTo(-6.86).Within(0.01));
            Assert.That(lower, Is.EqualTo(upper));
        });
    }

    /// <summary>
    /// GCGCG (most stable 5mer) produces exact Primer3 reference value.
    /// Source: Primer3 Manual PRIMER_MAX_END_STABILITY with SantaLucia (1998) parameters.
    /// GCGCG = GC(-2.24) + CG(-2.17) + GC(-2.24) + CG(-2.17) + Init(G·C)(+0.98) + Init(G·C)(+0.98) = -6.86
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_MostStable5mer_MatchesPrimer3()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("AAAAAGCGCG");

        // Primer3 Manual: "most stable 5mer duplex = 6.86 kcal/mol (GCGCG)"
        Assert.That(result, Is.EqualTo(-6.86).Within(0.01));
    }

    /// <summary>
    /// TATAT (least stable 5mer) produces exact Primer3 reference value.
    /// Source: Primer3 Manual PRIMER_MAX_END_STABILITY with SantaLucia (1998) parameters.
    /// TATAT = TA(-0.58) + AT(-0.88) + TA(-0.58) + AT(-0.88) + Init(A·T)(+1.03) + Init(A·T)(+1.03) = -0.86
    /// </summary>
    [Test]
    public void Calculate3PrimeStability_LeastStable5mer_MatchesPrimer3()
    {
        double result = PrimerDesigner.Calculate3PrimeStability("GGGGGTATAT");

        // Primer3 Manual: "most labile 5mer duplex = 0.86 kcal/mol (TATAT)"
        Assert.That(result, Is.EqualTo(-0.86).Within(0.01));
    }

    #endregion

    #region Integration Tests

    /// <summary>
    /// Well-designed primer produces expected exact metrics.
    /// Source: Primer3 primer evaluation workflow.
    /// ACGTACGTACGTACGTACGT: homopolymer=1, dinuc repeat=1 (no repeats),
    /// hairpin=true (ACGT pattern has reverse-complement matches),
    /// last 5 = TACGT → ΔG = -3.57 kcal/mol.
    /// </summary>
    [Test]
    public void PrimerStructureAnalysis_WellDesignedPrimer_ExactMetrics()
    {
        const string primer = "ACGTACGTACGTACGTACGT"; // 20 bp

        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindLongestHomopolymer(primer), Is.EqualTo(1));
            Assert.That(PrimerDesigner.FindLongestDinucleotideRepeat(primer), Is.EqualTo(1));
            Assert.That(PrimerDesigner.Calculate3PrimeStability(primer), Is.EqualTo(-3.57).Within(0.01));
        });
    }

    /// <summary>
    /// Problematic primer (20 G's) produces exact metrics showing all issues.
    /// Source: Primer3 PRIMER_MAX_POLY_X, SantaLucia (1998).
    /// GGGGG = GG(-1.84)×4 + Init(G·C)(+0.98)×2 = -5.40 kcal/mol.
    /// </summary>
    [Test]
    public void PrimerStructureAnalysis_ProblematicPrimer_ExactMetrics()
    {
        const string badPrimer = "GGGGGGGGGGGGGGGGGGGG"; // 20 G's

        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindLongestHomopolymer(badPrimer), Is.EqualTo(20));
            Assert.That(PrimerDesigner.Calculate3PrimeStability(badPrimer), Is.EqualTo(-5.40).Within(0.01));
        });
    }

    #endregion
}
