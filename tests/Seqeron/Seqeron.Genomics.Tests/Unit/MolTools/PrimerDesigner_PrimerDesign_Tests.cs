using System.Text;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for Primer Pair Design functionality (PRIMER-DESIGN-001).
/// Covers DesignPrimers, EvaluatePrimer, and GeneratePrimerCandidates methods.
/// 
/// Evidence sources:
/// - Wikipedia: Primer (molecular biology) - "typically between 18 and 24 bases in length" (length only;
///   fetched 2026-10-14 — the page states no GC, Tm or pair-Tm window)
/// - Addgene: How to Design a Primer (fetched 2026-10-14) - 18-24 bases, 40-60% G/C, Tm 50-60°C, pairs within 5°C
/// - Primer3 Manual (v2.6.1) / libprimer3.c pr_set_default_global_args_1 - PRIMER_MIN_TM=57, PRIMER_OPT_TM=60,
///   PRIMER_MAX_TM=63, PRIMER_MAX_POLY_X=5 (same defaults in EMBOSS 6.6.0 eprimer3.acd)
/// - MaxLength 25 is a library convention between Primer3 (27) and Addgene (24), not from either source
/// </summary>
[TestFixture]
public class PrimerDesigner_PrimerDesign_Tests
{

    private DnaSequence _standardTemplate = null!;
    private DnaSequence _gcRichTemplate = null!;
    private DnaSequence _atRichTemplate = null!;

    [SetUp]
    public void SetUp()
    {
        // Standard template with non-palindromic primer regions.
        // Repeating units are chosen so that NO 4-base window is a DNA palindrome
        // (i.e., no window where reverse complement equals itself), preventing
        // hairpin detection from rejecting all candidates.
        //
        // Forward unit: GAACTCGT (50% GC, no 4bp palindromes, max homopolymer=2)
        // Reverse unit: TCCGAAGT (50% GC, no 4bp palindromes, different from forward)
        //
        // Tm is the Primer3-default SantaLucia NN Tm (primer3-py calc_tm): the 20-mer
        // GAACTCGTGAACTCGTGAAC = 56.82°C (just below 57), longer windows fall in [57, 63].
        var sb = new StringBuilder();

        // Forward primer region (100bp, ~50% GC)
        while (sb.Length < 100) sb.Append("GAACTCGT");

        // Target region (50bp poly-T, clearly different from primer regions)
        sb.Append(new string('T', 50));

        // Reverse primer region (100bp, ~50% GC, different unit to avoid primer-dimer)
        int revStart = sb.Length;
        while (sb.Length - revStart < 100) sb.Append("TCCGAAGT");

        _standardTemplate = new DnaSequence(sb.ToString());

        // GC-rich template (70% GC)
        _gcRichTemplate = new DnaSequence(
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" + // 50bp
            "CCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCC" + // 50bp
            "ATATATATATATATATATATATATATATATATATATATATATATATATATAT" + // 50bp target
            "GGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGGGGCCCCGG" + // 50bp
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"   // 50bp
        );

        // AT-rich template (70% AT)
        _atRichTemplate = new DnaSequence(
            "ATATATATATATATATATATATATATATATATATATATATATATATATAT" + // 50bp
            "TTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATT" + // 50bp
            "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC" + // 50bp target
            "AAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTTAAAATTTTAA" + // 50bp
            "TATATATATATATATATATATATATATATATATATATATATATATATATAT"   // 50bp
        );
    }

    [Test]
    public void HasPrimerDimer_ComplementaryPrimers_ReturnsTrue()
    {
        // A 3'-3' primer-dimer needs the 3'-terminal bases of primer2 to be the REVERSE COMPLEMENT
        // of primer1's: primer1 ends ATCGATCG, primer2 ends CGATCGAT (= revcomp(ATCGATCG)).
        // primer3-py check_primers (alignment mode): PRIMER_PAIR_0_COMPL_END = 8.0 ("high end compl").
        // The former fixture (primer2 ending GCTAGCTA, the plain reverse) scores 1.0 in Primer3.
        string primer1 = "AACCGGTTAACCATCGATCG";
        string primer2 = "AACCGGTTAACGATCGAT";

        // Act
        bool hasDimer = PrimerDesigner.HasPrimerDimer(primer1, primer2);

        // Assert
        Assert.That(PrimerDesigner.CalculatePrimerDimerEndComplementarity(primer1, primer2), Is.EqualTo(8.0));
        Assert.That(PrimerDesigner.HasPrimerDimer(primer1, "AACCGGTTAAGCTAGCTA"), Is.False);
        Assert.That(hasDimer, Is.True,
            "Primers with fully complementary 3' ends should be detected as primer-dimer prone");
    }



    [Test]
    public void EvaluatePrimer_SelfComplementary_DetectsHairpin()
    {
        // Arrange - primer with clear hairpin potential:
        // GCGC (stem) + AAAA (loop) + GCGC (matching stem)
        // The reverse complement of GCGC is GCGC, forming a stable hairpin
        string primer = "GCGCAAAAGCGCATGCGATC"; // 20bp with hairpin structure

        // Act
        bool hasHairpin = PrimerDesigner.HasHairpinPotential(primer);

        // Assert
        Assert.That(hasHairpin, Is.True,
            "Primer with GCGC..GCGC stem-loop should be detected as hairpin-prone");
    }

    [Test]
    public void HasHairpinPotential_LongSelfComplementary_ReturnsTrue()
    {
        // Arrange - clearly self-complementary sequence
        // ACGT repeated is self-complementary after internal folding
        string primer = "ACGTACGTACGTACGTACGT"; // 20bp

        // Act
        bool hasHairpin = PrimerDesigner.HasHairpinPotential(primer);

        // Assert
        Assert.That(hasHairpin, Is.True,
            "Self-complementary sequence should be detected as hairpin-prone");
    }



    [Test]
    public void Calculate3PrimeStability_GCRich3Prime_MoreNegative()
    {
        // Arrange
        string gcRich3Prime = "ATATATATATATATGCGCGC"; // 3' = GCGCGC (GC-rich)
        string atRich3Prime = "GCGCGCGCGCGCGCATATAT"; // 3' = ATATAT (AT-rich)

        // Act
        double gcStability = PrimerDesigner.Calculate3PrimeStability(gcRich3Prime);
        double atStability = PrimerDesigner.Calculate3PrimeStability(atRich3Prime);

        // Assert - GC-rich 3' is more stable (more negative ΔG)
        Assert.That(gcStability, Is.LessThan(atStability),
            "GC-rich 3' end should have more negative (stable) ΔG");
    }

}
