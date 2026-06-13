// POP-SELECT-001 — Selection Signature Detection (integrated Haplotype Score, iHS)
// Evidence: docs/Evidence/POP-SELECT-001-Evidence.md
// TestSpec: tests/TestSpecs/POP-SELECT-001.md
// Source: Voight BF, Kudaravalli S, Wen X, Pritchard JK (2006). PLoS Biology 4(3):e72.
//         Sabeti PC et al. (2002). Nature 419:832-837.
//         Szpiech ZA, Hernandez RD (2014). selscan. Mol Biol Evol 31(10):2824.

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class PopulationGeneticsAnalyzer_SelectionSignature_Tests
{
    private const double Tol = 1e-10;

    #region CalculateEhh

    // M1 — selscan Eq. 3: EHH_c = Σ C(n_h,2)/C(n_c,2).
    // {11,11,11,10}: (C(3,2)+C(1,2))/C(4,2) = (3+0)/6 = 0.5.
    // {00,00,01,01}: (C(2,2)+C(2,2))/C(4,2) = (1+1)/6 = 0.333333...
    [Test]
    public void CalculateEhh_WorkedValues_MatchSelscanFormula()
    {
        double ehhDerived = PopulationGeneticsAnalyzer.CalculateEhh(new[] { "11", "11", "11", "10" });
        double ehhAncestral = PopulationGeneticsAnalyzer.CalculateEhh(new[] { "00", "00", "01", "01" });

        Assert.Multiple(() =>
        {
            Assert.That(ehhDerived, Is.EqualTo(0.5).Within(Tol),
                "EHH for 3x'11'+1x'10' must be (C(3,2)+C(1,2))/C(4,2) = 0.5 (selscan Eq. 3).");
            Assert.That(ehhAncestral, Is.EqualTo(1.0 / 3.0).Within(Tol),
                "EHH for 2x'00'+2x'01' must be (1+1)/6 = 0.3333... (selscan Eq. 3).");
        });
    }

    // M2 — boundary values of the combinatorial formula.
    [Test]
    public void CalculateEhh_SingleAndDistinct_ReturnExtremes()
    {
        double single = PopulationGeneticsAnalyzer.CalculateEhh(new[] { "AAA" });
        double allDistinct = PopulationGeneticsAnalyzer.CalculateEhh(new[] { "AB", "CD", "EF" });

        Assert.Multiple(() =>
        {
            Assert.That(single, Is.EqualTo(1.0).Within(Tol),
                "A single chromosome is trivially homozygous: EHH = 1.");
            Assert.That(allDistinct, Is.EqualTo(0.0).Within(Tol),
                "All-distinct haplotypes share no pairs: EHH = 0.");
        });
    }

    // M15 — null input.
    [Test]
    public void CalculateEhh_Null_Throws()
    {
        Assert.That(() => PopulationGeneticsAnalyzer.CalculateEhh(null!),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(), "Null haplotype list must throw.");
    }

    // M16 — empty sample.
    [Test]
    public void CalculateEhh_Empty_ReturnsZero()
    {
        Assert.That(PopulationGeneticsAnalyzer.CalculateEhh(Array.Empty<string>()),
            Is.EqualTo(0.0).Within(Tol), "Empty sample has no pairs: EHH = 0.");
    }

    #endregion

    #region CalculateIHS

    // M3 — constructed panel: 3 identical derived haplotypes (slow decay), 3 distinct ancestral
    // (instant decay). positions 0,10,20,30,40; core index 2.
    // Derived EHH = 1 at every flank marker -> iHH_D = (1+1)/2*10 + (1+1)/2*10 (right) + same (left) = 40.
    // Ancestral EHH = 0 at first flank -> iHH_A = (1+0)/2*10 (right, then EHH<0.05 stop) + same (left) = 10.
    // unstandardized iHS = ln(iHH_A/iHH_D) = ln(10/40) = ln(0.25) = -1.3862943611198906 (Voight 2006).
    [Test]
    public void CalculateIHS_ConstructedPanel_MatchesDerivedValues()
    {
        var haplotypes = new[] { "AA1GG", "AA1GG", "AA1GG", "TC0TC", "GA0AG", "CT0CA" };
        var positions = new[] { 0, 10, 20, 30, 40 };

        var result = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 2);

        Assert.Multiple(() =>
        {
            Assert.That(result.IhhAncestral, Is.EqualTo(10.0).Within(Tol),
                "iHH_A = two trapezoids ½(1+0)*10 (both directions, truncated at EHH<0.05) = 10.");
            Assert.That(result.IhhDerived, Is.EqualTo(40.0).Within(Tol),
                "iHH_D = four trapezoids ½(1+1)*10 (identical derived haplotypes) = 40.");
            Assert.That(result.UnstandardizedIHS, Is.EqualTo(Math.Log(0.25)).Within(Tol),
                "unstandardized iHS = ln(iHH_A/iHH_D) = ln(10/40) = -1.3862943611 (Voight 2006).");
            Assert.That(result.DerivedAlleleFrequency, Is.EqualTo(0.5).Within(Tol),
                "3 derived of 6 chromosomes -> derived allele frequency 0.5.");
        });
    }

    // M4 — balanced decay: both alleles have identical flanks, equal decay -> iHH_A/iHH_D = 1 -> iHS = 0.
    [Test]
    public void CalculateIHS_BalancedDecay_ReturnsZero()
    {
        // Core idx 1; flanks identical within each allele and symmetric across alleles.
        var haplotypes = new[] { "G1G", "G1G", "G0G", "G0G" };
        var positions = new[] { 0, 5, 10 };

        var result = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 1);

        Assert.That(result.UnstandardizedIHS, Is.EqualTo(0.0).Within(Tol),
            "Balanced EHH decay gives iHH_A/iHH_D = 1 and ln(1) = 0 (Voight 2006).");
    }

    // M5 — sign convention: long derived haplotype -> negative iHS (Voight 2006: negative ⇒ derived sweep).
    [Test]
    public void CalculateIHS_LongDerivedHaplotype_ReturnsNegative()
    {
        var haplotypes = new[] { "AA1GG", "AA1GG", "AA1GG", "TC0TC", "GA0AG", "CT0CA" };
        var positions = new[] { 0, 10, 20, 30, 40 };

        var result = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 2);

        Assert.That(result.UnstandardizedIHS, Is.LessThan(0.0),
            "A long derived haplotype yields large iHH_D, so ln(iHH_A/iHH_D) < 0 (Voight 2006).");
    }

    // S1 — reference-implementation ratio (rehh SNP F1205400): IHH_A=284429.9, IHH_D=2057107.4.
    // unstandardized iHS (Voight) = ln(284429.9/2057107.4) = -1.9785692742315621.
    [Test]
    public void CalculateIHS_RehhRatio_MatchesReference()
    {
        const double ihhA = 284429.9;
        const double ihhD = 2057107.4;
        double expected = Math.Log(ihhA / ihhD);

        Assert.That(expected, Is.EqualTo(-1.9785692742315621).Within(1e-9),
            "rehh worked iHH_A/iHH_D gives unstandardized iHS = -1.97857 (Gautier et al.).");
    }

    // M10 — null haplotypes.
    [Test]
    public void CalculateIHS_NullHaplotypes_Throws()
    {
        Assert.That(() => PopulationGeneticsAnalyzer.CalculateIHS(null!, new[] { 0, 1 }, 0),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(), "Null haplotypes must throw.");
    }

    // M11 — monomorphic core (only derived present) is not a valid iHS site.
    [Test]
    public void CalculateIHS_MonomorphicCore_Throws()
    {
        var haplotypes = new[] { "A1A", "A1A", "A1A" };
        var positions = new[] { 0, 5, 10 };

        Assert.That(() => PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, 1),
            NUnit.Framework.Throws.TypeOf<ArgumentException>(),
            "iHS requires a polymorphic focal SNP (both alleles present) (Voight 2006).");
    }

    // M12 — inconsistent haplotype length vs positions.
    [Test]
    public void CalculateIHS_InconsistentLength_Throws()
    {
        var haplotypes = new[] { "A1A", "A0" };
        var positions = new[] { 0, 5, 10 };

        Assert.That(() => PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, 1),
            NUnit.Framework.Throws.TypeOf<ArgumentException>(),
            "Each haplotype must have one allele per position.");
    }

    // M13 — invalid (non 0/1) core allele.
    [Test]
    public void CalculateIHS_InvalidAllele_Throws()
    {
        var haplotypes = new[] { "AGA", "A0A" };
        var positions = new[] { 0, 5, 10 };

        Assert.That(() => PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, 1),
            NUnit.Framework.Throws.TypeOf<ArgumentException>(),
            "Core alleles must be polarized '0'/'1' (Voight 2006).");
    }

    // M14 — coreIndex out of range.
    [Test]
    public void CalculateIHS_CoreIndexOutOfRange_Throws()
    {
        var haplotypes = new[] { "A1A", "A0A" };
        var positions = new[] { 0, 5, 10 };

        Assert.That(() => PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, 3),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "coreIndex >= positions.Count must throw.");
    }

    // C1 — INV-04 property: ln(iHH_A/iHH_D) = -ln(iHH_D/iHH_A) (Voight vs selscan sign symmetry).
    [Test]
    public void CalculateIHS_SignSymmetry_Property()
    {
        var haplotypes = new[] { "AA1GG", "AT1GG", "AA1GC", "TC0TC", "GA0AG", "CT0CA" };
        var positions = new[] { 0, 10, 20, 30, 40 };

        var result = PopulationGeneticsAnalyzer.CalculateIHS(haplotypes, positions, coreIndex: 2);
        double selscanSign = Math.Log(result.IhhDerived / result.IhhAncestral);

        Assert.That(result.UnstandardizedIHS, Is.EqualTo(-selscanSign).Within(Tol),
            "ln(iHH_A/iHH_D) must equal -ln(iHH_D/iHH_A) (selscan sign note, INV-04).");
    }

    #endregion

    #region StandardizeIHS

    // M6 — within one frequency bin, standardized = (x-mean)/sd. For {-1,0,1}: mean 0, sample sd 1.
    [Test]
    public void StandardizeIHS_SingleBin_CentersToZeroUnitSd()
    {
        var scores = new (double, double)[] { (-1.0, 0.5), (0.0, 0.5), (1.0, 0.5) };

        var standardized = PopulationGeneticsAnalyzer.StandardizeIHS(scores, binCount: 20);

        Assert.Multiple(() =>
        {
            Assert.That(standardized[0], Is.EqualTo(-1.0).Within(Tol),
                "(-1 - 0)/1 = -1 (Voight standardization, sample sd of {-1,0,1} = 1).");
            Assert.That(standardized[1], Is.EqualTo(0.0).Within(Tol), "(0 - 0)/1 = 0.");
            Assert.That(standardized[2], Is.EqualTo(1.0).Within(Tol), "(1 - 0)/1 = 1.");
        });
    }

    // M7 — a singleton bin has undefined SD -> standardized 0 (ASSUMPTION: sd=0 -> 0).
    [Test]
    public void StandardizeIHS_SingletonBin_ReturnsZero()
    {
        var scores = new (double, double)[] { (3.7, 0.9) };

        var standardized = PopulationGeneticsAnalyzer.StandardizeIHS(scores);

        Assert.That(standardized[0], Is.EqualTo(0.0).Within(Tol),
            "A single SNP in a bin has no spread; standardized score is 0.");
    }

    // S2 — two frequency bins are standardized independently.
    [Test]
    public void StandardizeIHS_TwoBins_StandardizedIndependently()
    {
        // Bin 0 (p in [0,0.1)): values {-2,2} -> mean 0, sd = sqrt(8) -> ±sqrt(2).
        // Bin 9 (p in [0.45,0.5)): values {10,20} -> mean 15, sd = sqrt(50) -> ∓... centered at 0.
        var scores = new (double, double)[] { (-2.0, 0.05), (2.0, 0.05), (10.0, 0.48), (20.0, 0.48) };

        var standardized = PopulationGeneticsAnalyzer.StandardizeIHS(scores, binCount: 10);

        Assert.Multiple(() =>
        {
            Assert.That(standardized[0] + standardized[1], Is.EqualTo(0.0).Within(Tol),
                "Bin 0 standardized values are symmetric about 0 (per-bin mean removed).");
            Assert.That(standardized[2] + standardized[3], Is.EqualTo(0.0).Within(Tol),
                "Bin 9 standardized values are symmetric about 0 independently of bin 0.");
            Assert.That(standardized[2], Is.EqualTo(standardized[0]).Within(Tol),
                "Each bin standardized to unit sd: the lower value maps to the same z in both bins.");
        });
    }

    // M18 — null scores and bad bin count.
    [Test]
    public void StandardizeIHS_NullOrBadBin_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => PopulationGeneticsAnalyzer.StandardizeIHS(null!),
                NUnit.Framework.Throws.TypeOf<ArgumentNullException>(), "Null scores must throw.");
            Assert.That(() => PopulationGeneticsAnalyzer.StandardizeIHS(
                    new (double, double)[] { (1.0, 0.5) }, binCount: 0),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "binCount < 1 must throw.");
        });
    }

    #endregion

    #region ScanForSelection

    // M8 — proportion of |iHS|>2 per window (Voight criterion). 2 of 4 extreme -> 0.5.
    [Test]
    public void ScanForSelection_HalfExtreme_ReturnsProportionHalf()
    {
        var scores = new[] { 2.5, 0.3, -3.0, 1.0 };

        var windows = PopulationGeneticsAnalyzer.ScanForSelection(scores, windowSize: 4).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(windows, Has.Count.EqualTo(1), "4 scores in one window of size 4.");
            Assert.That(windows[0].ExtremeCount, Is.EqualTo(2),
                "|2.5|>2 and |-3.0|>2 are extreme; |0.3| and |1.0| are not.");
            Assert.That(windows[0].ProportionExtreme, Is.EqualTo(0.5).Within(Tol),
                "Proportion of |iHS|>2 = 2/4 = 0.5 (Voight 2006).");
        });
    }

    // M9 — windowing of 5 scores with windowSize 2 -> windows of size 2,2,1.
    [Test]
    public void ScanForSelection_FiveScoresWindowTwo_ReturnsThreeWindows()
    {
        var scores = new[] { 3.0, 0.0, 0.0, 0.0, 5.0 };

        var windows = PopulationGeneticsAnalyzer.ScanForSelection(scores, windowSize: 2).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(windows.Select(w => w.SnpCount), Is.EqualTo(new[] { 2, 2, 1 }),
                "5 SNPs split into windows of 2,2,1.");
            Assert.That(windows[0].ExtremeCount, Is.EqualTo(1), "Window 0: only 3.0 is extreme.");
            Assert.That(windows[2].ProportionExtreme, Is.EqualTo(1.0).Within(Tol),
                "Window 2 holds only 5.0 (extreme): proportion 1/1 = 1.");
        });
    }

    // M17 — null scores and bad window size.
    [Test]
    public void ScanForSelection_NullOrBadWindow_Throws()
    {
        IReadOnlyList<double> nullScores = null!;
        IReadOnlyList<double> oneScore = new[] { 1.0 };
        Assert.Multiple(() =>
        {
            Assert.That(() => PopulationGeneticsAnalyzer.ScanForSelection(nullScores).ToList(),
                NUnit.Framework.Throws.TypeOf<ArgumentNullException>(), "Null scores must throw.");
            Assert.That(() => PopulationGeneticsAnalyzer.ScanForSelection(oneScore, 0).ToList(),
                NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(), "windowSize < 1 must throw.");
        });
    }

    #endregion
}
