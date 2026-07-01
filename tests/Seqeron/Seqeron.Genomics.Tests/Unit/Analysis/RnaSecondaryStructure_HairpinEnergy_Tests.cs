// RNA-HAIRPIN-001 — Hairpin Loop and Stem Free-Energy Calculation (Turner 2004)
// Evidence: docs/Evidence/RNA-HAIRPIN-001-Evidence.md
// TestSpec: tests/TestSpecs/RNA-HAIRPIN-001.md
// Source: Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH (2004).
//         Proc. Natl. Acad. Sci. USA 101:7287-7292. doi:10.1073/pnas.0401799101
//         Parameters/worked examples from NNDB Turner 2004 (rna.urmc.rochester.edu/NNDB/turner04).

using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.RnaSecondaryStructure;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class RnaSecondaryStructure_HairpinEnergy_Tests
{
    private const double Tol = 1e-9;

    #region CalculateHairpinLoopEnergy

    // M1 — NNDB Hairpin Example 1: closing A-U, 6-nt loop, first/last loop base A.
    // Loop component = initiation(6) 5.4 + terminal mismatch (A·A on A-U) -0.8 = 4.6.
    // Evidence: hairpin-example-1.html, loop.txt (init6=5.4), tstack.txt (AAAU=-0.8).
    [Test]
    public void CalculateHairpinLoopEnergy_Example1_SixNtLoop_Returns4_6()
    {
        double dg = CalculateHairpinLoopEnergy("AAAAAA", 'A', 'U');
        Assert.That(dg, Is.EqualTo(4.6).Within(Tol),
            "NNDB Example 1: hairpin loop ΔG°37 = init(6) 5.4 + terminal mismatch -0.8 = 4.6 kcal/mol");
    }

    // M2 — NNDB Hairpin Example 2: closing A-U, 5-nt loop, first/last loop base G.
    // Loop = initiation(5) 5.7 + terminal mismatch (G·G on A-U) -0.8 + GG first-mismatch bonus -0.8 = 4.1.
    // Evidence: hairpin-example-2.html, loop.txt (init5=5.7), tstack.txt (AGGU=-0.8), bonus GG=-0.8.
    [Test]
    public void CalculateHairpinLoopEnergy_Example2_GGFirstMismatch_Returns4_1()
    {
        double dg = CalculateHairpinLoopEnergy("GAAAG", 'A', 'U');
        Assert.That(dg, Is.EqualTo(4.1).Within(Tol),
            "NNDB Example 2: ΔG°37 = init(5) 5.7 + terminal mismatch -0.8 + GG bonus -0.8 = 4.1 kcal/mol");
    }

    // M3 — Special triloop CAACG: experimental total overrides the additive model.
    // Evidence: triloop.txt (CAACG = 6.8). Key = closing5' + loop + closing3' = C + AAC + G.
    [Test]
    public void CalculateHairpinLoopEnergy_SpecialTriloop_CAACG_Returns6_8()
    {
        double dg = CalculateHairpinLoopEnergy("AAC", 'C', 'G');
        Assert.That(dg, Is.EqualTo(6.8).Within(Tol),
            "NNDB special triloop CAACG total = 6.8 kcal/mol (overrides initiation+mismatch model)");
    }

    // M4 — Special tetraloop CCUCGG: experimental total overrides the model.
    // Evidence: tloop.txt (CCUCGG = 2.5). Key = C + CUCG + G.
    [Test]
    public void CalculateHairpinLoopEnergy_SpecialTetraloop_CCUCGG_Returns2_5()
    {
        double dg = CalculateHairpinLoopEnergy("CUCG", 'C', 'G');
        Assert.That(dg, Is.EqualTo(2.5).Within(Tol),
            "NNDB special tetraloop CCUCGG total = 2.5 kcal/mol (overrides model)");
    }

    // M5 — 3-nt loop receives NO sequence-dependent first-mismatch term: just initiation(3).
    // Evidence: hairpin.html (3-nt formula), loop.txt (init3=5.4). Closing G-C, loop AAA (not all-C, not special).
    [Test]
    public void CalculateHairpinLoopEnergy_ThreeNtLoop_NoMismatchTerm_ReturnsInit3()
    {
        double dg = CalculateHairpinLoopEnergy("AAA", 'G', 'C');
        Assert.That(dg, Is.EqualTo(5.4).Within(Tol),
            "NNDB: 3-nt hairpin = initiation(3) 5.4 only; no first-mismatch term for 3-nt loops");
    }

    // M6 — all-C 3-nt loop: initiation(3) 5.4 + C3 penalty +1.5 = 6.9.
    // Evidence: hairpin.html (3-nt all-C penalty), hairpin-mismatch-parameters.html (C3 = +1.5).
    [Test]
    public void CalculateHairpinLoopEnergy_AllC_ThreeNtLoop_AddsC3Penalty()
    {
        double dg = CalculateHairpinLoopEnergy("CCC", 'G', 'C');
        Assert.That(dg, Is.EqualTo(6.9).Within(Tol),
            "NNDB: all-C 3-nt hairpin = init(3) 5.4 + C3 penalty 1.5 = 6.9 kcal/mol");
    }

    // M7 — Loops with fewer than 3 nt are prohibited by the nearest-neighbor rules.
    // Evidence: hairpin.html ("nearest neighbor rules prohibit hairpin loops with fewer than 3 nucleotides").
    // The implementation returns a prohibitive energy so an optimizer never selects them.
    [Test]
    public void CalculateHairpinLoopEnergy_LoopShorterThanThree_ReturnsProhibitiveEnergy()
    {
        double dg = CalculateHairpinLoopEnergy("AA", 'G', 'C');
        // The source assigns NO thermodynamic value (loops <3 nt are prohibited); the implementation
        // returns an exact prohibitive sentinel of 100.0 (INV-02). Assert the exact sentinel so the
        // test cannot pass against a wrong-but-still-large value.
        Assert.That(dg, Is.EqualTo(100.0).Within(Tol),
            "NNDB: hairpin loops < 3 nt are prohibited; implementation returns the exact prohibitive sentinel 100.0 (INV-02)");
    }

    // S1 — special-GU closure bonus (-2.2) is applied only when the flag is set AND closing pair is G-U.
    // Evidence: hairpin-mismatch-parameters.html (special GU closure = -2.2); difference isolates the term.
    [Test]
    public void CalculateHairpinLoopEnergy_SpecialGUClosure_AppliesMinus2_2()
    {
        double withFlag = CalculateHairpinLoopEnergy("AAAA", 'G', 'U', specialGUClosure: true);
        double withoutFlag = CalculateHairpinLoopEnergy("AAAA", 'G', 'U', specialGUClosure: false);
        Assert.That(withFlag - withoutFlag, Is.EqualTo(-2.2).Within(Tol),
            "NNDB special GU closure adds exactly -2.2 kcal/mol when the G-U closing pair is preceded by two Gs");
    }

    // S2 — special-GU closure must NOT apply to a U-G closing pair, even with the flag set.
    // Evidence: hairpin.html ("a GU closing pair (not UG)").
    [Test]
    public void CalculateHairpinLoopEnergy_SpecialGUClosure_NotAppliedForUGClosing()
    {
        double withFlag = CalculateHairpinLoopEnergy("AAAA", 'U', 'G', specialGUClosure: true);
        double withoutFlag = CalculateHairpinLoopEnergy("AAAA", 'U', 'G', specialGUClosure: false);
        Assert.That(withFlag, Is.EqualTo(withoutFlag).Within(Tol),
            "NNDB: special GU closure applies to G-U closing pairs only, not U-G; no -2.2 for U-G");
    }

    // S3 — UU/GA first-mismatch bonus (-0.9) is added on top of the terminal mismatch.
    // Closing C-G, loop "UAAU" (first U, last U). Evidence: hairpin-mismatch-parameters.html (UU/GA = -0.9).
    // init(4)=5.6; terminal mismatch CUUG = -1.2 (tstack); UU bonus -0.9 => 5.6 -1.2 -0.9 = 3.5.
    [Test]
    public void CalculateHairpinLoopEnergy_UUFirstMismatch_AddsMinus0_9()
    {
        double dg = CalculateHairpinLoopEnergy("UAAU", 'C', 'G');
        Assert.That(dg, Is.EqualTo(3.5).Within(Tol),
            "NNDB: ΔG°37 = init(4) 5.6 + terminal mismatch CUUG -1.2 + UU first-mismatch bonus -0.9 = 3.5");
    }

    // S4 — all-C loop > 3 nt uses the linear penalty An+B (A=0.3, B=1.6).
    // Closing G-C, loop "CCCC": init(4) 5.6 + terminal mismatch GCCC -0.7 + (0.3*4 + 1.6 = 2.8) = 7.7.
    // Evidence: hairpin-mismatch-parameters.html (A=0.3,B=1.6), tstack.txt (GCCC=-0.7), loop.txt (init4=5.6).
    [Test]
    public void CalculateHairpinLoopEnergy_AllC_FourNtLoop_LinearPenalty()
    {
        double dg = CalculateHairpinLoopEnergy("CCCC", 'G', 'C');
        Assert.That(dg, Is.EqualTo(7.7).Within(Tol),
            "NNDB: all-C 4-nt hairpin = init(4) 5.6 + tm(GCCC) -0.7 + (0.3*4+1.6) 2.8 = 7.7 kcal/mol");
    }

    // C1 — Special hexaloop ACAGUGUU: experimental total overrides the model.
    // Evidence: hexaloop.txt (ACAGUGUU = 1.8). Key = A + CAGUGU + U.
    [Test]
    public void CalculateHairpinLoopEnergy_SpecialHexaloop_ACAGUGUU_Returns1_8()
    {
        double dg = CalculateHairpinLoopEnergy("CAGUGU", 'A', 'U');
        Assert.That(dg, Is.EqualTo(1.8).Within(Tol),
            "NNDB special hexaloop ACAGUGUU total = 1.8 kcal/mol (overrides model)");
    }

    // S6 — n>30 length extrapolation (Jacobson-Stockmayer): init(n>9) = init(9) + 1.75·R·T·ln(n/9).
    // Loop = 40×A, closing G-C. init(40) = 6.4 + 1.75·1.987·310.15/1000·ln(40/9) = 8.01;
    // terminal mismatch GAAC (closing G-C, first/last A) = -1.1 ⇒ 8.01 - 1.1 = 6.91.
    // Evidence: hairpin.html (init(n>9)=init(9)+1.75 RT ln(n/9)); loop.txt (init9=6.4); tstack.txt (GAAC=-1.1).
    [Test]
    public void CalculateHairpinLoopEnergy_LongLoop_UsesLogExtrapolation()
    {
        double dg = CalculateHairpinLoopEnergy(new string('A', 40), 'G', 'C');
        Assert.That(dg, Is.EqualTo(6.91).Within(1e-2),
            "NNDB: init(40)=init(9)+1.75 RT ln(40/9)=8.01; + terminal mismatch GAAC -1.1 = 6.91 kcal/mol");
    }

    // C2 — Determinism (INV-1): identical inputs yield identical output.
    [Test]
    public void CalculateHairpinLoopEnergy_SameInputs_IsDeterministic()
    {
        double a = CalculateHairpinLoopEnergy("GAAAG", 'A', 'U');
        double b = CalculateHairpinLoopEnergy("GAAAG", 'A', 'U');
        Assert.That(a, Is.EqualTo(b).Within(Tol), "Hairpin energy must be deterministic (INV-1)");
    }

    #endregion

    #region CalculateStemEnergy

    // M8 — Stem of NNDB Example 1 helix: pairs C-G, A-U, C-G, A-U (outer -> inner).
    // 3 stacks (CA/GU -2.11, AC/UG -2.24, CA/GU -2.11) + one terminal AU end penalty +0.45 = -6.01.
    // Evidence: wc-parameters.html (stacking + per-AU-end 0.45), hairpin-example-1.html (helix component -6.01).
    [Test]
    public void CalculateStemEnergy_Example1Helix_ReturnsMinus6_01()
    {
        var pairs = new List<BasePair>
        {
            new(0, 13, 'C', 'G', BasePairType.WatsonCrick),
            new(1, 12, 'A', 'U', BasePairType.WatsonCrick),
            new(2, 11, 'C', 'G', BasePairType.WatsonCrick),
            new(3, 10, 'A', 'U', BasePairType.WatsonCrick),
        };
        double dg = CalculateStemEnergy("CACAxxxxxxUGUG", pairs);
        Assert.That(dg, Is.EqualTo(-6.01).Within(Tol),
            "NNDB Example 1 helix: 3 stacks (-2.11,-2.24,-2.11) + one AU-end penalty +0.45 = -6.01 kcal/mol");
    }

    // M9 — Empty base-pair list contributes no stacking energy (P-1 stacks with P=0).
    // Evidence: wc-parameters.html ("For helices of P uninterrupted basepairs, there are P-1 stacks").
    [Test]
    public void CalculateStemEnergy_EmptyBasePairs_ReturnsZero()
    {
        double dg = CalculateStemEnergy("ACGU", new List<BasePair>());
        Assert.That(dg, Is.EqualTo(0.0).Within(Tol),
            "A stem of 0 base pairs has no stacking terms and no terminal pairs; energy = 0");
    }

    // S5 — A helix whose BOTH ends terminate in A-U receives the +0.45 per-AU-end penalty twice.
    // Pairs A-U, G-C, A-U: one stack-step sum + 2*0.45. Stacks: AU->GC "AG/UC"=-2.08, GC->AU "GA/CU"=-2.35.
    // Evidence: wc-parameters.html (AG/UC=-2.08, GA/CU=-2.35, per-AU-end=0.45).
    [Test]
    public void CalculateStemEnergy_BothEndsAU_AddsAuEndPenaltyTwice()
    {
        var pairs = new List<BasePair>
        {
            new(0, 5, 'A', 'U', BasePairType.WatsonCrick),
            new(1, 4, 'G', 'C', BasePairType.WatsonCrick),
            new(2, 3, 'A', 'U', BasePairType.WatsonCrick),
        };
        double dg = CalculateStemEnergy("AGAUCU", pairs);
        // -2.08 (AG/UC) + -2.35 (GA/CU) + 0.45 + 0.45 = -3.53
        Assert.That(dg, Is.EqualTo(-3.53).Within(Tol),
            "NNDB: stacks AG/UC -2.08 + GA/CU -2.35 + two AU-end penalties (+0.45 each) = -3.53 kcal/mol");
    }

    // S7 — A helix terminating in a G-U wobble pair receives the +0.45 "per GU end" penalty.
    // Pairs G-C (outer, WC), G-U (inner, wobble): one stack GG/CU = -1.53; G-C end no penalty,
    // G-U wobble end +0.45 ⇒ -1.53 + 0.45 = -1.08.
    // Evidence: gu-parameters.html ("Per GU end +0.45"; 5'GG3'/3'CU5' = -1.53).
    [Test]
    public void CalculateStemEnergy_WobbleGUEnd_AddsGuEndPenalty()
    {
        var pairs = new List<BasePair>
        {
            new(0, 3, 'G', 'C', BasePairType.WatsonCrick),
            new(1, 2, 'G', 'U', BasePairType.Wobble),
        };
        double dg = CalculateStemEnergy("GGUC", pairs);
        Assert.That(dg, Is.EqualTo(-1.08).Within(Tol),
            "NNDB: stack GG/CU -1.53 + one GU-end penalty +0.45 = -1.08 kcal/mol");
    }

    // S8 — Special GGUC/CUGG 3-stack context: the 3 individual stacks are replaced by -4.12 total.
    // Pairs G-C, G-U, U-G, C-G (outer→inner). Neither terminal pair (G-C / C-G) is AU/GU, so no end penalty.
    // Evidence: gu-parameters.html note b ("5'GGUC3'/3'CUGG5' = -4.12").
    [Test]
    public void CalculateStemEnergy_SpecialGGUC_CUGG_Returns_Minus4_12()
    {
        var pairs = new List<BasePair>
        {
            new(0, 7, 'G', 'C', BasePairType.WatsonCrick),
            new(1, 6, 'G', 'U', BasePairType.Wobble),
            new(2, 5, 'U', 'G', BasePairType.Wobble),
            new(3, 4, 'C', 'G', BasePairType.WatsonCrick),
        };
        double dg = CalculateStemEnergy("GGUCCUGG", pairs);
        Assert.That(dg, Is.EqualTo(-4.12).Within(Tol),
            "NNDB: special 5'GGUC3'/3'CUGG5' context = -4.12 kcal/mol (replaces 3 individual stacks)");
    }

    #endregion
}
