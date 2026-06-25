// PRIMER-TM-001 — Self-dimer / hetero-dimer (intermolecular) Tm via thermodynamic alignment (opt-in)
// Evidence: docs/Evidence/PRIMER-TM-001-DIMER-Evidence.md
// TestSpec: tests/TestSpecs/PRIMER-TM-001-DIMER.md
// Sources:
//   SantaLucia J, Hicks D (2004). Annu Rev Biophys 33:415-440 — unified NN ΔH°/ΔS° (Table 1),
//     bimolecular Tm Eq. 3 (x=1 self-comp / x=4 non), Eq. 5 entropy salt correction (0.368).
//   Untergasser A et al. (2012). Nucleic Acids Res 40:e115 — Primer3/ntthal thermodynamic alignment.
//   primer3-py 2.3.0 (calc_homodimer / calc_heterodimer; mv=50, dv=0, dntp=0, dna_conc=50 nM)
//     — reference ΔH/ΔS/Tm captured this session (ntthal engine).

using System;
using NUnit.Framework;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the opt-in self-/hetero-dimer Tm added under PRIMER-TM-001. The MUST cases use
/// values hand-derived from the SantaLucia &amp; Hicks (2004) Table 1 NN parameters
/// (independent of the implementation) and/or captured from primer3-py 2.3.0 (the ntthal
/// reference). Hand-derived ΔH°/ΔS°/Tm are asserted to 1e-9; primer3 parity to 1e-3.
/// </summary>
[TestFixture]
public class PrimerDesigner_DimerTm_Tests
{
    private const double ExactTol = 1e-9;   // hand-derived (matches both derivation and primer3)
    private const double ParityTol = 1e-3;  // primer3-py reference rounded to 4 dp

    // primer3-py default dimer conditions captured this session.
    private const double Na = 0.05;     // 50 mM monovalent
    private const double Ct = 50e-9;    // 50 nM total strand (ntthal dna_conc)

    #region FindMostStableDimer — thermodynamics (hand-derived)

    // M1 — GCGCGCGC self-dimer, full 8-bp Watson-Crick duplex (palindrome → x=1).
    // Hand-derived from SantaLucia & Hicks (2004) Table 1:
    //   stacks 4·GC(-9.8,-24.4) + 3·CG(-10.6,-27.2); init(+0.2,-5.7); no A·T end;
    //   ΔH° = -70.8 kcal/mol; ΔS° = -184.9 + 7·0.368·ln(0.05) = -192.61700633667505;
    //   ΔG°37 = ΔH° - 310.15·ΔS°/1000 = -11.059835484680235.
    [Test]
    public void FindMostStableDimer_GcgcgcgcSelfDimer_MatchesHandDerivedThermodynamics()
    {
        var result = PrimerDesigner.FindMostStableDimer("GCGCGCGC", "GCGCGCGC", Na, Ct);

        Assert.That(result, Is.Not.Null, "GCGCGCGC forms a fully complementary self-dimer.");
        var d = result!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(d.BasePairs, Is.EqualTo(8),
                "The most stable self-dimer is the full 8-bp Watson-Crick duplex.");
            Assert.That(d.DeltaH, Is.EqualTo(-70.8).Within(ExactTol),
                "ΔH° = init(+0.2) + 4·GC(-9.8) + 3·CG(-10.6) = -70.8 kcal/mol (salt-independent).");
            Assert.That(d.DeltaS, Is.EqualTo(-192.61700633667505).Within(ExactTol),
                "ΔS° = -5.7 + 4·(-24.4) + 3·(-27.2) + 7·0.368·ln(0.05) (Eq. 5 salt correction).");
            Assert.That(d.DeltaG37, Is.EqualTo(-11.059835484680235).Within(ExactTol),
                "ΔG°37 = ΔH° - 310.15·ΔS°/1000.");
        });
    }

    // M12 — alignment spans: the optimal self-dimer starts at index 0 on both strands.
    [Test]
    public void FindMostStableDimer_GcgcgcgcSelfDimer_ReportsFullAlignmentSpans()
    {
        var d = PrimerDesigner.FindMostStableDimer("GCGCGCGC", "GCGCGCGC", Na, Ct)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(d.Strand1Start, Is.EqualTo(0), "Duplex covers strand 1 from index 0.");
            Assert.That(d.Strand2Start, Is.EqualTo(0), "Duplex covers strand 2 from index 0.");
        });
    }

    // M3 — TGCATGCATG / CATGCATGCA: non-palindromic hetero-dimer (x=4).
    // Hand-derived: stacks TG,GC,CA,AT,TG,GC,CA,AT,TG; init(+0.2,-5.7); one A·T end (5'-T).
    //   ΔH° = -74.1 kcal/mol; ΔS° = -211.8218652900108; Tm(x=4) = 25.659587124835923 °C.
    [Test]
    public void FindMostStableDimer_NonPalindromicHeteroDimer_MatchesHandDerivedThermodynamics()
    {
        var result = PrimerDesigner.FindMostStableDimer("TGCATGCATG", "CATGCATGCA", Na, Ct);

        Assert.That(result, Is.Not.Null, "These strands form a 10-bp Watson-Crick hetero-dimer.");
        var d = result!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(d.BasePairs, Is.EqualTo(10), "Full 10-bp duplex.");
            Assert.That(d.DeltaH, Is.EqualTo(-74.1).Within(ExactTol),
                "ΔH° = init + Σ stacks + one terminal-A·T penalty (5'-T end).");
            Assert.That(d.DeltaS, Is.EqualTo(-211.8218652900108).Within(ExactTol),
                "ΔS° includes one A·T-penalty (+6.9) and 9·0.368·ln(0.05) salt correction.");
        });
    }

    #endregion

    #region Tm — hand-derived (exact)

    // M2 — GCGCGCGC self-dimer Tm at C_T=50 nM, x=1 (palindrome). Hand-derived 40.09064476882935 °C.
    [Test]
    public void CalculateSelfDimerMeltingTemperature_Gcgcgcgc_MatchesHandDerivedTm()
    {
        double tm = PrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC", Na, Ct);
        Assert.That(tm, Is.EqualTo(40.09064476882935).Within(ExactTol),
            "Tm = ΔH°·1000/(ΔS° + R·ln(C_T/1)) - 273.15 with the hand-derived ΔH°/ΔS° (palindrome → x=1).");
    }

    // M3 — TGCATGCATG/CATGCATGCA hetero-dimer Tm, x=4. Hand-derived 25.659587124835923 °C.
    [Test]
    public void CalculateDimerMeltingTemperature_NonPalindromicHeteroDimer_MatchesHandDerivedTm()
    {
        double tm = PrimerDesigner.CalculateDimerMeltingTemperature("TGCATGCATG", "CATGCATGCA", Na, Ct);
        Assert.That(tm, Is.EqualTo(25.659587124835923).Within(ExactTol),
            "Non-palindromic pair uses x=4 in the bimolecular Tm equation.");
    }

    #endregion

    #region Tm — primer3-py 2.3.0 parity

    // M4..M8 — parity with primer3-py 2.3.0 (ntthal) on pairs whose optimum is a contiguous WC duplex.
    [Test]
    [TestCase("GCGCGCGC", "GCGCGCGC", 40.0906, TestName = "GCGCGCGC self-dimer (x=1)")]
    [TestCase("ACGTACGTACGT", "ACGTACGTACGT", 37.6251, TestName = "ACGTACGTACGT self-dimer (x=1)")]
    [TestCase("ATCGATCGATCG", "CGATCGATCGAT", 32.6107, TestName = "ATCGATCGATCG/CGATCGATCGAT hetero (x=1)")]
    [TestCase("CGATCGATCG", "CGATCGATCG", 29.6600, TestName = "CGATCGATCG self-dimer palindrome (x=1)")]
    [TestCase("GCATGC", "GCATGC", 0.6859, TestName = "GCATGC self-dimer (x=1)")]
    [TestCase("GGGGCCCC", "GGGGCCCC", 29.0150, TestName = "GGGGCCCC dimer (x=1)")]
    [TestCase("TGCATGCATG", "CATGCATGCA", 25.6596, TestName = "TGCATGCATG/CATGCATGCA hetero (x=4)")]
    public void CalculateDimerMeltingTemperature_ContiguousWcDuplex_MatchesPrimer3(
        string a, string b, double primer3Tm)
    {
        double tm = PrimerDesigner.CalculateDimerMeltingTemperature(a, b, Na, Ct);
        Assert.That(tm, Is.EqualTo(primer3Tm).Within(ParityTol),
            $"ntthal (primer3-py 2.3.0) reports Tm={primer3Tm} °C for the {a}/{b} dimer.");
    }

    #endregion

    #region Full ntthal DP — non-contiguous optima (internal loop / bulge / overhang) parity

    // The remaining PRIMER-TM-001 limitation: dimers whose most stable structure is NOT a
    // contiguous Watson–Crick run (internal mismatch / internal loop, bulge, terminal overhang).
    // CalculateDimerThermodynamicsNtthal runs the full Primer3 ntthal DP and must reproduce
    // primer3-py 2.3.0 (mv=50, dv=0, dntp=0, dna_conc=50 nM) ΔH/ΔS/ΔG/Tm. Reference numbers were
    // captured this session directly from primer3-py 2.3.0 calc_homodimer/calc_heterodimer. These
    // are the cases the contiguous-only scorer (FindMostStableDimer) gets wrong.

    // N1 — 2×2 internal loop. GCGCATGCGC self-dimer: the AT centre cannot pair (A·G, T·C),
    // so the optimum is two 4-bp GCGC arms bridged by a symmetric internal loop, NOT a 4-bp run.
    [Test]
    public void NtthalDimer_InternalLoop2x2_MatchesPrimer3()
    {
        var d = PrimerDesigner.CalculateDimerThermodynamicsNtthal("GCGCATGCGC", "GCGCATGCGC", Na, Ct);
        Assert.That(d, Is.Not.Null, "An internal-loop dimer exists.");
        var v = d!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(v.TmCelsius, Is.EqualTo(43.1572).Within(ParityTol), "primer3 ntthal Tm.");
            Assert.That(v.DeltaH, Is.EqualTo(-84.40000).Within(ParityTol), "primer3 ntthal ΔH (kcal/mol).");
            Assert.That(v.DeltaS, Is.EqualTo(-233.42187).Within(ParityTol), "primer3 ntthal ΔS (cal/K/mol).");
            Assert.That(v.DeltaG37, Is.EqualTo(-12.00421).Within(ParityTol), "primer3 ntthal ΔG°37 (kcal/mol).");
        });
    }

    // N2 — 3×3 internal loop. GCGCAAAGCGC/GCGCTTTGCGC: the AAA/TTT centre forms an internal loop
    // (the contiguous scorer would split this into a 4-bp arm and miss the loop contribution).
    [Test]
    public void NtthalDimer_InternalLoop3x3_MatchesPrimer3()
    {
        var v = PrimerDesigner.CalculateDimerThermodynamicsNtthal("GCGCAAAGCGC", "GCGCTTTGCGC", Na, Ct)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(v.TmCelsius, Is.EqualTo(41.8816).Within(ParityTol), "primer3 ntthal Tm.");
            Assert.That(v.DeltaH, Is.EqualTo(-92.30000).Within(ParityTol), "primer3 ntthal ΔH (kcal/mol).");
            Assert.That(v.DeltaS, Is.EqualTo(-256.82429).Within(ParityTol), "primer3 ntthal ΔS (cal/K/mol).");
            Assert.That(v.DeltaG37, Is.EqualTo(-12.64594).Within(ParityTol), "primer3 ntthal ΔG°37 (kcal/mol).");
        });
    }

    // N3 — single-base bulge. GCGCGCGC/GCGCAGCGC: the extra A on strand 2 bulges out, leaving two
    // 4-bp GCGC arms bridged by a 1-base bulge (the contiguous scorer cannot span the bulge).
    [Test]
    public void NtthalDimer_SingleBaseBulge_MatchesPrimer3()
    {
        var v = PrimerDesigner.CalculateDimerThermodynamicsNtthal("GCGCGCGC", "GCGCAGCGC", Na, Ct)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(v.TmCelsius, Is.EqualTo(19.8125).Within(ParityTol), "primer3 ntthal Tm.");
            Assert.That(v.DeltaH, Is.EqualTo(-70.80000).Within(ParityTol), "primer3 ntthal ΔH (kcal/mol).");
            Assert.That(v.DeltaS, Is.EqualTo(-205.50701).Within(ParityTol), "primer3 ntthal ΔS (cal/K/mol).");
            Assert.That(v.DeltaG37, Is.EqualTo(-7.06200).Within(ParityTol), "primer3 ntthal ΔG°37 (kcal/mol).");
        });
    }

    // N4 — 2×2 internal loop with mixed AC/AT centre. GCGCACGCGC/GCGCTAGCGC.
    [Test]
    public void NtthalDimer_InternalLoopMixed_MatchesPrimer3()
    {
        var v = PrimerDesigner.CalculateDimerThermodynamicsNtthal("GCGCACGCGC", "GCGCTAGCGC", Na, Ct)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(v.TmCelsius, Is.EqualTo(18.5604).Within(ParityTol), "primer3 ntthal Tm.");
            Assert.That(v.DeltaH, Is.EqualTo(-68.40000).Within(ParityTol), "primer3 ntthal ΔH (kcal/mol).");
            Assert.That(v.DeltaG37, Is.EqualTo(-6.89198).Within(ParityTol), "primer3 ntthal ΔG°37 (kcal/mol).");
        });
    }

    // N5 — terminal overhang / dangling end. GCGCGCAAAA/AAAAGCGCGC: a 6-bp GCGCGC core with AAAA
    // overhangs on each side; ntthal scores the terminal dangling stacks (tstack2 / dangle tables).
    [Test]
    public void NtthalDimer_TerminalOverhang_MatchesPrimer3()
    {
        var v = PrimerDesigner.CalculateDimerThermodynamicsNtthal("GCGCGCAAAA", "AAAAGCGCGC", Na, Ct)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(v.TmCelsius, Is.EqualTo(24.6547).Within(ParityTol), "primer3 ntthal Tm.");
            Assert.That(v.DeltaH, Is.EqualTo(-60.00000).Within(ParityTol), "primer3 ntthal ΔH (kcal/mol).");
            Assert.That(v.DeltaG37, Is.EqualTo(-8.72844).Within(ParityTol), "primer3 ntthal ΔG°37 (kcal/mol).");
        });
    }

    // N6 — contiguous-WC regression through the full DP: the ntthal DP must STILL match primer3
    // (and the legacy contiguous result) for a fully complementary duplex. GCGCGCGC self-dimer.
    [Test]
    public void NtthalDimer_ContiguousDuplex_RegressionMatchesPrimer3()
    {
        var v = PrimerDesigner.CalculateDimerThermodynamicsNtthal("GCGCGCGC", "GCGCGCGC", Na, Ct)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(v.TmCelsius, Is.EqualTo(40.0906).Within(ParityTol), "Contiguous regression: ntthal Tm.");
            Assert.That(v.DeltaH, Is.EqualTo(-70.8).Within(ParityTol), "Contiguous regression: ΔH (kcal/mol).");
            Assert.That(v.BasePairs, Is.EqualTo(8), "Full 8-bp contiguous duplex.");
        });
    }

    // N7 — CalculateDimerMeltingTemperature now routes through the full ntthal DP, so on a
    // non-contiguous optimum (the 3×3 internal loop of GCGCAAAGCGC/GCGCTTTGCGC) it reports the
    // ntthal internal-loop Tm rather than the cold best-contiguous-run Tm. The Tm method and the
    // full thermodynamics method must agree, and both must match primer3-py (41.88 °C).
    [Test]
    public void NtthalDimer_DimerTmMethodDelegatesToFullDp()
    {
        double tm = PrimerDesigner.CalculateDimerMeltingTemperature("GCGCAAAGCGC", "GCGCTTTGCGC", Na, Ct);
        var full = PrimerDesigner.CalculateDimerThermodynamicsNtthal("GCGCAAAGCGC", "GCGCTTTGCGC", Na, Ct)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(tm, Is.EqualTo(full.TmCelsius).Within(ExactTol),
                "CalculateDimerMeltingTemperature delegates to the full ntthal DP.");
            Assert.That(tm, Is.EqualTo(41.8816).Within(ParityTol),
                "Full ntthal DP reproduces the primer3-py internal-loop dimer Tm.");
        });
    }

    // N8 — no stable dimer through the full DP: poly-A vs poly-A returns null (ntthal no_structure).
    [Test]
    public void NtthalDimer_NonComplementary_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateDimerThermodynamicsNtthal("AAAAAAAA", "AAAAAAAA", Na, Ct), Is.Null,
                "poly-A cannot pair with poly-A — ntthal structure_found=False.");
            Assert.That(PrimerDesigner.CalculateDimerThermodynamicsNtthal("ACGN", "ACGT", Na, Ct), Is.Null,
                "Non-ACGT input ⇒ null.");
            Assert.That(PrimerDesigner.CalculateDimerThermodynamicsNtthal(null!, "ACGT", Na, Ct), Is.Null,
                "Null input ⇒ null.");
        });
    }

    #endregion

    #region No-dimer and invalid input

    // M9 — poly-A self-dimer: no Watson-Crick duplex forms (primer3 structure_found=False).
    [Test]
    public void FindMostStableDimer_PolyASelfDimer_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindMostStableDimer("AAAAAAAA", "AAAAAAAA", Na, Ct), Is.Null,
                "A poly-A oligo cannot pair with a copy of itself — no self-dimer (primer3 structure_found=False).");
            Assert.That(PrimerDesigner.CalculateSelfDimerMeltingTemperature("AAAAAAAA", Na, Ct), Is.NaN,
                "No self-dimer ⇒ Tm is NaN.");
        });
    }

    // M10 — fully non-complementary pair: no duplex.
    [Test]
    public void FindMostStableDimer_NonComplementaryPair_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindMostStableDimer("GGGGGGGG", "AAAAAAAA", Na, Ct), Is.Null,
                "poly-G and poly-A share no Watson-Crick base pair.");
            Assert.That(PrimerDesigner.CalculateDimerMeltingTemperature("GGGGGGGG", "AAAAAAAA", Na, Ct), Is.NaN,
                "No duplex ⇒ Tm is NaN.");
        });
    }

    // M11 — invalid input: null, too short, non-ACGT.
    [Test]
    public void FindMostStableDimer_InvalidInput_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.FindMostStableDimer(null!, "ACGT", Na, Ct), Is.Null, "Null strand 1.");
            Assert.That(PrimerDesigner.FindMostStableDimer("ACGT", null!, Na, Ct), Is.Null, "Null strand 2.");
            Assert.That(PrimerDesigner.FindMostStableDimer("", "ACGT", Na, Ct), Is.Null, "Empty strand.");
            Assert.That(PrimerDesigner.FindMostStableDimer("A", "T", Na, Ct), Is.Null, "Single base < 2.");
            Assert.That(PrimerDesigner.FindMostStableDimer("ACGN", "NCGT", Na, Ct), Is.Null, "Non-ACGT base.");
            Assert.That(PrimerDesigner.CalculateDimerMeltingTemperature("ACGN", "ACGT", Na, Ct), Is.NaN,
                "Invalid input ⇒ Tm is NaN.");
        });
    }

    #endregion

    #region Delegation, salt, concentration, selection (SHOULD / COULD)

    // S1 — CalculateSelfDimerMeltingTemperature delegates to the two-argument method.
    [Test]
    public void CalculateSelfDimerMeltingTemperature_DelegatesToTwoArgument()
    {
        double self = PrimerDesigner.CalculateSelfDimerMeltingTemperature("ACGTACGTACGT", Na, Ct);
        double two = PrimerDesigner.CalculateDimerMeltingTemperature("ACGTACGTACGT", "ACGTACGTACGT", Na, Ct);
        Assert.That(self, Is.EqualTo(two).Within(ExactTol),
            "The self-dimer wrapper must equal the two-argument call with the sequence twice.");
    }

    // S2 — INV-04: lower [Na+] strictly lowers Tm for a fixed duplex (Eq. 5 monotonicity).
    [Test]
    public void CalculateSelfDimerMeltingTemperature_LowerSodium_LowersTm()
    {
        double tmLowSalt = PrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC", 0.01, Ct);
        double tmHighSalt = PrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC", 1.0, Ct);
        Assert.That(tmLowSalt, Is.LessThan(tmHighSalt),
            "0.368·N·ln[Na⁺] makes ΔS° more negative at lower [Na⁺], lowering Tm (Eq. 5).");
    }

    // S3 — higher total strand concentration raises the bimolecular Tm.
    [Test]
    public void CalculateSelfDimerMeltingTemperature_HigherConcentration_RaisesTm()
    {
        double tm50nM = PrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC", Na, 50e-9);
        double tm500nM = PrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC", Na, 500e-9);
        Assert.That(tm500nM, Is.GreaterThan(tm50nM),
            "R·ln(C_T/x) grows with C_T, raising the bimolecular Tm.");
    }

    // C1 — most-stable selection: a longer complementary core is preferred over a short one.
    // 'GCGCGCGCAA' / 'GCGCGCGCAA' can form a short 2-bp AA-edge pairing or the strong 8-bp
    // GCGCGCGC core; the 8-bp core has the higher Tm and must be the one returned.
    [Test]
    public void FindMostStableDimer_PrefersStrongerLongerDuplex()
    {
        var d = PrimerDesigner.FindMostStableDimer("GCGCGCGCAA", "GCGCGCGCAA", Na, Ct)!.Value;
        Assert.That(d.BasePairs, Is.EqualTo(8),
            "The 8-bp GCGCGCGC core is the highest-Tm duplex; the AA tail cannot self-pair.");
    }

    #endregion
}
