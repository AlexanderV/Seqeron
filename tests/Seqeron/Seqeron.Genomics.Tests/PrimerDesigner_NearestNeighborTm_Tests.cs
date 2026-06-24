// PRIMER-TM-001 — Nearest-neighbour salt-corrected melting temperature (opt-in)
// Evidence: docs/Evidence/PRIMER-TM-001-Evidence.md
// TestSpec: tests/TestSpecs/PRIMER-TM-001-NN.md
// Sources:
//   SantaLucia J (1998). PNAS 95(4):1460-65 — unified NN ΔH°/ΔS° (Table 1), Tm Eq. 3.
//   SantaLucia J, Hicks D (2004). Annu Rev Biophys 33:415-440 — Table 1 + Eq. 3 + Eq. 5 (cross-check).
//   Owczarzy R et al. (2004). Biochemistry 43:3537-54 — monovalent Na⁺ correction.
//   Owczarzy R et al. (2008). Biochemistry 47:5336-53 — divalent Mg²⁺ correction.
//   Biopython Bio.SeqUtils.MeltingTemp (DNA_NN4, salt_correction methods 6/7) — reference impl.

using NUnit.Framework;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the opt-in nearest-neighbour (SantaLucia 1998) salt-corrected melting
/// temperature added under PRIMER-TM-001. Expected values are hand-derived from the
/// unified NN ΔH°/ΔS° table and the bimolecular Tm equation, independent of the code.
/// </summary>
[TestFixture]
public class PrimerDesigner_NearestNeighborTm_Tests
{
    private const double Tol = 1e-9;

    #region CalculateNearestNeighborThermodynamics — ΔH°/ΔS° sums

    // M1 — Non-self-complementary oligo ATGCATGC.
    // NN stacks (SantaLucia&Hicks 2004 Table 1):
    //   AT(-7.2,-20.4) TG(-8.5,-22.7) GC(-9.8,-24.4) CA(-8.5,-22.7)
    //   AT(-7.2,-20.4) TG(-8.5,-22.7) GC(-9.8,-24.4)
    // init (+0.2,-5.7); both ends A/T → 2× terminal-AT (+2.2,+6.9); not self-comp → no symmetry.
    // ΔH° = 0.2 + (-7.2-8.5-9.8-8.5-7.2-8.5-9.8) + 2.2 + 2.2 = -57.1 kcal/mol
    // ΔS° = -5.7 + (-20.4-22.7-24.4-22.7-20.4-22.7-24.4) + 6.9 + 6.9 = -156.5 cal/(K·mol)
    [Test]
    public void CalculateNearestNeighborThermodynamics_NonSelfComp_MatchesTable1Sum()
    {
        var r = PrimerDesigner.CalculateNearestNeighborThermodynamics("ATGCATGC");

        Assert.That(r, Is.Not.Null, "ATGCATGC has only ACGT bases and length ≥ 2.");
        Assert.Multiple(() =>
        {
            Assert.That(r!.Value.DeltaH, Is.EqualTo(-57.1).Within(Tol),
                "ΔH° = init(+0.2) + Σ NN stacks + 2× terminal-AT(+2.2) per SantaLucia&Hicks 2004 Table 1.");
            Assert.That(r.Value.DeltaS, Is.EqualTo(-156.5).Within(Tol),
                "ΔS° = init(-5.7) + Σ NN stacks + 2× terminal-AT(+6.9); no symmetry (non-self-comp).");
            Assert.That(r.Value.IsSelfComplementary, Is.False,
                "revcomp(ATGCATGC)=GCATGCAT ≠ ATGCATGC → not self-complementary.");
        });
    }

    // M2 — Self-complementary oligo GCGCGC (revcomp = GCGCGC).
    // NN stacks: GC(-9.8,-24.4) CG(-10.6,-27.2) GC(-9.8,-24.4) CG(-10.6,-27.2) GC(-9.8,-24.4)
    // init (+0.2,-5.7); both ends G/C → no terminal-AT; self-comp → symmetry ΔS° -1.4.
    // ΔH° = 0.2 + (-9.8-10.6-9.8-10.6-9.8) = -50.4 kcal/mol
    // ΔS° = -5.7 + (-24.4-27.2-24.4-27.2-24.4) + (-1.4) = -134.7 cal/(K·mol)
    [Test]
    public void CalculateNearestNeighborThermodynamics_SelfComp_AppliesSymmetryNoTerminalAt()
    {
        var r = PrimerDesigner.CalculateNearestNeighborThermodynamics("GCGCGC");

        Assert.That(r, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(r!.Value.DeltaH, Is.EqualTo(-50.4).Within(Tol),
                "ΔH° = init(+0.2) + Σ NN stacks; symmetry ΔH°=0, no terminal-AT (G/C ends).");
            Assert.That(r.Value.DeltaS, Is.EqualTo(-134.7).Within(Tol),
                "ΔS° = init(-5.7) + Σ NN stacks + symmetry(-1.4).");
            Assert.That(r.Value.IsSelfComplementary, Is.True,
                "GCGCGC equals its own reverse complement.");
        });
    }

    // S1 — non-ACGT base makes NN lookup fail → null.
    [Test]
    public void CalculateNearestNeighborThermodynamics_NonAcgt_ReturnsNull()
    {
        Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamics("ATGN"), Is.Null,
            "An N base has no NN parameter; thermodynamics is not computable.");
    }

    // S2 — too short (< 2 bases) → null.
    [Test]
    public void CalculateNearestNeighborThermodynamics_SingleBase_ReturnsNull()
    {
        Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamics("A"), Is.Null,
            "A single base has no nearest-neighbour stack.");
    }

    #endregion

    #region CalculateMeltingTemperatureNN — Tm equation + salt corrections

    // M3 — Published worked example (SantaLucia&Hicks 2004, p.419):
    // ΔH°=-43.5 kcal/mol, ΔS°=-122.5 e.u., 0.2 mM each strand (C_T=0.0004, x=4) → Tm = 35.8 °C.
    // Reproduced via the same Tm equation the method uses; locks the equation/constants
    // independently of the NN table.
    [Test]
    public void TmEquation_PublishedWorkedExample_Gives35Point8()
    {
        // Tm = ΔH°·1000 / (ΔS° + R·ln(C_T/x)) − 273.15
        const double R = 1.9872;
        double tmK = (-43.5 * 1000.0) / (-122.5 + R * System.Math.Log(0.0004 / 4.0));
        double tmC = tmK - 273.15;

        Assert.That(tmC, Is.EqualTo(35.8).Within(0.05),
            "SantaLucia&Hicks 2004 p.419 worked example: ΔH°=-43.5, ΔS°=-122.5, 0.2 mM → 35.8 °C.");
    }

    // M4 — Full NN Tm, no salt correction (1 M NaCl reference state), self-comp GCGCGC.
    // C_T=0.5 µM, x=1. Tm = -50.4·1000/(-134.7 + R·ln(0.5e-6/1)) − 273.15 = 35.0473059911 °C.
    [Test]
    public void CalculateMeltingTemperatureNN_NoSalt_SelfComp_MatchesEquation()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNN(
            "GCGCGC", strandConcentrationMolar: 0.5e-6,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.That(tm, Is.EqualTo(35.0473059911).Within(1e-8),
            "Tm = ΔH°·1000/(ΔS° + R·ln(C_T/1)) − 273.15 with x=1 (self-complementary).");
    }

    // M5 — Full NN Tm, no salt correction, non-self-comp ATGCATGC, x=4.
    // Tm = -57.1·1000/(-156.5 + R·ln(0.5e-6/4)) − 273.15 = 30.4338060665 °C.
    [Test]
    public void CalculateMeltingTemperatureNN_NoSalt_NonSelfComp_UsesX4()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNN(
            "ATGCATGC", strandConcentrationMolar: 0.5e-6,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.That(tm, Is.EqualTo(30.4338060665).Within(1e-8),
            "Non-self-complementary duplex uses x=4 in the strand-concentration term (Eq. 3).");
    }

    // M6 — Owczarzy 2004 monovalent correction at 50 mM Na⁺, self-comp GCGCGC.
    // 1/Tm[Na] = 1/Tm[1M] + (4.29e-5·fGC − 3.95e-5)·ln[Na] + 9.40e-6·(ln[Na])²; fGC=1.0.
    // Result = 28.1593085080 °C (cooler than the 1 M value, 35.05 °C — physically correct).
    [Test]
    public void CalculateMeltingTemperatureNN_Owczarzy2004_50mM_LowersTm()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNN(
            "GCGCGC", strandConcentrationMolar: 0.5e-6, sodiumMolar: 0.05,
            saltMode: PrimerDesigner.SaltCorrectionMode.Owczarzy2004Monovalent);

        Assert.That(tm, Is.EqualTo(28.1593085080).Within(1e-8),
            "Owczarzy 2004 quadratic 1/Tm Na⁺ correction at 50 mM (Biochemistry 43:3537).");
    }

    // M7 — Owczarzy 2004 at 50 mM Na⁺, non-self-comp ATGCATGC (fGC=0.5).
    // Result = 18.1899960529 °C.
    [Test]
    public void CalculateMeltingTemperatureNN_Owczarzy2004_NonSelfComp_MatchesDerivation()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNN(
            "ATGCATGC", strandConcentrationMolar: 0.5e-6, sodiumMolar: 0.05,
            saltMode: PrimerDesigner.SaltCorrectionMode.Owczarzy2004Monovalent);

        Assert.That(tm, Is.EqualTo(18.1899960529).Within(1e-8),
            "Owczarzy 2004 correction with GC fraction 0.5 and x=4.");
    }

    // M8 — Owczarzy 2004 is the DEFAULT salt mode: omitting saltMode equals the explicit mode.
    [Test]
    public void CalculateMeltingTemperatureNN_DefaultSaltMode_IsOwczarzy2004()
    {
        double withDefault = PrimerDesigner.CalculateMeltingTemperatureNN(
            "ATGCATGC", strandConcentrationMolar: 0.5e-6, sodiumMolar: 0.05);
        double withExplicit = PrimerDesigner.CalculateMeltingTemperatureNN(
            "ATGCATGC", strandConcentrationMolar: 0.5e-6, sodiumMolar: 0.05,
            saltMode: PrimerDesigner.SaltCorrectionMode.Owczarzy2004Monovalent);

        Assert.That(withDefault, Is.EqualTo(withExplicit).Within(Tol),
            "Default saltMode must be Owczarzy2004Monovalent.");
    }

    // S3 — SantaLucia entropy correction (Eq. 5) at 50 mM, self-comp GCGCGC.
    // ΔS°[Na] = ΔS° + 0.368·(N/2)·ln[Na], N=2·(L−1)=10 → N/2=5. Result = 24.9976652723 °C.
    [Test]
    public void CalculateMeltingTemperatureNN_SantaLuciaEntropy_50mM_MatchesEq5()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNN(
            "GCGCGC", strandConcentrationMolar: 0.5e-6, sodiumMolar: 0.05,
            saltMode: PrimerDesigner.SaltCorrectionMode.SantaLuciaEntropy);

        Assert.That(tm, Is.EqualTo(24.9976652723).Within(1e-8),
            "SantaLucia&Hicks 2004 Eq. 5 entropy salt correction with N=10 phosphates.");
    }

    // S4 — EcoRI self-complementary 12-mer CGCGAATTCGCG at 1 M reference → 61.1452300219 °C.
    [Test]
    public void CalculateMeltingTemperatureNN_EcoRiSelfComp_NoSalt_MatchesEquation()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNN(
            "CGCGAATTCGCG", strandConcentrationMolar: 0.5e-6,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.That(tm, Is.EqualTo(61.1452300219).Within(1e-8),
            "Self-complementary EcoRI 12-mer Tm at the 1 M NaCl reference state.");
    }

    // S5 — Lowering [Na⁺] monotonically lowers the Owczarzy-2004-corrected Tm.
    [Test]
    public void CalculateMeltingTemperatureNN_LowerSodium_LowersTm()
    {
        double tmHigh = PrimerDesigner.CalculateMeltingTemperatureNN(
            "CGCGAATTCGCG", sodiumMolar: 1.0);
        double tmLow = PrimerDesigner.CalculateMeltingTemperatureNN(
            "CGCGAATTCGCG", sodiumMolar: 0.01);

        Assert.That(tmLow, Is.LessThan(tmHigh),
            "Lower monovalent salt destabilises the duplex → lower Tm (Owczarzy 2004 monotonicity).");
    }

    // C1 — Divalent mode with no Mg²⁺ falls back to the 2004 monovalent value.
    [Test]
    public void CalculateMeltingTemperatureNN_DivalentNoMg_EqualsMonovalent()
    {
        double divalent = PrimerDesigner.CalculateMeltingTemperatureNN(
            "ATGCATGC", strandConcentrationMolar: 0.5e-6, sodiumMolar: 0.05, magnesiumMolar: 0.0,
            saltMode: PrimerDesigner.SaltCorrectionMode.Owczarzy2008Divalent);
        double monovalent = PrimerDesigner.CalculateMeltingTemperatureNN(
            "ATGCATGC", strandConcentrationMolar: 0.5e-6, sodiumMolar: 0.05,
            saltMode: PrimerDesigner.SaltCorrectionMode.Owczarzy2004Monovalent);

        Assert.That(divalent, Is.EqualTo(monovalent).Within(1e-8),
            "With [Mg²⁺]=0 the Owczarzy 2008 divalent model reduces to the 2004 monovalent form.");
    }

    // C2 — Adding Mg²⁺ raises Tm relative to the same buffer with no Mg²⁺ (divalent stabilises).
    [Test]
    public void CalculateMeltingTemperatureNN_AddMagnesium_RaisesTm()
    {
        double noMg = PrimerDesigner.CalculateMeltingTemperatureNN(
            "CGCGAATTCGCG", sodiumMolar: 0.05, magnesiumMolar: 0.0,
            saltMode: PrimerDesigner.SaltCorrectionMode.Owczarzy2008Divalent);
        double withMg = PrimerDesigner.CalculateMeltingTemperatureNN(
            "CGCGAATTCGCG", sodiumMolar: 0.05, magnesiumMolar: 0.003,
            saltMode: PrimerDesigner.SaltCorrectionMode.Owczarzy2008Divalent);

        Assert.That(withMg, Is.GreaterThan(noMg),
            "Divalent Mg²⁺ stabilises the duplex → higher Tm than the Mg²⁺-free buffer.");
    }

    #endregion

    #region Edge cases and default-unchanged guarantee

    // E1 — empty / null / non-ACGT → NaN (no crash).
    [Test]
    public void CalculateMeltingTemperatureNN_InvalidInput_ReturnsNaN()
    {
        Assert.Multiple(() =>
        {
            Assert.That(double.IsNaN(PrimerDesigner.CalculateMeltingTemperatureNN("")), Is.True);
            Assert.That(double.IsNaN(PrimerDesigner.CalculateMeltingTemperatureNN(null!)), Is.True);
            Assert.That(double.IsNaN(PrimerDesigner.CalculateMeltingTemperatureNN("ATGN")), Is.True);
            Assert.That(double.IsNaN(PrimerDesigner.CalculateMeltingTemperatureNN("A")), Is.True);
        });
    }

    // E2 — The legacy default Tm method is UNCHANGED by this opt-in addition.
    // Wallace for ATATATAT = 2·8 = 16 °C (unchanged baseline).
    [Test]
    public void CalculateMeltingTemperature_DefaultMethod_Unchanged()
    {
        Assert.That(PrimerDesigner.CalculateMeltingTemperature("ATATATAT"), Is.EqualTo(16.0),
            "The default Wallace/Marmur-Doty Tm must be unchanged by the opt-in NN method.");
    }

    #endregion
}
