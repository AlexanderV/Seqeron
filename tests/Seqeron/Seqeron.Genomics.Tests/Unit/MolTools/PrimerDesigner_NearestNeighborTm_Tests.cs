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

namespace Seqeron.Genomics.Tests.Unit.MolTools;

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

    // S3 — exactly 2 bases is the SHORTEST valid oligo (one nearest-neighbour stack), so the result must be
    // non-null. This pins the "< 2" length guard against the off-by-one boundary mutants ("<= 2" would null a
    // legal dimer; "> 2" would null every oligo of length >= 3).
    [Test]
    public void CalculateNearestNeighborThermodynamics_TwoBaseAndLongerOligos_AreNonNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamics("AT"), Is.Not.Null,
                "a 2-base oligo has exactly one nearest-neighbour stack and is valid");
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamics("ATG"), Is.Not.Null,
                "a 3-base oligo is valid (guards against the '> 2' mutant)");
        });
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

    #region Internal mismatch + dangling-end NN (opt-in extension)
    // Convention (Biopython Tm_NN with imm_table=DNA_IMM, de_table=DNA_DE): top strand
    // 5'→3'; bottom strand written 3'→5' aligned base-for-base (complement direction, NOT
    // reverse complement); '.' marks a single unpaired dangling base.
    // Sources: Allawi & SantaLucia (1997/1998); Peyret et al. (1999) [internal mismatches];
    //          Bommarito, Peyret & SantaLucia (2000) NAR 28:1929 [dangling ends].

    // MM1 — Internal single T·G mismatch (G·T pair, Allawi & SantaLucia 1997).
    // Top 5'-CGTGAC-3' / bottom 3'-GCGCTG-5' (perfect complement of CGTGAC is GCACTG;
    // column 2 changed A→G → a T·G mismatch under the top T). Hand-derived from the
    // NN + DNA_IMM tables (bottom written 3'→5'):
    //   init (+0.2,-5.7); ends C/C → no terminal-AT; not self-comp → no symmetry.
    //   CG/GC = (-10.6,-27.2)  [WC]
    //   GT/CG = (-4.4,-12.3)   [IMM, G·T]
    //   TG/GC = CG/GT reversed = (-4.1,-11.7)  [IMM, T·G]
    //   GA/CT = (-8.2,-22.2)   [WC]
    //   AC/TG = GT/CA reversed = (-8.4,-22.4)  [WC]
    //   ΔH° = 0.2 - 10.6 - 4.4 - 4.1 - 8.2 - 8.4 = -35.5 kcal/mol
    //   ΔS° = -5.7 - 27.2 - 12.3 - 11.7 - 22.2 - 22.4 = -101.5 cal/(K·mol)
    [Test]
    public void CalculateNearestNeighborThermodynamicsMismatch_InternalMismatch_MatchesImmSum()
    {
        var r = PrimerDesigner.CalculateNearestNeighborThermodynamicsMismatch("CGTGAC", "GCGCTG");

        Assert.That(r, Is.Not.Null, "Single internal mismatch is computable (each stack has an NN parameter).");
        Assert.Multiple(() =>
        {
            Assert.That(r!.Value.DeltaH, Is.EqualTo(-35.5).Within(Tol),
                "ΔH° = init + Σ(WC stacks) + Σ(internal-mismatch stacks) per Allawi & SantaLucia.");
            Assert.That(r.Value.DeltaS, Is.EqualTo(-101.5).Within(Tol),
                "ΔS° = init + Σ(WC) + Σ(internal-mismatch) cal/(K·mol).");
            Assert.That(r.Value.IsSelfComplementary, Is.False,
                "A mismatched duplex is not self-complementary.");
        });
    }

    // MM1-Tm — Tm of the MM1 mismatched duplex, no salt correction, x=4 (non-self-comp).
    // Tm = -35.5·1000/(-101.5 + R·ln(0.5e-6/4)) − 273.15 = -6.4060879279 °C.
    [Test]
    public void CalculateMeltingTemperatureNNMismatch_InternalMismatch_MatchesEquation()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNNMismatch(
            "CGTGAC", "GCGCTG", strandConcentrationMolar: 0.5e-6,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.That(tm, Is.EqualTo(-6.4060879279).Within(1e-8),
            "NN Tm of a single-internal-mismatch duplex via the bimolecular Tm equation (x=4).");
    }

    // MM2 — A mismatch DESTABILISES: the same duplex with the mismatch repaired has a
    // higher ΔH°-magnitude and a higher Tm. Repaired bottom = GCACTG (perfect complement).
    [Test]
    public void CalculateMeltingTemperatureNNMismatch_Mismatch_LowersTmVsPerfect()
    {
        double mismatched = PrimerDesigner.CalculateMeltingTemperatureNNMismatch(
            "CGTGAC", "GCGCTG", saltMode: PrimerDesigner.SaltCorrectionMode.None);
        double perfect = PrimerDesigner.CalculateMeltingTemperatureNNMismatch(
            "CGTGAC", "GCACTG", saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.That(mismatched, Is.LessThan(perfect),
            "An internal mismatch destabilises the duplex → lower Tm than the perfectly paired duplex.");
    }

    // DE1 — 5'-dangling A on a self-complementary GCGCGC core (Bommarito 2000).
    // Top 5'-AGCGCGC-3' / bottom 3'-.CGCGCG-5' ('.' under the dangling A; the paired
    // region is GCGCGC/CGCGCG). Hand-derived:
    //   init (+0.2,-5.7); ends (top[0],top[-1]) = A,C → one terminal-AT (+2.2,+6.9);
    //   dangling present → not self-comp → no symmetry.
    //   left dangling end AG/.C = (-3.7,-10.0)  [DNA_DE, 5'-dangling A on a C·G pair]
    //   then strip → GCGCGC/CGCGCG, WC stacks:
    //     GC/CG ·3 = 3·(-9.8,-24.4); CG/GC ·2 = 2·(-10.6,-27.2)
    //   ΔH° = 0.2 + 2.2 - 3.7 + 3·(-9.8) + 2·(-10.6) = -51.9 kcal/mol
    //   ΔS° = -5.7 + 6.9 - 10.0 + 3·(-24.4) + 2·(-27.2) = -136.4 cal/(K·mol)
    [Test]
    public void CalculateNearestNeighborThermodynamicsMismatch_FivePrimeDanglingEnd_MatchesDeSum()
    {
        var r = PrimerDesigner.CalculateNearestNeighborThermodynamicsMismatch("AGCGCGC", ".CGCGCG");

        Assert.That(r, Is.Not.Null, "A single 5'-dangling end has a DNA_DE parameter.");
        Assert.Multiple(() =>
        {
            Assert.That(r!.Value.DeltaH, Is.EqualTo(-51.9).Within(Tol),
                "ΔH° = init + terminal-AT + Bommarito dangling-end term + Σ(WC stacks).");
            Assert.That(r.Value.DeltaS, Is.EqualTo(-136.4).Within(Tol),
                "ΔS° = init + terminal-AT + dangling-end + Σ(WC stacks).");
            Assert.That(r.Value.IsSelfComplementary, Is.False,
                "A duplex with a dangling end is not treated as self-complementary.");
        });
    }

    // DE1-Tm — Tm of the 5'-dangling duplex, no salt, x=4.
    // Tm = -51.9·1000/(-136.4 + R·ln(0.5e-6/4)) − 273.15 = 35.8034921829 °C.
    [Test]
    public void CalculateMeltingTemperatureNNMismatch_DanglingEnd_MatchesEquation()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNNMismatch(
            "AGCGCGC", ".CGCGCG", strandConcentrationMolar: 0.5e-6,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.That(tm, Is.EqualTo(35.8034921829).Within(1e-8),
            "NN Tm of a 5'-dangling-end duplex via the bimolecular Tm equation (x=4).");
    }

    // EQ1 — A perfectly complementary duplex through the mismatch/dangling path equals the
    // existing perfect-match CalculateMeltingTemperatureNN exactly (the extension is a
    // strict superset). GCGCGC top / CGCGCG bottom (3'→5') is fully Watson-Crick paired.
    [Test]
    public void CalculateMeltingTemperatureNNMismatch_PerfectDuplex_EqualsPerfectMatchPath()
    {
        double viaExtension = PrimerDesigner.CalculateMeltingTemperatureNNMismatch(
            "GCGCGC", "CGCGCG", strandConcentrationMolar: 0.5e-6,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);
        double viaPerfect = PrimerDesigner.CalculateMeltingTemperatureNN(
            "GCGCGC", strandConcentrationMolar: 0.5e-6,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.That(viaExtension, Is.EqualTo(viaPerfect).Within(1e-9),
            "A fully paired duplex through the extension must equal the perfect-match NN Tm.");
    }

    // EQ2 — The perfect-duplex thermodynamics also match the perfect-match path's
    // (ΔH°, ΔS°, IsSelfComplementary) tuple, including the self-complementary flag.
    [Test]
    public void CalculateNearestNeighborThermodynamicsMismatch_PerfectDuplex_EqualsPerfectMatchThermo()
    {
        var ext = PrimerDesigner.CalculateNearestNeighborThermodynamicsMismatch("GCGCGC", "CGCGCG");
        var perfect = PrimerDesigner.CalculateNearestNeighborThermodynamics("GCGCGC");

        Assert.That(ext, Is.Not.Null);
        Assert.That(perfect, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(ext!.Value.DeltaH, Is.EqualTo(perfect!.Value.DeltaH).Within(Tol),
                "Perfect-duplex ΔH° must match the perfect-match path.");
            Assert.That(ext.Value.DeltaS, Is.EqualTo(perfect.Value.DeltaS).Within(Tol),
                "Perfect-duplex ΔS° must match the perfect-match path.");
            Assert.That(ext.Value.IsSelfComplementary, Is.EqualTo(perfect.Value.IsSelfComplementary),
                "Self-complementary detection must agree (GCGCGC is self-complementary).");
        });
    }

    // C3 — Invalid inputs → NaN / null: null strands, unequal length, a tandem (two adjacent)
    // mismatch with no NN parameter.
    [Test]
    public void CalculateMeltingTemperatureNNMismatch_InvalidOrUncomputable_ReturnsNaN()
    {
        Assert.Multiple(() =>
        {
            Assert.That(double.IsNaN(PrimerDesigner.CalculateMeltingTemperatureNNMismatch(null!, "AT")), Is.True,
                "Null top strand → NaN.");
            Assert.That(double.IsNaN(PrimerDesigner.CalculateMeltingTemperatureNNMismatch("ATGC", "TAC")), Is.True,
                "Unequal-length strands → NaN.");
            // Tandem mismatch: top GG / bottom GG (3'→5') = both columns mismatched, no NN term.
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsMismatch("AGGT", "TGGA"), Is.Null,
                "A stack with two adjacent mismatches has no NN parameter → not computable.");
        });
    }

    #endregion
}
