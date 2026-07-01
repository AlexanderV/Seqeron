// SEQ-THERMO-001 — DNA Duplex Thermodynamics (Nearest-Neighbor ΔH°/ΔS°/ΔG°/Tm)
// Evidence: docs/Evidence/SEQ-THERMO-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-THERMO-001.md
// Source: Allawi HT, SantaLucia J Jr (1997). Biochemistry 36(34):10581-10594.
//         SantaLucia J Jr (1998). PNAS 95(4):1460-1465.
//         Biopython Bio.SeqUtils.MeltingTemp (DNA_NN3 table, Tm_NN).
//         MELTING 5 User Guide (Dumousseau et al. 2012) §4.2/§4.3.

using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_CalculateThermodynamics_Tests
{
    // Expected values are derived from the DNA_NN3 NN parameters and the
    // SantaLucia (1998) salt + Tm equations; rounding matches the implementation
    // (ΔH/ΔS/ΔG to 2 decimals, Tm to 1 decimal).
    private const double Tolerance = 1e-10;

    #region CalculateThermodynamics

    // M1 — Biopython worked example: Tm_NN('CGTTCCAAAGATGTGGGCATGAGCTTAC') = 60.32 °C
    //      at dnac1 = dnac2 = 25 nM (k = C_T/4 with C_T = 50 nM) and Na = 50 mM.
    // Evidence: Biopython Tm_NN docstring. Rounds to 60.3 at one decimal.
    [Test]
    public void CalculateThermodynamics_BiopythonWorkedExample_ReturnsTm60Point3()
    {
        var result = SequenceStatistics.CalculateThermodynamics(
            "CGTTCCAAAGATGTGGGCATGAGCTTAC",
            naConcentration: 0.05,
            primerConcentration: 0.00000005); // C_T = 50 nM => C_T/4 = 12.5 nM (Biopython k)

        Assert.That(result.MeltingTemperature, Is.EqualTo(60.3).Within(Tolerance),
            "Reproduces the Biopython Tm_NN reference value (60.32 -> 60.3) for the documented sequence (INV-03)");
    }

    // M2 — GCGC full result under repository defaults (Na = 0.05 M, C_T = 250 nM).
    // Evidence: init G+C twice = ΔH 0.2 / ΔS -5.6; NN GC,CG,GC = ΔH -30.2 / ΔS -76.0;
    //           salt ΔS = 0.368*3*ln(0.05) = -3.307; ΔH=-30.0, ΔS=-84.91, ΔG=-3.67, Tm=-18.6.
    [Test]
    public void CalculateThermodynamics_GcgcDefaults_ReturnsExactTuple()
    {
        var result = SequenceStatistics.CalculateThermodynamics("GCGC");

        Assert.Multiple(() =>
        {
            Assert.That(result.DeltaH, Is.EqualTo(-30.0).Within(Tolerance),
                "ΔH° = init(0.2) + NN(GC,CG,GC = -30.2) = -30.0 kcal/mol");
            Assert.That(result.DeltaS, Is.EqualTo(-84.91).Within(Tolerance),
                "ΔS° = init(-5.6) + NN(-76.0) + salt(-3.307) = -84.91 cal/(mol·K)");
            Assert.That(result.DeltaG, Is.EqualTo(-3.67).Within(Tolerance),
                "ΔG°37 = ΔH - 310.15·ΔS/1000 = -3.67 kcal/mol (INV-02)");
            Assert.That(result.MeltingTemperature, Is.EqualTo(-18.6).Within(Tolerance),
                "Tm = (1000·ΔH)/(ΔS + R·ln(C_T/4)) - 273.15 = -18.6 °C (INV-03)");
        });
    }

    // M3 — Two-end initiation with mixed termini: ATCG (A...G) gets init_A/T + init_G/C.
    // Evidence: Biopython two-end init (ends = seq[0]+seq[-1]). ΔH=-23.6, ΔS=-71.81.
    [Test]
    public void CalculateThermodynamics_MixedTermini_AppliesInitiationAtBothEnds()
    {
        var result = SequenceStatistics.CalculateThermodynamics("ATCG");

        Assert.Multiple(() =>
        {
            Assert.That(result.DeltaH, Is.EqualTo(-23.6).Within(Tolerance),
                "ΔH° includes init_A/T (2.3) at the A-end and init_G/C (0.1) at the G-end (INV-01)");
            Assert.That(result.DeltaS, Is.EqualTo(-71.81).Within(Tolerance),
                "ΔS° includes init_A/T (4.1) and init_G/C (-2.8) plus salt correction (INV-01)");
        });
    }

    // M4 — Both termini A/T: AATT applies init_A/T twice.
    // Evidence: DNA_NN3 init_A/T = (2.3, 4.1); derivation ΔH=-18.4, ΔS=-59.91, ΔG=0.18, Tm=-75.0.
    [Test]
    public void CalculateThermodynamics_BothAtTermini_ReturnsExactTuple()
    {
        var result = SequenceStatistics.CalculateThermodynamics("AATT");

        Assert.Multiple(() =>
        {
            Assert.That(result.DeltaH, Is.EqualTo(-18.4).Within(Tolerance),
                "ΔH° = init_A/T twice (4.6) + NN(AA,AT,TT = -23.0) = -18.4");
            Assert.That(result.DeltaS, Is.EqualTo(-59.91).Within(Tolerance),
                "ΔS° = init_A/T twice (8.2) + NN(-63.9) + salt(-4.21) = -59.91");
            Assert.That(result.DeltaG, Is.EqualTo(0.18).Within(Tolerance),
                "ΔG°37 = -18.4 - 310.15·(-59.91)/1000 = 0.18 (INV-02)");
            Assert.That(result.MeltingTemperature, Is.EqualTo(-75.0).Within(Tolerance),
                "Tm for the weakly-stable AT-rich duplex (INV-03)");
        });
    }

    // M5 — Empty input returns all-zero (NN model undefined for length < 2).
    [Test]
    public void CalculateThermodynamics_EmptyInput_ReturnsAllZero()
    {
        var result = SequenceStatistics.CalculateThermodynamics("");

        Assert.Multiple(() =>
        {
            Assert.That(result.DeltaH, Is.EqualTo(0.0), "Empty input has no NN step (INV-06)");
            Assert.That(result.DeltaS, Is.EqualTo(0.0), "Empty input has no NN step (INV-06)");
            Assert.That(result.DeltaG, Is.EqualTo(0.0), "Empty input has no NN step (INV-06)");
            Assert.That(result.MeltingTemperature, Is.EqualTo(0.0), "Empty input has no NN step (INV-06)");
        });
    }

    // M6 — Length-1 input returns all-zero (no dinucleotide exists).
    [Test]
    public void CalculateThermodynamics_SingleBase_ReturnsAllZero()
    {
        var result = SequenceStatistics.CalculateThermodynamics("A");

        Assert.Multiple(() =>
        {
            Assert.That(result.DeltaH, Is.EqualTo(0.0), "Length-1 has no dinucleotide (INV-06)");
            Assert.That(result.DeltaS, Is.EqualTo(0.0), "Length-1 has no dinucleotide (INV-06)");
            Assert.That(result.DeltaG, Is.EqualTo(0.0), "Length-1 has no dinucleotide (INV-06)");
            Assert.That(result.MeltingTemperature, Is.EqualTo(0.0), "Length-1 has no dinucleotide (INV-06)");
        });
    }

    // M6b — Null input returns all-zero (guarded by string.IsNullOrEmpty; no throw).
    // The NN model is undefined for length < 2; the API contract returns (0,0,0,0).
    [Test]
    public void CalculateThermodynamics_NullInput_ReturnsAllZero()
    {
        var result = SequenceStatistics.CalculateThermodynamics(null!);

        Assert.Multiple(() =>
        {
            Assert.That(result.DeltaH, Is.EqualTo(0.0), "Null input returns zero (INV-06)");
            Assert.That(result.DeltaS, Is.EqualTo(0.0), "Null input returns zero (INV-06)");
            Assert.That(result.DeltaG, Is.EqualTo(0.0), "Null input returns zero (INV-06)");
            Assert.That(result.MeltingTemperature, Is.EqualTo(0.0), "Null input returns zero (INV-06)");
        });
    }

    // M7 — Gibbs relation: ΔG°37 = ΔH° - 310.15·ΔS°/1000 must hold for the reported values.
    // Evidence: SantaLucia (1998) ΔG at 310.15 K.
    [Test]
    public void CalculateThermodynamics_DeltaG_SatisfiesGibbsRelation()
    {
        var result = SequenceStatistics.CalculateThermodynamics("GCGC");

        double expectedDeltaG = Math.Round(result.DeltaH - (310.15 * result.DeltaS / 1000.0), 2);

        Assert.That(result.DeltaG, Is.EqualTo(expectedDeltaG).Within(Tolerance),
            "ΔG°37 equals ΔH - 310.15·ΔS/1000 recomputed from the reported ΔH/ΔS (INV-02)");
    }

    // S1 — Case-insensitivity: lowercase equals uppercase.
    [Test]
    public void CalculateThermodynamics_LowercaseInput_EqualsUppercase()
    {
        var lower = SequenceStatistics.CalculateThermodynamics("gcgc");
        var upper = SequenceStatistics.CalculateThermodynamics("GCGC");

        Assert.Multiple(() =>
        {
            Assert.That(lower.DeltaH, Is.EqualTo(upper.DeltaH).Within(Tolerance), "Case-insensitive ΔH (INV-05)");
            Assert.That(lower.MeltingTemperature, Is.EqualTo(upper.MeltingTemperature).Within(Tolerance),
                "Case-insensitive Tm (INV-05)");
        });
    }

    // S2 — Higher [Na+] raises Tm (salt entropy term is monotonic).
    // Evidence: salt correction 0.368·(N-1)·ln[Na+]; GCGC Tm -18.6 (0.05 M) vs -11.3 (1.0 M).
    [Test]
    public void CalculateThermodynamics_HigherSalt_RaisesMeltingTemperature()
    {
        double tmLow = SequenceStatistics.CalculateThermodynamics("GCGC", naConcentration: 0.05).MeltingTemperature;
        double tmHigh = SequenceStatistics.CalculateThermodynamics("GCGC", naConcentration: 1.0).MeltingTemperature;

        Assert.Multiple(() =>
        {
            Assert.That(tmLow, Is.EqualTo(-18.6).Within(Tolerance), "Tm at 50 mM Na+");
            Assert.That(tmHigh, Is.EqualTo(-11.3).Within(Tolerance), "Tm at 1 M Na+");
            Assert.That(tmHigh, Is.GreaterThan(tmLow), "Higher monovalent salt stabilizes the duplex (raises Tm)");
        });
    }

    // S3 — NN table Watson-Crick symmetry: an AA run and its TT mirror give equal ΔH/ΔS.
    // Evidence: DNA_NN3 AA/TT identical; init_A/T identical at both ends for both.
    [Test]
    public void CalculateThermodynamics_WatsonCrickSymmetry_AaRunEqualsTtRun()
    {
        var aa = SequenceStatistics.CalculateThermodynamics("AAAA");
        var tt = SequenceStatistics.CalculateThermodynamics("TTTT");

        Assert.Multiple(() =>
        {
            Assert.That(aa.DeltaH, Is.EqualTo(tt.DeltaH).Within(Tolerance),
                "AA and TT NN steps share ΔH; A/T termini identical (INV-04)");
            Assert.That(aa.DeltaS, Is.EqualTo(tt.DeltaS).Within(Tolerance),
                "AA and TT NN steps share ΔS (INV-04)");
        });
    }

    #endregion

    #region CalculateMeltingTemperature (Delegate — smoke)

    // C1 — Wallace rule for short oligos (< 14 bp): Tm = 2(A+T) + 4(G+C).
    // Evidence: Wallace rule (ThermoConstants). ATGC: 2*2 + 4*2 = 12.
    [Test]
    public void CalculateMeltingTemperature_ShortOligoWallace_ReturnsTwoAtPlusFourGc()
    {
        double tm = SequenceStatistics.CalculateMeltingTemperature("ATGC", useWallaceRule: true);

        Assert.That(tm, Is.EqualTo(12.0).Within(Tolerance),
            "Wallace rule: 2*(A+T)=4 plus 4*(G+C)=8 equals 12 °C");
    }

    // C2 — Marmur-Doty GC formula (useWallaceRule=false): 64.9 + 41*(GC-16.4)/N.
    // Evidence: Marmur-Doty (ThermoConstants). 20-mer with 10 GC: 64.9 + 41*(10-16.4)/20 = 51.78.
    [Test]
    public void CalculateMeltingTemperature_MarmurDoty_ReturnsGcFormulaValue()
    {
        double tm = SequenceStatistics.CalculateMeltingTemperature(
            "GCGCGCGCGCATATATATAT", useWallaceRule: false); // 20 nt, 10 G/C

        Assert.That(tm, Is.EqualTo(51.78).Within(1e-9),
            "Marmur-Doty: 64.9 + 41*(10-16.4)/20 = 51.78 °C");
    }

    // C2b — Auto-switch branch: useWallaceRule=true but length >= 14 falls through to
    //       the Marmur-Doty GC formula (the Wallace rule of thumb is only applied below
    //       WallaceMaxLength = 14). 'ACGTTGCAATGCCGTA' is 16 nt, GC = 8.
    // Evidence: Marmur-Doty/GC form Tm = 64.9 + 41*(GC-16.4)/N (Marmur & Doty 1962;
    //           UGENE/Primer3 GC method). 64.9 + 41*(8-16.4)/16 = 43.375.
    //           NOT the Wallace value 48.0, because length 16 >= 14 disables Wallace.
    [Test]
    public void CalculateMeltingTemperature_WallaceRequestedButTooLong_UsesMarmurDoty()
    {
        double tm = SequenceStatistics.CalculateMeltingTemperature(
            "ACGTTGCAATGCCGTA", useWallaceRule: true); // 16 nt (>= 14) => Marmur-Doty

        Assert.That(tm, Is.EqualTo(43.375).Within(1e-9),
            "Length 16 (>= WallaceMaxLength 14) auto-switches to Marmur-Doty: 64.9 + 41*(8-16.4)/16 = 43.375");
    }

    // C3 — Null / empty input to the delegate returns 0 (guarded by string.IsNullOrEmpty).
    [Test]
    public void CalculateMeltingTemperature_NullOrEmptyInput_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateMeltingTemperature(null!), Is.EqualTo(0.0),
                "Null input returns 0");
            Assert.That(SequenceStatistics.CalculateMeltingTemperature(""), Is.EqualTo(0.0),
                "Empty input returns 0");
        });
    }

    #endregion
}
