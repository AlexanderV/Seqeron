// PROBE-DESIGN-001 — LNA-adjusted nearest-neighbour Tm + citable MGB design rules
// Evidence: docs/Evidence/PROBE-DESIGN-001-LNA-Evidence.md
// TestSpec: tests/TestSpecs/PROBE-DESIGN-001-LNA.md
// Sources:
//   McTigue PM, Peterson RJ, Kahn JD (2004). Biochemistry 43:5388-5405 — LNA-DNA NN increments
//     (DOI 10.1021/bi035976d); values transcribed verbatim from MELTING 5 McTigue2004lockedmn.xml.
//   rmelting tutorial worked example (MELTING mct04): CCATT(L)GCTACC → Tm 63.61426 °C.
//   SantaLucia J (1998). PNAS 95(4):1460-65 — base DNA NN model.
//   Kutyavin IV et al. (2000). Nucleic Acids Res 28(2):655-661 — 3'-MGB design rules.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the opt-in LNA (locked nucleic acid)-adjusted nearest-neighbour Tm added under
/// PROBE-DESIGN-001, and the citable 3'-MGB design-rule check. Expected ΔH°/ΔS°/Tm are
/// hand-derived from the McTigue (2004) increment table (verbatim from the MELTING data file)
/// added to the SantaLucia (1998) DNA NN model — independent of the implementation.
/// </summary>
[TestFixture]
public class ProbeDesigner_LnaTm_Tests
{
    private const double Tol = 1e-9;
    private const double TmTol = 1e-4;

    // Worked example duplex (rmelting tutorial CCATT(L)GCTACC → DNA CCATTGCTACC, LNA at index 4).
    private const string WorkedSeq = "CCATTGCTACC";
    private const int WorkedLnaIndex = 4;

    // Worked-example conditions: C_T = 1e-4 M, [Na+] = 1 M reference state, no salt correction.
    private const double WorkedConc = 1e-4;
    private const double WorkedNa = 1.0;

    #region CalculateNearestNeighborThermodynamicsLna — ΔH°/ΔS°

    // M1 — CCATTGCTACC, LNA at index 4 (the second T).
    // Base DNA NN (SantaLucia 1998 unified, library NnUnifiedParams):
    //   ends C,C → no terminal-AT; not self-comp → no symmetry.
    //   ΔH° = -80.8 kcal/mol, ΔS° = -221.7 cal/(mol·K)  (library CalculateNearestNeighborThermodynamics).
    // LNA increments for locked index 4 (McTigue 2004, verbatim XML):
    //   step i=3 "TT", 3' base (idx 4) locked → TTL/AA = (+2326 cal, +8.1) = (+2.326 kcal, +8.1)
    //   step i=4 "TG", 5' base (idx 4) locked → TLG/AC = (-1540 cal, -3.0) = (-1.540 kcal, -3.0)
    //   ΔH°_LNA = -80.8 + 2.326 - 1.540 = -80.014 kcal/mol
    //   ΔS°_LNA = -221.7 + 8.1 - 3.0    = -216.6  cal/(mol·K)
    [Test]
    public void CalculateNearestNeighborThermodynamicsLna_SingleInternalLna_MatchesMcTigueIncrement()
    {
        var lna = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, new[] { WorkedLnaIndex });

        Assert.That(lna, Is.Not.Null, "CCATTGCTACC is all-ACGT and the LNA at index 4 is internal.");
        Assert.Multiple(() =>
        {
            Assert.That(lna!.Value.DeltaH, Is.EqualTo(-80.014).Within(Tol),
                "ΔH° = base DNA -80.8 + TTL/AA(+2.326) + TLG/AC(-1.540) (McTigue 2004 verbatim).");
            Assert.That(lna.Value.DeltaS, Is.EqualTo(-216.6).Within(Tol),
                "ΔS° = base DNA -221.7 + TTL/AA(+8.1) + TLG/AC(-3.0) (McTigue 2004 verbatim).");
            Assert.That(lna.Value.IsSelfComplementary, Is.False,
                "CCATTGCTACC is not self-complementary (computed on the underlying DNA).");
        });
    }

    // M6 — Negative-ΔΔ increment applied with the correct sign.
    // Step "GG" with the 5' base (the first G) locked → GLG/CC = (-2844 cal, -6.7) = (-2.844 kcal, -6.7).
    // Sequence GGGCC, LNA at index 1: only step i=1 "GG" has the locked 5' base (index 1);
    //   step i=0 "GG" has the locked base at i+1 (index 1) → GGL/CC = (-943 cal, -0.9) = (-0.943, -0.9);
    //   step i=2 "GC", step i=3 "CC" no LNA.
    // So increments = GGL/CC (i=0, 3'-locked) + GLG/CC (i=1, 5'-locked).
    //   ΔΔH = -0.943 + -2.844 = -3.787 kcal/mol ; ΔΔS = -0.9 + -6.7 = -7.6 cal/(mol·K).
    [Test]
    public void CalculateNearestNeighborThermodynamicsLna_NegativeIncrement_LowersDeltaHByExactAmount()
    {
        var dna = PrimerDesigner.CalculateNearestNeighborThermodynamics("GGGCC");
        var lna = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("GGGCC", new[] { 1 });

        Assert.That(dna, Is.Not.Null);
        Assert.That(lna, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(lna!.Value.DeltaH - dna!.Value.DeltaH, Is.EqualTo(-3.787).Within(Tol),
                "ΔΔH = GGL/CC(-0.943) + GLG/CC(-2.844) applied with correct (negative) sign.");
            Assert.That(lna.Value.DeltaS - dna.Value.DeltaS, Is.EqualTo(-7.6).Within(Tol),
                "ΔΔS = GGL/CC(-0.9) + GLG/CC(-6.7).");
        });
    }

    // C1 — Two internal LNA positions: increments are additive.
    // CCATTGCTACC, LNA at indices 4 and 6 (T and C).
    //   index 4 (as M1): TTL/AA(i=3,3') + TLG/AC(i=4,5')  → already counted.
    //   index 6 (the C): step i=5 "GC", 3' base (idx 6) locked → GCL/CG = (-0.925, -1.1)
    //                    step i=6 "CT", 5' base (idx 6) locked → CLT/GA = (+0.708, +4.2)
    //   total over single-LNA(4): add GCL/CG + CLT/GA.
    [Test]
    public void CalculateNearestNeighborThermodynamicsLna_TwoInternalLna_AddsBothIncrements()
    {
        var one = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, new[] { 4 });
        var two = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, new[] { 4, 6 });

        Assert.That(one, Is.Not.Null);
        Assert.That(two, Is.Not.Null);
        // ΔΔ from adding the index-6 LNA = GCL/CG(-0.925,-1.1) + CLT/GA(+0.708,+4.2).
        Assert.Multiple(() =>
        {
            Assert.That(two!.Value.DeltaH - one!.Value.DeltaH, Is.EqualTo(-0.925 + 0.708).Within(Tol),
                "Second LNA adds GCL/CG(-0.925) + CLT/GA(+0.708) on top of the first.");
            Assert.That(two.Value.DeltaS - one.Value.DeltaS, Is.EqualTo(-1.1 + 4.2).Within(Tol),
                "Second LNA adds GCL/CG(-1.1) + CLT/GA(+4.2) on top of the first.");
        });
    }

    #endregion

    #region CalculateMeltingTemperatureNNLna — Tm

    // M2 — LNA-adjusted Tm of CCATT(L)GCTACC at C=1e-4, Na=1 (reference state, no salt correction).
    //   Tm = ΔH°·1000 / (ΔS° + R·ln(C/4)) - 273.15,  R = 1.9872, x = 4 (non-self-comp).
    //      = -80014 / (-216.6 + 1.9872·ln(1e-4/4)) - 273.15 = 63.527594 °C.
    //   MELTING mct04 reports 63.61426 °C; agreement within ~0.1 °C (different base DNA NN set).
    [Test]
    public void CalculateMeltingTemperatureNNLna_WorkedExample_MatchesHandDerivedAndMelting()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNNLna(
            WorkedSeq, new[] { WorkedLnaIndex },
            strandConcentrationMolar: WorkedConc, sodiumMolar: WorkedNa,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.Multiple(() =>
        {
            Assert.That(tm, Is.EqualTo(63.527594).Within(TmTol),
                "Hand-derived LNA Tm from -80.014 kcal / -216.6 cal·K⁻¹ via the bimolecular Tm equation.");
            Assert.That(Math.Abs(tm - 63.61426), Is.LessThan(0.1),
                "Within 0.1 °C of MELTING mct04 (63.61426 °C); residual is the base DNA NN model choice.");
        });
    }

    // M3 — Adding the internal LNA RAISES Tm vs the all-DNA duplex (McTigue 2004 stabilization).
    //   all-DNA Tm = -80800 / (-221.7 + R·ln(1e-4/4)) - 273.15 = 59.692264 °C.
    [Test]
    public void CalculateMeltingTemperatureNNLna_AddingLna_RaisesTmVsAllDna()
    {
        double tmDna = PrimerDesigner.CalculateMeltingTemperatureNN(
            WorkedSeq, strandConcentrationMolar: WorkedConc, sodiumMolar: WorkedNa,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);
        double tmLna = PrimerDesigner.CalculateMeltingTemperatureNNLna(
            WorkedSeq, new[] { WorkedLnaIndex },
            strandConcentrationMolar: WorkedConc, sodiumMolar: WorkedNa,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.Multiple(() =>
        {
            Assert.That(tmDna, Is.EqualTo(59.692264).Within(TmTol),
                "All-DNA CCATTGCTACC Tm (no LNA) by the same equation.");
            Assert.That(tmLna, Is.GreaterThan(tmDna),
                "An internal LNA substitution stabilises the duplex → higher Tm (McTigue 2004).");
            Assert.That(tmLna - tmDna, Is.EqualTo(63.527594 - 59.692264).Within(TmTol),
                "ΔTm = +3.835 °C from the McTigue increments for this duplex.");
        });
    }

    // M4 — With no LNA positions, the LNA Tm exactly equals the plain perfect-match NN Tm.
    [Test]
    public void CalculateMeltingTemperatureNNLna_NoLnaPositions_EqualsPlainNNTm()
    {
        double plain = PrimerDesigner.CalculateMeltingTemperatureNN(
            WorkedSeq, strandConcentrationMolar: WorkedConc, sodiumMolar: WorkedNa,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);
        double lnaEmpty = PrimerDesigner.CalculateMeltingTemperatureNNLna(
            WorkedSeq, Array.Empty<int>(),
            strandConcentrationMolar: WorkedConc, sodiumMolar: WorkedNa,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);

        Assert.That(lnaEmpty, Is.EqualTo(plain).Within(Tol),
            "Empty LNA-position set must reduce to the unchanged perfect-match NN Tm (opt-in additivity).");
    }

    // M5 — Terminal LNA (index 0 or last) is not parameterised by McTigue (2004) → not computable.
    [Test]
    public void CalculateMeltingTemperatureNNLna_TerminalLna_ReturnsNotComputable()
    {
        int last = WorkedSeq.Length - 1;
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, new[] { 0 }),
                Is.Null, "LNA at index 0 (terminal) has no McTigue parameter.");
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, new[] { last }),
                Is.Null, "LNA at the last index (terminal) has no McTigue parameter.");
            Assert.That(PrimerDesigner.CalculateMeltingTemperatureNNLna(WorkedSeq, new[] { 0 }),
                Is.NaN, "Terminal LNA → Tm NaN.");
            Assert.That(PrimerDesigner.CalculateMeltingTemperatureNNLna(WorkedSeq, new[] { last }),
                Is.NaN, "Terminal LNA → Tm NaN.");
        });
    }

    #endregion

    #region Edge cases / guards

    // S1 — null / empty / single-base sequence are not computable.
    [Test]
    public void CalculateNearestNeighborThermodynamicsLna_NullEmptyShort_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(null!, Array.Empty<int>()),
                Is.Null, "null sequence → null.");
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("", Array.Empty<int>()),
                Is.Null, "empty sequence → null.");
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("A", Array.Empty<int>()),
                Is.Null, "single base (< 2 nt) → null.");
            Assert.Throws<ArgumentNullException>(
                () => PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, null!),
                "null lnaPositions → ArgumentNullException.");
        });
    }

    // S2 — out-of-range LNA index → null; duplicates and unsorted order tolerated.
    [Test]
    public void CalculateNearestNeighborThermodynamicsLna_IndexRangeAndOrder()
    {
        var outOfRange = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("ACGTA", new[] { 7 });
        var dup = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, new[] { 4, 4 });
        var single = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, new[] { 4 });
        var unsorted = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, new[] { 6, 4 });
        var sorted = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(WorkedSeq, new[] { 4, 6 });

        Assert.Multiple(() =>
        {
            Assert.That(outOfRange, Is.Null, "Index beyond the sequence → null.");
            Assert.That(dup!.Value.DeltaH, Is.EqualTo(single!.Value.DeltaH).Within(Tol),
                "Duplicate LNA index counts once (set semantics).");
            Assert.That(unsorted!.Value.DeltaH, Is.EqualTo(sorted!.Value.DeltaH).Within(Tol),
                "LNA-position order does not change the result.");
        });
    }

    // S3 — non-ACGT base makes the underlying DNA NN lookup fail → not computable.
    [Test]
    public void CalculateMeltingTemperatureNNLna_NonAcgtBase_ReturnsNaN()
    {
        double tm = PrimerDesigner.CalculateMeltingTemperatureNNLna("CCANTGCTACC", new[] { 4 });
        Assert.That(tm, Is.NaN, "An 'N' makes the DNA NN lookup fail → NaN.");
    }

    // S4 — completeness: every NN step has both a 5'-locked and a 3'-locked increment so that any
    // internal LNA in an all-ACGT sequence is computable (32 = 16 steps × 2 locked positions).
    [Test]
    public void CalculateNearestNeighborThermodynamicsLna_AllInternalContexts_AreParameterised()
    {
        // A sequence visiting many contexts; any single internal LNA must be computable.
        const string seq = "ACGTACGTACGTACGT";
        Assert.Multiple(() =>
        {
            for (int i = 1; i < seq.Length - 1; i++)
            {
                var r = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { i });
                Assert.That(r, Is.Not.Null,
                    $"Internal LNA at index {i} must have a McTigue increment for both adjacent steps.");
            }
        });
    }

    #endregion

    #region MGB design rules (Kutyavin 2000)

    // M7 — MGB design rules: 12–20mer length window; 3'-MGB attachment guidance.
    [Test]
    public void EvaluateMgbProbeDesign_LengthWindowAndThreePrimePlacement()
    {
        var fifteen = ProbeDesigner.EvaluateMgbProbeDesign("ACGTACGTACGTACG");   // 15 nt — in range
        var twentyFive = ProbeDesigner.EvaluateMgbProbeDesign("ACGTACGTACGTACGTACGTACGTA"); // 25 nt — too long

        Assert.Multiple(() =>
        {
            Assert.That(fifteen.Length, Is.EqualTo(15));
            Assert.That(fifteen.LengthInMgbRange, Is.True,
                "15mer is within the Kutyavin (2000) MGB 12–20mer window.");
            Assert.That(fifteen.MgbAttachmentEnd, Is.EqualTo("3'"),
                "MGB is attached at the 3' end (Kutyavin 2000).");

            Assert.That(twentyFive.Length, Is.EqualTo(25));
            Assert.That(twentyFive.LengthInMgbRange, Is.False,
                "25mer exceeds the MGB 12–20mer window — MGB probes are designed shorter.");
            Assert.That(twentyFive.Guidance,
                Has.Some.Contains("12-20mer"),
                "Out-of-range length is flagged against the cited 12–20mer window.");
        });
    }

    // S-guard — null probe throws.
    [Test]
    public void EvaluateMgbProbeDesign_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => ProbeDesigner.EvaluateMgbProbeDesign(null!));
    }

    #endregion
}
