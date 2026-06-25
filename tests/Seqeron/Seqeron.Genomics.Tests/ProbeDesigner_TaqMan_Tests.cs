// PROBE-DESIGN-001 — TaqMan (hydrolysis-probe) design rules (opt-in)
// Evidence: docs/Evidence/PROBE-DESIGN-001-Evidence.md
// TestSpec: tests/TestSpecs/PROBE-DESIGN-001.md
// Source: Applied Biosystems / Thermo Fisher TaqMan probe-design guidelines;
//         PREMIER Biosoft "TaqMan probe design tips"
//         (http://www.premierbiosoft.com/tech_notes/TaqMan.html), retrieved 2026-06-24.
//   Rules: no G at the 5' end (a 5' G adjacent to the reporter dye quenches reporter
//   fluorescence even after cleavage); more Cs than Gs; no run of >=4 consecutive Gs;
//   G+C content 30-80%; probe Tm ~10 degC above the primer Tm; length 18-22 nt.
using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using System;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Tests for the opt-in TaqMan hydrolysis-probe rules on <see cref="ProbeDesigner"/>:
/// <c>EvaluateTaqManProbe</c> and <c>SelectTaqManStrand</c>. The generic designer is unchanged.
/// Test Unit: PROBE-DESIGN-001
/// </summary>
[TestFixture]
public class ProbeDesigner_TaqMan_Tests
{
    #region Test Data (hand-derived expected values)

    // A — satisfies every rule. len=18, 5'=C, C=10/G=0, maxGrun=0, gc=10/18=0.5556 (30-80%),
    // salt-adjusted Tm = 81.5 + 16.6*log10(0.05) + 41*0.5556 - 600/18 = 49.3473 degC.
    private const string PassAllProbe = "CCATCACCCTACATCACC";
    private const double PassAllTm = 49.3473; // hand-computed, see header

    // B — identical to A but a guanine at the 5' end -> violates the no-5'-G rule only.
    private const string FivePrimeGProbe = "GCATCACCCTACATCACC";

    // C — run of four consecutive Gs -> violates the no->=4-G-run rule only. 5'=A, C=8/G=4, gc=0.6667.
    private const string FourGRunProbe = "ACCCCGGGGACCCTACAT";

    // D — more G than C (C=1, G=9) -> violates the more-C-than-G rule. 5'=A, maxGrun=3.
    private const string MoreGThanCProbe = "ACGGGAGGTAGGTAGGTA";

    // E — 16 nt -> violates the 18-22 length rule.
    private const string ShortProbe = "CCATCACCCTACATCA";

    // F — G+C content = 100% -> violates the 30-80% GC rule.
    private const string AllGcProbe = "CCCGCCCCGCCCCGCCCC";

    #endregion

    #region EvaluateTaqManProbe — individual rules (Must)

    [Test]
    public void EvaluateTaqManProbe_FivePrimeGuanine_FlaggedAndRejected()
    {
        // TM1: a 5' G quenches the reporter dye even after cleavage -> NoGuanineAt5Prime is false,
        // and the probe must not pass overall.
        var e = ProbeDesigner.EvaluateTaqManProbe(FivePrimeGProbe);

        Assert.Multiple(() =>
        {
            Assert.That(e.NoGuanineAt5Prime, Is.False,
                "A 5' guanine must be flagged: it quenches the reporter dye even after cleavage.");
            Assert.That(e.PassesAll, Is.False,
                "A probe with a 5' G cannot pass the TaqMan rule set.");
            Assert.That(e.Violations, Has.Some.Contains("5' end"),
                "The 5'-G violation must be recorded.");
        });
    }

    [Test]
    public void EvaluateTaqManProbe_RunOfFourGuanines_Rejected()
    {
        // TM2: "no runs of identical nucleotides (especially four or more consecutive Gs)".
        // FourGRunProbe = ACCCCGGGGACCCTACAT has GGGG (run of 4) -> NoRunOfFourOrMoreG false.
        var e = ProbeDesigner.EvaluateTaqManProbe(FourGRunProbe);

        Assert.Multiple(() =>
        {
            Assert.That(e.NoRunOfFourOrMoreG, Is.False,
                "A run of four consecutive Gs must be flagged.");
            Assert.That(e.NoGuanineAt5Prime, Is.True, "5' base is A — the 5'-G rule is satisfied.");
            Assert.That(e.MoreCytosineThanGuanine, Is.True, "C=8 > G=4, so the C>G rule holds.");
            Assert.That(e.PassesAll, Is.False, "The G-run violation prevents passing.");
        });
    }

    [Test]
    public void EvaluateTaqManProbe_MoreGuanineThanCytosine_FlagsCgRule()
    {
        // TM3: probe must have more Cs than Gs. MoreGThanCProbe = ACGGGAGGTAGGTAGGTA: C=1, G=9.
        var e = ProbeDesigner.EvaluateTaqManProbe(MoreGThanCProbe);

        Assert.Multiple(() =>
        {
            Assert.That(e.CytosineCount, Is.EqualTo(1), "Exactly one C in ACGGGAGGTAGGTAGGTA.");
            Assert.That(e.GuanineCount, Is.EqualTo(9), "Exactly nine Gs in ACGGGAGGTAGGTAGGTA.");
            Assert.That(e.MoreCytosineThanGuanine, Is.False, "C=1 is not > G=9.");
            Assert.That(e.PassesAll, Is.False, "Failing the C>G rule prevents passing.");
        });
    }

    [Test]
    public void EvaluateTaqManProbe_LengthOutsideRange_FlagsLengthRule()
    {
        // TM4: probe length must be 18-22 nt. ShortProbe is 16 nt.
        var e = ProbeDesigner.EvaluateTaqManProbe(ShortProbe);

        Assert.Multiple(() =>
        {
            Assert.That(e.LengthInRange, Is.False, "16 nt is outside the 18-22 nt range.");
            Assert.That(e.PassesAll, Is.False, "Failing the length rule prevents passing.");
        });
    }

    [Test]
    public void EvaluateTaqManProbe_GcContentOutsideRange_FlagsGcRule()
    {
        // TM5: G+C content must be 30-80%. AllGcProbe has GC = 100%.
        var e = ProbeDesigner.EvaluateTaqManProbe(AllGcProbe);

        Assert.Multiple(() =>
        {
            Assert.That(e.GcContent, Is.EqualTo(1.0).Within(1e-10), "All G/C -> GC fraction 1.0.");
            Assert.That(e.GcContentInRange, Is.False, "100% GC is outside 30-80%.");
            Assert.That(e.PassesAll, Is.False, "Failing the GC rule prevents passing.");
        });
    }

    [Test]
    public void EvaluateTaqManProbe_ProbeTmNotTenAbovePrimer_FlagsTmGate()
    {
        // TM6: probe Tm must be >= primerTm + 10 degC. Probe A Tm = 49.3473.
        // primerTm = 45 -> 45 + 10 = 55 > 49.35 -> gate fails.
        var e = ProbeDesigner.EvaluateTaqManProbe(PassAllProbe, primerTm: 45.0);

        Assert.Multiple(() =>
        {
            Assert.That(e.Tm, Is.EqualTo(PassAllTm).Within(1e-3),
                "Salt-adjusted Tm of CCATCACCCTACATCACC must be 49.3473 degC.");
            Assert.That(e.ProbeTmAbovePrimer, Is.False,
                "Probe Tm 49.35 is below primer Tm 45 + 10 = 55.");
            Assert.That(e.PassesAll, Is.False, "Failing the Tm gate prevents passing.");
        });
    }

    [Test]
    public void EvaluateTaqManProbe_AllRulesSatisfied_Accepted()
    {
        // TM7: probe A satisfies every rule. With primerTm = 38, gate needs Tm >= 48; Tm = 49.35 -> pass.
        var e = ProbeDesigner.EvaluateTaqManProbe(PassAllProbe, primerTm: 38.0);

        Assert.Multiple(() =>
        {
            Assert.That(e.NoGuanineAt5Prime, Is.True, "5' base is C.");
            Assert.That(e.MoreCytosineThanGuanine, Is.True, "C=10 > G=0.");
            Assert.That(e.NoRunOfFourOrMoreG, Is.True, "No G run (G=0).");
            Assert.That(e.GcContentInRange, Is.True, "GC = 0.5556 is within 30-80%.");
            Assert.That(e.LengthInRange, Is.True, "Length 18 is within 18-22.");
            Assert.That(e.ProbeTmAbovePrimer, Is.True, "Tm 49.35 >= primer 38 + 10 = 48.");
            Assert.That(e.PassesAll, Is.True, "All six TaqMan rules are satisfied.");
            Assert.That(e.Violations, Is.Empty, "No violations recorded for a fully compliant probe.");
        });
    }

    [Test]
    public void EvaluateTaqManProbe_NoPrimerTm_TmGateReportedSatisfied()
    {
        // TM8: when no primer Tm is supplied, the probe-Tm gate is skipped (reported satisfied).
        var e = ProbeDesigner.EvaluateTaqManProbe(PassAllProbe, primerTm: null);

        Assert.Multiple(() =>
        {
            Assert.That(e.ProbeTmAbovePrimer, Is.True, "Gate skipped when primerTm is null.");
            Assert.That(e.PassesAll, Is.True, "All other rules pass and the gate is skipped.");
        });
    }

    [Test]
    public void EvaluateTaqManProbe_NullSequence_Throws()
    {
        // Edge: null input is a contract violation.
        Assert.Throws<ArgumentNullException>(() => ProbeDesigner.EvaluateTaqManProbe(null!));
    }

    #endregion

    #region SelectTaqManStrand (Must)

    [Test]
    public void SelectTaqManStrand_SenseHas5PrimeGAndMoreG_PicksAntisense()
    {
        // TM9: sense = GTTAGGGTTAGGGTTAGG (5'=G, C=0, G=9) fails the two hard reporter-dye rules.
        // Reverse complement = CCTAACCCTAACCCTAAC (5'=C, C=9, G=0) satisfies them.
        // Strand selection must choose the antisense strand.
        const string sense = "GTTAGGGTTAGGGTTAGG";
        const string expectedAntisense = "CCTAACCCTAACCCTAAC";

        var (probe, isRc, eval) = ProbeDesigner.SelectTaqManStrand(sense);

        Assert.Multiple(() =>
        {
            Assert.That(isRc, Is.True, "Antisense strand must be selected when the sense strand has a 5' G.");
            Assert.That(probe, Is.EqualTo(expectedAntisense), "Chosen probe is the reverse complement.");
            Assert.That(eval.NoGuanineAt5Prime, Is.True, "Chosen strand has no 5' G.");
            Assert.That(eval.MoreCytosineThanGuanine, Is.True, "Chosen strand has more C than G.");
        });
    }

    [Test]
    public void SelectTaqManStrand_SenseAlreadyCompliant_KeepsSense()
    {
        // TM10: probe A is already TaqMan-compliant -> sense strand kept (no RC).
        var (probe, isRc, eval) = ProbeDesigner.SelectTaqManStrand(PassAllProbe);

        Assert.Multiple(() =>
        {
            Assert.That(isRc, Is.False, "Compliant sense strand should be kept.");
            Assert.That(probe, Is.EqualTo(PassAllProbe), "Sense sequence is returned unchanged.");
            Assert.That(eval.NoGuanineAt5Prime, Is.True, "Sense strand has no 5' G.");
            Assert.That(eval.MoreCytosineThanGuanine, Is.True, "Sense strand has more C than G.");
        });
    }

    [Test]
    public void SelectTaqManStrand_NeitherStrandPasses_PicksHigherRankedByHardRules()
    {
        // TM11: ranking-fallback path (private RankTaqManStrand). Both strands fail (16 nt is
        // outside 18-22, so PassesAll is false on both), forcing the rank-based tie-break.
        // sense = GATCACCCTACATCAC : 5'=G (fails no-5'-G), C=7 > G=1 (more-C-than-G holds),
        //   no >=4-G run, GC=0.5 in range -> rank = 0 + 2 + 1 + 1 + 1(Tm gate skipped) = 5.
        // antisense = GTGATGTAGGGTGATC : 5'=G (fails), C=1 < G=7 (fails more-C), no >=4-G run,
        //   GC=0.5 in range -> rank = 0 + 0 + 1 + 1 + 1 = 3.
        // The more-C-than-G strand (sense, rank 5) must win even though neither passes outright.
        const string sense = "GATCACCCTACATCAC";

        var (probe, isRc, eval) = ProbeDesigner.SelectTaqManStrand(sense);

        Assert.Multiple(() =>
        {
            Assert.That(eval.PassesAll, Is.False, "Neither 16-nt strand passes the length rule.");
            Assert.That(isRc, Is.False, "Sense ranks higher (more C than G) and must be chosen.");
            Assert.That(probe, Is.EqualTo(sense), "Higher-ranked sense strand returned unchanged.");
            Assert.That(eval.MoreCytosineThanGuanine, Is.True, "Chosen strand has more C than G (C=7 > G=1).");
        });
    }

    #endregion

    #region EvaluateTaqManProbe — boundary pass sides (Should)

    [Test]
    public void EvaluateTaqManProbe_ThreeGuanineRun_PassesNoFourGRule()
    {
        // TM12: the no->=4-G-run rule is "< 4"; a run of exactly three Gs must PASS it.
        // GCC GGG ... -> wait, keep 5'!=G. ACCGGGACCCTACATCACC trimmed to 18: ACCGGGACCCTACATCAC.
        // 5'=A, run of 3 Gs (GGG), C=8 > G=3, len=18, GC=11/18=0.611 in range.
        const string threeGRun = "ACCGGGACCCTACATCAC";
        var e = ProbeDesigner.EvaluateTaqManProbe(threeGRun);

        Assert.Multiple(() =>
        {
            Assert.That(e.NoRunOfFourOrMoreG, Is.True, "A run of exactly three Gs is allowed (rule is < 4).");
            Assert.That(e.PassesAll, Is.True, "All other rules also hold for this probe.");
        });
    }

    [Test]
    public void EvaluateTaqManProbe_EmptySequence_FailsAllStructuralRules()
    {
        // TM13: defined behaviour for the empty-string edge (no throw; structural rules fail).
        var e = ProbeDesigner.EvaluateTaqManProbe("");

        Assert.Multiple(() =>
        {
            Assert.That(e.NoGuanineAt5Prime, Is.False, "Empty sequence has no valid 5' base -> rule not satisfied.");
            Assert.That(e.MoreCytosineThanGuanine, Is.False, "C=0 is not > G=0.");
            Assert.That(e.LengthInRange, Is.False, "Length 0 is outside 18-22.");
            Assert.That(e.GcContentInRange, Is.False, "GC=0 is below 30%.");
            Assert.That(e.PassesAll, Is.False, "An empty probe cannot pass.");
        });
    }

    #endregion
}
