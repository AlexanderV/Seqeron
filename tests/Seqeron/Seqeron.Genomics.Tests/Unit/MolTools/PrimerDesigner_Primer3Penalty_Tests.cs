// PRIMER-TM-001 — Primer3 weighted per-primer penalty (objective function)
// Evidence: docs/Evidence/PRIMER-TM-001-Evidence.md
// TestSpec: tests/TestSpecs/PRIMER-TM-001-Penalty.md
// Source: Primer3 source libprimer3.cc p_obj_fn / pr_set_default_global_args_2 (branch main);
//         Primer3 manual §19; Untergasser et al. (2012) NAR 40(15):e115; Koressaar & Remm (2007).
namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Canonical tests for the Primer3 weighted penalty objective (PRIMER-TM-001).
/// Expected values are computed by hand from the Primer3 formula
/// penalty = WT_TM·|Tm−OPT_TM| + WT_SIZE·|len−OPT_SIZE| + WT_GC·|GC−OPT_GC|
///           + WT_SELF_ANY·selfAny + WT_SELF_END·selfEnd + WT_NUM_NS·N
/// with the one-sided weights/optima documented in libprimer3.cc and the Primer3 manual.
/// </summary>
[TestFixture]
public class PrimerDesigner_Primer3Penalty_Tests
{
    private const double Tol = 1e-10;

    #region CalculatePrimer3Penalty — default weights/optima

    // M1 — At the optimum (Tm=60, len=20, GC=50, no N) every term is 0.
    [Test]
    public void CalculatePrimer3Penalty_AtOptimum_ReturnsZero()
    {
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 50.0));
        Assert.That(p, Is.EqualTo(0.0).Within(Tol),
            "All parameters at their Primer3 optima must yield penalty 0 (libprimer3.cc p_obj_fn).");
    }

    // M2 — Tm only above optimum: WT_TM_GT(1)·(63−60) = 3; size term 0.
    [Test]
    public void CalculatePrimer3Penalty_TmAboveOptimum_AddsTmGtTerm()
    {
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 63.0, Length: 20, GcPercent: 50.0));
        Assert.That(p, Is.EqualTo(3.0).Within(Tol),
            "Default WT_TM_GT=1, OPT_TM=60: penalty = 1*(63-60) = 3.0.");
    }

    // M3 — Tm below and size below optimum: 1*(60-57) + 1*(20-18) = 5.
    [Test]
    public void CalculatePrimer3Penalty_TmAndSizeBelowOptimum_SumsBothTerms()
    {
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 57.0, Length: 18, GcPercent: 50.0));
        Assert.That(p, Is.EqualTo(5.0).Within(Tol),
            "WT_TM_LT*(60-57)=3 plus WT_SIZE_LT*(20-18)=2 => 5.0.");
    }

    // M4 — Fractional Tm + size above optimum: 1*(62.5-60) + 1*(22-20) = 4.5.
    [Test]
    public void CalculatePrimer3Penalty_FractionalTmAndSizeAbove_SumsTerms()
    {
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 62.5, Length: 22, GcPercent: 50.0));
        Assert.That(p, Is.EqualTo(4.5).Within(Tol),
            "1*(62.5-60)=2.5 plus 1*(22-20)=2 => 4.5.");
    }

    #endregion

    #region CalculatePrimer3Penalty — GC, self-complementarity, N terms (non-default weights)

    // M5 — GC above optimum with WT_GC_PERCENT_GT=0.5: 0.5*(60-50) = 5.0.
    [Test]
    public void CalculatePrimer3Penalty_GcAboveOptimum_AddsGcGtTerm()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { GcGt = 0.5 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 60.0), w);
        Assert.That(p, Is.EqualTo(5.0).Within(Tol),
            "WT_GC_GT=0.5, OPT_GC=50: penalty = 0.5*(60-50) = 5.0.");
    }

    // M6 — GC below optimum with WT_GC_PERCENT_LT=0.5: 0.5*(50-40) = 5.0.
    [Test]
    public void CalculatePrimer3Penalty_GcBelowOptimum_AddsGcLtTerm()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { GcLt = 0.5 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 40.0), w);
        Assert.That(p, Is.EqualTo(5.0).Within(Tol),
            "WT_GC_LT=0.5, OPT_GC=50: penalty = 0.5*(50-40) = 5.0.");
    }

    // M7 — Self-any term scales linearly: WT_SELF_ANY(0.1)*4 = 0.4.
    [Test]
    public void CalculatePrimer3Penalty_SelfAny_AddsLinearTerm()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { SelfAny = 0.1 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 50.0, SelfAny: 4.0), w);
        Assert.That(p, Is.EqualTo(0.4).Within(Tol),
            "WT_SELF_ANY=0.1, selfAny=4: penalty = 0.1*4 = 0.4 (p_obj_fn compl_any).");
    }

    // M8 — Self-end term scales linearly: WT_SELF_END(0.2)*3 = 0.6.
    [Test]
    public void CalculatePrimer3Penalty_SelfEnd_AddsLinearTerm()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { SelfEnd = 0.2 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 50.0, SelfEnd: 3.0), w);
        Assert.That(p, Is.EqualTo(0.6).Within(Tol),
            "WT_SELF_END=0.2, selfEnd=3: penalty = 0.2*3 = 0.6 (p_obj_fn compl_end).");
    }

    // M9 — Num-Ns term: WT_NUM_NS(1)*2 = 2.0.
    [Test]
    public void CalculatePrimer3Penalty_NumNs_AddsLinearTerm()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { NumNs = 1.0 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 50.0, NumNs: 2), w);
        Assert.That(p, Is.EqualTo(2.0).Within(Tol),
            "WT_NUM_NS=1, N=2: penalty = 1*2 = 2.0 (p_obj_fn num_ns).");
    }

    // M10 — Combined: TM_GT*2 + SIZE_GT*2 + GC_GT(0.5)*5 + SELF_ANY(0.25)*2 + NUM_NS*1
    //               = 2 + 2 + 2.5 + 0.5 + 1 = 8.0.
    [Test]
    public void CalculatePrimer3Penalty_CombinedTerms_SumsAllContributions()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { GcGt = 0.5, SelfAny = 0.25, NumNs = 1.0 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 62.0, Length: 22, GcPercent: 55.0, SelfAny: 2.0, NumNs: 1), w);
        Assert.That(p, Is.EqualTo(8.0).Within(Tol),
            "1*(62-60)+1*(22-20)+0.5*(55-50)+0.25*2+1*1 = 2+2+2.5+0.5+1 = 8.0.");
    }

    #endregion

    #region Default weights/optima constants

    // M11 — Default weights and optima must match the Primer3 source / manual values.
    [Test]
    public void DefaultPrimer3WeightsAndOptima_MatchPrimer3Defaults()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights;
        var o = PrimerDesigner.DefaultPrimer3Optima;
        Assert.Multiple(() =>
        {
            Assert.That(w.TmGt, Is.EqualTo(1.0).Within(Tol), "PRIMER_WT_TM_GT default = 1 (libprimer3.cc temp_gt).");
            Assert.That(w.TmLt, Is.EqualTo(1.0).Within(Tol), "PRIMER_WT_TM_LT default = 1 (temp_lt).");
            Assert.That(w.SizeGt, Is.EqualTo(1.0).Within(Tol), "PRIMER_WT_SIZE_GT default = 1 (length_gt).");
            Assert.That(w.SizeLt, Is.EqualTo(1.0).Within(Tol), "PRIMER_WT_SIZE_LT default = 1 (length_lt).");
            Assert.That(w.GcGt, Is.EqualTo(0.0).Within(Tol), "PRIMER_WT_GC_PERCENT_GT default = 0 (gc_content_gt).");
            Assert.That(w.GcLt, Is.EqualTo(0.0).Within(Tol), "PRIMER_WT_GC_PERCENT_LT default = 0 (gc_content_lt).");
            Assert.That(w.SelfAny, Is.EqualTo(0.0).Within(Tol), "PRIMER_WT_SELF_ANY default = 0 (compl_any).");
            Assert.That(w.SelfEnd, Is.EqualTo(0.0).Within(Tol), "PRIMER_WT_SELF_END default = 0 (compl_end).");
            Assert.That(w.NumNs, Is.EqualTo(0.0).Within(Tol), "PRIMER_WT_NUM_NS default = 0 (num_ns).");
            Assert.That(o.OptTm, Is.EqualTo(60.0).Within(Tol), "PRIMER_OPT_TM default = 60.0 (opt_tm).");
            Assert.That(o.OptSize, Is.EqualTo(20), "PRIMER_OPT_SIZE default = 20 (opt_size).");
            Assert.That(o.OptGcPercent, Is.EqualTo(50.0).Within(Tol), "PRIMER_OPT_GC_PERCENT default = 50.0 (manual).");
        });
    }

    #endregion

    #region Invariants and edge cases

    // S1 — Tm exactly at optimum: neither _gt nor _lt fires (strict gates), so 0 even with large weights.
    [Test]
    public void CalculatePrimer3Penalty_TmExactlyAtOptimum_ContributesZero()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { TmGt = 5.0, TmLt = 5.0 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 50.0), w);
        Assert.That(p, Is.EqualTo(0.0).Within(Tol),
            "At Tm=OPT_TM the strict >/< gates exclude both terms => 0 (INV-03).");
    }

    // S2 — Asymmetric Tm weights: above optimum uses WT_TM_GT only (2*(62-60)=4), not WT_TM_LT.
    [Test]
    public void CalculatePrimer3Penalty_AsymmetricTmWeights_UsesGtSideOnly()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { TmGt = 2.0, TmLt = 1.0 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 62.0, Length: 20, GcPercent: 50.0), w);
        Assert.That(p, Is.EqualTo(4.0).Within(Tol),
            "Tm>OPT uses WT_TM_GT=2: 2*(62-60)=4.0; the _lt weight must not apply (one-sided).");
    }

    // S3 — Linearity in the weight: doubling WT_SELF_ANY doubles that term.
    [Test]
    public void CalculatePrimer3Penalty_DoublingWeight_DoublesTerm()
    {
        var inputs = new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 50.0, SelfAny: 4.0);
        var p1 = PrimerDesigner.CalculatePrimer3Penalty(inputs, PrimerDesigner.DefaultPrimer3Weights with { SelfAny = 0.1 });
        var p2 = PrimerDesigner.CalculatePrimer3Penalty(inputs, PrimerDesigner.DefaultPrimer3Weights with { SelfAny = 0.2 });
        Assert.Multiple(() =>
        {
            Assert.That(p1, Is.EqualTo(0.4).Within(Tol), "0.1*4 = 0.4.");
            Assert.That(p2, Is.EqualTo(0.8).Within(Tol), "0.2*4 = 0.8 = 2*p1 (INV-04 linearity).");
        });
    }

    // S4 — Non-negativity across a battery of deviating inputs (INV-01).
    [Test]
    public void CalculatePrimer3Penalty_VariousDeviations_AlwaysNonNegative()
    {
        var w = new Primer3PenaltyWeights(
            TmGt: 1, TmLt: 1, SizeGt: 1, SizeLt: 1, GcGt: 0.5, GcLt: 0.5,
            SelfAny: 0.1, SelfEnd: 0.1, NumNs: 1);
        var cases = new[]
        {
            new Primer3PenaltyInputs(40.0, 12, 10.0, 8.0, 6.0, 3),
            new Primer3PenaltyInputs(80.0, 30, 90.0, 0.0, 0.0, 0),
            new Primer3PenaltyInputs(60.0, 20, 50.0, 0.0, 0.0, 0),
            new Primer3PenaltyInputs(59.9, 19, 49.0, 1.0, 0.0, 1),
        };
        Assert.Multiple(() =>
        {
            foreach (var c in cases)
                Assert.That(PrimerDesigner.CalculatePrimer3Penalty(c, w), Is.GreaterThanOrEqualTo(0.0),
                    "Every penalty term is weight*non-negative deviation, so the total is >= 0 (INV-01).");
        });
    }

    // C1 — Lower penalty selects the better candidate: near-optimal primer scores below a far-off one.
    [Test]
    public void CalculatePrimer3Penalty_NearOptimalCandidate_HasLowerPenalty()
    {
        var near = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.5, Length: 20, GcPercent: 50.0)); // 0.5
        var far = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 66.0, Length: 25, GcPercent: 50.0)); // 6 + 5 = 11
        Assert.Multiple(() =>
        {
            Assert.That(near, Is.EqualTo(0.5).Within(Tol), "Near-optimal: 1*(60.5-60)=0.5.");
            Assert.That(far, Is.EqualTo(11.0).Within(Tol), "Far: 1*(66-60)+1*(25-20)=6+5=11.");
            Assert.That(near, Is.LessThan(far), "Primer3 selects the lowest-penalty primer (lower is better).");
        });
    }

    // C2 — GC is interpreted as a percentage 0-100: GC=50 gives 0 even with a non-zero GC weight.
    [Test]
    public void CalculatePrimer3Penalty_GcAsPercent_FiftyPercentIsOptimum()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { GcGt = 1.0, GcLt = 1.0 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 50.0), w);
        Assert.That(p, Is.EqualTo(0.0).Within(Tol),
            "GC% is a percentage (0-100); GC=50 equals OPT_GC=50 so the GC term is 0 (libprimer3.cc gc_content = 100*num_gc/num_gcat).");
    }

    #endregion

    #region Mutation killers — custom optima and one-sided weight gating

    // K1 — A CUSTOM optima must actually be used (not silently replaced by the defaults). With OptTm shifted
    // to 65, a 60 °C primer is now 5 below optimum, so WT_TM_LT*(65-60)=5; under the default optima (60) it
    // would be 0. Pins the "optima ?? DefaultPrimer3Optima" null-coalescing against dropping the left operand.
    [Test]
    public void CalculatePrimer3Penalty_CustomOptima_UsesProvidedOptimaNotDefaults()
    {
        var optima = new Primer3Optima(OptTm: 65.0, OptSize: 24, OptGcPercent: 50.0);
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 50.0),
            weights: null, optima: optima);

        // WT_TM_LT*(65-60)=5 ; WT_SIZE_LT*(24-20)=4 ; GC at optimum -> 0. Total 9.
        Assert.That(p, Is.EqualTo(9.0).Within(Tol),
            "Custom OPT_TM=65, OPT_SIZE=24: penalty = 1*(65-60) + 1*(24-20) = 9 (defaults would give 0).");
    }

    // K2 — A non-zero GC_GT weight must NOT penalise a BELOW-optimum GC: the gating is
    // "weight != 0 AND value > optimum" (logical AND, not OR). GC=40 < 50 with WT_GC_GT=0.5 and WT_GC_LT=0
    // contributes nothing; an OR-gating mutant would wrongly add 0.5*(40-50) = -5.
    [Test]
    public void CalculatePrimer3Penalty_GcGtWeight_DoesNotFireBelowOptimum()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { GcGt = 0.5, GcLt = 0.0 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 40.0), w);

        Assert.That(p, Is.EqualTo(0.0).Within(Tol),
            "GC_GT only penalises GC ABOVE optimum; GC=40 < 50 leaves the GC term at 0 (p_obj_fn gc_content_gt).");
    }

    // K3 — Symmetric check for the LT side: a non-zero GC_LT weight must NOT fire for an ABOVE-optimum GC.
    [Test]
    public void CalculatePrimer3Penalty_GcLtWeight_DoesNotFireAboveOptimum()
    {
        var w = PrimerDesigner.DefaultPrimer3Weights with { GcLt = 0.5, GcGt = 0.0 };
        var p = PrimerDesigner.CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(Tm: 60.0, Length: 20, GcPercent: 60.0), w);

        Assert.That(p, Is.EqualTo(0.0).Within(Tol),
            "GC_LT only penalises GC BELOW optimum; GC=60 > 50 leaves the GC term at 0 (p_obj_fn gc_content_lt).");
    }

    #endregion
}
