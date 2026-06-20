using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology MHC peptide–binding classification unit — ONCO-MHC-001.
/// The units under test are
/// <see cref="OncologyAnalyzer.ClassifyBindingAffinity(double)"/> (IC50 nM → BindingStrength),
/// <see cref="OncologyAnalyzer.ClassifyBindingRank(double, OncologyAnalyzer.MhcClass)"/> (%Rank → BindingStrength),
/// <see cref="OncologyAnalyzer.IsValidPeptideLength(int, OncologyAnalyzer.MhcClass)"/> (length range gate), and
/// <see cref="OncologyAnalyzer.ClassifyMhcBinding(int, double, OncologyAnalyzer.MhcClass)"/> (length gate + affinity),
/// all in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// This file is scoped STRICTLY to ONCO-MHC-001 (binding/affinity classification). It does NOT touch the
/// neoantigen peptide-window generator (ONCO-NEO-001, row 109, <see cref="OncologyAnalyzer.GenerateNeoantigenPeptides"/>)
/// nor HLA-allele nomenclature / allele-specific LOH (ONCO-HLA-001, row 117,
/// <see cref="OncologyAnalyzer.ParseHlaAllele"/> / <see cref="OncologyAnalyzer.DetectHlaLoh"/>).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts the code never fails in an
/// undisciplined way: no infinity / NaN leaking through a threshold comparison, no log(0)=−∞ in an affinity
/// transform, no matrix-index overflow on an out-of-range peptide length, no KeyNotFound / NullReference on an
/// unsupported MHC class. Every input must resolve to EITHER a well-defined, theory-correct
/// <see cref="OncologyAnalyzer.BindingStrength"/> OR a *documented, intentional* ArgumentOutOfRangeException.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-MHC-001 — MHC-peptide binding classification (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 110.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty).
///     Targets (checklist row 110): "IC50=0, IC50=∞, peptide too short/long, unknown allele".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Target mapping for this specification-driven (no trained model) unit:
///   • "IC50=0"   = a zero / negative / sub-normal concentration → documented ArgumentOutOfRangeException
///       (INV-01 "IC50 must be finite > 0"); a positive concentration may be arbitrarily small but the
///       extreme tiny IC50 still classifies Strong (no log(0)=−∞ leak — the unit does NO log transform).
///   • "IC50=∞"   = +∞ / −∞ / NaN affinity → documented ArgumentOutOfRangeException (INV-01); Inf must never
///       leak into the threshold comparison and silently become a NonBinder.
///   • "peptide too short / long" = length outside the MHC-class range (class I 8–11, class II 13–25),
///       including 0, −1, int.MinValue, int.MaxValue → IsValidPeptideLength is total and returns false;
///       ClassifyMhcBinding returns NonBinder (the length gate), never an IndexOutOfRange on a scoring matrix.
///   • "unknown allele" = the model has NO learned per-allele set (§5.2: prediction is the caller's job); the
///       only allele-axis the classifier exposes is the <see cref="OncologyAnalyzer.MhcClass"/> enum. An
///       out-of-domain / undefined enum value must not KeyNotFound / NullReference — it falls into the class II
///       branch deterministically (a defined, non-throwing outcome), and a known class behaves per spec.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// MHC_Peptide_Binding_Classification.md (docs/algorithms/Oncology/MHC_Peptide_Binding_Classification.md):
///   • Core model (§2.2): Strong iff v &lt; s; Weak iff s ≤ v &lt; w; NonBinder iff v ≥ w; cutoffs strict `&lt;`.
///   • IC50 cutoffs (§2.2, §4.2): s = 50 nM (Strong), w = 500 nM (Weak). 50 → Weak, 500 → NonBinder.
///   • %Rank cutoffs (§2.2, §4.2): class I s = 0.5, w = 2; class II s = 2, w = 10. Strict `&lt;`.
///   • Length ranges (§2.2, §4.2): class I 8–11 inclusive; class II 13–25 inclusive.
///   • INV-01 (§2.4): IC50 must be finite &gt; 0, else ArgumentOutOfRangeException.
///   • INV-02 (§2.4): %Rank must be in [0, 100], else ArgumentOutOfRangeException.
///   • INV-03 (§2.4): categories are monotone in the score (smaller v ⇒ equal-or-stronger category).
///   • INV-04 (§2.4): cutoffs strict `&lt;`; a value exactly at a cutoff falls in the weaker category.
///   • §3.3: IsValidPeptideLength is total (false for any out-of-range / non-positive length);
///       ClassifyMhcBinding returns NonBinder for invalid length BEFORE evaluating affinity.
///   • §6.1 edge cases: IC50 = 50 → Weak; IC50 = 500 → NonBinder; %Rank = 0.5 (class I) → Weak;
///       %Rank = 2.0 (class I) → NonBinder; IC50 ≤ 0 / NaN / ∞ → throw; length out of range → NonBinder.
///   • §7.1 worked example: ClassifyMhcBinding(9, 42.0, ClassI) == Strong; ClassifyBindingRank(0.5, ClassI) == Weak.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyMhcBindingFuzzTests
{
    // Documented cutoffs / ranges (mirrored from the algorithm doc §2.2/§4.2 and the public consts), used to
    // pin the business contract independently of the production constants where a hand-built example is clearer.
    private const double StrongIc50 = 50.0;   // strict `<` ⇒ Strong
    private const double WeakIc50 = 500.0;    // strict `<` ⇒ Weak; ≥ ⇒ NonBinder

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // A BindingStrength is a 3-valued closed enum; "well formed" means the call returned one of the three
    // defined members (never an undefined cast value, never a leaked Inf-driven default). This is the guard
    // that stops a fuzz test rubber-stamping an out-of-contract result.
    private static void AssertWellFormed(BindingStrength strength)
    {
        strength.Should().BeOneOf(
            BindingStrength.Strong, BindingStrength.Weak, BindingStrength.NonBinder);
    }

    // Independent reference oracle for the IC50 contract (§2.2/§4.2), used to cross-check the implementation
    // over random valid inputs without re-using its branch structure.
    private static BindingStrength ExpectedFromIc50(double ic50Nm)
    {
        if (ic50Nm < StrongIc50) return BindingStrength.Strong;
        if (ic50Nm < WeakIc50) return BindingStrength.Weak;
        return BindingStrength.NonBinder;
    }

    #region ONCO-MHC-001 — Positive sanity (documented thresholds and worked example)

    // §7.1 worked example: a 9-mer with caller-supplied IC50 42 nM on class I is a Strong binder.
    [Test]
    public void ClassifyMhcBinding_DocWorkedExample_9merIc50_42_IsStrong()
    {
        ClassifyMhcBinding(peptideLength: 9, ic50Nm: 42.0, MhcClass.ClassI)
            .Should().Be(BindingStrength.Strong);
    }

    // §7.1 worked example: %Rank exactly 0.5 is NOT < 0.5 ⇒ Weak (strict cutoff, INV-04).
    [Test]
    public void ClassifyBindingRank_DocWorkedExample_ClassIRankHalf_IsWeak()
    {
        ClassifyBindingRank(0.5, MhcClass.ClassI).Should().Be(BindingStrength.Weak);
    }

    // §2.2 IC50 model: a strong, a weak, and a non-binding hand-built affinity classify as documented.
    [Test]
    public void ClassifyBindingAffinity_HandBuiltStrongWeakNonBinder_MapToThresholds()
    {
        ClassifyBindingAffinity(10.0).Should().Be(BindingStrength.Strong);     // 10 < 50
        ClassifyBindingAffinity(200.0).Should().Be(BindingStrength.Weak);      // 50 ≤ 200 < 500
        ClassifyBindingAffinity(5000.0).Should().Be(BindingStrength.NonBinder); // ≥ 500
    }

    // §6.1 / INV-04: the two IC50 cutoffs are strict `<` — boundary values fall to the weaker category.
    [Test]
    public void ClassifyBindingAffinity_AtCutoffs_FallToWeakerCategory()
    {
        ClassifyBindingAffinity(50.0).Should().Be(BindingStrength.Weak);       // not < 50
        ClassifyBindingAffinity(500.0).Should().Be(BindingStrength.NonBinder); // not < 500
    }

    // §4.2 %Rank cutoffs differ per class; verify both class I (0.5/2) and class II (2/10) boundaries.
    [Test]
    public void ClassifyBindingRank_PerClassCutoffs_AreStrict()
    {
        // class I: strong < 0.5, weak < 2
        ClassifyBindingRank(0.4, MhcClass.ClassI).Should().Be(BindingStrength.Strong);
        ClassifyBindingRank(0.5, MhcClass.ClassI).Should().Be(BindingStrength.Weak);
        ClassifyBindingRank(2.0, MhcClass.ClassI).Should().Be(BindingStrength.NonBinder);
        // class II: strong < 2, weak < 10
        ClassifyBindingRank(1.9, MhcClass.ClassII).Should().Be(BindingStrength.Strong);
        ClassifyBindingRank(2.0, MhcClass.ClassII).Should().Be(BindingStrength.Weak);
        ClassifyBindingRank(10.0, MhcClass.ClassII).Should().Be(BindingStrength.NonBinder);
    }

    // §4.2 documented length ranges: the inclusive endpoints are valid, the just-outside lengths are not.
    [Test]
    public void IsValidPeptideLength_DocumentedRangeEndpoints()
    {
        IsValidPeptideLength(8, MhcClass.ClassI).Should().BeTrue();
        IsValidPeptideLength(11, MhcClass.ClassI).Should().BeTrue();
        IsValidPeptideLength(7, MhcClass.ClassI).Should().BeFalse();
        IsValidPeptideLength(12, MhcClass.ClassI).Should().BeFalse();

        IsValidPeptideLength(13, MhcClass.ClassII).Should().BeTrue();
        IsValidPeptideLength(25, MhcClass.ClassII).Should().BeTrue();
        IsValidPeptideLength(12, MhcClass.ClassII).Should().BeFalse();
        IsValidPeptideLength(26, MhcClass.ClassII).Should().BeFalse();
    }

    #endregion

    #region ONCO-MHC-001 — BE: IC50 = 0 / negative / sub-normal (documented reject, no log(0) leak)

    // §3.3 / INV-01: IC50 = 0 is rejected (a concentration must be strictly positive).
    [Test]
    public void ClassifyBindingAffinity_Ic50Zero_Throws()
    {
        Action act = () => ClassifyBindingAffinity(0.0);
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    // Negative "concentrations" (including −0.0 and the most negative double) are rejected, not silently Strong.
    [Test]
    public void ClassifyBindingAffinity_NegativeIc50_Throws()
    {
        foreach (double bad in new[] { -1e-300, -1.0, -50.0, -double.MaxValue, -0.0 })
        {
            Action act = () => ClassifyBindingAffinity(bad);
            act.Should().Throw<ArgumentOutOfRangeException>($"IC50 {bad} is not a positive concentration");
        }
    }

    // A perfect/extreme affinity is a *positive* but vanishingly small IC50 (not literally 0). It must classify
    // Strong with NO log(0)=−∞ leak — the unit performs NO log transform, so even the smallest sub-normal is
    // a finite, in-contract Strong binder.
    [Test]
    public void ClassifyBindingAffinity_ExtremeTinyPositiveIc50_IsStrongNoInfLeak()
    {
        foreach (double tiny in new[] { double.Epsilon, 1e-300, 1e-12, 0.001 })
        {
            BindingStrength s = ClassifyBindingAffinity(tiny);
            AssertWellFormed(s);
            s.Should().Be(BindingStrength.Strong, $"a near-zero IC50 ({tiny}) is the strongest possible binder");
        }
    }

    // The end-to-end wrapper validates affinity (after a valid length) — IC50 = 0 still throws, never NonBinder.
    [Test]
    public void ClassifyMhcBinding_ValidLengthIc50Zero_Throws()
    {
        Action act = () => ClassifyMhcBinding(peptideLength: 9, ic50Nm: 0.0, MhcClass.ClassI);
        act.Should().Throw<ArgumentOutOfRangeException>();
    }

    #endregion

    #region ONCO-MHC-001 — BE: IC50 = ∞ / NaN (no Inf leaking into the threshold ⇒ documented reject)

    // §3.3 / INV-01: +∞, −∞ and NaN are non-finite and must throw — an Inf must NEVER pass the `< 500`
    // comparison and silently become a NonBinder, and NaN must never make every comparison false.
    [Test]
    public void ClassifyBindingAffinity_InfinityAndNaN_Throw()
    {
        foreach (double bad in new[] { double.PositiveInfinity, double.NegativeInfinity, double.NaN })
        {
            Action act = () => ClassifyBindingAffinity(bad);
            act.Should().Throw<ArgumentOutOfRangeException>($"IC50 {bad} is not finite");
        }
    }

    [Test]
    public void ClassifyMhcBinding_ValidLengthInfiniteIc50_Throws()
    {
        foreach (double bad in new[] { double.PositiveInfinity, double.NaN })
        {
            Action act = () => ClassifyMhcBinding(peptideLength: 10, ic50Nm: bad, MhcClass.ClassI);
            act.Should().Throw<ArgumentOutOfRangeException>();
        }
    }

    // §3.3 / INV-02: %Rank guards its own boundary — ∞ / NaN / out-of-[0,100] throw, so no Inf reaches the
    // percentile cutoff comparison.
    [Test]
    public void ClassifyBindingRank_NonFiniteOrOutOfRange_Throw()
    {
        foreach (double bad in new[]
                 {
                     double.PositiveInfinity, double.NegativeInfinity, double.NaN,
                     -0.0001, -1.0, 100.0001, 1e9
                 })
        {
            Action actI = () => ClassifyBindingRank(bad, MhcClass.ClassI);
            Action actII = () => ClassifyBindingRank(bad, MhcClass.ClassII);
            actI.Should().Throw<ArgumentOutOfRangeException>($"%Rank {bad} is not a percentile");
            actII.Should().Throw<ArgumentOutOfRangeException>($"%Rank {bad} is not a percentile");
        }
    }

    // %Rank boundaries [0, 100] are themselves valid (0 = strongest, 100 = weakest), no false reject.
    [Test]
    public void ClassifyBindingRank_ZeroAndHundred_AreValid()
    {
        ClassifyBindingRank(0.0, MhcClass.ClassI).Should().Be(BindingStrength.Strong);
        ClassifyBindingRank(100.0, MhcClass.ClassI).Should().Be(BindingStrength.NonBinder);
        ClassifyBindingRank(100.0, MhcClass.ClassII).Should().Be(BindingStrength.NonBinder);
    }

    #endregion

    #region ONCO-MHC-001 — BE: peptide too short / too long (length gate, no matrix index overflow)

    // §3.3: IsValidPeptideLength is TOTAL — even the integer extremes return false, never overflow / throw.
    [Test]
    public void IsValidPeptideLength_ExtremeLengths_ReturnFalseNoThrow()
    {
        foreach (int len in new[] { int.MinValue, -1, 0, 1, 7, 12, 100, int.MaxValue })
        {
            // outside class I 8–11 for all these except none (7,12 just outside; 1,0,neg,huge far outside)
            Action actI = () => IsValidPeptideLength(len, MhcClass.ClassI);
            actI.Should().NotThrow();
            IsValidPeptideLength(len, MhcClass.ClassI).Should().BeFalse();
        }
    }

    // §3.3: a too-short / too-long peptide is gated to NonBinder BEFORE the affinity is even looked at — so an
    // out-of-range length never indexes a scoring matrix and a perfectly Strong IC50 is still reported NonBinder.
    [Test]
    public void ClassifyMhcBinding_OutOfRangeLengthWithStrongIc50_IsNonBinder()
    {
        foreach (int len in new[] { int.MinValue, -1, 0, 7, 12, 25, int.MaxValue })
        {
            // 25 is valid for class II but NOT class I; all the rest are invalid for class I.
            ClassifyMhcBinding(len, ic50Nm: 1.0, MhcClass.ClassI)
                .Should().Be(BindingStrength.NonBinder, $"class I rejects length {len}");
        }
    }

    // The length gate runs first: an invalid length must NOT throw even with a degenerate (would-throw) IC50,
    // because the method short-circuits to NonBinder before validating affinity.
    [Test]
    public void ClassifyMhcBinding_InvalidLengthWithBadIc50_DoesNotThrow_NonBinder()
    {
        foreach (double ic50 in new[] { 0.0, -1.0, double.PositiveInfinity, double.NaN })
        {
            BindingStrength s = ClassifyMhcBinding(peptideLength: 0, ic50Nm: ic50, MhcClass.ClassI);
            s.Should().Be(BindingStrength.NonBinder, "invalid length short-circuits before affinity validation");
        }
    }

    // Random length fuzz: for any int length and any class, ClassifyMhcBinding with a valid IC50 returns a
    // well-formed result and is NonBinder exactly when the length is out of the documented range.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifyMhcBinding_RandomLengths_GateMatchesRange()
    {
        var rng = new Random(110_0001);
        for (int i = 0; i < 20_000; i++)
        {
            int len = rng.Next(-50, 60);
            MhcClass cls = rng.Next(2) == 0 ? MhcClass.ClassI : MhcClass.ClassII;
            double ic50 = rng.NextDouble() * 1000.0 + 0.0001; // valid (> 0, finite)

            BindingStrength s = ClassifyMhcBinding(len, ic50, cls);
            AssertWellFormed(s);

            bool valid = IsValidPeptideLength(len, cls);
            if (!valid)
            {
                s.Should().Be(BindingStrength.NonBinder, $"len {len} cls {cls} is out of range");
            }
            else
            {
                s.Should().Be(ExpectedFromIc50(ic50), $"valid len {len} cls {cls} ⇒ affinity rules");
            }
        }
    }

    #endregion

    #region ONCO-MHC-001 — BE: unknown allele / undefined MHC class (no KeyNotFound / NullReference)

    // §5.2: the unit models NO learned per-allele set; the only allele axis is the MhcClass enum. An
    // out-of-domain / undefined enum value (the "unknown allele" analogue) must resolve deterministically with
    // NO KeyNotFoundException / NullReferenceException. The implementation's `== ClassI ? ... : ...` makes any
    // non-ClassI value behave as the class II branch — a defined, non-throwing fallback.
    [Test]
    public void ClassifyBindingRank_UndefinedMhcClass_FallsToClassIIBranchNoThrow()
    {
        var unknown = (MhcClass)999;

        BindingStrength s = default;
        Action act = () => s = ClassifyBindingRank(3.0, unknown);
        act.Should().NotThrow();
        AssertWellFormed(s);
        // class II cutoffs (2/10): 3.0 ⇒ Weak. Confirms the deterministic class-II fallback, not a crash.
        s.Should().Be(BindingStrength.Weak);
    }

    [Test]
    public void IsValidPeptideLength_UndefinedMhcClass_UsesClassIIRangeNoThrow()
    {
        var unknown = (MhcClass)(-7);
        Action act = () => IsValidPeptideLength(20, unknown);
        act.Should().NotThrow();
        // class II range 13–25 ⇒ 20 valid.
        IsValidPeptideLength(20, unknown).Should().BeTrue();
        IsValidPeptideLength(9, unknown).Should().BeFalse();
    }

    [Test]
    public void ClassifyMhcBinding_UndefinedMhcClass_NoThrowWellFormed()
    {
        var unknown = (MhcClass)int.MaxValue;
        BindingStrength s = default;
        Action act = () => s = ClassifyMhcBinding(20, 30.0, unknown);
        act.Should().NotThrow();
        AssertWellFormed(s);
        // class II range admits 20; IC50 30 < 50 ⇒ Strong.
        s.Should().Be(BindingStrength.Strong);
    }

    #endregion

    #region ONCO-MHC-001 — Invariants: monotonicity & well-formedness over broad random fuzz

    // INV-03 monotonicity: over many random valid IC50 pairs, the smaller IC50 yields an equal-or-stronger
    // category. Encodes Strong < Weak < NonBinder as 0 < 1 < 2 and asserts the ordering is preserved.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifyBindingAffinity_RandomValidIc50_MonotoneAndWellFormed()
    {
        var rng = new Random(110_0002);
        for (int i = 0; i < 30_000; i++)
        {
            // valid IC50s spanning sub-nM to far above the weak cutoff
            double a = Math.Pow(10, rng.NextDouble() * 6 - 2); // 1e-2 .. 1e4 nM
            double b = Math.Pow(10, rng.NextDouble() * 6 - 2);

            BindingStrength sa = ClassifyBindingAffinity(a);
            BindingStrength sb = ClassifyBindingAffinity(b);
            AssertWellFormed(sa);
            AssertWellFormed(sb);

            // cross-check against the independent oracle
            sa.Should().Be(ExpectedFromIc50(a));
            sb.Should().Be(ExpectedFromIc50(b));

            if (a <= b)
            {
                ((int)sa).Should().BeLessThanOrEqualTo((int)sb,
                    $"smaller IC50 {a} must be equal-or-stronger than {b}");
            }
        }
    }

    // INV-03 monotonicity for %Rank, per class: smaller rank ⇒ equal-or-stronger category.
    [Test]
    [CancelAfter(20_000)]
    public void ClassifyBindingRank_RandomValidRank_MonotonePerClass()
    {
        var rng = new Random(110_0003);
        for (int i = 0; i < 30_000; i++)
        {
            double a = rng.NextDouble() * 100.0;
            double b = rng.NextDouble() * 100.0;
            MhcClass cls = rng.Next(2) == 0 ? MhcClass.ClassI : MhcClass.ClassII;

            BindingStrength sa = ClassifyBindingRank(a, cls);
            BindingStrength sb = ClassifyBindingRank(b, cls);
            AssertWellFormed(sa);
            AssertWellFormed(sb);

            if (a <= b)
            {
                ((int)sa).Should().BeLessThanOrEqualTo((int)sb,
                    $"smaller %Rank {a} must be equal-or-stronger than {b} ({cls})");
            }
        }
    }

    #endregion
}
