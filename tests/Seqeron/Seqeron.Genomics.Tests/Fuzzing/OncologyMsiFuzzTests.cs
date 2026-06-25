using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology microsatellite-instability area — ONCO-MSI-001.
/// The units under test are the deterministic scoring-and-classification entry
/// points implemented in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs:
/// <see cref="OncologyAnalyzer.CalculateMSIScore(int, int)"/> (fraction
/// unstable/total), <see cref="OncologyAnalyzer.ClassifyMSIStatus(double)"/>
/// (continuous MSIsensor2 ≥20% cutoff), <see cref="OncologyAnalyzer.ClassifyBethesdaPanel(int, int)"/>
/// (categorical NCI/Bethesda count rule), and the end-to-end
/// <see cref="OncologyAnalyzer.DetectMSI(System.Collections.Generic.IEnumerable{bool})"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / extreme inputs to a unit and asserts
/// that the code NEVER fails in an undisciplined way: no hang, no nonsense
/// output, and no *unhandled* runtime fault. Every input must resolve to EITHER
/// a well-defined, theory-correct value OR a *documented, intentional* outcome
/// (here, an <see cref="ArgumentOutOfRangeException"/> for ≤0 loci, an
/// out-of-range count, or a non-finite / out-of-[0,1] score; an
/// <see cref="ArgumentNullException"/> for a null flag sequence).
/// For MSI the headline hazards (Microsatellite_Instability_Detection.md §3.3,
/// §6.1) are:
///   • ZERO loci → the denominator of u/n is 0; the contract is a documented
///     ArgumentOutOfRangeException ("MSI score undefined with no evaluable loci"),
///     NEVER a DivideByZero, a +Infinity / NaN score leaking out, and NEVER a
///     silent MSS/MSI-H call on an undefined sample (§3.3, §6.1, edge "0 valid loci");
///   • ALL-STABLE (0 unstable) → score = 0.0, classified MSS, no throw, no NaN (§6.1);
///   • ALL-UNSTABLE (u = n) → score = 1.0, classified MSI-High (1.0 ≥ 0.20) (INV-01/03);
///   • SINGLE read / minimal evidence (n = 1) → still a valid sample of exactly
///     one locus: score is 0.0 (stable) or 1.0 (unstable), classified accordingly,
///     never a crash and never DivideByZero (§4.1, INV-02);
///   • out-of-range counts (unstable<0, unstable>total) and out-of-[0,1] / NaN /
///     ±Infinity scores → documented ArgumentOutOfRangeException, never a score
///     outside [0,1] (INV-01).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-MSI-001 — Microsatellite instability detection (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 93.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, 1, n=u, empty/null.
///     Targets (checklist row 93): "zero loci, all-stable, all-unstable,
///     single read".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// Microsatellite_Instability_Detection.md (docs/algorithms/Oncology/...):
///   • MSI score = u / n (unstable / valid loci), a fraction in [0,1]   (§2.2, INV-01/02)
///   • continuous status = MSI-High iff score ≥ 0.20 (inclusive), else MSS
///                                                                       (§2.2, INV-03)
///   • Bethesda categorical: ≥2 unstable markers → MSI-H, exactly 1 → MSI-L,
///       0 → MSS                                                         (§2.2, INV-04)
///   • DetectMSI: count → score → continuous status, one entry per valid locus (§4.1)
///   • 0 valid loci / unstable>valid ⇒ ArgumentOutOfRangeException        (§3.3, §6.1)
///   • non-finite or out-of-[0,1] score ⇒ ArgumentOutOfRangeException     (§3.3)
///   • empty flags ⇒ ArgumentOutOfRangeException; null ⇒ ArgumentNullException (§3.3, §6.1)
///   • score = 0.20 ⇒ MSI-H (inclusive); score = 0.0 ⇒ MSS               (§6.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyMsiFuzzTests
{
    private const double Cutoff = MsiHighScoreThreshold; // 0.20, MSIsensor2 inclusive

    // ── Well-formed-score assertion helper ───────────────────────────────────
    // Pins the documented numeric contract on EVERY returned MSI score: it must
    // be FINITE (no DivideByZero +Infinity, no NaN leaking through) and inside
    // the closed unit interval [0,1] (INV-01). This is what stops a fuzz test
    // from rubber-stamping an Infinity / NaN / out-of-range score green.
    private static void AssertWellFormedScore(double score)
    {
        double.IsNaN(score).Should().BeFalse("MSI score must never be NaN");
        double.IsInfinity(score).Should().BeFalse("MSI score = u/n must be finite (n > 0)");
        score.Should().BeInRange(0.0, 1.0, "0 ≤ MSI score ≤ 1 (INV-01)");
    }

    // Builds a flag sequence with the requested number of unstable/stable loci.
    private static bool[] Flags(int unstable, int stable) =>
        Enumerable.Repeat(true, unstable).Concat(Enumerable.Repeat(false, stable)).ToArray();

    #region ONCO-MSI-001 — Positive sanity (hand-computed fraction + correct call)

    // ── POSITIVE sanity: documented unstable fraction → exact score + class ───
    [Test]
    public void DetectMSI_KnownUnstableFraction_HandComputedScoreAndMsiHigh()
    {
        // 8 unstable of 20 valid loci = 0.40 ≥ 0.20 ⇒ MSI-High (MSIsensor2).
        MsiResult r = DetectMSI(Flags(unstable: 8, stable: 12));

        r.UnstableLoci.Should().Be(8);
        r.TotalLoci.Should().Be(20);
        AssertWellFormedScore(r.Score);
        r.Score.Should().BeApproximately(0.40, 1e-12, "8 / 20 = 0.40 (§2.2)");
        r.Status.Should().Be(MsiStatus.MSI_High, "0.40 ≥ 0.20 ⇒ MSI-High (INV-03)");
    }

    [Test]
    public void DetectMSI_DocWorkedExample_SixOfTwentyIsMsiHigh()
    {
        // Docs §7.1 worked example: 6 unstable of 20 → 0.30 ≥ 0.20 → MSI-High.
        MsiResult r = DetectMSI(Flags(unstable: 6, stable: 14));

        AssertWellFormedScore(r.Score);
        r.Score.Should().BeApproximately(0.30, 1e-12, "6 / 20 = 0.30 (§7.1)");
        r.Status.Should().Be(MsiStatus.MSI_High);
    }

    [Test]
    public void DetectMSI_SubThresholdFraction_ClassifiesMss()
    {
        // 3 unstable of 20 = 0.15 < 0.20 ⇒ MSS (no MSI-L band on the score, §5.2).
        MsiResult r = DetectMSI(Flags(unstable: 3, stable: 17));

        AssertWellFormedScore(r.Score);
        r.Score.Should().BeApproximately(0.15, 1e-12, "3 / 20 = 0.15");
        r.Status.Should().Be(MsiStatus.MSS, "0.15 < 0.20 ⇒ MSS (INV-03)");
    }

    [Test]
    public void ClassifyMSIStatus_ExactlyAtCutoff_IsMsiHigh()
    {
        // §6.1: score = 0.20 ⇒ MSI-H because the cutoff is inclusive (≥).
        ClassifyMSIStatus(Cutoff).Should().Be(MsiStatus.MSI_High, "cutoff ≥ 0.20 inclusive");
    }

    [Test]
    public void ClassifyMSIStatus_JustBelowCutoff_IsMss()
    {
        ClassifyMSIStatus(Math.BitDecrement(Cutoff)).Should().Be(MsiStatus.MSS,
            "the largest double < 0.20 must be MSS — the boundary is strict from below");
    }

    [Test]
    public void ClassifyBethesdaPanel_TwoOfFive_IsMsiHigh_OneIsLow_NoneIsMss()
    {
        // §2.2 / INV-04: ≥2→MSI-H, exactly 1→MSI-L, 0→MSS over the 5-marker panel.
        ClassifyBethesdaPanel(2, 5).Should().Be(MsiStatus.MSI_High);
        ClassifyBethesdaPanel(1, 5).Should().Be(MsiStatus.MSI_Low);
        ClassifyBethesdaPanel(0, 5).Should().Be(MsiStatus.MSS);
    }

    #endregion

    #region ONCO-MSI-001 — BE: zero loci (no DivideByZero, documented throw)

    [Test]
    public void CalculateMSIScore_ZeroTotalLoci_ThrowsNoDivideByZero()
    {
        // §3.3 / §6.1 edge "0 valid loci": u/n is undefined with n = 0. The
        // contract is a documented ArgumentOutOfRangeException — NEVER a
        // DivideByZeroException, a +Infinity / NaN score, nor a silent call.
        Action act = () => CalculateMSIScore(unstableLoci: 0, totalLoci: 0);

        act.Should().Throw<ArgumentOutOfRangeException>()
            .Which.ParamName.Should().Be("totalLoci");
    }

    [Test]
    public void DetectMSI_EmptyFlagSequence_ThrowsNoDivideByZero()
    {
        // §3.3 / §6.1 edge "empty flags": an empty locus set has no valid loci.
        Action act = () => DetectMSI(Array.Empty<bool>());

        act.Should().Throw<ArgumentOutOfRangeException>(
            "an empty flag set has 0 valid loci ⇒ undefined score (§6.1)");
    }

    [Test]
    public void DetectMSI_NullFlagSequence_ThrowsArgumentNull()
    {
        Action act = () => DetectMSI(null!);

        act.Should().Throw<ArgumentNullException>("null is guarded (§3.3)");
    }

    [Test]
    public void ClassifyBethesdaPanel_ZeroTotalMarkers_Throws()
    {
        Action act = () => ClassifyBethesdaPanel(unstableMarkers: 0, totalMarkers: 0);

        act.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("totalMarkers");
    }

    #endregion

    #region ONCO-MSI-001 — BE: all-stable (0% unstable → MSS, not flagged)

    [Test]
    public void DetectMSI_AllStable_ScoreZeroAndMss()
    {
        // §6.1 edge "score = 0.0": every locus stable ⇒ u = 0 ⇒ score 0.0 ⇒ MSS.
        // The hazard guarded: a 0 numerator must be a valid finite 0, not a NaN
        // or a tripped guard, and must NOT be mis-flagged MSI-H.
        MsiResult r = DetectMSI(Flags(unstable: 0, stable: 30));

        r.UnstableLoci.Should().Be(0);
        r.TotalLoci.Should().Be(30);
        AssertWellFormedScore(r.Score);
        r.Score.Should().Be(0.0, "0 / 30 = 0 (§6.1)");
        r.Status.Should().Be(MsiStatus.MSS, "0 < 0.20 ⇒ MSS, never flagged MSI-H");
    }

    [Test]
    public void CalculateMSIScore_ZeroUnstable_AnyValidLoci_IsExactlyZero()
    {
        var rng = new Random(93_001);
        for (int i = 0; i < 200; i++)
        {
            int total = rng.Next(1, 5_000);
            double score = CalculateMSIScore(unstableLoci: 0, totalLoci: total);

            AssertWellFormedScore(score);
            score.Should().Be(0.0, "0 unstable ⇒ score 0 regardless of n (INV-02)");
            ClassifyMSIStatus(score).Should().Be(MsiStatus.MSS);
        }
    }

    #endregion

    #region ONCO-MSI-001 — BE: all-unstable (100% unstable → MSI-High)

    [Test]
    public void DetectMSI_AllUnstable_ScoreOneAndMsiHigh()
    {
        // INV-01/03: every locus unstable ⇒ u = n ⇒ score = 1.0 ⇒ MSI-High.
        MsiResult r = DetectMSI(Flags(unstable: 25, stable: 0));

        r.UnstableLoci.Should().Be(25);
        r.TotalLoci.Should().Be(25);
        AssertWellFormedScore(r.Score);
        r.Score.Should().Be(1.0, "25 / 25 = 1.0 (INV-01)");
        r.Status.Should().Be(MsiStatus.MSI_High, "1.0 ≥ 0.20 ⇒ MSI-High (INV-03)");
    }

    [Test]
    public void CalculateMSIScore_FullUnstable_AnyValidLoci_IsExactlyOne()
    {
        var rng = new Random(93_002);
        for (int i = 0; i < 200; i++)
        {
            int total = rng.Next(1, 5_000);
            double score = CalculateMSIScore(unstableLoci: total, totalLoci: total);

            AssertWellFormedScore(score);
            score.Should().Be(1.0, "u = n ⇒ score 1.0 (INV-01)");
            ClassifyMSIStatus(score).Should().Be(MsiStatus.MSI_High);
        }
    }

    #endregion

    #region ONCO-MSI-001 — BE: single read / minimal evidence (n = 1)

    [Test]
    public void DetectMSI_SingleStableLocus_ScoreZeroMss_NoDivideByZero()
    {
        // §4.1 / INV-02: a sample of exactly one valid locus is well-defined.
        // One stable locus ⇒ score 0.0 ⇒ MSS, never a crash, never DivideByZero.
        MsiResult r = DetectMSI(new[] { false });

        r.TotalLoci.Should().Be(1);
        r.UnstableLoci.Should().Be(0);
        AssertWellFormedScore(r.Score);
        r.Score.Should().Be(0.0);
        r.Status.Should().Be(MsiStatus.MSS);
    }

    [Test]
    public void DetectMSI_SingleUnstableLocus_ScoreOneMsiHigh()
    {
        // One unstable locus ⇒ 1/1 = 1.0 ⇒ MSI-High (the minimal MSI-H sample).
        MsiResult r = DetectMSI(new[] { true });

        r.TotalLoci.Should().Be(1);
        r.UnstableLoci.Should().Be(1);
        AssertWellFormedScore(r.Score);
        r.Score.Should().Be(1.0);
        r.Status.Should().Be(MsiStatus.MSI_High, "1/1 = 1.0 ≥ 0.20 ⇒ MSI-High");
    }

    [Test]
    public void ClassifyBethesdaPanel_SingleMarkerPanel_OneUnstableIsMsiHigh()
    {
        // Minimal Bethesda panel of one marker: 1 unstable ≥ 2-rule? No — the
        // rule is ≥2 → MSI-H, so on a 1-marker panel 1 unstable is MSI-L
        // (exactly one), and 0 is MSS. No DivideByZero / no crash on n = 1.
        ClassifyBethesdaPanel(unstableMarkers: 1, totalMarkers: 1).Should().Be(MsiStatus.MSI_Low);
        ClassifyBethesdaPanel(unstableMarkers: 0, totalMarkers: 1).Should().Be(MsiStatus.MSS);
    }

    #endregion

    #region ONCO-MSI-001 — BE: out-of-range counts and scores (documented throws)

    [Test]
    public void CalculateMSIScore_NegativeUnstable_Throws()
    {
        Action act = () => CalculateMSIScore(unstableLoci: -1, totalLoci: 10);

        act.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("unstableLoci");
    }

    [Test]
    public void CalculateMSIScore_UnstableExceedsTotal_Throws()
    {
        // 0 ≤ u ≤ n (INV-01 precondition); u > n is impossible by definition.
        Action act = () => CalculateMSIScore(unstableLoci: 11, totalLoci: 10);

        act.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("unstableLoci");
    }

    [Test]
    public void CalculateMSIScore_NegativeTotal_Throws()
    {
        Action act = () => CalculateMSIScore(unstableLoci: 0, totalLoci: -5);

        act.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("totalLoci");
    }

    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NegativeInfinity)]
    [TestCase(-0.0001)]
    [TestCase(1.0001)]
    public void ClassifyMSIStatus_NonFiniteOrOutOfRange_Throws(double badScore)
    {
        Action act = () => ClassifyMSIStatus(badScore);

        act.Should().Throw<ArgumentOutOfRangeException>(
            "a score must be a finite value in [0,1] (§3.3)");
    }

    #endregion

    #region ONCO-MSI-001 — BE: randomized fraction invariants (no leak, exact u/n)

    [Test]
    [CancelAfter(20_000)]
    public void CalculateMSIScore_RandomCounts_ScoreInUnitIntervalAndMatchesRatio()
    {
        // Sweep arbitrary (u, n) with 0 ≤ u ≤ n: the score must stay finite and
        // in [0,1] (INV-01) and equal u/n exactly (INV-02); the continuous call
        // must agree with the documented inclusive 0.20 cutoff (INV-03). This is
        // the anti-rubber-stamp net: it would fail on any leaked NaN/Infinity, a
        // score outside [0,1], or a flipped cutoff direction.
        var rng = new Random(93_777);
        for (int i = 0; i < 5_000; i++)
        {
            int total = rng.Next(1, 10_000);
            int unstable = rng.Next(0, total + 1);

            double score = CalculateMSIScore(unstable, total);

            AssertWellFormedScore(score);
            score.Should().BeApproximately((double)unstable / total, 1e-12, "score = u/n (INV-02)");

            MsiStatus expected = score >= Cutoff ? MsiStatus.MSI_High : MsiStatus.MSS;
            ClassifyMSIStatus(score).Should().Be(expected, "continuous call follows the ≥0.20 cutoff (INV-03)");
        }
    }

    [Test]
    [CancelAfter(20_000)]
    public void DetectMSI_RandomFlagSets_CountsScoreAndStatusConsistent()
    {
        // End-to-end consistency over random per-locus flag sets: the reported
        // counts, the score (= u/n), and the status must all agree internally,
        // and no degenerate-but-non-empty set may crash or leak a bad score.
        var rng = new Random(93_888);
        for (int i = 0; i < 2_000; i++)
        {
            int total = rng.Next(1, 500);
            var flags = new bool[total];
            int unstable = 0;
            for (int j = 0; j < total; j++)
            {
                bool u = rng.NextDouble() < 0.5;
                flags[j] = u;
                if (u) unstable++;
            }

            MsiResult r = DetectMSI(flags);

            r.TotalLoci.Should().Be(total);
            r.UnstableLoci.Should().Be(unstable);
            AssertWellFormedScore(r.Score);
            r.Score.Should().BeApproximately((double)unstable / total, 1e-12);
            r.Status.Should().Be(r.Score >= Cutoff ? MsiStatus.MSI_High : MsiStatus.MSS);
        }
    }

    #endregion
}
