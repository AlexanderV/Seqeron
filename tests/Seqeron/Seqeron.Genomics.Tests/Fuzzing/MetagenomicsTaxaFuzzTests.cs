using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Metagenomics significant-taxa unit — the two-group
/// differential-abundance test via the Mann–Whitney U (Wilcoxon rank-sum)
/// statistic: <see cref="MetagenomicsAnalyzer.MannWhitneyU"/> and the per-taxon
/// driver <see cref="MetagenomicsAnalyzer.FindSignificantTaxa"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// no NaN/Infinity leaking into a z-score or p-value, and no *unhandled* runtime
/// exception (IndexOutOfRangeException, NullReferenceException,
/// DivideByZeroException, OverflowException, …). Every input must produce EITHER
/// a well-defined, theory-correct result, OR a *documented, intentional*
/// validation exception (ArgumentNullException / ArgumentException). A raw
/// runtime exception, a hang, a p-value outside [0, 1], or a *fabricated*
/// significance on degenerate input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-TAXA-001 — significant taxa detection (Metagenomics)
/// Checklist: docs/checklists/03_FUZZING.md, row 197.
/// Fuzz strategy for THIS unit: BE = Boundary Exploitation (0, -1, MaxInt, empty)
///   — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Fuzz targets (checklist row 197): "identical samples, single taxon, empty".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// MannWhitneyU pools n1 + n2 = n observations, assigns midranks (ties share the
/// average of their positions), and computes
///   U1 = R1 − n1(n1+1)/2,   U2 = n1·n2 − U1,
///   m_U = n1·n2/2,
///   σ_U = sqrt( n1·n2·(n+1)/12 − n1·n2·Σ(t³−t) / (12·n·(n−1)) ),
///   z   = (|max(U1,U2) − m_U| − cc) / σ_U,   cc ∈ {0, 0.5},
///   p   = clamp(2·(1 − Φ(z)), 0, 1);   σ_U ≤ 0 ⇒ z = 0, p = 1.
/// FindSignificantTaxa runs that test per taxon over the two label-defined groups
/// (absent taxon ⇒ abundance 0), flags PValue &lt; pThreshold, and returns the taxa
/// ascending by p-value.
///   — docs/algorithms/Metagenomics/Significant_Taxa_Detection.md §2.2, §4;
///     src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs
///     (MannWhitneyU lines 1493–1567, FindSignificantTaxa 1581–1629).
///
/// Documented invariants this fixture pins (Significant_Taxa_Detection.md §2.4):
///   • INV-01: U1 + U2 = n1·n2.
///   • INV-02: 0 ≤ U ≤ n1·n2 (both U ≥ 0).
///   • INV-03: p-value ∈ [0, 1].
///   • INV-05: all-tied groups → p = 1, not significant (σ_U → 0).
///   • INV-06: swapping group1/group2 leaves p (and the larger U) unchanged.
///
/// Boundary / malformed-input handling fixed by the doc (§3.3, §6.1) and pinned
/// here so the contract can never silently drift:
///   • IDENTICAL SAMPLES (BE): two groups with the SAME multiset of abundances
///     (the pooled values are all tied across groups) ⇒ σ_U → 0 ⇒ z = 0, p = 1,
///     NOT significant. Significance is NEVER fabricated when the groups match.
///     — INV-05, §2.4, §6.1 ("All observations identical → z=0, p=1").
///   • SINGLE TAXON (BE): a profile set with exactly one taxon ⇒ exactly one
///     SignificantTaxon result, no crash. — §3.2, §4.1 step 6.
///   • EMPTY (BE):
///       – empty profiles list ⇒ empty result, no crash. — §3.3, §6.1.
///       – empty group to MannWhitneyU ⇒ ArgumentException (test undefined at n=0).
///         — §3.3, §6.1 ("Empty group → ArgumentException").
///       – a taxon absent from a profile counts as abundance 0 (NOT missing).
///         — ASM-03, §3.3, §6.1.
///   • Argument validation: null group/profiles/groups ⇒ ArgumentNullException;
///     length mismatch, a label ∉ {1,2}, or a missing group ⇒ ArgumentException.
///     — §3.3, §6.1, source 1498–1604.
///
/// Positive sanity (worked example, derived INDEPENDENTLY from the U formula and
/// the standard-normal SF, NOT echoed off the implementation):
///   • §7.1 SciPy example: x=[19,22,16,29,24] (n1=5), y=[20,11,17,12] (n2=4)
///     ⇒ U1 = 17, U2 = 3 (U1+U2 = 20 = n1·n2), m_U = 10, σ_U = sqrt(200/12) = 4.0825.
///       with continuity:  z = (17−10−0.5)/4.0825 = 1.5922 → p ≈ 0.1113.
///       without:          z = (17−10)/4.0825     = 1.7146 → p ≈ 0.0864.
/// A genuine separation must therefore yield U1 = 17 / U2 = 3 with a small,
/// strictly-positive, &lt; 1 p-value, so a passing "no crash" result cannot be a
/// degenerate analyzer that returns p = 1.0 (or U = 0) for everything.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Determinism
/// ───────────────────────────────────────────────────────────────────────────
/// All inputs are hand-built or generated from a LOCALLY fixed-seed
/// `new Random(seed)` (never a shared static Rng), so every run is reproducible.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MetagenomicsTaxaFuzzTests
{
    // ── Independently-derived constants for the §7.1 SciPy worked example ──
    //   x=[19,22,16,29,24], y=[20,11,17,12]: U1 = 17, U2 = 3, m_U = 10,
    //   σ_U = sqrt(200/12) = 4.08248..., p_cc ≈ 0.1113, p_nocc ≈ 0.0864.
    // These are NOT echoed off the implementation; they are the doc/SciPy values.
    private static readonly double[] ScipyX = { 19, 22, 16, 29, 24 };
    private static readonly double[] ScipyY = { 20, 11, 17, 12 };
    private const double ScipyU1 = 17.0;
    private const double ScipyU2 = 3.0;
    private const double ScipyZContinuity = 1.5921683328090657;   // (17−10−0.5)/sqrt(200/12)
    private const double ScipyZNoContinuity = 1.7146428199482247; // (17−10)/sqrt(200/12)
    private const double ScipyPContinuity = 0.1113468865331404;   // 2·(1−Φ(1.5922))
    private const double ScipyPNoContinuity = 0.08641073297370006; // 2·(1−Φ(1.7146))

    #region META-TAXA-001 — Mann–Whitney significant-taxa detection

    // ════════════════════════════════════════════════════════════════════════
    //  Positive sanity — the §7.1 worked example must be reproduced EXACTLY.
    //  Guards against a degenerate analyzer (constant p = 1.0, or U = 0) that
    //  would otherwise pass every boundary test below.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void MannWhitneyU_ScipyWorkedExample_ReturnsU17And3()
    {
        // §7.1: pooled midranks ⇒ U1 = 17, U2 = 3 (independently of the code).
        var r = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY);

        r.U1.Should().BeApproximately(ScipyU1, 1e-12, "§7.1: U1 = R1 − n1(n1+1)/2 = 17");
        r.U2.Should().BeApproximately(ScipyU2, 1e-12, "§7.1: U2 = n1·n2 − U1 = 20 − 17 = 3");
    }

    [Test]
    public void MannWhitneyU_ScipyWorkedExample_WithContinuity_MatchesZAndP()
    {
        // §7.1 with the default SciPy continuity correction (cc = 0.5).
        var r = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY); // cc = true by default

        r.Z.Should().BeApproximately(ScipyZContinuity, 1e-6,
            "§7.1: z = (17−10−0.5)/sqrt(200/12) = 1.5922");
        r.PValue.Should().BeApproximately(ScipyPContinuity, 1e-6,
            "§7.1: p = 2·(1−Φ(1.5922)) ≈ 0.1113 (SciPy reference)");
    }

    [Test]
    public void MannWhitneyU_ScipyWorkedExample_NoContinuity_MatchesZAndP()
    {
        // §7.1 without the continuity correction (cc = 0).
        var r = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY, useContinuityCorrection: false);

        r.Z.Should().BeApproximately(ScipyZNoContinuity, 1e-6,
            "§7.1: z = (17−10)/sqrt(200/12) = 1.7146");
        r.PValue.Should().BeApproximately(ScipyPNoContinuity, 1e-6,
            "§7.1: p = 2·(1−Φ(1.7146)) ≈ 0.0864 (SciPy reference)");
    }

    // ───────────────────────────────────────────────────────────────────────
    // INV-01 / INV-02: U1 + U2 = n1·n2 with both U ≥ 0 ⇒ 0 ≤ U ≤ n1·n2.
    // Pinned over a randomized sweep so the rank arithmetic can never drift.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void MannWhitneyU_RandomGroups_USumEqualsProductWithBothNonNegative(CancellationToken token)
    {
        var rng = new Random(197_001);

        for (int trial = 0; trial < 4000 && !token.IsCancellationRequested; trial++)
        {
            int n1 = rng.Next(1, 12);
            int n2 = rng.Next(1, 12);
            // Deliberately a small value range so ties (zeros, repeats) occur often.
            var g1 = Enumerable.Range(0, n1).Select(_ => (double)rng.Next(0, 5)).ToList();
            var g2 = Enumerable.Range(0, n2).Select(_ => (double)rng.Next(0, 5)).ToList();

            var r = MetagenomicsAnalyzer.MannWhitneyU(g1, g2);

            double product = (double)n1 * n2;
            r.U1.Should().BeGreaterThanOrEqualTo(0.0, "INV-02: U1 ≥ 0");
            r.U2.Should().BeGreaterThanOrEqualTo(0.0, "INV-02: U2 ≥ 0");
            (r.U1 + r.U2).Should().BeApproximately(product, 1e-9, "INV-01: U1 + U2 = n1·n2");
            Math.Max(r.U1, r.U2).Should().BeLessThanOrEqualTo(product + 1e-9, "INV-02: U ≤ n1·n2");

            double.IsNaN(r.Z).Should().BeFalse("z must never be NaN");
            double.IsInfinity(r.Z).Should().BeFalse("z must never be Infinity");
            r.PValue.Should().BeInRange(0.0, 1.0, "INV-03: p-value ∈ [0,1]");
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // INV-06: swapping group1/group2 leaves the larger U and the p-value
    // unchanged (z uses max(U1,U2) / |U − m_U|, symmetric in the groups).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void MannWhitneyU_SwapGroups_LargerUAndPValueUnchanged()
    {
        var forward = MetagenomicsAnalyzer.MannWhitneyU(ScipyX, ScipyY);
        var swapped = MetagenomicsAnalyzer.MannWhitneyU(ScipyY, ScipyX);

        Math.Max(swapped.U1, swapped.U2).Should().BeApproximately(
            Math.Max(forward.U1, forward.U2), 1e-12, "INV-06: max(U1,U2) is symmetric");
        swapped.PValue.Should().BeApproximately(forward.PValue, 1e-12,
            "INV-06: swapping the groups leaves p unchanged");
        // Sanity: the swap actually exchanges U1 and U2 (not a no-op).
        swapped.U1.Should().BeApproximately(forward.U2, 1e-12);
        swapped.U2.Should().BeApproximately(forward.U1, 1e-12);
    }

    #endregion

    #region META-TAXA-001 — BE boundary: IDENTICAL SAMPLES (all-tied → p = 1)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: IDENTICAL SAMPLES (BE). When the two groups carry the SAME
    // multiset of abundances, U1 = U2 = m_U, so |U − m_U| = 0 ⇒ (with the default
    // continuity correction, max(0, 0 − 0.5)) z = 0. The two-tailed
    // p = 2·(1 − Φ(0)) ≈ 1, and the taxon is NOT significant. Significance must
    // NEVER be fabricated when the samples match.
    //   — INV-05/INV-06, §2.4, §6.1. (p is ≈ 1, not EXACTLY 1, because Φ here is
    //   the A&S 7.1.26 approximation, |ε| ≤ 1.5×10⁻⁷ at z=0 — §5.3; EXACTLY 1.0
    //   requires the genuinely degenerate σ=0 branch, covered separately below.)
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void MannWhitneyU_IdenticalGroups_GivesZeroZAndPNearOne()
    {
        var g1 = new double[] { 3, 7, 7, 2, 9 };
        var g2 = new double[] { 3, 7, 7, 2, 9 }; // same multiset

        var r = MetagenomicsAnalyzer.MannWhitneyU(g1, g2);

        // U1 = U2 = m_U = n1·n2/2 = 12.5 ⇒ |U − m_U| = 0.
        r.U1.Should().BeApproximately(r.U2, 1e-12, "identical multiset ⇒ U1 = U2 = m_U");
        r.Z.Should().Be(0.0, "INV-05/INV-06: |U − m_U| = 0, cc clamps to 0 ⇒ z = 0");
        r.PValue.Should().BeApproximately(1.0, 1e-6, "z = 0 ⇒ p = 2·(1 − Φ(0)) ≈ 1");
    }

    [Test]
    public void MannWhitneyU_AllObservationsConstant_GivesExactlyZeroZAndPOne()
    {
        // Every pooled value equal to ONE value (within and across groups) ⇒
        // genuinely degenerate σ_U ≤ 0 ⇒ the σ=0 branch fires: z = 0, p = 1 EXACTLY.
        var g1 = Enumerable.Repeat(4.0, 6).ToArray();
        var g2 = Enumerable.Repeat(4.0, 5).ToArray();

        var r = MetagenomicsAnalyzer.MannWhitneyU(g1, g2);

        r.Z.Should().Be(0.0, "INV-05: all-tied to one value ⇒ σ=0 branch ⇒ z = 0");
        r.PValue.Should().Be(1.0, "INV-05: σ=0 ⇒ p = 1 exactly (no Φ involved)");
    }

    [Test]
    public void FindSignificantTaxa_IdenticalGroupProfiles_NoTaxonSignificant()
    {
        // Each taxon has the SAME per-sample abundance vector in both groups, so
        // every per-taxon test is all-tied ⇒ p = 1, never flagged significant.
        var profiles = new IReadOnlyDictionary<string, double>[]
        {
            new Dictionary<string, double> { ["A"] = 5, ["B"] = 2 }, // group 1
            new Dictionary<string, double> { ["A"] = 8, ["B"] = 1 }, // group 1
            new Dictionary<string, double> { ["A"] = 5, ["B"] = 2 }, // group 2 (mirror)
            new Dictionary<string, double> { ["A"] = 8, ["B"] = 1 }, // group 2 (mirror)
        };
        var groups = new[] { 1, 1, 2, 2 };

        var results = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, pThreshold: 0.05);

        results.Should().HaveCount(2, "two taxa A, B");
        // Each taxon's group vectors are the same multiset ⇒ U1 = U2 = m_U ⇒ z = 0
        // ⇒ p = 2·(1 − Φ(0)) ≈ 1 (A&S approximation; ≈, not exactly, 1 — §5.3).
        results.Should().OnlyContain(t => Math.Abs(t.PValue - 1.0) < 1e-6,
            "INV-05/INV-06: matched group abundances ⇒ z = 0 per taxon ⇒ p ≈ 1");
        results.Should().OnlyContain(t => !t.Significant,
            "INV-04/INV-05: nothing is significant when the groups are identical");
    }

    // Randomized boundary sweep: identical-multiset groups give U1 = U2 = m_U, so
    // (with the continuity correction) z = 0 always, and p ≈ 1 — regardless of the
    // random shared abundance vector. (Without cc and with distinct pooled values,
    // |U − m_U| = 0 still ⇒ z = 0; the all-one-value case hits the exact σ=0 branch.)
    [Test]
    [CancelAfter(30000)]
    public void MannWhitneyU_RandomIdenticalGroups_AlwaysZeroZAndPNearOne(CancellationToken token)
    {
        var rng = new Random(197_002);

        for (int trial = 0; trial < 3000 && !token.IsCancellationRequested; trial++)
        {
            int n = rng.Next(1, 12);
            var shared = Enumerable.Range(0, n).Select(_ => (double)rng.Next(-3, 8)).ToArray();

            // group2 is a shuffled copy of group1 ⇒ same multiset ⇒ U1 = U2 = m_U.
            var g2 = shared.OrderBy(_ => rng.Next()).ToArray();

            // cc = true here so |U − m_U| = 0 always clamps to z = 0 (cc = false can
            // leave a tiny non-zero z when 0 − 0 is divided by σ — still p ≈ 1, but we
            // pin the strict z = 0 contract under the documented default cc).
            var r = MetagenomicsAnalyzer.MannWhitneyU(shared, g2, useContinuityCorrection: true);

            r.U1.Should().BeApproximately(r.U2, 1e-9, "identical multiset ⇒ U1 = U2 = m_U");
            r.Z.Should().Be(0.0, "INV-05/INV-06: identical multiset ⇒ z = 0");
            r.PValue.Should().BeApproximately(1.0, 1e-6, "INV-05/INV-06: identical multiset ⇒ p ≈ 1");
        }
    }

    #endregion

    #region META-TAXA-001 — BE boundary: SINGLE TAXON

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: SINGLE TAXON (BE). A profile set mentioning exactly one taxon
    // must yield exactly one SignificantTaxon, with a well-formed p-value, no
    // crash. The separated case must also produce the §7.1 statistics, proving
    // the single-taxon path is the real Mann–Whitney test, not a stub.
    //   — §3.2, §4.1 step 6.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindSignificantTaxa_SingleSeparatedTaxon_MatchesScipyExample()
    {
        // Encode the §7.1 example as one taxon "T": group 1 = ScipyX, group 2 = ScipyY.
        var profiles = new List<IReadOnlyDictionary<string, double>>();
        var groups = new List<int>();
        foreach (var v in ScipyX) { profiles.Add(new Dictionary<string, double> { ["T"] = v }); groups.Add(1); }
        foreach (var v in ScipyY) { profiles.Add(new Dictionary<string, double> { ["T"] = v }); groups.Add(2); }

        var results = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, pThreshold: 0.05);

        var t = results.Should().ContainSingle("exactly one taxon was supplied").Subject;
        t.Taxon.Should().Be("T");
        // FindSignificantTaxa reports U = max(U1,U2) = max(17,3) = 17 (§3.2).
        t.U.Should().BeApproximately(ScipyU1, 1e-9, "§3.2: reported U = max(U1,U2) = 17");
        t.Z.Should().BeApproximately(ScipyZContinuity, 1e-6, "§7.1: z ≈ 1.5922 (cc on by default)");
        t.PValue.Should().BeApproximately(ScipyPContinuity, 1e-6, "§7.1: p ≈ 0.1113");
        t.Significant.Should().BeFalse("p ≈ 0.1113 ≥ 0.05 ⇒ not significant — INV-04");
    }

    [Test]
    public void FindSignificantTaxa_SingleTaxonClearlySeparated_IsFlaggedSignificant()
    {
        // Group 1 abundances all strictly below group 2's ⇒ complete separation
        // ⇒ U = n1·n2, the most extreme rank-sum ⇒ small p, flagged significant.
        var profiles = new List<IReadOnlyDictionary<string, double>>();
        var groups = new List<int>();
        foreach (var v in new[] { 1.0, 2.0, 3.0, 4.0, 5.0 })
        { profiles.Add(new Dictionary<string, double> { ["X"] = v }); groups.Add(1); }
        foreach (var v in new[] { 10.0, 11.0, 12.0, 13.0, 14.0 })
        { profiles.Add(new Dictionary<string, double> { ["X"] = v }); groups.Add(2); }

        var t = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, pThreshold: 0.05)
            .Should().ContainSingle().Subject;

        // Complete separation: U2 = R2 − n2(n2+1)/2; with n1=n2=5, U = n1·n2 = 25.
        t.U.Should().BeApproximately(25.0, 1e-9, "complete separation ⇒ U = n1·n2 = 25");
        t.PValue.Should().BeGreaterThan(0.0).And.BeLessThan(0.05,
            "fully separated groups are significant at α = 0.05");
        t.Significant.Should().BeTrue("p < 0.05 ⇒ significant — INV-04");
    }

    #endregion

    #region META-TAXA-001 — BE boundary: EMPTY + argument validation

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY (BE). An empty profiles list ⇒ empty result, no crash.
    //   — §3.3, §6.1 ("Empty profiles list → empty result").
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void FindSignificantTaxa_EmptyProfiles_ReturnsEmptyResult()
    {
        var results = MetagenomicsAnalyzer.FindSignificantTaxa(
            Array.Empty<IReadOnlyDictionary<string, double>>(),
            Array.Empty<int>());

        results.Should().BeEmpty("§6.1: no taxa to test ⇒ empty result");
    }

    // EMPTY (BE): an empty group to the core test is undefined ⇒ ArgumentException.
    //   — §3.3, §6.1 ("Empty group → ArgumentException").
    [Test]
    public void MannWhitneyU_EmptyGroup_ThrowsArgumentException()
    {
        Action g1Empty = () => MetagenomicsAnalyzer.MannWhitneyU(Array.Empty<double>(), ScipyY);
        Action g2Empty = () => MetagenomicsAnalyzer.MannWhitneyU(ScipyX, Array.Empty<double>());

        g1Empty.Should().Throw<ArgumentException>("§6.1: n=0 group ⇒ test undefined");
        g2Empty.Should().Throw<ArgumentException>("§6.1: n=0 group ⇒ test undefined");
    }

    // EMPTY (BE): a profile with NO taxa keys is valid — a taxon absent from a
    // profile counts as abundance 0 (ASM-03), so the present taxon is still
    // tested against zeros, never crashing on the missing key.
    [Test]
    public void FindSignificantTaxa_AbsentTaxonTreatedAsZeroAbundance()
    {
        // Taxon "Z" appears only in group 1; in group 2 it is absent ⇒ 0.
        // Group 1 has high Z (10..14), group 2's Z is 0,0 ⇒ separation ⇒ significant.
        var profiles = new IReadOnlyDictionary<string, double>[]
        {
            new Dictionary<string, double> { ["Z"] = 10 },                 // g1
            new Dictionary<string, double> { ["Z"] = 12 },                 // g1
            new Dictionary<string, double> { ["Z"] = 14 },                 // g1
            new Dictionary<string, double>(),                               // g2: Z absent ⇒ 0
            new Dictionary<string, double> { ["other"] = 1 },              // g2: Z absent ⇒ 0
            new Dictionary<string, double>(),                               // g2: Z absent ⇒ 0
        };
        var groups = new[] { 1, 1, 1, 2, 2, 2 };

        var results = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, pThreshold: 0.05);

        var z = results.Single(t => t.Taxon == "Z");
        // group1 Z = {10,12,14} all > group2 Z = {0,0,0} ⇒ complete separation ⇒ U = 9.
        z.U.Should().BeApproximately(9.0, 1e-9, "ASM-03: absent ⇒ 0 ⇒ separation ⇒ U = n1·n2 = 9");
        z.PValue.Should().BeInRange(0.0, 1.0, "INV-03");
    }

    // Argument validation (BE / INJ surface): nulls and malformed labels must be
    // rejected eagerly and cleanly — never a raw NullReference / index error.
    [Test]
    public void MannWhitneyU_NullGroup_ThrowsArgumentNullException()
    {
        Action g1Null = () => MetagenomicsAnalyzer.MannWhitneyU(null!, ScipyY);
        Action g2Null = () => MetagenomicsAnalyzer.MannWhitneyU(ScipyX, null!);

        g1Null.Should().Throw<ArgumentNullException>("§6.1: null group ⇒ ArgumentNullException");
        g2Null.Should().Throw<ArgumentNullException>("§6.1: null group ⇒ ArgumentNullException");
    }

    [Test]
    public void FindSignificantTaxa_NullArguments_ThrowArgumentNullException()
    {
        var profiles = new IReadOnlyDictionary<string, double>[]
        {
            new Dictionary<string, double> { ["A"] = 1 },
        };
        Action nullProfiles = () => MetagenomicsAnalyzer.FindSignificantTaxa(null!, new[] { 1 });
        Action nullGroups = () => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, null!);

        nullProfiles.Should().Throw<ArgumentNullException>("§6.1: null profiles");
        nullGroups.Should().Throw<ArgumentNullException>("§6.1: null groups");
    }

    [Test]
    public void FindSignificantTaxa_LengthMismatch_ThrowsArgumentException()
    {
        var profiles = new IReadOnlyDictionary<string, double>[]
        {
            new Dictionary<string, double> { ["A"] = 1 },
            new Dictionary<string, double> { ["A"] = 2 },
        };
        // groups shorter than profiles
        Action act = () => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, new[] { 1 });

        act.Should().Throw<ArgumentException>("§3.3: profiles and groups must have equal length");
    }

    [Test]
    public void FindSignificantTaxa_InvalidGroupLabel_ThrowsArgumentException()
    {
        var profiles = new IReadOnlyDictionary<string, double>[]
        {
            new Dictionary<string, double> { ["A"] = 1 },
            new Dictionary<string, double> { ["A"] = 2 },
        };
        // 0 and -1 are out-of-range labels (BE: 0, -1) — only {1,2} allowed.
        Action zeroLabel = () => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, new[] { 0, 2 });
        Action negLabel = () => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, new[] { 1, -1 });
        Action bigLabel = () => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, new[] { 1, int.MaxValue });

        zeroLabel.Should().Throw<ArgumentException>("§6.1: label 0 ∉ {1,2}");
        negLabel.Should().Throw<ArgumentException>("§6.1: label -1 ∉ {1,2}");
        bigLabel.Should().Throw<ArgumentException>("§6.1: label MaxInt ∉ {1,2}");
    }

    [Test]
    public void FindSignificantTaxa_MissingOneGroup_ThrowsArgumentException()
    {
        var profiles = new IReadOnlyDictionary<string, double>[]
        {
            new Dictionary<string, double> { ["A"] = 1 },
            new Dictionary<string, double> { ["A"] = 2 },
        };
        // All profiles in group 1; group 2 is empty ⇒ no two-group comparison.
        Action act = () => MetagenomicsAnalyzer.FindSignificantTaxa(profiles, new[] { 1, 1 });

        act.Should().Throw<ArgumentException>("§6.1: both groups must contain ≥ 1 profile");
    }

    #endregion

    #region META-TAXA-001 — randomized robustness sweep (FindSignificantTaxa)

    // ───────────────────────────────────────────────────────────────────────
    // Broad boundary sweep over FindSignificantTaxa: random taxa, random sparse
    // profiles (so taxa are frequently absent ⇒ 0), random {1,2} labels with both
    // groups guaranteed non-empty. The driver must never crash, must produce one
    // result per observed taxon, p ∈ [0,1] (INV-03), finite z, ascending p-value
    // ordering, and Significant ⇔ p < threshold (INV-04).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void FindSignificantTaxa_RandomSparseProfiles_NeverCrashAndContractHolds(CancellationToken token)
    {
        var rng = new Random(197_003);
        string[] taxonPool = { "A", "B", "C", "D", "E", "F" };

        for (int trial = 0; trial < 1500 && !token.IsCancellationRequested; trial++)
        {
            int nSamples = rng.Next(2, 9);
            var profiles = new List<IReadOnlyDictionary<string, double>>(nSamples);
            var expectedTaxa = new HashSet<string>(StringComparer.Ordinal);

            for (int s = 0; s < nSamples; s++)
            {
                var map = new Dictionary<string, double>();
                foreach (var taxon in taxonPool)
                {
                    if (rng.Next(2) == 0) // sparse: present roughly half the time
                    {
                        double abundance = rng.Next(0, 6); // includes 0 abundances
                        map[taxon] = abundance;
                        expectedTaxa.Add(taxon);
                    }
                }
                profiles.Add(map);
            }

            // Build labels guaranteeing both groups are non-empty: first → 1, second → 2.
            var groups = new int[nSamples];
            groups[0] = 1;
            groups[1] = 2;
            for (int i = 2; i < nSamples; i++) groups[i] = rng.Next(1, 3);

            double pThreshold = rng.NextDouble(); // [0,1)

            var results = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, pThreshold);

            results.Should().HaveCount(expectedTaxa.Count,
                "one result per taxon observed in any profile");
            results.Select(r => r.PValue).Should().BeInAscendingOrder("§4.1 step 6: ascending p-value");

            foreach (var r in results)
            {
                r.PValue.Should().BeInRange(0.0, 1.0, "INV-03: p-value ∈ [0,1]");
                double.IsNaN(r.Z).Should().BeFalse("z must never be NaN");
                double.IsInfinity(r.Z).Should().BeFalse("z must never be Infinity");
                r.U.Should().BeGreaterThanOrEqualTo(0.0, "INV-02: reported U = max(U1,U2) ≥ 0");
                r.Significant.Should().Be(r.PValue < pThreshold, "INV-04: Significant ⇔ p < threshold");
            }
        }
    }

    #endregion
}
