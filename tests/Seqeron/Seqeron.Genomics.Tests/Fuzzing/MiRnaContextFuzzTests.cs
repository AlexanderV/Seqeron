using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the MiRNA area — the opt-in TargetScan <b>context++</b> target-site scorer
/// <see cref="MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(string, MiRnaAnalyzer.MiRna, MiRnaAnalyzer.TargetSite, MiRnaAnalyzer.ContextPlusPlusInputs)"/>
/// (Agarwal et al. 2015, eLife 4:e05005).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-CONTEXT-001 — context++ target-site scoring
/// Checklist: docs/checklists/03_FUZZING.md, row 252.
/// Source doc: docs/algorithms/MiRNA/Target_Site_Prediction.md (§5.1-§5.2, §7.1).
/// Source: src/.../Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs
///   • ScoreTargetSiteContextPlusPlus(...) (~764), ContextPlusPlusScore (~518),
///     ContextPlusPlusInputs (~552), and the per-feature contribution helpers.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// context++ is a multiple-linear-regression model: the score is the SUM of per-feature
/// contributions, each `coeff × scaledFeature`, with a distinct coefficient set per site
/// type (Agarwal 2015 [4][6]; targetscan_70_context_scores.pl [7]). The features fold in
/// `log10` transforms (Min_dist, Len_3UTR), min-max scaling, fixed-width window extraction
/// (Local_AU 30-nt flanks, SA 14-nt window, 3P_score subseq), and integer counts (Off6m) —
/// every one a hazard for a degenerate / boundary input to produce a NaN (log10 of a
/// non-positive number, divide-by-(max−min)), an Inf, or an IndexOutOfRange off a window
/// that walks past a too-short UTR. Fuzzing feeds those degenerate inputs and asserts that
/// the routine NEVER fails undisciplined: every input resolves to EITHER a finite,
/// theory-correct context++ score OR a documented validation exception (a non-canonical
/// site type ⇒ ArgumentException). The bar: the returned `ContextScorePartial` and every
/// component contribution are FINITE for any valid site, and equal the documented linear
/// combination of feature × coefficient.
///
/// Fuzz strategies exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate corners of the regression:
///        empty / null 3'UTR (every UTR-derived feature defaults to 0 with no crash),
///        a miRNA LONGER than the target (the 3P / Local_AU / window extractions must clamp,
///        never IndexOutOfRange), and a site placed at the VERY START / VERY END of the UTR
///        (the Min_dist edge: nearest-end distance 0 ⇒ log10 floored to 0, never log10(0) = −∞).
///        — docs/checklists/03_FUZZING.md §Description (code BE).
///   • MC = Malformed Content — a non-canonical / out-of-contract site type (Offset6mer,
///        Supplementary, Centered) ⇒ the documented ArgumentException, not a wrong score;
///        a "no seed match" UTR carrying no real site is scored only via a manufactured site
///        and must still produce a finite, defined score (the scorer trusts the supplied
///        site type and never re-discovers the seed). — §Description (code MC).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The context++ contract under test (Target_Site_Prediction.md §5.1-§5.2, §7.1)
/// ───────────────────────────────────────────────────────────────────────────
/// ScoreTargetSiteContextPlusPlus(mRnaSequence, miRna, site, inputs) → ContextPlusPlusScore.
///   • The site MUST be one of the four canonical seed-match types {Seed8mer, Seed7merM8,
///     Seed7merA1, Seed6mer}; any other type (Offset6mer / Supplementary / Centered) throws
///     ArgumentException (§5.1 / source ~770). This is the documented validation outcome.
///   • ContextScorePartial = Intercept + Local_AU + sRNA1 + sRNA8 + Site8 + SA + 3P_score
///       + Min_dist + Len_3UTR + Off6m + SPS + TA + Len_ORF + ORF8m + PCT.
///     Each feature is `coeff × (min-max scaled raw)` (or raw, for the Off6m/ORF8m counts and
///     the nucleotide-identity indicators), with the VERBATIM Agarwal_2015_parameters.txt
///     coefficients per site type [6].
///   • Optional features (SPS, TA, Len_ORF, ORF8m, PCT) default to 0 when the caller does NOT
///     supply them, and each missing one is named in OmittedFeatures (source BuildOmittedFeatures,
///     ~991). SA is omitted (contributes 0) when its 14-nt window does not fit the local UTR.
///   • More-negative context++ ⇒ stronger predicted repression. The site-type intercepts are
///     calibrated 8mer (−0.589) ≪ 7mer-m8 (−0.224) < 7mer-A1 (−0.195) < 6mer (−0.079) [6]:
///     in an otherwise feature-free context an 8mer scores more negative than a 6mer.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The hand-derived linear-combination pin (Target_Site_Prediction.md §7.1)
/// ───────────────────────────────────────────────────────────────────────────
/// The doc's published numerical walk-through: let-7a (seed GAGGUAG; nt1 = U, nt8 = G) against
/// mRNA = "GGGGG" + "CUACCUCA" + "GGGGG" (length 18) presents an 8mer at Start = 5, End = 12.
/// Both 30-nt flanks are all-G ⇒ local-AU fraction 0. Re-derived INDEPENDENTLY from the
/// Agarwal_2015_parameters.txt coefficients (NOT read off the code):
///   Intercept(8mer)            = -0.589
///   Local_AU = -0.254 × ((0 − 0.308)/(0.814 − 0.308))              = +0.15460869565217392
///   sRNA8 = G (8mer)                                                = +0.015
///       (nt1 = U ⇒ sRNA1 unscored; Site8 undefined for 8mer)
///   3P_score: raw = 0 ⇒ -0.040 × ((0 − 1)/(3.5 − 1))               = +0.016
///   Min_dist: nearest end = min(5, 18−13) = 5 ⇒
///             0.118 × ((log10 5 − 1.415)/(3.113 − 1.415))          = -0.049759446106213065
///   Len_3UTR: 0.310 × ((log10 18 − 2.392)/(3.637 − 2.392))         = -0.2830405810586145
///   Off6m: pattern CUACCU occurs once ⇒ -0.020 × 1                 = -0.020
///   SA: window [−1..12] does not fit (windowStart0 = −1 &lt; 0) ⇒ omitted, 0.
///   Sum = ContextScorePartial = -0.7561913315126536.
/// The test asserts EACH contribution field individually AND the partial sum, so a wrong
/// coefficient, a missing log10, a divide-by-(max−min) slip, or a feature dropped/double-counted
/// is caught — a test that passes against a wrong implementation is rejected.
///
/// All inputs are fixed / deterministically generated; the random helper uses a LOCALLY seeded
/// `new Random(seed)` (no shared static Rng), so every fuzz input is reproducible. The assembly
/// runs under LimitationMode.Permissive (LimitationPolicyTestBootstrap), so the partial-score
/// path (which enforces MIRNA-TARGET-001 when features are omitted) does not throw.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MiRnaContextFuzzTests
{
    #region Helpers

    // hsa-let-7a-5p — seed GAGGUAG (nt1 = U, nt8 = G). The doc's §7.1 worked-example miRNA.
    private const string Let7aSequence = "UGAGGUAGUAGGUUGUAUAGUU";

    // §7.1 worked example: all-G flanks around the 8mer site CUACCUCA → length-18 3'UTR.
    private const string WorkedExampleUtr = "GGGGGCUACCUCAGGGGG";
    private const int WorkedExampleStart = 5;   // 8mer occupies [5..12]
    private const int WorkedExampleEnd = 12;
    private const double WorkedExamplePartial = -0.7561913315126536;
    private const double Tol = 1e-12;

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomRna(int length, int seed)
    {
        const string bases = "ACGU";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal context++ finiteness invariant: the partial score and EVERY component
    /// contribution are finite (never NaN, never ±Inf) — the core fuzz bar for a regression model
    /// full of log10 transforms, min-max divisions and fixed-width window extractions.
    /// </summary>
    private static void AssertAllFinite(in ContextPlusPlusScore s)
    {
        foreach (var (value, name) in new[]
                 {
                     (s.Intercept, nameof(s.Intercept)),
                     (s.LocalAuContribution, nameof(s.LocalAuContribution)),
                     (s.SRna1Contribution, nameof(s.SRna1Contribution)),
                     (s.SRna8Contribution, nameof(s.SRna8Contribution)),
                     (s.Site8Contribution, nameof(s.Site8Contribution)),
                     (s.SaContribution, nameof(s.SaContribution)),
                     (s.ThreePrimePairingContribution, nameof(s.ThreePrimePairingContribution)),
                     (s.MinDistContribution, nameof(s.MinDistContribution)),
                     (s.Len3UtrContribution, nameof(s.Len3UtrContribution)),
                     (s.Off6mContribution, nameof(s.Off6mContribution)),
                     (s.SpsContribution, nameof(s.SpsContribution)),
                     (s.TaContribution, nameof(s.TaContribution)),
                     (s.LenOrfContribution, nameof(s.LenOrfContribution)),
                     (s.Orf8mContribution, nameof(s.Orf8mContribution)),
                     (s.PctContribution, nameof(s.PctContribution)),
                     (s.ContextScorePartial, nameof(s.ContextScorePartial)),
                 })
        {
            double.IsNaN(value).Should().BeFalse($"{name} must be finite (not NaN) for a valid site");
            double.IsInfinity(value).Should().BeFalse($"{name} must be finite (not ±Inf) for a valid site");
        }

        // The partial sum equals the sum of its parts — the regression is additive (no hidden term).
        double recomputed = s.Intercept + s.LocalAuContribution + s.SRna1Contribution
            + s.SRna8Contribution + s.Site8Contribution + s.SaContribution
            + s.ThreePrimePairingContribution + s.MinDistContribution + s.Len3UtrContribution
            + s.Off6mContribution + s.SpsContribution + s.TaContribution
            + s.LenOrfContribution + s.Orf8mContribution + s.PctContribution;
        s.ContextScorePartial.Should().BeApproximately(recomputed, Tol,
            "ContextScorePartial is the additive sum of all per-feature contributions");
    }

    /// <summary>Builds a canonical seed-match site of a given type at the supplied span (the scorer
    /// trusts the type — it does not re-discover the seed — so a manufactured site is a legal probe).</summary>
    private static TargetSite Site(TargetSiteType type, int start, int end, int seedLen) =>
        new TargetSite(
            Start: start, End: end, TargetSequence: "", MiRnaName: "let-7a",
            Type: type, SeedMatchLength: seedLen, Score: 0.0, FreeEnergy: 0.0, Alignment: "");

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  MIRNA-CONTEXT-001 — context++ target-site scoring : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region MIRNA-CONTEXT-001 — context++ scoring

    #region Pin — the §7.1 hand-derived linear combination

    /// <summary>
    /// THE anchor test: the doc's §7.1 published numerical walk-through, re-derived INDEPENDENTLY
    /// from the Agarwal_2015_parameters.txt coefficients (see the class-header derivation — NOT read
    /// off the code). let-7a vs the all-G-flanked 8mer UTR yields the per-feature contributions and
    /// the partial sum below. Asserting EACH field plus the sum catches a wrong coefficient, a
    /// dropped log10, a min-max divide slip, or a feature double-counted — proving the harness pins
    /// the real regression, not a tautology.
    /// </summary>
    [Test]
    public void ScoreContextPlusPlus_WorkedExample_MatchesHandDerivedLinearCombination()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);
        let7a.SeedSequence.Should().Be("GAGGUAG", "the pin starts from the canonical let-7a seed (nt1=U, nt8=G)");

        // Sanity: the real finder discovers exactly this 8mer at the documented span.
        var found = FindTargetSites(WorkedExampleUtr, let7a, minScore: 0.1).ToList();
        var site = found.Should().ContainSingle(s => s.Type == TargetSiteType.Seed8mer).Subject;
        site.Start.Should().Be(WorkedExampleStart, "the worked-example 8mer occupies [5..12] (§7.1)");
        site.End.Should().Be(WorkedExampleEnd);

        var s = ScoreTargetSiteContextPlusPlus(WorkedExampleUtr, let7a, site);

        AssertAllFinite(s);
        s.SiteType.Should().Be(TargetSiteType.Seed8mer);

        // Each contribution, hand-derived from Agarwal_2015_parameters.txt (§7.1).
        s.Intercept.Should().BeApproximately(-0.589, Tol, "8mer site-type intercept");
        s.LocalAuContribution.Should().BeApproximately(0.15460869565217392, Tol,
            "Local_AU = -0.254 × ((0 − 0.308)/(0.814 − 0.308)); all-G flanks ⇒ AU fraction 0");
        s.SRna1Contribution.Should().Be(0.0, "miRNA nt1 = U ⇒ sRNA1 indicators unscored");
        s.SRna8Contribution.Should().BeApproximately(0.015, Tol, "miRNA nt8 = G ⇒ sRNA8G(8mer) = +0.015");
        s.Site8Contribution.Should().Be(0.0, "Site8 is defined only for 7mer-A1 / 6mer, not 8mer");
        s.SaContribution.Should().Be(0.0, "the 14-nt SA window [−1..12] does not fit ⇒ SA omitted (0)");
        s.ThreePrimePairingContribution.Should().BeApproximately(0.016, Tol,
            "3P_score raw 0 ⇒ -0.040 × ((0 − 1)/(3.5 − 1))");
        s.MinDistContribution.Should().BeApproximately(-0.049759446106213065, Tol,
            "Min_dist nearest end = min(5, 18−13) = 5 ⇒ 0.118 × ((log10 5 − 1.415)/(3.113 − 1.415))");
        s.Len3UtrContribution.Should().BeApproximately(-0.2830405810586145, Tol,
            "Len_3UTR = 0.310 × ((log10 18 − 2.392)/(3.637 − 2.392))");
        s.Off6mContribution.Should().BeApproximately(-0.020, Tol, "Off6m: pattern CUACCU occurs once ⇒ -0.020 × 1");

        // Optional, caller-only features default to 0 and stay residual.
        s.SpsContribution.Should().Be(0.0, "SPS not supplied ⇒ 0");
        s.TaContribution.Should().Be(0.0, "TA not supplied ⇒ 0");
        s.LenOrfContribution.Should().Be(0.0, "Len_ORF not supplied ⇒ 0");
        s.Orf8mContribution.Should().Be(0.0, "ORF8m not supplied ⇒ 0");
        s.PctContribution.Should().Be(0.0, "no Conservation input ⇒ PCT residual ⇒ 0");

        // The published partial sum.
        s.ContextScorePartial.Should().BeApproximately(WorkedExamplePartial, Tol,
            "the partial context++ score equals the documented §7.1 sum -0.7561913315126536");
    }

    /// <summary>
    /// The omitted-feature contract (BuildOmittedFeatures, source ~991): for the worked-example call
    /// (no optional inputs) PCT, SPS, TA_3UTR, Len_ORF and ORF8m are ALL residual and named in
    /// OmittedFeatures, while the realised features (Local_AU, 3P, Min_dist, Len_3UTR, Off6m) are
    /// folded in. SA is also omitted here (its window does not fit). Missing optional features
    /// default to a 0 CONTRIBUTION — asserted in the pin above — never to a guessed value.
    /// </summary>
    [Test]
    public void ScoreContextPlusPlus_NoOptionalInputs_ListsResidualFeaturesAsOmitted()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);
        var site = Site(TargetSiteType.Seed8mer, WorkedExampleStart, WorkedExampleEnd, 8);

        var s = ScoreTargetSiteContextPlusPlus(WorkedExampleUtr, let7a, site);

        s.OmittedFeatures.Should().NotBeNull();
        // Each data-/parameter-blocked feature is honestly reported as omitted.
        s.OmittedFeatures.Should().Contain(f => f.StartsWith("PCT"), "no conservation ⇒ PCT omitted");
        s.OmittedFeatures.Should().Contain(f => f.StartsWith("SPS"), "no SPS supplied ⇒ omitted");
        s.OmittedFeatures.Should().Contain(f => f.StartsWith("TA_3UTR"), "no TA supplied ⇒ omitted");
        s.OmittedFeatures.Should().Contain(f => f.StartsWith("Len_ORF"), "no ORF length ⇒ omitted");
        s.OmittedFeatures.Should().Contain(f => f.StartsWith("ORF8m"), "no ORF-8mer count ⇒ omitted");
    }

    /// <summary>
    /// Missing-optional-feature DEFAULT vs SUPPLIED: supplying an optional feature must add EXACTLY
    /// its documented `coeff × scaled(raw)` to the partial score and remove it from OmittedFeatures;
    /// omitting it must contribute exactly 0. We supply TA = log10(N) for an 8mer and assert the
    /// delta equals the hand-derived TA contribution 0.222 × ((TA − 3.113)/(3.865 − 3.113)).
    /// </summary>
    [Test]
    public void ScoreContextPlusPlus_SuppliedOptionalFeature_AddsExactlyItsDocumentedContribution()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);
        var site = Site(TargetSiteType.Seed8mer, WorkedExampleStart, WorkedExampleEnd, 8);

        var baseline = ScoreTargetSiteContextPlusPlus(WorkedExampleUtr, let7a, site);
        baseline.TaContribution.Should().Be(0.0, "TA omitted ⇒ 0 contribution (the documented default)");

        const double ta = 3.5; // an in-range log10(site-count) value
        var withTa = ScoreTargetSiteContextPlusPlus(
            WorkedExampleUtr, let7a, site, new ContextPlusPlusInputs(Ta: ta));

        // Hand-derived TA(8mer) = 0.222 × ((3.5 − 3.113)/(3.865 − 3.113)).
        double expectedTa = 0.222 * ((ta - 3.113) / (3.865 - 3.113));
        withTa.TaContribution.Should().BeApproximately(expectedTa, Tol,
            "supplied TA folds in as 0.222 × ((TA − 3.113)/(3.865 − 3.113)) for an 8mer");
        withTa.OmittedFeatures.Should().NotContain(f => f.StartsWith("TA_3UTR"),
            "a supplied TA is no longer residual");

        // The ONLY change to the partial score is the TA term (every other feature is unchanged).
        (withTa.ContextScorePartial - baseline.ContextScorePartial)
            .Should().BeApproximately(expectedTa, Tol,
            "supplying TA shifts the partial score by exactly its contribution, nothing else");
        AssertAllFinite(withTa);
    }

    /// <summary>
    /// Calibration ordering (Agarwal 2015 intercepts [6]): in an otherwise feature-neutral context a
    /// STRONGER site type scores MORE NEGATIVE (stronger predicted repression) than a weaker one.
    /// We score the SAME span/UTR as an 8mer and as a 6mer; the 8mer's partial score must be the
    /// more negative of the two (its intercept −0.589 ≪ the 6mer's −0.079). This pins the doc's
    /// site-type calibration, not just the arithmetic of one site.
    /// </summary>
    [Test]
    public void ScoreContextPlusPlus_StrongerSiteType_ScoresMoreNegativeThanWeaker()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);

        // A neutral middle-of-UTR span in a long all-G UTR (no AU, no Off6m pattern, symmetric ends).
        string utr = new string('G', 81);
        const int start = 37, end8 = 44, end6 = 42;

        var eightMer = ScoreTargetSiteContextPlusPlus(utr, let7a, Site(TargetSiteType.Seed8mer, start, end8, 8));
        var sixMer = ScoreTargetSiteContextPlusPlus(utr, let7a, Site(TargetSiteType.Seed6mer, start, end6, 6));

        AssertAllFinite(eightMer);
        AssertAllFinite(sixMer);
        eightMer.Intercept.Should().BeApproximately(-0.589, Tol, "8mer intercept");
        sixMer.Intercept.Should().BeApproximately(-0.079, Tol, "6mer intercept");
        eightMer.ContextScorePartial.Should().BeLessThan(sixMer.ContextScorePartial,
            "an 8mer is calibrated to stronger (more negative) repression than a 6mer (§5.2, [6])");
    }

    #endregion

    #region MC — Malformed Content: non-canonical site type ⇒ documented ArgumentException

    /// <summary>
    /// MC — the documented validation exception. context++ is defined only for the four canonical
    /// seed-match types; any other site type (Offset6mer, Supplementary, Centered) must throw
    /// ArgumentException naming the `site` argument (source ~770) — a defined refusal, NOT a wrong
    /// score, NaN, or crash. We sweep all three out-of-contract types.
    /// </summary>
    [Test]
    public void ScoreContextPlusPlus_NonCanonicalSiteType_ThrowsArgumentException()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);

        foreach (var bad in new[]
                 {
                     TargetSiteType.Offset6mer,
                     TargetSiteType.Supplementary,
                     TargetSiteType.Centered,
                 })
        {
            var site = Site(bad, WorkedExampleStart, WorkedExampleEnd, 6);
            var act = () => ScoreTargetSiteContextPlusPlus(WorkedExampleUtr, let7a, site);
            act.Should().Throw<ArgumentException>(
                    $"{bad} is not a canonical context++ site type ⇒ documented validation exception")
               .Which.ParamName.Should().Be("site", "the rejected argument is the site");
        }

        // Sanity: a canonical type at the same span is accepted (the guard is type-specific, not blanket).
        var ok = () => ScoreTargetSiteContextPlusPlus(
            WorkedExampleUtr, let7a, Site(TargetSiteType.Seed6mer, WorkedExampleStart, WorkedExampleEnd, 6));
        ok.Should().NotThrow("a canonical 6mer site is a valid context++ input");
    }

    #endregion

    #region BE — Boundary: empty / null 3'UTR

    /// <summary>
    /// BE — "empty 3'UTR". When mRnaSequence is empty (or null) the scorer normalises it to "" so
    /// every UTR-DERIVED feature (Local_AU, Site8, SA, 3P_score, Min_dist, Off6m) defaults to 0, and
    /// only the miRNA-derived terms (Intercept, sRNA1/8) and Len_3UTR (log10 of length 0 ⇒ floored
    /// to log10 0 = 0 input) remain. The score must be FINITE — no log10(0) = −∞, no divide-by-zero,
    /// no IndexOutOfRange — for both the empty and null UTR. We assert finiteness and that the
    /// UTR-derived contributions are exactly 0.
    /// </summary>
    [Test]
    public void ScoreContextPlusPlus_EmptyOrNullUtr_AllUtrFeaturesZeroAndScoreFinite()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);
        var site = Site(TargetSiteType.Seed8mer, 0, 7, 8);

        foreach (string? utr in new[] { "", (string?)null })
        {
            var s = ScoreTargetSiteContextPlusPlus(utr!, let7a, site);
            AssertAllFinite(s);

            s.LocalAuContribution.Should().Be(0.0, "empty UTR ⇒ Local_AU flank weight 0 ⇒ 0");
            s.Site8Contribution.Should().Be(0.0, "empty UTR ⇒ Site8 base out of range ⇒ 0");
            s.SaContribution.Should().Be(0.0, "empty UTR ⇒ SA window cannot fit ⇒ 0");
            // 3P_score RAW is 0 for an empty UTR, but the min-max scaled CONTRIBUTION of raw 0 is
            // -0.040 × ((0 − 1)/(3.5 − 1)) = +0.016 (Ctx3PMin = 1 ≠ 0): the contribution is not 0.
            s.ThreePrimePairingContribution.Should().BeApproximately(0.016, Tol,
                "empty UTR ⇒ 3P_score raw 0 ⇒ scaled contribution -0.040 × ((0 − 1)/(3.5 − 1)) = +0.016");
            s.MinDistContribution.Should().Be(0.0, "empty UTR ⇒ Min_dist 0 (mrna.Length 0 short-circuit)");
            s.Off6mContribution.Should().Be(0.0, "empty UTR ⇒ Off6m count 0 ⇒ 0");
            // Intercept and sRNA8 (G) still apply: the score is not trivially 0.
            s.Intercept.Should().BeApproximately(-0.589, Tol, "the 8mer intercept is miRNA/type-derived");
            s.SRna8Contribution.Should().BeApproximately(0.015, Tol, "sRNA8 G is miRNA-derived ⇒ unaffected");
        }
    }

    #endregion

    #region BE — Boundary: miRNA longer than the target

    /// <summary>
    /// BE — "miRNA longer than target". The 30-nt Local_AU flanks, the 14-nt SA window and the
    /// 3P_score subseq are all wider than a tiny UTR; a naive extraction would walk off the end. The
    /// scorer must clamp every window to the available UTR and never IndexOutOfRange / NaN. We sweep
    /// UTRs FAR SHORTER than the 22-nt let-7a (down to length 1, with a 6mer site clamped inside),
    /// asserting a finite, additive score every time.
    /// </summary>
    [Test]
    [CancelAfter(30_000)]
    public void ScoreContextPlusPlus_MiRnaLongerThanTarget_FiniteNeverIndexOverflow()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence); // 22 nt
        let7a.Sequence.Length.Should().BeGreaterThan(7, "the probe miRNA is longer than the tiny UTRs below");

        foreach (int len in new[] { 1, 2, 3, 4, 6, 8, 10 })
        {
            string tinyUtr = new string('C', len);
            // A 6mer site fully inside the tiny UTR (End ≤ len-1).
            int end = Math.Min(len - 1, len >= 6 ? 5 : len - 1);
            int start = Math.Max(0, end - 5);
            var site = Site(TargetSiteType.Seed6mer, start, end, 6);

            var act = () => ScoreTargetSiteContextPlusPlus(tinyUtr, let7a, site);
            var s = act.Should().NotThrow(
                $"a {len}-nt UTR far shorter than the 22-nt miRNA must not crash the scorer").Subject;
            AssertAllFinite(s);
        }
    }

    #endregion

    #region BE — Boundary: site at the very start / very end of the UTR (Min_dist edge)

    /// <summary>
    /// BE — site at the VERY START / VERY END of the 3'UTR: the Min_dist log10 edge. Min_dist takes
    /// log10(min(distTo5′, distTo3′)). When the site abuts an end the nearest distance is 0, and a
    /// naive log10(0) = −∞ would make the whole score −∞ / NaN. The source floors the log10 input to
    /// 0 for a nearest distance of 0 (source ~1429), so the contribution stays finite. We place a
    /// site flush against each end and assert: Min_dist is finite, equals the documented value for a
    /// nearest distance of 0 (log10 floored to 0), and the whole score is finite.
    /// </summary>
    [Test]
    public void ScoreContextPlusPlus_SiteAtUtrEnd_MinDistLog10FlooredFinite()
    {
        var let7a = CreateMiRna("let-7a", Let7aSequence);
        string utr = new string('G', 40);

        // 8mer flush against the 5' end: distTo5′ = 0 ⇒ nearest = 0 ⇒ log10 input floored to 0.
        var atStart = ScoreTargetSiteContextPlusPlus(utr, let7a, Site(TargetSiteType.Seed8mer, 0, 7, 8));
        // 8mer flush against the 3' end: distTo3′ = 0 ⇒ nearest = 0 likewise.
        var atEnd = ScoreTargetSiteContextPlusPlus(utr, let7a, Site(TargetSiteType.Seed8mer, 32, 39, 8));

        // Documented Min_dist(8mer) at nearest = 0: 0.118 × ((0 − 1.415)/(3.113 − 1.415)).
        double expectedMinDistAtEnd = 0.118 * ((0 - 1.415) / (3.113 - 1.415));

        foreach (var (s, where) in new[] { (atStart, "5' end"), (atEnd, "3' end") })
        {
            AssertAllFinite(s);
            s.MinDistContribution.Should().BeApproximately(expectedMinDistAtEnd, Tol,
                $"site flush against the {where}: nearest distance 0 ⇒ log10 floored to 0, finite Min_dist");
            double.IsNegativeInfinity(s.ContextScorePartial).Should().BeFalse(
                $"the {where} site must not produce a −∞ score from log10(0)");
        }
    }

    #endregion

    #region MC/BE — robustness sweep: random miRNA / UTR / type / span never NaN, Inf or crash

    /// <summary>
    /// Broad robustness fuzz: across fixed seeds, random RNA miRNAs and UTRs, all four canonical site
    /// types, and random in-range spans, ScoreTargetSiteContextPlusPlus must NEVER throw (other than
    /// the documented type guard, which we avoid by only using canonical types) and ALWAYS return a
    /// finite, additive score. This is the §8 "no undisciplined failure" bar over a wide input space.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void ScoreContextPlusPlus_RandomInputs_AlwaysFiniteAndCrashFree()
    {
        var types = new[]
        {
            TargetSiteType.Seed8mer, TargetSiteType.Seed7merM8,
            TargetSiteType.Seed7merA1, TargetSiteType.Seed6mer,
        };

        foreach (int seed in new[] { 1, 7, 42, 2026 })
        {
            var rng = new Random(seed);
            foreach (int utrLen in new[] { 6, 20, 50, 120 })
            {
                string utr = RandomRna(utrLen, seed + utrLen);
                string mirnaSeq = RandomRna(rng.Next(8, 25), seed * 31 + utrLen);
                var mirna = CreateMiRna($"rnd-{seed}-{utrLen}", mirnaSeq);

                foreach (var type in types)
                {
                    int siteLen = type == TargetSiteType.Seed6mer || type == TargetSiteType.Seed7merA1 ? 7 : 8;
                    if (utrLen < siteLen) continue;
                    int start = rng.Next(0, utrLen - siteLen + 1);
                    int end = start + siteLen - 1;
                    var site = Site(type, start, end, type == TargetSiteType.Seed8mer ? 8
                                                     : type == TargetSiteType.Seed6mer ? 6 : 7);

                    var act = () => ScoreTargetSiteContextPlusPlus(utr, mirna, site);
                    var s = act.Should().NotThrow(
                        $"random valid input must not crash (seed {seed}, utrLen {utrLen}, type {type})").Subject;
                    AssertAllFinite(s);
                }
            }
        }
    }

    #endregion

    #endregion
}
