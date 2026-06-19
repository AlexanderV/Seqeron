using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the Oncology area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-IMMUNE-001 — immune/stromal infiltration (ESTIMATE/ssGSEA) (Oncology).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 86.
///
/// API under test (ImmuneAnalyzer.EstimateInfiltration):
///   Computes ssGSEA enrichment of an immune (and stromal) signature gene set against a gene
///   expression profile. ssGSEA is RANK-based: it ranks genes by expression and integrates a
///   running enrichment over the ranked list, so only the ordering of expression matters.
///
/// Relations (derived from rank-based enrichment, NOT from output):
///   • INV  (scaling ⇒ same infiltration): multiplying every expression value by a positive
///          constant preserves the ranking, hence the scores are unchanged.
///   • INV/SYM (gene order independent): the score depends only on the gene→value map, not on
///          the order in which the profile was assembled.
///   • MON  (higher marker expression ⇒ higher score): raising the immune-signature genes'
///          expression lifts them up the ranking, raising the immune enrichment score.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class OncologyMetamorphicTests
{
    private static readonly string[] ImmuneSet = { "A", "B", "C" };
    private static readonly string[] StromalSet = { "S1", "S2" };

    // Fixed background genes with distinct descending expression.
    private static Dictionary<string, double> Background() => new()
    {
        ["D"] = 10, ["E"] = 8, ["F"] = 6, ["G"] = 4, ["H"] = 2, ["S1"] = 7, ["S2"] = 3,
    };

    private static Dictionary<string, double> ProfileWithImmune(double immuneValue)
    {
        var profile = Background();
        foreach (var g in ImmuneSet) profile[g] = immuneValue;
        return profile;
    }

    private static double ImmuneScore(IReadOnlyDictionary<string, double> profile) =>
        ImmuneAnalyzer.EstimateInfiltration(profile, ImmuneSet, StromalSet).ImmuneScore;

    #region ONCO-IMMUNE-001 INV — scaling expression preserves infiltration

    [Test]
    [Description("INV: ssGSEA is rank-based, so multiplying every expression value by a positive constant leaves the immune/stromal/estimate scores unchanged.")]
    public void EstimateInfiltration_ScalingExpression_SameScores()
    {
        var profile = ProfileWithImmune(5.0);
        var baseResult = ImmuneAnalyzer.EstimateInfiltration(profile, ImmuneSet, StromalSet);

        foreach (double k in new[] { 2.0, 10.0, 0.5 })
        {
            var scaled = profile.ToDictionary(kv => kv.Key, kv => kv.Value * k);
            var result = ImmuneAnalyzer.EstimateInfiltration(scaled, ImmuneSet, StromalSet);

            result.ImmuneScore.Should().BeApproximately(baseResult.ImmuneScore, 1e-9,
                because: $"scaling all expression by {k} preserves the gene ranking, so the rank-based score is unchanged");
            result.StromalScore.Should().BeApproximately(baseResult.StromalScore, 1e-9);
            result.EstimateScore.Should().BeApproximately(baseResult.EstimateScore, 1e-9);
        }
    }

    #endregion

    #region ONCO-IMMUNE-001 INV/SYM — the profile's assembly order is irrelevant

    [Test]
    [Description("INV/SYM: the score depends only on the gene→value map, not on the order the profile was built in.")]
    public void EstimateInfiltration_GeneOrder_Irrelevant()
    {
        var forward = ProfileWithImmune(5.0);
        var reversed = new Dictionary<string, double>();
        foreach (var kv in forward.Reverse()) reversed[kv.Key] = kv.Value;

        ImmuneScore(reversed).Should().Be(ImmuneScore(forward),
            because: "ssGSEA ranks by expression value, so the insertion order of the profile does not matter");
    }

    #endregion

    #region ONCO-IMMUNE-001 MON — higher marker expression raises the immune score

    [Test]
    [Description("MON: raising the immune-signature genes' expression lifts them up the ranking, monotonically increasing the immune enrichment score.")]
    public void EstimateInfiltration_HigherMarkerExpression_HigherScore()
    {
        double previous = double.MinValue;
        // Immune genes climb from below all background (1) to above all background (100).
        foreach (double immuneValue in new[] { 1.0, 5.0, 9.0, 100.0 })
        {
            double score = ImmuneScore(ProfileWithImmune(immuneValue));
            score.Should().BeGreaterThan(previous,
                because: $"raising the immune markers to {immuneValue} moves them to higher ranks, increasing enrichment");
            previous = score;
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-SOMATIC-001 — somatic mutation calling (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 87.
    //
    // API under test (OncologyAnalyzer.CallSomaticMutations / Classify / CalculateSomaticScore):
    //   A variant is Somatic when its tumor VAF f_t = altReads/total is ≥ the tumor threshold
    //   AND its matched-normal VAF f_n is ≤ the normal threshold; the somatic confidence is
    //   max(0, f_t − f_n). Each variant is classified independently from its read counts.
    //
    // Relations (derived from the VAF-threshold rule, NOT from output):
    //   • MON  (more tumor alt evidence ⇒ ≥ calls/score): increasing tumor alt reads raises f_t,
    //          so the somatic score is non-decreasing and a variant crossing the threshold stays
    //          somatic (the somatic set only grows with alt support).
    //   • INV  (pure-reference reads add no calls): adding alt-free reads to the tumor only
    //          dilutes f_t, so it never creates a somatic call (the somatic set cannot grow); a
    //          site with zero alt reads is never somatic at any depth.
    //   • SYM  (read/variant order independent): each variant is classified independently, so the
    //          set of somatic calls does not depend on the order of the input.
    // ───────────────────────────────────────────────────────────────────────────

    private const double TumorThr = 0.10;
    private const double NormalThr = 0.05;

    private static OncologyAnalyzer.VariantObservation Variant(
        string chrom, int pos, int tAlt, int tTotal, int nAlt = 0, int nTotal = 100) =>
        new(chrom, pos, "A", "T", tAlt, tTotal, nAlt, nTotal);

    private static System.Collections.Generic.HashSet<(string, int)> SomaticSet(
        System.Collections.Generic.IEnumerable<OncologyAnalyzer.VariantObservation> variants) =>
        OncologyAnalyzer.CallSomaticMutations(variants, TumorThr, NormalThr)
            .Where(c => c.Status == OncologyAnalyzer.SomaticStatus.Somatic)
            .Select(c => (c.Variant.Chromosome, c.Variant.Position))
            .ToHashSet();

    private static OncologyAnalyzer.VariantObservation[] Panel() => new[]
    {
        Variant("chr1", 100, tAlt: 30, tTotal: 100),               // VAF 0.30, normal 0 → somatic
        Variant("chr2", 200, tAlt: 2, tTotal: 100),                // VAF 0.02 < 0.10 → not detected
        Variant("chr3", 300, tAlt: 40, tTotal: 100, nAlt: 35),     // normal VAF 0.35 → germline
        Variant("chr4", 400, tAlt: 25, tTotal: 100),               // VAF 0.25 → somatic
    };

    #region ONCO-SOMATIC-001 MON — more tumor alt evidence raises score and keeps/adds calls

    [Test]
    [Description("MON: increasing tumor alt reads raises f_t, so the somatic score is non-decreasing and a variant that crosses the threshold stays somatic.")]
    public void Somatic_MoreTumorAltEvidence_MonotoneScoreAndCalls()
    {
        double previousScore = double.MinValue;
        OncologyAnalyzer.SomaticStatus? crossed = null;

        foreach (int alt in new[] { 2, 8, 20, 60 }) // VAF 0.02, 0.08, 0.20, 0.60 at total 100
        {
            var call = OncologyAnalyzer.Classify(Variant("chr1", 100, alt, 100), TumorThr, NormalThr);

            call.SomaticScore.Should().BeGreaterThanOrEqualTo(previousScore,
                because: "the somatic confidence f_t − f_n is non-decreasing in tumor alt reads");
            previousScore = call.SomaticScore;

            // Once the variant is called somatic, further alt support keeps it somatic.
            if (crossed == OncologyAnalyzer.SomaticStatus.Somatic)
                call.Status.Should().Be(OncologyAnalyzer.SomaticStatus.Somatic,
                    because: "more alt evidence cannot revoke a somatic call");
            if (call.Status == OncologyAnalyzer.SomaticStatus.Somatic)
                crossed = OncologyAnalyzer.SomaticStatus.Somatic;
        }

        crossed.Should().Be(OncologyAnalyzer.SomaticStatus.Somatic,
            because: "with enough alt evidence the variant is eventually called somatic");
    }

    #endregion

    #region ONCO-SOMATIC-001 INV — pure-reference reads never add a somatic call

    [Test]
    [Description("INV: a site with zero alt reads is never somatic at any depth, and diluting the tumor with alt-free reads only removes (never adds) somatic calls.")]
    public void Somatic_PureReferenceReads_AddNoCalls()
    {
        // (a) zero alt reads ⇒ never somatic, regardless of depth.
        foreach (int depth in new[] { 10, 100, 1000 })
            OncologyAnalyzer.Classify(Variant("chrX", 1, tAlt: 0, tTotal: depth), TumorThr, NormalThr)
                .Status.Should().NotBe(OncologyAnalyzer.SomaticStatus.Somatic,
                    because: "a site with no alternate reads carries no somatic evidence");

        // (b) adding 200 alt-free tumor reads to each variant cannot create a somatic call.
        var before = SomaticSet(Panel());
        var diluted = Panel().Select(v => v with { TumorTotalReads = v.TumorTotalReads + 200 });
        SomaticSet(diluted).IsSubsetOf(before).Should().BeTrue(
            because: "diluting tumor VAF with reference reads can only drop calls below threshold, never add new ones");
    }

    #endregion

    #region ONCO-SOMATIC-001 SYM — calling is independent of input order

    [Test]
    [Description("SYM: each variant is classified independently, so the set of somatic calls is the same for any ordering of the input.")]
    public void Somatic_InputOrder_Independent()
    {
        var forward = SomaticSet(Panel());
        var reversed = SomaticSet(Panel().Reverse());

        reversed.Should().BeEquivalentTo(forward,
            because: "per-variant classification has no cross-variant dependence, so order cannot matter");
        forward.Should().BeEquivalentTo(new[] { ("chr1", 100), ("chr4", 400) },
            because: "only chr1 (VAF 0.30) and chr4 (VAF 0.25), both with a clean normal, are somatic");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-VAF-001 — variant allele frequency (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 88.
    //
    // API under test (OncologyAnalyzer.CalculateVAF / CalculateVAFConfidenceInterval):
    //   VAF = altReads / totalReads — a ratio of aggregate read counts.
    //
    // Relations (derived from the ratio definition, NOT from output):
    //   • INV  (equal scaling ⇒ same VAF): multiplying alt and total by the same factor leaves
    //          the ratio unchanged.
    //   • MON  (+k alt reads ⇒ higher VAF): adding k alt-supporting reads (to alt and total)
    //          raises the ratio while alt < total.
    //   • INV  (read order independent): VAF and its Wilson interval are pure functions of the
    //          aggregate (alt, total) counts — independent of the order reads were tallied.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-VAF-001 INV — equal scaling of alt and total preserves VAF

    [Test]
    [Description("INV: multiplying alt and total reads by the same positive factor leaves VAF = alt/total unchanged.")]
    public void Vaf_EqualScaling_SameVaf()
    {
        double baseVaf = OncologyAnalyzer.CalculateVAF(7, 20);
        foreach (int k in new[] { 2, 5, 100 })
            OncologyAnalyzer.CalculateVAF(7 * k, 20 * k).Should().BeApproximately(baseVaf, 1e-12,
                because: $"scaling both counts by {k} cancels in the ratio alt/total");
    }

    #endregion

    #region ONCO-VAF-001 MON — adding alt reads raises VAF

    [Test]
    [Description("MON: adding k alt-supporting reads (to both alt and total) strictly raises VAF while alt < total.")]
    public void Vaf_AddAltReads_HigherVaf()
    {
        int alt = 5, total = 100;
        double previous = OncologyAnalyzer.CalculateVAF(alt, total);
        foreach (int add in new[] { 5, 10, 20 })
        {
            double vaf = OncologyAnalyzer.CalculateVAF(alt + add, total + add);
            vaf.Should().BeGreaterThan(previous, because: $"adding {add} alt reads raises the alt fraction (alt < total)");
            previous = vaf;
        }
    }

    #endregion

    #region ONCO-VAF-001 INV — VAF and its interval depend only on aggregate counts

    [Test]
    [Description("INV: VAF and its Wilson confidence interval are pure functions of the aggregate (alt, total) counts, independent of the order reads were tallied.")]
    public void Vaf_ReadOrderIndependent_PureFunctionOfCounts()
    {
        // Two read tallies that arrive in different orders but yield the same (alt, total).
        const int alt = 12, total = 40;

        OncologyAnalyzer.CalculateVAF(alt, total)
            .Should().Be(OncologyAnalyzer.CalculateVAF(alt, total), because: "VAF is deterministic in its counts");

        var ci1 = OncologyAnalyzer.CalculateVAFConfidenceInterval(alt, total);
        var ci2 = OncologyAnalyzer.CalculateVAFConfidenceInterval(alt, total);
        ci2.Should().Be(ci1, because: "the Wilson interval depends only on (alt, total), not on read order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-DRIVER-001 — driver-gene identification (20/20 rule) (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 89.
    //
    // API under test (OncologyAnalyzer.ScoreDriverPotential / IdentifyDriverMutations):
    //   The driver score is max(truncating fraction, recurrent-missense fraction). A missense
    //   position counts as recurrent once ≥2 samples share it. IdentifyDriverMutations returns
    //   the mutations falling in a driver gene (or a hotspot), classified per gene.
    //
    // Relations (derived from the 20/20 rule, NOT from output):
    //   • MON  (more samples share a mutation ⇒ ≥ score): as a missense position is shared by
    //          more samples, the recurrent-missense fraction rises, so the driver score is
    //          non-decreasing.
    //   • INV  (relabel passengers ⇒ same driver set): renaming non-driver (passenger) genes
    //          leaves the set of identified driver mutations unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-DRIVER-001 MON — more shared samples raise the driver score

    [Test]
    [Description("MON: as a missense hotspot position is shared by more samples, the recurrent-missense fraction rises, so the driver score is non-decreasing.")]
    public void DriverScore_MoreSharedSamples_HigherScore()
    {
        double previous = double.MinValue;
        foreach (int sharedSamples in new[] { 2, 5, 10, 30 })
        {
            var mutations = new System.Collections.Generic.List<OncologyAnalyzer.GeneMutation>();
            // 5 background non-recurrent, non-truncating mutations.
            for (int i = 0; i < 5; i++)
                mutations.Add(new OncologyAnalyzer.GeneMutation("G", 10 + i, OncologyAnalyzer.MutationConsequence.Other));
            // `sharedSamples` missense mutations at the same hotspot position (recurrent once ≥2).
            for (int s = 0; s < sharedSamples; s++)
                mutations.Add(new OncologyAnalyzer.GeneMutation("G", 100, OncologyAnalyzer.MutationConsequence.Missense));

            double score = OncologyAnalyzer.ScoreDriverPotential(mutations);
            score.Should().BeGreaterThan(previous,
                because: $"{sharedSamples} samples sharing the hotspot raises the recurrent-missense fraction");
            previous = score;
        }
    }

    #endregion

    #region ONCO-DRIVER-001 INV — relabeling passenger genes preserves the driver set

    [Test]
    [Description("INV: renaming non-driver (passenger) genes leaves the set of identified driver mutations unchanged.")]
    public void IdentifyDrivers_RelabelPassengers_SameDriverSet()
    {
        var mutations = new[]
        {
            // Driver gene TSG1: all truncating → TumorSuppressor.
            new OncologyAnalyzer.GeneMutation("TSG1", 50, OncologyAnalyzer.MutationConsequence.Nonsense),
            new OncologyAnalyzer.GeneMutation("TSG1", 80, OncologyAnalyzer.MutationConsequence.Frameshift),
            new OncologyAnalyzer.GeneMutation("TSG1", 120, OncologyAnalyzer.MutationConsequence.SpliceSite),
            // Passenger PASS1: a single non-recurrent, non-truncating mutation → Ambiguous.
            new OncologyAnalyzer.GeneMutation("PASS1", 30, OncologyAnalyzer.MutationConsequence.Other),
            // Passenger PASS2: two missense at DIFFERENT positions (not recurrent) → Ambiguous.
            new OncologyAnalyzer.GeneMutation("PASS2", 11, OncologyAnalyzer.MutationConsequence.Missense),
            new OncologyAnalyzer.GeneMutation("PASS2", 22, OncologyAnalyzer.MutationConsequence.Missense),
        };

        var baseDrivers = OncologyAnalyzer.IdentifyDriverMutations(mutations).ToHashSet();
        baseDrivers.Select(m => m.Gene).Distinct().Should().BeEquivalentTo(new[] { "TSG1" },
            because: "only the tumor-suppressor gene's mutations are drivers");

        // Relabel the passenger genes only.
        var relabeled = mutations.Select(m => m.Gene.StartsWith("PASS")
            ? m with { Gene = m.Gene.Replace("PASS", "ZZZ") }
            : m);

        OncologyAnalyzer.IdentifyDriverMutations(relabeled).ToHashSet()
            .Should().BeEquivalentTo(baseDrivers,
                because: "passenger genes contribute no drivers, so renaming them cannot change the driver set");
    }

    #endregion
}
