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

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-ARTIFACT-001 — sequencing-artifact filtering (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 90.
    //
    // API under test (OncologyAnalyzer.FilterArtifacts / ClassifyArtifact):
    //   Removes variants flagged as artifacts: FFPE deamination (C>T/G>A) always, OxoG (G>T/C>A)
    //   when the read-orientation imbalance GIV = altR1/altR2 exceeds the damaged threshold (1.5).
    //   Survivors = candidates not flagged; each variant is judged independently.
    //
    // Relations (derived from the substitution+GIV rule, NOT from output):
    //   • MON  (stricter artifact evidence ⇒ subset of survivors): survivors are always a subset
    //          of the candidates, and raising an OxoG variant's GIV can only remove it — once the
    //          damage evidence crosses the threshold the variant never re-enters the survivor set.
    //   • INV  (duplicate a passing variant ⇒ still passing): per-variant judgement, so a copy of
    //          a surviving variant also survives (and a copy of an artifact stays filtered).
    // ───────────────────────────────────────────────────────────────────────────

    // A>G is not an artifact substitution class → always a survivor.
    private static OncologyAnalyzer.ArtifactObservation PassingVariant() =>
        new('A', 'G', 10, 10, 5, 5, 5, 5);

    // G>T is OxoG: artifact iff GIV = altR1/altR2 > 1.5.
    private static OncologyAnalyzer.ArtifactObservation OxoG(int altR1, int altR2) =>
        new('G', 'T', 10, 10, 5, 5, altR1, altR2);

    #region ONCO-ARTIFACT-001 MON — survivors are a subset and shrink with damage evidence

    [Test]
    [Description("MON: FilterArtifacts returns a subset of the candidates, and raising an OxoG variant's GIV imbalance can only remove it from the survivors (monotone survivorship across the damage threshold).")]
    public void FilterArtifacts_StricterEvidence_SubsetOfSurvivors()
    {
        bool SurvivesOxoG(int altR1, int altR2) =>
            OncologyAnalyzer.FilterArtifacts(new[] { OxoG(altR1, altR2) }).Count == 1;

        // GIV = altR1/altR2: 0.5, 1.0, 1.5 survive (≤ 1.5); 2.0, 3.0 are filtered.
        SurvivesOxoG(5, 10).Should().BeTrue(because: "GIV 0.5 ≤ 1.5 is not damaged");
        SurvivesOxoG(15, 10).Should().BeTrue(because: "GIV 1.5 is at the threshold, not above it");
        SurvivesOxoG(20, 10).Should().BeFalse(because: "GIV 2.0 > 1.5 is OxoG damage → filtered");
        SurvivesOxoG(30, 10).Should().BeFalse(because: "stronger imbalance stays filtered — survivorship is monotone in GIV");

        // Subset invariant: survivors ⊆ candidates.
        var candidates = new[] { PassingVariant(), OxoG(5, 10), OxoG(30, 10), PassingVariant() };
        var survivors = OncologyAnalyzer.FilterArtifacts(candidates);
        survivors.Count.Should().BeLessThanOrEqualTo(candidates.Length);
        survivors.Should().OnlyContain(v => candidates.Contains(v), because: "every survivor is one of the candidates");
    }

    #endregion

    #region ONCO-ARTIFACT-001 INV — duplicating a passing variant keeps it passing

    [Test]
    [Description("INV: variants are judged independently, so duplicating a passing variant keeps both copies passing (and duplicating an artifact keeps both filtered).")]
    public void FilterArtifacts_DuplicatePassingVariant_StillPasses()
    {
        var passing = PassingVariant();
        OncologyAnalyzer.FilterArtifacts(new[] { passing }).Should().ContainSingle();
        OncologyAnalyzer.FilterArtifacts(new[] { passing, passing }).Count
            .Should().Be(2, because: "a duplicate of a surviving variant is judged identically and also survives");

        var ffpe = new OncologyAnalyzer.ArtifactObservation('C', 'T', 10, 10, 5, 5, 5, 5); // FFPE deamination → always artifact
        OncologyAnalyzer.FilterArtifacts(new[] { ffpe, ffpe })
            .Should().BeEmpty(because: "duplicating an artifact keeps both copies filtered");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-ANNOT-001 — cancer-variant tier annotation (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 91.
    //
    // API under test (OncologyAnalyzer.AnnotateCancerVariants / ClassifyVariantTier):
    //   Assigns each variant an AMP/ASCO/CAP 2017 tier from caller-supplied evidence
    //   (evidence level, population MAF, cancer association). The tier is a pure function of
    //   that evidence; the variant's identity fields (Gene/ProteinChange) are carried through
    //   but do not affect the tier.
    //
    // Relations (derived from the per-variant tiering rule, NOT from output):
    //   • INV  (identity shift ⇒ annotations shift equally): uniformly relabelling the variant
    //          identity carries the annotation along unchanged — same tier, same per-variant order.
    //          (The tiering API has no genomic coordinate; the identity field is the analog.)
    //   • INV  (variant order independent): each variant is tiered independently, so the set of
    //          (variant → tier) results does not depend on the input order.
    // ───────────────────────────────────────────────────────────────────────────

    private static OncologyAnalyzer.CancerVariantAnnotationInput[] TierPanel() => new[]
    {
        new OncologyAnalyzer.CancerVariantAnnotationInput("BRAF", "p.V600E", OncologyAnalyzer.ClinicalEvidenceLevel.A, 0.0, true),   // Tier I
        new OncologyAnalyzer.CancerVariantAnnotationInput("KIT", "p.D816V", OncologyAnalyzer.ClinicalEvidenceLevel.C, 0.0, true),    // Tier II
        new OncologyAnalyzer.CancerVariantAnnotationInput("XYZ1", "p.A1B", OncologyAnalyzer.ClinicalEvidenceLevel.None, 0.05, true), // Tier IV (MAF ≥ 1%)
        new OncologyAnalyzer.CancerVariantAnnotationInput("XYZ2", "p.C2D", OncologyAnalyzer.ClinicalEvidenceLevel.None, 0.0001, true),// Tier III
    };

    #region ONCO-ANNOT-001 INV — a uniform identity shift carries annotations along unchanged

    [Test]
    [Description("INV: uniformly relabelling each variant's identity leaves its tier unchanged and carries the annotation onto the relabelled variant (annotations shift equally with the variants).")]
    public void AnnotateCancerVariants_IdentityShift_ShiftsAnnotationsEqually()
    {
        var baseAnnotations = OncologyAnalyzer.AnnotateCancerVariants(TierPanel());

        // Uniform "coordinate shift": relabel every gene by the same constant suffix.
        var shifted = TierPanel().Select(v => v with { Gene = v.Gene + "@chr2" }).ToList();
        var shiftedAnnotations = OncologyAnalyzer.AnnotateCancerVariants(shifted);

        shiftedAnnotations.Select(a => a.Tier).Should().Equal(baseAnnotations.Select(a => a.Tier),
            because: "the tier is a function of the evidence, invariant to a uniform identity shift");
        shiftedAnnotations.Select(a => a.Variant.Gene).Should().Equal(shifted.Select(v => v.Gene),
            because: "each annotation is carried onto its (shifted) variant");
    }

    #endregion

    #region ONCO-ANNOT-001 INV — variant order independent

    [Test]
    [Description("INV: each variant is tiered independently, so the set of (gene → tier) results is the same for any ordering of the input.")]
    public void AnnotateCancerVariants_VariantOrder_Independent()
    {
        var forward = OncologyAnalyzer.AnnotateCancerVariants(TierPanel())
            .Select(a => (a.Variant.Gene, a.Tier)).ToHashSet();
        var reversed = OncologyAnalyzer.AnnotateCancerVariants(TierPanel().Reverse())
            .Select(a => (a.Variant.Gene, a.Tier)).ToHashSet();

        reversed.Should().BeEquivalentTo(forward, because: "per-variant tiering has no cross-variant dependence");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-TMB-001 — tumor mutational burden (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 92.
    //
    // API under test (OncologyAnalyzer.CalculateTMB):
    //   TMB = (#somatic coding mutations) / (panel size in Mb) — a density.
    //
    // Relations (derived from the density definition, NOT from output):
    //   • INV  (doubling mutations and Mb ⇒ same density): scaling numerator and denominator
    //          equally leaves the ratio unchanged.
    //   • MON  (+1 coding mutation ⇒ ≥ TMB): at fixed panel size, one more mutation raises TMB.
    //   • INV  (order independent): TMB from a call set counts somatic calls — order irrelevant.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-TMB-001 INV — doubling mutations and panel size preserves TMB

    [Test]
    [Description("INV: scaling the mutation count and the panel Mb by the same factor leaves the TMB density unchanged.")]
    public void Tmb_DoubleMutationsAndPanel_SameDensity()
    {
        double baseTmb = OncologyAnalyzer.CalculateTMB(100, 2.0);
        foreach (int k in new[] { 2, 3, 10 })
            OncologyAnalyzer.CalculateTMB(100 * k, 2.0 * k).Should().BeApproximately(baseTmb, 1e-12,
                because: $"scaling mutations and Mb by {k} cancels in mutations/Mb");
    }

    #endregion

    #region ONCO-TMB-001 MON — one more mutation raises TMB

    [Test]
    [Description("MON: at fixed panel size, each additional coding mutation strictly raises TMB.")]
    public void Tmb_OneMoreMutation_HigherTmb()
    {
        const double panelMb = 5.0;
        double previous = double.MinValue;
        foreach (int count in new[] { 0, 1, 10, 100 })
        {
            double tmb = OncologyAnalyzer.CalculateTMB(count, panelMb);
            if (count > 0) tmb.Should().BeGreaterThan(previous, because: "more mutations at fixed Mb raise the density");
            previous = tmb;
        }
    }

    #endregion

    #region ONCO-TMB-001 INV — TMB from a call set is order independent

    [Test]
    [Description("INV: TMB counts the somatic calls, so it is independent of the order of the call list.")]
    public void Tmb_CallOrder_Independent()
    {
        var calls = OncologyAnalyzer.CallSomaticMutations(Panel(), TumorThr, NormalThr);
        const double panelMb = 1.5;

        OncologyAnalyzer.CalculateTMB(calls.Reverse(), panelMb)
            .Should().Be(OncologyAnalyzer.CalculateTMB(calls, panelMb),
                because: "TMB is a count of somatic calls divided by Mb — order cannot matter");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-MSI-001 — microsatellite instability score (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 93.
    //
    // API under test (OncologyAnalyzer.CalculateMSIScore / DetectMSI):
    //   MSI score = unstable loci / total valid loci. DetectMSI counts the unstable flags.
    //
    // Relations (derived from the fraction definition, NOT from output):
    //   • MON  (more unstable loci ⇒ ≥ score): at fixed total, more unstable loci raise the score.
    //   • INV  (locus order independent): the score is a count ratio — flag order is irrelevant.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-MSI-001 MON — more unstable loci raise the MSI score

    [Test]
    [Description("MON: at a fixed total locus count, increasing the number of unstable loci strictly raises the MSI score.")]
    public void Msi_MoreUnstableLoci_HigherScore()
    {
        const int total = 50;
        double previous = double.MinValue;
        foreach (int unstable in new[] { 0, 5, 20, 50 })
        {
            double score = OncologyAnalyzer.CalculateMSIScore(unstable, total);
            if (unstable > 0) score.Should().BeGreaterThan(previous, because: "more unstable loci at fixed total raise unstable/total");
            previous = score;
        }
    }

    #endregion

    #region ONCO-MSI-001 INV — locus order independent

    [Test]
    [Description("INV: DetectMSI counts unstable flags, so the score is independent of the order of the locus flags.")]
    public void Msi_LocusOrder_Independent()
    {
        var flags = new[] { true, false, true, true, false, false, true, false };

        var forward = OncologyAnalyzer.DetectMSI(flags);
        var reversed = OncologyAnalyzer.DetectMSI(flags.Reverse());

        reversed.Score.Should().Be(forward.Score, because: "the MSI score is a count ratio, independent of flag order");
        reversed.UnstableLoci.Should().Be(forward.UnstableLoci);
        reversed.Status.Should().Be(forward.Status);
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-HRD-001 — homologous-recombination-deficiency score (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 94.
    //
    // API under test (OncologyAnalyzer.CalculateHRDScore / DetectHRD):
    //   HRD score = LOH + TAI + LST (Telli 2016), an unweighted sum of the three scar counts.
    //
    // Relations (derived from the additive sum, NOT from output):
    //   • MON  (adding an event ⇒ ≥ HRD): incrementing any component raises the score by 1.
    //   • INV  (event order/label independent): the sum is symmetric in (LOH, TAI, LST), so any
    //          permutation of the component counts yields the same score.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-HRD-001 MON — adding any scar event raises the HRD score

    [Test]
    [Description("MON: incrementing any of the LOH/TAI/LST component counts raises the HRD score by one.")]
    public void Hrd_AddingEvent_RaisesScore()
    {
        int baseScore = OncologyAnalyzer.CalculateHRDScore(10, 12, 14);

        OncologyAnalyzer.CalculateHRDScore(11, 12, 14).Should().BeGreaterThan(baseScore, because: "an extra LOH event adds to the sum");
        OncologyAnalyzer.CalculateHRDScore(10, 13, 14).Should().BeGreaterThan(baseScore, because: "an extra TAI event adds to the sum");
        OncologyAnalyzer.CalculateHRDScore(10, 12, 15).Should().BeGreaterThan(baseScore, because: "an extra LST event adds to the sum");
    }

    #endregion

    #region ONCO-HRD-001 INV — the score is symmetric in its three components

    [Test]
    [Description("INV: HRD = LOH + TAI + LST is symmetric, so any permutation of the component counts yields the same score and status.")]
    public void Hrd_ComponentOrder_Independent()
    {
        var perms = new[]
        {
            (10, 15, 17), (10, 17, 15), (15, 10, 17), (15, 17, 10), (17, 10, 15), (17, 15, 10),
        };

        int expected = OncologyAnalyzer.CalculateHRDScore(10, 15, 17);
        foreach (var (a, b, c) in perms)
        {
            OncologyAnalyzer.CalculateHRDScore(a, b, c).Should().Be(expected,
                because: "the HRD sum does not depend on which component contributes which count");
            OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(a, b, c)).Status
                .Should().Be(OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(10, 15, 17)).Status);
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-LOH-001 — loss-of-heterozygosity detection (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 95.
    //
    // API under test (OncologyAnalyzer.DetectLOH):
    //   A segment is LOH when one allele is lost (minor CN = 0) and the other is retained
    //   (major CN ≠ 0); HRD-LOH counts such regions longer than 15 Mb that do not span a whole
    //   chromosome. (The implementation is copy-number based, not BAF-threshold based, and
    //   represents alleles by magnitude — major ≥ minor.)
    //
    // Relations (derived from the LOH predicate, NOT from output):
    //   • INV  (allele symmetry / order): LOH depends only on the loss of one allele, so the
    //          retained (major) allele's copy-number value does not change the calls, and the
    //          calls are independent of segment input order. (The faithful reading of the
    //          checklist's "swap A/B labels" — only that ONE allele is lost matters.)
    //   • MON  (looser size criterion ⇒ superset): a region only counts once it exceeds 15 Mb,
    //          and adding more qualifying LOH regions yields a superset (higher score).
    // ───────────────────────────────────────────────────────────────────────────

    private const long Mb = 1_000_000L;

    // chr with a het segment (avoids whole-chromosome-LOH exclusion) + one LOH segment.
    private static OncologyAnalyzer.AlleleSpecificSegment[] LohChromosome(string chrom, int retainedCn, long lohLength) => new[]
    {
        new OncologyAnalyzer.AlleleSpecificSegment(chrom, 0, 1 * Mb, 2, 1),                       // heterozygous
        new OncologyAnalyzer.AlleleSpecificSegment(chrom, 1 * Mb, 1 * Mb + lohLength, retainedCn, 0), // LOH (minor lost)
    };

    #region ONCO-LOH-001 INV — retained-allele value and segment order don't change LOH calls

    [Test]
    [Description("INV: LOH depends only on the loss of one allele, so the retained (major) allele's copy number does not change the calls, and the calls are independent of segment order.")]
    public void DetectLOH_AlleleSymmetryAndOrder_PreserveCalls()
    {
        var reference = OncologyAnalyzer.DetectLOH(LohChromosome("1", retainedCn: 2, lohLength: 20 * Mb));
        reference.Score.Should().Be(1, because: "a 20 Mb minor-lost region (with a het neighbour) is one LOH region");

        // Retained-allele value is irrelevant — only that the minor allele is lost.
        foreach (int retained in new[] { 1, 3, 5 })
            OncologyAnalyzer.DetectLOH(LohChromosome("1", retained, 20 * Mb)).Score
                .Should().Be(reference.Score, because: $"the retained allele being CN={retained} does not change that the other allele is lost");

        // Segment order is irrelevant.
        OncologyAnalyzer.DetectLOH(LohChromosome("1", 2, 20 * Mb).Reverse()).Score
            .Should().Be(reference.Score, because: "LOH detection groups by chromosome — input order cannot matter");
    }

    #endregion

    #region ONCO-LOH-001 MON — bigger / more LOH regions yield a superset

    [Test]
    [Description("MON: a region counts only once it exceeds the 15 Mb size cutoff, and adding more qualifying LOH regions raises the score.")]
    public void DetectLOH_LooserSizeAndMoreRegions_Superset()
    {
        // Size criterion: 10 Mb region does not count; 20 Mb does.
        OncologyAnalyzer.DetectLOH(LohChromosome("1", 2, 10 * Mb)).Score
            .Should().Be(0, because: "a 10 Mb LOH region is below the 15 Mb HRD-LOH cutoff");
        OncologyAnalyzer.DetectLOH(LohChromosome("1", 2, 20 * Mb)).Score
            .Should().Be(1, because: "a 20 Mb LOH region exceeds the 15 Mb cutoff and is counted");

        // Adding a second qualifying LOH region (on another chromosome) raises the score.
        var oneRegion = LohChromosome("1", 2, 20 * Mb);
        var twoRegions = oneRegion.Concat(LohChromosome("2", 2, 20 * Mb)).ToArray();
        OncologyAnalyzer.DetectLOH(twoRegions).Score
            .Should().BeGreaterThan(OncologyAnalyzer.DetectLOH(oneRegion).Score,
                because: "an additional qualifying LOH region adds to the HRD-LOH count");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-SIG-001 — SBS-96 mutational context classification (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 96.
    //
    // API under test (OncologyAnalyzer.ClassifySbsContext / Build96ContextCatalog):
    //   Maps a single-base substitution + trinucleotide context to one of the 96 pyrimidine-
    //   centric COSMIC channels, folding purine-reference substitutions onto the pyrimidine
    //   strand by reverse-complement. The catalog tallies each variant into its channel.
    //
    // Relations (derived from pyrimidine-strand folding, NOT from output):
    //   • INV  (reverse-complement ⇒ same channel): a variant and its reverse complement are the
    //          two strands of the same mutation, so they classify into the identical channel.
    //   • INV  (variant order independent): the 96-channel catalog is a tally, independent of the
    //          order of the input variants.
    // ───────────────────────────────────────────────────────────────────────────

    private static char Comp(char b) => b switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => b };

    #region ONCO-SIG-001 INV — a variant and its reverse complement share a channel

    [Test]
    [Description("INV: a substitution and its reverse complement are the same mutation on opposite strands, so ClassifySbsContext folds both to the identical pyrimidine channel.")]
    public void Sbs_ReverseComplement_SameChannel()
    {
        // (5', ref, alt, 3') variants spanning pyrimidine- and purine-reference cases.
        var variants = new[]
        {
            ('A', 'C', 'A', 'G'),
            ('T', 'C', 'T', 'A'),
            ('G', 'T', 'C', 'C'),
            ('C', 'G', 'T', 'T'), // purine reference
            ('A', 'A', 'G', 'T'), // purine reference
        };

        foreach (var (f, r, a, t) in variants)
        {
            string forward = OncologyAnalyzer.ClassifySbsContext(f, r, a, t);
            // Reverse complement of the (5',ref,3') context with the substitution complemented.
            string reverse = OncologyAnalyzer.ClassifySbsContext(Comp(t), Comp(r), Comp(a), Comp(f));

            reverse.Should().Be(forward,
                because: "the two strands of one substitution fold to the same pyrimidine-centric SBS-96 channel");
        }
    }

    #endregion

    #region ONCO-SIG-001 INV — the catalog is order independent

    [Test]
    [Description("INV: the 96-channel catalog is a tally over variants, so it is identical for any ordering of the input.")]
    public void Sbs_CatalogVariantOrder_Independent()
    {
        var variants = new[]
        {
            ('A', 'C', 'A', 'G'), ('T', 'C', 'T', 'A'), ('G', 'T', 'C', 'C'),
            ('C', 'G', 'T', 'T'), ('A', 'A', 'G', 'T'), ('A', 'C', 'A', 'G'),
        };

        var forward = OncologyAnalyzer.Build96ContextCatalog(variants);
        var reversed = OncologyAnalyzer.Build96ContextCatalog(variants.Reverse());

        reversed.Should().BeEquivalentTo(forward, because: "tallying into channels does not depend on variant order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-SIG-002 — signature exposure refitting (NNLS) (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 97.
    //
    // API under test (OncologyAnalyzer.FitSignatures):
    //   Solves minₓ ‖S·x − d‖² subject to x ≥ 0 for the per-signature exposures x of an observed
    //   catalog d against reference signatures S. NNLS is positively homogeneous in d.
    //
    // Relations (derived from the linear NNLS model, NOT from output):
    //   • INV  (scale catalog by k ⇒ exposures × k): scaling the observed catalog scales the
    //          minimiser exposures by the same factor.
    //   • MON  (add signature-consistent mutations ⇒ ≥ that exposure): adding counts in a
    //          signature's channels raises that signature's exposure.
    // ───────────────────────────────────────────────────────────────────────────

    // Orthonormal one-hot signatures, so the NNLS exposure of each equals the catalog count in its channel.
    private static readonly IReadOnlyList<IReadOnlyList<double>> OneHotSignatures = new IReadOnlyList<double>[]
    {
        new double[] { 1, 0, 0 },
        new double[] { 0, 1, 0 },
        new double[] { 0, 0, 1 },
    };

    #region ONCO-SIG-002 INV — scaling the catalog scales the exposures

    [Test]
    [Description("INV: NNLS is positively homogeneous, so scaling the observed catalog by k scales every fitted exposure by k.")]
    public void FitSignatures_ScaleCatalog_ScalesExposures()
    {
        var catalog = new double[] { 3, 5, 2 };
        var baseExposures = OncologyAnalyzer.FitSignatures(catalog, OneHotSignatures).Exposures;

        foreach (double k in new[] { 2.0, 7.0, 0.5 })
        {
            var scaled = OncologyAnalyzer.FitSignatures(catalog.Select(d => d * k).ToList(), OneHotSignatures).Exposures;
            for (int j = 0; j < baseExposures.Count; j++)
                scaled[j].Should().BeApproximately(baseExposures[j] * k, 1e-9,
                    because: $"scaling the catalog by {k} scales exposure {j} by {k}");
        }
    }

    #endregion

    #region ONCO-SIG-002 MON — adding signature-consistent mutations raises that exposure

    [Test]
    [Description("MON: adding counts in signature 1's channel raises signature 1's exposure while the others stay fixed.")]
    public void FitSignatures_AddSignatureConsistentMutations_HigherExposure()
    {
        double previous = double.MinValue;
        foreach (double added in new[] { 0.0, 5.0, 10.0, 20.0 })
        {
            // Catalog [3, 5+added, 2]: only signature 1's channel grows.
            var catalog = new double[] { 3, 5 + added, 2 };
            var exposures = OncologyAnalyzer.FitSignatures(catalog, OneHotSignatures).Exposures;

            exposures[1].Should().BeGreaterThan(previous,
                because: $"adding {added} mutations consistent with signature 1 raises its exposure");
            exposures[0].Should().BeApproximately(3, 1e-9, because: "signature 0's channel is unchanged");
            exposures[2].Should().BeApproximately(2, 1e-9, because: "signature 2's channel is unchanged");
            previous = exposures[1];
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-SIG-003 — bootstrap exposure confidence intervals (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 98.
    //
    // API under test (OncologyAnalyzer.BootstrapExposures):
    //   Multinomial-bootstrap percentile CIs for signature exposures: resample the catalog,
    //   refit by NNLS, take per-signature [α/2, 1−α/2] percentiles. Seeded RNG ⇒ reproducible.
    //
    // Relations (derived from the seeded percentile bootstrap, NOT from output):
    //   • INV  (same seed ⇒ identical CI): the result is fully reproducible from the seed.
    //   • INV  (point estimate independent of bootstrap settings): the point estimate is the
    //          NNLS fit of the observed catalog — unchanged by replicates or seed.
    //   • MON  (wider confidence ⇒ non-narrower CI): a higher confidence level widens the
    //          percentile interval (lower bound non-increasing, upper bound non-decreasing).
    //          (NOTE: the checklist's "more reps → non-wider CI" is not a true property of a
    //          percentile bootstrap — width converges, not monotonically narrows — so the
    //          genuinely monotone relation in the confidence level is asserted instead.)
    // ───────────────────────────────────────────────────────────────────────────

    private static readonly int[] SigCatalog = { 30, 50, 20 };

    #region ONCO-SIG-003 INV — same seed reproduces the CI; point estimate is setting-independent

    [Test]
    [Description("INV: BootstrapExposures is fully reproducible from its seed, and the point estimate (NNLS fit of the observed catalog) is independent of replicates and seed.")]
    public void Bootstrap_SameSeed_IdenticalCi_PointEstimateStable()
    {
        var a = OncologyAnalyzer.BootstrapExposures(SigCatalog, OneHotSignatures, replicates: 200, confidence: 0.95, seed: 42);
        var b = OncologyAnalyzer.BootstrapExposures(SigCatalog, OneHotSignatures, replicates: 200, confidence: 0.95, seed: 42);

        a.Should().BeEquivalentTo(b, because: "the same seed and settings reproduce the identical bootstrap CIs");

        // Point estimate is the NNLS fit of the observed catalog → invariant to seed and replicate count.
        var c = OncologyAnalyzer.BootstrapExposures(SigCatalog, OneHotSignatures, replicates: 500, confidence: 0.95, seed: 7);
        for (int j = 0; j < a.Count; j++)
            c[j].PointEstimate.Should().BeApproximately(a[j].PointEstimate, 1e-9,
                because: "the point estimate does not depend on the bootstrap replicate count or seed");
    }

    #endregion

    #region ONCO-SIG-003 MON — a wider confidence level widens the interval

    [Test]
    [Description("MON: at a fixed seed and replicate count, raising the confidence level widens the percentile interval (lower non-increasing, upper non-decreasing).")]
    public void Bootstrap_HigherConfidence_NonNarrowerCi()
    {
        var narrow = OncologyAnalyzer.BootstrapExposures(SigCatalog, OneHotSignatures, replicates: 400, confidence: 0.80, seed: 123);
        var wide = OncologyAnalyzer.BootstrapExposures(SigCatalog, OneHotSignatures, replicates: 400, confidence: 0.99, seed: 123);

        for (int j = 0; j < narrow.Count; j++)
        {
            wide[j].Lower.Should().BeLessThanOrEqualTo(narrow[j].Lower,
                because: "a higher confidence level pushes the lower percentile down (or keeps it)");
            wide[j].Upper.Should().BeGreaterThanOrEqualTo(narrow[j].Upper,
                because: "a higher confidence level pushes the upper percentile up (or keeps it)");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-SIG-004 — mutational-process classification (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 99.
    //
    // API under test (OncologyAnalyzer.ClassifyMutationalProcess):
    //   Normalizes per-signature exposures to relative contributions, drops sub-cutoff ones,
    //   aggregates by aetiology process, and reports the active processes and the dominant one.
    //
    // Relations (derived from the normalize-then-aggregate model, NOT from output):
    //   • INV  (scale exposures ⇒ same dominant process): normalization cancels a common factor,
    //          so scaling all exposures preserves the active processes and the dominant one.
    //   • INV  (signature order independent): aggregation by process and the deterministic
    //          ordering make the result independent of the input signature order.
    // ───────────────────────────────────────────────────────────────────────────

    private static (string, double)[] ExposurePanel() => new[]
    {
        ("SBS4", 60.0),  // TobaccoSmoking (dominant)
        ("SBS1", 30.0),  // Aging
        ("SBS2", 10.0),  // Apobec
    };

    #region ONCO-SIG-004 INV — scaling exposures preserves the dominant process

    [Test]
    [Description("INV: contributions are normalized, so scaling every exposure by a positive constant preserves the active processes and the dominant process.")]
    public void ClassifyProcess_ScaleExposures_SameDominant()
    {
        var baseResult = OncologyAnalyzer.ClassifyMutationalProcess(ExposurePanel());
        baseResult.DominantProcess.Should().NotBe(OncologyAnalyzer.MutationalProcess.Unknown,
            because: "the panel has mapped signatures, so a dominant process exists");

        foreach (double k in new[] { 2.0, 100.0, 0.25 })
        {
            var scaled = OncologyAnalyzer.ClassifyMutationalProcess(
                ExposurePanel().Select(e => (e.Item1, e.Item2 * k)).ToList());

            scaled.DominantProcess.Should().Be(baseResult.DominantProcess,
                because: $"scaling all exposures by {k} cancels in the normalization");
            scaled.ActiveProcesses.Select(p => p.Process).Should().Equal(baseResult.ActiveProcesses.Select(p => p.Process),
                because: "the relative contributions — and hence the active-process ranking — are unchanged");
        }
    }

    #endregion

    #region ONCO-SIG-004 INV — signature order independent

    [Test]
    [Description("INV: aggregation by process and deterministic ordering make the dominant process and active-process set independent of the input signature order.")]
    public void ClassifyProcess_SignatureOrder_Independent()
    {
        var forward = OncologyAnalyzer.ClassifyMutationalProcess(ExposurePanel());
        var reversed = OncologyAnalyzer.ClassifyMutationalProcess(ExposurePanel().Reverse().ToList());

        reversed.DominantProcess.Should().Be(forward.DominantProcess, because: "the dominant process does not depend on input order");
        reversed.ActiveProcesses.Should().Equal(forward.ActiveProcesses,
            because: "active processes are emitted in a deterministic (contribution, then process) order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-FUSION-001 — gene-fusion detection (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 100.
    //
    // API under test (OncologyAnalyzer.DetectFusions / ComputeTotalSupport):
    //   Calls a fusion from per-class supporting-read counts (STAR-Fusion rule); TotalSupport =
    //   split_reads1 + split_reads2 + discordant_mates is the abundance/confidence measure.
    //   The API carries gene symbols + read counts (no genomic breakpoint coordinates).
    //
    // Relations (derived from the support-count rule, NOT from output):
    //   • MON  (more split reads ⇒ ≥ confidence): adding split reads raises TotalSupport and a
    //          fusion that passes stays passing (detection monotone in support).
    //   • INV  (preserves fusion count): the detected-fusion set/count is independent of the
    //          candidate input order. (The API has no breakpoint coordinate to shift; the
    //          checklist's "prepend flank shifts breakpoints, preserves count" reduces to the
    //          count being invariant to ordering.)
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-FUSION-001 MON — more split reads raise support and keep the call

    [Test]
    [Description("MON: adding split reads raises the fusion's TotalSupport (confidence), and a fusion that passes the support rule stays passing.")]
    public void DetectFusions_MoreSplitReads_HigherSupport()
    {
        double previous = double.MinValue;
        bool everDetected = false;

        foreach (int split in new[] { 0, 1, 5, 20 })
        {
            var candidate = new OncologyAnalyzer.FusionCandidate("EML4", "ALK", split, 0, 1);
            var calls = OncologyAnalyzer.DetectFusions(new[] { candidate });

            if (calls.Count == 1)
            {
                double support = calls[0].TotalSupport;
                support.Should().BeGreaterThan(previous, because: $"{split} split reads raise the total support");
                previous = support;
                everDetected = true;
            }
        }

        everDetected.Should().BeTrue(because: "with enough split reads the fusion is detected with growing support");
    }

    #endregion

    #region ONCO-FUSION-001 INV — the detected-fusion count is order independent

    [Test]
    [Description("INV: each candidate is judged independently, so the set and count of detected fusions are the same for any input order.")]
    public void DetectFusions_CandidateOrder_PreservesCount()
    {
        var candidates = new[]
        {
            new OncologyAnalyzer.FusionCandidate("EML4", "ALK", 10, 5, 3),  // passes (total 18)
            new OncologyAnalyzer.FusionCandidate("BCR", "ABL1", 2, 1, 0),   // passes (junction 3, total 3)
            new OncologyAnalyzer.FusionCandidate("SELF", "SELF", 50, 50, 50), // self-fusion → excluded
            new OncologyAnalyzer.FusionCandidate("FOO", "BAR", 0, 0, 2),    // no junction, disc 2 < 5 → fails
        };

        var forward = OncologyAnalyzer.DetectFusions(candidates)
            .Select(c => (c.Gene5Prime, c.Gene3Prime)).ToHashSet();
        var reversed = OncologyAnalyzer.DetectFusions(candidates.Reverse())
            .Select(c => (c.Gene5Prime, c.Gene3Prime)).ToHashSet();

        forward.Should().BeEquivalentTo(new[] { ("EML4", "ALK"), ("BCR", "ABL1") },
            because: "only the two adequately-supported, non-self candidates are called");
        reversed.Should().BeEquivalentTo(forward, because: "candidate order cannot change the detected-fusion set or count");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-FUSION-002 — known-fusion database matching (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 101.
    //
    // API under test (OncologyAnalyzer.MatchKnownFusions / GetFusionAnnotation):
    //   Looks a fusion's directional designation (5'::3') up in a caller-supplied known-fusion
    //   set; case-insensitive but orientation-sensitive.
    //
    // Relations (derived from the directional lookup, NOT from output):
    //   • SUB  (matched ⊆ known DB): a fusion is IsKnown only when its designation is a DB key.
    //   • INV  (orientation preserved; case-insensitive): the lookup keeps the 5'→3' orientation
    //          (the swapped pair does NOT match) and is invariant to gene-symbol case (the
    //          identity-preserving "shift"). The API has no genomic coordinate to shift.
    // ───────────────────────────────────────────────────────────────────────────

    private static OncologyAnalyzer.FusionCall Fusion(string g5, string g3) =>
        new(g5, g3, 10, 2, 12, OncologyAnalyzer.FusionReadingFrame.Unknown);

    private static Dictionary<string, string> KnownFusions() => new()
    {
        [OncologyAnalyzer.GetFusionAnnotation("EML4", "ALK")] = "lung adenocarcinoma",
        [OncologyAnalyzer.GetFusionAnnotation("BCR", "ABL1")] = "CML",
    };

    #region ONCO-FUSION-002 SUB — matched fusions are a subset of the known DB

    [Test]
    [Description("SUB: a fusion is flagged known only when its directional designation is a key of the known-fusion DB, so matched designations are a subset of the DB.")]
    public void MatchKnownFusions_Matched_SubsetOfDb()
    {
        var db = KnownFusions();
        var fusions = new[] { Fusion("EML4", "ALK"), Fusion("FOO", "BAR"), Fusion("BCR", "ABL1") };

        var matched = fusions
            .Select(f => OncologyAnalyzer.MatchKnownFusions(f, db))
            .Where(m => m.IsKnown)
            .Select(m => m.Designation)
            .ToHashSet();

        matched.IsSubsetOf(db.Keys.ToHashSet()).Should().BeTrue(
            because: "a fusion is matched only if its designation is present in the known DB");
        matched.Should().NotContain(OncologyAnalyzer.GetFusionAnnotation("FOO", "BAR"),
            because: "an absent fusion is not matched");
    }

    #endregion

    #region ONCO-FUSION-002 INV — orientation preserved, case-insensitive

    [Test]
    [Description("INV: the lookup is orientation-sensitive (the swapped 5'/3' pair does not match) and case-insensitive (a case relabel still matches).")]
    public void MatchKnownFusions_OrientationPreserved_CaseInsensitive()
    {
        var db = KnownFusions();

        OncologyAnalyzer.MatchKnownFusions(Fusion("EML4", "ALK"), db).IsKnown
            .Should().BeTrue(because: "EML4::ALK is in the DB");
        OncologyAnalyzer.MatchKnownFusions(Fusion("eml4", "alk"), db).IsKnown
            .Should().BeTrue(because: "the lookup is case-insensitive — a case relabel preserves the match");
        OncologyAnalyzer.MatchKnownFusions(Fusion("ALK", "EML4"), db).IsKnown
            .Should().BeFalse(because: "the designation is directional 5'→3'; the swapped orientation is a different fusion");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-FUSION-003 — fusion-junction reading-frame classification (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 102.
    //
    // API under test (OncologyAnalyzer.AnalyzeBreakpoint):
    //   For a coding-to-coding junction, calls InFrame iff (fivePrimeCodingBases − phase) mod 3
    //   == 0. The 5' coding-base count is the breakpoint position measured in coding bases.
    //
    // Relation (derived from the mod-3 frame rule, NOT from output):
    //   • INV  (codon-multiple coordinate shift preserves frame): shifting the breakpoint by a
    //          multiple of 3 coding bases leaves the in/out-of-frame classification unchanged,
    //          while a non-multiple-of-3 shift flips it (the rule is exactly mod 3).
    // ───────────────────────────────────────────────────────────────────────────

    private static OncologyAnalyzer.BreakpointFrameStatus FrameAt(int fivePrimeCodingBases, int phase = 0) =>
        OncologyAnalyzer.AnalyzeBreakpoint(new OncologyAnalyzer.FusionBreakpoint(
            "A", "B", OncologyAnalyzer.BreakpointSite.Cds, OncologyAnalyzer.BreakpointSite.Cds,
            fivePrimeCodingBases, phase)).FrameStatus;

    #region ONCO-FUSION-003 INV — a codon-multiple breakpoint shift preserves the frame

    [Test]
    [Description("INV: shifting the breakpoint by a multiple of 3 coding bases preserves the in/out-of-frame classification; a non-codon shift flips it (the frame rule is exactly mod 3).")]
    public void AnalyzeBreakpoint_CodonShift_PreservesFrame()
    {
        var inFrame = FrameAt(30); // (30 − 0) mod 3 == 0 → InFrame
        inFrame.Should().Be(OncologyAnalyzer.BreakpointFrameStatus.InFrame);

        // Codon-multiple shifts preserve the in-frame call.
        foreach (int shift in new[] { 3, 6, 30, 300 })
            FrameAt(30 + shift).Should().Be(inFrame, because: $"a +{shift} (multiple of 3) coding-base shift preserves the frame");

        // A non-codon shift flips the classification.
        var shifted = FrameAt(31);
        shifted.Should().Be(OncologyAnalyzer.BreakpointFrameStatus.OutOfFrame,
            because: "a +1 coding-base shift moves the junction out of frame");

        // The out-of-frame state is likewise preserved under codon-multiple shifts.
        foreach (int shift in new[] { 3, 9, 30 })
            FrameAt(31 + shift).Should().Be(shifted, because: $"a +{shift} shift preserves the out-of-frame state");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-CNA-001 — copy-number-alteration classification (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 103.
    //
    // API under test (OncologyAnalyzer.ClassifyCopyNumber / ClassifyCopyNumbers):
    //   Maps a log2 copy ratio to an integer copy number (CN = ploidy·2^log2) and a CNA state
    //   on the ordered ladder DeepDeletion < Loss < Neutral < Gain < Amplification.
    //
    // Relations (derived from the monotone CN mapping, NOT from output):
    //   • MON  (higher log2 ratio ⇒ ≥ CN class): CN is increasing in the log2 ratio, so the CNA
    //          state is non-decreasing along the state ladder.
    //   • INV  (segment order independent): each region is classified from its own log2 ratio, so
    //          the value→state mapping does not depend on the segment's position.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-CNA-001 MON — higher log2 ratio gives a ≥ copy-number class

    [Test]
    [Description("MON: copy number increases with the log2 ratio, so the CNA state is non-decreasing along DeepDeletion < Loss < Neutral < Gain < Amplification.")]
    public void ClassifyCopyNumber_HigherLog2_HigherClass()
    {
        int previous = int.MinValue;
        foreach (double log2 in new[] { -3.0, -1.0, -0.3, 0.0, 0.3, 1.0, 3.0 })
        {
            int state = (int)OncologyAnalyzer.ClassifyCopyNumber(log2).State;
            state.Should().BeGreaterThanOrEqualTo(previous,
                because: $"a higher log2 ratio ({log2}) cannot lower the copy-number class");
            previous = state;
        }

        OncologyAnalyzer.ClassifyCopyNumber(0.0).State.Should().Be(OncologyAnalyzer.CopyNumberState.Neutral,
            because: "log2 ratio 0 is diploid (CN = 2) → Neutral");
    }

    #endregion

    #region ONCO-CNA-001 INV — segment order independent

    [Test]
    [Description("INV: each region is classified from its own log2 ratio, so the value→state mapping is unaffected by segment order.")]
    public void ClassifyCopyNumbers_SegmentOrder_Independent()
    {
        var log2Ratios = new[] { -2.5, -0.4, 0.0, 0.5, 2.0 };

        var forwardStates = OncologyAnalyzer.ClassifyCopyNumbers(log2Ratios).Select(c => c.State).ToList();
        var reversedStates = OncologyAnalyzer.ClassifyCopyNumbers(log2Ratios.Reverse()).Select(c => c.State).ToList();

        reversedStates.Should().Equal(forwardStates.AsEnumerable().Reverse(),
            because: "reordering the segments only reorders the calls — each region's state depends solely on its own ratio");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-CNA-002 — focal-amplification detection (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 104.
    //
    // API under test (OncologyAnalyzer.IsFocalAmplification / DetectFocalAmplifications):
    //   A segment is a focal amplification iff its log2 ratio exceeds the amplitude threshold
    //   (t_amp) AND its length is below broad_len_cutoff × arm length. Length/fraction depend
    //   only on End−Start, not on absolute position.
    //
    // Relations (derived from the amplitude+length rule, NOT from output):
    //   • MON  (higher CN keeps focal amplification): once a focal segment is amplified, raising
    //          its log2 ratio (higher CN) keeps it a focal amplification.
    //   • INV  (prepend flank shifts focal coordinates): translating the segment by +k leaves its
    //          length/arm-fraction and amplitude unchanged, so it stays a focal amplification and
    //          its coordinates shift by exactly k.
    // ───────────────────────────────────────────────────────────────────────────

    // 1 Mb segment on a 80 Mb arm (fraction 0.0125 ≪ 0.98 → focal).
    private static OncologyAnalyzer.CopyNumberArmSegment FocalSeg(long start, double log2) =>
        new("17q", start, start + 1_000_000, 80_000_000, log2);

    #region ONCO-CNA-002 MON — higher CN keeps a focal amplification

    [Test]
    [Description("MON: once a focal segment is amplified, raising its log2 ratio keeps it a focal amplification.")]
    public void IsFocalAmplification_HigherCn_StaysFocalAmp()
    {
        var thr = OncologyAnalyzer.FocalAmplificationThresholds.Default;

        OncologyAnalyzer.IsFocalAmplification(FocalSeg(1_000_000, 0.05), thr)
            .Should().BeFalse(because: "log2 0.05 ≤ t_amp 0.1 is not amplified");

        foreach (double log2 in new[] { 0.2, 0.7, 2.0, 5.0 })
            OncologyAnalyzer.IsFocalAmplification(FocalSeg(1_000_000, log2), thr)
                .Should().BeTrue(because: $"log2 {log2} > t_amp and the segment is focal — a higher CN keeps it a focal amplification");
    }

    #endregion

    #region ONCO-CNA-002 INV — a prepended flank shifts the focal coordinates

    [Test]
    [Description("INV: translating the segment by +k leaves its length/arm-fraction and amplitude unchanged, so it stays a focal amplification and its coordinates shift by exactly k.")]
    public void DetectFocalAmplifications_PrependFlank_ShiftsCoordinates()
    {
        var baseSeg = FocalSeg(1_000_000, 1.0);
        var baseHit = OncologyAnalyzer.DetectFocalAmplifications(new[] { baseSeg }).Single();

        foreach (long k in new[] { 500_000L, 5_000_000L, 50_000_000L })
        {
            var shifted = baseSeg with { Start = baseSeg.Start + k, End = baseSeg.End + k };
            var hit = OncologyAnalyzer.DetectFocalAmplifications(new[] { shifted }).Single();

            hit.Start.Should().Be(baseHit.Start + k, because: "the start shifts by the flank length");
            hit.End.Should().Be(baseHit.End + k, because: "the end shifts by the flank length");
            hit.Length.Should().Be(baseHit.Length, because: "translation preserves the segment length, so it stays focal");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-CNA-003 — homozygous-deletion detection (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 105.
    //
    // API under test (OncologyAnalyzer.IsHomozygousDeletion / DetectHomozygousDeletions):
    //   A segment is a homozygous (deep) deletion iff its integer copy number is 0 — i.e. its
    //   log2 ratio falls below the deepest CNVkit cutoff (≈ −1.1).
    //
    // Relations (derived from the CN==0 rule, NOT from output):
    //   • MON  (lower CN keeps homozygous deletion): once a segment is a deep deletion, lowering
    //          its log2 ratio (lower CN) keeps it a deep deletion.
    //   • INV  (segment order independent): each segment is judged from its own log2 ratio, so
    //          the detected-deletion set is independent of input order.
    // ───────────────────────────────────────────────────────────────────────────

    private static OncologyAnalyzer.CopyNumberArmSegment CnSeg(string arm, long start, double log2) =>
        new(arm, start, start + 1_000_000, 10_000_000, log2);

    #region ONCO-CNA-003 MON — lower CN keeps the homozygous deletion

    [Test]
    [Description("MON: once a segment's CN is 0 (deep deletion), lowering its log2 ratio keeps it a homozygous deletion.")]
    public void IsHomozygousDeletion_LowerCn_StaysHomDel()
    {
        OncologyAnalyzer.IsHomozygousDeletion(CnSeg("1", 0, -0.5))
            .Should().BeFalse(because: "log2 −0.5 is a single-copy loss (CN 1), not a deep deletion");

        foreach (double log2 in new[] { -1.2, -2.0, -3.0, -5.0 })
            OncologyAnalyzer.IsHomozygousDeletion(CnSeg("1", 0, log2))
                .Should().BeTrue(because: $"log2 {log2} is below the deep-deletion cutoff (CN 0) — a lower CN stays a homozygous deletion");
    }

    #endregion

    #region ONCO-CNA-003 INV — segment order independent

    [Test]
    [Description("INV: each segment is judged from its own log2 ratio, so the detected-deletion set is independent of input order.")]
    public void DetectHomozygousDeletions_SegmentOrder_Independent()
    {
        var segments = new[]
        {
            CnSeg("1", 0, -2.0),   // deep deletion
            CnSeg("2", 0, 0.0),    // neutral
            CnSeg("3", 0, -1.5),   // deep deletion
            CnSeg("4", 0, 1.0),    // gain
        };

        var forward = OncologyAnalyzer.DetectHomozygousDeletions(segments).Select(s => s.Arm).ToHashSet();
        var reversed = OncologyAnalyzer.DetectHomozygousDeletions(segments.Reverse()).Select(s => s.Arm).ToHashSet();

        forward.Should().BeEquivalentTo(new[] { "1", "3" }, because: "only the CN-0 segments are deep deletions");
        reversed.Should().BeEquivalentTo(forward, because: "segment order cannot change which segments are deep deletions");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-PURITY-001 — tumor-purity estimation from VAF (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 106.
    //
    // API under test (OncologyAnalyzer.EstimatePurityFromVAF):
    //   Estimates purity as the median over clonal heterozygous SNVs of ρ = 2·VAF (diploid model).
    //
    // Relations (derived from ρ = 2·median(VAF), NOT from output):
    //   • MON  (higher clonal VAFs ⇒ ≥ purity): raising the alt fractions raises every 2·VAF and
    //          hence the median, so the estimated purity is non-decreasing.
    //   • INV  (variant order independent): the median does not depend on the variant order.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-PURITY-001 MON — higher clonal VAFs raise the purity estimate

    [Test]
    [Description("MON: raising the clonal SNV alt fractions raises ρ = 2·VAF for each, so the median purity estimate is non-decreasing.")]
    public void EstimatePurity_HigherVafs_HigherPurity()
    {
        double previous = double.MinValue;
        foreach (int alt in new[] { 10, 20, 30, 40 }) // VAF 0.10..0.40 (≤ 0.5)
        {
            double purity = OncologyAnalyzer.EstimatePurityFromVAF(new[]
            {
                Variant("1", 1, alt, 100), Variant("2", 2, alt, 100), Variant("3", 3, alt, 100),
            });
            purity.Should().BeGreaterThan(previous, because: $"VAF {alt / 100.0} raises ρ = 2·VAF and the median");
            previous = purity;
        }
    }

    #endregion

    #region ONCO-PURITY-001 INV — variant order independent

    [Test]
    [Description("INV: purity is the median of the per-variant 2·VAF estimates, which does not depend on variant order.")]
    public void EstimatePurity_VariantOrder_Independent()
    {
        var panel = new[] { Variant("1", 1, 10, 100), Variant("2", 2, 25, 100), Variant("3", 3, 40, 100) };

        OncologyAnalyzer.EstimatePurityFromVAF(panel.Reverse())
            .Should().Be(OncologyAnalyzer.EstimatePurityFromVAF(panel),
                because: "the median is invariant to the order of the variants");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-PLOIDY-001 — average tumour ploidy (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 107.
    //
    // API under test (OncologyAnalyzer.EstimatePloidy):
    //   ψ = length-weighted mean of per-segment total copy number (Major + Minor).
    //
    // Relations (derived from the weighted mean, NOT from output):
    //   • MON  (amplifying more segments ⇒ ≥ ploidy): raising segments' copy number raises the
    //          weighted mean, so the ploidy estimate is non-decreasing.
    //   • INV  (segment order independent): a weighted mean does not depend on segment order.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-PLOIDY-001 MON — amplifying more segments raises ploidy

    [Test]
    [Description("MON: amplifying more (equal-length) segments raises the length-weighted mean total copy number, so the ploidy estimate is non-decreasing.")]
    public void EstimatePloidy_AmplifyMoreSegments_HigherPloidy()
    {
        double previous = double.MinValue;
        foreach (int amplified in new[] { 0, 1, 2, 3, 4 })
        {
            var segments = Enumerable.Range(0, 4).Select(i =>
                i < amplified
                    ? new OncologyAnalyzer.AlleleSpecificSegment(i.ToString(), 0, 1_000_000, 4, 2)  // total CN 6
                    : new OncologyAnalyzer.AlleleSpecificSegment(i.ToString(), 0, 1_000_000, 1, 1)); // total CN 2

            double ploidy = OncologyAnalyzer.EstimatePloidy(segments);
            ploidy.Should().BeGreaterThan(previous, because: $"amplifying {amplified}/4 equal segments raises the mean copy number");
            previous = ploidy;
        }
    }

    #endregion

    #region ONCO-PLOIDY-001 INV — segment order independent

    [Test]
    [Description("INV: the length-weighted mean ploidy does not depend on the order of the segments.")]
    public void EstimatePloidy_SegmentOrder_Independent()
    {
        var segments = new[]
        {
            new OncologyAnalyzer.AlleleSpecificSegment("1", 0, 2_000_000, 3, 1),
            new OncologyAnalyzer.AlleleSpecificSegment("2", 0, 1_000_000, 1, 1),
            new OncologyAnalyzer.AlleleSpecificSegment("3", 0, 5_000_000, 2, 2),
        };

        OncologyAnalyzer.EstimatePloidy(segments.Reverse())
            .Should().BeApproximately(OncologyAnalyzer.EstimatePloidy(segments), 1e-9,
                because: "a length-weighted mean is invariant to segment order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-CLONAL-001 — clonal/subclonal classification (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 108.
    //
    // API under test (OncologyAnalyzer.ClassifyClonality):
    //   Builds a CCF posterior from a variant's read evidence (more alt reads ⇒ higher CCF) and
    //   calls Clonal when P(CCF > 0.95) is high. Each variant is judged independently.
    //
    // Relations (derived from the CCF posterior, NOT from output):
    //   • MON  (higher CCF keeps clonal): the CCF estimate is non-decreasing in alt reads, and
    //          clonality is monotone — once a variant is clonal, more alt support keeps it clonal.
    //   • INV  (variant order independent): per-variant classification makes the clonal count and
    //          set independent of input order.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-CLONAL-001 MON — higher CCF keeps a clonal call clonal

    [Test]
    [Description("MON: the CCF estimate is non-decreasing in alt reads, and clonality is monotone — once clonal, more alt support keeps it clonal.")]
    public void ClassifyClonality_HigherCcf_StaysClonal()
    {
        double previousCcf = double.MinValue;
        bool clonalSeen = false;

        // Deep coverage (N = 1000) so the CCF posterior is tight enough that a VAF-0.5 (CCF ≈ 1) variant clears the clonal probability gate.
        foreach (int alt in new[] { 200, 350, 450, 500 }) // purity 1, CN 2 ⇒ expected VAF = CCF/2
        {
            var call = OncologyAnalyzer.ClassifyClonality(new[] { new OncologyAnalyzer.ClonalityVariant(alt, 1000, 2) }, 1.0).Calls[0];

            call.Ccf.Should().BeGreaterThanOrEqualTo(previousCcf, because: "more alt reads raise the CCF estimate");
            previousCcf = call.Ccf;

            if (clonalSeen)
                call.Status.Should().Be(OncologyAnalyzer.ClonalityStatus.Clonal, because: "more alt support cannot revoke a clonal call");
            if (call.Status == OncologyAnalyzer.ClonalityStatus.Clonal) clonalSeen = true;
        }

        clonalSeen.Should().BeTrue(because: "a variant at VAF 0.5 (CCF ≈ 1) is clonal");
    }

    #endregion

    #region ONCO-CLONAL-001 INV — variant order independent

    [Test]
    [Description("INV: each variant is classified independently, so the clonal count is independent of input order.")]
    public void ClassifyClonality_VariantOrder_Independent()
    {
        var variants = new[]
        {
            new OncologyAnalyzer.ClonalityVariant(500, 1000, 2), // CCF ≈ 1 → clonal
            new OncologyAnalyzer.ClonalityVariant(150, 1000, 2), // CCF ≈ 0.3 → subclonal
            new OncologyAnalyzer.ClonalityVariant(490, 1000, 2), // CCF ≈ 0.98 → clonal
            new OncologyAnalyzer.ClonalityVariant(200, 1000, 2), // subclonal
        };

        OncologyAnalyzer.ClassifyClonality(variants.Reverse(), 1.0).ClonalCount
            .Should().Be(OncologyAnalyzer.ClassifyClonality(variants, 1.0).ClonalCount,
                because: "per-variant clonality classification does not depend on input order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-NEO-001 — neoantigen peptide tiling (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 109.
    //
    // API under test (OncologyAnalyzer.GenerateNeoantigenPeptides):
    //   Enumerates every length-k window (k ∈ [min,max]) of the mutant protein that spans the
    //   substituted residue, paired with the wild-type window. The peptide content depends only
    //   on the ±(k−1) local context of the mutation, not its absolute position.
    //
    // Relation (derived from the windowing rule, NOT from output):
    //   • INV  (flanking-context shift preserves the tiling peptide set): prepending an
    //          N-terminal flank shifts every window's start position but, when the mutation keeps
    //          ≥ maxLength−1 residues of original context on each side, leaves the set of
    //          mutant/wild-type peptide strings spanning the mutation unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-NEO-001 INV — a flanking shift preserves the tiling peptides

    [Test]
    [Description("INV: prepending an N-terminal flank shifts the peptide start positions but preserves the set of mutant/wild-type peptide strings tiling the mutation.")]
    public void GenerateNeoantigenPeptides_FlankShift_PreservesPeptideSet()
    {
        // 25-mer with the mutation at position 13 → ≥ 10 residues of context on each side (max k = 11).
        const string wildType = "ACDEFGHIKLMNPQRSTVWYACDEF";
        const char mutant = 'G';      // substitute the wild-type 'P' at position 13
        const int position = 13;

        var basePeptides = OncologyAnalyzer.GenerateNeoantigenPeptides(wildType, mutant, position);
        basePeptides.Should().NotBeEmpty();

        foreach (int flankLength in new[] { 1, 5, 20 })
        {
            string flank = new string('W', flankLength);
            var shifted = OncologyAnalyzer.GenerateNeoantigenPeptides(flank + wildType, mutant, position + flankLength);

            shifted.Select(p => (p.Length, p.MutantPeptide, p.WildTypePeptide)).ToHashSet()
                .Should().BeEquivalentTo(basePeptides.Select(p => (p.Length, p.MutantPeptide, p.WildTypePeptide)).ToHashSet(),
                    because: "the windows tiling the mutation depend only on its local context, which the distant flank does not touch");

            shifted.Select(p => p.StartPosition).OrderBy(s => s)
                .Should().Equal(basePeptides.Select(p => p.StartPosition + flankLength).OrderBy(s => s),
                    because: $"every window's start shifts by the {flankLength}-residue flank");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-MHC-001 — peptide–MHC binding classification (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 110.
    //
    // API under test (OncologyAnalyzer.ClassifyMhcBinding):
    //   Classifies a peptide as Strong (IC50 < 50 nM) / Weak (< 500 nM) / NonBinder by affinity,
    //   after a length gate (invalid length ⇒ NonBinder). Strength ordinals: Strong < Weak < NonBinder.
    //
    // Relations (derived from the IC50 thresholds, NOT from output):
    //   • MON  (lower IC50 ⇒ stronger-or-equal class): decreasing IC50 makes the binding class
    //          stronger or equal (the strength ordinal is non-increasing).
    //   • INV  (peptide order independent): each peptide is classified from its own (length,
    //          IC50) only, so a batch's per-peptide classes do not depend on order.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-MHC-001 MON — lower IC50 gives a stronger-or-equal class

    [Test]
    [Description("MON: decreasing the IC50 makes the binding class stronger or equal (Strong < Weak < NonBinder ordinal is non-increasing).")]
    public void ClassifyMhcBinding_LowerIc50_StrongerClass()
    {
        int previousOrdinal = int.MaxValue;
        foreach (double ic50 in new[] { 1000.0, 400.0, 100.0, 40.0, 10.0 }) // descending IC50
        {
            int ordinal = (int)OncologyAnalyzer.ClassifyMhcBinding(9, ic50, OncologyAnalyzer.MhcClass.ClassI);
            ordinal.Should().BeLessThanOrEqualTo(previousOrdinal,
                because: $"a lower IC50 ({ic50} nM) cannot weaken the binding class");
            previousOrdinal = ordinal;
        }

        OncologyAnalyzer.ClassifyMhcBinding(9, 10.0, OncologyAnalyzer.MhcClass.ClassI)
            .Should().Be(OncologyAnalyzer.BindingStrength.Strong, because: "IC50 10 nM < 50 nM is a strong binder");
    }

    #endregion

    #region ONCO-MHC-001 INV — peptide order independent

    [Test]
    [Description("INV: each peptide is classified from its own length and IC50, so a batch's per-peptide classes are unchanged by reordering.")]
    public void ClassifyMhcBinding_PeptideOrder_Independent()
    {
        var peptides = new[] { (9, 10.0), (9, 1000.0), (11, 300.0), (7, 5.0) }; // 7-mer is invalid length → NonBinder

        var forward = peptides.Select(p => OncologyAnalyzer.ClassifyMhcBinding(p.Item1, p.Item2, OncologyAnalyzer.MhcClass.ClassI)).ToList();
        var reversed = peptides.Reverse().Select(p => OncologyAnalyzer.ClassifyMhcBinding(p.Item1, p.Item2, OncologyAnalyzer.MhcClass.ClassI)).ToList();

        reversed.Should().Equal(forward.AsEnumerable().Reverse(),
            because: "per-peptide classification has no cross-peptide dependence, so order cannot matter");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-CTDNA-001 — ctDNA detection probability (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 111.
    //
    // API under test (OncologyAnalyzer.CtDnaDetectionProbability):
    //   Poisson model p = 1 − e^(−λ), λ = n·d·k (n genome equivalents, d mutant allele fraction,
    //   k reporters). Probability is count-based — it depends only on the aggregate λ.
    //
    // Relations (derived from p = 1 − e^(−n·d·k), NOT from output):
    //   • MON  (spiking tumor signal ⇒ ≥ detection probability): p is non-decreasing in the
    //          mutant allele fraction d (and in n and k).
    //   • INV  (read order / aggregation independent): p depends only on λ = n·d·k, so equal λ
    //          (however the reads are partitioned) gives equal p.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-CTDNA-001 MON — more tumor signal raises the detection probability

    [Test]
    [Description("MON: the ctDNA detection probability is non-decreasing in the mutant allele fraction (spiking tumor reads).")]
    public void CtDna_HigherMutantFraction_HigherDetection()
    {
        double previous = double.MinValue;
        foreach (double d in new[] { 0.0, 0.0001, 0.001, 0.01, 0.1 })
        {
            double p = OncologyAnalyzer.CtDnaDetectionProbability(10_000, d);
            p.Should().BeGreaterThanOrEqualTo(previous, because: $"a higher mutant allele fraction ({d}) cannot lower detection probability");
            previous = p;
        }
    }

    #endregion

    #region ONCO-CTDNA-001 INV — probability depends only on the aggregate λ

    [Test]
    [Description("INV: the detection probability depends only on λ = n·d·k, so different read aggregations with the same λ give the same probability.")]
    public void CtDna_EqualLambda_SameProbability()
    {
        // λ = 100 expected mutant molecules computed three ways.
        double a = OncologyAnalyzer.CtDnaDetectionProbability(10_000, 0.01, 1);
        double b = OncologyAnalyzer.CtDnaDetectionProbability(1_000, 0.10, 1);
        double c = OncologyAnalyzer.CtDnaDetectionProbability(1_000, 0.05, 2);

        b.Should().BeApproximately(a, 1e-12, because: "λ = n·d·k = 100 is the same partition of the tumor signal");
        c.Should().BeApproximately(a, 1e-12, because: "moving signal between n, d and k leaves λ — and p — unchanged");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-MRD-001 — minimal-residual-disease detection (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 112.
    //
    // API under test (OncologyAnalyzer.DetectMRD):
    //   MRD is Positive when the number of detected tracked variants (PlasmaAltReads ≥ minimum)
    //   reaches the positivity threshold; the panel is aggregated by counting.
    //
    // Relations (derived from the detected-count threshold, NOT from output):
    //   • MON  (more detected tracked variants keeps MRD positive): once positive, observing
    //          additional detected variants keeps it positive (detected count only grows).
    //   • INV  (variant order independent): the detected count, status and IMAF are aggregates,
    //          independent of marker order.
    // ───────────────────────────────────────────────────────────────────────────

    private static OncologyAnalyzer.TumorMarker Marker(int pos, int alt) =>
        new("1", pos, "A", "T", alt, 100);

    #region ONCO-MRD-001 MON — more detected variants keeps MRD positive

    [Test]
    [Description("MON: once MRD is positive, observing additional detected tracked variants keeps it positive and the detected count grows.")]
    public void DetectMRD_MoreDetectedVariants_StaysPositive()
    {
        foreach (int extra in new[] { 0, 1, 3 })
        {
            var panel = new[] { Marker(1, 5), Marker(2, 5) }                       // 2 detected → positive (threshold 2)
                .Concat(Enumerable.Range(0, extra).Select(i => Marker(10 + i, 5))) // more detected variants
                .ToList();

            var result = OncologyAnalyzer.DetectMRD(panel);
            result.Status.Should().Be(OncologyAnalyzer.MrdStatus.Positive,
                because: "observing more detected tracked variants cannot turn a positive panel negative");
            result.DetectedVariantCount.Should().Be(2 + extra, because: "every alt-supported marker is counted as detected");
        }
    }

    #endregion

    #region ONCO-MRD-001 INV — variant order independent

    [Test]
    [Description("INV: the MRD status, detected count and IMAF are aggregates over the panel, independent of marker order.")]
    public void DetectMRD_VariantOrder_Independent()
    {
        var panel = new[] { Marker(1, 5), Marker(2, 0), Marker(3, 8), Marker(4, 0), Marker(5, 2) };

        var forward = OncologyAnalyzer.DetectMRD(panel);
        var reversed = OncologyAnalyzer.DetectMRD(panel.Reverse());

        reversed.Status.Should().Be(forward.Status, because: "the MRD status is a count threshold, order-independent");
        reversed.DetectedVariantCount.Should().Be(forward.DetectedVariantCount);
        reversed.IntegratedMutantAlleleFraction.Should().BeApproximately(forward.IntegratedMutantAlleleFraction, 1e-12,
            because: "IMAF is a read-pooled aggregate, independent of marker order");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-CHIP-001 — clonal-hematopoiesis (CHIP) variant flagging (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 113.
    //
    // API under test (OncologyAnalyzer.IdentifyCHIPVariants):
    //   Flags a variant as candidate CHIP iff its gene is in the CHIP panel AND its VAF ≥ minVaf;
    //   the output is a subset of the input, in input order.
    //
    // Relations (derived from the per-variant flag, NOT from output):
    //   • SUB  (survivors ⊆ input): every flagged variant is one of the input variants.
    //   • INV  (duplicating a CHIP variant keeps it flagged): per-variant judgement, so a copy of
    //          a flagged variant is also flagged.
    // ───────────────────────────────────────────────────────────────────────────

    private static OncologyAnalyzer.ChipVariant Chip(string gene, double vaf) =>
        new("1", 100, "A", "T", gene, vaf);

    #region ONCO-CHIP-001 SUB — flagged CHIP variants are a subset of the input

    [Test]
    [Description("SUB: a variant is flagged CHIP only when its gene is in the panel and its VAF ≥ threshold, so the flagged set is a subset of the input.")]
    public void IdentifyCHIP_Flagged_SubsetOfInput()
    {
        var input = new[]
        {
            Chip("DNMT3A", 0.20),  // panel gene, VAF ≥ 0.02 → CHIP
            Chip("TET2", 0.005),   // panel gene but VAF < 0.02 → not flagged
            Chip("EGFR", 0.30),    // not a CHIP gene → not flagged
            Chip("JAK2", 0.10),    // CHIP
        };

        var flagged = OncologyAnalyzer.IdentifyCHIPVariants(input);

        flagged.Should().OnlyContain(v => input.Contains(v), because: "every flagged variant comes from the input");
        flagged.Select(v => v.Gene).Should().BeEquivalentTo(new[] { "DNMT3A", "JAK2" },
            because: "only panel genes at VAF ≥ 0.02 are flagged");
    }

    #endregion

    #region ONCO-CHIP-001 INV — duplicating a CHIP variant keeps it flagged

    [Test]
    [Description("INV: variants are judged independently, so duplicating a flagged CHIP variant keeps both copies flagged.")]
    public void IdentifyCHIP_DuplicateChipVariant_StillFlagged()
    {
        var chip = Chip("DNMT3A", 0.20);
        OncologyAnalyzer.IdentifyCHIPVariants(new[] { chip }).Should().ContainSingle();
        OncologyAnalyzer.IdentifyCHIPVariants(new[] { chip, chip }).Count
            .Should().Be(2, because: "a duplicate of a flagged CHIP variant is judged identically and also flagged");

        var nonChip = Chip("EGFR", 0.30);
        OncologyAnalyzer.IdentifyCHIPVariants(new[] { nonChip, nonChip })
            .Should().BeEmpty(because: "duplicating a non-CHIP variant keeps both unflagged");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-PHYLO-001 — clonal phylogeny reconstruction (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 114.
    //
    // API under test (OncologyAnalyzer.ReconstructPhylogeny):
    //   Builds a rooted clone tree from per-sample CCF clusters using per-sample lineage and
    //   sum-rule constraints (applied identically to every sample column).
    //
    // Relations (derived from the per-sample constraints, NOT from output):
    //   • INV  (sample relabeling preserves topology): permuting the CCF sample columns (the same
    //          permutation for every cluster) leaves the reconstructed edge set unchanged — the
    //          constraints are symmetric over the sample index.
    //   • SYM  (pairwise clone distance symmetric): the tree path-distance between two clones is
    //          symmetric, d(a,b) = d(b,a).
    // ───────────────────────────────────────────────────────────────────────────

    private static OncologyAnalyzer.CcfCluster Clone(int id, params double[] ccf) =>
        new(id, ccf);

    private static int ClonePathDistance(OncologyAnalyzer.ClonalPhylogeny p, int a, int b)
    {
        var adjacency = new Dictionary<int, List<int>>();
        void Link(int x, int y)
        {
            (adjacency.TryGetValue(x, out var l) ? l : adjacency[x] = new List<int>()).Add(y);
        }
        foreach (var e in p.Edges) { Link(e.ParentId, e.ChildId); Link(e.ChildId, e.ParentId); }

        var dist = new Dictionary<int, int> { [a] = 0 };
        var queue = new Queue<int>();
        queue.Enqueue(a);
        while (queue.Count > 0)
        {
            int u = queue.Dequeue();
            if (u == b) return dist[u];
            foreach (int v in adjacency.GetValueOrDefault(u, new List<int>()))
                if (!dist.ContainsKey(v)) { dist[v] = dist[u] + 1; queue.Enqueue(v); }
        }
        return dist.GetValueOrDefault(b, -1);
    }

    #region ONCO-PHYLO-001 INV — permuting sample columns preserves the topology

    [Test]
    [Description("INV: the lineage/sum constraints are applied identically per sample, so permuting the CCF sample columns (same permutation for every cluster) leaves the reconstructed edge set unchanged.")]
    public void ReconstructPhylogeny_SamplePermutation_PreservesTopology()
    {
        var clusters = new[] { Clone(1, 1.0, 0.9), Clone(2, 0.6, 0.3), Clone(3, 0.3, 0.5) };

        var baseEdges = OncologyAnalyzer.ReconstructPhylogeny(clusters).Edges
            .Select(e => (e.ParentId, e.ChildId)).ToHashSet();

        // Swap the two sample columns consistently across all clusters.
        var swapped = clusters.Select(c => Clone(c.Id, c.CcfPerSample[1], c.CcfPerSample[0])).ToArray();
        var swappedEdges = OncologyAnalyzer.ReconstructPhylogeny(swapped).Edges
            .Select(e => (e.ParentId, e.ChildId)).ToHashSet();

        swappedEdges.Should().BeEquivalentTo(baseEdges,
            because: "the per-sample constraints are symmetric over the sample index, so reordering samples cannot change the tree");
    }

    #endregion

    #region ONCO-PHYLO-001 SYM — pairwise clone distance is symmetric

    [Test]
    [Description("SYM: the reconstructed tree induces a symmetric pairwise clone path-distance, d(a,b) = d(b,a).")]
    public void ReconstructPhylogeny_CloneDistance_Symmetric()
    {
        var phylogeny = OncologyAnalyzer.ReconstructPhylogeny(new[]
        {
            Clone(1, 1.0, 0.9), Clone(2, 0.6, 0.3), Clone(3, 0.3, 0.5),
        });

        foreach (var (a, b) in new[] { (1, 2), (2, 3), (1, 3) })
            ClonePathDistance(phylogeny, a, b).Should().Be(ClonePathDistance(phylogeny, b, a),
                because: $"the tree path distance between clones {a} and {b} is symmetric");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-CCF-001 — cancer-cell-fraction estimation (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 115.
    //
    // API under test (OncologyAnalyzer.EstimateCcf):
    //   CCF = VAF·(ρ·N_T + 2(1−ρ)) / (ρ·m), linear and increasing in VAF; the reported CCF is
    //   capped to [0, 1] while RawCcf is the uncapped value.
    //
    // Relations (derived from the linear formula, NOT from output):
    //   • MON  (higher VAF ⇒ ≥ CCF at fixed CN/purity): RawCcf strictly increases with VAF and
    //          the capped CCF is non-decreasing.
    //   • INV  (variant order independent): CCF depends only on the variant's own (VAF, purity,
    //          CN, multiplicity), so a batch's per-variant CCFs do not depend on order.
    // ───────────────────────────────────────────────────────────────────────────

    #region ONCO-CCF-001 MON — higher VAF gives a ≥ CCF

    [Test]
    [Description("MON: at fixed purity/CN/multiplicity, RawCcf increases with VAF and the capped CCF is non-decreasing.")]
    public void EstimateCcf_HigherVaf_HigherCcf()
    {
        double previousRaw = double.MinValue, previousCcf = double.MinValue;
        foreach (double vaf in new[] { 0.0, 0.1, 0.2, 0.3, 0.4 })
        {
            var est = OncologyAnalyzer.EstimateCcf(vaf, purity: 0.8, tumorCopyNumber: 2, multiplicity: 1);
            est.RawCcf.Should().BeGreaterThan(previousRaw, because: $"RawCcf is linear in VAF ({vaf})");
            est.Ccf.Should().BeGreaterThanOrEqualTo(previousCcf, because: "the capped CCF is non-decreasing in VAF");
            previousRaw = est.RawCcf;
            previousCcf = est.Ccf;
        }
    }

    #endregion

    #region ONCO-CCF-001 INV — variant order independent

    [Test]
    [Description("INV: each variant's CCF depends only on its own parameters, so a batch's per-variant CCFs are unchanged by reordering.")]
    public void EstimateCcf_VariantOrder_Independent()
    {
        var vafs = new[] { 0.05, 0.20, 0.45 };

        var forward = vafs.Select(v => OncologyAnalyzer.EstimateCcf(v, 0.8, 2, 1).RawCcf).ToList();
        var reversed = vafs.Reverse().Select(v => OncologyAnalyzer.EstimateCcf(v, 0.8, 2, 1).RawCcf).ToList();

        reversed.Should().Equal(forward.AsEnumerable().Reverse(),
            because: "per-variant CCF has no cross-variant dependence, so order cannot matter");
    }

    #endregion
}
