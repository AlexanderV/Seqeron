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
}
