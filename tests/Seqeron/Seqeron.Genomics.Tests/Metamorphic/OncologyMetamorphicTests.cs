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
}
