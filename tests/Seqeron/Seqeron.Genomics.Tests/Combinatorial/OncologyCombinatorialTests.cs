namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Oncology area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids (3³ = 27) use the exhaustive <c>[Combinatorial]</c> product,
/// a strict superset of pairwise.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Oncology")]
public class OncologyCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-IMMUNE-001 — Immune infiltration estimation (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 86.
    // Spec: tests/TestSpecs/ONCO-IMMUNE-001.md (ImmuneAnalyzer.EstimateInfiltration / DeconvoluteImmuneCells).
    // ADVANCED_TESTING_CHECKLIST.md §10 (Combinatorial / Pairwise).
    //
    // Sources:
    //   Yoshihara et al. (2013), Nat. Commun. 4:2612 — ESTIMATE immune/stromal scoring + purity.
    //   Barbie et al. (2009), Nature 462:108-112 — ssGSEA enrichment (rank-weighted, τ=0.25).
    //   Hänzelmann et al. (2013), BMC Bioinformatics 14:7 — GSVA ssGSEA integral form.
    //   Lawson & Hanson (1995), SIAM, Ch.23 — NNLS. Abbas et al. (2009), PLoS One 4:e6098 — deconvolution.
    //   Newman et al. (2015), Nat. Methods 12:453-457 — LM22 cell-type framework.
    //
    // Checklist axes for this row are geneSetSize(3) × normalization(3) × nPermutations(3).
    // Both ImmuneAnalyzer methods are CLOSED-FORM and DETERMINISTIC: there is no permutation
    // resampling step and no tunable normalization "mode" knob in the API. Per the campaign
    // convention we map the nominal axes onto the real, implemented knobs and document it here:
    //
    //   • geneSetSize  → realised in BOTH grids: ssGSEA grid varies the number of signature
    //                    (hit) genes nHits ∈ {1,2,3}; NNLS grid varies the cell-type signatures
    //                    that compose the mixture.
    //   • normalization→ the two normalization regimes that genuinely exist: the ssGSEA
    //                    rank-weight normalization (Σ rank^τ; grid 1) and the NNLS Σf = 1
    //                    simplex normalization (grid 2, asserted in every cell, INV-2).
    //   • nPermutations→ no permutation step exists; replaced by the structural axes that DO
    //                    drive the output — profile size N and hit placement (grid 1, ssGSEA)
    //                    and the mixture weight (grid 2, NNLS).
    //
    // Each grid is 3 × 3 × 3 = 27 cells = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: the enrichment / deconvolution output is a JOINT function of all
    // three axes at once — no single axis explains a cell, and flipping one axis (hit placement,
    // or mixture weight) flips the result. The interaction-witness tests prove the suite genuinely
    // needs the combination.
    // ═══════════════════════════════════════════════════════════════════════

    private const double SsGseaTau = 0.25; // Barbie et al. (2009); GSVA default.

    /// <summary>Where the signature (hit) genes sit in the descending-expression ranking.</summary>
    public enum HitPlacement
    {
        /// <summary>Hits occupy the highest-expression ranks → positive enrichment.</summary>
        Top,

        /// <summary>Hits occupy the lowest-expression ranks → negative enrichment.</summary>
        Bottom,

        /// <summary>Hits spread evenly across the ranking.</summary>
        Spread,
    }

    // ───────────────────────────────────────────────────────────────────────
    // Grid 1 — ssGSEA enrichment (EstimateInfiltration)
    // Axes: profileSize N(3) × geneSetSize nHits(3) × hitPlacement(3) = 27.
    // ───────────────────────────────────────────────────────────────────────

    /// <summary>
    /// For every (N, nHits, placement) the production ssGSEA immune score must equal the
    /// rank-based enrichment integral computed independently from its mathematical definition —
    /// ES = Σ_p RS(p) over the descending-ranked walk, hits weighted by rank^τ / Σrank^τ and
    /// misses by 1/N_miss (Barbie et al. 2009; Hänzelmann et al. 2013). Stromal genes are empty,
    /// so EstimateScore collapses to the immune score (definition cross-check).
    /// </summary>
    [Test, Combinatorial]
    public void EstimateInfiltration_SsGseaGrid_MatchesRankBasedIntegral(
        [Values(4, 6, 8)] int profileSize,
        [Values(1, 2, 3)] int geneSetSize,
        [Values] HitPlacement placement)
    {
        var (profile, hitGenes, values, hitIdx) = BuildEnrichmentProfile(profileSize, geneSetSize, placement);

        double expected = GroundTruthSsGseaIntegral(values, hitIdx);

        var result = ImmuneAnalyzer.EstimateInfiltration(profile, hitGenes, Array.Empty<string>());

        result.ImmuneScore.Should().BeApproximately(expected, 1e-9,
            "production ssGSEA must equal the rank-based integral defined by Barbie/Hänzelmann");
        result.StromalScore.Should().Be(0.0, "no stromal genes overlap the profile");
        result.EstimateScore.Should().BeApproximately(expected, 1e-9,
            "ESTIMATE score = immune + stromal collapses to the immune score when stromal = 0");
        result.OverlappingImmuneGenes.Should().Be(geneSetSize, "all hit genes are present in the profile");
    }

    /// <summary>
    /// Interaction witness (placement axis): for a fixed profile and gene-set size, concentrating
    /// the signature among the highest-expression genes yields positive enrichment, among the
    /// lowest-expression genes yields negative enrichment. This is the biological meaning of
    /// ssGSEA — a result that flips purely on the placement axis, proving the grid needs it.
    /// </summary>
    [Test]
    public void EstimateInfiltration_HitPlacement_FlipsEnrichmentSign()
    {
        var top = BuildEnrichmentProfile(8, 2, HitPlacement.Top);
        var bottom = BuildEnrichmentProfile(8, 2, HitPlacement.Bottom);

        double topScore = ImmuneAnalyzer.EstimateInfiltration(top.profile, top.hitGenes, Array.Empty<string>()).ImmuneScore;
        double bottomScore = ImmuneAnalyzer.EstimateInfiltration(bottom.profile, bottom.hitGenes, Array.Empty<string>()).ImmuneScore;

        topScore.Should().BeGreaterThan(0.0, "a set concentrated in highly-expressed genes is positively enriched");
        bottomScore.Should().BeLessThan(0.0, "a set concentrated in lowly-expressed genes is negatively enriched");
        topScore.Should().BeGreaterThan(bottomScore, "enrichment is monotone in the hit ranks");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Grid 2 — NNLS deconvolution (DeconvoluteImmuneCells)
    // Axes: cellTypeA(3) × cellTypeB(3) × mixtureWeight(3) = 27.
    // ───────────────────────────────────────────────────────────────────────

    /// <summary>
    /// For every two-component mixture m = wₐ·Sₐ + w_b·S_b (a point inside the column space of the
    /// signature matrix), NNLS must satisfy its defining guarantees across the whole grid:
    /// non-negativity (f ≥ 0, INV-1), simplex normalization (Σf = 1, INV-2), and a zero-residual
    /// reconstruction (RMSE → 0, correlation → 1) because the global least-squares minimum is 0.
    /// Source: Lawson & Hanson (1995); Abbas et al. (2009).
    /// </summary>
    [Test, Combinatorial]
    public void DeconvoluteImmuneCells_TwoComponentMixtureGrid_SatisfiesNnlsInvariants(
        [Values("T_cells_CD8", "Monocytes", "Plasma_cells")] string cellTypeA,
        [Values("B_cells_naive", "Neutrophils", "Mast_cells_resting")] string cellTypeB,
        [Values(0.25, 0.5, 0.75)] double weightA)
    {
        var profile = BuildMixtureProfile(cellTypeA, weightA, cellTypeB, 1.0 - weightA);

        var result = ImmuneAnalyzer.DeconvoluteImmuneCells(profile);

        result.CellFractions.Values.Should().OnlyContain(f => f >= -1e-9,
            "[INV-1] NNLS enforces non-negative fractions");
        result.CellFractions.Values.Sum().Should().BeApproximately(1.0, 1e-6,
            "[INV-2] fractions are normalized onto the simplex Σf = 1");
        result.Rmse.Should().BeLessThan(1e-6,
            "the mixture lies in the column space, so the least-squares residual is 0");
        result.Correlation.Should().BeApproximately(1.0, 1e-6,
            "zero-residual reconstruction means observed and reconstructed profiles are identical");
        result.OverlappingGenes.Should().BeGreaterThan(0, "the mixture genes overlap the signature matrix");
    }

    /// <summary>
    /// Interaction witness (cell-type × weight): when both components carry marker genes unique to
    /// their column, NNLS recovers the exact mixing weights, and flipping the weight axis flips the
    /// recovered fractions. CD8A/CD8B are unique to T_cells_CD8; CD79A/PAX5/BANK1 are unique to
    /// B_cells_naive — so the zero-residual solution is unique. Source: Abbas et al. (2009) linearity.
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_UniqueMarkerMixture_RecoversWeights_FlipsOnWeightAxis()
    {
        var mix25 = BuildMixtureProfile("T_cells_CD8", 0.25, "B_cells_naive", 0.75);
        var mix75 = BuildMixtureProfile("T_cells_CD8", 0.75, "B_cells_naive", 0.25);

        var r25 = ImmuneAnalyzer.DeconvoluteImmuneCells(mix25);
        var r75 = ImmuneAnalyzer.DeconvoluteImmuneCells(mix75);

        r25.CellFractions["T_cells_CD8"].Should().BeApproximately(0.25, 1e-6);
        r25.CellFractions["B_cells_naive"].Should().BeApproximately(0.75, 1e-6);
        r75.CellFractions["T_cells_CD8"].Should().BeApproximately(0.75, 1e-6);
        r75.CellFractions["B_cells_naive"].Should().BeApproximately(0.25, 1e-6);

        foreach (var (cellType, fraction) in r25.CellFractions)
        {
            if (cellType is not ("T_cells_CD8" or "B_cells_naive"))
                fraction.Should().BeApproximately(0.0, 1e-6, "non-mixture components are absent");
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-SOMATIC-001 — Somatic mutation calling (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 87.
    // Spec: tests/TestSpecs/ONCO-SOMATIC-001.md (OncologyAnalyzer.CallSomaticMutations / Classify).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Saunders et al. (2012) Strelka; Kim et al. (2018) Strelka2; Benjamin et al. (2019)
    // Mutect2; Yan et al. (2021) Sci. Rep. 11:11640 (5% WES limit of detection).
    //
    // The caller is a deterministic two-threshold rule over (f_t, f_n):
    //   f_t < tumorThreshold              → NotDetected   (below the tumor LoD; Yan 2021)
    //   f_t ≥ tumorThreshold, f_n ≤ 0.01  → Somatic       (present in tumor, ref/ref normal; Saunders 2012)
    //   f_t ≥ tumorThreshold, f_n > 0.01  → Germline      (present in matched normal; Benjamin 2019)
    //   score = max(0, f_t − f_n) for somatic, else 0.
    //
    // Checklist axes minVaf(3) × minDepth(3) × strandBias(2) — the caller has no depth or strand
    // parameter (those belong to artifact filtering, ONCO-ARTIFACT-001 row 90). Per the campaign
    // convention we map the nominal axes onto the real knobs and document it:
    //   • minVaf      → tumorVafThreshold ∈ {0.02, 0.05, 0.10} (the real parameter).
    //   • minDepth    → tumor sequencing depth ∈ {50, 100, 200}. With the alt-supporting count held
    //                   FIXED (6 reads), depth genuinely drives f_t = 6/depth, so depth × threshold
    //                   jointly decide presence — the real "depth matters for low-count calls" effect.
    //   • strandBias  → the binary normal-channel artifact axis the somatic rule actually responds to:
    //                   {Clean normal (f_n = 0, ref/ref) → somatic, Contaminated normal (f_n = 0.20,
    //                   CHIP/germline leakage) → germline}.
    // Grid = 3 × 3 × 2 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: status is a JOINT function of all three axes — no single axis fixes
    // the outcome (each one appears in at least one status flip). Every cell asserts the full
    // SomaticCall against the rule re-derived independently from the realized read counts.
    // ═══════════════════════════════════════════════════════════════════════

    private const int SomaticFixedAltReads = 6;        // alt-supporting reads held constant across depths
    private const double NormalNoiseCeiling = 0.01;    // DefaultNormalVafThreshold
    private const double ContaminatedNormalVaf = 0.20; // CHIP/germline leakage well above the noise ceiling

    /// <summary>Matched-normal channel state — the binary "artifact" axis the somatic rule responds to.</summary>
    public enum NormalState
    {
        /// <summary>Homozygous-reference normal (f_n = 0): a true somatic site.</summary>
        Clean,

        /// <summary>Normal carries the allele above the noise ceiling (f_n = 0.20): germline/CHIP leakage.</summary>
        Contaminated,
    }

    /// <summary>
    /// For every (tumorThreshold, depth, normalState) the production classification, VAFs and somatic
    /// score must equal the two-threshold rule re-derived from the realized read counts. With the alt
    /// count fixed at 6, depth changes f_t = 6/depth so the depth and threshold axes jointly decide
    /// presence; the normal axis decides somatic-vs-germline once present.
    /// </summary>
    [Test, Combinatorial]
    public void CallSomaticMutations_VafDepthNormalGrid_MatchesTwoThresholdRule(
        [Values(0.02, 0.05, 0.10)] double tumorThreshold,
        [Values(50, 100, 200)] int depth,
        [Values] NormalState normalState)
    {
        int normalAlt = normalState == NormalState.Clean ? 0 : (int)Math.Round(ContaminatedNormalVaf * depth);
        var variant = new OncologyAnalyzer.VariantObservation(
            "chr1", 1000, "C", "T",
            TumorAltReads: SomaticFixedAltReads, TumorTotalReads: depth,
            NormalAltReads: normalAlt, NormalTotalReads: depth);

        // Independent ground truth from the realized read counts.
        double tumorVaf = (double)SomaticFixedAltReads / depth;
        double normalVaf = (double)normalAlt / depth;
        OncologyAnalyzer.SomaticStatus expectedStatus;
        if (tumorVaf < tumorThreshold)
            expectedStatus = OncologyAnalyzer.SomaticStatus.NotDetected;
        else if (normalVaf <= NormalNoiseCeiling)
            expectedStatus = OncologyAnalyzer.SomaticStatus.Somatic;
        else
            expectedStatus = OncologyAnalyzer.SomaticStatus.Germline;
        double expectedScore = expectedStatus == OncologyAnalyzer.SomaticStatus.Somatic ? tumorVaf - normalVaf : 0.0;

        var call = OncologyAnalyzer.Classify(variant, tumorThreshold);

        call.Status.Should().Be(expectedStatus, "the two-threshold rule decides status jointly from f_t and f_n");
        call.TumorVaf.Should().BeApproximately(tumorVaf, 1e-12);
        call.NormalVaf.Should().BeApproximately(normalVaf, 1e-12);
        call.SomaticScore.Should().BeApproximately(expectedScore, 1e-12, "somatic score = max(0, f_t − f_n)");
    }

    /// <summary>
    /// Interaction witness (depth × threshold): with the alt count fixed at 6 and a strict 10% tumor
    /// threshold, the very same allele is called Somatic at 50× (f_t = 0.12) but drops to NotDetected
    /// at 200× (f_t = 0.03). The status flips purely on the depth axis — the grid genuinely needs it.
    /// </summary>
    [Test]
    public void CallSomaticMutations_DepthAxis_FlipsPresence()
    {
        var shallow = new OncologyAnalyzer.VariantObservation("chr1", 1, "C", "T", 6, 50, 0, 50);
        var deep = new OncologyAnalyzer.VariantObservation("chr1", 1, "C", "T", 6, 200, 0, 200);

        OncologyAnalyzer.Classify(deep, tumorVafThreshold: 0.10).Status
            .Should().Be(OncologyAnalyzer.SomaticStatus.NotDetected, "6/200 = 3% is below the 10% LoD");
        OncologyAnalyzer.Classify(shallow, tumorVafThreshold: 0.10).Status
            .Should().Be(OncologyAnalyzer.SomaticStatus.Somatic, "6/50 = 12% clears the 10% LoD with a clean normal");
    }

    /// <summary>
    /// Interaction witness (normal axis): a present tumor allele is Somatic against a clean normal but
    /// Germline once the matched normal carries the allele above the noise ceiling. Source: Mutect2
    /// germline filtering (Benjamin et al. 2019).
    /// </summary>
    [Test]
    public void CallSomaticMutations_NormalAxis_FlipsSomaticVsGermline()
    {
        var clean = new OncologyAnalyzer.VariantObservation("chr1", 1, "C", "T", 30, 100, 0, 100);
        var contaminated = new OncologyAnalyzer.VariantObservation("chr1", 1, "C", "T", 30, 100, 20, 100);

        OncologyAnalyzer.Classify(clean).Status
            .Should().Be(OncologyAnalyzer.SomaticStatus.Somatic, "ref/ref normal → somatic");
        OncologyAnalyzer.Classify(contaminated).Status
            .Should().Be(OncologyAnalyzer.SomaticStatus.Germline, "normal carries the allele at 20% → germline");
    }

    /// <summary>
    /// Worked example (spec M1): a clear somatic SNV, tumor 25/100 against a clean normal 0/100, is
    /// classified Somatic with score = f_t − f_n = 0.25. Source: Saunders et al. (2012).
    /// </summary>
    [Test]
    public void CallSomaticMutations_ClearSomatic_WorkedExample()
    {
        var variant = new OncologyAnalyzer.VariantObservation("chr7", 140_453_136, "A", "T", 25, 100, 0, 100);

        var call = OncologyAnalyzer.Classify(variant);

        call.Status.Should().Be(OncologyAnalyzer.SomaticStatus.Somatic);
        call.TumorVaf.Should().BeApproximately(0.25, 1e-12);
        call.SomaticScore.Should().BeApproximately(0.25, 1e-12);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-DRIVER-001 — Driver mutation detection, Vogelstein 20/20 rule (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 89.
    // Spec: tests/TestSpecs/ONCO-DRIVER-001.md (OncologyAnalyzer.IdentifyDriverMutations / ClassifyGene).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Vogelstein et al. (2013) Science 339:1546-1558 (20/20 rule, strict ">20%"); Tokheim &
    // Karchin (2020); Schroeder et al. (2014) (truncating categories); Miller et al. (2017) (a recurrent
    // position needs ≥2 mutations at the same codon).
    //
    // A mutation is reported by IdentifyDriverMutations iff its gene classifies as a driver
    // (Oncogene: recurrent-missense fraction > 0.20; TumorSuppressor: truncating fraction > 0.20) OR its
    // (gene, position) is in the caller-supplied hotspot set.
    //
    // Checklist axes recurrence(3) × hotspot(2) × geneList(2) map onto the real predicate inputs:
    //   • recurrence → copies of the query missense at the SAME codon, k ∈ {1, 2, 4}: 1 = not recurrent
    //     (needs ≥2); 2 = recurrent but fraction exactly 0.20 (NOT > 0.20, strict); 4 = fraction 0.40 →
    //     Oncogene. Exercises the strict-inequality OG criterion (INV-2, INV-5).
    //   • geneList → whether the SAME gene also carries a truncating block making it a TumorSuppressor
    //     independent of recurrence (the data-derived "driver gene list"): {Absent, Present}.
    //   • hotspot → whether the query (gene, position) is in the caller hotspot set: {Absent, Present}
    //     — the rescue path that reports a mutation even when its gene is Ambiguous.
    // Grid = 3 × 2 × 2 = 12 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: each axis independently can turn a non-driver into a reported driver —
    // recurrence (k=4), the TSG background, and the hotspot rescue are three separate routes. The role
    // and the reported-membership are checked against the 20/20 rule re-derived from the construction.
    // ═══════════════════════════════════════════════════════════════════════

    private const string FocalGene = "GENEX";
    private const int FocalPosition = 100;
    private const int DriverSpectrumTotal = 10; // fixed denominator → predictable 20/20 fractions

    /// <summary>Recurrence of the query missense codon — drives the oncogene (recurrent-missense) criterion.</summary>
    public enum Recurrence
    {
        /// <summary>k = 1 copy: a singleton, not a recurrent position (Miller 2017 needs ≥2).</summary>
        NotRecurrent,

        /// <summary>k = 2 copies: recurrent, fraction exactly 0.20 — NOT &gt; 0.20 (strict threshold).</summary>
        AtThreshold,

        /// <summary>k = 4 copies: fraction 0.40 &gt; 0.20 → Oncogene.</summary>
        AboveThreshold,
    }

    private static int RecurrenceCopies(Recurrence r) => r switch
    {
        Recurrence.NotRecurrent => 1,
        Recurrence.AtThreshold => 2,
        Recurrence.AboveThreshold => 4,
        _ => throw new ArgumentOutOfRangeException(nameof(r)),
    };

    /// <summary>
    /// For every (recurrence, TSG background, hotspot) the gene's 20/20 role and whether the query
    /// missense is reported as a driver must match the rule re-derived from the construction. The query
    /// is reported iff its gene is a driver OR its position is a known hotspot; the result is always a
    /// subset of the input.
    /// </summary>
    [Test, Combinatorial]
    public void IdentifyDriverMutations_RecurrenceTsgHotspotGrid_MatchesTwentyTwentyRule(
        [Values] Recurrence recurrence,
        [Values(false, true)] bool tsgBackground,
        [Values(false, true)] bool hotspot)
    {
        int k = RecurrenceCopies(recurrence);
        var spectrum = BuildDriverGeneSpectrum(k, tsgBackground);
        var query = new OncologyAnalyzer.GeneMutation(FocalGene, FocalPosition, OncologyAnalyzer.MutationConsequence.Missense);
        var hotspots = hotspot
            ? new HashSet<(string, int)> { (FocalGene, FocalPosition) }
            : null;

        // Independent ground truth (20/20 rule applied to the known construction; fixed denominator 10).
        double recurrentMissenseFraction = (k >= 2 ? k : 0) / (double)DriverSpectrumTotal;
        double truncatingFraction = (tsgBackground ? 5 : 0) / (double)DriverSpectrumTotal;
        bool isOg = recurrentMissenseFraction > 0.20;
        bool isTsg = truncatingFraction > 0.20;
        OncologyAnalyzer.DriverGeneRole expectedRole =
            isTsg && isOg ? (truncatingFraction > recurrentMissenseFraction ? OncologyAnalyzer.DriverGeneRole.TumorSuppressor
                           : recurrentMissenseFraction > truncatingFraction ? OncologyAnalyzer.DriverGeneRole.Oncogene
                           : OncologyAnalyzer.DriverGeneRole.Ambiguous)
            : isTsg ? OncologyAnalyzer.DriverGeneRole.TumorSuppressor
            : isOg ? OncologyAnalyzer.DriverGeneRole.Oncogene
            : OncologyAnalyzer.DriverGeneRole.Ambiguous;
        bool expectedReported = expectedRole != OncologyAnalyzer.DriverGeneRole.Ambiguous || hotspot;

        var classification = OncologyAnalyzer.ClassifyGene(spectrum);
        var drivers = OncologyAnalyzer.IdentifyDriverMutations(spectrum, hotspots);

        classification.Role.Should().Be(expectedRole, "the 20/20 rule decides the gene role");
        classification.RecurrentMissenseFraction.Should().BeApproximately(recurrentMissenseFraction, 1e-12);
        classification.TruncatingFraction.Should().BeApproximately(truncatingFraction, 1e-12);
        drivers.Contains(query).Should().Be(expectedReported,
            "a mutation is reported iff its gene is a driver or its position is a hotspot");
        drivers.Should().OnlyContain(m => spectrum.Contains(m), "[INV-1] drivers ⊆ input");
    }

    /// <summary>
    /// Interaction witness: an otherwise-Ambiguous gene (a singleton missense, no truncating block) is
    /// NOT reported on its own, yet each of the three axes independently rescues it — bumping recurrence
    /// to k=4 (Oncogene), adding a truncating block (TumorSuppressor), or listing the position as a
    /// hotspot. A result that flips on every axis proves the grid genuinely needs the combination.
    /// </summary>
    [Test]
    public void IdentifyDriverMutations_EachAxisIndependentlyRescuesAmbiguousGene()
    {
        var query = new OncologyAnalyzer.GeneMutation(FocalGene, FocalPosition, OncologyAnalyzer.MutationConsequence.Missense);

        var ambiguous = BuildDriverGeneSpectrum(1, tsgBackground: false);
        OncologyAnalyzer.IdentifyDriverMutations(ambiguous).Should().NotContain(query,
            "a singleton missense in an Ambiguous gene is not a driver");

        OncologyAnalyzer.IdentifyDriverMutations(BuildDriverGeneSpectrum(4, tsgBackground: false))
            .Should().Contain(query, "recurrence axis: k=4 → Oncogene");
        OncologyAnalyzer.IdentifyDriverMutations(BuildDriverGeneSpectrum(1, tsgBackground: true))
            .Should().Contain(query, "geneList axis: truncating block → TumorSuppressor");
        OncologyAnalyzer.IdentifyDriverMutations(ambiguous, new HashSet<(string, int)> { (FocalGene, FocalPosition) })
            .Should().Contain(query, "hotspot axis: caller hotspot rescues the Ambiguous gene");
    }

    /// <summary>
    /// Worked example (spec M1): IDH1-like gene, all 10 mutations missense at codon 132 → recurrent-
    /// missense fraction 1.00, Oncogene. Source: Vogelstein et al. (2013); Miller et al. (2017).
    /// </summary>
    [Test]
    public void ClassifyGene_AllRecurrentMissense_OncogeneWorkedExample()
    {
        var mutations = Enumerable.Range(0, 10)
            .Select(_ => new OncologyAnalyzer.GeneMutation("IDH1", 132, OncologyAnalyzer.MutationConsequence.Missense))
            .ToList();

        var c = OncologyAnalyzer.ClassifyGene(mutations);

        c.Role.Should().Be(OncologyAnalyzer.DriverGeneRole.Oncogene);
        c.RecurrentMissenseFraction.Should().BeApproximately(1.0, 1e-12);
        c.TruncatingFraction.Should().BeApproximately(0.0, 1e-12);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-ARTIFACT-001 — Sequencing artifact detection (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 90.
    // Spec: tests/TestSpecs/ONCO-ARTIFACT-001.md (OncologyAnalyzer.ClassifyArtifact / FilterArtifacts).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Chen et al. (2017) Science 355:752-756 (OxoG G>T excess, GIV); Nature Methods (2017)
    // (GIV=1 undamaged, >1.5 damaged); Do & Dobrovic (2015) (FFPE C>T/G>A deamination, disjoint from
    // OxoG); GATK FisherStrand (FS = −10·log10 two-sided Fisher p).
    //
    // IsArtifact rule:  FFPE (C>T/G>A) → always artifact;  OxoG (G>T/C>A) → artifact iff GIV > 1.5;
    // any other substitution → not an artifact. The Fisher strand-bias FS is computed and reported as a
    // separate output (it is monotone in strand segregation; it does NOT gate IsArtifact here).
    //
    // Checklist axes strandBias(3) × baseQual(3) × position(2) — the model has no base-quality or read-
    // position field. Per the campaign convention we map them onto the real knobs and document it:
    //   • strandBias → the 2×2 strand contingency table {Balanced [20,20,20,20] → FS=0, Mild [20,20,16,4],
    //     Strong [30,0,0,30]}: the real FisherStrand input; drives FS monotonically (INV-3, INV-5).
    //   • baseQual   → the GIV read-orientation imbalance {Undamaged 100/100=1.0, Boundary 150/100=1.5,
    //     Damaged 200/100=2.0}: the per-read-pair evidence that gates the OxoG call (strict > 1.5).
    //   • position   → the artifact substitution FAMILY {FFPE deamination C>T, OxoG oxidation G>T}: the two
    //     disjoint artifact classes (INV-4). The non-artifact class (e.g. A>G) is covered by a witness.
    // Grid = 3 × 3 × 2 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: IsArtifact is a JOINT function of family and GIV — FFPE ignores GIV while
    // OxoG flips on the GIV axis at the strict 1.5 boundary — while FS depends only on the strand table.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>Artifact substitution family — the two disjoint artifact classes (Do & Dobrovic; Chen).</summary>
    public enum ArtifactFamily
    {
        /// <summary>FFPE cytosine-deamination C&gt;T: always an artifact by substitution class.</summary>
        FfpeDeamination,

        /// <summary>OxoG oxidation G&gt;T: an artifact only when the GIV imbalance exceeds 1.5.</summary>
        OxoG,
    }

    /// <summary>GIV read-orientation imbalance level (the OxoG-gating evidence).</summary>
    public enum GivLevel
    {
        /// <summary>R1 = R2 → GIV 1.0 (balanced, undamaged).</summary>
        Undamaged,

        /// <summary>150/100 → GIV exactly 1.5: NOT &gt; 1.5 (strict damaged threshold).</summary>
        Boundary,

        /// <summary>200/100 → GIV 2.0 &gt; 1.5 (damaged).</summary>
        Damaged,
    }

    /// <summary>Strand contingency-table segregation level (drives the FisherStrand FS).</summary>
    public enum StrandSegregation
    {
        /// <summary>[20,20,20,20]: perfectly balanced, Fisher p = 1 → FS = 0.</summary>
        Balanced,

        /// <summary>[20,20,16,4]: moderate strand skew → FS &gt; 0.</summary>
        Mild,

        /// <summary>[30,0,0,30]: complete strand segregation → large FS.</summary>
        Strong,
    }

    /// <summary>
    /// For every (family, GIV level, strand segregation) the artifact classification must match the rule:
    /// FFPE is always an artifact, OxoG is an artifact only at GIV &gt; 1.5, and the GIV score / FS are the
    /// expected values. The grid checks that FFPE ignores GIV while OxoG flips on it, and that FS = 0 only
    /// for the balanced strand table; FilterArtifacts keeps a variant iff it is not flagged (INV-1 subset).
    /// </summary>
    [Test, Combinatorial]
    public void ClassifyArtifact_StrandGivFamilyGrid_MatchesArtifactRule(
        [Values] StrandSegregation strand,
        [Values] GivLevel giv,
        [Values] ArtifactFamily family)
    {
        var obs = BuildArtifactObservation(family, giv, strand);

        double expectedGiv = giv switch { GivLevel.Undamaged => 1.0, GivLevel.Boundary => 1.5, GivLevel.Damaged => 2.0, _ => double.NaN };
        var expectedType = family == ArtifactFamily.FfpeDeamination
            ? OncologyAnalyzer.ArtifactType.FfpeDeamination
            : OncologyAnalyzer.ArtifactType.OxoG;
        bool expectedIsArtifact = family == ArtifactFamily.FfpeDeamination || expectedGiv > 1.5;

        var call = OncologyAnalyzer.ClassifyArtifact(obs);

        call.Type.Should().Be(expectedType, "substitution maps to its disjoint artifact class (INV-4)");
        call.GivScore.Should().BeApproximately(expectedGiv, 1e-12);
        call.IsArtifact.Should().Be(expectedIsArtifact,
            "FFPE is always an artifact; OxoG only when GIV > 1.5 (strict)");
        call.StrandBiasPhred.Should().BeGreaterThanOrEqualTo(0.0, "[INV-3] FS ≥ 0");
        if (strand == StrandSegregation.Balanced)
            call.StrandBiasPhred.Should().Be(0.0, "[INV-3] a balanced table has Fisher p = 1 → FS = 0");

        var kept = OncologyAnalyzer.FilterArtifacts(new[] { obs });
        kept.Should().BeSubsetOf(new[] { obs }, "[INV-1] FilterArtifacts ⊆ input");
        kept.Should().HaveCount(expectedIsArtifact ? 0 : 1, "a variant is kept iff it is not flagged");
    }

    /// <summary>
    /// Interaction witness (GIV axis flips OxoG, but not FFPE): the same OxoG G>T variant is kept at the
    /// strict 1.5 boundary (GIV = 1.5 is NOT damaged) yet removed at GIV = 2.0, whereas an FFPE C>T variant
    /// is removed at both — proving the GIV axis only matters in combination with the OxoG family.
    /// </summary>
    [Test]
    public void FilterArtifacts_GivAxis_FlipsOxoGButNotFfpe()
    {
        var oxoBoundary = BuildArtifactObservation(ArtifactFamily.OxoG, GivLevel.Boundary, StrandSegregation.Balanced);
        var oxoDamaged = BuildArtifactObservation(ArtifactFamily.OxoG, GivLevel.Damaged, StrandSegregation.Balanced);
        var ffpeBoundary = BuildArtifactObservation(ArtifactFamily.FfpeDeamination, GivLevel.Boundary, StrandSegregation.Balanced);

        OncologyAnalyzer.FilterArtifacts(new[] { oxoBoundary }).Should().ContainSingle("GIV 1.5 is not > 1.5 → OxoG kept");
        OncologyAnalyzer.FilterArtifacts(new[] { oxoDamaged }).Should().BeEmpty("GIV 2.0 > 1.5 → OxoG removed");
        OncologyAnalyzer.FilterArtifacts(new[] { ffpeBoundary }).Should().BeEmpty("FFPE is removed regardless of GIV");
    }

    /// <summary>
    /// Interaction witness (INV-5 monotonicity): the Fisher strand-bias FS is non-decreasing as the strand
    /// table segregates — 0 for a balanced table, larger for moderate skew, largest for complete
    /// segregation. Source: GATK FisherStrand (the Fisher p decreases as alleles segregate by strand).
    /// </summary>
    [Test]
    public void CalculateStrandBias_IncreasingSegregation_FsMonotone()
    {
        double balanced = OncologyAnalyzer.CalculateStrandBias(20, 20, 20, 20);
        double mild = OncologyAnalyzer.CalculateStrandBias(20, 20, 16, 4);
        double strong = OncologyAnalyzer.CalculateStrandBias(30, 0, 0, 30);

        balanced.Should().Be(0.0, "balanced table → p = 1 → FS = 0");
        mild.Should().BeGreaterThan(balanced);
        strong.Should().BeGreaterThan(mild);
    }

    /// <summary>
    /// Witness (non-artifact class): a transition substitution outside the artifact classes (A>G) is never
    /// flagged and is always kept, regardless of strand/GIV evidence. Source: artifact classes are specific
    /// (Do & Dobrovic 2015; Chen et al. 2017).
    /// </summary>
    [Test]
    public void FilterArtifacts_NonArtifactSubstitution_AlwaysKept()
    {
        var realVariant = new OncologyAnalyzer.ArtifactObservation('A', 'G', 30, 0, 0, 30, 200, 100);

        OncologyAnalyzer.ClassifyArtifact(realVariant).Type.Should().Be(OncologyAnalyzer.ArtifactType.None);
        OncologyAnalyzer.FilterArtifacts(new[] { realVariant }).Should().ContainSingle("A>G is a candidate true variant");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-TMB-001 — Tumor mutational burden (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 92.
    // Spec: tests/TestSpecs/ONCO-TMB-001.md (OncologyAnalyzer.CalculateTMB / ClassifyTMB).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Chalmers et al. (2017) Genome Medicine 9:34 (TMB = somatic coding mutations / Mb;
    // FoundationOne 315-gene panel = 1.1 Mb); Marcus et al. (2021) FDA pembrolizumab approval
    // (TMB-High ⟺ TMB ≥ 10 mut/Mb, inclusive).
    //
    // TMB = somaticCount / targetRegionMb;  ClassifyTMB → High ⟺ TMB ≥ 10. The SomaticCall overload
    // counts ONLY Somatic-status calls (Germline / NotDetected are excluded).
    //
    // Checklist axes panelSize(3) × mutationCount(3) × includeSilent(2) map onto the real knobs:
    //   • panelSize     → targetRegionMb ∈ {1.1 (FoundationOne 315-gene panel), 10, 30 (≈WES)} — the
    //     real denominator.
    //   • mutationCount → number of somatic calls ∈ {3, 11, 150} — the real numerator.
    //   • includeSilent → whether non-counted "silent" calls (Germline + NotDetected) pad the call set.
    //     There is no includeSilent parameter; the business property is that such calls must NOT inflate
    //     TMB — so TMB is INVARIANT to this axis (the silent calls are correctly excluded).
    // Grid = 3 × 3 × 2 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: TMB and the High/Low call are a JOINT function of numerator and
    // denominator (a small panel with the same count flips Low→High), while the includeSilent axis must
    // leave the result unchanged — a genuine three-way interaction with a built-in invariance.
    // ═══════════════════════════════════════════════════════════════════════

    private const int SilentPadGermline = 7;     // non-somatic calls that must not be counted
    private const int SilentPadNotDetected = 5;

    /// <summary>
    /// For every (panelSize, somatic count, includeSilent) the TMB equals somaticCount / targetRegionMb
    /// and the High/Low classification matches the FDA ≥10 cutoff — and both are invariant to padding the
    /// call set with non-somatic "silent" calls, which must be excluded from the count.
    /// </summary>
    [Test, Combinatorial]
    public void CalculateTMB_PanelCountSilentGrid_CountsOnlySomaticPerMegabase(
        [Values(1.1, 10.0, 30.0)] double targetRegionMb,
        [Values(3, 11, 150)] int somaticCount,
        [Values(false, true)] bool includeSilent)
    {
        var calls = new List<OncologyAnalyzer.SomaticCall>();
        for (int i = 0; i < somaticCount; i++)
            calls.Add(MakeSomaticCall(OncologyAnalyzer.SomaticStatus.Somatic));
        if (includeSilent)
        {
            for (int i = 0; i < SilentPadGermline; i++)
                calls.Add(MakeSomaticCall(OncologyAnalyzer.SomaticStatus.Germline));
            for (int i = 0; i < SilentPadNotDetected; i++)
                calls.Add(MakeSomaticCall(OncologyAnalyzer.SomaticStatus.NotDetected));
        }

        // Independent ground truth (counts only somatic, divides by Mb).
        double expectedTmb = somaticCount / targetRegionMb;
        var expectedStatus = expectedTmb >= 10.0 ? OncologyAnalyzer.TmbStatus.High : OncologyAnalyzer.TmbStatus.Low;

        double tmb = OncologyAnalyzer.CalculateTMB(calls, targetRegionMb);

        tmb.Should().BeApproximately(expectedTmb, 1e-12, "TMB = somatic count / Mb, silent calls excluded");
        OncologyAnalyzer.ClassifyTMB(tmb).Should().Be(expectedStatus, "High ⟺ TMB ≥ 10 mut/Mb (FDA, inclusive)");
    }

    /// <summary>
    /// Interaction witness (panel × count flips the call): 11 somatic mutations are TMB-High on the 1.1 Mb
    /// FoundationOne panel (exactly 10.0 mut/Mb) but TMB-Low across a 30 Mb exome (≈0.37 mut/Mb). The
    /// classification flips purely on the panel-size axis at fixed count — the grid needs the combination.
    /// </summary>
    [Test]
    public void CalculateTMB_PanelSizeAxis_FlipsHighLowAtFixedCount()
    {
        OncologyAnalyzer.ClassifyTMB(OncologyAnalyzer.CalculateTMB(11, 1.1))
            .Should().Be(OncologyAnalyzer.TmbStatus.High, "11 / 1.1 = 10.0 ≥ 10 → High");
        OncologyAnalyzer.ClassifyTMB(OncologyAnalyzer.CalculateTMB(11, 30.0))
            .Should().Be(OncologyAnalyzer.TmbStatus.Low, "11 / 30 ≈ 0.37 < 10 → Low");
    }

    /// <summary>
    /// Interaction witness (includeSilent invariance): adding germline / not-detected calls to a fixed
    /// somatic set leaves TMB unchanged — silent variants do not inflate the burden (Chalmers et al. 2017
    /// count only somatic mutations).
    /// </summary>
    [Test]
    public void CalculateTMB_SilentCalls_DoNotInflateBurden()
    {
        var somaticOnly = Enumerable.Range(0, 11)
            .Select(_ => MakeSomaticCall(OncologyAnalyzer.SomaticStatus.Somatic)).ToList();
        var padded = new List<OncologyAnalyzer.SomaticCall>(somaticOnly);
        padded.AddRange(Enumerable.Range(0, 50).Select(_ => MakeSomaticCall(OncologyAnalyzer.SomaticStatus.Germline)));
        padded.AddRange(Enumerable.Range(0, 50).Select(_ => MakeSomaticCall(OncologyAnalyzer.SomaticStatus.NotDetected)));

        OncologyAnalyzer.CalculateTMB(padded, 1.1)
            .Should().Be(OncologyAnalyzer.CalculateTMB(somaticOnly, 1.1), "only somatic calls count toward TMB");
    }

    /// <summary>
    /// Boundary witness (FDA inclusive cutoff): classification flips only at exactly 10 mut/Mb — 9.9 is
    /// Low, 10.0 is High. Source: Marcus et al. (2021).
    /// </summary>
    [Test]
    public void ClassifyTMB_FdaCutoff_IsInclusiveAtTen()
    {
        OncologyAnalyzer.ClassifyTMB(9.9).Should().Be(OncologyAnalyzer.TmbStatus.Low);
        OncologyAnalyzer.ClassifyTMB(10.0).Should().Be(OncologyAnalyzer.TmbStatus.High);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-MSI-001 — Microsatellite instability detection (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 93.
    // Spec: tests/TestSpecs/ONCO-MSI-001.md (OncologyAnalyzer.DetectMSI / CalculateMSIScore / ClassifyMSIStatus).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Niu et al. (2014) MSIsensor, Bioinformatics 30(7):1015 (MSI score = unstable sites / valid
    // sites); niu-lab/msisensor2 README (MSI-High ⟺ score ≥ 20%, inclusive); Boland et al. (1998).
    //
    // MSI score = unstableLoci / totalLoci;  ClassifyMSIStatus → MSI-High ⟺ score ≥ 0.20.
    //
    // Checklist axes nLoci(3) × instabilityThreshold(3) map onto the real knobs:
    //   • nLoci → totalLoci (panel size) ∈ {10, 20, 50} valid evaluated loci.
    //   • instabilityThreshold → the UNSTABLE FRACTION relative to the fixed 20% MSIsensor2 cutoff:
    //     {BelowCutoff 0.10, AtCutoff 0.20 (inclusive boundary), AboveCutoff 0.40}. The 20% cutoff is a
    //     fixed constant, so this axis is realised as where the sample's unstable fraction sits w.r.t. it.
    // Grid = 3 × 3 = 9 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: the MSI-High call is a function of the FRACTION (it flips MSS↔MSI-High at
    // exactly 0.20) and is INVARIANT to the panel size — the same fraction gives the same status whether
    // measured over 10 or 50 loci. The grid exercises both the flip and the invariance jointly.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>Where the sample's unstable fraction sits relative to the fixed 20% MSI-High cutoff.</summary>
    public enum InstabilityLevel
    {
        /// <summary>Unstable fraction 0.10 &lt; 0.20 → MSS.</summary>
        BelowCutoff,

        /// <summary>Unstable fraction exactly 0.20 ≥ 0.20 → MSI-High (inclusive boundary).</summary>
        AtCutoff,

        /// <summary>Unstable fraction 0.40 &gt; 0.20 → MSI-High.</summary>
        AboveCutoff,
    }

    private static double InstabilityFraction(InstabilityLevel level) => level switch
    {
        InstabilityLevel.BelowCutoff => 0.10,
        InstabilityLevel.AtCutoff => 0.20,
        InstabilityLevel.AboveCutoff => 0.40,
        _ => throw new ArgumentOutOfRangeException(nameof(level)),
    };

    /// <summary>
    /// For every (panel size, instability level) the MSI score equals unstable/total and the MSI status
    /// matches the 20% cutoff — and the status depends only on the fraction, not the panel size (the same
    /// fraction yields the same MSI-High/MSS call at 10, 20 or 50 loci).
    /// </summary>
    [Test, Combinatorial]
    public void DetectMSI_LociFractionGrid_MatchesFractionAndCutoff(
        [Values(10, 20, 50)] int totalLoci,
        [Values] InstabilityLevel instability)
    {
        int unstable = (int)Math.Round(InstabilityFraction(instability) * totalLoci);
        var flags = BuildLocusFlags(totalLoci, unstable);

        // Independent ground truth.
        double expectedScore = (double)unstable / totalLoci;
        var expectedStatus = expectedScore >= 0.20 ? OncologyAnalyzer.MsiStatus.MSI_High : OncologyAnalyzer.MsiStatus.MSS;

        var result = OncologyAnalyzer.DetectMSI(flags);

        result.UnstableLoci.Should().Be(unstable);
        result.TotalLoci.Should().Be(totalLoci);
        result.Score.Should().BeApproximately(expectedScore, 1e-12, "MSI score = unstable / valid loci");
        result.Status.Should().Be(expectedStatus, "MSI-High ⟺ score ≥ 0.20 (MSIsensor2, inclusive)");
        OncologyAnalyzer.CalculateMSIScore(unstable, totalLoci).Should().Be(result.Score, "[INV-5] DetectMSI composes CalculateMSIScore");
    }

    /// <summary>
    /// Interaction witness (instability axis, inclusive boundary): the MSI call flips MSS→MSI-High exactly
    /// at the 20% cutoff — 4/25 = 16% is MSS, 5/25 = 20% is MSI-High. Source: niu-lab/msisensor2 (≥20%).
    /// </summary>
    [Test]
    public void ClassifyMSIStatus_TwentyPercentCutoff_IsInclusive()
    {
        OncologyAnalyzer.ClassifyMSIStatus(OncologyAnalyzer.CalculateMSIScore(4, 25))
            .Should().Be(OncologyAnalyzer.MsiStatus.MSS, "16% < 20% → MSS");
        OncologyAnalyzer.ClassifyMSIStatus(OncologyAnalyzer.CalculateMSIScore(5, 25))
            .Should().Be(OncologyAnalyzer.MsiStatus.MSI_High, "20% ≥ 20% → MSI-High");
    }

    /// <summary>
    /// Interaction witness (panel-size invariance): the same 20% unstable fraction gives an identical
    /// MSI-High call and identical score whether measured over 10 loci (2/10) or 50 loci (10/50) — the
    /// status depends on the fraction, not the absolute locus count.
    /// </summary>
    [Test]
    public void DetectMSI_SameFractionDifferentPanelSize_SameStatus()
    {
        var small = OncologyAnalyzer.DetectMSI(BuildLocusFlags(10, 2));
        var large = OncologyAnalyzer.DetectMSI(BuildLocusFlags(50, 10));

        small.Score.Should().Be(large.Score, "2/10 = 10/50 = 0.20");
        small.Status.Should().Be(large.Status).And.Be(OncologyAnalyzer.MsiStatus.MSI_High);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-HRD-001 — Homologous-recombination-deficiency genomic-scar score (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 94.
    // Spec: tests/TestSpecs/ONCO-HRD-001.md (OncologyAnalyzer.CalculateHRDScore / ClassifyHRDStatus / DetectHRD).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Telli et al. (2016) Clin Cancer Res 22:3764 (HRD = unweighted sum LOH+TAI+LST; HRD-High ⟺
    // score ≥ 42, inclusive); Birkbak et al. (2012) (TAI); Popova et al. (2012) (LST); Stewart et al. (2022).
    //
    // HRD score = LOH + TAI + LST;  ClassifyHRDStatus → HrdHigh ⟺ score ≥ 42.
    //
    // Checklist axes LOH(3) × TAI(3) × LST(3) map DIRECTLY onto the three real integer component inputs
    // (no deviation): LOH ∈ {5,14,25}, TAI ∈ {4,14,20}, LST ∈ {3,14,18}. Grid = 3³ = 27 = the checklist's
    // "Full Combos" for this row.
    //
    // The combinatorial point: HrdHigh depends on the SUM of all three scars — no single component decides
    // it. Each cell checks the unweighted sum and the inclusive 42 cutoff; the end-to-end DetectHRD result
    // preserves the components and composes the sum and classification.
    // ═══════════════════════════════════════════════════════════════════════

    private const int HrdHighCutoff = 42;

    /// <summary>
    /// For every (LOH, TAI, LST) the HRD score equals the unweighted sum and the status matches the
    /// inclusive 42 cutoff; DetectHRD preserves the components and composes the same score/status.
    /// </summary>
    [Test, Combinatorial]
    public void CalculateHRDScore_LohTaiLstGrid_UnweightedSumAndCutoff(
        [Values(5, 14, 25)] int loh,
        [Values(4, 14, 20)] int tai,
        [Values(3, 14, 18)] int lst)
    {
        int expectedScore = loh + tai + lst;
        var expectedStatus = expectedScore >= HrdHighCutoff ? OncologyAnalyzer.HrdStatus.HrdHigh : OncologyAnalyzer.HrdStatus.HrdNegative;

        OncologyAnalyzer.CalculateHRDScore(loh, tai, lst).Should().Be(expectedScore, "[INV-1] HRD = LOH + TAI + LST (unweighted)");
        OncologyAnalyzer.ClassifyHRDStatus(expectedScore).Should().Be(expectedStatus, "[INV-3] HrdHigh ⟺ score ≥ 42");

        var result = OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(loh, tai, lst));
        result.Score.Should().Be(expectedScore);
        result.Status.Should().Be(expectedStatus);
        result.Components.Should().Be(new OncologyAnalyzer.HrdComponents(loh, tai, lst), "[INV-4] DetectHRD preserves components");
    }

    /// <summary>
    /// Interaction witness: starting from a sub-threshold sum of 41, bumping ANY single scar component by
    /// one to reach 42 flips HRD-negative → HRD-high — no component is privileged; the call is on the sum.
    /// Source: Telli et al. (2016) unweighted sum, inclusive ≥42 cutoff.
    /// </summary>
    [Test]
    public void DetectHRD_AnyComponentCrossingFortyTwo_FlipsToHrdHigh()
    {
        OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(13, 14, 14)).Status
            .Should().Be(OncologyAnalyzer.HrdStatus.HrdNegative, "13+14+14 = 41 < 42");
        OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(14, 14, 14)).Status
            .Should().Be(OncologyAnalyzer.HrdStatus.HrdHigh, "LOH axis: 14+14+14 = 42");
        OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(13, 15, 14)).Status
            .Should().Be(OncologyAnalyzer.HrdStatus.HrdHigh, "TAI axis: 13+15+14 = 42");
        OncologyAnalyzer.DetectHRD(new OncologyAnalyzer.HrdComponents(13, 14, 15)).Status
            .Should().Be(OncologyAnalyzer.HrdStatus.HrdHigh, "LST axis: 13+14+15 = 42");
    }

    /// <summary>
    /// Interaction witness (INV-2 commutativity): the unweighted sum is order-independent — every
    /// permutation of (20, 15, 12) gives 47. Source: Telli et al. (2016) "unweighted sum".
    /// </summary>
    [Test]
    public void CalculateHRDScore_ComponentOrderPermuted_SameSum()
    {
        OncologyAnalyzer.CalculateHRDScore(20, 15, 12).Should().Be(47);
        OncologyAnalyzer.CalculateHRDScore(12, 20, 15).Should().Be(47);
        OncologyAnalyzer.CalculateHRDScore(15, 12, 20).Should().Be(47);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-SIG-002 — Mutational signature fitting / refitting (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 97.
    // Spec: tests/TestSpecs/ONCO-SIG-002.md (OncologyAnalyzer.FitSignatures / CosineSimilarity / ReconstructCatalog).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Blokzijl et al. (2018) MutationalPatterns (NNLS refit min‖Sx−d‖², x≥0; cosine quality);
    // Rosenthal et al. (2016) deconstructSigs (R = S·W, normalised weights); Lawson & Hanson (1974) NNLS.
    //
    // FitSignatures solves x = argmin‖S·x − d‖² s.t. x ≥ 0, then reports normalised exposures, the
    // reconstruction S·x, and cos(d, S·x).
    //
    // Checklist axes nSignatures(3) × nMutations(3) × solver(2). There is a single NNLS solver
    // (Lawson-Hanson); per the campaign convention we map the axes onto the real knobs and document it:
    //   • nSignatures → number of reference signatures k ∈ {1,2,3} (here built with DISJOINT channel
    //     support, so the NNLS reduces to an exact orthogonal projection with closed-form exposures).
    //   • nMutations  → the per-signature intensity v ∈ {50,500,2500} (∝ total mutation count).
    //   • solver      → the fit REGIME the single solver lands in: {ExactFit — catalog is an exact non-
    //     negative mix → reconstruction cosine = 1; ResidualBearing — catalog carries an unmodeled
    //     orthogonal residual block → exposures unchanged but cosine = √(k/(k+1)) < 1}. The active-set
    //     clamp branch (negative unconstrained coefficient) is covered by a dedicated witness below.
    // Grid = 3 × 3 × 2 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: the NNLS projection recovers the true per-signature exposures (each = v,
    // INV-4/INV-7) regardless of intensity or an orthogonal residual, while the reconstruction cosine is a
    // JOINT function of k and the fit regime (1 for an exact mix, √(k/(k+1)) when a residual is present).
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>Which fit regime the single NNLS solver lands in for the constructed catalog.</summary>
    public enum FitRegime
    {
        /// <summary>Catalog is an exact non-negative mix of the signatures → reconstruction cosine = 1.</summary>
        ExactFit,

        /// <summary>Catalog carries an unmodeled orthogonal residual block → cosine = √(k/(k+1)) &lt; 1.</summary>
        ResidualBearing,
    }

    /// <summary>
    /// For every (k signatures, intensity v, fit regime) the NNLS projection recovers exposure v for each
    /// disjoint-support signature (all ≥ 0, normalised to 1/k summing to 1), and the reconstruction cosine
    /// is 1 for an exact mix or √(k/(k+1)) when an unmodeled orthogonal residual is present — the exposures
    /// being unchanged by the residual (the projection ignores the unmodeled block).
    /// </summary>
    [Test, Combinatorial]
    public void FitSignatures_SignaturesIntensityRegimeGrid_ProjectsAndScoresReconstruction(
        [Values(1, 2, 3)] int nSignatures,
        [Values(50.0, 500.0, 2500.0)] double intensity,
        [Values] FitRegime regime)
    {
        var (catalog, signatures) = BuildOrthogonalFit(nSignatures, intensity, regime == FitRegime.ResidualBearing);

        double expectedCosine = regime == FitRegime.ExactFit ? 1.0 : Math.Sqrt((double)nSignatures / (nSignatures + 1));

        var result = OncologyAnalyzer.FitSignatures(catalog, signatures);

        result.Exposures.Should().HaveCount(nSignatures);
        result.Exposures.Should().OnlyContain(e => e >= 0.0, "[INV-4] NNLS exposures are non-negative");
        result.Exposures.Should().OnlyContain(e => Math.Abs(e - intensity) < 1e-9,
            "[INV-7] orthogonal projection recovers exposure v for each signature, residual ignored");
        result.NormalizedExposures.Sum().Should().BeApproximately(1.0, 1e-9, "[INV-6] proportions sum to 1");
        result.NormalizedExposures.Should().OnlyContain(p => Math.Abs(p - 1.0 / nSignatures) < 1e-9, "equal exposures → equal proportions");
        result.ReconstructionCosineSimilarity.Should().BeApproximately(expectedCosine, 1e-9,
            "reconstruction cosine is 1 for an exact mix, √(k/(k+1)) with an unmodeled residual");
    }

    /// <summary>
    /// Interaction witness (regime axis flips reconstruction quality, not exposures): for fixed k=2 and
    /// intensity, the exact-mix catalog reconstructs perfectly (cosine 1) while the residual-bearing
    /// catalog scores √(2/3) ≈ 0.8165 — yet both recover the same exposures. The cosine flips purely on the
    /// regime axis; the projection is invariant to the unmodeled residual.
    /// </summary>
    [Test]
    public void FitSignatures_RegimeAxis_FlipsCosineNotExposures()
    {
        var (exactCatalog, sigs) = BuildOrthogonalFit(2, 500.0, residual: false);
        var (residualCatalog, _) = BuildOrthogonalFit(2, 500.0, residual: true);

        var exact = OncologyAnalyzer.FitSignatures(exactCatalog, sigs);
        var residual = OncologyAnalyzer.FitSignatures(residualCatalog, sigs);

        exact.ReconstructionCosineSimilarity.Should().BeApproximately(1.0, 1e-9);
        residual.ReconstructionCosineSimilarity.Should().BeApproximately(Math.Sqrt(2.0 / 3.0), 1e-9);
        residual.Exposures.Should().BeEquivalentTo(exact.Exposures, "the unmodeled residual does not change the projection");
    }

    /// <summary>
    /// Witness for the active-set clamp branch (solver constraint binds): S columns [1,0] and [1,1] with
    /// d = [0,1] have a negative unconstrained coefficient for the first signature, so NNLS clamps it to 0
    /// and refits the second to 0.5 → exposures [0, 0.5]. Source: Lawson & Hanson (1974) active-set NNLS;
    /// spec M6.
    /// </summary>
    [Test]
    public void FitSignatures_NegativeUnconstrainedCoefficient_ClampsToZero()
    {
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1.0, 0.0 }, // signature 0
            new double[] { 1.0, 1.0 }, // signature 1
        };
        var catalog = new double[] { 0.0, 1.0 };

        var result = OncologyAnalyzer.FitSignatures(catalog, signatures);

        result.Exposures[0].Should().BeApproximately(0.0, 1e-9, "negative LS coefficient clamps to 0 (NNLS x ≥ 0)");
        result.Exposures[1].Should().BeApproximately(0.5, 1e-9, "second signature refits to 0.5");
    }

    /// <summary>
    /// Witness (cosine formula, spec M1-M4): identical vectors → 1, orthogonal → 0, and positive scaling is
    /// invariant. Source: Blokzijl et al. (2018) cosine similarity.
    /// </summary>
    [Test]
    public void CosineSimilarity_IdentityOrthogonalScale_MatchFormula()
    {
        OncologyAnalyzer.CosineSimilarity(new double[] { 1, 2, 3 }, new double[] { 1, 2, 3 }).Should().BeApproximately(1.0, 1e-12);
        OncologyAnalyzer.CosineSimilarity(new double[] { 1, 0 }, new double[] { 0, 1 }).Should().BeApproximately(0.0, 1e-12);
        OncologyAnalyzer.CosineSimilarity(new double[] { 3, 4 }, new double[] { 6, 8 }).Should().BeApproximately(1.0, 1e-12);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-SIG-003 — Bootstrap confidence intervals on signature exposures (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 98.
    // Spec: tests/TestSpecs/ONCO-SIG-003.md (OncologyAnalyzer.BootstrapExposures).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Senkin (2021) MSA (multinomial resample of the count catalog + NNLS refit + [2.5%,97.5%]
    // percentile CI); Huang et al. (2018); Efron (1979) percentile method; Hyndman & Fan (1996) type-7.
    //
    // BootstrapExposures resamples the catalog as a multinomial draw of N = Σcatalog, refits each draw by
    // NNLS, and reports a per-signature percentile CI; the point estimate is the NNLS fit of the observed
    // catalog. All three checklist axes are REAL parameters: nBootstrap → replicates, seed → seed,
    // nMutations → catalog total N.
    //
    // The grid uses a DEGENERATE catalog with a single non-zero channel: the multinomial draw is then
    // deterministic (every resample reproduces the observed catalog), so the bootstrap distribution
    // collapses to the point estimate exactly — INDEPENDENT of replicates and seed (Senkin 2021 corner
    // case). This makes every cell exactly verifiable: the active signature's interval is [N,N,N,N], the
    // absent signature's is [0,0,0,0].
    //   • nBootstrap → replicates ∈ {1, 10, 100}  • seed → ∈ {7, 42}  • nMutations → N ∈ {0, 100, 1000}.
    // Grid = 3 × 2 × 3 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: under the degenerate resample the CI collapses to the point estimate for
    // ALL (replicates, seed) — an invariance across two axes — while the nMutations axis scales the point
    // estimate (N). Seed-sensitivity on a NON-degenerate catalog (where the axis is genuinely live) and the
    // confidence-width / determinism properties are covered by witnesses.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// With a single-non-zero-channel catalog the multinomial resample is deterministic, so for every
    /// (replicates, seed, N) the bootstrap CI collapses exactly to the point estimate: the active signature
    /// reports [N,N,N,N] and the absent signature [0,0,0,0]. Replicates and seed are provably inert here;
    /// N scales the point estimate.
    /// </summary>
    [Test, Combinatorial]
    public void BootstrapExposures_ReplicatesSeedMutationsGrid_DegenerateResampleCollapsesToPointEstimate(
        [Values(1, 10, 100)] int replicates,
        [Values(7, 42)] int seed,
        [Values(0, 100, 1000)] int totalMutations)
    {
        var catalog = new[] { totalMutations, 0 }; // all mass on channel 0 → deterministic resample
        var signatures = new IReadOnlyList<double>[]
        {
            new double[] { 1.0, 0.0 }, // active signature (channel 0)
            new double[] { 0.0, 1.0 }, // absent signature (channel 1)
        };

        var intervals = OncologyAnalyzer.BootstrapExposures(catalog, signatures, replicates, 0.95, seed);

        intervals.Should().HaveCount(2, "[INV-5] one interval per signature, in order");

        var active = intervals[0];
        active.PointEstimate.Should().BeApproximately(totalMutations, 1e-9, "[INV-5] point estimate = NNLS exposure of observed catalog = N");
        active.Mean.Should().BeApproximately(totalMutations, 1e-9, "degenerate resample → every replicate = N");
        active.Lower.Should().BeApproximately(totalMutations, 1e-9);
        active.Upper.Should().BeApproximately(totalMutations, 1e-9);

        var absent = intervals[1];
        absent.PointEstimate.Should().Be(0.0);
        absent.Mean.Should().Be(0.0);
        absent.Lower.Should().Be(0.0);
        absent.Upper.Should().Be(0.0);

        foreach (var ci in intervals)
        {
            ci.Lower.Should().BeGreaterThanOrEqualTo(0.0, "[INV-1] bounds ≥ 0");
            ci.Lower.Should().BeLessThanOrEqualTo(ci.Upper, "[INV-2] lower ≤ upper");
            ci.Mean.Should().BeInRange(ci.Lower, ci.Upper, "[INV-3] lower ≤ mean ≤ upper");
        }
    }

    /// <summary>
    /// Interaction witness (determinism, INV-4): the same (catalog, signatures, replicates, confidence,
    /// seed) yields an element-wise identical result on a NON-degenerate catalog where the resample is
    /// genuinely random. Source: fixed RNG seed (Senkin 2021).
    /// </summary>
    [Test]
    public void BootstrapExposures_SameSeed_IsDeterministic()
    {
        var catalog = new[] { 13, 7 };
        var signatures = new IReadOnlyList<double>[] { new double[] { 1.0, 0.0 }, new double[] { 0.0, 1.0 } };

        var first = OncologyAnalyzer.BootstrapExposures(catalog, signatures, 200, 0.95, seed: 42);
        var second = OncologyAnalyzer.BootstrapExposures(catalog, signatures, 200, 0.95, seed: 42);

        second.Should().BeEquivalentTo(first, "identical inputs + seed → identical bootstrap result");
    }

    /// <summary>
    /// Interaction witness (seed axis is live on a non-degenerate catalog): two different seeds produce
    /// different bootstrap distributions — the seed genuinely matters when the multinomial draw is not
    /// degenerate (unlike the collapsed grid cells where it is provably inert).
    /// </summary>
    [Test]
    public void BootstrapExposures_DifferentSeeds_ProduceDifferentDistributions()
    {
        var catalog = new[] { 13, 7 };
        var signatures = new IReadOnlyList<double>[] { new double[] { 1.0, 0.0 }, new double[] { 0.0, 1.0 } };

        var seedA = OncologyAnalyzer.BootstrapExposures(catalog, signatures, 200, 0.95, seed: 42);
        var seedB = OncologyAnalyzer.BootstrapExposures(catalog, signatures, 200, 0.95, seed: 7);

        seedB.Should().NotBeEquivalentTo(seedA, "different seeds explore different resamples");
    }

    /// <summary>
    /// Interaction witness (confidence-width monotonicity, S1): a wider confidence level gives a
    /// wider-or-equal percentile interval — upper grows (or holds) and lower shrinks (or holds). Source:
    /// Efron (1979) percentile method.
    /// </summary>
    [Test]
    public void BootstrapExposures_WiderConfidence_GivesWiderInterval()
    {
        var catalog = new[] { 13, 7 };
        var signatures = new IReadOnlyList<double>[] { new double[] { 1.0, 0.0 }, new double[] { 0.0, 1.0 } };

        var narrow = OncologyAnalyzer.BootstrapExposures(catalog, signatures, 500, 0.50, seed: 42)[0];
        var wide = OncologyAnalyzer.BootstrapExposures(catalog, signatures, 500, 0.99, seed: 42)[0];

        wide.Upper.Should().BeGreaterThanOrEqualTo(narrow.Upper, "wider confidence → upper bound does not decrease");
        wide.Lower.Should().BeLessThanOrEqualTo(narrow.Lower, "wider confidence → lower bound does not increase");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-FUSION-001 — Gene-fusion detection from supporting reads (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 100.
    // Spec: tests/TestSpecs/ONCO-FUSION-001.md (OncologyAnalyzer.DetectFusions / IsInFrame / ComputeTotalSupport).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Haas et al. (2019); STAR-Fusion (min_junction_reads=1, min_sum_frags=2,
    // min_spanning_frags_only=5); Uhrig et al. (2021) Arriba (split_reads1/2 + discordant_mates).
    //
    // A candidate (distinct genes) is reported iff: junction (split) reads ≥ 1 AND total support ≥ 2; OR,
    // when there are no junction reads, discordant mates ≥ min_spanning_frags_only. Total support =
    // split5p + split3p + discordant.
    //
    // Checklist axes splitReads(3) × spanningReads(3) × minMapQ(2) map onto the real knobs:
    //   • splitReads   → junction (split) read count ∈ {0, 1, 3}.
    //   • spanningReads → discordant (spanning) mate count ∈ {0, 4, 5}.
    //   • minMapQ      → there is no read-level mapping-quality parameter; the analogous binary stringency
    //     knob is the threshold set: {Default (min_spanning_frags_only=5) vs Relaxed (=3)}. Relaxing it
    //     moves the spanning-only floor, so a 4-discordant spanning-only candidate flips reject→detect.
    // Grid = 3 × 3 × 2 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: detection is a JOINT function of junction count, spanning count and the
    // stringency — the spanning-only rule applies only when junction = 0, total-support gates the junction
    // path, and the stringency axis moves the spanning floor. Each cell checks the rule and total-support
    // (INV-2) re-derived from the counts.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>Detection stringency (the minMapQ-proxy threshold axis).</summary>
    public enum FusionStringency
    {
        /// <summary>STAR-Fusion defaults (min_spanning_frags_only = 5).</summary>
        Default,

        /// <summary>Relaxed spanning floor (min_spanning_frags_only = 3).</summary>
        Relaxed,
    }

    /// <summary>
    /// For every (junction reads, discordant mates, stringency) the detection decision and the reported
    /// total support match the STAR-Fusion min-support rule re-derived from the counts.
    /// </summary>
    [Test, Combinatorial]
    public void DetectFusions_SplitSpanningStringencyGrid_MatchesMinSupportRule(
        [Values(0, 1, 3)] int splitReads,
        [Values(0, 4, 5)] int discordantMates,
        [Values] FusionStringency stringency)
    {
        var candidate = new OncologyAnalyzer.FusionCandidate("EML4", "ALK", splitReads, 0, discordantMates);
        OncologyAnalyzer.FusionDetectionThresholds? thresholds = stringency == FusionStringency.Relaxed
            ? new OncologyAnalyzer.FusionDetectionThresholds(MinSpanningFragsOnly: 3)
            : null; // STAR-Fusion defaults

        // Independent ground truth.
        int spanningOnlyFloor = stringency == FusionStringency.Relaxed ? 3 : 5;
        int junction = splitReads;
        int totalSupport = splitReads + discordantMates;
        bool expectedDetected = junction >= 1 ? totalSupport >= 2 : discordantMates >= spanningOnlyFloor;

        var calls = OncologyAnalyzer.DetectFusions(new[] { candidate }, thresholds);

        calls.Should().HaveCount(expectedDetected ? 1 : 0,
            "detection = junction-path (total ≥ 2) when junction ≥ 1, else spanning-only (discordant ≥ floor)");
        if (expectedDetected)
        {
            var call = calls[0];
            call.JunctionReads.Should().Be(junction);
            call.DiscordantMates.Should().Be(discordantMates);
            call.TotalSupport.Should().Be(totalSupport, "[INV-2] total = split5p + split3p + discordant");
        }
    }

    /// <summary>
    /// Interaction witness (stringency × spanning, junction=0): a spanning-only candidate with 4 discordant
    /// mates is rejected under the STAR-Fusion default floor (5) but detected once the floor is relaxed to
    /// 3 — the call flips purely on the stringency axis. Source: STAR-Fusion min_spanning_frags_only.
    /// </summary>
    [Test]
    public void DetectFusions_StringencyAxis_FlipsSpanningOnlyCandidate()
    {
        var candidate = new OncologyAnalyzer.FusionCandidate("CD74", "ROS1", 0, 0, 4);

        OncologyAnalyzer.DetectFusions(new[] { candidate }).Should().BeEmpty("4 < default spanning floor 5");
        OncologyAnalyzer.DetectFusions(new[] { candidate }, new OncologyAnalyzer.FusionDetectionThresholds(MinSpanningFragsOnly: 3))
            .Should().ContainSingle("4 ≥ relaxed spanning floor 3");
    }

    /// <summary>
    /// Interaction witness (INV-1 same-gene, INV-4 ordering): a gene is never fused with itself regardless
    /// of support, and reported fusions are ordered by descending total support. Source: Registry invariant;
    /// STAR-Fusion scoring by abundance.
    /// </summary>
    [Test]
    public void DetectFusions_RejectsSameGeneAndOrdersBySupport()
    {
        var candidates = new[]
        {
            new OncologyAnalyzer.FusionCandidate("ALK", "ALK", 20, 20, 20),   // same gene → rejected
            new OncologyAnalyzer.FusionCandidate("EML4", "ALK", 3, 2, 4),     // total 9
            new OncologyAnalyzer.FusionCandidate("TMPRSS2", "ERG", 1, 0, 1),  // total 2
        };

        var calls = OncologyAnalyzer.DetectFusions(candidates);

        calls.Should().HaveCount(2, "[INV-1] the ALK-ALK self-fusion is excluded");
        calls.Select(c => c.TotalSupport).Should().BeInDescendingOrder("[INV-4] ordered by descending support");
        calls[0].TotalSupport.Should().Be(9);
    }

    /// <summary>
    /// Witness (INV-5 reading frame): a junction is in-frame iff (fivePrimeCodingBases − threePrimeStartPhase)
    /// mod 3 == 0. Source: Genomics England exon-phase rule + reading-frame modulo-3.
    /// </summary>
    [Test]
    public void IsInFrame_CodonPhaseRule_MatchesModuloThree()
    {
        OncologyAnalyzer.IsInFrame(300, 0).Should().BeTrue("300 mod 3 = 0");
        OncologyAnalyzer.IsInFrame(301, 0).Should().BeFalse("301 mod 3 = 1");
        OncologyAnalyzer.IsInFrame(301, 1).Should().BeTrue("(301 − 1) mod 3 = 0");
        OncologyAnalyzer.IsInFrame(302, 0).Should().BeFalse("302 mod 3 = 2");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-CNA-001 — Copy-number alteration classification (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 103.
    // Spec: tests/TestSpecs/ONCO-CNA-001.md (OncologyAnalyzer.ClassifyCopyNumber / CallCopyNumber / Log2RatioToCopyNumber).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: CNVkit absolute_threshold (hard-threshold integer CN; default cutoffs −1.1,−0.25,0.2,0.7;
    // inclusive log2 ≤ thresh; amp = ceil(ploidy·2^log2)); Mermel et al. (2011) GISTIC2 amplitude semantics.
    //
    // ClassifyCopyNumber: integer CN = index of the first cutoff log2 ≤, else ceil(ploidy·2^log2); state
    // 0→DeepDeletion,1→Loss,2→Neutral,3→Gain,≥4→Amplification; absolute CN = ploidy·2^log2.
    //
    // Checklist axes log2Range(3) × binSize(3) × ploidy(2) map onto the real knobs:
    //   • log2Range → the log2 copy ratio ∈ {−0.3, 0.25, 1.0}, chosen NEAR cutoff boundaries so the bin
    //     partition matters.
    //   • binSize   → the cutoff set that defines the CN bins {Default (−1.1,−0.25,0.2,0.7), Germline
    //     (−1.1,−0.4,0.3,0.7), Narrow (−0.5,−0.1,0.1,0.5)} — the real thresholds parameter.
    //   • ploidy    → reference ploidy ∈ {2, 4} (the real parameter); scales the absolute CN and the
    //     amplification integer CN.
    // Grid = 3 × 3 × 2 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: the state is a JOINT function of log2 and the bin cutoffs (the same −0.3 is
    // Loss under default/narrow but Neutral under germline cutoffs; 0.25 is Gain vs Neutral), while ploidy
    // scales the absolute CN and the amplification integer CN (CN 4 at ploidy 2 vs 8 at ploidy 4). Each
    // cell is checked against the CNVkit rule re-derived independently from the inputs.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>The cutoff set defining the copy-number bins (the binSize-proxy axis).</summary>
    public enum CnThresholdSet
    {
        /// <summary>CNVkit tumor defaults (−1.1, −0.25, 0.2, 0.7).</summary>
        Default,

        /// <summary>Germline-tuned cutoffs (−1.1, −0.4, 0.3, 0.7).</summary>
        Germline,

        /// <summary>Narrow neutral band (−0.5, −0.1, 0.1, 0.5).</summary>
        Narrow,
    }

    private static double[] CnCutoffs(CnThresholdSet set) => set switch
    {
        CnThresholdSet.Default => new[] { -1.1, -0.25, 0.2, 0.7 },
        CnThresholdSet.Germline => new[] { -1.1, -0.4, 0.3, 0.7 },
        CnThresholdSet.Narrow => new[] { -0.5, -0.1, 0.1, 0.5 },
        _ => throw new ArgumentOutOfRangeException(nameof(set)),
    };

    /// <summary>
    /// For every (log2, cutoff set, ploidy) the integer CN, CNA state and absolute copy number match the
    /// CNVkit hard-threshold rule re-derived independently from the inputs.
    /// </summary>
    [Test, Combinatorial]
    public void ClassifyCopyNumber_Log2BinPloidyGrid_MatchesCnvkitThresholdRule(
        [Values(-0.3, 0.25, 1.0)] double log2Ratio,
        [Values] CnThresholdSet thresholdSet,
        [Values(2.0, 4.0)] double ploidy)
    {
        double[] cutoffs = CnCutoffs(thresholdSet);

        // Independent ground truth (CNVkit absolute_threshold).
        int expectedCn = ExpectedIntegerCopyNumber(log2Ratio, cutoffs, ploidy);
        var expectedState = ExpectedCopyNumberState(expectedCn);
        double expectedAbsolute = ploidy * Math.Pow(2.0, log2Ratio);

        var call = OncologyAnalyzer.ClassifyCopyNumber(log2Ratio, cutoffs, ploidy);

        call.IntegerCopyNumber.Should().Be(expectedCn, "CN = first cutoff index log2 ≤, else ceil(ploidy·2^log2)");
        call.State.Should().Be(expectedState, "[INV-4] state ↔ integer CN mapping");
        call.AbsoluteCopyNumber.Should().BeApproximately(expectedAbsolute, 1e-9, "absolute CN = ploidy·2^log2");
        call.Log2Ratio.Should().Be(log2Ratio);
    }

    /// <summary>
    /// Interaction witness (binSize axis flips state): the same log2 = 0.25 is a Gain under the default
    /// cutoffs (&gt; 0.2) but Neutral under the germline cutoffs (≤ 0.3) — the state flips purely on the
    /// cutoff-set axis. Source: CNVkit threshold bins.
    /// </summary>
    [Test]
    public void ClassifyCopyNumber_ThresholdSetAxis_FlipsGainVsNeutral()
    {
        OncologyAnalyzer.ClassifyCopyNumber(0.25, CnCutoffs(CnThresholdSet.Default)).State
            .Should().Be(OncologyAnalyzer.CopyNumberState.Gain, "0.25 > 0.2 default gain cutoff");
        OncologyAnalyzer.ClassifyCopyNumber(0.25, CnCutoffs(CnThresholdSet.Germline)).State
            .Should().Be(OncologyAnalyzer.CopyNumberState.Neutral, "0.25 ≤ 0.3 germline neutral cutoff");
    }

    /// <summary>
    /// Interaction witness (ploidy axis scales amplification CN): an amplified log2 = 1.0 is CN 4 at diploid
    /// reference but CN 8 at tetraploid — the amplification integer CN = ceil(ploidy·2^log2) tracks ploidy,
    /// while both remain the Amplification state. Source: CNVkit amp ceil.
    /// </summary>
    [Test]
    public void ClassifyCopyNumber_PloidyAxis_ScalesAmplificationCopyNumber()
    {
        var diploid = OncologyAnalyzer.ClassifyCopyNumber(1.0, ploidy: 2.0);
        var tetraploid = OncologyAnalyzer.ClassifyCopyNumber(1.0, ploidy: 4.0);

        diploid.IntegerCopyNumber.Should().Be(4, "ceil(2·2) = 4");
        tetraploid.IntegerCopyNumber.Should().Be(8, "ceil(4·2) = 8");
        diploid.State.Should().Be(OncologyAnalyzer.CopyNumberState.Amplification);
        tetraploid.State.Should().Be(OncologyAnalyzer.CopyNumberState.Amplification);
    }

    /// <summary>
    /// Witness (inclusive boundaries, spec M7-M10): a log2 exactly on a default cutoff takes the LOWER CN
    /// state of that bin (log2 ≤ thresh). Source: CNVkit binning loop.
    /// </summary>
    [Test]
    public void ClassifyCopyNumber_ExactCutoffs_AreInclusiveLowerState()
    {
        OncologyAnalyzer.ClassifyCopyNumber(-1.1).State.Should().Be(OncologyAnalyzer.CopyNumberState.DeepDeletion);
        OncologyAnalyzer.ClassifyCopyNumber(-0.25).State.Should().Be(OncologyAnalyzer.CopyNumberState.Loss);
        OncologyAnalyzer.ClassifyCopyNumber(0.2).State.Should().Be(OncologyAnalyzer.CopyNumberState.Neutral);
        OncologyAnalyzer.ClassifyCopyNumber(0.7).State.Should().Be(OncologyAnalyzer.CopyNumberState.Gain);
    }

    private static int ExpectedIntegerCopyNumber(double log2Ratio, double[] cutoffs, double ploidy)
    {
        if (double.IsNaN(log2Ratio))
            return (int)Math.Round(ploidy, MidpointRounding.AwayFromZero);
        for (int cn = 0; cn < cutoffs.Length; cn++)
        {
            if (log2Ratio <= cutoffs[cn])
                return cn;
        }
        return (int)Math.Ceiling(ploidy * Math.Pow(2.0, log2Ratio));
    }

    private static OncologyAnalyzer.CopyNumberState ExpectedCopyNumberState(int cn) => cn switch
    {
        0 => OncologyAnalyzer.CopyNumberState.DeepDeletion,
        1 => OncologyAnalyzer.CopyNumberState.Loss,
        2 => OncologyAnalyzer.CopyNumberState.Neutral,
        3 => OncologyAnalyzer.CopyNumberState.Gain,
        _ => OncologyAnalyzer.CopyNumberState.Amplification,
    };

    // ───────────────────────────────────────────────────────────────────────
    // Helpers — engineered constructs + independent ground truth
    // ───────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Builds a profile of <paramref name="n"/> genes with strictly-distinct descending expression
    /// (G0 highest … G(n-1) lowest), designating <paramref name="nHits"/> of them as signature genes
    /// at the requested ranking <paramref name="placement"/>. Distinct values make the descending
    /// rank order unambiguous; nHits &lt; n guarantees at least one miss (required by ssGSEA).
    /// </summary>
    private static (Dictionary<string, double> profile, string[] hitGenes, double[] values, HashSet<int> hitIdx)
        BuildEnrichmentProfile(int n, int nHits, HitPlacement placement)
    {
        var values = new double[n];
        for (int idx = 0; idx < n; idx++)
            values[idx] = n - idx; // G0 = n (highest) … G(n-1) = 1 (lowest)

        var hitIdx = new HashSet<int>();
        switch (placement)
        {
            case HitPlacement.Top:
                for (int k = 0; k < nHits; k++) hitIdx.Add(k);
                break;
            case HitPlacement.Bottom:
                for (int k = 0; k < nHits; k++) hitIdx.Add(n - 1 - k);
                break;
            case HitPlacement.Spread:
                int stride = n / nHits;
                for (int k = 0; k < nHits; k++) hitIdx.Add(k * stride);
                break;
        }

        var profile = new Dictionary<string, double>();
        for (int idx = 0; idx < n; idx++)
            profile[$"G{idx}"] = values[idx];

        string[] hitGenes = hitIdx.Select(i => $"G{i}").ToArray();
        return (profile, hitGenes, values, hitIdx);
    }

    /// <summary>
    /// Independent ground-truth ssGSEA enrichment score computed directly from the mathematical
    /// definition: rank genes descending, accumulate a running sum that gains rank^τ / Σrank^τ on
    /// each hit and loses 1/N_miss on each miss, and integrate (sum) the running sum over all
    /// positions. Source: Barbie et al. (2009); Hänzelmann et al. (2013) GSVA integral form.
    /// </summary>
    private static double GroundTruthSsGseaIntegral(IReadOnlyList<double> values, HashSet<int> hitIdx)
    {
        int n = values.Count;
        var order = Enumerable.Range(0, n).OrderByDescending(i => values[i]).ToList();

        int nHits = hitIdx.Count;
        int nMiss = n - nHits;
        if (nHits == 0 || nMiss == 0) return 0.0;

        double totalHitWeight = 0.0;
        for (int pos = 0; pos < n; pos++)
        {
            if (hitIdx.Contains(order[pos]))
                totalHitWeight += Math.Pow(n - pos, SsGseaTau);
        }
        if (totalHitWeight == 0.0) return 0.0;

        double runningSum = 0.0, integral = 0.0, missStep = 1.0 / nMiss;
        for (int pos = 0; pos < n; pos++)
        {
            if (hitIdx.Contains(order[pos]))
                runningSum += Math.Pow(n - pos, SsGseaTau) / totalHitWeight;
            else
                runningSum -= missStep;
            integral += runningSum;
        }
        return integral;
    }

    /// <summary>
    /// Builds <paramref name="k"/> reference signatures with DISJOINT channel support (signature j occupies
    /// channels [2j, 2j+1]) and a catalog that is their exact non-negative mix with per-signature intensity
    /// <paramref name="v"/>. With orthogonal signatures the NNLS reduces to a closed-form projection giving
    /// exposure v for each. When <paramref name="residual"/> is set, an extra unmodeled block [2k, 2k+1] —
    /// covered by NO signature — is filled with v, so the reconstruction cannot represent it and the cosine
    /// drops to √(k/(k+1)) while the exposures are unchanged.
    /// </summary>
    private static (List<double> catalog, IReadOnlyList<double>[] signatures) BuildOrthogonalFit(int k, double v, bool residual)
    {
        int channelCount = 2 * (k + 1);

        var signatures = new IReadOnlyList<double>[k];
        for (int j = 0; j < k; j++)
        {
            var sig = new double[channelCount];
            sig[2 * j] = 1.0;
            sig[2 * j + 1] = 1.0;
            signatures[j] = sig;
        }

        var catalog = new List<double>(new double[channelCount]);
        for (int j = 0; j < k; j++)
        {
            catalog[2 * j] = v;
            catalog[2 * j + 1] = v;
        }
        if (residual)
        {
            catalog[2 * k] = v;
            catalog[2 * k + 1] = v;
        }

        return (catalog, signatures);
    }

    /// <summary>
    /// Builds <paramref name="totalLoci"/> per-locus stability flags with <paramref name="unstable"/> of them
    /// marked unstable (true) and the rest stable (false).
    /// </summary>
    private static List<bool> BuildLocusFlags(int totalLoci, int unstable)
    {
        var flags = new List<bool>(totalLoci);
        for (int i = 0; i < unstable; i++) flags.Add(true);
        for (int i = unstable; i < totalLoci; i++) flags.Add(false);
        return flags;
    }

    /// <summary>
    /// Builds a <see cref="OncologyAnalyzer.SomaticCall"/> carrying the requested status; only the status
    /// matters for the TMB count, so the variant and VAF fields are placeholders.
    /// </summary>
    private static OncologyAnalyzer.SomaticCall MakeSomaticCall(OncologyAnalyzer.SomaticStatus status)
    {
        var variant = new OncologyAnalyzer.VariantObservation("chr1", 1, "C", "T", 0, 0, 0, 0);
        return new OncologyAnalyzer.SomaticCall(variant, 0.0, 0.0, status, 0.0);
    }

    /// <summary>
    /// Builds an artifact observation for the requested substitution family (FFPE C&gt;T or OxoG G&gt;T),
    /// GIV read-orientation imbalance (R1/R2 ∈ {1.0, 1.5, 2.0}) and strand contingency table
    /// (Balanced/Mild/Strong). The three knobs populate disjoint fields of the record.
    /// </summary>
    private static OncologyAnalyzer.ArtifactObservation BuildArtifactObservation(
        ArtifactFamily family, GivLevel giv, StrandSegregation strand)
    {
        (char reference, char alternate) = family switch
        {
            ArtifactFamily.FfpeDeamination => ('C', 'T'), // C>T deamination
            ArtifactFamily.OxoG => ('G', 'T'),            // G>T oxidation
            _ => throw new ArgumentOutOfRangeException(nameof(family)),
        };

        (int r1, int r2) = giv switch
        {
            GivLevel.Undamaged => (100, 100), // 1.0
            GivLevel.Boundary => (150, 100),  // 1.5 (not > 1.5)
            GivLevel.Damaged => (200, 100),   // 2.0
            _ => throw new ArgumentOutOfRangeException(nameof(giv)),
        };

        (int rf, int rr, int af, int ar) = strand switch
        {
            StrandSegregation.Balanced => (20, 20, 20, 20),
            StrandSegregation.Mild => (20, 20, 16, 4),
            StrandSegregation.Strong => (30, 0, 0, 30),
            _ => throw new ArgumentOutOfRangeException(nameof(strand)),
        };

        return new OncologyAnalyzer.ArtifactObservation(reference, alternate, rf, rr, af, ar, r1, r2);
    }

    /// <summary>
    /// Builds a 10-mutation spectrum for the focal driver gene: <paramref name="recurrentCopies"/> copies
    /// of the query missense at the shared focal codon, an optional 5-mutation truncating (nonsense) block
    /// at distinct codons when <paramref name="tsgBackground"/> is set, and distinct-codon missense filler
    /// padding the total to a fixed 10 so the 20/20 fractions are exactly predictable. Filler codons never
    /// collide with the focal or truncating codons, so they add to the denominator only (not recurrent, not
    /// truncating).
    /// </summary>
    private static List<OncologyAnalyzer.GeneMutation> BuildDriverGeneSpectrum(int recurrentCopies, bool tsgBackground)
    {
        var mutations = new List<OncologyAnalyzer.GeneMutation>();

        for (int i = 0; i < recurrentCopies; i++)
            mutations.Add(new OncologyAnalyzer.GeneMutation(FocalGene, FocalPosition, OncologyAnalyzer.MutationConsequence.Missense));

        int truncating = tsgBackground ? 5 : 0;
        for (int i = 0; i < truncating; i++)
            mutations.Add(new OncologyAnalyzer.GeneMutation(FocalGene, 200 + i, OncologyAnalyzer.MutationConsequence.Nonsense));

        for (int p = 0; mutations.Count < DriverSpectrumTotal; p++)
            mutations.Add(new OncologyAnalyzer.GeneMutation(FocalGene, 300 + p, OncologyAnalyzer.MutationConsequence.Missense));

        return mutations;
    }

    /// <summary>
    /// Builds a bulk expression profile as the weighted sum of two cell-type signatures from the
    /// default LM22-style signature matrix: m = wₐ·Sₐ + w_b·S_b. By construction m lies in the
    /// column space of the signature matrix.
    /// </summary>
    private static Dictionary<string, double> BuildMixtureProfile(string cellTypeA, double weightA, string cellTypeB, double weightB)
    {
        var sigA = ImmuneAnalyzer.DefaultSignatureMatrix[cellTypeA];
        var sigB = ImmuneAnalyzer.DefaultSignatureMatrix[cellTypeB];

        var profile = new Dictionary<string, double>();
        foreach (var gene in sigA.Keys.Union(sigB.Keys))
        {
            double a = sigA.TryGetValue(gene, out double va) ? va : 0.0;
            double b = sigB.TryGetValue(gene, out double vb) ? vb : 0.0;
            profile[gene] = a * weightA + b * weightB;
        }
        return profile;
    }
}
