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
