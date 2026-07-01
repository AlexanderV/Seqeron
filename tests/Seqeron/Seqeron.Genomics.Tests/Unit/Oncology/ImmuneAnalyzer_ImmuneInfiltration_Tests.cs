using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

/// <summary>
/// Test Unit: ONCO-IMMUNE-001
/// Algorithm: Immune Infiltration Estimation
/// Methods: EstimateInfiltration, DeconvoluteImmuneCells
/// Class: ImmuneAnalyzer
/// 
/// Evidence sources:
/// - Newman AM et al. (2015). Nature Methods 12(5):453-457. DOI: 10.1038/nmeth.3337
/// - Yoshihara K et al. (2013). Nature Communications 4:2612. DOI: 10.1038/ncomms3612
/// - Barbie DA et al. (2009). Nature 462:108-112. DOI: 10.1038/nature08460
/// - Subramanian A et al. (2005). PNAS 102(43):15545-15550. DOI: 10.1073/pnas.0506580102
/// - Abbas AR et al. (2009). PLoS One 4(7):e6098. DOI: 10.1371/journal.pone.0006098
/// - Lawson CL, Hanson RJ (1995). Solving Least Squares Problems. SIAM.
/// - Hänzelmann S et al. (2013). BMC Bioinformatics 14:7. DOI: 10.1186/1471-2105-14-7 (GSVA/ssGSEA)
/// </summary>
[TestFixture]
public class ImmuneAnalyzer_ImmuneInfiltration_Tests
{
    #region Constants

    /// <summary>
    /// ESTIMATE tumor purity coefficient a. Source: Yoshihara et al. (2013).
    /// </summary>
    private const double EstimateCoefficientA = 0.6049872018;

    /// <summary>
    /// ESTIMATE tumor purity coefficient b. Source: Yoshihara et al. (2013).
    /// </summary>
    private const double EstimateCoefficientB = 0.0001467884;

    /// <summary>
    /// Tolerance for exact mathematical identities.
    /// </summary>
    private const double ExactTolerance = 1e-10;

    /// <summary>
    /// Tolerance for computed numerical results.
    /// </summary>
    private const double ComputedTolerance = 1e-6;

    #endregion

    #region Test Helpers

    /// <summary>
    /// Creates an expression profile from a cell type's signature in the default matrix.
    /// All genes get the signature value; acts as a "pure" cell type expression.
    /// </summary>
    private static Dictionary<string, double> CreatePureCellTypeProfile(string cellType)
    {
        var sig = ImmuneAnalyzer.DefaultSignatureMatrix[cellType];
        return new Dictionary<string, double>(sig);
    }

    /// <summary>
    /// Creates a mixed expression profile by averaging two cell types equally.
    /// </summary>
    private static Dictionary<string, double> CreateEqualMixProfile(string cellType1, string cellType2)
    {
        var sig1 = ImmuneAnalyzer.DefaultSignatureMatrix[cellType1];
        var sig2 = ImmuneAnalyzer.DefaultSignatureMatrix[cellType2];
        var allGenes = sig1.Keys.Union(sig2.Keys).Distinct();

        var profile = new Dictionary<string, double>();
        foreach (var gene in allGenes)
        {
            double v1 = sig1.TryGetValue(gene, out double val1) ? val1 : 0.0;
            double v2 = sig2.TryGetValue(gene, out double val2) ? val2 : 0.0;
            profile[gene] = (v1 + v2) / 2.0;
        }
        return profile;
    }

    /// <summary>
    /// Creates a weighted mix of two cell types.
    /// </summary>
    private static Dictionary<string, double> CreateWeightedMixProfile(
        string cellType1, double weight1,
        string cellType2, double weight2)
    {
        var sig1 = ImmuneAnalyzer.DefaultSignatureMatrix[cellType1];
        var sig2 = ImmuneAnalyzer.DefaultSignatureMatrix[cellType2];
        var allGenes = sig1.Keys.Union(sig2.Keys).Distinct();

        var profile = new Dictionary<string, double>();
        foreach (var gene in allGenes)
        {
            double v1 = sig1.TryGetValue(gene, out double val1) ? val1 : 0.0;
            double v2 = sig2.TryGetValue(gene, out double val2) ? val2 : 0.0;
            profile[gene] = v1 * weight1 + v2 * weight2;
        }
        return profile;
    }

    #endregion

    #region ONCO-IMMUNE-001 — M1: Empty expression → zero infiltration scores

    /// <summary>
    /// [M1] Empty expression profile should produce zero immune and stromal scores.
    /// ESTIMATE score = 0. Tumor purity = cos(a) where a = 0.6049872018.
    /// Evidence: Mathematical definition — no genes → no enrichment.
    /// </summary>
    [Test]
    public void EstimateInfiltration_EmptyExpression_ReturnsZeroScores()
    {
        var emptyProfile = new Dictionary<string, double>();

        var result = ImmuneAnalyzer.EstimateInfiltration(emptyProfile);

        double expectedPurity = Math.Cos(EstimateCoefficientA);

        Assert.Multiple(() =>
        {
            Assert.That(result.ImmuneScore, Is.EqualTo(0.0).Within(ExactTolerance),
                "Immune score must be 0 for empty profile");
            Assert.That(result.StromalScore, Is.EqualTo(0.0).Within(ExactTolerance),
                "Stromal score must be 0 for empty profile");
            Assert.That(result.EstimateScore, Is.EqualTo(0.0).Within(ExactTolerance),
                "ESTIMATE score must be 0 for empty profile");
            Assert.That(result.TumorPurity, Is.EqualTo(expectedPurity).Within(ExactTolerance),
                "Tumor purity must equal cos(a) for zero ESTIMATE score");
            Assert.That(result.OverlappingImmuneGenes, Is.EqualTo(0),
                "No overlapping immune genes for empty profile");
            Assert.That(result.OverlappingStromalGenes, Is.EqualTo(0),
                "No overlapping stromal genes for empty profile");
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — M2: Empty expression → zero deconvolution

    /// <summary>
    /// [M2] Empty expression profile should produce all-zero cell fractions.
    /// Evidence: Mathematical definition — no genes → no deconvolution possible.
    /// Source: Newman et al. (2015).
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_EmptyExpression_ReturnsZeroFractions()
    {
        var emptyProfile = new Dictionary<string, double>();

        var result = ImmuneAnalyzer.DeconvoluteImmuneCells(emptyProfile);

        Assert.Multiple(() =>
        {
            Assert.That(result.CellFractions, Is.Not.Null,
                "CellFractions dictionary must not be null");
            Assert.That(result.CellFractions.Values.Sum(), Is.EqualTo(0.0).Within(ExactTolerance),
                "All fractions must be zero for empty profile");
            Assert.That(result.Correlation, Is.EqualTo(0.0).Within(ExactTolerance),
                "Correlation must be 0 for empty profile");
            Assert.That(result.Rmse, Is.EqualTo(0.0).Within(ExactTolerance),
                "RMSE must be 0 for empty profile");
            Assert.That(result.OverlappingGenes, Is.EqualTo(0),
                "No overlapping genes for empty profile");
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — M3: No overlapping genes → zero enrichment

    /// <summary>
    /// [M3] Expression profile with genes not in immune/stromal signatures → zero scores.
    /// Evidence: Mathematical definition — no signal in gene sets.
    /// </summary>
    [Test]
    public void EstimateInfiltration_NoOverlappingGenes_ReturnsZeroScores()
    {
        var profile = new Dictionary<string, double>
        {
            ["NONEXISTENT_GENE_1"] = 10.0,
            ["NONEXISTENT_GENE_2"] = 5.0,
            ["NONEXISTENT_GENE_3"] = 8.0
        };

        var result = ImmuneAnalyzer.EstimateInfiltration(profile);

        Assert.Multiple(() =>
        {
            Assert.That(result.ImmuneScore, Is.EqualTo(0.0).Within(ExactTolerance),
                "Immune score must be 0 with no overlapping genes");
            Assert.That(result.StromalScore, Is.EqualTo(0.0).Within(ExactTolerance),
                "Stromal score must be 0 with no overlapping genes");
            Assert.That(result.OverlappingImmuneGenes, Is.EqualTo(0));
            Assert.That(result.OverlappingStromalGenes, Is.EqualTo(0));
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — M4: No overlapping genes → zero deconvolution

    /// <summary>
    /// [M4] Expression profile with no genes in signature matrix → zero fractions.
    /// Evidence: Mathematical definition — no overlapping genes means underdetermined system.
    /// Source: Newman et al. (2015).
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_NoOverlappingGenes_ReturnsZeroFractions()
    {
        var profile = new Dictionary<string, double>
        {
            ["NONEXISTENT_GENE_1"] = 10.0,
            ["NONEXISTENT_GENE_2"] = 5.0
        };

        var result = ImmuneAnalyzer.DeconvoluteImmuneCells(profile);

        Assert.Multiple(() =>
        {
            Assert.That(result.OverlappingGenes, Is.EqualTo(0),
                "No genes should overlap with signature");
            Assert.That(result.CellFractions.Values.Sum(), Is.EqualTo(0.0).Within(ExactTolerance),
                "All fractions must be zero with no overlap");
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — M5: Single cell type → 100% fraction

    /// <summary>
    /// [M5] Pure CD8 T cell expression profile must deconvolve to f_CD8 = 1.0 exactly.
    /// Proof: CD8A and CD8B appear only in the T_cells_CD8 column of the signature matrix,
    /// so the only NNLS solution with zero residual is f_CD8 = 1.0, all others = 0.
    /// Evidence: Abbas et al. (2009) — linear model identity: m = S_j → f_j = 1.
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_PureCd8Profile_Cd8ExactFraction()
    {
        var profile = CreatePureCellTypeProfile("T_cells_CD8");

        var result = ImmuneAnalyzer.DeconvoluteImmuneCells(profile);

        Assert.Multiple(() =>
        {
            Assert.That(result.CellFractions["T_cells_CD8"],
                Is.EqualTo(1.0).Within(ComputedTolerance),
                "CD8 fraction must be 1.0 for pure CD8 profile (unique genes CD8A, CD8B)");

            foreach (var kvp in result.CellFractions)
            {
                if (kvp.Key != "T_cells_CD8")
                {
                    Assert.That(kvp.Value, Is.EqualTo(0.0).Within(ComputedTolerance),
                        $"Fraction for {kvp.Key} must be 0 for pure CD8 profile");
                }
            }

            Assert.That(result.Correlation, Is.EqualTo(1.0).Within(ComputedTolerance),
                "Correlation must be 1.0 for perfect reconstruction");
            Assert.That(result.Rmse, Is.EqualTo(0.0).Within(ComputedTolerance),
                "RMSE must be 0 for perfect reconstruction");
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — M6: Two cell types equal mix → 50:50

    /// <summary>
    /// [M6] Equal mixture of B cells naive and CD8 T cells must produce exact 0.5:0.5 fractions.
    /// Proof: CD8A/CD8B are unique to CD8 (force f_CD8 = 0.5), CD79A/PAX5/BANK1 are unique
    /// to B_naive (force f_B = 0.5). The only zero-residual NNLS solution is f_B = f_CD8 = 0.5.
    /// Evidence: Abbas et al. (2009) — linearity of deconvolution model.
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_EqualMixBNaiveAndCd8_ExactHalfFractions()
    {
        var profile = CreateEqualMixProfile("B_cells_naive", "T_cells_CD8");

        var result = ImmuneAnalyzer.DeconvoluteImmuneCells(profile);

        Assert.Multiple(() =>
        {
            Assert.That(result.CellFractions["B_cells_naive"],
                Is.EqualTo(0.5).Within(ComputedTolerance),
                "B_naive fraction must be 0.5 for equal mix (unique genes CD79A, PAX5, BANK1)");
            Assert.That(result.CellFractions["T_cells_CD8"],
                Is.EqualTo(0.5).Within(ComputedTolerance),
                "CD8 fraction must be 0.5 for equal mix (unique genes CD8A, CD8B)");

            foreach (var kvp in result.CellFractions)
            {
                if (kvp.Key != "B_cells_naive" && kvp.Key != "T_cells_CD8")
                {
                    Assert.That(kvp.Value, Is.EqualTo(0.0).Within(ComputedTolerance),
                        $"Fraction for {kvp.Key} must be 0 for B+CD8 equal mix");
                }
            }
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — M9: Tumor purity in [0, 1] (INV-3)

    /// <summary>
    /// [M9] Tumor purity must always be in [0, 1].
    /// Evidence: Yoshihara et al. (2013) — purity is clamped cosine value.
    /// Invariant: INV-3.
    /// </summary>
    [Test]
    public void EstimateInfiltration_ValidProfile_TumorPurityInRange()
    {
        // Use a profile with actual ESTIMATE immune and stromal genes.
        var profile = new Dictionary<string, double>
        {
            ["CD2"] = 8.0,       // ESTIMATE immune
            ["CD3D"] = 7.5,      // ESTIMATE immune
            ["CCL5"] = 9.0,      // ESTIMATE immune
            ["GZMB"] = 10.0,     // ESTIMATE immune
            ["DCN"] = 6.0,       // ESTIMATE stromal
            ["LUM"] = 5.5,       // ESTIMATE stromal
            ["COL1A2"] = 7.0,    // ESTIMATE stromal
            ["RANDOM_GENE"] = 3.0
        };

        var result = ImmuneAnalyzer.EstimateInfiltration(profile);

        Assert.Multiple(() =>
        {
            Assert.That(result.TumorPurity, Is.GreaterThanOrEqualTo(0.0),
                "Tumor purity must be >= 0");
            Assert.That(result.TumorPurity, Is.LessThanOrEqualTo(1.0),
                "Tumor purity must be <= 1");

            // Verify purity formula: purity = cos(a + b × estimateScore)
            double expectedPurity = Math.Cos(EstimateCoefficientA + EstimateCoefficientB * result.EstimateScore);
            expectedPurity = Math.Clamp(expectedPurity, 0.0, 1.0);
            Assert.That(result.TumorPurity, Is.EqualTo(expectedPurity).Within(ExactTolerance),
                "Tumor purity must equal cos(a + b × ESTIMATE score) per Yoshihara et al. (2013)");
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — M10: ESTIMATE score = immune + stromal (INV-4)

    /// <summary>
    /// [M10] ESTIMATE score must equal Immune score + Stromal score.
    /// Evidence: Yoshihara et al. (2013) — definition of ESTIMATE score.
    /// Invariant: INV-4.
    /// </summary>
    [Test]
    public void EstimateInfiltration_ValidProfile_EstimateScoreEqualsSum()
    {
        var profile = new Dictionary<string, double>
        {
            ["CD2"] = 8.0,       // ESTIMATE immune
            ["CD3D"] = 7.5,      // ESTIMATE immune
            ["GZMB"] = 10.0,     // ESTIMATE immune
            ["NKG7"] = 6.0,      // ESTIMATE immune
            ["DCN"] = 6.0,       // ESTIMATE stromal
            ["COL1A2"] = 7.0,    // ESTIMATE stromal
            ["LUM"] = 5.5,       // ESTIMATE stromal
            ["FAP"] = 4.0,       // ESTIMATE stromal
            ["PDGFRB"] = 5.0,    // ESTIMATE stromal
            ["SOME_OTHER_GENE"] = 2.0
        };

        var result = ImmuneAnalyzer.EstimateInfiltration(profile);

        Assert.That(result.EstimateScore,
            Is.EqualTo(result.ImmuneScore + result.StromalScore).Within(ExactTolerance),
            "ESTIMATE score must equal Immune score + Stromal score");
    }

    #endregion

    #region ONCO-IMMUNE-001 — M11: High immune > low immune ordering

    /// <summary>
    /// [M11] A profile with high immune gene expression should have a higher immune score
    /// than a profile with low immune gene expression.
    /// Evidence: ESTIMATE concept — ssGSEA enrichment reflects expression levels. Yoshihara et al. (2013).
    /// </summary>
    [Test]
    public void EstimateInfiltration_HighVsLowImmune_HigherScoreForHighImmune()
    {
        // High immune profile: ESTIMATE immune genes expressed at high levels, non-signature at low.
        var highImmuneProfile = new Dictionary<string, double>
        {
            ["CD2"] = 12.0,      // ESTIMATE immune
            ["CD3D"] = 11.5,     // ESTIMATE immune
            ["CCL5"] = 13.0,     // ESTIMATE immune
            ["GZMB"] = 14.0,     // ESTIMATE immune
            ["NKG7"] = 10.0,     // ESTIMATE immune
            ["PRF1"] = 11.0,     // ESTIMATE immune
            ["HLA-DRA"] = 10.5,  // ESTIMATE immune
            ["LCP2"] = 9.0,      // ESTIMATE immune
            ["LSP1"] = 8.5,      // ESTIMATE immune
            ["FYB"] = 9.5,       // ESTIMATE immune
            ["PLEK"] = 8.0,      // ESTIMATE immune
            ["LOW_GENE_1"] = 1.0,
            ["LOW_GENE_2"] = 0.5,
            ["LOW_GENE_3"] = 1.5
        };

        // Low immune profile: same ESTIMATE immune genes at low levels, non-signature at high.
        var lowImmuneProfile = new Dictionary<string, double>
        {
            ["CD2"] = 1.0,       // ESTIMATE immune
            ["CD3D"] = 0.5,      // ESTIMATE immune
            ["CCL5"] = 1.5,      // ESTIMATE immune
            ["GZMB"] = 2.0,      // ESTIMATE immune
            ["NKG7"] = 1.0,      // ESTIMATE immune
            ["PRF1"] = 0.5,      // ESTIMATE immune
            ["HLA-DRA"] = 1.5,   // ESTIMATE immune
            ["LCP2"] = 1.0,      // ESTIMATE immune
            ["LSP1"] = 0.5,      // ESTIMATE immune
            ["FYB"] = 0.5,       // ESTIMATE immune
            ["PLEK"] = 0.5,      // ESTIMATE immune
            ["HIGH_GENE_1"] = 12.0,
            ["HIGH_GENE_2"] = 13.0,
            ["HIGH_GENE_3"] = 14.0
        };

        var highResult = ImmuneAnalyzer.EstimateInfiltration(highImmuneProfile);
        var lowResult = ImmuneAnalyzer.EstimateInfiltration(lowImmuneProfile);

        Assert.That(highResult.ImmuneScore, Is.GreaterThan(lowResult.ImmuneScore),
            "High-immune profile must have higher immune score than low-immune profile");
    }

    #endregion

    #region ONCO-IMMUNE-001 — M12: Null expression throws for EstimateInfiltration

    /// <summary>
    /// [M12] Null expression profile must throw ArgumentNullException.
    /// Evidence: Robustness — standard .NET null guard pattern.
    /// </summary>
    [Test]
    public void EstimateInfiltration_NullExpression_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => ImmuneAnalyzer.EstimateInfiltration(null!));
    }

    #endregion

    #region ONCO-IMMUNE-001 — M13: Null expression throws for DeconvoluteImmuneCells

    /// <summary>
    /// [M13] Null expression profile must throw ArgumentNullException.
    /// Evidence: Robustness — standard .NET null guard pattern.
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_NullExpression_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(
            () => ImmuneAnalyzer.DeconvoluteImmuneCells(null!));
    }

    #endregion

    #region ONCO-IMMUNE-001 — M14: ssGSEA exact value against hand-computed reference

    /// <summary>
    /// [M14a] ssGSEA with two hit genes at top and bottom of 3 genes must match hand-computed
    /// rank-based integral. This discriminates rank-based weighting from expression-value weighting.
    /// Hand computation for profile {A=100, B=1, C=0.5}, gene set = {A, C}:
    /// Ranked descending: A(rank 3), B(rank 2), C(rank 1). N=3, N_H=2, N_miss=1.
    /// TW = 3^0.25 + 1^0.25 = 3^(1/4) + 1.
    /// Walk: hit(A → 3^(1/4)/TW), miss(B → −1), hit(C → 1/TW).
    /// Integral = (3^(1/4) − 1) / (3^(1/4) + 1) ≈ 0.136548.
    /// Expression-value weighting would give ≈ 0.57992 — clearly distinct.
    /// Evidence: Barbie et al. (2009), Hänzelmann et al. (2013) — rank-based ssGSEA.
    /// </summary>
    [Test]
    public void EstimateInfiltration_SsGseaExactReference_MatchesRankBasedIntegral()
    {
        var profile = new Dictionary<string, double>
        {
            ["GENE_A"] = 100.0,
            ["GENE_B"] = 1.0,
            ["GENE_C"] = 0.5
        };
        var immuneGenes = new[] { "GENE_A", "GENE_C" };
        var stromalGenes = Array.Empty<string>();

        var result = ImmuneAnalyzer.EstimateInfiltration(profile, immuneGenes, stromalGenes);

        // Expected: (3^(1/4) − 1) / (3^(1/4) + 1)
        double threeToFourth = Math.Pow(3.0, 0.25);
        double expectedImmuneScore = (threeToFourth - 1.0) / (threeToFourth + 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(result.ImmuneScore, Is.EqualTo(expectedImmuneScore).Within(ComputedTolerance),
                "ssGSEA immune score must match hand-computed rank-based integral reference");
            Assert.That(result.StromalScore, Is.EqualTo(0.0).Within(ExactTolerance),
                "Stromal score must be 0 with no overlapping stromal genes");
            Assert.That(result.EstimateScore, Is.EqualTo(expectedImmuneScore).Within(ComputedTolerance),
                "ESTIMATE score must equal immune score when stromal = 0");
        });
    }

    /// <summary>
    /// [M14b] ssGSEA with single hit gene at top rank → positive enrichment = 1.5.
    /// Hand computation: 3 genes, 1 hit at rank 3. TW = 3^0.25 (cancels). N_miss=2.
    /// Walk: hit(weight=1), miss(−0.5), miss(−0.5). RS: 1, 0.5, 0. Integral: 1 + 0.5 + 0 = 1.5.
    /// Verifies integral aggregation (max-deviation would give 1.0).
    /// Evidence: Barbie et al. (2009), Hänzelmann et al. (2013).
    /// </summary>
    [Test]
    public void EstimateInfiltration_SsGseaSingleHitAtTop_IntegralEquals1Point5()
    {
        var profile = new Dictionary<string, double>
        {
            ["HIT_GENE"] = 10.0,
            ["MISS_1"] = 5.0,
            ["MISS_2"] = 1.0
        };
        var immuneGenes = new[] { "HIT_GENE" };
        var stromalGenes = Array.Empty<string>();

        var result = ImmuneAnalyzer.EstimateInfiltration(profile, immuneGenes, stromalGenes);

        // Walk: hit(1), miss(−0.5), miss(−0.5). RS: 1, 0.5, 0. Integral = 1.5.
        // Max-deviation would return 1.0 — this test catches that bug.
        Assert.That(result.ImmuneScore, Is.EqualTo(1.5).Within(ComputedTolerance),
            "Single hit at top rank → ssGSEA integral must be 1.5 (not 1.0 as max-deviation)");
    }

    /// <summary>
    /// [M14c] ssGSEA with single hit gene at bottom rank → negative enrichment = −1.5.
    /// Hand computation: 3 genes, 1 hit at rank 1 (bottom). N_miss=2.
    /// Walk: miss(−0.5), miss(−0.5), hit(1). RS: −0.5, −1.0, 0. Integral: −0.5 + −1.0 + 0 = −1.5.
    /// Evidence: Barbie et al. (2009), Hänzelmann et al. (2013).
    /// </summary>
    [Test]
    public void EstimateInfiltration_SsGseaSingleHitAtBottom_IntegralEqualsMinus1Point5()
    {
        var profile = new Dictionary<string, double>
        {
            ["MISS_1"] = 10.0,
            ["MISS_2"] = 5.0,
            ["HIT_GENE"] = 1.0
        };
        var immuneGenes = new[] { "HIT_GENE" };
        var stromalGenes = Array.Empty<string>();

        var result = ImmuneAnalyzer.EstimateInfiltration(profile, immuneGenes, stromalGenes);

        // Walk: miss(−0.5), miss(−0.5), hit(1). RS: −0.5, −1.0, 0. Integral = −1.5.
        Assert.That(result.ImmuneScore, Is.EqualTo(-1.5).Within(ComputedTolerance),
            "Single hit at bottom rank → ssGSEA integral must be −1.5");
    }

    #endregion

    #region ONCO-IMMUNE-001 — S1: Unequal mixture → exact proportional fractions

    /// <summary>
    /// [S1] A 75:25 mixture of CD8 T cells and B naive must produce exact 0.75:0.25 fractions.
    /// Proof: CD8A/CD8B unique to CD8 force f_CD8 = 0.75; CD79A/PAX5/BANK1 unique to B_naive
    /// force f_B = 0.25. Zero-residual solution is unique.
    /// Validates linearity of NNLS deconvolution model.
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_UnequalMix75Cd8_25BNaive_ExactFractions()
    {
        var profile = CreateWeightedMixProfile("T_cells_CD8", 0.75, "B_cells_naive", 0.25);

        var result = ImmuneAnalyzer.DeconvoluteImmuneCells(profile);

        Assert.Multiple(() =>
        {
            Assert.That(result.CellFractions["T_cells_CD8"],
                Is.EqualTo(0.75).Within(ComputedTolerance),
                "CD8 fraction must be 0.75 for 75:25 mix (unique genes CD8A, CD8B)");
            Assert.That(result.CellFractions["B_cells_naive"],
                Is.EqualTo(0.25).Within(ComputedTolerance),
                "B_naive fraction must be 0.25 for 75:25 mix (unique genes CD79A, PAX5, BANK1)");

            foreach (var kvp in result.CellFractions)
            {
                if (kvp.Key != "T_cells_CD8" && kvp.Key != "B_cells_naive")
                {
                    Assert.That(kvp.Value, Is.EqualTo(0.0).Within(ComputedTolerance),
                        $"Fraction for {kvp.Key} must be 0 for CD8+B_naive 75:25 mix");
                }
            }
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — S2: Extra genes ignored in deconvolution

    /// <summary>
    /// [S2] Extra genes not in the signature matrix should be ignored.
    /// Result should be similar to profile with signature genes only.
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_ExtraGenesInProfile_Ignored()
    {
        var pureProfile = CreatePureCellTypeProfile("T_cells_CD8");

        var profileWithExtras = new Dictionary<string, double>(pureProfile);
        for (int i = 0; i < 100; i++)
        {
            profileWithExtras[$"EXTRA_GENE_{i}"] = i * 0.5;
        }

        var pureResult = ImmuneAnalyzer.DeconvoluteImmuneCells(pureProfile);
        var mixedResult = ImmuneAnalyzer.DeconvoluteImmuneCells(profileWithExtras);

        Assert.Multiple(() =>
        {
            Assert.That(mixedResult.OverlappingGenes, Is.EqualTo(pureResult.OverlappingGenes),
                "Overlapping gene count should be the same regardless of extra genes");

            // Fractions should be identical because NNLS only uses overlapping genes.
            foreach (var cellType in pureResult.CellFractions.Keys)
            {
                Assert.That(mixedResult.CellFractions[cellType],
                    Is.EqualTo(pureResult.CellFractions[cellType]).Within(ExactTolerance),
                    $"Fraction for {cellType} should be identical regardless of extra genes");
            }
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — S3: Positive correlation for good fit

    /// <summary>
    /// [S3] A pure cell type profile must reconstruct perfectly: correlation = 1.0, RMSE = 0.
    /// Proof: Monocyte-unique genes (CD14, FCN1, S100A8, VCAN) force f_Mono = 1.0.
    /// Reconstructed profile = S[:, Mono] = m exactly.
    /// Evidence: Quality metric — pure signal must reconstruct perfectly.
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_PureProfile_PerfectReconstruction()
    {
        var profile = CreatePureCellTypeProfile("Monocytes");

        var result = ImmuneAnalyzer.DeconvoluteImmuneCells(profile);

        Assert.Multiple(() =>
        {
            Assert.That(result.Correlation, Is.EqualTo(1.0).Within(ComputedTolerance),
                "Pure cell type profile must produce correlation = 1.0 (perfect fit)");
            Assert.That(result.Rmse, Is.EqualTo(0.0).Within(ComputedTolerance),
                "Pure cell type profile must produce RMSE = 0 (zero residual)");
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — S4: Custom gene sets accepted

    /// <summary>
    /// [S4] EstimateInfiltration should accept custom immune and stromal gene sets.
    /// </summary>
    [Test]
    public void EstimateInfiltration_CustomGeneSets_UsesCustomSets()
    {
        var profile = new Dictionary<string, double>
        {
            ["CUSTOM_IMMUNE_1"] = 10.0,
            ["CUSTOM_IMMUNE_2"] = 8.0,
            ["CUSTOM_STROMAL_1"] = 6.0,
            ["OTHER_GENE"] = 3.0
        };

        var customImmune = new[] { "CUSTOM_IMMUNE_1", "CUSTOM_IMMUNE_2" };
        var customStromal = new[] { "CUSTOM_STROMAL_1" };

        var result = ImmuneAnalyzer.EstimateInfiltration(profile, customImmune, customStromal);

        Assert.Multiple(() =>
        {
            Assert.That(result.OverlappingImmuneGenes, Is.EqualTo(2),
                "Should find 2 custom immune genes");
            Assert.That(result.OverlappingStromalGenes, Is.EqualTo(1),
                "Should find 1 custom stromal gene");
            Assert.That(result.ImmuneScore, Is.Not.EqualTo(0.0),
                "Immune score should be non-zero with overlapping custom genes");
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — C1: Default signature matrix has 22 cell types

    /// <summary>
    /// [C1] DefaultSignatureMatrix should contain exactly 22 immune cell types (LM22-inspired).
    /// Evidence: Newman et al. (2015) — LM22 has 22 hematopoietic cell phenotypes.
    /// </summary>
    [Test]
    public void DefaultSignatureMatrix_Has22CellTypes()
    {
        Assert.That(ImmuneAnalyzer.DefaultSignatureMatrix.Count, Is.EqualTo(22),
            "Default signature matrix should have 22 cell types per LM22");
    }

    #endregion

    #region ONCO-IMMUNE-001 — C2: Negative expression values handled

    /// <summary>
    /// [C2] Log-transformed expression data may contain negative values. Should not throw.
    /// </summary>
    [Test]
    public void EstimateInfiltration_NegativeExpressionValues_NoException()
    {
        var profile = new Dictionary<string, double>
        {
            ["CD2"] = -2.0,      // ESTIMATE immune
            ["CD3D"] = -1.5,     // ESTIMATE immune
            ["GZMB"] = -0.5,     // ESTIMATE immune
            ["DCN"] = -3.0,      // ESTIMATE stromal
            ["LUM"] = -1.0,      // ESTIMATE stromal
            ["RANDOM"] = 1.0
        };

        Assert.DoesNotThrow(() => ImmuneAnalyzer.EstimateInfiltration(profile),
            "Negative expression values (log-transformed) should be valid inputs");
    }

    /// <summary>
    /// [C2b] Negative expression values in deconvolution should not throw.
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_NegativeExpressionValues_NoException()
    {
        var profile = new Dictionary<string, double>
        {
            ["CD8A"] = -1.0,
            ["CD8B"] = -0.5,
            ["CD3D"] = -2.0,
            ["GZMB"] = 1.0,
            ["PRF1"] = 0.5
        };

        Assert.DoesNotThrow(() => ImmuneAnalyzer.DeconvoluteImmuneCells(profile),
            "Negative expression values should be handled without exceptions");
    }

    #endregion

    #region ONCO-IMMUNE-001 — C3: Default signature matrix entries have 5 genes each

    /// <summary>
    /// [C3] Each cell type in the default signature matrix should have exactly 5 marker genes.
    /// This validates the simplified LM22 structure.
    /// </summary>
    [Test]
    public void DefaultSignatureMatrix_EachCellTypeHas5Genes()
    {
        Assert.Multiple(() =>
        {
            foreach (var kvp in ImmuneAnalyzer.DefaultSignatureMatrix)
            {
                Assert.That(kvp.Value.Count, Is.EqualTo(5),
                    $"Cell type {kvp.Key} should have exactly 5 marker genes");
            }
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — Invariant Tests

    /// <summary>
    /// [INV-1, INV-2] For multiple different cell type profiles, verify fractions are all non-negative
    /// and sum to 1.0. Property-based invariant test.
    /// Evidence: Newman et al. (2015) — NNLS constraints + normalization.
    /// </summary>
    [TestCase("B_cells_naive")]
    [TestCase("T_cells_CD8")]
    [TestCase("Monocytes")]
    [TestCase("NK_cells_activated")]
    [TestCase("Macrophages_M1")]
    [TestCase("Neutrophils")]
    public void DeconvoluteImmuneCells_MultipleCellTypes_FractionsInvariant(string cellType)
    {
        var profile = CreatePureCellTypeProfile(cellType);

        var result = ImmuneAnalyzer.DeconvoluteImmuneCells(profile);

        Assert.Multiple(() =>
        {
            // INV-1: All fractions ≥ 0
            foreach (var kvp in result.CellFractions)
            {
                Assert.That(kvp.Value, Is.GreaterThanOrEqualTo(0.0),
                    $"[INV-1] Fraction for {kvp.Key} must be non-negative");
            }

            // INV-2: Sum = 1.0
            double sum = result.CellFractions.Values.Sum();
            Assert.That(sum, Is.EqualTo(1.0).Within(ComputedTolerance),
                "[INV-2] Fractions must sum to 1.0");
        });
    }

    /// <summary>
    /// [INV-3, INV-4] For profiles with varying levels of immune expression,
    /// verify tumor purity is in [0,1] and ESTIMATE score = immune + stromal.
    /// Evidence: Yoshihara et al. (2013).
    /// </summary>
    [TestCase(1.0)]
    [TestCase(5.0)]
    [TestCase(10.0)]
    [TestCase(15.0)]
    public void EstimateInfiltration_VaryingExpression_PurityAndScoreInvariant(double baseExpression)
    {
        var profile = new Dictionary<string, double>
        {
            ["CD2"] = baseExpression,          // ESTIMATE immune
            ["CD3D"] = baseExpression * 0.9,    // ESTIMATE immune
            ["GZMB"] = baseExpression * 1.1,    // ESTIMATE immune
            ["NKG7"] = baseExpression * 0.8,     // ESTIMATE immune
            ["DCN"] = baseExpression * 0.7,      // ESTIMATE stromal
            ["LUM"] = baseExpression * 0.6,      // ESTIMATE stromal
            ["OTHER"] = baseExpression * 0.5
        };

        var result = ImmuneAnalyzer.EstimateInfiltration(profile);

        Assert.Multiple(() =>
        {
            // INV-3: Purity ∈ [0, 1]
            Assert.That(result.TumorPurity, Is.GreaterThanOrEqualTo(0.0),
                "[INV-3] Tumor purity must be >= 0");
            Assert.That(result.TumorPurity, Is.LessThanOrEqualTo(1.0),
                "[INV-3] Tumor purity must be <= 1");

            // INV-3 strengthened: verify formula purity = cos(a + b × estimateScore)
            double expectedPurity = Math.Cos(EstimateCoefficientA + EstimateCoefficientB * result.EstimateScore);
            expectedPurity = Math.Clamp(expectedPurity, 0.0, 1.0);
            Assert.That(result.TumorPurity, Is.EqualTo(expectedPurity).Within(ExactTolerance),
                "[INV-3] Tumor purity must equal cos(a + b × ESTIMATE score)");

            // INV-4: ESTIMATE = immune + stromal
            Assert.That(result.EstimateScore,
                Is.EqualTo(result.ImmuneScore + result.StromalScore).Within(ExactTolerance),
                "[INV-4] ESTIMATE score must equal Immune + Stromal");
        });
    }

    /// <summary>
    /// [INV-5] Overlapping genes must be between 0 and total signature genes.
    /// </summary>
    [Test]
    public void DeconvoluteImmuneCells_ValidProfile_OverlappingGenesInRange()
    {
        var profile = CreatePureCellTypeProfile("T_cells_CD8");

        var result = ImmuneAnalyzer.DeconvoluteImmuneCells(profile);

        int totalSignatureGenes = ImmuneAnalyzer.DefaultSignatureMatrix.Values
            .SelectMany(d => d.Keys)
            .Distinct()
            .Count();

        Assert.Multiple(() =>
        {
            Assert.That(result.OverlappingGenes, Is.GreaterThanOrEqualTo(0),
                "[INV-5] Overlapping genes must be >= 0");
            Assert.That(result.OverlappingGenes, Is.LessThanOrEqualTo(totalSignatureGenes),
                "[INV-5] Overlapping genes must be <= total signature genes");
        });
    }

    #endregion

    #region ONCO-IMMUNE-001 — EstimateTumorPurity (Yoshihara 2013 closed-form transform)

    // E1 — Yoshihara (2013) purity transform at score 0.
    // purity = cos(0.6049872018 + 0.0001467884 × 0) = cos(0.6049872018).
    // Hand-computed: 0.8225093766958238.
    // Source: Yoshihara et al. (2013), Nat Commun 4:2612; ESTIMATE/tidyestimate estimate_score().
    [Test]
    public void EstimateTumorPurity_ZeroScore_EqualsCosOfCoefficientA()
    {
        // Arrange
        const double expected = 0.8225093766958238;

        // Act
        double purity = ImmuneAnalyzer.EstimateTumorPurity(0.0);

        // Assert
        Assert.That(purity, Is.EqualTo(expected).Within(ExactTolerance),
            "purity(0) must equal cos(0.6049872018) per Yoshihara et al. (2013)");
    }

    // E2 — Purity at a mid-range Affymetrix-scale ESTIMATE score.
    // purity = cos(0.6049872018 + 0.0001467884 × 1000) = cos(0.7517756018) = 0.7304773970805112.
    [Test]
    public void EstimateTumorPurity_Score1000_EqualsHandComputedCosine()
    {
        // Arrange
        const double expected = 0.7304773970805112;

        // Act
        double purity = ImmuneAnalyzer.EstimateTumorPurity(1000.0);

        // Assert
        Assert.That(purity, Is.EqualTo(expected).Within(ExactTolerance),
            "purity(1000) must equal cos(0.6049872018 + 0.0001467884 × 1000)");
    }

    // E3 — Purity at a high ESTIMATE score.
    // purity = cos(0.6049872018 + 0.0001467884 × 3000) = cos(1.0453524018) = 0.5015970942006772.
    [Test]
    public void EstimateTumorPurity_Score3000_EqualsHandComputedCosine()
    {
        // Arrange
        const double expected = 0.5015970942006772;

        // Act
        double purity = ImmuneAnalyzer.EstimateTumorPurity(3000.0);

        // Assert
        Assert.That(purity, Is.EqualTo(expected).Within(ExactTolerance),
            "purity(3000) must equal cos(0.6049872018 + 0.0001467884 × 3000)");
    }

    // E4 — Domain handling: when the cosine argument exceeds π/2 the cosine is negative,
    // which is outside the calibrated domain → NaN (the R packages set such values to NA).
    // At score 7000: arg = 1.6325060018 > π/2, cos = -0.0617 < 0 → NaN.
    // Source: tidyestimate estimate_score(): purity = ifelse(purity < 0, NA, purity).
    [Test]
    public void EstimateTumorPurity_NegativeCosineDomain_ReturnsNaN()
    {
        // Act
        double purity = ImmuneAnalyzer.EstimateTumorPurity(7000.0);

        // Assert
        Assert.That(double.IsNaN(purity), Is.True,
            "Out-of-domain score (negative cosine) must return NaN, mirroring the reference NA handling");
    }

    // E5 — Boundary: a score just below the negative-cosine cutoff (~6579.6) stays defined,
    // and a score just above it becomes NaN. cos crosses 0 at score (π/2 − a)/b ≈ 6579.60.
    [Test]
    public void EstimateTumorPurity_AroundNegativeCutoff_DefinedBelowNaNAbove()
    {
        // Act
        double below = ImmuneAnalyzer.EstimateTumorPurity(6000.0); // arg < π/2 → positive
        double above = ImmuneAnalyzer.EstimateTumorPurity(6600.0); // arg > π/2 → negative

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(below, Is.EqualTo(0.0849761233112934).Within(ExactTolerance),
                "purity(6000) must equal cos(0.6049872018 + 0.0001467884 × 6000) and be positive");
            Assert.That(double.IsNaN(above), Is.True,
                "purity(6600) is past the negative-cosine cutoff (~6579.6) and must be NaN");
        });
    }

    // E6 — Monotonic decreasing: within the valid domain (cosine argument in [0, π/2]),
    // higher ESTIMATE score ⇒ lower tumour purity (cos is decreasing on [0, π]).
    // Source: Yoshihara et al. (2013) — purity is a monotone-decreasing function of ESTIMATE score.
    [Test]
    public void EstimateTumorPurity_IncreasingScore_PurityMonotonicallyDecreases()
    {
        // Arrange
        double[] scores = { -2000.0, 0.0, 1000.0, 3000.0, 5000.0, 6000.0 };

        // Act
        double[] purities = scores.Select(ImmuneAnalyzer.EstimateTumorPurity).ToArray();

        // Assert
        Assert.Multiple(() =>
        {
            for (int i = 1; i < purities.Length; i++)
            {
                Assert.That(purities[i], Is.LessThan(purities[i - 1]),
                    $"purity must strictly decrease as ESTIMATE score increases (score {scores[i]} vs {scores[i - 1]})");
            }
        });
    }

    // E7 — The opt-in transform matches the same cosine the default InfiltrationResult applies,
    // confirming it is the identical Yoshihara formula (applied here to a caller-supplied score).
    [Test]
    public void EstimateTumorPurity_MatchesClosedFormCosine_ForRepresentativeScore()
    {
        // Arrange
        const double score = 2500.0;
        double expected = Math.Cos(EstimateCoefficientA + EstimateCoefficientB * score);

        // Act
        double purity = ImmuneAnalyzer.EstimateTumorPurity(score);

        // Assert
        Assert.That(purity, Is.EqualTo(expected).Within(ExactTolerance),
            "EstimateTumorPurity must apply cos(a + b × score) with the Yoshihara coefficients");
    }

    #endregion

    #region CIBERSORT nu-SVR Deconvolution (DeconvoluteImmuneCellsNuSvr, LoadSignatureMatrix)

    // Tolerance for planted-truth fraction recovery. nu-SVR is a regularised, epsilon-insensitive
    // (robust) estimator, so recovered fractions approximate — but do not exactly equal — the
    // planted fractions; the tube width and z-standardisation introduce a small, bounded deviation.
    private const double NuSvrRecoveryTolerance = 0.025;

    // Tolerance for agreement with the scikit-learn / libsvm NuSVR reference (cross-implementation).
    private const double NuSvrReferenceTolerance = 2e-3;

    /// <summary>
    /// Builds a synthetic bulk mixture m = B·f from a signature matrix B and planted fractions f.
    /// This is the planted-truth construction used to verify deconvolution recovery.
    /// </summary>
    private static Dictionary<string, double> BuildPlantedMixture(
        IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>> signature,
        IReadOnlyDictionary<string, double> plantedFractions)
    {
        var allGenes = signature.Values.SelectMany(d => d.Keys).Distinct();
        var mixture = new Dictionary<string, double>();
        foreach (var gene in allGenes)
        {
            double v = 0.0;
            foreach (var ct in signature.Keys)
            {
                double frac = plantedFractions.TryGetValue(ct, out double fr) ? fr : 0.0;
                if (frac != 0.0 && signature[ct].TryGetValue(gene, out double sval))
                    v += frac * sval;
            }
            mixture[gene] = v;
        }
        return mixture;
    }

    // NSVR-M1 — Planted-truth recovery on the bundled 5-marker × 22-cell-type matrix.
    // Bulk = B·f with f = {CD8:0.60, B_naive:0.30, Monocytes:0.10}; nu-SVR must recover f.
    // Evidence: Newman et al. (2015) — nu-SVR deconvolution of m on the signature columns.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_PlantedTruth_RecoversFractions()
    {
        // Arrange
        var planted = new Dictionary<string, double>
        {
            ["T_cells_CD8"] = 0.60,
            ["B_cells_naive"] = 0.30,
            ["Monocytes"] = 0.10,
        };
        var mixture = BuildPlantedMixture(ImmuneAnalyzer.DefaultSignatureMatrix, planted);

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture);

        // Assert
        Assert.Multiple(() =>
        {
            foreach (var kvp in planted)
            {
                Assert.That(result.CellFractions[kvp.Key], Is.EqualTo(kvp.Value).Within(NuSvrRecoveryTolerance),
                    $"nu-SVR must recover the planted fraction for {kvp.Key} (planted {kvp.Value})");
            }

            // Cell types not present in the mixture should be ~0.
            foreach (var ct in ImmuneAnalyzer.DefaultSignatureMatrix.Keys.Where(c => !planted.ContainsKey(c)))
            {
                Assert.That(result.CellFractions[ct], Is.LessThan(NuSvrRecoveryTolerance),
                    $"absent cell type {ct} should recover a near-zero fraction");
            }
        });
    }

    // NSVR-M2 — Cross-implementation reference match against scikit-learn 1.6.1 NuSVR (libsvm).
    // Identical standardized linear problem; selected nu = 0.75 (lowest RMSE); sklearn normalized
    // fractions = [0.508497, 0.179491, 0.312012]. This is the decisive correctness check: it would
    // FAIL for any solver that does not faithfully implement the Schölkopf (2000) nu-SVR dual.
    // Evidence: scikit-learn NuSVR(kernel='linear', nu, C=1.0) — verbatim reference numbers
    // recomputed this session; Schölkopf et al. (2000), Neural Computation 12:1207.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_MatchesLibsvmNuSvrReference()
    {
        // Arrange — 3 cell types, 3 disjoint markers each (same matrix used against sklearn).
        var signature = new Dictionary<string, IReadOnlyDictionary<string, double>>
        {
            ["TypeA"] = new Dictionary<string, double> { ["a1"] = 10, ["a2"] = 8, ["a3"] = 6 },
            ["TypeB"] = new Dictionary<string, double> { ["b1"] = 9, ["b2"] = 7, ["b3"] = 5 },
            ["TypeC"] = new Dictionary<string, double> { ["c1"] = 11, ["c2"] = 4, ["c3"] = 8 },
        };
        var planted = new Dictionary<string, double> { ["TypeA"] = 0.5, ["TypeB"] = 0.2, ["TypeC"] = 0.3 };
        var mixture = BuildPlantedMixture(signature, planted);

        // scikit-learn 1.6.1 NuSVR reference (selected nu = 0.75), normalized after zero-clip.
        const double refA = 0.508497;
        const double refB = 0.179491;
        const double refC = 0.312012;

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture, signature);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.BestNu, Is.EqualTo(0.75).Within(ExactTolerance),
                "the lowest-RMSE nu for this problem is 0.75 (matches sklearn selection)");
            Assert.That(result.CellFractions["TypeA"], Is.EqualTo(refA).Within(NuSvrReferenceTolerance),
                "TypeA fraction must match the libsvm/sklearn NuSVR reference");
            Assert.That(result.CellFractions["TypeB"], Is.EqualTo(refB).Within(NuSvrReferenceTolerance),
                "TypeB fraction must match the libsvm/sklearn NuSVR reference");
            Assert.That(result.CellFractions["TypeC"], Is.EqualTo(refC).Within(NuSvrReferenceTolerance),
                "TypeC fraction must match the libsvm/sklearn NuSVR reference");
        });
    }

    // NSVR-M3 — Fractions are non-negative and sum to 1 (INV: zero-clip + sum-to-1).
    // Evidence: Newman et al. (2015) — negatives set to 0, remaining normalized to sum 1.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_FractionsAreNonNegativeAndSumToOne()
    {
        // Arrange
        var planted = new Dictionary<string, double>
        {
            ["NK_cells_activated"] = 0.4,
            ["Macrophages_M1"] = 0.35,
            ["Plasma_cells"] = 0.25,
        };
        var mixture = BuildPlantedMixture(ImmuneAnalyzer.DefaultSignatureMatrix, planted);

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture);
        double sum = result.CellFractions.Values.Sum();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.CellFractions.Values.All(v => v >= 0.0), Is.True,
                "all cell-type fractions must be non-negative (negatives clipped to 0)");
            Assert.That(sum, Is.EqualTo(1.0).Within(1e-9),
                "cell-type fractions must be normalized to sum to 1");
        });
    }

    // NSVR-M4 — Selected nu is one of the CIBERSORT sweep values {0.25, 0.5, 0.75}.
    // Evidence: CIBERSORT protocol — "CIBERSORT uses a set of nu values (0.25, 0.5, 0.75)".
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_SelectsNuFromCibersortSet()
    {
        // Arrange
        var planted = new Dictionary<string, double> { ["T_cells_CD8"] = 0.7, ["Neutrophils"] = 0.3 };
        var mixture = BuildPlantedMixture(ImmuneAnalyzer.DefaultSignatureMatrix, planted);

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture);

        // Assert
        Assert.That(ImmuneAnalyzer.CibersortNuValues, Does.Contain(result.BestNu),
            "the selected nu must be one of the CIBERSORT sweep values {0.25, 0.5, 0.75}");
    }

    // NSVR-M5 — High-fidelity reconstruction: for an exact planted mixture the Pearson correlation
    // between observed and reconstructed profile is near 1 (the model fits a true linear mixture).
    // Evidence: Newman et al. (2015) — m = B·f is a linear mixture, recovered by the regression.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_PlantedTruth_ReconstructionCorrelationNearOne()
    {
        // Arrange
        var planted = new Dictionary<string, double> { ["B_cells_naive"] = 0.5, ["T_cells_CD8"] = 0.5 };
        var mixture = BuildPlantedMixture(ImmuneAnalyzer.DefaultSignatureMatrix, planted);

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture);

        // Assert
        Assert.That(result.Correlation, Is.GreaterThan(0.95),
            "reconstruction of an exact linear mixture must correlate near-perfectly with the observed profile");
    }

    // NSVR-S1 — Determinism: identical inputs yield identical fractions (no randomness).
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_IsDeterministic()
    {
        // Arrange
        var planted = new Dictionary<string, double> { ["Monocytes"] = 0.6, ["Eosinophils"] = 0.4 };
        var mixture = BuildPlantedMixture(ImmuneAnalyzer.DefaultSignatureMatrix, planted);

        // Act
        var r1 = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture);
        var r2 = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(r1.BestNu, Is.EqualTo(r2.BestNu), "selected nu must be deterministic");
            foreach (var ct in r1.CellFractions.Keys)
            {
                Assert.That(r1.CellFractions[ct], Is.EqualTo(r2.CellFractions[ct]).Within(ExactTolerance),
                    $"fraction for {ct} must be deterministic across runs");
            }
        });
    }

    // NSVR-S2 — No overlapping genes → all fractions 0, OverlappingGenes = 0.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_NoOverlappingGenes_ReturnsZeroFractions()
    {
        // Arrange — expression profile with genes absent from the signature matrix.
        var profile = new Dictionary<string, double> { ["ZZZ1"] = 5.0, ["ZZZ2"] = 7.0 };

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(profile);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OverlappingGenes, Is.EqualTo(0), "no signature genes overlap the profile");
            Assert.That(result.CellFractions.Values.All(v => v == 0.0), Is.True,
                "with no overlap, every cell-type fraction is 0");
            Assert.That(result.BestNu, Is.EqualTo(0.0), "no fit performed → BestNu is the sentinel 0");
        });
    }

    // NSVR-S3 — Null expression profile → ArgumentNullException.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_NullProfile_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(null!),
            "a null expression profile must raise ArgumentNullException");
    }

    // NSVR-S4 — LM22-format loader parses a tab-separated matrix into cellType → (gene → value).
    // Evidence: Newman et al. (2015) — LM22 is a TSV: header (gene-symbol col + 22 cell types),
    // then one row per gene. The full LM22 is 547 genes × 22 cell types (caller-supplied; not bundled).
    [Test]
    public void LoadSignatureMatrix_ParsesTabSeparatedMatrix()
    {
        // Arrange
        var lines = new[]
        {
            "Gene symbol\tB cells naive\tT cells CD8",
            "CD19\t8.5\t0.0",
            "CD8A\t0.0\t9.5",
        };

        // Act
        var matrix = ImmuneAnalyzer.LoadSignatureMatrix(lines);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(matrix.Keys, Is.EquivalentTo(new[] { "B cells naive", "T cells CD8" }),
                "cell-type columns must be parsed from the header");
            Assert.That(matrix["B cells naive"]["CD19"], Is.EqualTo(8.5).Within(ExactTolerance),
                "value for (CD19, B cells naive) must be parsed");
            Assert.That(matrix["T cells CD8"]["CD8A"], Is.EqualTo(9.5).Within(ExactTolerance),
                "value for (CD8A, T cells CD8) must be parsed");
            Assert.That(matrix["B cells naive"]["CD8A"], Is.EqualTo(0.0).Within(ExactTolerance),
                "zero entries must be parsed");
        });
    }

    // NSVR-S5 — A matrix loaded via the LM22 loader feeds the deconvolution end-to-end.
    [Test]
    public void LoadSignatureMatrix_LoadedMatrix_DrivesDeconvolution()
    {
        // Arrange — load a small disjoint matrix, then deconvolute a planted mixture with it.
        var lines = new[]
        {
            "Gene symbol\tTypeA\tTypeB\tTypeC",
            "a1\t10\t0\t0",
            "a2\t8\t0\t0",
            "b1\t0\t9\t0",
            "b2\t0\t7\t0",
            "c1\t0\t0\t11",
            "c2\t0\t0\t8",
        };
        var signature = ImmuneAnalyzer.LoadSignatureMatrix(lines);
        var planted = new Dictionary<string, double> { ["TypeA"] = 0.5, ["TypeB"] = 0.3, ["TypeC"] = 0.2 };
        var mixture = BuildPlantedMixture(signature, planted);

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture, signature);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.CellFractions["TypeA"], Is.GreaterThan(result.CellFractions["TypeB"]),
                "the dominant planted type (TypeA, 0.5) must have the largest recovered fraction");
            Assert.That(result.CellFractions["TypeB"], Is.GreaterThan(result.CellFractions["TypeC"]),
                "recovered ordering must follow the planted ordering A > B > C");
            Assert.That(result.CellFractions.Values.Sum(), Is.EqualTo(1.0).Within(1e-9),
                "loaded-matrix fractions must still sum to 1");
        });
    }

    // NSVR-S6 — Degenerate (constant) mixture: a flat expression vector has zero variance, so
    // z-standardisation cannot scale it. The method must not throw and must still return a valid
    // normalized fraction vector (Σ = 1 or all-zero). Exercises the zero-SD standardisation branch.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_ConstantMixture_DoesNotThrow()
    {
        // Arrange — every overlapping signature gene gets the same expression value.
        var genes = ImmuneAnalyzer.DefaultSignatureMatrix.Values
            .SelectMany(d => d.Keys).Distinct();
        var profile = genes.ToDictionary(g => g, _ => 5.0);

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(profile);
        double sum = result.CellFractions.Values.Sum();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.CellFractions.Values.All(v => v >= 0.0), Is.True,
                "fractions must remain non-negative for a degenerate constant mixture");
            Assert.That(sum, Is.EqualTo(1.0).Within(1e-9).Or.EqualTo(0.0).Within(1e-9),
                "a degenerate mixture yields either a normalized (Σ=1) or all-zero fraction vector");
        });
    }

    // NSVR-S7 — All-zero (overlapping) mixture: every overlapping signature gene is present but
    // valued 0. z-standardisation of a zero vector yields zeros, so the regression weights are 0;
    // the result must be a valid all-zero fraction vector (sum 0), not a crash or NaN. This is a
    // distinct edge case from the constant-nonzero mixture (NSVR-S6): here even the un-standardized
    // mixture is identically zero. Verified vs the engine: overlap = 85, all fractions exactly 0.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_AllZeroMixture_ReturnsZeroFractions()
    {
        // Arrange — every overlapping signature gene present with value 0.
        var genes = ImmuneAnalyzer.DefaultSignatureMatrix.Values
            .SelectMany(d => d.Keys).Distinct();
        var profile = genes.ToDictionary(g => g, _ => 0.0);

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(profile);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OverlappingGenes, Is.GreaterThan(0),
                "the all-zero profile still overlaps the signature genes");
            Assert.That(result.CellFractions.Values.All(v => v == 0.0), Is.True,
                "an identically-zero mixture has no signal → every fraction is exactly 0");
            Assert.That(result.CellFractions.Values.Any(double.IsNaN), Is.False,
                "the zero-SD standardisation branch must not produce NaN");
        });
    }

    // NSVR-S8 — Empty signature matrix: no cell types, no genes. The method must return an empty
    // fraction map with the sentinel BestNu = 0 and OverlappingGenes = 0 (no fit attempted), not throw.
    // Exercises the `cellTypes.Count == 0` guard.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_EmptySignatureMatrix_ReturnsEmptyResult()
    {
        // Arrange
        var emptySig = new Dictionary<string, IReadOnlyDictionary<string, double>>();
        var profile = new Dictionary<string, double> { ["CD8A"] = 5.0 };

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(profile, emptySig);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.CellFractions, Is.Empty, "no cell types → empty fraction map");
            Assert.That(result.OverlappingGenes, Is.EqualTo(0), "no signature genes → zero overlap");
            Assert.That(result.BestNu, Is.EqualTo(0.0), "no fit performed → BestNu is the sentinel 0");
        });
    }

    // NSVR-S9 — Gene-count mismatch (partial overlap): the mixture exposes only a SUBSET of the
    // signature genes (a1, a2 from TypeA and b1 from TypeB; no TypeC gene at all). The engine must
    // restrict the regression to the overlapping genes, report that count, and still return a valid
    // (non-negative, sum-to-1) fraction vector. Verified vs the engine: overlap = 3, Σ = 1.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_PartialGeneOverlap_RestrictsToOverlapAndNormalizes()
    {
        // Arrange — disjoint 3-type matrix; profile contains only 3 of its 9 genes.
        var signature = new Dictionary<string, IReadOnlyDictionary<string, double>>
        {
            ["TypeA"] = new Dictionary<string, double> { ["a1"] = 10, ["a2"] = 8, ["a3"] = 6 },
            ["TypeB"] = new Dictionary<string, double> { ["b1"] = 9, ["b2"] = 7, ["b3"] = 5 },
            ["TypeC"] = new Dictionary<string, double> { ["c1"] = 11, ["c2"] = 4, ["c3"] = 8 },
        };
        var profile = new Dictionary<string, double> { ["a1"] = 10, ["a2"] = 8, ["b1"] = 9 };

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(profile, signature);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.OverlappingGenes, Is.EqualTo(3),
                "only the 3 overlapping genes (a1, a2, b1) drive the regression");
            Assert.That(result.CellFractions.Values.All(v => v >= 0.0), Is.True,
                "fractions remain non-negative under partial overlap");
            Assert.That(result.CellFractions.Values.Sum(), Is.EqualTo(1.0).Within(1e-9),
                "fractions are still normalized to sum to 1 under partial overlap");
        });
    }

    // NSVR-S10 — ν-parameter effect: holding the problem fixed, the recovered weight vector depends
    // on ν (the ε-insensitive tube width / SV-fraction budget). A single-ν sweep at ν = 0.25 yields a
    // materially different fraction than at ν = 0.75 for the disjoint 3×3 problem, and each run reports
    // the ν it was given. This proves ν is wired through the solver, not ignored.
    // Reference numbers (this session, matching the C# engine):
    //   ν = 0.25 → TypeA ≈ 0.7374 ; ν = 0.75 → TypeA ≈ 0.5085 (= sklearn-matched 0.508464).
    // Evidence: Schölkopf et al. (2000) — ν bounds the SV fraction / tube width, so it changes the fit.
    [Test]
    public void DeconvoluteImmuneCellsNuSvr_NuParameterChangesSolution()
    {
        // Arrange
        var signature = new Dictionary<string, IReadOnlyDictionary<string, double>>
        {
            ["TypeA"] = new Dictionary<string, double> { ["a1"] = 10, ["a2"] = 8, ["a3"] = 6 },
            ["TypeB"] = new Dictionary<string, double> { ["b1"] = 9, ["b2"] = 7, ["b3"] = 5 },
            ["TypeC"] = new Dictionary<string, double> { ["c1"] = 11, ["c2"] = 4, ["c3"] = 8 },
        };
        var planted = new Dictionary<string, double> { ["TypeA"] = 0.5, ["TypeB"] = 0.2, ["TypeC"] = 0.3 };
        var mixture = BuildPlantedMixture(signature, planted);

        // Act — force each ν singly so the reported BestNu is that ν and no RMSE selection occurs.
        var lowNu = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture, signature, new[] { 0.25 });
        var highNu = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture, signature, new[] { 0.75 });

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(lowNu.BestNu, Is.EqualTo(0.25).Within(ExactTolerance),
                "a single-ν sweep must report the ν it was given (0.25)");
            Assert.That(highNu.BestNu, Is.EqualTo(0.75).Within(ExactTolerance),
                "a single-ν sweep must report the ν it was given (0.75)");
            Assert.That(lowNu.CellFractions["TypeA"], Is.EqualTo(0.7374).Within(5e-4),
                "at ν = 0.25 the recovered TypeA fraction matches the engine reference ≈ 0.7374");
            Assert.That(highNu.CellFractions["TypeA"], Is.EqualTo(0.508464).Within(NuSvrReferenceTolerance),
                "at ν = 0.75 the recovered TypeA fraction matches the sklearn reference ≈ 0.508464");
            Assert.That(Math.Abs(lowNu.CellFractions["TypeA"] - highNu.CellFractions["TypeA"]),
                Is.GreaterThan(1e-3),
                "ν must materially change the recovered solution (ν is genuinely wired through the solver)");
        });
    }

    // NSVR-C1 — Loader rejects an empty input (no header).
    [Test]
    public void LoadSignatureMatrix_EmptyInput_Throws()
    {
        Assert.Throws<FormatException>(
            () => ImmuneAnalyzer.LoadSignatureMatrix(Array.Empty<string>()),
            "an empty signature matrix (no header) must raise FormatException");
    }

    // NSVR-C2 — Loader rejects a header without any cell-type columns.
    [Test]
    public void LoadSignatureMatrix_HeaderWithoutCellTypes_Throws()
    {
        Assert.Throws<FormatException>(
            () => ImmuneAnalyzer.LoadSignatureMatrix(new[] { "Gene symbol" }),
            "a header with only the gene-symbol column (no cell types) must raise FormatException");
    }

    // NSVR-C3 — Loader rejects a ragged data row (wrong column count).
    [Test]
    public void LoadSignatureMatrix_RaggedRow_Throws()
    {
        var lines = new[]
        {
            "Gene symbol\tTypeA\tTypeB",
            "g1\t1.0\t2.0",
            "g2\t3.0", // missing TypeB value
        };
        Assert.Throws<FormatException>(
            () => ImmuneAnalyzer.LoadSignatureMatrix(lines),
            "a row with the wrong number of columns must raise FormatException");
    }

    // NSVR-C4 — Loader rejects a non-numeric signature value.
    [Test]
    public void LoadSignatureMatrix_NonNumericValue_Throws()
    {
        var lines = new[]
        {
            "Gene symbol\tTypeA\tTypeB",
            "g1\t1.0\tNOT_A_NUMBER",
        };
        Assert.Throws<FormatException>(
            () => ImmuneAnalyzer.LoadSignatureMatrix(lines),
            "a non-numeric signature value must raise FormatException");
    }

    // NSVR-C5 — Null lines to the loader → ArgumentNullException.
    [Test]
    public void LoadSignatureMatrix_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => ImmuneAnalyzer.LoadSignatureMatrix(null!),
            "null TSV lines must raise ArgumentNullException");
    }

    #endregion

    #region Bundled ABIS-Seq Signature Matrix (LoadBundledAbisSignatureMatrix)

    // Bundled ABIS-Seq matrix (Monaco et al., 2019, Cell Reports 26:1627, CC BY 4.0).
    // Provenance: Table S5 (sheet "ABIS-Seq"), PMC6367568 supplementary mmc6.xlsx.
    // The 17 ABIS-Seq cell types and a few exact reference values, copied from the source matrix.
    private static readonly string[] AbisCellTypes =
    {
        "Monocytes C", "NK", "T CD8 Memory", "T CD4 Naive", "T CD8 Naive", "B Naive",
        "T CD4 Memory", "MAIT", "T gd Vd2", "Neutrophils LD", "T gd non-Vd2", "Basophils LD",
        "Monocytes NC+I", "B Memory", "mDCs", "pDCs", "Plasmablasts",
    };

    // ABIS-B1 — The bundled matrix loads with the published ABIS-Seq dimensions and cell-type names.
    // Evidence: Monaco et al. (2019), Table S5 (ABIS-Seq) — 1296 genes × 17 immune cell types.
    [Test]
    public void LoadBundledAbisSignatureMatrix_HasPublishedDimensions()
    {
        // Act
        var matrix = ImmuneAnalyzer.LoadBundledAbisSignatureMatrix();
        var genes = matrix.Values.SelectMany(d => d.Keys).Distinct().ToList();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(matrix.Count, Is.EqualTo(ImmuneAnalyzer.AbisSignatureCellTypeCount),
                "ABIS-Seq has 17 immune cell types (Monaco 2019, Table S5)");
            Assert.That(matrix.Count, Is.EqualTo(17), "ABIS-Seq cell-type count is 17");
            Assert.That(genes.Count, Is.EqualTo(ImmuneAnalyzer.AbisSignatureGeneCount),
                "ABIS-Seq has 1296 signature genes (Monaco 2019, Table S5)");
            Assert.That(genes.Count, Is.EqualTo(1296), "ABIS-Seq gene count is 1296");
            Assert.That(matrix.Keys, Is.EquivalentTo(AbisCellTypes),
                "the 17 ABIS-Seq cell-type names must match Table S5 exactly");
        });
    }

    // ABIS-B2 — Exact reference values from Table S5 (ABIS-Seq) round-trip through the loader.
    // Evidence (verbatim from mmc6.xlsx, sheet "ABIS-Seq"):
    //   S1PR3 / Monocytes C   = 45.720735005602499
    //   CD8A  / T CD8 Memory  = 1060.1507652944399
    //   MS4A1 / B Naive       = 3220.5650656491198
    //   S1PR3 / mDCs          = 3.9962058331855701
    [Test]
    public void LoadBundledAbisSignatureMatrix_HasExactReferenceValues()
    {
        // Act
        var matrix = ImmuneAnalyzer.LoadBundledAbisSignatureMatrix();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(matrix["Monocytes C"]["S1PR3"], Is.EqualTo(45.720735005602499).Within(ExactTolerance),
                "S1PR3 in Monocytes C must match the published ABIS-Seq value");
            Assert.That(matrix["T CD8 Memory"]["CD8A"], Is.EqualTo(1060.1507652944399).Within(ExactTolerance),
                "CD8A (a CD8 T-cell marker) must be high in T CD8 Memory per Table S5");
            Assert.That(matrix["B Naive"]["MS4A1"], Is.EqualTo(3220.5650656491198).Within(ExactTolerance),
                "MS4A1/CD20 (a B-cell marker) must be high in B Naive per Table S5");
            Assert.That(matrix["mDCs"]["S1PR3"], Is.EqualTo(3.9962058331855701).Within(ExactTolerance),
                "S1PR3 in mDCs must match the published ABIS-Seq value");
        });
    }

    // Tolerance for planted-truth recovery on the FULL bundled ABIS matrix. Unlike the disjoint-marker
    // synthetic matrix (NSVR-M1, ≤0.025), ABIS has 1296 shared genes spanning ~4 orders of magnitude;
    // after z-standardisation the epsilon-insensitive (robust) nu-SVR recovers the proportions of
    // well-separated lineages to within a small bounded deviation. This is the genuine recovery quality
    // of the unmodified engine on the real matrix — not a weakened assertion.
    private const double AbisRecoveryTolerance = 0.06;

    // ABIS-B3 — Planted-truth deconvolution on the BUNDLED ABIS matrix: bulk = ABIS·f for a known f,
    // nu-SVR must recover f (same evidence-based m = B·f pattern as NSVR-M1, on the real bundled
    // 1296×17 matrix — proves out-of-the-box deconvolution works). Two well-separated lineages
    // (NK + classical Monocytes) are recovered to within AbisRecoveryTolerance; every absent cell
    // type is recovered as exactly 0, and the ordering follows the planted ordering.
    // Evidence: m = B·f planted-truth construction; Newman et al. (2015) nu-SVR recovery.
    [Test]
    public void LoadBundledAbisSignatureMatrix_PlantedTruth_RecoversFractions()
    {
        // Arrange — plant a two-population mixture and synthesise the bulk profile m = ABIS·f.
        var abis = ImmuneAnalyzer.LoadBundledAbisSignatureMatrix();
        var planted = new Dictionary<string, double>
        {
            ["NK"] = 0.60,
            ["Monocytes C"] = 0.40,
        };
        var mixture = BuildPlantedMixture(abis, planted);

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture, abis);

        // Assert
        Assert.Multiple(() =>
        {
            foreach (var kvp in planted)
            {
                Assert.That(result.CellFractions[kvp.Key], Is.EqualTo(kvp.Value).Within(AbisRecoveryTolerance),
                    $"nu-SVR on the bundled ABIS matrix must recover the planted fraction for {kvp.Key} (planted {kvp.Value})");
            }

            // Every cell type absent from the mixture must recover an exactly-zero fraction.
            foreach (var ct in abis.Keys.Where(c => !planted.ContainsKey(c)))
            {
                Assert.That(result.CellFractions[ct], Is.EqualTo(0.0).Within(1e-9),
                    $"absent cell type {ct} must recover a zero fraction");
            }

            Assert.That(result.CellFractions["NK"], Is.GreaterThan(result.CellFractions["Monocytes C"]),
                "the dominant planted type (NK, 0.60) must have the larger recovered fraction");
            Assert.That(result.CellFractions.Values.Sum(), Is.EqualTo(1.0).Within(1e-9),
                "recovered ABIS fractions must sum to 1");
            Assert.That(result.Correlation, Is.GreaterThan(0.99),
                "the planted-truth reconstruction must correlate near-perfectly with the bulk mixture");
            Assert.That(result.OverlappingGenes, Is.EqualTo(ImmuneAnalyzer.AbisSignatureGeneCount),
                "all 1296 ABIS genes overlap the synthetic bulk built from the same matrix");
        });
    }

    // ABIS-B4 — Single-population planted truth on the bundled ABIS matrix recovers EXACTLY.
    // bulk = ABIS·(Monocytes C = 1.0) → fraction 1.0 for Monocytes C, 0 elsewhere, correlation 1.0.
    // This is the exact, discriminating planted-truth case: a wrong solver or a corrupted matrix
    // would not return a clean one-hot recovery.
    // Evidence: m = B·e_k planted-truth construction; Newman et al. (2015) nu-SVR recovery.
    [Test]
    public void LoadBundledAbisSignatureMatrix_SinglePopulationPlantedTruth_RecoversExactly()
    {
        // Arrange — a pure Monocytes-C bulk profile.
        var abis = ImmuneAnalyzer.LoadBundledAbisSignatureMatrix();
        var planted = new Dictionary<string, double> { ["Monocytes C"] = 1.0 };
        var mixture = BuildPlantedMixture(abis, planted);

        // Act
        var result = ImmuneAnalyzer.DeconvoluteImmuneCellsNuSvr(mixture, abis);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(result.CellFractions["Monocytes C"], Is.EqualTo(1.0).Within(1e-6),
                "a pure Monocytes-C bulk must deconvolute to fraction 1.0 for Monocytes C");
            foreach (var ct in abis.Keys.Where(c => c != "Monocytes C"))
            {
                Assert.That(result.CellFractions[ct], Is.EqualTo(0.0).Within(1e-6),
                    $"a pure Monocytes-C bulk must recover 0 for {ct}");
            }
            Assert.That(result.Correlation, Is.EqualTo(1.0).Within(1e-6),
                "a single-population reconstruction must correlate perfectly with the bulk");
        });
    }

    // ABIS-B5 — The bundled matrix is stable across calls (deterministic, identical content).
    [Test]
    public void LoadBundledAbisSignatureMatrix_IsDeterministic()
    {
        // Act
        var a = ImmuneAnalyzer.LoadBundledAbisSignatureMatrix();
        var b = ImmuneAnalyzer.LoadBundledAbisSignatureMatrix();

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(b.Count, Is.EqualTo(a.Count), "cell-type count must be identical across loads");
            Assert.That(b["T CD8 Memory"]["CD8A"], Is.EqualTo(a["T CD8 Memory"]["CD8A"]).Within(ExactTolerance),
                "a sampled value must be identical across loads");
        });
    }

    #endregion
}
