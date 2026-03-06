using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

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
}
