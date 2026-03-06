using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Oncology;

/// <summary>
/// Provides immune infiltration estimation algorithms for tumor samples.
/// <list type="bullet">
/// <item><description>
/// <b>EstimateInfiltration</b>: ssGSEA-based immune/stromal enrichment scoring and
/// tumor purity estimation per Yoshihara et al. (2013) and Barbie et al. (2009).
/// </description></item>
/// <item><description>
/// <b>DeconvoluteImmuneCells</b>: NNLS-based immune cell type deconvolution per
/// Lawson &amp; Hanson (1995) and Abbas et al. (2009, PLoS One 4:e6098).
/// </description></item>
/// </list>
/// Both methods accept configurable gene sets and signature matrices.
/// </summary>
public static class ImmuneAnalyzer
{
    #region Constants

    /// <summary>
    /// ssGSEA weighting exponent τ for hit accumulation in the enrichment score walk.
    /// Source: Barbie et al. (2009), Nature 462:108–112; GSVA default (Hänzelmann et al., 2013, BMC Bioinformatics 14:7).
    /// </summary>
    public const double SsGseaTau = 0.25;

    /// <summary>
    /// Coefficient for ESTIMATE tumor purity formula: cos(a + b × ESTIMATE_score).
    /// Source: Yoshihara et al. (2013), Nature Communications 4:2612.
    /// </summary>
    public const double EstimatePurityCoefficientA = 0.6049872018;

    /// <summary>
    /// Coefficient for ESTIMATE tumor purity formula: cos(a + b × ESTIMATE_score).
    /// Source: Yoshihara et al. (2013), Nature Communications 4:2612.
    /// </summary>
    public const double EstimatePurityCoefficientB = 0.0001467884;

    /// <summary>
    /// Default immune signature genes: the complete 141-gene ESTIMATE immune signature.
    /// Source: Yoshihara et al. (2013), Nature Communications 4:2612, Supplementary Data 1.
    /// Extracted from: ESTIMATE R package v1.0.11, inst/extdata/SI_geneset.gmt (ImmuneSignature).
    /// Custom gene sets can be supplied via the <c>immuneGenes</c> parameter of <see cref="EstimateInfiltration"/>.
    /// </summary>
    public static readonly IReadOnlyList<string> DefaultImmuneSignatureGenes = new[]
    {
        "LCP2", "LSP1", "FYB", "PLEK", "HCK", "IL10RA", "LILRB1", "NCKAP1L",
        "LAIR1", "NCF2", "CYBB", "PTPRC", "IL7R", "LAPTM5", "CD53", "EVI2B",
        "SLA", "ITGB2", "GIMAP4", "MYO1F", "HCLS1", "MNDA", "IL2RG", "CD48",
        "AOAH", "CCL5", "LTB", "GMFG", "GIMAP6", "GZMK", "LST1", "GPR65",
        "LILRB2", "WIPF1", "CD37", "BIN2", "FCER1G", "IKZF1", "TYROBP", "FGL2",
        "FLI1", "IRF8", "ARHGAP15", "SH2B3", "TNFRSF1B", "DOCK2", "CD2", "ARHGEF6",
        "CORO1A", "LY96", "LYZ", "ITGAL", "TNFAIP3", "RNASE6", "TGFB1", "PSTPIP1",
        "CST7", "RGS1", "FGR", "SELL", "MICAL1", "TRAF3IP3", "ITGA4", "MAFB",
        "ARHGDIB", "IL4R", "RHOH", "HLA-DPA1", "NKG7", "NCF4", "LPXN", "ITK",
        "SELPLG", "HLA-DPB1", "CD3D", "CD300A", "IL2RB", "ADCY7", "PTGER4", "SRGN",
        "CD247", "CCR7", "MSN", "ALOX5AP", "PTGER2", "RAC2", "GBP2", "VAV1",
        "CLEC2B", "P2RY14", "NFKBIA", "S100A9", "IFI30", "MFSD1", "RASSF2", "TPP1",
        "RHOG", "CLEC4A", "GZMB", "PVRIG", "S100A8", "CASP1", "BCL2A1", "HLA-E",
        "KLRB1", "GNLY", "RAB27A", "IL18RAP", "TPST2", "EMP3", "GMIP", "LCK",
        "IL32", "PTPRCAP", "LGALS9", "CCDC69", "SAMHD1", "TAP1", "GBP1", "CTSS",
        "GZMH", "ADAM8", "GLRX", "PRF1", "CD69", "HLA-B", "HLA-DMA", "CD74",
        "KLRK1", "PTPRE", "HLA-DRA", "VNN2", "TCIRG1", "RABGAP1L", "CSTA", "ZAP70",
        "HLA-F", "HLA-G", "CD52", "CD302", "CD27"
    };

    /// <summary>
    /// Default stromal signature genes: the complete 141-gene ESTIMATE stromal signature.
    /// Source: Yoshihara et al. (2013), Nature Communications 4:2612, Supplementary Data 1.
    /// Extracted from: ESTIMATE R package v1.0.11, inst/extdata/SI_geneset.gmt (StromalSignature).
    /// Custom gene sets can be supplied via the <c>stromalGenes</c> parameter of <see cref="EstimateInfiltration"/>.
    /// </summary>
    public static readonly IReadOnlyList<string> DefaultStromalSignatureGenes = new[]
    {
        "DCN", "PAPPA", "SFRP4", "THBS2", "LY86", "CXCL14", "FOXF1", "COL10A1",
        "ACTG2", "APBB1IP", "SH2D1A", "SULF1", "MSR1", "C3AR1", "FAP", "PTGIS",
        "ITGBL1", "BGN", "CXCL12", "ECM2", "FCGR2A", "MS4A4A", "WISP1", "COL1A2",
        "MS4A6A", "EDNRA", "VCAM1", "GPR124", "SCUBE2", "AIF1", "HEPH", "LUM",
        "PTGER3", "RUNX1T1", "CDH5", "PIK3R5", "RAMP3", "LDB2", "COX7A1", "EDIL3",
        "DDR2", "FCGR2B", "LPPR4", "COL15A1", "AOC3", "ITIH3", "FMO1", "PRKG1",
        "PLXDC1", "VSIG4", "COL6A3", "SGCD", "COL3A1", "F13A1", "OLFML1", "IGSF6",
        "COMP", "HGF", "GIMAP5", "ABCA6", "ITGAM", "MAF", "ITM2A", "CLEC7A",
        "ASPN", "LRRC15", "ERG", "CD86", "TRAT1", "COL8A2", "TCF21", "CD93",
        "CD163", "GREM1", "LMOD1", "TLR2", "ZEB2", "C1QB", "KCNJ8", "KDR",
        "CD33", "RASGRP3", "TNFSF4", "CCR1", "CSF1R", "BTK", "MFAP5", "MXRA5",
        "ISLR", "ARHGAP28", "ZFPM2", "TLR7", "ADAM12", "OLFML2B", "ENPP2", "CILP",
        "SIGLEC1", "SPON2", "PLXNC1", "ADAMTS5", "SAMSN1", "CH25H", "COL14A1", "EMCN",
        "RGS4", "PCDH12", "RARRES2", "CD248", "PDGFRB", "C1QA", "COL5A3", "IGF1",
        "SP140", "TFEC", "TNN", "ATP8B4", "ZNF423", "FRZB", "SERPING1", "ENPEP",
        "CD14", "DIO2", "FPR1", "IL18R1", "HDC", "TXNDC3", "PDE2A", "RSAD2",
        "ITIH5", "FASLG", "MMP3", "NOX4", "WNT2", "LRRC32", "CXCL9", "ODZ4",
        "FBLN2", "EGFL6", "IL1B", "SPON1", "CD200"
    };

    /// <summary>
    /// Default signature matrix for 22 immune cell types with curated marker gene expression weights.
    /// The 22 cell phenotypes correspond to the LM22 reference (Newman et al., 2015).
    /// Per-cell-type genes are representative markers from literature, not the full 547-gene LM22 matrix.
    /// A custom signature matrix (including the full LM22) can be supplied via the
    /// <c>signatureMatrix</c> parameter of <see cref="DeconvoluteImmuneCells"/>.
    /// Source: Newman et al. (2015), Nature Methods 12(5):453-457.
    /// </summary>
    public static readonly IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>> DefaultSignatureMatrix =
        new Dictionary<string, IReadOnlyDictionary<string, double>>
        {
            ["B_cells_naive"] = new Dictionary<string, double>
            {
                ["CD19"] = 8.5,
                ["MS4A1"] = 9.2,
                ["CD79A"] = 8.8,
                ["PAX5"] = 7.1,
                ["BANK1"] = 6.5
            },
            ["B_cells_memory"] = new Dictionary<string, double>
            {
                ["CD27"] = 7.4,
                ["CD38"] = 5.2,
                ["MS4A1"] = 8.1,
                ["AIM2"] = 5.8,
                ["CD19"] = 7.9
            },
            ["Plasma_cells"] = new Dictionary<string, double>
            {
                ["CD38"] = 10.2,
                ["SDC1"] = 9.1,
                ["XBP1"] = 8.7,
                ["MZB1"] = 9.5,
                ["IGKC"] = 11.0
            },
            ["T_cells_CD8"] = new Dictionary<string, double>
            {
                ["CD8A"] = 9.5,
                ["CD8B"] = 9.0,
                ["CD3D"] = 7.2,
                ["GZMB"] = 7.8,
                ["PRF1"] = 6.5
            },
            ["T_cells_CD4_naive"] = new Dictionary<string, double>
            {
                ["CD4"] = 7.1,
                ["CCR7"] = 8.2,
                ["SELL"] = 7.5,
                ["LEF1"] = 8.0,
                ["TCF7"] = 7.8
            },
            ["T_cells_CD4_memory_resting"] = new Dictionary<string, double>
            {
                ["CD4"] = 7.0,
                ["IL7R"] = 7.5,
                ["CD27"] = 6.0,
                ["LDHB"] = 6.5,
                ["CD3D"] = 7.0
            },
            ["T_cells_CD4_memory_activated"] = new Dictionary<string, double>
            {
                ["CD4"] = 7.3,
                ["NKG7"] = 6.0,
                ["GZMA"] = 5.5,
                ["CCL5"] = 7.2,
                ["CD3D"] = 7.1
            },
            ["T_cells_follicular_helper"] = new Dictionary<string, double>
            {
                ["CXCR5"] = 8.5,
                ["ICOS"] = 7.0,
                ["BCL6"] = 6.5,
                ["SH2D1A"] = 7.2,
                ["CD200"] = 6.0
            },
            ["T_cells_regulatory"] = new Dictionary<string, double>
            {
                ["FOXP3"] = 9.0,
                ["IL2RA"] = 8.5,
                ["CTLA4"] = 7.8,
                ["TNFRSF18"] = 7.0,
                ["IKZF2"] = 6.5
            },
            ["T_cells_gamma_delta"] = new Dictionary<string, double>
            {
                ["TRDC"] = 9.5,
                ["TRGC1"] = 9.0,
                ["TRGC2"] = 8.5,
                ["NKG7"] = 6.5,
                ["GZMB"] = 6.0
            },
            ["NK_cells_resting"] = new Dictionary<string, double>
            {
                ["NCAM1"] = 8.0,
                ["KLRD1"] = 7.5,
                ["KLRF1"] = 7.0,
                ["KIR2DL3"] = 6.0,
                ["KIR3DL1"] = 5.5
            },
            ["NK_cells_activated"] = new Dictionary<string, double>
            {
                ["NCAM1"] = 8.5,
                ["GZMB"] = 9.0,
                ["PRF1"] = 8.5,
                ["KLRD1"] = 7.5,
                ["IFNG"] = 7.0
            },
            ["Monocytes"] = new Dictionary<string, double>
            {
                ["CD14"] = 9.5,
                ["CSF1R"] = 8.0,
                ["FCN1"] = 8.5,
                ["S100A8"] = 9.0,
                ["VCAN"] = 7.0
            },
            ["Macrophages_M0"] = new Dictionary<string, double>
            {
                ["CD68"] = 8.5,
                ["CD163"] = 5.0,
                ["MSR1"] = 6.5,
                ["CSF1R"] = 7.0,
                ["MARCO"] = 5.5
            },
            ["Macrophages_M1"] = new Dictionary<string, double>
            {
                ["CD68"] = 8.0,
                ["NOS2"] = 7.5,
                ["IL1B"] = 8.0,
                ["TNF"] = 7.0,
                ["CXCL9"] = 7.5
            },
            ["Macrophages_M2"] = new Dictionary<string, double>
            {
                ["CD68"] = 8.0,
                ["CD163"] = 9.0,
                ["MRC1"] = 8.5,
                ["MSR1"] = 7.5,
                ["IL10"] = 7.0
            },
            ["Dendritic_cells_resting"] = new Dictionary<string, double>
            {
                ["ITGAX"] = 7.5,
                ["HLA-DRA"] = 8.0,
                ["CD1C"] = 7.0,
                ["FCER1A"] = 6.5,
                ["CLEC10A"] = 7.0
            },
            ["Dendritic_cells_activated"] = new Dictionary<string, double>
            {
                ["ITGAX"] = 8.0,
                ["HLA-DRA"] = 9.5,
                ["CD83"] = 8.0,
                ["CCR7"] = 7.5,
                ["LAMP3"] = 8.5
            },
            ["Mast_cells_resting"] = new Dictionary<string, double>
            {
                ["KIT"] = 8.5,
                ["CPA3"] = 9.0,
                ["MS4A2"] = 8.0,
                ["TPSAB1"] = 7.5,
                ["HDC"] = 7.0
            },
            ["Mast_cells_activated"] = new Dictionary<string, double>
            {
                ["KIT"] = 8.0,
                ["CPA3"] = 8.5,
                ["TPSAB1"] = 9.0,
                ["IL4"] = 6.5,
                ["IL13"] = 6.0
            },
            ["Eosinophils"] = new Dictionary<string, double>
            {
                ["SIGLEC8"] = 9.0,
                ["CCR3"] = 8.5,
                ["IL5RA"] = 7.5,
                ["PRG2"] = 8.0,
                ["EPX"] = 7.0
            },
            ["Neutrophils"] = new Dictionary<string, double>
            {
                ["CSF3R"] = 9.0,
                ["FCGR3B"] = 8.5,
                ["CEACAM8"] = 8.0,
                ["S100A12"] = 7.5,
                ["CXCR2"] = 7.0
            }
        };

    #endregion

    #region Records

    /// <summary>
    /// Result of ESTIMATE-style immune infiltration estimation.
    /// Source: Yoshihara et al. (2013), Nature Communications 4:2612.
    /// </summary>
    /// <param name="ImmuneScore">Enrichment score for immune signature genes. Higher values indicate greater immune infiltration.</param>
    /// <param name="StromalScore">Enrichment score for stromal signature genes. Higher values indicate greater stromal content.</param>
    /// <param name="EstimateScore">Combined score: ImmuneScore + StromalScore.</param>
    /// <param name="TumorPurity">Estimated tumor purity in [0, 1]. Computed as cos(a + b × EstimateScore).</param>
    /// <param name="OverlappingImmuneGenes">Number of immune signature genes found in the expression profile.</param>
    /// <param name="OverlappingStromalGenes">Number of stromal signature genes found in the expression profile.</param>
    public readonly record struct InfiltrationResult(
        double ImmuneScore,
        double StromalScore,
        double EstimateScore,
        double TumorPurity,
        int OverlappingImmuneGenes,
        int OverlappingStromalGenes);

    /// <summary>
    /// Result of NNLS-based immune cell type deconvolution.
    /// NNLS: Lawson &amp; Hanson (1995). Cell type framework: Newman et al. (2015).
    /// </summary>
    /// <param name="CellFractions">Dictionary mapping cell type names to their estimated fractions in [0, 1].</param>
    /// <param name="Correlation">Pearson correlation between observed and reconstructed expression profiles.</param>
    /// <param name="Rmse">Root mean square error of the deconvolution fit.</param>
    /// <param name="OverlappingGenes">Number of signature genes found in the expression profile.</param>
    public readonly record struct DeconvolutionResult(
        IReadOnlyDictionary<string, double> CellFractions,
        double Correlation,
        double Rmse,
        int OverlappingGenes);

    #endregion

    #region Immune Infiltration Estimation (ESTIMATE)

    /// <summary>
    /// Estimates immune and stromal infiltration from a gene expression profile
    /// using ssGSEA-based enrichment scoring per the ESTIMATE algorithm.
    /// <para>
    /// Algorithm: For each gene set (immune, stromal), computes a
    /// single-sample Gene Set Enrichment Analysis (ssGSEA) score:
    /// 1. Rank all genes in the expression profile by expression level (descending).
    /// 2. Walk through the ranked list; accumulate +hit (weighted by rank^τ, τ=0.25, where rank = N−i for position i in descending sort) for signature genes, −miss otherwise.
    /// 3. The enrichment score is the integral (sum) of the running sum across all positions.
    /// </para>
    /// <para>
    /// Source: Yoshihara et al. (2013), Nature Communications 4:2612.
    /// Source: Barbie et al. (2009), Nature 462:108-112 (ssGSEA method).
    /// Source: Hänzelmann et al. (2013), BMC Bioinformatics 14:7 (GSVA package, ssGSEA integral form).
    /// </para>
    /// </summary>
    /// <param name="expressionProfile">Dictionary mapping gene names to expression values (e.g., log2 TPM).</param>
    /// <param name="immuneGenes">Custom immune signature gene set. If null, uses <see cref="DefaultImmuneSignatureGenes"/>.</param>
    /// <param name="stromalGenes">Custom stromal signature gene set. If null, uses <see cref="DefaultStromalSignatureGenes"/>.</param>
    /// <returns>An <see cref="InfiltrationResult"/> with immune score, stromal score, ESTIMATE score, and tumor purity.</returns>
    public static InfiltrationResult EstimateInfiltration(
        IReadOnlyDictionary<string, double> expressionProfile,
        IReadOnlyList<string>? immuneGenes = null,
        IReadOnlyList<string>? stromalGenes = null)
    {
        ArgumentNullException.ThrowIfNull(expressionProfile);

        immuneGenes ??= DefaultImmuneSignatureGenes;
        stromalGenes ??= DefaultStromalSignatureGenes;

        if (expressionProfile.Count == 0)
        {
            return new InfiltrationResult(0.0, 0.0, 0.0, ComputeTumorPurity(0.0), 0, 0);
        }

        var immuneOverlap = immuneGenes.Where(g => expressionProfile.ContainsKey(g)).ToHashSet();
        var stromalOverlap = stromalGenes.Where(g => expressionProfile.ContainsKey(g)).ToHashSet();

        double immuneScore = ComputeSsGseaScore(expressionProfile, immuneOverlap);
        double stromalScore = ComputeSsGseaScore(expressionProfile, stromalOverlap);
        double estimateScore = immuneScore + stromalScore;
        double tumorPurity = ComputeTumorPurity(estimateScore);

        return new InfiltrationResult(
            immuneScore,
            stromalScore,
            estimateScore,
            tumorPurity,
            immuneOverlap.Count,
            stromalOverlap.Count);
    }

    #endregion

    #region Immune Cell Deconvolution (NNLS)

    /// <summary>
    /// Deconvolutes immune cell type proportions from a bulk expression profile
    /// using Non-Negative Least Squares (NNLS) regression.
    /// <para>
    /// Algorithm: Solves the linear system m = S × f where:
    /// - m is the mixture expression vector (observed gene expression),
    /// - S is the signature matrix (reference expression of each cell type),
    /// - f is the fraction vector (unknown cell type proportions).
    /// The NNLS algorithm finds f that minimizes ||m - S·f||² subject to f ≥ 0,
    /// then normalizes f so that sum(f) = 1.
    /// </para>
    /// <para>
    /// This is the NNLS/LLSR deconvolution approach as used by Abbas et al. (2009,
    /// PLoS One 4:e6098) — distinct from CIBERSORT which uses ν-SVR (Newman et al., 2015).
    /// NNLS: Lawson &amp; Hanson (1995), Solving Least Squares Problems, SIAM, Chapter 23.
    /// Cell type framework: Newman et al. (2015), Nature Methods 12(5):453-457.
    /// </para>
    /// </summary>
    /// <param name="expressionProfile">Dictionary mapping gene names to expression values.</param>
    /// <param name="signatureMatrix">
    /// Signature matrix: outer key = cell type name, inner dictionary = gene name → reference expression.
    /// If null, uses <see cref="DefaultSignatureMatrix"/>.
    /// </param>
    /// <param name="maxIterations">Maximum number of NNLS iterations (default: 1000).</param>
    /// <returns>A <see cref="DeconvolutionResult"/> with cell type fractions, correlation, RMSE, and overlap count.</returns>
    public static DeconvolutionResult DeconvoluteImmuneCells(
        IReadOnlyDictionary<string, double> expressionProfile,
        IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>>? signatureMatrix = null,
        int maxIterations = 1000)
    {
        ArgumentNullException.ThrowIfNull(expressionProfile);

        signatureMatrix ??= DefaultSignatureMatrix;

        // Collect all genes present in both signature and expression profile.
        var allSignatureGenes = signatureMatrix.Values
            .SelectMany(cellGenes => cellGenes.Keys)
            .Distinct()
            .ToList();

        var overlappingGenes = allSignatureGenes
            .Where(g => expressionProfile.ContainsKey(g))
            .ToList();

        var cellTypes = signatureMatrix.Keys.ToList();

        if (overlappingGenes.Count == 0 || cellTypes.Count == 0)
        {
            var emptyFractions = cellTypes.ToDictionary(ct => ct, _ => 0.0);
            return new DeconvolutionResult(
                emptyFractions,
                Correlation: 0.0,
                Rmse: 0.0,
                OverlappingGenes: 0);
        }

        // Build the mixture vector m and signature matrix S.
        int nGenes = overlappingGenes.Count;
        int nCellTypes = cellTypes.Count;

        double[] m = new double[nGenes];
        double[,] s = new double[nGenes, nCellTypes];

        for (int i = 0; i < nGenes; i++)
        {
            m[i] = expressionProfile[overlappingGenes[i]];
            for (int j = 0; j < nCellTypes; j++)
            {
                if (signatureMatrix[cellTypes[j]].TryGetValue(overlappingGenes[i], out double val))
                {
                    s[i, j] = val;
                }
            }
        }

        // Solve using NNLS (Lawson-Hanson algorithm).
        double[] fractions = SolveNnls(s, m, nGenes, nCellTypes, maxIterations);

        // Normalize fractions to sum to 1.
        double sum = fractions.Sum();
        if (sum > 0)
        {
            for (int j = 0; j < nCellTypes; j++)
            {
                fractions[j] /= sum;
            }
        }

        // Compute reconstructed expression and quality metrics.
        double[] reconstructed = new double[nGenes];
        for (int i = 0; i < nGenes; i++)
        {
            for (int j = 0; j < nCellTypes; j++)
            {
                reconstructed[i] += s[i, j] * fractions[j];
            }
        }

        double correlation = ComputePearsonCorrelation(m, reconstructed);
        double rmse = ComputeRmse(m, reconstructed);

        var cellFractions = new Dictionary<string, double>();
        for (int j = 0; j < nCellTypes; j++)
        {
            cellFractions[cellTypes[j]] = fractions[j];
        }

        return new DeconvolutionResult(
            cellFractions,
            correlation,
            rmse,
            overlappingGenes.Count);
    }

    #endregion

    #region Private Helpers

    /// <summary>
    /// Computes ssGSEA enrichment score for a gene set in an expression profile.
    /// Returns the integral (sum) of the weighted running sum across all ranked positions.
    /// Hit weights use rank^τ where τ = <see cref="SsGseaTau"/> (0.25) and rank = N − i
    /// for gene at sorted position i (0-indexed, descending by expression).
    /// <para>
    /// Source: Barbie et al. (2009), Nature 462:108-112.
    /// Source: Hänzelmann et al. (2013), BMC Bioinformatics 14:7 (GSVA ssGSEA implementation).
    /// The GSVA ssGSEA function weights hits by rank (integer position), not by expression value.
    /// </para>
    /// </summary>
    private static double ComputeSsGseaScore(
        IReadOnlyDictionary<string, double> expressionProfile,
        HashSet<string> geneSet)
    {
        if (geneSet.Count == 0 || expressionProfile.Count == 0)
            return 0.0;

        // Rank genes by expression (descending).
        var rankedGenes = expressionProfile
            .OrderByDescending(kvp => kvp.Value)
            .Select(kvp => kvp.Key)
            .ToList();

        int n = rankedGenes.Count;
        int nHits = geneSet.Count;
        int nMiss = n - nHits;

        if (nHits == 0 || nMiss == 0)
            return 0.0;

        // Compute total hit weight: sum of rank^τ for genes in the set.
        // rank = N − i where i is the gene's zero-indexed position in descending-expression order.
        // Highest expression → highest rank = N. Lowest expression → rank = 1.
        // τ = 0.25 per Barbie et al. (2009) and GSVA default (Hänzelmann et al., 2013).
        double totalHitWeight = 0.0;
        for (int i = 0; i < n; i++)
        {
            if (geneSet.Contains(rankedGenes[i]))
            {
                int rank = n - i;
                totalHitWeight += Math.Pow(rank, SsGseaTau);
            }
        }

        if (totalHitWeight == 0.0)
            return 0.0;

        // Compute ssGSEA enrichment score as the integral (sum) of the running sum.
        // This is the area under the enrichment curve — distinct from classic GSEA which uses max deviation.
        // Source: Hänzelmann et al. (2013), GSVA package ssGSEA implementation.
        double runningSum = 0.0;
        double integral = 0.0;
        double missStep = 1.0 / nMiss;

        for (int i = 0; i < n; i++)
        {
            if (geneSet.Contains(rankedGenes[i]))
            {
                int rank = n - i;
                double hitWeight = Math.Pow(rank, SsGseaTau) / totalHitWeight;
                runningSum += hitWeight;
            }
            else
            {
                runningSum -= missStep;
            }

            integral += runningSum;
        }

        return integral;
    }

    /// <summary>
    /// Computes tumor purity using the ESTIMATE formula.
    /// Source: Yoshihara et al. (2013), purity = cos(a + b × estimateScore), clamped to [0, 1].
    /// </summary>
    private static double ComputeTumorPurity(double estimateScore)
    {
        double purity = Math.Cos(EstimatePurityCoefficientA + EstimatePurityCoefficientB * estimateScore);
        return Math.Clamp(purity, 0.0, 1.0);
    }

    /// <summary>
    /// Non-negative least squares (NNLS) solver using the Lawson-Hanson active set method.
    /// Finds x that minimizes ||Ax - b||² subject to x ≥ 0.
    /// Source: Lawson &amp; Hanson (1995), Solving Least Squares Problems, SIAM, Chapter 23.
    /// </summary>
    private static double[] SolveNnls(double[,] a, double[] b, int nRows, int nCols, int maxIterations)
    {
        double[] x = new double[nCols];
        double[] w = new double[nCols]; // gradient
        bool[] passiveSet = new bool[nCols]; // true = variable is in the passive (free) set

        // Compute initial gradient: w = A^T * (b - A*x)
        ComputeGradient(a, b, x, w, nRows, nCols);

        int iteration = 0;
        while (iteration < maxIterations)
        {
            // Find the maximum gradient among variables in the active (zero) set.
            int maxIdx = -1;
            double maxW = 0.0;
            for (int j = 0; j < nCols; j++)
            {
                if (!passiveSet[j] && w[j] > maxW)
                {
                    maxW = w[j];
                    maxIdx = j;
                }
            }

            // If no positive gradient in active set, solution is optimal.
            if (maxIdx == -1 || maxW <= 1e-15)
                break;

            // Move variable to passive set.
            passiveSet[maxIdx] = true;

            // Inner loop: solve least squares on passive set, fix infeasible variables.
            bool feasible = false;
            while (!feasible)
            {
                // Solve least squares for passive variables only.
                double[] z = SolvePassiveSetLeastSquares(a, b, passiveSet, nRows, nCols);

                // Check feasibility.
                feasible = true;
                for (int j = 0; j < nCols; j++)
                {
                    if (passiveSet[j] && z[j] <= 0)
                    {
                        feasible = false;
                        break;
                    }
                }

                if (feasible)
                {
                    Array.Copy(z, x, nCols);
                }
                else
                {
                    // Find alpha to stay feasible.
                    double alpha = double.MaxValue;
                    for (int j = 0; j < nCols; j++)
                    {
                        if (passiveSet[j] && z[j] <= 0)
                        {
                            double ratio = x[j] / (x[j] - z[j]);
                            if (ratio < alpha)
                                alpha = ratio;
                        }
                    }

                    // Update x.
                    for (int j = 0; j < nCols; j++)
                    {
                        x[j] += alpha * (z[j] - x[j]);
                        if (passiveSet[j] && Math.Abs(x[j]) < 1e-15)
                        {
                            x[j] = 0.0;
                            passiveSet[j] = false;
                        }
                    }
                }

                iteration++;
                if (iteration >= maxIterations)
                    break;
            }

            ComputeGradient(a, b, x, w, nRows, nCols);
        }

        // Ensure all values are non-negative (numerical cleanup).
        for (int j = 0; j < nCols; j++)
        {
            if (x[j] < 0)
                x[j] = 0.0;
        }

        return x;
    }

    /// <summary>
    /// Computes gradient w = A^T * (b - A*x).
    /// </summary>
    private static void ComputeGradient(double[,] a, double[] b, double[] x, double[] w, int nRows, int nCols)
    {
        // residual = b - A*x
        double[] residual = new double[nRows];
        for (int i = 0; i < nRows; i++)
        {
            residual[i] = b[i];
            for (int j = 0; j < nCols; j++)
            {
                residual[i] -= a[i, j] * x[j];
            }
        }

        // w = A^T * residual
        for (int j = 0; j < nCols; j++)
        {
            w[j] = 0.0;
            for (int i = 0; i < nRows; i++)
            {
                w[j] += a[i, j] * residual[i];
            }
        }
    }

    /// <summary>
    /// Solves least squares for variables in the passive set only.
    /// Uses normal equations: (A_p^T * A_p) * z_p = A_p^T * b.
    /// </summary>
    private static double[] SolvePassiveSetLeastSquares(double[,] a, double[] b, bool[] passiveSet, int nRows, int nCols)
    {
        double[] z = new double[nCols];
        var passiveIndices = new List<int>();
        for (int j = 0; j < nCols; j++)
        {
            if (passiveSet[j])
                passiveIndices.Add(j);
        }

        int nPassive = passiveIndices.Count;
        if (nPassive == 0)
            return z;

        // Build A_p^T * A_p and A_p^T * b
        double[,] ata = new double[nPassive, nPassive];
        double[] atb = new double[nPassive];

        for (int pi = 0; pi < nPassive; pi++)
        {
            int ji = passiveIndices[pi];
            for (int pk = 0; pk < nPassive; pk++)
            {
                int jk = passiveIndices[pk];
                double val = 0;
                for (int i = 0; i < nRows; i++)
                    val += a[i, ji] * a[i, jk];
                ata[pi, pk] = val;
            }

            double bval = 0;
            for (int i = 0; i < nRows; i++)
                bval += a[i, ji] * b[i];
            atb[pi] = bval;
        }

        // Solve via Gaussian elimination with partial pivoting.
        double[] solution = SolveLinearSystem(ata, atb, nPassive);

        for (int pi = 0; pi < nPassive; pi++)
        {
            z[passiveIndices[pi]] = solution[pi];
        }

        return z;
    }

    /// <summary>
    /// Solves a linear system Ax = b using Gaussian elimination with partial pivoting.
    /// </summary>
    private static double[] SolveLinearSystem(double[,] a, double[] b, int n)
    {
        // Create augmented matrix.
        double[,] aug = new double[n, n + 1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
                aug[i, j] = a[i, j];
            aug[i, n] = b[i];
        }

        // Forward elimination with partial pivoting.
        for (int k = 0; k < n; k++)
        {
            // Find pivot.
            int maxRow = k;
            double maxVal = Math.Abs(aug[k, k]);
            for (int i = k + 1; i < n; i++)
            {
                if (Math.Abs(aug[i, k]) > maxVal)
                {
                    maxVal = Math.Abs(aug[i, k]);
                    maxRow = i;
                }
            }

            // Swap rows.
            if (maxRow != k)
            {
                for (int j = k; j <= n; j++)
                {
                    (aug[k, j], aug[maxRow, j]) = (aug[maxRow, j], aug[k, j]);
                }
            }

            double pivot = aug[k, k];
            if (Math.Abs(pivot) < 1e-15)
                continue; // Singular or near-singular — skip.

            // Eliminate below.
            for (int i = k + 1; i < n; i++)
            {
                double factor = aug[i, k] / pivot;
                for (int j = k; j <= n; j++)
                {
                    aug[i, j] -= factor * aug[k, j];
                }
            }
        }

        // Back substitution.
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            if (Math.Abs(aug[i, i]) < 1e-15)
            {
                x[i] = 0;
                continue;
            }
            double sum = aug[i, n];
            for (int j = i + 1; j < n; j++)
            {
                sum -= aug[i, j] * x[j];
            }
            x[i] = sum / aug[i, i];
        }

        return x;
    }

    /// <summary>
    /// Computes Pearson correlation coefficient between two vectors.
    /// </summary>
    private static double ComputePearsonCorrelation(double[] x, double[] y)
    {
        int n = x.Length;
        if (n < 2)
            return 0.0;

        double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;
        for (int i = 0; i < n; i++)
        {
            sumX += x[i];
            sumY += y[i];
            sumXY += x[i] * y[i];
            sumX2 += x[i] * x[i];
            sumY2 += y[i] * y[i];
        }

        double numerator = n * sumXY - sumX * sumY;
        double denominator = Math.Sqrt((n * sumX2 - sumX * sumX) * (n * sumY2 - sumY * sumY));

        if (denominator < 1e-15)
            return 0.0;

        return numerator / denominator;
    }

    /// <summary>
    /// Computes root mean square error between two vectors.
    /// </summary>
    private static double ComputeRmse(double[] observed, double[] predicted)
    {
        int n = observed.Length;
        if (n == 0)
            return 0.0;

        double sumSqErr = 0;
        for (int i = 0; i < n; i++)
        {
            double diff = observed[i] - predicted[i];
            sumSqErr += diff * diff;
        }

        return Math.Sqrt(sumSqErr / n);
    }

    #endregion
}
