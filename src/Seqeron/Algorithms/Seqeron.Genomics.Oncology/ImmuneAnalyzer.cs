using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
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
    /// The set of ν (nu) values CIBERSORT sweeps for ν-support-vector-regression deconvolution.
    /// For each ν the linear-kernel ν-SVR is solved; the ν giving the lowest RMSE between the
    /// observed mixture <c>m</c> and the reconstruction <c>B·f</c> is selected.
    /// Source: Newman et al. (2015), Nature Methods 12(5):453-457; CIBERSORT protocol
    /// (Chen et al., 2018, Methods Mol Biol 1711:243-259) — "CIBERSORT uses a set of ν values
    /// (0.25, 0.5, 0.75) and chooses the value producing the best [lowest-RMSE] result."
    /// </summary>
    public static readonly IReadOnlyList<double> CibersortNuValues = new[] { 0.25, 0.5, 0.75 };

    /// <summary>
    /// Regularization constant C of the ν-SVR primal/dual (box bound on the dual variables,
    /// α_i, α_i* ∈ [0, C]). CIBERSORT/libsvm use the default C = 1. Source: Schölkopf et al. (2000),
    /// "New Support Vector Algorithms", Neural Computation 12(5):1207-1245 (ν-SVR formulation);
    /// libsvm default cost = 1 (Chang &amp; Lin, 2011).
    /// </summary>
    public const double NuSvrCost = 1.0;

    /// <summary>
    /// Number of immune cell types in the bundled ABIS-Seq signature matrix (Monaco et al., 2019).
    /// Source: Monaco G et al. (2019), Cell Reports 26(6):1627-1640.e7, Table S5 (ABIS-Seq sheet).
    /// </summary>
    public const int AbisSignatureCellTypeCount = 17;

    /// <summary>
    /// Number of genes (rows) in the bundled ABIS-Seq signature matrix (Monaco et al., 2019).
    /// Source: Monaco G et al. (2019), Cell Reports 26(6):1627-1640.e7, Table S5 (ABIS-Seq sheet).
    /// </summary>
    public const int AbisSignatureGeneCount = 1296;

    /// <summary>
    /// Embedded-resource name of the bundled ABIS-Seq signature matrix TSV.
    /// </summary>
    private const string AbisResourceName = "Seqeron.Genomics.Oncology.Resources.ABIS_sigmatrixRNAseq.tsv";

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

    /// <summary>
    /// Result of CIBERSORT-style ν-support-vector-regression (ν-SVR) immune cell deconvolution.
    /// Source: Newman et al. (2015), Nature Methods 12(5):453-457 (CIBERSORT);
    /// Schölkopf et al. (2000), Neural Computation 12(5):1207-1245 (ν-SVR).
    /// </summary>
    /// <param name="CellFractions">Cell-type name → estimated fraction in [0, 1], summing to 1 (after zero-clip and renormalisation).</param>
    /// <param name="BestNu">The ν value (from <see cref="CibersortNuValues"/>) selected as giving the lowest RMSE.</param>
    /// <param name="Correlation">Pearson correlation between observed mixture and reconstructed profile at the selected ν.</param>
    /// <param name="Rmse">Root mean square error of the fit at the selected ν.</param>
    /// <param name="OverlappingGenes">Number of signature genes found in the expression profile.</param>
    public readonly record struct NuSvrDeconvolutionResult(
        IReadOnlyDictionary<string, double> CellFractions,
        double BestNu,
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

    /// <summary>
    /// Converts a cohort-scaled ESTIMATE score into an <b>absolute</b> tumour-purity estimate
    /// using the closed-form transform published in the ESTIMATE paper.
    /// <para>
    /// <c>purity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)</c>
    /// </para>
    /// <para>
    /// This is an <b>opt-in</b> alternative to the relative <see cref="InfiltrationResult.TumorPurity"/>
    /// returned by <see cref="EstimateInfiltration"/>: that field applies the same cosine to the
    /// library's <i>single-sample, un-normalised</i> ssGSEA integral and is therefore a relative
    /// value, whereas this method applies the transform to a caller-supplied ESTIMATE score that is
    /// on the original ESTIMATE numeric scale (the sum of the cohort-/rank-normalised immune and
    /// stromal ssGSEA scores produced by the ESTIMATE R package).
    /// </para>
    /// <para>
    /// <b>Domain.</b> The regression curve was fit by nonlinear least squares against ABSOLUTE
    /// purity on TCGA <b>Affymetrix</b> data; it is calibrated for, and should only be applied to,
    /// Affymetrix-derived ESTIMATE scores (it is not valid for RNA-seq-derived scores). Following the
    /// reference implementation, when the cosine evaluates to a <b>negative</b> value the purity is
    /// undefined and this method returns <see cref="double.NaN"/> (the R <c>estimate</c>/
    /// <c>tidyestimate</c> packages set such values to <c>NA</c>); otherwise the returned purity lies
    /// in (0, 1].
    /// </para>
    /// <para>
    /// Source: Yoshihara et al. (2013), Nature Communications 4:2612, doi:10.1038/ncomms3612
    /// ("Tumor purity was estimated …  Tumor_purity = cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)";
    /// "regression curve … based on ABSOLUTE … by applying the nonlinear least squares method").
    /// Reference implementation: ESTIMATE R package; tidyestimate <c>estimate_score()</c>
    /// (<c>purity = cos(0.6049872018 + 0.0001467884 * estimate); purity = ifelse(purity &lt; 0, NA, purity)</c>).
    /// </para>
    /// </summary>
    /// <param name="estimateScore">
    /// The ESTIMATE score (immune + stromal enrichment) on the original ESTIMATE numeric scale,
    /// derived from Affymetrix expression data.
    /// </param>
    /// <returns>
    /// The absolute tumour purity in (0, 1], or <see cref="double.NaN"/> when the score lies outside
    /// the domain where the cosine model yields a non-negative purity.
    /// </returns>
    public static double EstimateTumorPurity(double estimateScore)
    {
        // ESTIMATE->ABSOLUTE transform calibrated on Affymetrix data (Yoshihara 2013) — platform-contingent
        // off-array. Strict mode throws; Moderate/Permissive allow it with that caveat.
        Seqeron.Genomics.Core.LimitationPolicy.Enforce("ONCO-IMMUNE-001");

        double purity = Math.Cos(EstimatePurityCoefficientA + EstimatePurityCoefficientB * estimateScore);

        // Reference implementation (ESTIMATE / tidyestimate): negative cosine values are out of the
        // calibrated domain and are reported as NA. We mirror that as NaN rather than clamping,
        // because a clamped 0 would falsely imply a valid (fully non-tumour) estimate.
        return purity < 0.0 ? double.NaN : purity;
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

        // Deconvolution against the bundled ABIS (or caller) matrix — NOT CIBERSORT-LM22-identical.
        // Strict mode throws; Moderate/Permissive allow it.
        Seqeron.Genomics.Core.LimitationPolicy.Enforce("ONCO-IMMUNE-001");

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

    #region Immune Cell Deconvolution (CIBERSORT ν-SVR)

    /// <summary>
    /// Deconvolutes immune cell-type proportions from a bulk expression profile using
    /// <b>ν-support-vector regression</b> (ν-SVR) with a linear kernel — the regression engine of
    /// CIBERSORT (Newman et al., 2015).
    /// <para>
    /// For the mixture vector <c>m</c> (genes) and a signature matrix <c>B</c> (genes × cell types),
    /// a linear ν-SVR of <c>m</c> on the columns of <c>B</c> is solved at each ν in
    /// <see cref="CibersortNuValues"/> (0.25, 0.5, 0.75). For each ν the regression weight vector
    /// <c>f</c> = <c>w</c> is recovered; the ν giving the lowest RMSE between <c>m</c> and
    /// <c>B·f</c> is selected. Negative weights are set to 0 and the remaining weights are
    /// normalised to sum to 1, giving cell-type fractions.
    /// </para>
    /// <para>
    /// The ν-SVR optimisation follows Schölkopf et al. (2000): minimise
    /// <c>½‖w‖² + C(Σ(ξ_i+ξ_i*)/ℓ + νε)</c> subject to the ε-insensitive tube constraints, whose
    /// linear-kernel dual is maximise
    /// <c>−½ Σ(α_i−α_i*)(α_j−α_j*)⟨x_i,x_j⟩ + Σ y_i(α_i−α_i*)</c> subject to
    /// <c>Σ(α_i−α_i*)=0</c>, <c>Σ(α_i+α_i*) ≤ C·ν·ℓ</c>, <c>α_i, α_i* ∈ [0, C]</c>; here
    /// <c>w = Σ(α_i−α_i*)x_i</c>. Per Theorem 9 of Schölkopf et al. (2000), ν is an upper bound on
    /// the fraction of points outside the tube and a lower bound on the fraction of support vectors.
    /// </para>
    /// <para>
    /// <b>LM22.</b> CIBERSORT is normally run with the LM22 signature matrix (547 genes × 22 cell
    /// types). LM22 is distributed by Stanford under a non-commercial licence that forbids
    /// redistribution ("RECIPIENT shall not distribute the Program or transfer it to any other
    /// person or organization without prior written permission from STANFORD"), so it is <b>not</b>
    /// bundled with this library. Obtain <c>LM22.txt</c> from https://cibersort.stanford.edu under
    /// your own licence and supply it via <see cref="LoadSignatureMatrix"/>. If
    /// <paramref name="signatureMatrix"/> is null, the bundled representative 5-marker
    /// <see cref="DefaultSignatureMatrix"/> is used (a small synthetic matrix, not LM22).
    /// </para>
    /// <para>
    /// Source: Newman et al. (2015), Nature Methods 12(5):453-457, doi:10.1038/nmeth.3337.
    /// Source: Schölkopf, Smola, Williamson &amp; Bartlett (2000), Neural Computation 12(5):1207-1245,
    /// doi:10.1162/089976600300015565; Smola &amp; Schölkopf (2004) tutorial, eqs (60)-(62).
    /// </para>
    /// </summary>
    /// <param name="expressionProfile">Dictionary mapping gene names to expression values (mixture vector m).</param>
    /// <param name="signatureMatrix">
    /// Signature matrix B: outer key = cell type, inner dictionary = gene → reference expression.
    /// If null, uses <see cref="DefaultSignatureMatrix"/> (NOT LM22 — see remarks). Supply LM22 via
    /// <see cref="LoadSignatureMatrix"/>.
    /// </param>
    /// <param name="nuValues">ν values to sweep. If null, uses <see cref="CibersortNuValues"/> (0.25, 0.5, 0.75).</param>
    /// <returns>A <see cref="NuSvrDeconvolutionResult"/> with cell-type fractions, the selected ν, correlation, RMSE, and overlap count.</returns>
    public static NuSvrDeconvolutionResult DeconvoluteImmuneCellsNuSvr(
        IReadOnlyDictionary<string, double> expressionProfile,
        IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>>? signatureMatrix = null,
        IReadOnlyList<double>? nuValues = null)
    {
        ArgumentNullException.ThrowIfNull(expressionProfile);

        // Deconvolution against the bundled ABIS (or caller) matrix — NOT CIBERSORT-LM22-identical.
        // Strict mode throws; Moderate/Permissive allow it.
        Seqeron.Genomics.Core.LimitationPolicy.Enforce("ONCO-IMMUNE-001");

        signatureMatrix ??= DefaultSignatureMatrix;
        nuValues ??= CibersortNuValues;

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
            return new NuSvrDeconvolutionResult(emptyFractions, BestNu: 0.0, Correlation: 0.0, Rmse: 0.0, OverlappingGenes: 0);
        }

        int nGenes = overlappingGenes.Count;
        int nCellTypes = cellTypes.Count;

        // Build mixture vector m and signature matrix B (genes × cell types).
        double[] m = new double[nGenes];
        double[,] b = new double[nGenes, nCellTypes];
        for (int i = 0; i < nGenes; i++)
        {
            m[i] = expressionProfile[overlappingGenes[i]];
            for (int j = 0; j < nCellTypes; j++)
            {
                if (signatureMatrix[cellTypes[j]].TryGetValue(overlappingGenes[i], out double val))
                    b[i, j] = val;
            }
        }

        // CIBERSORT normalises both the mixture and the signature to zero mean / unit variance
        // before regression so the SVR scale is comparable across features and samples.
        // Source: Newman et al. (2015) — "the mixture and signature gene expression values are
        // standardized (z-score)". The recovered weights w of the standardized regression are the
        // cell-type coefficients; only their sign and relative magnitude (post-normalisation) matter.
        double[] mz = Standardize(m);
        double[,] bz = StandardizeColumns(b, nGenes, nCellTypes);

        // Solve the linear ν-SVR for each ν; keep the weight vector with the lowest reconstruction RMSE.
        double[]? bestWeights = null;
        double bestNu = 0.0;
        double bestRmse = double.PositiveInfinity;

        foreach (double nu in nuValues)
        {
            double[] w = SolveNuSvrLinear(bz, mz, nGenes, nCellTypes, nu, NuSvrCost);

            // Reconstruct on the standardized scale and measure RMSE (CIBERSORT's selection metric).
            double[] reconstructed = new double[nGenes];
            for (int i = 0; i < nGenes; i++)
            {
                double acc = 0.0;
                for (int j = 0; j < nCellTypes; j++)
                    acc += bz[i, j] * w[j];
                reconstructed[i] = acc;
            }

            double rmse = ComputeRmse(mz, reconstructed);
            if (rmse < bestRmse)
            {
                bestRmse = rmse;
                bestNu = nu;
                bestWeights = w;
            }
        }

        double[] fractions = bestWeights ?? new double[nCellTypes];

        // Zero-clip negative coefficients and normalise to sum 1.
        // Source: Newman et al. (2015) — negative weights are set to zero and the remaining
        // coefficients are normalised to sum to 1 to yield cell-type fractions.
        for (int j = 0; j < nCellTypes; j++)
        {
            if (fractions[j] < 0.0)
                fractions[j] = 0.0;
        }

        double sum = fractions.Sum();
        if (sum > 0.0)
        {
            for (int j = 0; j < nCellTypes; j++)
                fractions[j] /= sum;
        }

        // Quality metrics on the ORIGINAL (non-standardized) scale, against the final fractions.
        double[] recon = new double[nGenes];
        for (int i = 0; i < nGenes; i++)
        {
            double acc = 0.0;
            for (int j = 0; j < nCellTypes; j++)
                acc += b[i, j] * fractions[j];
            recon[i] = acc;
        }

        double correlation = ComputePearsonCorrelation(m, recon);
        double finalRmse = ComputeRmse(m, recon);

        var cellFractions = new Dictionary<string, double>();
        for (int j = 0; j < nCellTypes; j++)
            cellFractions[cellTypes[j]] = fractions[j];

        return new NuSvrDeconvolutionResult(cellFractions, bestNu, correlation, finalRmse, overlappingGenes.Count);
    }

    /// <summary>
    /// Parses a CIBERSORT LM22-format signature matrix from tab-separated text into the
    /// <c>cellType → (gene → value)</c> shape consumed by <see cref="DeconvoluteImmuneCellsNuSvr"/>
    /// and <see cref="DeconvoluteImmuneCells"/>.
    /// <para>
    /// Format (per the LM22 distribution): a tab-separated table whose first row is a header —
    /// the first column label (e.g. "Gene symbol") followed by the cell-type column names — and
    /// whose subsequent rows each give a gene symbol followed by one numeric value per cell type.
    /// The full LM22 is 547 genes × 22 cell types.
    /// </para>
    /// <para>
    /// <b>Licence.</b> The LM22 matrix itself is NOT bundled: Stanford distributes it under a
    /// non-commercial licence forbidding redistribution. Obtain <c>LM22.txt</c> from
    /// https://cibersort.stanford.edu under your own licence and pass its contents here.
    /// Source (format): Newman et al. (2015), Nature Methods 12(5):453-457, Supplementary Table.
    /// </para>
    /// </summary>
    /// <param name="tsvLines">The lines of an LM22-format TSV (first line = header).</param>
    /// <returns>A signature matrix: cell-type name → (gene symbol → reference expression value).</returns>
    /// <exception cref="ArgumentNullException">If <paramref name="tsvLines"/> is null.</exception>
    /// <exception cref="FormatException">If the header or any data row is malformed (no columns, ragged rows, or non-numeric value).</exception>
    public static IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>> LoadSignatureMatrix(
        IEnumerable<string> tsvLines)
    {
        ArgumentNullException.ThrowIfNull(tsvLines);

        using var e = tsvLines.GetEnumerator();

        // Skip leading blank lines, then read the header.
        string? header = null;
        while (e.MoveNext())
        {
            if (!string.IsNullOrWhiteSpace(e.Current))
            {
                header = e.Current;
                break;
            }
        }

        if (header is null)
            throw new FormatException("Signature matrix is empty: no header row found.");

        string[] headerCols = header.Split('\t');
        if (headerCols.Length < 2)
            throw new FormatException("Signature matrix header must have a gene-symbol column followed by at least one cell-type column.");

        int nCellTypes = headerCols.Length - 1;
        var cellTypeNames = new string[nCellTypes];
        var matrix = new Dictionary<string, Dictionary<string, double>>(nCellTypes);
        for (int j = 0; j < nCellTypes; j++)
        {
            string name = headerCols[j + 1].Trim();
            cellTypeNames[j] = name;
            matrix[name] = new Dictionary<string, double>();
        }

        while (e.MoveNext())
        {
            string line = e.Current;
            if (string.IsNullOrWhiteSpace(line))
                continue;

            string[] cols = line.Split('\t');
            if (cols.Length != headerCols.Length)
                throw new FormatException(
                    $"Ragged signature-matrix row: expected {headerCols.Length} columns, got {cols.Length}.");

            string gene = cols[0].Trim();
            for (int j = 0; j < nCellTypes; j++)
            {
                if (!double.TryParse(cols[j + 1].Trim(), NumberStyles.Float, CultureInfo.InvariantCulture, out double val))
                    throw new FormatException($"Non-numeric signature value '{cols[j + 1]}' for gene '{gene}', cell type '{cellTypeNames[j]}'.");
                matrix[cellTypeNames[j]][gene] = val;
            }
        }

        return matrix.ToDictionary(
            kvp => kvp.Key,
            kvp => (IReadOnlyDictionary<string, double>)kvp.Value);
    }

    /// <summary>
    /// Loads the <b>bundled ABIS-Seq immune-cell signature matrix</b> (Monaco et al., 2019) so that
    /// <see cref="DeconvoluteImmuneCellsNuSvr"/> works out-of-the-box without a caller-supplied matrix.
    /// <para>
    /// The matrix is the RNA-seq (ABIS-Seq) "well-conditioned signature matrix" of Monaco et al. (2019):
    /// <see cref="AbisSignatureGeneCount"/> genes × <see cref="AbisSignatureCellTypeCount"/> immune cell
    /// types (Monocytes C, NK, T CD8 Memory, T CD4 Naive, T CD8 Naive, B Naive, T CD4 Memory, MAIT,
    /// T gd Vd2, Neutrophils LD, T gd non-Vd2, Basophils LD, Monocytes NC+I, B Memory, mDCs, pDCs,
    /// Plasmablasts). Values are mRNA-abundance-normalised expression on the original ABIS scale.
    /// Pass the result as the <c>signatureMatrix</c> argument of <see cref="DeconvoluteImmuneCellsNuSvr"/>
    /// or <see cref="DeconvoluteImmuneCells"/>.
    /// </para>
    /// <para>
    /// <b>Provenance and licence.</b> The matrix is Table S5 (sheet "ABIS-Seq") of the open-access paper,
    /// retrieved from PMC6367568 supplementary file <c>mmc6.xlsx</c>. The article is published under the
    /// Creative Commons Attribution 4.0 (CC BY 4.0) licence ("© 2019 The Authors. This is an open access
    /// article under the CC BY license"), which permits redistribution with attribution; the bundled
    /// resource carries this provenance/licence in its header. Unlike CIBERSORT LM22 (Stanford,
    /// no-redistribution — caller-supplied via <see cref="LoadSignatureMatrix"/>), ABIS may be bundled.
    /// Source: Monaco G, Lee B, Xu W, et al. "RNA-Seq Signatures Normalized by mRNA Abundance Allow
    /// Absolute Deconvolution of Human Immune Cell Types." Cell Reports 26(6):1627-1640.e7, 2019,
    /// doi:10.1016/j.celrep.2019.01.041; PMID:30726743.
    /// </para>
    /// </summary>
    /// <returns>The bundled ABIS-Seq signature matrix: cell-type name → (gene symbol → reference expression value).</returns>
    public static IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>> LoadBundledAbisSignatureMatrix()
    {
        var asm = typeof(ImmuneAnalyzer).Assembly;
        using Stream stream = asm.GetManifestResourceStream(AbisResourceName)
            ?? throw new InvalidOperationException($"Bundled ABIS signature matrix resource '{AbisResourceName}' was not found.");
        using var reader = new StreamReader(stream);

        // The bundled TSV carries a provenance/licence header in '#'-prefixed comment lines; strip
        // them so LoadSignatureMatrix sees the gene/cell-type header as the first non-blank line.
        return LoadSignatureMatrix(ReadNonCommentLines(reader));

        static IEnumerable<string> ReadNonCommentLines(StreamReader r)
        {
            string? line;
            while ((line = r.ReadLine()) is not null)
            {
                if (line.StartsWith('#'))
                    continue;
                yield return line;
            }
        }
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
    /// Standardizes a vector to zero mean and unit (population) standard deviation (z-score).
    /// If the standard deviation is ~0, returns the mean-centred vector (all zeros) to avoid div-by-0.
    /// Source: Newman et al. (2015) — mixture/signature z-score standardisation prior to ν-SVR.
    /// </summary>
    private static double[] Standardize(double[] v)
    {
        int n = v.Length;
        if (n == 0)
            return Array.Empty<double>();

        double mean = v.Average();
        double sumSq = 0.0;
        for (int i = 0; i < n; i++)
        {
            double d = v[i] - mean;
            sumSq += d * d;
        }

        double sd = Math.Sqrt(sumSq / n);
        double[] result = new double[n];
        if (sd < 1e-15)
        {
            for (int i = 0; i < n; i++)
                result[i] = 0.0;
            return result;
        }

        for (int i = 0; i < n; i++)
            result[i] = (v[i] - mean) / sd;

        return result;
    }

    /// <summary>
    /// Standardizes each column of an n×k matrix to zero mean and unit (population) standard deviation.
    /// </summary>
    private static double[,] StandardizeColumns(double[,] matrix, int nRows, int nCols)
    {
        double[,] result = new double[nRows, nCols];
        for (int j = 0; j < nCols; j++)
        {
            double mean = 0.0;
            for (int i = 0; i < nRows; i++)
                mean += matrix[i, j];
            mean /= nRows;

            double sumSq = 0.0;
            for (int i = 0; i < nRows; i++)
            {
                double d = matrix[i, j] - mean;
                sumSq += d * d;
            }

            double sd = Math.Sqrt(sumSq / nRows);
            if (sd < 1e-15)
            {
                for (int i = 0; i < nRows; i++)
                    result[i, j] = 0.0;
            }
            else
            {
                for (int i = 0; i < nRows; i++)
                    result[i, j] = (matrix[i, j] - mean) / sd;
            }
        }

        return result;
    }

    /// <summary>
    /// Solves a linear-kernel ν-support-vector regression of the target vector <paramref name="y"/>
    /// on the columns of the design matrix <paramref name="x"/> (rows = samples/genes,
    /// columns = features/cell types) and returns the primal weight vector <c>w</c> (length nCols).
    /// <para>
    /// Implements the ν-SVR dual of Schölkopf et al. (2000) (Smola &amp; Schölkopf tutorial eqs (60)-(62)):
    /// maximise <c>−½ Σ(α_i−α_i*)(α_j−α_j*)⟨x_i,x_j⟩ + Σ y_i(α_i−α_i*)</c> subject to
    /// <c>Σ(α_i−α_i*)=0</c>, <c>Σ(α_i+α_i*) ≤ C·ν·ℓ</c> and <c>α_i, α_i* ∈ [0, C]</c>, where here a
    /// "sample" is a gene (row of the signature matrix) whose feature vector is that gene's row.
    /// The primal weights are recovered as <c>w = Σ(α_i−α_i*)·rowFeature_i</c> — these are the
    /// cell-type regression coefficients <c>f</c>.
    /// </para>
    /// <para>
    /// The dual QP is convex; it is solved by an SMO-style pairwise coordinate ascent over the dual
    /// variables β_i = α_i − α_i* that exactly maintains the equality constraint Σβ_i = 0 (each step
    /// moves a pair (β_p, β_q) by (+δ, −δ)) and enforces the ν budget Σ|β_i| ≤ C·ν·ℓ and the box
    /// |β_i| ≤ C. At optimum the KKT conditions of the dual hold; see the accompanying tests.
    /// </para>
    /// <para>Source: Schölkopf et al. (2000), Neural Computation 12(5):1207-1245; Newman et al. (2015), Nature Methods 12:453.</para>
    /// </summary>
    /// <param name="x">Design matrix, nRows × nCols (gene rows × cell-type features).</param>
    /// <param name="y">Target vector (mixture), length nRows.</param>
    /// <param name="nRows">Number of samples (genes).</param>
    /// <param name="nCols">Number of features (cell types) = length of the returned weight vector.</param>
    /// <param name="nu">ν parameter in (0, 1].</param>
    /// <param name="cost">Regularisation constant C (box bound on dual variables).</param>
    private static double[] SolveNuSvrLinear(double[,] x, double[] y, int nRows, int nCols, double nu, double cost)
    {
        // Work in the "sample = gene row" space: each sample i has feature vector x[i, :].
        // Dual variable per sample: beta_i = alpha_i - alpha_i* in [-C, C].
        // Linear kernel K_ij = <x_i, x_j>.
        int l = nRows;
        double[,] k = new double[l, l];
        for (int i = 0; i < l; i++)
        {
            for (int j = i; j < l; j++)
            {
                double dot = 0.0;
                for (int c = 0; c < nCols; c++)
                    dot += x[i, c] * x[j, c];
                k[i, j] = dot;
                k[j, i] = dot;
            }
        }

        // ν budget on Σ(α_i + α_i*) = Σ|β_i|.
        double nuBudget = cost * nu * l;

        double[] beta = new double[l];           // β_i = α_i − α_i*
        // Gradient of the dual objective W(β) = −½ βᵀKβ + Σ y_i β_i (subject to constraints)
        // with respect to β_i: g_i = y_i − Σ_j K_ij β_j. (β starts at 0 ⇒ g_i = y_i.)
        double[] grad = new double[l];
        Array.Copy(y, grad, l);

        const double Tol = 1e-10;
        const double Eps = 1e-12;
        int maxIterations = 200 * Math.Max(l, 1);

        for (int iter = 0; iter < maxIterations; iter++)
        {
            // Select a working pair (p, q) by maximum-violating-pair on the equality constraint.
            // Moving (β_p, β_q) by (+δ, −δ) keeps Σβ_i = 0. The directional gradient is g_p − g_q;
            // pick p maximising g (subject to being allowed to increase) and q minimising g.
            double sumAbs = 0.0;
            for (int i = 0; i < l; i++)
                sumAbs += Math.Abs(beta[i]);

            int p = -1, q = -1;
            double gMaxUp = double.NegativeInfinity;   // best candidate to increase β
            double gMaxDown = double.PositiveInfinity;  // best candidate to decrease β

            for (int i = 0; i < l; i++)
            {
                bool canIncrease = beta[i] < cost - Eps && (sumAbs < nuBudget - Eps || beta[i] < -Eps);
                bool canDecrease = beta[i] > -cost + Eps && (sumAbs < nuBudget - Eps || beta[i] > Eps);

                if (canIncrease && grad[i] > gMaxUp)
                {
                    gMaxUp = grad[i];
                    p = i;
                }
                if (canDecrease && grad[i] < gMaxDown)
                {
                    gMaxDown = grad[i];
                    q = i;
                }
            }

            if (p == -1 || q == -1 || p == q)
                break;

            double violation = gMaxUp - gMaxDown;
            if (violation < Tol)
                break;

            // Unconstrained optimal step for moving (β_p += δ, β_q -= δ):
            // δ* = (g_p − g_q) / (K_pp + K_qq − 2 K_pq).
            double curvature = k[p, p] + k[q, q] - 2.0 * k[p, q];
            if (curvature < Eps)
                curvature = Eps;

            double delta = (grad[p] - grad[q]) / curvature;
            if (delta <= 0.0)
                break;

            // Clip δ against the box on β_p and β_q.
            double maxUpP = cost - beta[p];      // β_p may rise to +C
            double maxDownQ = beta[q] + cost;    // β_q may fall to −C
            delta = Math.Min(delta, Math.Min(maxUpP, maxDownQ));

            // Clip δ against the ν budget Σ|β_i| ≤ nuBudget.
            // The change in Σ|β| from this (+δ, −δ) move depends on the signs of β_p, β_q.
            delta = ClipDeltaForNuBudget(beta[p], beta[q], delta, sumAbs, nuBudget, Eps);

            if (delta <= Eps)
                break;

            beta[p] += delta;
            beta[q] -= delta;

            // Update gradient: g_i -= δ (K_ip − K_iq).
            for (int i = 0; i < l; i++)
                grad[i] -= delta * (k[i, p] - k[i, q]);
        }

        // Recover primal weights w = Σ_i β_i x_i (length nCols).
        double[] w = new double[nCols];
        for (int i = 0; i < l; i++)
        {
            double bi = beta[i];
            if (bi == 0.0)
                continue;
            for (int c = 0; c < nCols; c++)
                w[c] += bi * x[i, c];
        }

        return w;
    }

    /// <summary>
    /// Clips the SMO step δ for a (+δ, −δ) move on (β_p, β_q) so that Σ|β_i| does not exceed the
    /// ν budget <paramref name="nuBudget"/>. Computes the largest feasible δ ∈ (0, <paramref name="delta"/>]
    /// such that |β_p+δ| + |β_q−δ| − |β_p| − |β_q| keeps the total ≤ budget.
    /// </summary>
    private static double ClipDeltaForNuBudget(double betaP, double betaQ, double delta, double sumAbs, double nuBudget, double eps)
    {
        // ΔsumAbs(δ) = (|β_p+δ| − |β_p|) + (|β_q−δ| − |β_q|). Piecewise-linear, non-decreasing slope
        // in regions; the worst case is when both moves add magnitude (β_p ≥ 0, β_q ≤ 0), giving
        // ΔsumAbs = 2δ. Solve sumAbs + ΔsumAbs(δ) ≤ nuBudget for δ.
        double headroom = nuBudget - sumAbs;
        if (headroom <= eps)
        {
            // No budget headroom: only moves that do not increase Σ|β| are allowed, i.e. δ limited so
            // that β_p moves toward 0 (β_p < 0) and/or β_q moves toward 0 (β_q > 0).
            // Slope of ΔsumAbs at δ→0+ is sign-based: +1 if β_p ≥ 0 else −1, plus +1 if β_q ≤ 0 else −1.
            double slope = (betaP >= 0 ? 1.0 : -1.0) + (betaQ <= 0 ? 1.0 : -1.0);
            if (slope > 0)
                return 0.0;
            // slope ≤ 0: the move (weakly) frees budget up to the kink where a sign flips.
            double kink = double.PositiveInfinity;
            if (betaP < 0) kink = Math.Min(kink, -betaP);       // β_p reaches 0
            if (betaQ > 0) kink = Math.Min(kink, betaQ);        // β_q reaches 0
            return Math.Min(delta, kink);
        }

        // Compute ΔsumAbs as a function of δ and cap δ so it does not exceed headroom.
        // Use the conservative worst-case slope of +2 only where both add magnitude; otherwise the
        // exact piecewise value. Walk to the first sign kink, accumulate.
        double remaining = delta;
        double budgetLeft = headroom;
        double bp = betaP, bq = betaQ;
        double allowed = 0.0;

        while (remaining > eps)
        {
            double slope = (bp >= 0 ? 1.0 : -1.0) + (bq <= 0 ? 1.0 : -1.0);

            // Distance to the next kink (where bp hits 0 going up, or bq hits 0 going down).
            double kink = remaining;
            if (bp < 0) kink = Math.Min(kink, -bp);
            if (bq > 0) kink = Math.Min(kink, bq);

            if (slope <= 0)
            {
                // This segment does not consume budget; take it fully.
                allowed += kink;
                bp += kink;
                bq -= kink;
                remaining -= kink;
                continue;
            }

            // slope > 0: this segment consumes budget at the given slope.
            double maxBySegment = budgetLeft / slope;
            double step = Math.Min(kink, maxBySegment);
            allowed += step;
            budgetLeft -= step * slope;
            bp += step;
            bq -= step;
            remaining -= step;

            if (budgetLeft <= eps)
                break;
        }

        return allowed;
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
