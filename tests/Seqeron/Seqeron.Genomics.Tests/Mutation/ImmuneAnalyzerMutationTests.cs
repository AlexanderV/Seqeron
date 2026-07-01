using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.ImmuneAnalyzer;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// ONCO-IMMUNE-001 mutation killers: deterministic NNLS deconvolution cases with known
/// exact solutions that exercise the full Lawson-Hanson active-set solver (including the
/// infeasibility inner loop that the canonical tests never reach) and pin the Pearson /
/// RMSE fit metrics.
///
/// Evidence: NNLS deconvolution m = S·f, f ≥ 0, Σf = 1 (Abbas et al. 2009; Lawson &amp;
/// Hanson 1995, Ch. 23); ESTIMATE purity cos(a + b·score) (Yoshihara et al. 2013).
/// </summary>
[TestFixture]
public class ImmuneAnalyzerMutationTests
{
    private const double Tol = 1e-6;

    private static IReadOnlyDictionary<string, IReadOnlyDictionary<string, double>> Sig(
        params (string Cell, (string Gene, double Val)[] Genes)[] cells)
    {
        var outer = new Dictionary<string, IReadOnlyDictionary<string, double>>();
        foreach (var (cell, genes) in cells)
        {
            var inner = new Dictionary<string, double>();
            foreach (var (g, v) in genes) inner[g] = v;
            outer[cell] = inner;
        }
        return outer;
    }

    #region NNLS deconvolution — exact recoverable mixtures

    [Test]
    public void DeconvoluteImmuneCells_ExactTwoTypeMixture_RecoversFractions()
    {
        // Disjoint marker genes ⇒ the mixture 0.3·TypeA + 0.7·TypeB is exactly recoverable.
        var sig = Sig(
            ("TypeA", new[] { ("GA1", 2.0), ("GA2", 4.0) }),
            ("TypeB", new[] { ("GB1", 3.0), ("GB2", 6.0) }));
        var expr = new Dictionary<string, double>
        {
            ["GA1"] = 0.6, ["GA2"] = 1.2, ["GB1"] = 2.1, ["GB2"] = 4.2,
        };

        var r = DeconvoluteImmuneCells(expr, sig);

        Assert.That(r.CellFractions["TypeA"], Is.EqualTo(0.3).Within(Tol));
        Assert.That(r.CellFractions["TypeB"], Is.EqualTo(0.7).Within(Tol));
        Assert.That(r.OverlappingGenes, Is.EqualTo(4));
        Assert.That(r.Correlation, Is.EqualTo(1.0).Within(Tol)); // perfect reconstruction
        Assert.That(r.Rmse, Is.EqualTo(0.0).Within(Tol));
    }

    [Test]
    public void DeconvoluteImmuneCells_InfeasibleJointFit_ActiveSetZeroesOneType()
    {
        // S = [[1,1],[1,2],[1,3]], m = [1,1,1]. The unconstrained joint fit drives TypeB to 0,
        // forcing the Lawson-Hanson infeasibility inner loop; the NNLS solution is f = (1, 0).
        var sig = Sig(
            ("TypeA", new[] { ("G1", 1.0), ("G2", 1.0), ("G3", 1.0) }),
            ("TypeB", new[] { ("G1", 1.0), ("G2", 2.0), ("G3", 3.0) }));
        var expr = new Dictionary<string, double> { ["G1"] = 1.0, ["G2"] = 1.0, ["G3"] = 1.0 };

        var r = DeconvoluteImmuneCells(expr, sig);

        Assert.That(r.CellFractions["TypeA"], Is.EqualTo(1.0).Within(Tol));
        Assert.That(r.CellFractions["TypeB"], Is.EqualTo(0.0).Within(Tol));
    }

    [Test]
    public void DeconvoluteImmuneCells_SingleTypeImperfectFit_ExactRmseAndCorrelation()
    {
        // One cell type {G1:2, G2:4}; observed [1,3] is not proportional ⇒ NNLS gives the
        // least-squares scalar, normalised to fraction 1.0, reconstructing [2,4].
        // RMSE = sqrt(((1-2)² + (3-4)²)/2) = 1.0; Pearson(observed, reconstructed) = 1.0.
        var sig = Sig(("TypeA", new[] { ("G1", 2.0), ("G2", 4.0) }));
        var expr = new Dictionary<string, double> { ["G1"] = 1.0, ["G2"] = 3.0 };

        var r = DeconvoluteImmuneCells(expr, sig);

        Assert.That(r.CellFractions["TypeA"], Is.EqualTo(1.0).Within(Tol));
        Assert.That(r.Rmse, Is.EqualTo(1.0).Within(Tol));         // kills sumSqErr*n vs /n
        Assert.That(r.Correlation, Is.EqualTo(1.0).Within(Tol));
    }

    #endregion

    #region ESTIMATE tumor purity

    [Test]
    public void EstimateInfiltration_EmptyProfile_PurityIsCosOfInterceptA()
    {
        // Empty profile ⇒ scores 0; purity = cos(a + b·0) = cos(a), overlaps 0.
        var r = EstimateInfiltration(new Dictionary<string, double>());

        Assert.That(r.ImmuneScore, Is.EqualTo(0.0).Within(Tol));
        Assert.That(r.StromalScore, Is.EqualTo(0.0).Within(Tol));
        Assert.That(r.OverlappingImmuneGenes, Is.EqualTo(0));
        Assert.That(r.TumorPurity, Is.EqualTo(Math.Cos(EstimatePurityCoefficientA)).Within(Tol));
    }

    #endregion

    #region LoadSignatureMatrix — header/row parsing (CIBERSORT-style TSV)

    // A 2-column header (gene-symbol column + exactly ONE cell-type column) is the SMALLEST valid
    // signature matrix and must parse, pinning the "headerCols.Length < 2" guard against the "<= 2"
    // off-by-one (which would wrongly reject a single-cell-type matrix). Values are pinned exactly.
    [Test]
    public void LoadSignatureMatrix_TwoColumnHeader_IsSmallestValidMatrix()
    {
        var lines = new[] { "Gene\tTcell", "CD8A\t5.5", "CD4\t2.0" };

        var matrix = ImmuneAnalyzer.LoadSignatureMatrix(lines);

        Assert.That(matrix.Keys, Is.EquivalentTo(new[] { "Tcell" }), "one cell-type column");
        Assert.That(matrix["Tcell"]["CD8A"], Is.EqualTo(5.5).Within(Tol));
        Assert.That(matrix["Tcell"]["CD4"], Is.EqualTo(2.0).Within(Tol));
    }

    // Exact per-cell-type value placement across MULTIPLE columns pins the column index arithmetic
    // (headerCols[j+1] / cols[j+1]): an off-by-one would map values to the wrong cell type.
    [Test]
    public void LoadSignatureMatrix_MultipleColumns_PlacesEachValueInItsCellType()
    {
        var lines = new[] { "Gene\tBcell\tNK\tMono", "MS4A1\t9.0\t1.0\t0.5", "NCAM1\t0.2\t8.0\t0.3" };

        var matrix = ImmuneAnalyzer.LoadSignatureMatrix(lines);

        Assert.Multiple(() =>
        {
            Assert.That(matrix["Bcell"]["MS4A1"], Is.EqualTo(9.0).Within(Tol));
            Assert.That(matrix["NK"]["MS4A1"], Is.EqualTo(1.0).Within(Tol));
            Assert.That(matrix["Mono"]["MS4A1"], Is.EqualTo(0.5).Within(Tol));
            Assert.That(matrix["NK"]["NCAM1"], Is.EqualTo(8.0).Within(Tol), "NCAM1 is the NK marker here");
            Assert.That(matrix["Bcell"]["NCAM1"], Is.EqualTo(0.2).Within(Tol));
        });
    }

    // A non-numeric signature value must raise a FormatException (NOT an index error): this drives the
    // error path whose message references cols[j+1] / cellTypeNames[j], pinning that column index.
    [Test]
    public void LoadSignatureMatrix_NonNumericValueInFirstCellType_ThrowsFormatException()
    {
        var lines = new[] { "Gene\tTcell\tBcell", "CD8A\tNaNish\t2.0" };

        Assert.Throws<FormatException>(() => ImmuneAnalyzer.LoadSignatureMatrix(lines),
            "a non-numeric value in the first cell-type column is a format error, not an index error");
    }

    // When NONE of the signature genes are present in the expression profile, the deconvolution must
    // early-out to all-zero fractions. The guard is "overlappingGenes == 0 OR cellTypes == 0"; the
    // OR-vs-AND mutant would (with cell types present) skip the early-out and try to deconvolve with an
    // empty gene set. Pins the disjunction.
    [Test]
    public void DeconvoluteImmuneCells_NoOverlappingGenes_EarlyOutsToZeroFractions()
    {
        var sig = Sig(("Tcell", new[] { ("CD8A", 5.0), ("CD3D", 4.0) }));
        var expr = new Dictionary<string, double> { ["FOXP3"] = 7.0, ["MKI67"] = 3.0 }; // disjoint genes

        var r = DeconvoluteImmuneCells(expr, sig);

        Assert.That(r.OverlappingGenes, Is.EqualTo(0), "no signature gene is present in the profile");
        Assert.That(r.CellFractions["Tcell"], Is.EqualTo(0.0).Within(Tol),
            "zero overlapping genes ⇒ all-zero fractions via the early-out guard");
    }

    #endregion
}
