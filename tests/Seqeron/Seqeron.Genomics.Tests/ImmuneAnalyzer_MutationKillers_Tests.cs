using System;
using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.ImmuneAnalyzer;

namespace Seqeron.Genomics.Tests;

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
public class ImmuneAnalyzer_MutationKillers_Tests
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
}
