// IMMUNE-NUSVR-001 / ONCO-IMMUNE-001 (08_DIFFERENTIAL strategy REF) — differential killers for the
// pure numeric kernels of the CIBERSORT/ABIS deconvolution pipeline. These kernels (Gaussian
// elimination, Pearson correlation, population z-score standardization, RMSE) are private and only
// observable through the robust NNLS/ν-SVR solvers, which RE-CONVERGE to the same fractions under an
// internal arithmetic mutation — so their index/operator mutants survive at the public API. They are
// exposed internal (Seqeron.Genomics.Tests via IVT) and asserted against CLOSED-FORM expected values
// derived by hand (independent of the implementation), so a mutated index/operator diverges -> killed.

using System;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

[TestFixture]
[Category("IMMUNE-NUSVR-001")]
public class ImmuneAnalyzer_NumericKernels_Tests
{
    private const double Tol = 1e-10;

    // ---- SolveLinearSystem: Gaussian elimination with partial pivoting ----

    [Test]
    public void SolveLinearSystem_2x2_ExactSolution()
    {
        // 2x + y = 3 ; x + 3y = 5  =>  x = 0.8, y = 1.4 (hand-solved).
        var a = new double[,] { { 2, 1 }, { 1, 3 } };
        var x = ImmuneAnalyzer.SolveLinearSystem(a, new double[] { 3, 5 }, 2);
        Assert.That(x[0], Is.EqualTo(0.8).Within(Tol));
        Assert.That(x[1], Is.EqualTo(1.4).Within(Tol));
    }

    [Test]
    public void SolveLinearSystem_ZeroLeadingPivot_RequiresRowSwap()
    {
        // [[0,1],[1,0]] x = [2,3]. The first pivot a[0,0]=0 forces a partial-pivot row swap;
        // solution is x0 = 3, x1 = 2. Kills the pivot-search / row-swap mutants.
        var a = new double[,] { { 0, 1 }, { 1, 0 } };
        var x = ImmuneAnalyzer.SolveLinearSystem(a, new double[] { 2, 3 }, 2);
        Assert.That(x[0], Is.EqualTo(3.0).Within(Tol));
        Assert.That(x[1], Is.EqualTo(2.0).Within(Tol));
    }

    [Test]
    public void SolveLinearSystem_3x3_ExactSolution()
    {
        // Classic system with solution (2, 3, -1).
        var a = new double[,] { { 2, 1, -1 }, { -3, -1, 2 }, { -2, 1, 2 } };
        var b = new double[] { 8, -11, -3 };
        var x = ImmuneAnalyzer.SolveLinearSystem(a, b, 3);
        Assert.That(x[0], Is.EqualTo(2.0).Within(1e-9));
        Assert.That(x[1], Is.EqualTo(3.0).Within(1e-9));
        Assert.That(x[2], Is.EqualTo(-1.0).Within(1e-9));

        // Reconstruction invariant A·x == b (kills back-substitution / elimination index mutants).
        for (int i = 0; i < 3; i++)
        {
            double row = a[i, 0] * x[0] + a[i, 1] * x[1] + a[i, 2] * x[2];
            Assert.That(row, Is.EqualTo(b[i]).Within(1e-9), $"row {i} reconstruction");
        }
    }

    // ---- ComputePearsonCorrelation ----

    [Test]
    public void Pearson_PerfectPositive_IsOne()
    {
        var r = ImmuneAnalyzer.ComputePearsonCorrelation(new double[] { 1, 2, 3 }, new double[] { 2, 4, 6 });
        Assert.That(r, Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void Pearson_PerfectNegative_IsMinusOne()
    {
        var r = ImmuneAnalyzer.ComputePearsonCorrelation(new double[] { 1, 2, 3 }, new double[] { 6, 4, 2 });
        Assert.That(r, Is.EqualTo(-1.0).Within(Tol));
    }

    [Test]
    public void Pearson_KnownIntermediate_IsExact()
    {
        // x=[1,2,3,4,5], y=[2,1,4,3,5]: n=5, num=5*53-15*15=40, varX=varY=50, r=40/50=0.8.
        var r = ImmuneAnalyzer.ComputePearsonCorrelation(
            new double[] { 1, 2, 3, 4, 5 }, new double[] { 2, 1, 4, 3, 5 });
        Assert.That(r, Is.EqualTo(0.8).Within(Tol));
    }

    [Test]
    public void Pearson_ZeroVarianceOrTooShort_IsZero()
    {
        Assert.That(ImmuneAnalyzer.ComputePearsonCorrelation(new double[] { 5, 5, 5 }, new double[] { 1, 2, 3 }),
            Is.EqualTo(0.0), "constant x → undefined → 0");
        Assert.That(ImmuneAnalyzer.ComputePearsonCorrelation(new double[] { 7 }, new double[] { 9 }),
            Is.EqualTo(0.0), "n < 2 → 0");
    }

    // ---- Standardize / StandardizeColumns: population z-score (÷ n) ----

    [Test]
    public void Standardize_KnownVector_IsExactPopulationZScore()
    {
        // [1,2,3]: mean=2, population sd = sqrt(2/3); z = (-1,0,1)/sd = (-sqrt(1.5), 0, sqrt(1.5)).
        var z = ImmuneAnalyzer.Standardize(new double[] { 1, 2, 3 });
        double s = Math.Sqrt(1.5);
        Assert.That(z[0], Is.EqualTo(-s).Within(Tol));
        Assert.That(z[1], Is.EqualTo(0.0).Within(Tol));
        Assert.That(z[2], Is.EqualTo(s).Within(Tol));
    }

    [Test]
    public void Standardize_ConstantVector_IsAllZero()
    {
        var z = ImmuneAnalyzer.Standardize(new double[] { 4, 4, 4 });
        Assert.That(z, Is.EqualTo(new double[] { 0, 0, 0 }));
    }

    [Test]
    public void StandardizeColumns_PerColumnPopulationZScore()
    {
        // col0=[1,2,3], col1=[10,20,30] → both standardize to (-sqrt(1.5), 0, sqrt(1.5)).
        var m = new double[,] { { 1, 10 }, { 2, 20 }, { 3, 30 } };
        var r = ImmuneAnalyzer.StandardizeColumns(m, 3, 2);
        double s = Math.Sqrt(1.5);
        Assert.That(r[0, 0], Is.EqualTo(-s).Within(Tol));
        Assert.That(r[1, 0], Is.EqualTo(0.0).Within(Tol));
        Assert.That(r[2, 0], Is.EqualTo(s).Within(Tol));
        Assert.That(r[0, 1], Is.EqualTo(-s).Within(Tol));
        Assert.That(r[2, 1], Is.EqualTo(s).Within(Tol));
    }

    // ---- ComputeRmse ----

    [Test]
    public void Rmse_KnownVectors_IsExact()
    {
        Assert.That(ImmuneAnalyzer.ComputeRmse(new double[] { 1, 2, 3 }, new double[] { 1, 2, 3 }),
            Is.EqualTo(0.0).Within(Tol), "identical → 0");
        // observed=[0,0], predicted=[3,4]: sqrt((9+16)/2) = sqrt(12.5) = 3.53553390593…
        Assert.That(ImmuneAnalyzer.ComputeRmse(new double[] { 0, 0 }, new double[] { 3, 4 }),
            Is.EqualTo(Math.Sqrt(12.5)).Within(Tol));
    }
}
