// StatisticsHelper.Sigmoid / SolveSymmetricPositiveDefinite — shared helpers extracted from MiRnaAnalyzer (review 2026-09, B13 WP-R2-6,
// DUP_MAP §22). Reference values: scipy.special.expit and numpy.linalg.solve (numpy 2.x / scipy 1.x).

using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class StatisticsHelper_SigmoidCholesky_Tests
{
    [TestCase(0.0, 0.5)]
    [TestCase(2.5, 0.9241418199787566)]
    [TestCase(-40.0, 4.248354255291589e-18)]
    public void Sigmoid_MatchesScipyExpit(double z, double expected) =>
        Assert.That(StatisticsHelper.Sigmoid(z), Is.EqualTo(expected).Within(1e-15 * Math.Max(1, expected)).Or.EqualTo(expected).Within(1e-30));

    [Test]
    public void SolveSymmetricPositiveDefinite_MatchesNumpySolve()
    {
        // Classic Cholesky example A = L·Lᵀ with L = [[2,0,0],[6,1,0],[-8,5,3]].
        var a = new double[,] { { 4, 12, -16 }, { 12, 37, -43 }, { -16, -43, 98 } };
        var x = StatisticsHelper.SolveSymmetricPositiveDefinite(a, new[] { 1.0, 2, 3 });
        Assert.That(x, Is.EqualTo(new[] { 28.583333333333332, -7.666666666666667, 1.3333333333333333 }).Within(1e-12));
        var y = StatisticsHelper.SolveSymmetricPositiveDefinite(new double[,] { { 2, -1 }, { -1, 2 } }, new[] { 1.0, 0 });
        Assert.That(y, Is.EqualTo(new[] { 2.0 / 3, 1.0 / 3 }).Within(1e-15));
    }

    [Test]
    public void SolveSymmetricPositiveDefinite_ContractAndNoMutation()
    {
        var a = new double[,] { { 2, -1 }, { -1, 2 } };
        var b = new[] { 1.0, 0 };
        StatisticsHelper.SolveSymmetricPositiveDefinite(a, b);
        Assert.That(a, Is.EqualTo(new double[,] { { 2, -1 }, { -1, 2 } }));
        Assert.That(b, Is.EqualTo(new[] { 1.0, 0 }));
        Assert.Throws<ArgumentNullException>(() => StatisticsHelper.SolveSymmetricPositiveDefinite(null!, b));
        Assert.Throws<ArgumentNullException>(() => StatisticsHelper.SolveSymmetricPositiveDefinite(a, null!));
        Assert.Throws<ArgumentException>(() => StatisticsHelper.SolveSymmetricPositiveDefinite(a, new[] { 1.0, 2, 3 }));
        // Not positive definite (negative pivot) => NaN, no exception.
        var bad = StatisticsHelper.SolveSymmetricPositiveDefinite(new double[,] { { -1, 0 }, { 0, 1 } }, new[] { 1.0, 1 });
        Assert.That(double.IsNaN(bad[0]), Is.True);
    }

    [TestCase("", "")]
    [TestCase("A", "A")]
    [TestCase("AUGC", "CGUA")]
    [TestCase("acgu N-", "-N ugca")]
    public void ReverseSequence_ReversesCharacters(string s, string expected) =>
        Assert.That(SequenceExtensions.ReverseSequence(s), Is.EqualTo(expected));

    [Test]
    public void ReverseSequence_NullThrows_InputUnchanged()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceExtensions.ReverseSequence(null!));
        string s = "AUGC";
        _ = SequenceExtensions.ReverseSequence(s);
        Assert.That(s, Is.EqualTo("AUGC"));
    }
}
