using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for <see cref="FitSubclonalCopyNumber"/> — the port of Battenberg <c>determine_copynumber</c>
/// (Wedge-lab/battenberg R/fitcopynumber.R; Nik-Zainal et al. 2012) added by the 2026-09 review (B24 F14), which had no
/// heavy-tier coverage.
///
/// Documented contract (XML doc + Battenberg code):
/// <list type="bullet">
/// <item>Valid input (finite logR, BAF ∈ [0, 1], ρ ∈ (0, 1], ψ &gt; 0, γ &gt; 0) never throws; one fit per segment, in order.</item>
/// <item>Every state has integer major/minor copy numbers ≥ 0 (a negative minor is raised to 0.01 before flooring;
/// Battenberg <c>cn_upper_limit</c> = 1000 at BAF 1).</item>
/// <item>Clonal ⇔ no second state, fraction exactly 1; sub-clonal ⇔ a second state with fraction 1 − τ.</item>
/// <item>Invalid input ⇒ ArgumentNullException / ArgumentException (segment signal) / ArgumentOutOfRangeException (ρ, ψ, γ).</item>
/// </list>
/// All randomness is locally seeded.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologySubclonalCopyNumberFuzzTests
{
    private static AlleleSpecificSegmentSummary Seg(double logR, double baf) => new("1", 0, 1000, logR, baf, 25);

    [Test]
    [CancelAfter(20_000)]
    public void FitSubclonalCopyNumber_ExtremeValidInputs_WellFormedStates()
    {
        double[] logRs = { -20.0, -5.0, -1.0, -1e-9, 0.0, 1e-9, 0.585, 1.0, 3.0, 8.0 };
        double[] bafs = { 0.0, 1e-12, 0.01, 0.25, 0.5, 0.5 + 1e-15, 0.75, 0.99, 1.0 - 1e-12, 1.0 };
        double[] purities = { 1e-6, 0.01, 0.2, 0.5, 0.9, 1.0 };
        double[] ploidies = { 1e-3, 1.0, 2.0, 3.7, 8.0 };
        double[] gammas = { 0.55, 1.0 };

        foreach (double rho in purities)
        foreach (double psi in ploidies)
        foreach (double gamma in gammas)
        {
            var segs = logRs.SelectMany(r => bafs.Select(b => Seg(r, b))).ToArray();
            var fits = FitSubclonalCopyNumber(segs, rho, psi, gamma);

            fits.Should().HaveCount(segs.Length);
            for (int i = 0; i < fits.Count; i++)
            {
                var f = fits[i];
                f.Segment.Should().Be(segs[i], "fits are in input order");
                f.PrimaryState.MajorCopyNumber.Should().BeGreaterThanOrEqualTo(0);
                f.PrimaryState.MinorCopyNumber.Should().BeGreaterThanOrEqualTo(0);
                if (f.IsSubclonal)
                {
                    f.SecondaryState.Should().NotBeNull("a sub-clonal segment has two states (ρ={0}, ψ={1}, i={2})", rho, psi, i);
                    f.SecondaryState!.Value.MajorCopyNumber.Should().BeGreaterThanOrEqualTo(0);
                    f.SecondaryState!.Value.MinorCopyNumber.Should().BeGreaterThanOrEqualTo(0);
                }
                else
                {
                    f.SecondaryState.Should().BeNull();
                    f.PrimaryState.CellFraction.Should().Be(1.0, "a clonal state covers every tumour cell");
                }
            }
        }
    }

    [Test]
    [CancelAfter(20_000)]
    public void FitSubclonalCopyNumber_RandomValidInputs_NeverThrowAndStatesNonNegative()
    {
        for (int seed = 0; seed < 400; seed++)
        {
            var rng = new Random(seed);
            double rho = rng.Next(4) == 0 ? 1.0 : Math.Max(1e-9, rng.NextDouble());
            double psi = 0.1 + 7.9 * rng.NextDouble();
            var segs = Enumerable.Range(0, 1 + rng.Next(20))
                .Select(_ => Seg((rng.NextDouble() - 0.5) * 12.0, rng.Next(6) == 0 ? rng.Next(2) : rng.NextDouble()))
                .ToArray();

            var fits = FitSubclonalCopyNumber(segs, rho, psi);

            fits.Should().HaveCount(segs.Length, "seed {0}", seed);
            fits.All(f => f.PrimaryState.MajorCopyNumber >= 0 && f.PrimaryState.MinorCopyNumber >= 0
                          && (f.IsSubclonal == f.SecondaryState.HasValue)).Should().BeTrue("seed {0}", seed);
        }
    }

    [Test]
    public void FitSubclonalCopyNumber_EmptyInput_ReturnsEmpty() =>
        FitSubclonalCopyNumber(Array.Empty<AlleleSpecificSegmentSummary>(), 0.7, 2.0).Should().BeEmpty();

    [Test]
    public void FitSubclonalCopyNumber_InvalidArguments_ThrowDocumentedExceptions()
    {
        var ok = new[] { Seg(0.0, 0.5) };
        ((Action)(() => FitSubclonalCopyNumber(null!, 0.7, 2.0))).Should().Throw<ArgumentNullException>();

        foreach (double rho in new[] { 0.0, -0.1, 1.0 + 1e-12, double.NaN, double.PositiveInfinity })
            ((Action)(() => FitSubclonalCopyNumber(ok, rho, 2.0))).Should().Throw<ArgumentOutOfRangeException>("ρ = {0}", rho);

        foreach (double psi in new[] { 0.0, -2.0, double.NaN, double.PositiveInfinity })
            ((Action)(() => FitSubclonalCopyNumber(ok, 0.7, psi))).Should().Throw<ArgumentOutOfRangeException>("ψ = {0}", psi);

        foreach (double gamma in new[] { 0.0, -1.0, double.NaN, double.PositiveInfinity })
            ((Action)(() => FitSubclonalCopyNumber(ok, 0.7, 2.0, gamma))).Should().Throw<ArgumentOutOfRangeException>("γ = {0}", gamma);

        foreach (var bad in new[] { Seg(double.NaN, 0.5), Seg(double.PositiveInfinity, 0.5), Seg(0.0, -0.01), Seg(0.0, 1.01), Seg(0.0, double.NaN) })
            ((Action)(() => FitSubclonalCopyNumber(new[] { bad }, 0.7, 2.0))).Should().Throw<ArgumentException>("{0}", bad);
    }
}
