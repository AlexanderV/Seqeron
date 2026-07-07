// ONCO-CLONAL-001 — Clonal vs Subclonal Mutation Classification
// Evidence: docs/Evidence/ONCO-CLONAL-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-CLONAL-001.md
// Source: Landau DA et al. (2013). Evolution and Impact of Subclonal Mutations in CLL.
//         Cell 152(4):714-726. https://doi.org/10.1016/j.cell.2013.01.019
//         Satas G et al. (2021). DeCiFering the CCF. Cell Systems 12(10):1004-1018.
//         https://doi.org/10.1016/j.cels.2021.07.006
//
// Expected CCF means and P(CCF>0.95) below are derived independently from the Landau (2013)
// model — f(c)=alpha*M*c/(2(1-alpha)+alpha*q); P(c) prop Binom(a|N,f(c)); uniform prior;
// 100-point grid c in [0.01,1]; clonal iff P(CCF>0.95) > 0.5 — NOT copied from the implementation.

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_ClassifyClonality_Tests
{
    private static OncologyAnalyzer.ClonalityVariant V(int alt, int total, int cn, int mult) =>
        new(alt, total, cn, mult);

    private static OncologyAnalyzer.ClonalityVariant V(int alt, int total, int cn) =>
        new(alt, total, cn);

    #region ClassifyClonality

    // M1 — Clonal, pure het deep coverage: a=300,N=300,q=2,M=1,rho=1.0 -> f(c)=c, a=N. Landau (2013).
    [Test]
    public void ClassifyClonality_PureHetAllAltReads_IsClonal()
    {
        OncologyAnalyzer.ClonalityResult result =
            OncologyAnalyzer.ClassifyClonality(new[] { V(300, 300, 2) }, purity: 1.0);

        OncologyAnalyzer.ClonalityCall call = result.Calls[0];
        Assert.Multiple(() =>
        {
            Assert.That(call.Status, Is.EqualTo(OncologyAnalyzer.ClonalityStatus.Clonal),
                "P(CCF>0.95)=1.0 > 0.5 ⇒ clonal (Landau 2013)");
            Assert.That(call.Ccf, Is.EqualTo(0.999486).Within(1e-6),
                "posterior-mean CCF for a=N at f(c)=c concentrates at the grid top");
            Assert.That(call.ProbabilityClonal, Is.EqualTo(1.0).Within(1e-6),
                "all posterior mass lies above CCF 0.95");
        });
    }

    // M2 — Clonal, impure het deep coverage: a=400,N=1000,q=2,M=1,rho=0.8 (f at c=1 is 0.4). Landau (2013).
    [Test]
    public void ClassifyClonality_ImpureHetClonal_IsClonal()
    {
        OncologyAnalyzer.ClonalityResult result =
            OncologyAnalyzer.ClassifyClonality(new[] { V(400, 1000, 2) }, purity: 0.8);

        OncologyAnalyzer.ClonalityCall call = result.Calls[0];
        Assert.Multiple(() =>
        {
            Assert.That(call.Status, Is.EqualTo(OncologyAnalyzer.ClonalityStatus.Clonal),
                "P(CCF>0.95)≈0.864 > 0.5 ⇒ clonal (Landau 2013)");
            Assert.That(call.Ccf, Is.EqualTo(0.972455).Within(1e-6), "posterior-mean CCF");
            Assert.That(call.ProbabilityClonal, Is.EqualTo(0.864167).Within(1e-6),
                "posterior mass above CCF 0.95");
        });
    }

    // M3 — Subclonal, CCF~0.6: a=240,N=1000,q=2,M=1,rho=0.8 (f=0.6*0.4=0.24). Landau (2013).
    [Test]
    public void ClassifyClonality_Ccf06_IsSubclonal()
    {
        OncologyAnalyzer.ClonalityResult result =
            OncologyAnalyzer.ClassifyClonality(new[] { V(240, 1000, 2) }, purity: 0.8);

        OncologyAnalyzer.ClonalityCall call = result.Calls[0];
        Assert.Multiple(() =>
        {
            Assert.That(call.Status, Is.EqualTo(OncologyAnalyzer.ClonalityStatus.Subclonal),
                "P(CCF>0.95)≈0 < 0.5 ⇒ subclonal (Landau 2013)");
            Assert.That(call.Ccf, Is.EqualTo(0.601297).Within(1e-6), "posterior-mean CCF ≈ 0.6");
            Assert.That(call.ProbabilityClonal, Is.EqualTo(0.0).Within(1e-5),
                "essentially no mass above CCF 0.95");
        });
    }

    // M4 — Subclonal, CCF~0.4: a=200,N=1000,q=2,M=1,rho=1.0 (VAF 0.2 -> CCF 0.4). Landau (2013).
    [Test]
    public void ClassifyClonality_Ccf04_IsSubclonal()
    {
        OncologyAnalyzer.ClonalityResult result =
            OncologyAnalyzer.ClassifyClonality(new[] { V(200, 1000, 2) }, purity: 1.0);

        OncologyAnalyzer.ClonalityCall call = result.Calls[0];
        Assert.Multiple(() =>
        {
            Assert.That(call.Status, Is.EqualTo(OncologyAnalyzer.ClonalityStatus.Subclonal),
                "low CCF ⇒ subclonal (Landau 2013)");
            Assert.That(call.Ccf, Is.EqualTo(0.401198).Within(1e-6), "posterior-mean CCF ≈ 0.4");
            Assert.That(call.ProbabilityClonal, Is.EqualTo(0.0).Within(1e-9),
                "no mass above CCF 0.95");
        });
    }

    // M5 — Multiplicity M=2 raises CCF: a=100,N=100,q=2,M=2,rho=1.0 -> f(c)=c. Satas (2021) Eq. 1.
    [Test]
    public void ClassifyClonality_MultiplicityTwo_RaisesCcfToClonal()
    {
        OncologyAnalyzer.ClonalityResult result =
            OncologyAnalyzer.ClassifyClonality(new[] { V(100, 100, 2, 2) }, purity: 1.0);

        OncologyAnalyzer.ClonalityCall call = result.Calls[0];
        Assert.Multiple(() =>
        {
            Assert.That(call.Status, Is.EqualTo(OncologyAnalyzer.ClonalityStatus.Clonal),
                "M=2 doubles the per-CCF allele fraction ⇒ clonal (Satas 2021)");
            Assert.That(call.Ccf, Is.EqualTo(0.994330).Within(1e-6), "posterior-mean CCF");
            Assert.That(call.ProbabilityClonal, Is.EqualTo(0.998016).Within(1e-6),
                "posterior mass above CCF 0.95");
        });
    }

    // M5b — same VAF (alt/total) at M=1 would NOT be clonal: confirms multiplicity is load-bearing.
    [Test]
    public void ClassifyClonality_SameVafMultiplicityOne_IsSubclonal()
    {
        // a=100,N=100 with M=1, rho=1.0: f(c)=c so a=N would be clonal; use a lower count to isolate M.
        // a=50,N=100,q=2,M=1,rho=1.0 (VAF 0.5 -> CCF 1 point est, but shallow N) vs M=2 case above.
        OncologyAnalyzer.ClonalityResult result =
            OncologyAnalyzer.ClassifyClonality(new[] { V(50, 100, 2, 1) }, purity: 1.0);

        Assert.That(result.Calls[0].Status, Is.EqualTo(OncologyAnalyzer.ClonalityStatus.Subclonal),
            "VAF 0.5 at N=100, M=1 has P(CCF>0.95)≈0.443 < 0.5 ⇒ subclonal (Landau 2013)");
    }

    // M6 — Counts partition (INV-1): one clonal (M1) + one subclonal (M3).
    [Test]
    public void ClassifyClonality_MixedSet_CountsPartitionTotal()
    {
        var variants = new[] { V(300, 300, 2), V(240, 1000, 2) };

        OncologyAnalyzer.ClonalityResult result = OncologyAnalyzer.ClassifyClonality(variants, purity: 0.8);

        Assert.Multiple(() =>
        {
            Assert.That(result.ClonalCount, Is.EqualTo(1), "exactly one clonal variant");
            Assert.That(result.SubclonalCount, Is.EqualTo(1), "exactly one subclonal variant");
            Assert.That(result.ClonalCount + result.SubclonalCount, Is.EqualTo(variants.Length),
                "INV-1: clonal_count + subclonal_count = total_variants");
        });
    }

    // M7 — ClonalFraction (INV-2): 3 clonal + 1 subclonal -> 0.75.
    [Test]
    public void ClassifyClonality_ThreeClonalOneSubclonal_ClonalFractionIsThreeQuarters()
    {
        // All variants evaluated at purity 1.0 (f(c)=M*c/2*... ); 3 clonal + 1 subclonal -> fraction 0.75.
        var clonalSet = new[]
        {
            V(300, 300, 2),    // clonal
            V(990, 1000, 2),   // clonal (VAF 0.99, f(c)=c at rho=1.0)
            V(100, 100, 2, 2), // clonal (M=2)
            V(200, 1000, 2),   // subclonal (VAF 0.2 -> CCF 0.4)
        };

        OncologyAnalyzer.ClonalityResult result = OncologyAnalyzer.ClassifyClonality(clonalSet, purity: 1.0);

        Assert.Multiple(() =>
        {
            Assert.That(result.ClonalCount, Is.EqualTo(3), "three clonal variants");
            Assert.That(result.SubclonalCount, Is.EqualTo(1), "one subclonal variant");
            Assert.That(result.ClonalFraction, Is.EqualTo(0.75).Within(1e-12),
                "INV-2: ClonalFraction = 3/4");
        });
    }

    // INV-3/INV-4 bounds on a representative call.
    [Test]
    public void ClassifyClonality_AnyCall_CcfAndProbabilityWithinBounds()
    {
        OncologyAnalyzer.ClonalityResult result =
            OncologyAnalyzer.ClassifyClonality(new[] { V(240, 1000, 2) }, purity: 0.8);

        OncologyAnalyzer.ClonalityCall call = result.Calls[0];
        Assert.Multiple(() =>
        {
            Assert.That(call.Ccf, Is.InRange(0.01, 1.0), "INV-3: CCF estimate ∈ [0.01, 1]");
            Assert.That(call.ProbabilityClonal, Is.InRange(0.0, 1.0), "INV-4: probability ∈ [0, 1]");
        });
    }

    #endregion

    #region ClassifyClonality — validation and edge cases

    // S1 — null variants.
    [Test]
    public void ClassifyClonality_NullVariants_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.ClassifyClonality(null!, 0.8));
    }

    // S2 — purity out of range.
    [Test]
    public void ClassifyClonality_PurityOutOfRange_Throws()
    {
        var v = new[] { V(50, 100, 2) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyClonality(v, 0.0),
                "purity 0 is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyClonality(v, 1.5),
                "purity > 1 is invalid");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.ClassifyClonality(v, double.NaN),
                "purity NaN is invalid");
        });
    }

    // S3 — invalid read counts.
    [Test]
    public void ClassifyClonality_InvalidReadCounts_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClassifyClonality(new[] { V(1, 0, 2) }, 0.8),
                "TotalReads < 1 is invalid");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClassifyClonality(new[] { V(120, 100, 2) }, 0.8),
                "AltReads > TotalReads is invalid");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClassifyClonality(new[] { V(-1, 100, 2) }, 0.8),
                "AltReads < 0 is invalid (lower bound of [0, N])");
        });
    }

    // S4 — invalid copy number / multiplicity.
    [Test]
    public void ClassifyClonality_InvalidCopyNumberOrMultiplicity_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClassifyClonality(new[] { V(50, 100, 0) }, 0.8),
                "LocalCopyNumber < 1 is invalid");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClassifyClonality(new[] { V(50, 100, 2, 3) }, 0.8),
                "Multiplicity > LocalCopyNumber is invalid");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.ClassifyClonality(new[] { V(50, 100, 2, 0) }, 0.8),
                "Multiplicity < 1 is invalid (lower bound of [1, q])");
        });
    }

    // C1 — empty variant set.
    [Test]
    public void ClassifyClonality_EmptySet_ReturnsEmptyWithZeroCounts()
    {
        OncologyAnalyzer.ClonalityResult result =
            OncologyAnalyzer.ClassifyClonality(Array.Empty<OncologyAnalyzer.ClonalityVariant>(), 0.8);

        Assert.Multiple(() =>
        {
            Assert.That(result.Calls, Is.Empty, "no calls for an empty set");
            Assert.That(result.ClonalCount, Is.EqualTo(0), "zero clonal");
            Assert.That(result.SubclonalCount, Is.EqualTo(0), "zero subclonal");
            Assert.That(result.ClonalFraction, Is.EqualTo(0.0), "ClonalFraction = 0 for empty set");
        });
    }

    #endregion

    #region IdentifyClonalMutations

    // M8 — strict > 0.95 threshold; boundary 0.95 excluded. Landau (2013).
    [Test]
    public void IdentifyClonalMutations_StrictThreshold_ReturnsCorrectIndices()
    {
        var ccf = new[] { 0.96, 0.95, 1.00, 0.50, 0.951 };

        IReadOnlyList<int> clonal = OncologyAnalyzer.IdentifyClonalMutations(ccf);

        Assert.That(clonal, Is.EqualTo(new[] { 0, 2, 4 }),
            "clonal iff CCF > 0.95 (strict): 0.95 boundary excluded, 0.951 included (Landau 2013)");
    }

    // C2 — empty CCF list.
    [Test]
    public void IdentifyClonalMutations_Empty_ReturnsEmpty()
    {
        Assert.That(OncologyAnalyzer.IdentifyClonalMutations(Array.Empty<double>()), Is.Empty,
            "no clonal indices for an empty input");
    }

    // S5 — null CCF list.
    [Test]
    public void IdentifyClonalMutations_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.IdentifyClonalMutations(null!));
    }

    // S6 — CCF out of range.
    [Test]
    public void IdentifyClonalMutations_CcfOutOfRange_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.IdentifyClonalMutations(new[] { 1.2 }),
                "CCF > 1 is invalid");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.IdentifyClonalMutations(new[] { double.NaN }),
                "CCF NaN is invalid");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.IdentifyClonalMutations(new[] { -0.1 }),
                "CCF < 0 is invalid (lower bound of [0, 1])");
        });
    }

    #endregion
}
