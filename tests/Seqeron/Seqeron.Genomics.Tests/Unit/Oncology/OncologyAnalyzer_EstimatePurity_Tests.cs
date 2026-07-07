// ONCO-PURITY-001 — Tumor Purity Estimation
// Evidence: docs/Evidence/ONCO-PURITY-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-PURITY-001.md
// Source: Antonello et al. (2024). CNAqc. Genome Biology 25(1):38. https://doi.org/10.1186/s13059-024-03170-5
//         CNAqc vignette. https://caravagnalab.github.io/CNAqc/articles/CNAqc.html
//         Carter et al. (2012). ABSOLUTE. Nat Biotechnol 30(5):413-421. https://doi.org/10.1038/nbt.2203
//         Shen & Seshan (2016). FACETS. NAR 44(16):e131. https://doi.org/10.1093/nar/gkw520

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_EstimatePurity_Tests
{
    // Helper: a heterozygous somatic SNV observation with a given tumor VAF = alt/total.
    private static OncologyAnalyzer.VariantObservation Het(int alt, int total) =>
        new("1", 100, "A", "T", alt, total, 0, total);

    #region EstimatePurityFromVaf (single-value delegate, closed form ρ = 2·VAF)

    // M1 — CNAqc worked example: purity 60% ⇔ VAF 30% (m=1, n_tot=2 ⇒ v=π/2 ⇒ ρ=2·VAF).
    [Test]
    public void EstimatePurityFromVaf_DiploidHetVaf30Percent_ReturnsPurity60Percent()
    {
        Assert.That(OncologyAnalyzer.EstimatePurityFromVaf(0.30), Is.EqualTo(0.60).Within(1e-10),
            "CNAqc: a clonal het diploid SNV at VAF 0.30 implies purity 0.60 (ρ = 2·VAF)");
    }

    // M2/M3 — boundaries: VAF 0.5 ⇒ purity 1.0 (fully pure); VAF 0 ⇒ purity 0.0.
    [TestCase(0.50, 1.00)]
    [TestCase(0.00, 0.00)]
    [TestCase(0.275, 0.55)] // M5: CNAqc 55% purity ⇔ 27.5% VAF band edge.
    public void EstimatePurityFromVaf_ClosedForm_ReturnsTwiceVaf(double vaf, double expected)
    {
        Assert.That(OncologyAnalyzer.EstimatePurityFromVaf(vaf), Is.EqualTo(expected).Within(1e-10),
            $"ρ = 2·VAF = 2·{vaf} = {expected}");
    }

    // M10 — VAF > 0.5 implies purity > 1 under the diploid het model: rejected.
    [Test]
    public void EstimatePurityFromVaf_VafAbovePointFive_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimatePurityFromVaf(0.60),
            "VAF 0.6 would give purity 1.2 > 1, which is impossible for a het diploid locus");
    }

    // M13 — VAF outside [0,1] is invalid.
    [TestCase(-0.1)]
    [TestCase(1.1)]
    public void EstimatePurityFromVaf_VafOutOfRange_Throws(double vaf)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimatePurityFromVaf(vaf),
            "VAF must be in [0, 1]");
    }

    #endregion

    #region EstimatePurityFromVAF (read-count collection, median of ρ = 2·VAF)

    // M1 — single het diploid SNV with VAF 30/100 = 0.30 ⇒ purity 0.60.
    [Test]
    public void EstimatePurityFromVAF_SingleDiploidHet_ReturnsTwiceVaf()
    {
        double purity = OncologyAnalyzer.EstimatePurityFromVAF(new[] { Het(30, 100) });

        Assert.That(purity, Is.EqualTo(0.60).Within(1e-10),
            "CNAqc: VAF 0.30 ⇒ purity 0.60 (ρ = 2·VAF)");
    }

    // M4 — median aggregation: VAFs {0.10, 0.15, 0.30} ⇒ purities {0.20, 0.30, 0.60}, median 0.30.
    [Test]
    public void EstimatePurityFromVAF_MultipleVariants_ReturnsMedianPurity()
    {
        var variants = new[] { Het(10, 100), Het(15, 100), Het(30, 100) };

        double purity = OncologyAnalyzer.EstimatePurityFromVAF(variants);

        Assert.That(purity, Is.EqualTo(0.30).Within(1e-10),
            "Per-variant purities {0.20, 0.30, 0.60}; robust median = 0.30");
    }

    // S1 — low purity below detection (VAF 0.02 ⇒ purity 0.04) returns without error.
    [Test]
    public void EstimatePurityFromVAF_LowPurity_ReturnsSmallPurityNoError()
    {
        double purity = OncologyAnalyzer.EstimatePurityFromVAF(new[] { Het(2, 100) });

        Assert.That(purity, Is.EqualTo(0.04).Within(1e-10),
            "VAF 0.02 ⇒ purity 0.04 (ρ = 2·VAF); below detection but not an error");
    }

    // M11 — empty input: purity undefined.
    [Test]
    public void EstimatePurityFromVAF_Empty_Throws()
    {
        Assert.Throws<ArgumentException>(
            () => OncologyAnalyzer.EstimatePurityFromVAF(Array.Empty<OncologyAnalyzer.VariantObservation>()),
            "Purity is undefined for an empty variant set");
    }

    // M14 — null input.
    [Test]
    public void EstimatePurityFromVAF_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.EstimatePurityFromVAF(null!),
            "Null variant collection is invalid");
    }

    #endregion

    #region EstimatePurity (allele-specific, π = 2v / [m + v(2 − n_tot)])

    // M6 — diploid het (m=1, n_tot=2), VAF 0.30 ⇒ purity 0.60 (must agree with the VAF-only estimator).
    [Test]
    public void EstimatePurity_DiploidHet_AgreesWithVafOnly()
    {
        double purity = OncologyAnalyzer.EstimatePurity(new[]
        {
            new OncologyAnalyzer.PurityVariant(0.30, Multiplicity: 1, TumorTotalCopyNumber: 2)
        });

        Assert.That(purity, Is.EqualTo(0.60).Within(1e-10),
            "Allele-specific inversion at m=1,n=2 reduces to ρ=2·VAF=0.60");
    }

    // M7/M8 — CNAqc 2:1 segment (n_tot=3) at purity 1.0: m=1 peak at VAF 1/3, m=2 peak at VAF 2/3.
    [Test]
    public void EstimatePurity_TwoToOneSegment_PurePeaks_ReturnOne()
    {
        Assert.Multiple(() =>
        {
            double p1 = OncologyAnalyzer.EstimatePurity(new[]
            {
                new OncologyAnalyzer.PurityVariant(1.0 / 3.0, Multiplicity: 1, TumorTotalCopyNumber: 3)
            });
            double p2 = OncologyAnalyzer.EstimatePurity(new[]
            {
                new OncologyAnalyzer.PurityVariant(2.0 / 3.0, Multiplicity: 2, TumorTotalCopyNumber: 3)
            });

            Assert.That(p1, Is.EqualTo(1.0).Within(1e-10),
                "CNAqc 2:1 peak: m=1, VAF 1/3 ⇒ purity 1.0");
            Assert.That(p2, Is.EqualTo(1.0).Within(1e-10),
                "CNAqc 2:1 peak: m=2, VAF 2/3 ⇒ purity 1.0");
        });
    }

    // M9 — general non-trivial inversion: π=0.5, m=1, n_tot=4 ⇒ v = 0.5/3 = 1/6; inverting recovers 0.5.
    [Test]
    public void EstimatePurity_GeneralAlleleSpecific_RecoversPurity()
    {
        double purity = OncologyAnalyzer.EstimatePurity(new[]
        {
            new OncologyAnalyzer.PurityVariant(1.0 / 6.0, Multiplicity: 1, TumorTotalCopyNumber: 4)
        });

        Assert.That(purity, Is.EqualTo(0.5).Within(1e-10),
            "π = 2·(1/6) / [1 + (1/6)(2−4)] = (1/3)/(2/3) = 0.5");
    }

    // S2 — median across allele-specific variants (purities 0.60, 1.0, 0.5 ⇒ median 0.60).
    [Test]
    public void EstimatePurity_MultipleVariants_ReturnsMedian()
    {
        var variants = new[]
        {
            new OncologyAnalyzer.PurityVariant(0.30, 1, 2),       // 0.60
            new OncologyAnalyzer.PurityVariant(1.0 / 3.0, 1, 3),  // 1.0
            new OncologyAnalyzer.PurityVariant(1.0 / 6.0, 1, 4),  // 0.5
        };

        double purity = OncologyAnalyzer.EstimatePurity(variants);

        Assert.That(purity, Is.EqualTo(0.60).Within(1e-10),
            "Per-variant purities {0.60, 1.0, 0.5}; median = 0.60");
    }

    // M12 — invalid copy number / multiplicity.
    [Test]
    public void EstimatePurity_InvalidCopyNumberOrMultiplicity_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimatePurity(new[]
            {
                new OncologyAnalyzer.PurityVariant(0.30, 1, 0)
            }), "n_tot must be ≥ 1");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimatePurity(new[]
            {
                new OncologyAnalyzer.PurityVariant(0.30, 0, 2)
            }), "multiplicity must be ≥ 1");
        });
    }

    // M13 — VAF out of range for the allele-specific overload.
    [Test]
    public void EstimatePurity_VafOutOfRange_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimatePurity(new[]
        {
            new OncologyAnalyzer.PurityVariant(1.5, 1, 2)
        }), "VAF must be in [0, 1]");
    }

    // M10 (allele-specific) — a (VAF, m, n_tot) combination implying purity > 1 is rejected.
    [Test]
    public void EstimatePurity_CombinationImplyingPurityAboveOne_Throws()
    {
        // m=1, n_tot=2, VAF 0.6 ⇒ π = 1.2 > 1.
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimatePurity(new[]
        {
            new OncologyAnalyzer.PurityVariant(0.60, 1, 2)
        }), "VAF 0.6 at m=1,n=2 implies purity 1.2 > 1");
    }

    // M10b — a (VAF, m, n_tot) combination whose VAF is unreachable by any purity in [0,1] (non-positive
    // denominator) is rejected. For m=1, n_tot=4 the forward relation v = π/[2(1−π)+4π] maps π∈[0,1] onto
    // VAF∈[0, 0.25]; VAF 0.9 is impossible, so the inverse denominator m+v(2−n_tot)=1+0.9·(−2)=−0.8 ≤ 0.
    // Evidence: CNAqc/ABSOLUTE forward range of v (Antonello 2024; Carter 2012, F=αs_q/[αq+2(1−α)]).
    [Test]
    public void EstimatePurity_VafUnreachableForCopyState_NonPositiveDenominator_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EstimatePurity(new[]
        {
            new OncologyAnalyzer.PurityVariant(0.90, Multiplicity: 1, TumorTotalCopyNumber: 4)
        }), "VAF 0.9 at m=1,n_tot=4 is unreachable (denominator 1+0.9·(2−4)=−0.8 ≤ 0)");
    }

    // M11/M14 — empty and null.
    [Test]
    public void EstimatePurity_EmptyAndNull_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.EstimatePurity(Array.Empty<OncologyAnalyzer.PurityVariant>()),
                "Purity is undefined for an empty variant set");
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.EstimatePurity(null!),
                "Null variant collection is invalid");
        });
    }

    // C1 — determinism: same input twice yields identical results.
    [Test]
    public void EstimatePurity_RepeatedCalls_AreDeterministic()
    {
        var first = new[]
        {
            new OncologyAnalyzer.PurityVariant(0.30, 1, 2),
            new OncologyAnalyzer.PurityVariant(1.0 / 3.0, 1, 3),
        };
        var second = new[]
        {
            new OncologyAnalyzer.PurityVariant(0.30, 1, 2),
            new OncologyAnalyzer.PurityVariant(1.0 / 3.0, 1, 3),
        };

        double firstResult = OncologyAnalyzer.EstimatePurity(first);
        double secondResult = OncologyAnalyzer.EstimatePurity(second);

        Assert.That(secondResult, Is.EqualTo(firstResult).Within(1e-12),
            "Estimation is deterministic and order/state-independent");
    }

    #endregion
}
