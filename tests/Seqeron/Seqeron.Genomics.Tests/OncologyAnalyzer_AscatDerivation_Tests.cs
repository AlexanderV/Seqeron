// ONCO-ASCAT-001 — Upstream allele-specific derivation (segmentation, purity/ploidy fit, multiplicity)
// Evidence: docs/Evidence/ONCO-ASCAT-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-ASCAT-001.md
// Source: Van Loo P et al. (2010). PNAS 107(39):16910-16915. https://doi.org/10.1073/pnas.1009843107
//         VanLoo-lab/ascat, ASCAT/R/ascat.runAscat.R (nA/nB equations + goodness of fit).
//         McGranahan N et al. (2016). Science 351(6280):1463-1469. https://doi.org/10.1126/science.aaf1490
//         Zheng L et al. (2022). Bioinformatics 38(15):3677-3683. https://doi.org/10.1093/bioinformatics/btac440
//
// Planted-truth inputs are synthesised in-test by inverting the ASCAT forward model:
//   denom = rho*n + 2*(1-rho);  D = rho*psi + 2*(1-rho)
//   logR r = log2(denom / D);   BAF b = (rho*nB + (1-rho)) / denom
// These are the exact algebraic inverse of the two cited nA/nB equations (gamma=1). Expected
// purity/ploidy/CN/multiplicity/CCF values are computed independently of the implementation.

using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_AscatDerivation_Tests
{
    private const double Gamma = 1.0;

    // ---- Planted-truth forward model (ASCAT inverse) ----

    private static (double LogR, double Baf) Forward(int nA, int nB, double rho, double psi)
    {
        int n = nA + nB;
        double denom = rho * n + 2.0 * (1.0 - rho);
        double d = rho * psi + 2.0 * (1.0 - rho);
        double r = Math.Log2(denom / d);
        double b = (rho * nB + (1.0 - rho)) / denom;
        return (r, b);
    }

    // Builds per-locus measurements for a list of (chrom, nA, nB) integer segments, each replicated
    // into `lociPerSegment` adjacent loci 1000 bp apart, given the planted rho/psi.
    private static List<OncologyAnalyzer.AlleleSpecificLocus> SynthesiseLoci(
        IReadOnlyList<(string Chrom, int NA, int NB)> segments, double rho, double psi, int lociPerSegment = 5)
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        long pos = 1000;
        foreach (var (chrom, nA, nB) in segments)
        {
            (double r, double b) = Forward(nA, nB, rho, psi);
            for (int i = 0; i < lociPerSegment; i++)
            {
                loci.Add(new OncologyAnalyzer.AlleleSpecificLocus(chrom, pos, r, b));
                pos += 1000;
            }
        }

        return loci;
    }

    // Diploid planted genome (length-weighted mean total CN = 2.2). rho0 = 0.80.
    private const double PlantedPurity = 0.80;
    private const double PlantedPloidy = 2.2;

    private static readonly (string Chrom, int NA, int NB)[] PlantedSegments =
    {
        ("1", 1, 1), // balanced diploid, b=0.5
        ("1", 2, 0), // copy-neutral LOH, b=0.1
        ("1", 1, 1),
        ("1", 2, 1), // gain
        ("1", 1, 1),
    };

    #region SegmentAlleleSpecific

    // M1 — two clear logR levels on chr1 then a chromosome change: 3 segments at the planted boundaries.
    [Test]
    public void SegmentAlleleSpecific_TwoLevelsAndChromChange_RecoversThreeSegments()
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        // chr1: level 0.0 (5 loci), then level 1.0 (5 loci) -> mean-shift split
        for (int i = 0; i < 5; i++) loci.Add(new OncologyAnalyzer.AlleleSpecificLocus("1", 1000 + i * 1000, 0.0, 0.5));
        for (int i = 0; i < 5; i++) loci.Add(new OncologyAnalyzer.AlleleSpecificLocus("1", 6000 + i * 1000, 1.0, 0.5));
        // chr2: level 0.0 -> chromosome change forces a new segment
        for (int i = 0; i < 5; i++) loci.Add(new OncologyAnalyzer.AlleleSpecificLocus("2", 1000 + i * 1000, 0.0, 0.5));

        IReadOnlyList<OncologyAnalyzer.AlleleSpecificSegmentSummary> segs =
            OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.5, minLociPerSegment: 1);

        Assert.Multiple(() =>
        {
            Assert.That(segs.Count, Is.EqualTo(3), "Two logR levels on chr1 plus a chromosome change yield 3 segments.");
            Assert.That(segs[0].Chromosome, Is.EqualTo("1"), "First segment is the chr1 low level.");
            Assert.That(segs[0].MeanLogR, Is.EqualTo(0.0).Within(1e-12), "First segment mean logR is the planted 0.0.");
            Assert.That(segs[1].MeanLogR, Is.EqualTo(1.0).Within(1e-12), "Second segment mean logR is the planted 1.0.");
            Assert.That(segs[2].Chromosome, Is.EqualTo("2"), "Third segment is chr2 after the chromosome change.");
        });
    }

    // C1 — one locus per chromosome -> one segment per chromosome, LocusCount = 1.
    [Test]
    public void SegmentAlleleSpecific_SingleLocusPerChromosome_OneSegmentEach()
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>
        {
            new("1", 1000, 0.0, 0.5),
            new("2", 1000, 0.0, 0.5),
        };

        IReadOnlyList<OncologyAnalyzer.AlleleSpecificSegmentSummary> segs =
            OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2);

        Assert.Multiple(() =>
        {
            Assert.That(segs.Count, Is.EqualTo(2), "One locus per chromosome gives one segment per chromosome.");
            Assert.That(segs[0].LocusCount, Is.EqualTo(1), "Single-locus segment has LocusCount 1.");
            Assert.That(segs[0].Start, Is.EqualTo(segs[0].End), "A single-locus segment has Start == End.");
        });
    }

    // M11 — invalid arguments throw.
    [Test]
    public void SegmentAlleleSpecific_InvalidArguments_Throw()
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus> { new("1", 1000, 0.0, 0.5) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.SegmentAlleleSpecific(null!, 0.2), "Null loci must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.SegmentAlleleSpecific(loci, 0.0), "Non-positive threshold must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.SegmentAlleleSpecific(loci, 0.2, minLociPerSegment: 0), "minLoci < 1 must throw.");
        });
    }

    #endregion

    #region FitPurityPloidy

    // M2/M3/M4 — grid fit recovers rho0=0.80, psi0=2.2, the integer (nA,nB), and GoF ~= 100%.
    [Test]
    public void FitPurityPloidy_PlantedDiploid_RecoversPurityPloidyAndIntegerCopyNumbers()
    {
        List<OncologyAnalyzer.AlleleSpecificLocus> loci = SynthesiseLoci(PlantedSegments, PlantedPurity, PlantedPloidy);
        IReadOnlyList<OncologyAnalyzer.AlleleSpecificSegmentSummary> summaries =
            OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);

        OncologyAnalyzer.PurityPloidyFit fit = OncologyAnalyzer.FitPurityPloidy(
            summaries, purityMin: 0.05, purityMax: 1.0, purityStep: 0.01,
            ploidyMin: 1.5, ploidyMax: 5.0, ploidyStep: 0.05, gamma: Gamma);

        Assert.Multiple(() =>
        {
            // M2 — recovered within one grid step of the planted values.
            Assert.That(fit.Purity, Is.EqualTo(PlantedPurity).Within(0.01), "Recovers planted purity rho0 = 0.80.");
            Assert.That(fit.Ploidy, Is.EqualTo(PlantedPloidy).Within(0.05), "Recovers planted ploidy psi0 = 2.2.");
            // M4 — exact integer fit at the truth => distance ~ 0 => GoF ~ 100%.
            Assert.That(fit.GoodnessOfFit, Is.EqualTo(100.0).Within(1e-6), "GoF is ~100% at the integer-CN truth.");
        });

        // M3 — integer copy numbers per segment (major >= minor). Planted totals: 2,2,2,3,2.
        var expected = new (int Major, int Minor)[] { (1, 1), (2, 0), (1, 1), (2, 1), (1, 1) };
        Assert.That(fit.Segments.Count, Is.EqualTo(expected.Length), "One emitted segment per planted segment.");
        for (int i = 0; i < expected.Length; i++)
        {
            Assert.Multiple(() =>
            {
                Assert.That(fit.Segments[i].MajorCopyNumber, Is.EqualTo(expected[i].Major),
                    $"Segment {i} recovers planted major CN {expected[i].Major}.");
                Assert.That(fit.Segments[i].MinorCopyNumber, Is.EqualTo(expected[i].Minor),
                    $"Segment {i} recovers planted minor CN {expected[i].Minor}.");
            });
        }
    }

    // S3 — triploid planted genome (psi0 = 3.0) is recovered (aneuploidy).
    [Test]
    public void FitPurityPloidy_PlantedTriploid_RecoversPloidyThree()
    {
        var triploid = new (string Chrom, int NA, int NB)[]
        {
            ("1", 2, 1), // total 3
            ("1", 2, 1),
            ("1", 3, 0), // total 3, LOH
            ("1", 2, 1),
        };
        double psi0 = triploid.Average(s => s.NA + s.NB); // 3.0
        List<OncologyAnalyzer.AlleleSpecificLocus> loci = SynthesiseLoci(triploid, PlantedPurity, psi0);
        IReadOnlyList<OncologyAnalyzer.AlleleSpecificSegmentSummary> summaries =
            OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);

        OncologyAnalyzer.PurityPloidyFit fit = OncologyAnalyzer.FitPurityPloidy(summaries, gamma: Gamma);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Purity, Is.EqualTo(PlantedPurity).Within(0.01), "Recovers planted purity 0.80.");
            Assert.That(fit.Ploidy, Is.EqualTo(3.0).Within(0.05), "Recovers planted triploid ploidy 3.0.");
            Assert.That(fit.GoodnessOfFit, Is.EqualTo(100.0).Within(1e-6), "GoF ~100% at the integer-CN triploid truth.");
        });
    }

    // S1 — the goodness-of-fit discriminates: GoF at the true (rho,psi) exceeds GoF at a wrong (rho,psi).
    [Test]
    public void FitPurityPloidy_GoodnessOfFit_DiscriminatesTrueFromWrong()
    {
        List<OncologyAnalyzer.AlleleSpecificLocus> loci = SynthesiseLoci(PlantedSegments, PlantedPurity, PlantedPloidy);
        IReadOnlyList<OncologyAnalyzer.AlleleSpecificSegmentSummary> summaries =
            OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);

        // Pin the grid to the true point, then to a deliberately wrong purity, comparing the resulting GoF.
        OncologyAnalyzer.PurityPloidyFit atTruth = OncologyAnalyzer.FitPurityPloidy(
            summaries, purityMin: PlantedPurity, purityMax: PlantedPurity, purityStep: 1.0,
            ploidyMin: PlantedPloidy, ploidyMax: PlantedPloidy, ploidyStep: 1.0);
        OncologyAnalyzer.PurityPloidyFit atWrong = OncologyAnalyzer.FitPurityPloidy(
            summaries, purityMin: 0.30, purityMax: 0.30, purityStep: 1.0,
            ploidyMin: 1.6, ploidyMax: 1.6, ploidyStep: 1.0);

        Assert.That(atTruth.GoodnessOfFit, Is.GreaterThan(atWrong.GoodnessOfFit),
            "GoF must be higher at the true (rho,psi) than at a wrong (rho,psi) — the objective discriminates.");
    }

    // S2 — balanced-only genome (all b=0.5): fit completes and segments fold to BAF = 0.5.
    [Test]
    public void FitPurityPloidy_BalancedOnlyGenome_CompletesWithBalancedSegments()
    {
        var balanced = new (string Chrom, int NA, int NB)[] { ("1", 1, 1), ("1", 2, 2) };
        double psi0 = balanced.Average(s => s.NA + s.NB); // 3.0
        List<OncologyAnalyzer.AlleleSpecificLocus> loci = SynthesiseLoci(balanced, PlantedPurity, psi0);
        IReadOnlyList<OncologyAnalyzer.AlleleSpecificSegmentSummary> summaries =
            OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);

        Assert.That(summaries.All(s => Math.Abs(s.MeanBAF - 0.5) < 1e-9), Is.True,
            "All balanced (1:1, 2:2) segments fold to mean BAF 0.5.");
        OncologyAnalyzer.PurityPloidyFit fit = OncologyAnalyzer.FitPurityPloidy(summaries, gamma: Gamma);
        Assert.That(fit.GoodnessOfFit, Is.LessThanOrEqualTo(100.0 + 1e-9), "GoF percentage never exceeds 100%.");
    }

    // M12 — invalid arguments throw.
    [Test]
    public void FitPurityPloidy_InvalidArguments_Throw()
    {
        var one = new List<OncologyAnalyzer.AlleleSpecificSegmentSummary>
        {
            new("1", 1000, 2000, 0.0, 0.5, 5),
        };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.FitPurityPloidy(null!), "Null segments must throw.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.FitPurityPloidy(new List<OncologyAnalyzer.AlleleSpecificSegmentSummary>()),
                "Empty segments must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.FitPurityPloidy(one, purityMin: 0.0), "purityMin <= 0 must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.FitPurityPloidy(one, purityMin: 0.5, purityMax: 0.4), "purityMax < min must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.FitPurityPloidy(one, ploidyMin: 0.0), "ploidyMin <= 0 must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.FitPurityPloidy(one, gamma: 0.0), "gamma <= 0 must throw.");
        });
    }

    #endregion

    #region DeriveMultiplicity

    // M5 — VAF=0.40 from m=1 on CN=2 (1+1), rho=0.80: n_mut = 0.40*(0.8*2+2*0.2)/0.8 = 1.0 -> m=1.
    [Test]
    public void DeriveMultiplicity_ClonalDiploidSingleCopy_ReturnsOne()
    {
        int m = OncologyAnalyzer.DeriveMultiplicity(vaf: 0.40, purity: 0.80, totalCopyNumber: 2, majorCopyNumber: 1);
        Assert.That(m, Is.EqualTo(1), "n_mut = 0.40·(0.8·2+2·0.2)/0.8 = 1.0, rounds to multiplicity 1.");
    }

    // M6 — VAF=4/7 from m=2 on CN=3 (major=2), rho=0.80: n_mut = (4/7)*(0.8*3+0.4)/0.8 = 2.0 -> m=2.
    [Test]
    public void DeriveMultiplicity_ClonalGainTwoCopies_ReturnsTwo()
    {
        double vaf = 2.0 * 0.80 / (3.0 * 0.80 + 2.0 * 0.20); // synth VAF for m=2, CN=3, ccf=1 => 4/7
        int m = OncologyAnalyzer.DeriveMultiplicity(vaf, purity: 0.80, totalCopyNumber: 3, majorCopyNumber: 2);
        Assert.That(m, Is.EqualTo(2), "n_mut = (4/7)·(0.8·3+2·0.2)/0.8 = 2.0, rounds to multiplicity 2.");
    }

    // M7 — high VAF whose n_mut exceeds the major CN is clamped down to majorCopyNumber.
    [Test]
    public void DeriveMultiplicity_AboveMajorCopyNumber_ClampsToMajor()
    {
        int m = OncologyAnalyzer.DeriveMultiplicity(vaf: 1.0, purity: 1.0, totalCopyNumber: 4, majorCopyNumber: 2);
        Assert.That(m, Is.EqualTo(2), "n_mut = 1.0·4/1.0 = 4 > major CN 2, clamped to 2.");
    }

    // M8 — tiny VAF whose n_mut rounds to 0 is clamped up to 1 (an observed variant has >= 1 copy).
    [Test]
    public void DeriveMultiplicity_RoundsToZero_ClampsToOne()
    {
        int m = OncologyAnalyzer.DeriveMultiplicity(vaf: 0.0, purity: 0.80, totalCopyNumber: 2, majorCopyNumber: 1);
        Assert.That(m, Is.EqualTo(1), "n_mut = 0 rounds to 0 but is clamped up to 1.");
    }

    // M10 — invalid arguments throw.
    [Test]
    public void DeriveMultiplicity_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DeriveMultiplicity(1.5, 0.8, 2, 1), "VAF > 1 must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DeriveMultiplicity(0.4, 0.0, 2, 1), "purity <= 0 must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DeriveMultiplicity(0.4, 0.8, 0, 1), "totalCopyNumber < 1 must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DeriveMultiplicity(0.4, 0.8, 2, 3), "majorCopyNumber > total must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.DeriveMultiplicity(0.4, 0.8, 2, 0), "majorCopyNumber < 1 must throw.");
        });
    }

    #endregion

    #region End-to-end CCF (M9)

    // M9 — fit -> derive CN + multiplicity -> EstimateCcf on a planted clonal mutation yields CCF = 1.0.
    [Test]
    public void EndToEnd_PlantedClonalMutation_CcfEqualsOne()
    {
        // Recover purity and the integer copy-number segments from the planted per-locus signal.
        List<OncologyAnalyzer.AlleleSpecificLocus> loci = SynthesiseLoci(PlantedSegments, PlantedPurity, PlantedPloidy);
        IReadOnlyList<OncologyAnalyzer.AlleleSpecificSegmentSummary> summaries =
            OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);
        OncologyAnalyzer.PurityPloidyFit fit = OncologyAnalyzer.FitPurityPloidy(summaries, gamma: Gamma);

        // Take the balanced diploid segment (1+1, total CN 2). A clonal m=1 mutation there:
        OncologyAnalyzer.AlleleSpecificSegment seg = fit.Segments[0];
        int total = seg.MajorCopyNumber + seg.MinorCopyNumber; // 2
        // Synthesise the clonal VAF: VAF = m·CCF·rho / (CN·rho + 2(1-rho)) with m=1, CCF=1.
        double vaf = 1.0 * 1.0 * fit.Purity / (total * fit.Purity + 2.0 * (1.0 - fit.Purity));

        int m = OncologyAnalyzer.DeriveMultiplicity(vaf, fit.Purity, total, seg.MajorCopyNumber);
        OncologyAnalyzer.CcfEstimate ccf = OncologyAnalyzer.EstimateCcf(vaf, fit.Purity, total, m);

        Assert.Multiple(() =>
        {
            Assert.That(seg.MajorCopyNumber + seg.MinorCopyNumber, Is.EqualTo(2), "Recovered segment is diploid (total CN 2).");
            Assert.That(m, Is.EqualTo(1), "Derived multiplicity of the clonal single-copy mutation is 1.");
            Assert.That(ccf.RawCcf, Is.EqualTo(1.0).Within(1e-9), "End-to-end CCF of the planted clonal mutation is 1.0.");
        });
    }

    #endregion
}
