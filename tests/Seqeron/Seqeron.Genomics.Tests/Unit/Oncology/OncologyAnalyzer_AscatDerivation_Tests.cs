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

namespace Seqeron.Genomics.Tests.Unit.Oncology;

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

    // Each planted segment sits on its own chromosome (5 loci < ASCAT kmin = 6), so ascat.aspcf — which
    // SegmentAlleleSpecific runs (F35) — emits exactly one segment per planted segment: a chromosome with fewer than
    // 6 loci is one segment (S-ASPCF-3), and noise-free data has MAD sd = 0 so no in-chromosome breakpoint is placed
    // anyway (S-ASPCF-2).
    private static readonly (string Chrom, int NA, int NB)[] PlantedSegments =
    {
        ("1", 1, 1), // balanced diploid, b=0.5
        ("2", 2, 0), // copy-neutral LOH, b=0.1
        ("3", 1, 1),
        ("4", 2, 1), // gain
        ("5", 1, 1),
    };

    #region SegmentAlleleSpecific

    // M1 (F35) — SegmentAlleleSpecific runs ascat.aspcf with ASCAT's default penalty 70: on the R-locked noisy step
    // track it returns the ascat.aspcf output (same values as M-ASPCF-1), whatever the ignored legacy thresholds.
    [Test]
    public void SegmentAlleleSpecific_NoisyStep_MatchesAscatAspcfDefaultPenalty()
    {
        IReadOnlyList<OncologyAnalyzer.AlleleSpecificSegmentSummary> segs =
            OncologyAnalyzer.SegmentAlleleSpecific(AspcfStepTrack(), logRChangeThreshold: 0.5, minLociPerSegment: 1);

        Assert.Multiple(() =>
        {
            Assert.That(segs.Count, Is.EqualTo(2), "ascat.aspcf (penalty 70): 2 segments.");
            Assert.That((segs[0].Start, segs[0].End, segs[0].LocusCount), Is.EqualTo((1000L, 40000L, 40)), "Segment 1 = loci 1–40.");
            Assert.That((segs[1].Start, segs[1].End, segs[1].LocusCount), Is.EqualTo((41000L, 80000L, 40)), "Segment 2 = loci 41–80.");
            Assert.That(segs[0].MeanLogR, Is.EqualTo(-0.023564999999999999), "ascat.aspcf logR level 1.");
            Assert.That(segs[1].MeanLogR, Is.EqualTo(0.60321000000000002), "ascat.aspcf logR level 2.");
            Assert.That(segs[0].MeanBAF, Is.EqualTo(0.5), "Balanced segment BAF shrunk to exactly 0.5 (ASCAT).");
            Assert.That(segs[1].MeanBAF, Is.EqualTo(0.75129407874999998), "ascat.aspcf BAF level 2 (0.5 + μ).");
            Assert.That(segs, Is.EqualTo(OncologyAnalyzer.SegmentAlleleSpecificAspcf(AspcfStepTrack(), OncologyAnalyzer.AspcfDefaultPenalty)),
                "Identical to SegmentAlleleSpecificAspcf at the ASCAT default penalty.");
        });
    }

    // M1b (F35) — the legacy greedy thresholds have no ASPCF meaning and are ignored: very different values give the
    // identical ascat.aspcf segmentation.
    [Test]
    public void SegmentAlleleSpecific_LegacyThresholds_AreIgnored()
    {
        var a = OncologyAnalyzer.SegmentAlleleSpecific(AspcfStepTrack(), logRChangeThreshold: 0.05, bafChangeThreshold: 0.01, minLociPerSegment: 1);
        var b = OncologyAnalyzer.SegmentAlleleSpecific(AspcfStepTrack(), logRChangeThreshold: 10.0, bafChangeThreshold: 5.0, minLociPerSegment: 50);

        Assert.That(a, Is.EqualTo(b), "Thresholds / minLociPerSegment do not change the ASPCF result.");
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
            // F35: the ASPCF input contract (IsValidAlleleSignal) — finite logR, BAF in [0, 1].
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.SegmentAlleleSpecific(new[] { new OncologyAnalyzer.AlleleSpecificLocus("1", 1, double.NaN, 0.5) }, 0.2),
                "NaN logR must throw.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.SegmentAlleleSpecific(new[] { new OncologyAnalyzer.AlleleSpecificLocus("1", 1, double.PositiveInfinity, 0.5) }, 0.2),
                "Infinite logR must throw.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.SegmentAlleleSpecific(new[] { new OncologyAnalyzer.AlleleSpecificLocus("1", 1, 0.0, -0.3) }, 0.2),
                "BAF < 0 must throw.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.SegmentAlleleSpecific(new[] { new OncologyAnalyzer.AlleleSpecificLocus("1", 1, 0.0, 1.5) }, 0.2),
                "BAF > 1 must throw.");
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
            ("1", 2, 1), // total 3 (one chromosome per planted segment, see PlantedSegments)
            ("2", 2, 1),
            ("3", 3, 0), // total 3, LOH
            ("4", 2, 1),
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
    // Evaluated through the ASCAT rho_manual/psi_manual path (ascat.runAscat.R), which scores a fixed (ρ, ψ).
    [Test]
    public void FitPurityPloidy_GoodnessOfFit_DiscriminatesTrueFromWrong()
    {
        List<OncologyAnalyzer.AlleleSpecificLocus> loci = SynthesiseLoci(PlantedSegments, PlantedPurity, PlantedPloidy);
        IReadOnlyList<OncologyAnalyzer.AlleleSpecificSegmentSummary> summaries =
            OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);

        OncologyAnalyzer.PurityPloidyFit atTruth = OncologyAnalyzer.EvaluatePurityPloidy(summaries, PlantedPurity, PlantedPloidy);
        OncologyAnalyzer.PurityPloidyFit atWrong = OncologyAnalyzer.EvaluatePurityPloidy(summaries, 0.30, 1.6);

        Assert.Multiple(() =>
        {
            Assert.That(atTruth.GoodnessOfFit, Is.EqualTo(100.0).Within(1e-9), "d = 0 at the integer-CN truth ⇒ GoF 100 %.");
            Assert.That(atTruth.GoodnessOfFit, Is.GreaterThan(atWrong.GoodnessOfFit),
                "GoF must be higher at the true (rho,psi) than at a wrong (rho,psi) — the objective discriminates.");
            Assert.That(atTruth.Psi, Is.EqualTo(PlantedPloidy), "The manual path reports the supplied ψ.");
        });
    }

    // ---- ASCAT reference cases: values produced by the ORIGINAL ascat.runAscat.R (VanLoo-lab/ascat master) ----
    // runASCAT(lrr, baf, lrrsegmented, bafsegmented, "XX", …, gamma = 1, min_ploidy 1.5, max_ploidy 5.5,
    // min_purity 0.1, max_purity 1.05) run under R 4.3 on one probe per heterozygous locus (segment = LocusCount
    // probes, segmented BAF = 1 − mirrored BAF); segments from the seg_raw rounding block. Cross-check script:
    // docs/Evidence/ONCO-ASCAT-001-Evidence.md §"ASCAT R cross-check" (150/150 random genomes identical).

    private static OncologyAnalyzer.AlleleSpecificSegmentSummary Summary(string chrom, double r, double b, int n) =>
        new(chrom, 0, 1000, r, b, n);

    // ASCAT R: rho = 1, psi = 2.7, goodnessOfFit = 99.781420571107006, nonaberrant = FALSE, seg 2:1 ×3.
    // The chrX segment is excluded from the fit (sexchromosomes) but emitted with the diploid model.
    [Test]
    public void FitPurityPloidy_AscatReferenceCase_WithSexChromosome_MatchesRunAscat()
    {
        var segs = new[]
        {
            Summary("X", -0.19799990120548305, 0.7800586872805495, 55),
            Summary("3", 0.097918462573404322, 0.66269643393554578, 13),
            Summary("1", 0.065882791496869667, 0.6385382413628512, 15),
        };

        OncologyAnalyzer.PurityPloidyFit fit = OncologyAnalyzer.FitPurityPloidy(segs);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Purity, Is.EqualTo(1.0), "ASCAT rho (grid value 1.00).");
            Assert.That(fit.Psi, Is.EqualTo(2.7).Within(1e-12), "ASCAT psi.");
            Assert.That(fit.GoodnessOfFit, Is.EqualTo(99.781420571107006).Within(1e-9), "ASCAT goodnessOfFit.");
            Assert.That(fit.IsNonAberrant, Is.False, "ASCAT nonaberrant.");
            Assert.That(fit.Segments.Select(s => (s.MajorCopyNumber, s.MinorCopyNumber)),
                Is.EqualTo(new[] { (2, 1), (2, 1), (2, 1) }), "ASCAT seg_raw nMajor:nMinor.");
            Assert.That(fit.Ploidy, Is.EqualTo(3.0).Within(1e-12), "ASCAT ploidy = mean integer total CN over probes.");
        });
    }

    // ASCAT R: rho = 0.85, psi = 2.2, goodnessOfFit = 99.999772627448223, seg 2:0 2:1 1:1.
    // The pre-fix global-minimum grid search returned rho = 0.72, psi = 4.45 (3:0 4:2 2:2).
    [Test]
    public void FitPurityPloidy_AscatReferenceCase_LocalMinimumAndFilters_MatchRunAscat()
    {
        var segs = new[]
        {
            Summary("3", -0.44374838250380438, 0.90652646528359304, 6),
            Summary("3", 0.36009080318741721, 0.64110353482264759, 51),
            Summary("1", -0.11480994884775791, 0.5, 60),
        };

        OncologyAnalyzer.PurityPloidyFit fit = OncologyAnalyzer.FitPurityPloidy(segs);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Purity, Is.EqualTo(0.85).Within(1e-12), "ASCAT rho.");
            Assert.That(fit.Psi, Is.EqualTo(2.2).Within(1e-12), "ASCAT psi.");
            Assert.That(fit.GoodnessOfFit, Is.EqualTo(99.999772627448223).Within(1e-9), "ASCAT goodnessOfFit.");
            Assert.That(fit.Segments.Select(s => (s.MajorCopyNumber, s.MinorCopyNumber)),
                Is.EqualTo(new[] { (2, 0), (2, 1), (1, 1) }), "ASCAT seg_raw nMajor:nMinor.");
            Assert.That(fit.Ploidy, Is.EqualTo(285.0 / 117.0).Within(1e-12),
                "Probe-weighted mean integer total CN: (2·6 + 3·51 + 2·60)/117.");
        });
    }

    // ASCAT R returns rho = NA ("ASCAT could not find an optimal ploidy and purity value") for both genomes.
    [Test]
    public void FitPurityPloidy_NoAscatOptimum_TryReturnsFalseAndFitThrows()
    {
        var withX = new[]
        {
            Summary("X", 0.2316983774009044, 0.5, 17),
            Summary("3", 0.164318402357189, 0.89674532869091261, 27),
            Summary("1", -0.48381153881768874, 0.8116658364580579, 56),
        };
        var autosomal = new[]
        {
            Summary("1", 0.23763864436633977, 0.88277974701356954, 6),
            Summary("1", 0.2573573024164994, 0.86513206216626082, 28),
            Summary("1", -0.7611160923354886, 0.7525369149308675, 38),
        };

        Assert.Multiple(() =>
        {
            Assert.That(OncologyAnalyzer.TryFitPurityPloidy(withX, out _), Is.False, "ASCAT: rho = NA.");
            Assert.That(OncologyAnalyzer.TryFitPurityPloidy(autosomal, out _), Is.False, "ASCAT: rho = NA.");
            Assert.Throws<InvalidOperationException>(() => OncologyAnalyzer.FitPurityPloidy(autosomal),
                "FitPurityPloidy reports ASCAT's failure as InvalidOperationException.");
        });
    }

    // Sex chromosomes are excluded from the fit (ASCAT autoprobes): adding a wild chrX segment changes nothing
    // but the emitted segment list.
    [Test]
    public void FitPurityPloidy_SexChromosomeSegment_DoesNotAffectFit()
    {
        var auto = new List<OncologyAnalyzer.AlleleSpecificSegmentSummary>
        {
            Summary("3", -0.44374838250380438, 0.90652646528359304, 6),
            Summary("3", 0.36009080318741721, 0.64110353482264759, 51),
            Summary("1", -0.11480994884775791, 0.5, 60),
        };
        var withSex = new List<OncologyAnalyzer.AlleleSpecificSegmentSummary>(auto)
        {
            Summary("chrX", 1.7, 0.97, 500),
            Summary("Y", -2.0, 0.5, 300),
        };

        OncologyAnalyzer.PurityPloidyFit a = OncologyAnalyzer.FitPurityPloidy(auto);
        OncologyAnalyzer.PurityPloidyFit b = OncologyAnalyzer.FitPurityPloidy(withSex);

        Assert.Multiple(() =>
        {
            Assert.That(b.Purity, Is.EqualTo(a.Purity), "ρ unchanged by sex-chromosome segments.");
            Assert.That(b.Psi, Is.EqualTo(a.Psi), "ψ unchanged by sex-chromosome segments.");
            Assert.That(b.GoodnessOfFit, Is.EqualTo(a.GoodnessOfFit), "GoF is computed over autosomes only.");
            Assert.That(b.Segments.Count, Is.EqualTo(5), "Every summary is still emitted.");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.FitPurityPloidy(new[] { Summary("X", 0.0, 0.7, 10) }),
                "No autosomal segment ⇒ nothing to fit.");
        });
    }

    // ASCAT seg_raw rounding (ascat.runAscat.R, verified in R with the verbatim block):
    //  • balanced (BAF 0.5) odd total: ρ=1, ψ=2, r=log2(3/2) ⇒ nA = nB = 1.5 ⇒ limitround rule ⇒ 2:1 (pre-fix 2:2);
    //  • negative allele folded into the other: ρ=0.5, ψ=2, r=0, mirrored BAF 0.9 ⇒ nA=2.6, nB=−0.6 ⇒ 2:0 (pre-fix 3:0).
    [Test]
    public void EvaluatePurityPloidy_AscatSegmentRounding_BalancedOddTotalAndNegativeCorrection()
    {
        var balancedOdd = OncologyAnalyzer.EvaluatePurityPloidy(new[] { Summary("1", Math.Log2(1.5), 0.5, 10) }, 1.0, 2.0);
        var negative = OncologyAnalyzer.EvaluatePurityPloidy(new[] { Summary("1", 0.0, 0.9, 10) }, 0.5, 2.0);

        Assert.Multiple(() =>
        {
            Assert.That((balancedOdd.Segments[0].MajorCopyNumber, balancedOdd.Segments[0].MinorCopyNumber), Is.EqualTo((2, 1)),
                "Balanced segment with odd total 3 ⇒ ASCAT 2:1.");
            Assert.That((negative.Segments[0].MajorCopyNumber, negative.Segments[0].MinorCopyNumber), Is.EqualTo((2, 0)),
                "nB = −0.6 is folded into nA (2.6 − 0.6 = 2) ⇒ ASCAT 2:0.");
        });
    }

    // Segments must carry ≥ 1 probe (ASCAT's length is the probe count) and finite logR / BAF in [0,1].
    [Test]
    public void FitPurityPloidy_MalformedSegments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.FitPurityPloidy(new[] { Summary("1", 0.0, 0.7, 0) }),
                "LocusCount 0 ⇒ no probe weight.");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.FitPurityPloidy(new[] { Summary("1", double.NaN, 0.7, 5) }),
                "NaN logR rejected.");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.FitPurityPloidy(new[] { Summary("1", 0.0, 1.2, 5) }),
                "BAF > 1 rejected.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.EvaluatePurityPloidy(new[] { Summary("1", 0.0, 0.7, 5) }, 1.2, 2.0),
                "Manual purity > 1 rejected.");
        });
    }

    // S2 — balanced-only genome (all b=0.5): fit completes and segments fold to BAF = 0.5.
    [Test]
    public void FitPurityPloidy_BalancedOnlyGenome_CompletesWithBalancedSegments()
    {
        var balanced = new (string Chrom, int NA, int NB)[] { ("1", 1, 1), ("2", 2, 2) };
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

    // M11 — exact .5 ties round half-to-even, as facets-suite expected_mutant_copies (R round = IEC 60559).
    // All inputs are dyadic, so n_mut is the exact double shown. R 4.x reference output:
    //   expected_mutant_copies(0.625, 4, 1.0) = 2   (mu = 2.5)
    //   expected_mutant_copies(0.375, 2, 0.5) = 2   (mu = 1.5)
    //   expected_mutant_copies(0.875, 4, 1.0) = 4   (mu = 3.5)
    //   expected_mutant_copies(0.5625, 8, 1.0) = 4  (mu = 4.5)
    //   expected_mutant_copies(0.125, 4, 1.0) = 1   (mu = 0.5 < 1 -> 1)
    [TestCase(0.625, 1.0, 4, 4, 2)]
    [TestCase(0.375, 0.5, 2, 2, 2)]
    [TestCase(0.875, 1.0, 4, 4, 4)]
    [TestCase(0.5625, 1.0, 8, 8, 4)]
    [TestCase(0.125, 1.0, 4, 4, 1)]
    public void DeriveMultiplicity_ExactHalfTie_RoundsHalfToEven_MatchesFacetsSuite(
        double vaf, double purity, int totalCopyNumber, int majorCopyNumber, int expected)
    {
        int m = OncologyAnalyzer.DeriveMultiplicity(vaf, purity, totalCopyNumber, majorCopyNumber);
        Assert.That(m, Is.EqualTo(expected),
            "facets-suite expected_mutant_copies: mu < 1 -> 1, then R round() (half-to-even).");
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

    #region SegmentAlleleSpecificAspcf (ASCAT ascat.aspcf / fastAspcf / aspcfpart)

    // Deterministic noisy tracks (80 loci on chr1; 4-decimal literals shared with the R cross-check).
    // Expected values are the output of the ORIGINAL ascat.aspcf (VanLoo-lab/ascat master, R 4.3) on the same
    // loci with Germline_BAF = 0.5 (all heterozygous) — see Evidence §"ASCAT R cross-check" (60/60 random
    // genomes, 681 segments identical to ≤ 5e-16).
    private static readonly double[] AspcfLogRNoise = { -0.0453, -0.0731, 0.1381, -0.0817, -0.0632, 0.1167, -0.0624, -0.1293, 0.0513, -0.0016, -0.0042, -0.039, -0.0982, -0.1934, 0.1133, -0.0558, 0.0104, 0.0412, -0.0503, -0.16, -0.2202, 0.1674, -0.0273, -0.1205, 0.0077, -0.0957, -0.0398, -0.0175, -0.0372, -0.0352, 0.0519, 0.0563, 0.0991, -0.12, -0.1088, -0.0052, 0.0256, 0.0508, -0.0206, 0.0331, -0.0893, -0.1283, -0.1238, -0.0852, -0.1075, -0.0466, -0.0233, 0.0351, 0.0796, -0.0397, 0.1226, -0.1288, -0.0112, 0.0451, -0.1349, 0.0402, 0.1099, 0.1532, -0.0676, 0.0155, 0.0646, 0.1728, -0.1561, 0.0349, -0.1101, -0.0692, -0.0067, 0.0384, 0.2818, -0.0153, 0.1269, -0.0723, -0.0398, 0.0733, -0.0551, 0.104, -0.0577, 0.1, 0.0129, 0.0861 };
    private static readonly double[] AspcfBafNoise = { 0.0551, -0.043, 0.0313, 0.0063, 0.0189, -0.0093, 0.0282, -0.0104, 0.0366, 0.0327, -0.0089, 0.0162, 0.0422, -0.0015, -0.0042, 0.0176, 0.0547, 0.0701, -0.0727, 0.0183, -0.0305, -0.0791, -0.0044, -0.0028, -0.0073, -0.0053, -0.0095, -0.0097, 0.0231, 0.0169, 0.0052, -0.024, -0.0048, -0.0044, -0.0019, -0.0138, 0.0083, -0.0535, -0.0493, 0.0487, -0.0337, -0.0589, 0.0278, -0.044, -0.0299, -0.0028, 0.0535, 0.0257, -0.0473, -0.0817, 0.0147, 0.0468, -0.0068, -0.0375, -0.0206, 0.0489, -0.0407, -0.0614, -0.0107, -0.0145, -0.0128, -0.0009, 0.0529, -0.0075, 0.0448, 0.0199, -0.0101, -0.0299, -0.0394, 0.0183, 0.0099, 0.0669, 0.004, -0.0395, -0.0178, -0.0445, 0.0081, 0.0025, 0.0197, 0.0157 };
    private static readonly int[] AspcfBafSign = { 1, -1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, 1, -1 };

    // Dataset 1: logR 0 → 0.6 and BAF balanced → imbalanced (0.5 ± 0.25) at locus 41.
    internal static List<OncologyAnalyzer.AlleleSpecificLocus> AspcfStepTrack(double logRShift = 0.0)
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        for (int i = 0; i < 80; i++)
        {
            double r = Math.Round(AspcfLogRNoise[i] + (i >= 40 ? 0.6 : 0.0), 4) + logRShift;
            double b = Math.Round(0.5 + (i >= 40 ? AspcfBafSign[i] * 0.25 : 0.0) + AspcfBafNoise[i], 4);
            loci.Add(new OncologyAnalyzer.AlleleSpecificLocus("1", (i + 1) * 1000L, r, b));
        }

        return loci;
    }

    // Dataset 2: identical logR everywhere; BAF balanced → copy-neutral LOH (0.5 ± 0.47) at locus 41.
    internal static List<OncologyAnalyzer.AlleleSpecificLocus> AspcfLohTrack()
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        for (int i = 0; i < 80; i++)
        {
            double b = i < 40 ? 0.5 + AspcfBafNoise[i] : 0.5 + AspcfBafSign[i] * 0.47 + AspcfBafNoise[i];
            loci.Add(new OncologyAnalyzer.AlleleSpecificLocus("1", (i + 1) * 1000L, AspcfLogRNoise[i], Math.Round(Math.Clamp(b, 0.0, 1.0), 4)));
        }

        return loci;
    }

    // M-ASPCF-1 — ASCAT (penalty 70): one breakpoint after locus 40; logR level = mean raw logR; the balanced half's
    // BAF is shrunk to exactly 0.5 (sqrt(sd2² + μ²) < 2·sd2), the imbalanced half's is 0.5 + mean|b − 0.5|.
    [Test]
    public void SegmentAlleleSpecificAspcf_NoisyStep_MatchesAscatAspcf()
    {
        var segs = OncologyAnalyzer.SegmentAlleleSpecificAspcf(AspcfStepTrack(), penalty: 70.0);

        Assert.Multiple(() =>
        {
            Assert.That(segs.Count, Is.EqualTo(2), "ascat.aspcf: 2 segments.");
            Assert.That((segs[0].Start, segs[0].End, segs[0].LocusCount), Is.EqualTo((1000L, 40000L, 40)), "Segment 1 = loci 1–40.");
            Assert.That((segs[1].Start, segs[1].End, segs[1].LocusCount), Is.EqualTo((41000L, 80000L, 40)), "Segment 2 = loci 41–80.");
            Assert.That(segs[0].MeanLogR, Is.EqualTo(-0.023564999999999999), "ascat.aspcf logR level 1.");
            Assert.That(segs[1].MeanLogR, Is.EqualTo(0.60321000000000002), "ascat.aspcf logR level 2.");
            Assert.That(segs[0].MeanBAF, Is.EqualTo(0.5), "Balanced segment BAF shrunk to exactly 0.5 (ASCAT).");
            Assert.That(segs[1].MeanBAF, Is.EqualTo(0.75129407874999998), "ascat.aspcf BAF level 2 (0.5 + μ).");
        });
    }

    // M-ASPCF-2 — copy-neutral LOH vs balanced with the same logR: the joint (logR + BAF) cost still splits.
    [Test]
    public void SegmentAlleleSpecificAspcf_SameLogRDifferentBaf_SplitsOnBaf()
    {
        var segs = OncologyAnalyzer.SegmentAlleleSpecificAspcf(AspcfLohTrack(), penalty: 70.0);

        Assert.Multiple(() =>
        {
            Assert.That(segs.Count, Is.EqualTo(2), "ascat.aspcf: 2 segments despite identical logR.");
            Assert.That(segs[0].LocusCount, Is.EqualTo(40), "Breakpoint after locus 40.");
            Assert.That(segs[0].MeanBAF, Is.EqualTo(0.5), "Balanced half: BAF 0.5.");
            Assert.That(segs[1].MeanBAF, Is.EqualTo(0.96697500000000003), "LOH half: ascat.aspcf BAF 0.966975.");
            Assert.That(segs[1].MeanLogR, Is.EqualTo(0.003210000000000001), "ascat.aspcf logR level of the LOH half.");
        });
    }

    // S-ASPCF-1 — a huge penalty collapses to one segment (ascat.aspcf, penalty 1e6: logR 0.2898225, BAF 0.637905789375).
    [Test]
    public void SegmentAlleleSpecificAspcf_PenaltyControlsSegmentCount()
    {
        var big = OncologyAnalyzer.SegmentAlleleSpecificAspcf(AspcfStepTrack(), penalty: 1e6);
        var def = OncologyAnalyzer.SegmentAlleleSpecificAspcf(AspcfStepTrack());

        Assert.Multiple(() =>
        {
            Assert.That(big.Count, Is.EqualTo(1), "A very large penalty forces a single segment (no breakpoints).");
            Assert.That(big[0].MeanLogR, Is.EqualTo(0.28982249999999998), "ascat.aspcf single-segment logR.");
            Assert.That(big[0].MeanBAF, Is.EqualTo(0.63790578937499998), "ascat.aspcf single-segment BAF.");
            Assert.That(def.Count, Is.EqualTo(2), "The ASCAT default penalty (70) recovers the two levels.");
            Assert.That(OncologyAnalyzer.AspcfDefaultPenalty, Is.EqualTo(70.0), "ascat.aspcf(..., penalty = 70).");
        });
    }

    // S-ASPCF-2 — noise-free data: every window has MAD sd = 0, so ASCAT's fastAspcf skips it (sd.valid fails)
    // and places no breakpoint. ascat.aspcf on 10 × 0.0 then 10 × 1.0 (BAF 0.5), penalty 0.5: one segment,
    // logR 0.5, BAF 0.5.
    [Test]
    public void SegmentAlleleSpecificAspcf_NoiseFreeTrack_NoBreakpointAsAscat()
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        for (int i = 0; i < 10; i++) loci.Add(new OncologyAnalyzer.AlleleSpecificLocus("1", 1000 + i * 1000, 0.0, 0.5));
        for (int i = 0; i < 10; i++) loci.Add(new OncologyAnalyzer.AlleleSpecificLocus("1", 11000 + i * 1000, 1.0, 0.5));

        var segs = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, penalty: 0.5);

        Assert.Multiple(() =>
        {
            Assert.That(segs.Count, Is.EqualTo(1), "MAD = 0 ⇒ window skipped ⇒ no breakpoint (ascat.aspcf).");
            Assert.That(segs[0].MeanLogR, Is.EqualTo(0.5), "Mean raw logR.");
            Assert.That(segs[0].MeanBAF, Is.EqualTo(0.5), "Balanced BAF.");
        });
    }

    // S-ASPCF-3 — a chromosome with fewer than kmin = 6 loci is one segment with the mean winsorised mirrored BAF
    // (no shrinkage). ascat.aspcf on logR (0.1, 0.3, −0.2, 0.05, 0.4), BAF (0.8, 0.3, 0.75, 0.28, 0.9): logR 0.13, BAF 0.774.
    [Test]
    public void SegmentAlleleSpecificAspcf_FewerThanSixLoci_SingleSegmentAsAscat()
    {
        double[] r = { 0.1, 0.3, -0.2, 0.05, 0.4 };
        double[] b = { 0.8, 0.3, 0.75, 0.28, 0.9 };
        var loci = Enumerable.Range(0, 5).Select(i => new OncologyAnalyzer.AlleleSpecificLocus("1", (i + 1) * 1000L, r[i], b[i])).ToList();

        var segs = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci);

        Assert.Multiple(() =>
        {
            Assert.That(segs.Count, Is.EqualTo(1), "n < 6 ⇒ one segment.");
            Assert.That(segs[0].MeanLogR, Is.EqualTo(0.13), "ascat.aspcf logR.");
            Assert.That(segs[0].MeanBAF, Is.EqualTo(0.77400000000000002), "ascat.aspcf BAF (mean mirrored).");
        });
    }

    // S-ASPCF-4 — chromosome boundaries are never crossed by a segment.
    [Test]
    public void SegmentAlleleSpecificAspcf_ChromosomeBoundary_NeverCrossed()
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        for (int i = 0; i < 5; i++) loci.Add(new OncologyAnalyzer.AlleleSpecificLocus("1", 1000 + i * 1000, 0.0, 0.5));
        for (int i = 0; i < 5; i++) loci.Add(new OncologyAnalyzer.AlleleSpecificLocus("2", 1000 + i * 1000, 0.0, 0.5));

        // Same flat value across both chromosomes; only the contig change can split them.
        var segs = OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, penalty: 100.0);

        Assert.Multiple(() =>
        {
            Assert.That(segs.Count, Is.EqualTo(2), "A segment may not span two chromosomes even at a high penalty.");
            Assert.That(segs[0].Chromosome, Is.EqualTo("1"), "First segment is chr1.");
            Assert.That(segs[1].Chromosome, Is.EqualTo("2"), "Second segment is chr2.");
        });
    }

    // C-ASPCF-1 — invalid arguments throw.
    [Test]
    public void SegmentAlleleSpecificAspcf_InvalidArguments_Throw()
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus> { new("1", 1000, 0.0, 0.5) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.SegmentAlleleSpecificAspcf(null!), "Null loci must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.SegmentAlleleSpecificAspcf(loci, penalty: 0.0), "Non-positive penalty must throw.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.SegmentAlleleSpecificAspcf(new[] { new OncologyAnalyzer.AlleleSpecificLocus("1", 1, double.NaN, 0.5) }),
                "NaN logR must throw.");
            Assert.Throws<ArgumentException>(
                () => OncologyAnalyzer.SegmentAlleleSpecificAspcf(new[] { new OncologyAnalyzer.AlleleSpecificLocus("1", 1, 0.0, 1.5) }),
                "BAF outside [0,1] must throw.");
        });
    }

    #endregion

    #region FitSubclonalCopyNumber (Battenberg determine_copynumber)

    // Expected values: Battenberg determine_copynumber (Wedge-lab/battenberg R/fitcopynumber.R + R/orderEdges.R,
    // master) run in R 4.3 on the same segment (constant SNP BAFs ⇒ sd = 0 ⇒ pval = 0, maxdist = 0.01);
    // 402/402 random segments identical (Evidence §"Battenberg R cross-check").

    // M-SUB-1 — planted mixture along the nearest edge: 0.3·(2,1) + 0.7·(1,1) at ρ = 0.8, ψ = 2.5
    // (logR −0.099535673550914569, BAF 0.55357142857142849) ⇒ Battenberg state 1 = (1,1) at τ = 0.70000000000000051,
    // state 2 = (2,1) at 0.29999999999999949.
    [Test]
    public void FitSubclonalCopyNumber_PlantedEdgeMixture_MatchesBattenberg()
    {
        var seg = new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, -0.099535673550914569, 0.55357142857142849, 5);

        var fit = OncologyAnalyzer.FitSubclonalCopyNumber(new[] { seg }, 0.8, 2.5)[0];

        Assert.Multiple(() =>
        {
            Assert.That(fit.IsSubclonal, Is.True, "Battenberg: pval ≤ 0.05 ⇒ sub-clonal.");
            Assert.That((fit.PrimaryState.MajorCopyNumber, fit.PrimaryState.MinorCopyNumber), Is.EqualTo((1, 1)), "nMaj1_A:nMin1_A.");
            Assert.That(fit.PrimaryState.CellFraction, Is.EqualTo(0.70000000000000051).Within(1e-12), "frac1_A = τ.");
            Assert.That((fit.SecondaryState!.Value.MajorCopyNumber, fit.SecondaryState!.Value.MinorCopyNumber), Is.EqualTo((2, 1)), "nMaj2_A:nMin2_A.");
            Assert.That(fit.SecondaryState!.Value.CellFraction, Is.EqualTo(0.29999999999999949).Within(1e-12), "frac2_A = 1 − τ.");
        });
    }

    // M-SUB-2 — Battenberg mixes states that differ in ONE allele (an edge of the copy-number square). The pre-fix
    // code decomposed nA = 1.4, nB = 0.6 (ρ = 1, ψ = 2, logR 0, BAF 0.7) as 0.4·(2,0) + 0.6·(1,1); Battenberg's
    // nearest edge is (2,0)–(2,1) with τ = 0.14285714285714241 (the BAF-mixture fraction).
    [Test]
    public void FitSubclonalCopyNumber_BothAllelesFractional_UsesNearestEdgeAsBattenberg()
    {
        var seg = new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, 0.0, 0.7, 5);

        var fit = OncologyAnalyzer.FitSubclonalCopyNumber(new[] { seg }, 1.0, 2.0)[0];

        Assert.Multiple(() =>
        {
            Assert.That(fit.IsSubclonal, Is.True);
            Assert.That((fit.PrimaryState.MajorCopyNumber, fit.PrimaryState.MinorCopyNumber), Is.EqualTo((2, 0)), "nMaj1_A:nMin1_A.");
            Assert.That(fit.PrimaryState.CellFraction, Is.EqualTo(0.14285714285714241).Within(1e-12), "frac1_A.");
            Assert.That((fit.SecondaryState!.Value.MajorCopyNumber, fit.SecondaryState!.Value.MinorCopyNumber), Is.EqualTo((2, 1)), "nMaj2_A:nMin2_A.");
            Assert.That(fit.PrimaryState.CellFraction + fit.SecondaryState!.Value.CellFraction, Is.EqualTo(1.0).Within(1e-12),
                "frac1 + frac2 = 1.");
        });
    }

    // M-SUB-3 — clonal segments (|l − corner level| < maxdist 0.01) collapse to the corner with fraction 1.
    // Battenberg on the ASCAT forward values at ρ = 0.8, ψ = 2.5: (2,1) → C 2 1; (2,0) → C 2 0; (3,1) → C 3 1.
    [TestCase(0.22239242133644802, 0.64285714285714279, 2, 1)]
    [TestCase(-0.26303440583379378, 0.90000000000000002, 2, 0)]
    [TestCase(0.58496250072115619, 0.72222222222222232, 3, 1)]
    public void FitSubclonalCopyNumber_IntegerSegment_CollapsesToSingleClonalState(double logR, double baf, int major, int minor)
    {
        var seg = new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, logR, baf, 5);

        var fit = OncologyAnalyzer.FitSubclonalCopyNumber(new[] { seg }, 0.8, 2.5)[0];

        Assert.Multiple(() =>
        {
            Assert.That(fit.IsSubclonal, Is.False, "Battenberg maxdist rule ⇒ clonal.");
            Assert.That(fit.SecondaryState, Is.Null, "A clonal segment has no second state.");
            Assert.That((fit.PrimaryState.MajorCopyNumber, fit.PrimaryState.MinorCopyNumber), Is.EqualTo((major, minor)), "Clonal corner.");
            Assert.That(fit.PrimaryState.CellFraction, Is.EqualTo(1.0), "Clonal state is present in all tumour cells (f=1).");
        });
    }

    // C-SUB-1 — invalid arguments throw.
    [Test]
    public void FitSubclonalCopyNumber_InvalidArguments_Throw()
    {
        var seg = new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, 0.0, 0.5, 5);
        var segs = new[] { seg };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(
                () => OncologyAnalyzer.FitSubclonalCopyNumber(null!, 0.8, 2.0), "Null segments must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.FitSubclonalCopyNumber(segs, 0.0, 2.0), "purity ≤ 0 must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.FitSubclonalCopyNumber(segs, 0.8, 0.0), "ploidy ≤ 0 must throw.");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => OncologyAnalyzer.FitSubclonalCopyNumber(segs, 0.8, 2.0, gamma: 0.0), "gamma ≤ 0 must throw.");
        });
    }

    #endregion

    #region FitSubclonalCopyNumberWithSnpTest (Battenberg per-SNP t-test, B24 F39)

    // Expected values: Battenberg determine_copynumber (Wedge-lab/battenberg master R/fitcopynumber.R, sourced verbatim
    // with R/orderEdges.R; R 4.3.3, maxdist 0.01, siglevel 0.05, noperms 1000) on four genomes (ρ/ψ 0.8/2.5, 1/2,
    // 0.55/3.1, 0.35/1.9), one chromosome per segment, BAFphased = the listed SNP BAFs, BAFseg = the segment BAF,
    // LogR = the segment logR at every SNP. 180/180 segments identical (Evidence ONCO-ASCAT-001 §F39). Rows:
    // S41 noisy SNPs 0.03 off the (2,1) level (maxdist says sub-clonal, t-test p ≈ 0.56 ⇒ clonal); S42 the same offset
    // with tight SNPs (p ≈ 2e−5 ⇒ sub-clonal); S43 constant SNPs and S44 a single SNP (sd 0 / NA ⇒ pval 0 ⇒ sub-clonal);
    // S45 within maxdist (pval 1); G1S28/G3S31/G2S15 p ∈ (0.05, 0.1) ⇒ clonal; G1S40/G3S40/G4S27 p ∈ (0.02, 0.05) ⇒ sub-clonal.
    [TestCase(0.8, 2.5, -0.057038, 0.780347, new[] { 0.8864, 0.7955, 0.8486, 0.7765, 0.617 }, 0.0673250064382435, 2, 0, 1.0, -1, -1, double.NaN, TestName = "G1S28")]
    [TestCase(0.8, 2.5, 0.331898, 0.854267, new[] { 0.8622, 0.777, 0.8651, 0.9944, 0.9311, 0.8201, 0.844, 0.8369 }, 0.0354791289104208, 3, 0, 0.695568832695164, 3, 1, 0.304431167304836, TestName = "G1S40")]
    [TestCase(0.8, 2.5, 0.222392, 0.672857, new[] { 0.5529, 0.7629, 0.5729, 0.7829, 0.6929 }, 0.560904254225102, 2, 1, 1.0, -1, -1, double.NaN, TestName = "G1S41")]
    [TestCase(0.8, 2.5, 0.222392, 0.672857, new[] { 0.6689, 0.6759, 0.6739, 0.6709, 0.6749 }, 2.10208586143735e-05, 2, 0, 0.156050245445912, 2, 1, 0.843949754554088, TestName = "G1S42")]
    [TestCase(0.8, 2.5, 0.222392, 0.672857, new[] { 0.6729, 0.6729, 0.6729, 0.6729 }, 0.0, 2, 0, 0.156050245445912, 2, 1, 0.843949754554088, TestName = "G1S43")]
    [TestCase(0.8, 2.5, 0.222392, 0.672857, new[] { 0.6729 }, 0.0, 2, 0, 0.156050245445912, 2, 1, 0.843949754554088, TestName = "G1S44")]
    [TestCase(0.8, 2.5, 0.222392, 0.647857, new[] { 0.6379, 0.6579, 0.6484 }, 1.0, 2, 1, 1.0, -1, -1, double.NaN, TestName = "G1S45")]
    [TestCase(1.0, 2.0, 0.71718, 0.985, new[] { 1.0, 1.0, 0.9148, 0.9128, 1.0, 1.0, 1.0, 0.8588 }, 0.0914471543031092, 3, 0, 1.0, -1, -1, double.NaN, TestName = "G2S15")]
    [TestCase(1.0, 2.0, 0.584963, 0.696667, new[] { 0.5767, 0.7867, 0.5967, 0.8067, 0.7167 }, 0.561022844257267, 2, 1, 1.0, -1, -1, double.NaN, TestName = "G2S41")]
    [TestCase(1.0, 2.0, 0.584963, 0.696667, new[] { 0.6927, 0.6997, 0.6977, 0.6947, 0.6987 }, 2.10473683882307e-05, 2, 1, 0.703293080541847, 3, 1, 0.296706919458153, TestName = "G2S42")]
    [TestCase(0.55, 3.1, -0.787152, 0.702461, new[] { 0.7483, 0.7243, 0.7425, 0.6835, 0.7018 }, 0.0675778888092366, 1, 0, 1.0, -1, -1, double.NaN, TestName = "G3S31")]
    [TestCase(0.55, 3.1, 0.01946, 0.796306, new[] { 0.7938, 0.7866, 0.8057 }, 0.0369496702675139, 3, 0, 0.841496062614577, 3, 1, 0.158503937385423, TestName = "G3S40")]
    [TestCase(0.55, 3.1, -0.030786, 0.637843, new[] { 0.5178, 0.7278, 0.5378, 0.7478, 0.6578 }, 0.561975673235374, 2, 1, 1.0, -1, -1, double.NaN, TestName = "G3S41")]
    [TestCase(0.55, 3.1, -0.030786, 0.637843, new[] { 0.6338, 0.6408, 0.6388, 0.6358, 0.6398 }, 2.1261746475365e-05, 2, 1, 0.615939295037034, 3, 1, 0.384060704962966, TestName = "G3S42")]
    [TestCase(0.35, 1.9, 0.103462, 0.698502, new[] { 0.693, 0.7098, 0.7167, 0.7037, 0.6761 }, 0.0248852309205578, 2, 0, 0.554567052328231, 3, 0, 0.445432947671769, TestName = "G4S27")]
    [TestCase(0.35, 1.9, 0.258131, 0.604468, new[] { 0.4845, 0.6945, 0.5045, 0.7145, 0.6245 }, 0.561040508068104, 2, 1, 1.0, -1, -1, double.NaN, TestName = "G4S41")]
    [TestCase(0.35, 1.9, 0.258131, 0.604468, new[] { 0.6005, 0.6075, 0.6055, 0.6025, 0.6065 }, 2.10513202226578e-05, 2, 0, 0.333231866699312, 2, 1, 0.666768133300688, TestName = "G4S42")]
    [TestCase(0.35, 1.9, 0.258131, 0.604468, new[] { 0.6045 }, 0.0, 2, 0, 0.333231866699312, 2, 1, 0.666768133300688, TestName = "G4S44")]
    public void FitSubclonalCopyNumberWithSnpTest_MatchesBattenberg(
        double purity, double psit, double logR, double baf, double[] snps, double pval,
        int major1, int minor1, double frac1, int major2, int minor2, double frac2)
    {
        var seg = new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 1000L * snps.Length, logR, baf, snps.Length);

        var r = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(
            new[] { new OncologyAnalyzer.SubclonalSegmentSnpBafs(seg, snps) }, purity, psit)[0];

        Assert.Multiple(() =>
        {
            Assert.That(r.PValue, Is.EqualTo(pval).Within(1e-12 * pval), "Battenberg pval.");
            Assert.That(r.Fit.IsSubclonal, Is.EqualTo(major2 >= 0), "pval ≤ siglevel ⇔ sub-clonal.");
            Assert.That((r.Fit.PrimaryState.MajorCopyNumber, r.Fit.PrimaryState.MinorCopyNumber), Is.EqualTo((major1, minor1)), "nMaj1_A:nMin1_A.");
            Assert.That(r.Fit.PrimaryState.CellFraction, Is.EqualTo(frac1).Within(1e-12), "frac1_A.");
            if (major2 >= 0)
            {
                Assert.That((r.Fit.SecondaryState!.Value.MajorCopyNumber, r.Fit.SecondaryState!.Value.MinorCopyNumber),
                    Is.EqualTo((major2, minor2)), "nMaj2_A:nMin2_A.");
                Assert.That(r.Fit.SecondaryState!.Value.CellFraction, Is.EqualTo(frac2).Within(1e-12), "frac2_A.");
            }
            else
            {
                Assert.That(r.Fit.SecondaryState, Is.Null, "Clonal ⇒ one state.");
            }
        });
    }

    // M-SUBT-2 — the t-test overrides maxdist: G1S41 is sub-clonal under the maxdist-only rule but clonal in Battenberg.
    [Test]
    public void FitSubclonalCopyNumberWithSnpTest_NoisySnps_ClonalWhereMaxdistOnlyCallsSubclonal()
    {
        var seg = new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, 0.222392, 0.672857, 5);
        double[] snps = { 0.5529, 0.7629, 0.5729, 0.7829, 0.6929 };

        var maxdistOnly = OncologyAnalyzer.FitSubclonalCopyNumber(new[] { seg }, 0.8, 2.5)[0];
        var tested = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(
            new[] { new OncologyAnalyzer.SubclonalSegmentSnpBafs(seg, snps) }, 0.8, 2.5)[0];

        Assert.Multiple(() =>
        {
            Assert.That(maxdistOnly.IsSubclonal, Is.True, "|l − level| = 0.03 ≥ maxdist.");
            Assert.That(tested.Fit.IsSubclonal, Is.False, "Battenberg: p = 0.5609 > 0.05 ⇒ clonal (2,1).");
        });
    }

    // M-SUBT-3 — siglevel and maxdist are Battenberg callSubclones parameters (R: siglevel 0.1 ⇒ G1S28 sub-clonal
    // (2,0)@0.61666732876527997 + (2,1)@0.38333267123471998; maxdist 0.001 or 0 ⇒ G1S45 t-tested, p 0.46228196278297101, clonal (2,1)).
    [Test]
    public void FitSubclonalCopyNumberWithSnpTest_SiglevelAndMaxdist_MatchBattenberg()
    {
        var s28 = new OncologyAnalyzer.SubclonalSegmentSnpBafs(
            new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, -0.057038, 0.780347, 5),
            new[] { 0.8864, 0.7955, 0.8486, 0.7765, 0.617 });
        var s45 = new OncologyAnalyzer.SubclonalSegmentSnpBafs(
            new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 3000, 0.222392, 0.647857, 3),
            new[] { 0.6379, 0.6579, 0.6484 });

        var sig = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(new[] { s28 }, 0.8, 2.5, significanceLevel: 0.1)[0];
        var md1 = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(new[] { s45 }, 0.8, 2.5, maxBafDistance: 0.001)[0];
        var md0 = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(new[] { s45 }, 0.8, 2.5, maxBafDistance: 0.0)[0];

        Assert.Multiple(() =>
        {
            Assert.That(sig.Fit.IsSubclonal, Is.True, "p 0.0673 ≤ 0.1.");
            Assert.That(sig.Fit.PrimaryState.CellFraction, Is.EqualTo(0.61666732876527997).Within(1e-12), "frac1_A.");
            Assert.That(sig.Fit.SecondaryState!.Value.CellFraction, Is.EqualTo(0.38333267123471998).Within(1e-12), "frac2_A.");
            Assert.That(md1.PValue, Is.EqualTo(0.46228196278297101).Within(1e-12), "maxdist 0.001 ⇒ t-test p.");
            Assert.That(md0.PValue, Is.EqualTo(0.46228196278297101).Within(1e-12), "maxdist 0 ⇒ t-test p.");
            Assert.That(md1.Fit.IsSubclonal || md0.Fit.IsSubclonal, Is.False, "p > 0.05 ⇒ clonal.");
        });
    }

    // M-SUBT-4 — Battenberg mirrors only the segment level (l = max(BAFseg, 1 − BAFseg)), not the SNP BAFs: with
    // BAFseg 1 − 0.672857 and SNPs 1 − G1S42 the t-test compares ≈ 0.327 against the 0.643 level (R pval
    // 1.7441626943676701e-09 ⇒ sub-clonal (2,0)@0.15605024544591201).
    [Test]
    public void FitSubclonalCopyNumberWithSnpTest_SnpBafsNotMirrored_AsBattenberg()
    {
        var seg = new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, 0.222392, 1 - 0.672857, 5);
        double[] snps = { 1 - 0.6689, 1 - 0.6759, 1 - 0.6739, 1 - 0.6709, 1 - 0.6749 };

        var r = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(
            new[] { new OncologyAnalyzer.SubclonalSegmentSnpBafs(seg, snps) }, 0.8, 2.5)[0];

        Assert.Multiple(() =>
        {
            Assert.That(r.PValue, Is.EqualTo(1.7441626943676701e-09).Within(1e-12 * 1.7441626943676701e-09));
            Assert.That(r.Fit.IsSubclonal, Is.True);
            Assert.That(r.Fit.PrimaryState.CellFraction, Is.EqualTo(0.15605024544591201).Within(1e-12));
        });
    }

    // M-SUBT-5 — constant / single-SNP segments reproduce the maxdist-only overload exactly (Battenberg pval 0).
    [Test]
    public void FitSubclonalCopyNumberWithSnpTest_ConstantSnps_EqualsMaxdistOnlyOverload()
    {
        var segs = new[]
        {
            new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, -0.099535673550914569, 0.55357142857142849, 5),
            new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 6000, 9000, 0.22239242133644802, 0.64285714285714279, 5),
            new OncologyAnalyzer.AlleleSpecificSegmentSummary("2", 1000, 5000, 0.0, 0.7, 5),
        };
        var old = OncologyAnalyzer.FitSubclonalCopyNumber(segs, 0.8, 2.5);
        var tested = OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(
            segs.Select((s, i) => new OncologyAnalyzer.SubclonalSegmentSnpBafs(s, i == 2 ? new[] { 0.7 } : new[] { s.MeanBAF, s.MeanBAF })).ToArray(),
            0.8, 2.5);

        Assert.That(tested.Select(t => t.Fit), Is.EqualTo(old));
    }

    // C-SUBT-1 — invalid arguments throw.
    [Test]
    public void FitSubclonalCopyNumberWithSnpTest_InvalidArguments_Throw()
    {
        var seg = new OncologyAnalyzer.AlleleSpecificSegmentSummary("1", 1000, 5000, 0.0, 0.6, 5);
        var ok = new[] { new OncologyAnalyzer.SubclonalSegmentSnpBafs(seg, new[] { 0.6, 0.62 }) };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(null!, 0.8, 2.0));
            Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(
                new[] { new OncologyAnalyzer.SubclonalSegmentSnpBafs(seg, null!) }, 0.8, 2.0));
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(
                new[] { new OncologyAnalyzer.SubclonalSegmentSnpBafs(seg, new[] { 0.6, 1.2 }) }, 0.8, 2.0), "SNP BAF > 1.");
            Assert.Throws<ArgumentException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(
                new[] { new OncologyAnalyzer.SubclonalSegmentSnpBafs(seg, new[] { double.NaN }) }, 0.8, 2.0), "NaN SNP BAF.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(ok, 0.0, 2.0));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(ok, 0.8, 2.0, significanceLevel: 1.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.FitSubclonalCopyNumberWithSnpTest(ok, 0.8, 2.0, maxBafDistance: -0.1));
        });
    }

    #endregion
}
