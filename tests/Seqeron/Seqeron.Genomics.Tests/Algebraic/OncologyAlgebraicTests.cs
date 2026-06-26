using System.Collections.Generic;
using System.Linq;
using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Algebraic;

using VariantObservation = OncologyAnalyzer.VariantObservation;

/// <summary>
/// Algebraic-law tests for the Oncology area (VAF, TMB, signature exposures,
/// tumour fraction, cancer-cell fraction).
///
/// Algebraic testing pins the identity values on null/zero input and the
/// homogeneity (scale-invariance / linearity) laws of these quantitative
/// oncology estimators.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 88, 92, 97, 111, 115.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("Oncology")]
public class OncologyAlgebraicTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-VAF-001 — Variant allele fraction (Oncology), row 88.
    // ID — alt=0 → VAF=0.  HOMO — VAF(k·alt, k·total) = VAF(alt, total).
    //   — OncologyAnalyzer.CalculateVAF.
    // ═══════════════════════════════════════════════════════════════════════

    [FsCheck.NUnit.Property]
    public Property Vaf_Identity_ZeroAltIsZero()
    {
        return Prop.ForAll(Gen.Choose(1, 10_000).ToArbitrary(), total =>
            (OncologyAnalyzer.CalculateVAF(0, total) == 0.0).Label($"VAF(0,{total})"));
    }

    [FsCheck.NUnit.Property]
    public Property Vaf_Homogeneous_ScaleInvariant()
    {
        var gen = (from total in Gen.Choose(1, 1000)
                   from alt in Gen.Choose(0, total)
                   from k in Gen.Choose(1, 100)
                   select (alt, total, k)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            double baseVaf = OncologyAnalyzer.CalculateVAF(t.alt, t.total);
            double scaled = OncologyAnalyzer.CalculateVAF(t.k * t.alt, t.k * t.total);
            return (System.Math.Abs(baseVaf - scaled) < 1e-12).Label($"VAF scale variant: {baseVaf} vs {scaled}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-TMB-001 — Tumour mutational burden (Oncology), row 92.
    // ID — zero mutations → TMB=0.  HOMO — TMB scales inversely with panel-Mb.
    //   — OncologyAnalyzer.CalculateTMB.
    // ═══════════════════════════════════════════════════════════════════════

    [FsCheck.NUnit.Property]
    public Property Tmb_Identity_ZeroMutationsIsZero()
    {
        return Prop.ForAll(Gen.Choose(1, 100).Select(x => x / 10.0).ToArbitrary(), mb =>
            (OncologyAnalyzer.CalculateTMB(0, mb) == 0.0).Label($"TMB(0,{mb})"));
    }

    [FsCheck.NUnit.Property]
    public Property Tmb_Homogeneous_InverseInPanelSize()
    {
        var gen = (from count in Gen.Choose(0, 10_000)
                   from mbTimes10 in Gen.Choose(1, 500)
                   from k in Gen.Choose(2, 50)
                   select (count, mb: mbTimes10 / 10.0, k)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            double baseTmb = OncologyAnalyzer.CalculateTMB(t.count, t.mb);
            double scaled = OncologyAnalyzer.CalculateTMB(t.count, t.k * t.mb);
            // TMB(count, k·Mb) = TMB(count, Mb) / k.
            return (System.Math.Abs(scaled * t.k - baseTmb) < 1e-9).Label($"TMB inverse: {scaled}*{t.k} vs {baseTmb}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-SIG-002 — Signature exposure refit (Oncology), row 97.
    // ID — zero catalogue → zero exposures.  HOMO — exposure(k·catalogue)=k·exposure.
    //   — OncologyAnalyzer.FitSignatures (non-negative least squares).
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly IReadOnlyList<IReadOnlyList<double>> Signatures = new IReadOnlyList<double>[]
    {
        new double[] { 0.7, 0.2, 0.1 },
        new double[] { 0.1, 0.3, 0.6 },
    };

    [Test]
    public void SignatureExposure_Identity_ZeroCatalogIsZeroExposures()
    {
        var fit = OncologyAnalyzer.FitSignatures(new double[] { 0, 0, 0 }, Signatures);
        fit.Exposures.Should().OnlyContain(e => e == 0.0);
    }

    [FsCheck.NUnit.Property]
    public Property SignatureExposure_Homogeneous_LinearInCatalog()
    {
        var gen = (from a in Gen.Choose(0, 500)
                   from b in Gen.Choose(0, 500)
                   from c in Gen.Choose(0, 500)
                   from k in Gen.Choose(2, 20)
                   select (cat: new double[] { a, b, c }, k)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            var baseExp = OncologyAnalyzer.FitSignatures(t.cat, Signatures).Exposures;
            var scaledCat = t.cat.Select(x => x * t.k).ToArray();
            var scaledExp = OncologyAnalyzer.FitSignatures(scaledCat, Signatures).Exposures;
            bool ok = baseExp.Zip(scaledExp, (b2, s) => System.Math.Abs(s - t.k * b2) < 1e-6).All(x => x);
            return ok.Label($"exposures not linear: base=[{string.Join(",", baseExp)}] k={t.k}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-CTDNA-001 — Tumour (ctDNA) fraction (Oncology), row 111.
    // ID — no tumour reads → fraction 0.  HOMO — scale-invariant to total depth.
    //   — OncologyAnalyzer.CalculateTumorFraction.
    // ═══════════════════════════════════════════════════════════════════════

    private static VariantObservation Variant(int alt, int total) =>
        new("chr1", 100, "A", "T", alt, total, 0, 100);

    [Test]
    public void TumorFraction_Identity_NoTumorReadsIsZero()
    {
        var variants = new[] { Variant(0, 100), Variant(0, 200), Variant(0, 50) };
        OncologyAnalyzer.CalculateTumorFraction(variants).Should().Be(0.0);
    }

    [FsCheck.NUnit.Property]
    public Property TumorFraction_Homogeneous_ScaleInvariantToDepth()
    {
        // VAF kept <= 0.5 (clonal-het constraint); scale read depth by k.
        var gen = (from n in Gen.Choose(1, 5)
                   from altHalf in Gen.Choose(0, 25).ArrayOf(n)   // alt <= total/2 with total=100
                   from k in Gen.Choose(2, 20)
                   select (alts: altHalf, k)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            var baseVariants = t.alts.Select(a => Variant(a, 100)).ToList();
            var scaledVariants = t.alts.Select(a => Variant(a * t.k, 100 * t.k)).ToList();
            double baseFrac = OncologyAnalyzer.CalculateTumorFraction(baseVariants);
            double scaledFrac = OncologyAnalyzer.CalculateTumorFraction(scaledVariants);
            return (System.Math.Abs(baseFrac - scaledFrac) < 1e-12).Label($"{baseFrac} vs {scaledFrac}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-CCF-001 — Cancer cell fraction (Oncology), row 115.
    // ID — VAF=0 → CCF=0.  HOMO — CCF linear in VAF at fixed CN/purity/multiplicity.
    //   — OncologyAnalyzer.EstimateCcf.
    // ═══════════════════════════════════════════════════════════════════════

    [Test]
    public void Ccf_Identity_ZeroVafIsZero()
    {
        var ccf = OncologyAnalyzer.EstimateCcf(vaf: 0.0, purity: 0.8, tumorCopyNumber: 2, multiplicity: 1);
        ccf.Ccf.Should().Be(0.0);
        ccf.RawCcf.Should().Be(0.0);
    }

    [FsCheck.NUnit.Property]
    public Property Ccf_Homogeneous_LinearInVaf()
    {
        // RawCcf (uncapped) is exactly proportional to VAF at fixed CN/purity/multiplicity.
        var gen = (from vafCenti in Gen.Choose(1, 25)        // vaf in (0, 0.25]
                   from k in Gen.Choose(2, 4)                 // keep k·vaf <= 1
                   select (vaf: vafCenti / 100.0, k)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            double baseRaw = OncologyAnalyzer.EstimateCcf(t.vaf, 0.8, 2, 1).RawCcf;
            double scaledRaw = OncologyAnalyzer.EstimateCcf(t.vaf * t.k, 0.8, 2, 1).RawCcf;
            return (System.Math.Abs(scaledRaw - t.k * baseRaw) < 1e-9).Label($"{scaledRaw} vs {t.k}*{baseRaw}");
        });
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-ASCAT-001 — Upstream allele-specific derivation (Oncology), row 235.
    // ID — no loci → no segments (the empty input maps to the empty segmentation).
    // IDEMP — the segmentation and the ASCAT grid optimum are pure, deterministic
    //         functions: f(x) = f(x) for the same input and grid.
    //   — OncologyAnalyzer.SegmentAlleleSpecific / FitPurityPloidy.
    //   — Van Loo et al. (2010), PNAS 107:16910 (grid over ploidy × aberrant fraction);
    //     ascat.runAscat.R nA/nB equations + sunrise goodness of fit.
    // ═══════════════════════════════════════════════════════════════════════

    private const double AscatGamma = 1.0; // sequencing data (ASCAT README).

    // Exact algebraic inverse of the two cited ASCAT nA/nB equations (gamma = 1):
    //   denom = rho·n + 2(1−rho); D = rho·psi + 2(1−rho); r = log2(denom/D); b = (rho·nB + (1−rho))/denom.
    private static (double LogR, double Baf) AscatForward(int nA, int nB, double rho, double psi)
    {
        int n = nA + nB;
        double denom = rho * n + 2.0 * (1.0 - rho);
        double d = rho * psi + 2.0 * (1.0 - rho);
        return (System.Math.Log2(denom / d), (rho * nB + (1.0 - rho)) / denom);
    }

    private static List<OncologyAnalyzer.AlleleSpecificLocus> SynthesiseAscatLoci(
        IReadOnlyList<(string Chrom, int NA, int NB)> segments, double rho, double psi, int lociPerSegment = 5)
    {
        var loci = new List<OncologyAnalyzer.AlleleSpecificLocus>();
        long pos = 1000;
        foreach (var (chrom, nA, nB) in segments)
        {
            (double r, double b) = AscatForward(nA, nB, rho, psi);
            for (int i = 0; i < lociPerSegment; i++)
            {
                loci.Add(new OncologyAnalyzer.AlleleSpecificLocus(chrom, pos, r, b));
                pos += 1000;
            }
        }

        return loci;
    }

    private static readonly (string Chrom, int NA, int NB)[] AscatPlantedSegments =
    {
        ("1", 1, 1), // balanced diploid, b = 0.5
        ("1", 2, 0), // copy-neutral LOH, b ≈ 0.1
        ("1", 1, 1),
        ("1", 2, 1), // gain
        ("1", 1, 1),
    };

    [Test]
    public void Ascat_Identity_NoLociIsNoSegments()
    {
        // ID: the empty locus set is the neutral input; its segmentation is the empty list.
        var segments = OncologyAnalyzer.SegmentAlleleSpecific(
            System.Array.Empty<OncologyAnalyzer.AlleleSpecificLocus>(), logRChangeThreshold: 0.2);
        segments.Should().BeEmpty();
    }

    [FsCheck.NUnit.Property]
    public Property Ascat_Idempotent_SegmentationIsDeterministic()
    {
        // IDEMP: SegmentAlleleSpecific is a pure function — the same loci always yield the same segments.
        var gen = (from n in Gen.Choose(1, 30)
                   from logR in Gen.Choose(-300, 300).Select(x => x / 100.0).ArrayOf(n)
                   from baf in Gen.Choose(0, 100).Select(x => x / 100.0).ArrayOf(n)
                   select (logR, baf)).ToArbitrary();
        return Prop.ForAll(gen, t =>
        {
            int n = System.Math.Min(t.logR.Length, t.baf.Length);
            var loci = Enumerable.Range(0, n)
                .Select(i => new OncologyAnalyzer.AlleleSpecificLocus("1", 1000 + i * 1000, t.logR[i], t.baf[i]))
                .ToList();
            var a = OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.3, bafChangeThreshold: 0.1);
            var b = OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.3, bafChangeThreshold: 0.1);
            return a.SequenceEqual(b).Label($"segmentation non-deterministic for {n} loci");
        });
    }

    [Test]
    public void Ascat_Idempotent_GridOptimumIsDeterministic()
    {
        // IDEMP: the ASCAT grid search is a pure function — re-running it over the same segments and
        // the same (ρ × ψ) grid recovers the identical optimum (purity, ploidy, GoF, integer segments).
        var loci = SynthesiseAscatLoci(AscatPlantedSegments, rho: 0.80, psi: 2.2);
        var summaries = OncologyAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);

        OncologyAnalyzer.PurityPloidyFit fit1 = OncologyAnalyzer.FitPurityPloidy(summaries, gamma: AscatGamma);
        OncologyAnalyzer.PurityPloidyFit fit2 = OncologyAnalyzer.FitPurityPloidy(summaries, gamma: AscatGamma);

        fit2.Purity.Should().Be(fit1.Purity);
        fit2.Ploidy.Should().Be(fit1.Ploidy);
        fit2.GoodnessOfFit.Should().Be(fit1.GoodnessOfFit);
        fit2.Segments.Should().Equal(fit1.Segments);
    }
}
