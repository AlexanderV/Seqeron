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
}
