// SPDX-License-Identifier: GPL-3.0-only
// Tests of the Port of ASCAT 3.2.0 and Battenberg 3.0.0; part of Seqeron.Genomics.AlleleSpecificCopyNumber.Tests, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the production project).
using FsCheck;
using FsCheck.Fluent;

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

    private static List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> SynthesiseAscatLoci(
        IReadOnlyList<(string Chrom, int NA, int NB)> segments, double rho, double psi, int lociPerSegment = 5)
    {
        var loci = new List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus>();
        long pos = 1000;
        foreach (var (chrom, nA, nB) in segments)
        {
            (double r, double b) = AscatForward(nA, nB, rho, psi);
            for (int i = 0; i < lociPerSegment; i++)
            {
                loci.Add(new AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus(chrom, pos, r, b));
                pos += 1000;
            }
        }

        return loci;
    }

    [Test]
    public void Ascat_Identity_NoLociIsNoSegments()
    {
        // ID: the empty locus set is the neutral input; its segmentation is the empty list.
        var segments = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(
            System.Array.Empty<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus>(), logRChangeThreshold: 0.2);
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
                .Select(i => new AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus("1", 1000 + i * 1000, t.logR[i], t.baf[i]))
                .ToList();
            var a = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.3, bafChangeThreshold: 0.1);
            var b = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.3, bafChangeThreshold: 0.1);
            return a.SequenceEqual(b).Label($"segmentation non-deterministic for {n} loci");
        });
    }

    [Test]
    public void Ascat_Idempotent_GridOptimumIsDeterministic()
    {
        // IDEMP: the ASCAT grid search is a pure function — re-running it over the same segments and
        // the same (ρ × ψ) grid recovers the identical optimum (purity, ploidy, GoF, integer segments).
        var loci = SynthesiseAscatLoci(AscatPlantedSegments, rho: 0.80, psi: 2.2);
        var summaries = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);

        AlleleSpecificCopyNumberAnalyzer.PurityPloidyFit fit1 = AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy(summaries, gamma: AscatGamma);
        AlleleSpecificCopyNumberAnalyzer.PurityPloidyFit fit2 = AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy(summaries, gamma: AscatGamma);

        fit2.Purity.Should().Be(fit1.Purity);
        fit2.Ploidy.Should().Be(fit1.Ploidy);
        fit2.GoodnessOfFit.Should().Be(fit1.GoodnessOfFit);
        fit2.Segments.Should().Equal(fit1.Segments);
    }


    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-ASCAT-001 — Upstream allele-specific derivation (Oncology), row 235.
    // ID — no loci → no segments (the empty input maps to the empty segmentation).
    // IDEMP — the segmentation and the ASCAT grid optimum are pure, deterministic
    //         functions: f(x) = f(x) for the same input and grid.
    //   — AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific / FitPurityPloidy.
    //   — Van Loo et al. (2010), PNAS 107:16910 (grid over ploidy × aberrant fraction);
    //     ascat.runAscat.R nA/nB equations + sunrise goodness of fit.
    // ═══════════════════════════════════════════════════════════════════════

    private const double AscatGamma = 1.0;
 // sequencing data (ASCAT README).

    // Exact algebraic inverse of the two cited ASCAT nA/nB equations (gamma = 1):
    //   denom = rho·n + 2(1−rho); D = rho·psi + 2(1−rho); r = log2(denom/D); b = (rho·nB + (1−rho))/denom.
    private static (double LogR, double Baf) AscatForward(int nA, int nB, double rho, double psi)
    {
        int n = nA + nB;
        double denom = rho * n + 2.0 * (1.0 - rho);
        double d = rho * psi + 2.0 * (1.0 - rho);
        return (System.Math.Log2(denom / d), (rho * nB + (1.0 - rho)) / denom);
    }

    // One chromosome per planted segment (5 loci < ASCAT kmin = 6): SegmentAlleleSpecific runs ascat.aspcf (B24 F35),
    // which emits one segment per such chromosome (noise-free data has MAD 0 ⇒ no in-chromosome breakpoint).
    private static readonly (string Chrom, int NA, int NB)[] AscatPlantedSegments =
    {
        ("1", 1, 1), // balanced diploid, b = 0.5
        ("2", 2, 0), // copy-neutral LOH, b ≈ 0.1
        ("3", 1, 1),
        ("4", 2, 1), // gain
        ("5", 1, 1),
    };
}
