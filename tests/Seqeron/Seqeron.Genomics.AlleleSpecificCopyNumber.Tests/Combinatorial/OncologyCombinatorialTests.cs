// SPDX-License-Identifier: GPL-3.0-only
// Tests of the Port of ASCAT 3.2.0 and Battenberg 3.0.0; part of Seqeron.Genomics.AlleleSpecificCopyNumber.Tests, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Oncology area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids (3³ = 27) use the exhaustive <c>[Combinatorial]</c> product,
/// a strict superset of pairwise.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Oncology")]
public class OncologyCombinatorialTests
{

    /// <summary>Replicates each integer segment into adjacent loci carrying its planted (logR, BAF), one chromosome
    /// per segment: with 5 loci (&lt; ASCAT kmin = 6) ascat.aspcf — run by SegmentAlleleSpecific (B24 F35) — emits
    /// exactly one segment per chromosome (noise-free data has MAD 0 ⇒ no in-chromosome breakpoint).</summary>
    private static List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> AscatLoci(
        IReadOnlyList<(int NA, int NB)> segments, double rho, double psi, int lociPerSegment = 5)
    {
        var loci = new List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus>();
        long pos = 1000;
        int chromosome = 0;
        foreach (var (nA, nB) in segments)
        {
            (double r, double b) = AscatForward(nA, nB, rho, psi);
            string chrom = (++chromosome).ToString(System.Globalization.CultureInfo.InvariantCulture);
            for (int i = 0; i < lociPerSegment; i++)
            {
                loci.Add(new AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus(chrom, pos, r, b));
                pos += 1000;
            }
        }
        return loci;
    }

    /// <summary>
    /// Pairwise grid over (purity × ploidy): from data synthesised at the planted (ρ,ψ) the joint
    /// ASCAT grid fit recovers both parameters within one grid step and the exact integer (major,
    /// minor) copy numbers per segment, with goodness-of-fit ≈ 100%.
    /// </summary>
    [Test, Combinatorial]
    public void AscatFit_RecoversPlantedPurityPloidyAndIntegerCopyNumbers(
        [Values(0.50, 0.70, 0.90)] double purity,
        [Values(2.0, 2.5, 3.0)] double ploidy)
    {
        var segments = AscatSegmentsForPloidy(ploidy);
        var loci = AscatLoci(segments, purity, ploidy);

        var summaries = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(
            loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);

        AlleleSpecificCopyNumberAnalyzer.PurityPloidyFit fit = AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy(
            summaries, purityMin: 0.05, purityMax: 1.0, purityStep: 0.01,
            ploidyMin: 1.5, ploidyMax: 5.0, ploidyStep: 0.05, gamma: AscatGamma);

        fit.Purity.Should().BeApproximately(purity, 0.01 + 1e-9, "the ρ axis is recovered within one grid step");
        fit.Ploidy.Should().BeApproximately(ploidy, 0.05 + 1e-9, "the ψ axis is recovered within one grid step");
        fit.GoodnessOfFit.Should().BeApproximately(100.0, 1e-6, "integer CN are recovered exactly at the truth ⇒ GoF 100%");

        fit.Segments.Should().HaveCount(segments.Length, "one emitted run per planted (adjacent-distinct) segment");
        for (int i = 0; i < segments.Length; i++)
        {
            int major = Math.Max(segments[i].NA, segments[i].NB);
            int minor = Math.Min(segments[i].NA, segments[i].NB);
            fit.Segments[i].MajorCopyNumber.Should().Be(major, $"segment {i} recovers the planted major CN {major}");
            fit.Segments[i].MinorCopyNumber.Should().Be(minor, $"segment {i} recovers the planted minor CN {minor}");
        }
    }

    /// <summary>
    /// Interaction witness: the two axes genuinely move the answer. Holding ploidy fixed, data
    /// synthesised at purity 0.9 recovers a higher ρ than data synthesised at 0.5; holding purity
    /// fixed, data synthesised at ψ=3.0 recovers a higher ψ than data at ψ=2.0. A purity/ploidy
    /// inversion that ignored either axis could not satisfy both.
    /// </summary>
    [Test]
    public void AscatFit_PurityAndPloidyAxes_EachDriveTheFit()
    {
        AlleleSpecificCopyNumberAnalyzer.PurityPloidyFit FitAt(double rho, double psi)
        {
            var loci = AscatLoci(AscatSegmentsForPloidy(psi), rho, psi);
            var summaries = AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2, minLociPerSegment: 1);
            return AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy(summaries, gamma: AscatGamma);
        }

        FitAt(0.90, 2.0).Purity.Should().BeGreaterThan(FitAt(0.50, 2.0).Purity, "the purity axis drives the recovered ρ");
        FitAt(0.70, 3.0).Ploidy.Should().BeGreaterThan(FitAt(0.70, 2.0).Ploidy, "the ploidy axis drives the recovered ψ");
    }


    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ONCO-ASCAT-001 — Allele-specific purity/ploidy fit (Oncology)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 235.
    // Spec: tests/TestSpecs/ONCO-ASCAT-001.md (AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy). ADVANCED §10.
    //
    // Sources: Van Loo et al. (2010), ASCAT, PNAS 107:16910; VanLoo-lab/ascat runAscat.R (γ=1).
    //
    // Model: per-locus logR/BAF are generated from integer allele copy numbers (nA,nB) at an
    // aberrant cell fraction ρ (purity) and average sample ploidy ψ by the ASCAT forward model
    //   r = log2[(ρ·n + 2(1−ρ)) / (ρ·ψ + 2(1−ρ))],   b = (ρ·nB + (1−ρ)) / (ρ·n + 2(1−ρ)),  n = nA+nB
    // (the exact algebraic inverse of the cited nA/nB equations). FitPurityPloidy inverts this
    // with a (ρ × ψ) grid search, selecting the (ρ,ψ) whose inferred copy numbers are closest to
    // nonnegative integers (goodness-of-fit → 100%).
    //
    // Dimensions: purity(3) × ploidy(3). Full grid 3×3 = 9 cells (exhaustive [Combinatorial]).
    //
    // The combinatorial point: ρ and ψ INTERACT in both forward terms — the BAF scales with ρ,
    // the logR baseline with ρ·ψ — so a naive read of one axis is confounded by the other. Each
    // cell is synthesised from a segment set whose length-weighted mean total CN equals the
    // planted ψ, and the joint fit must recover BOTH planted parameters (within one grid step)
    // AND the exact integer (major,minor) CN per segment with GoF ≈ 100%, i.e. it disentangles
    // the two coupled axes in every cell.
    // ═══════════════════════════════════════════════════════════════════════

    private const double AscatGamma = 1.0;


    /// <summary>ASCAT forward model (γ=1): integer (nA,nB) at (ρ,ψ) ⇒ the locus (logR, BAF).</summary>
    private static (double LogR, double Baf) AscatForward(int nA, int nB, double rho, double psi)
    {
        int n = nA + nB;
        double denom = rho * n + 2.0 * (1.0 - rho);
        double d = rho * psi + 2.0 * (1.0 - rho);
        return (Math.Log2(denom / d), (rho * nB + (1.0 - rho)) / denom);
    }

    // Integer-CN segment sets whose mean total CN equals each planted ploidy. No two ADJACENT
    // segments are equal (each entry is also its own chromosome, so segmentation emits one run per entry) and each set carries an LOH
    // segment (minor = 0, BAF far from 0.5) which breaks the 2n-vs-4n sunrise symmetry.
    private static (int NA, int NB)[] AscatSegmentsForPloidy(double ploidy) => ploidy switch
    {
        2.0 => new[] { (1, 1), (2, 0), (1, 1) },         // mean (2+2+2)/3 = 2.0
        2.5 => new[] { (1, 1), (2, 1), (2, 0), (2, 1) }, // mean (2+3+2+3)/4 = 2.5
        _ => new[] { (2, 1), (3, 0), (2, 1) },           // ploidy 3.0: mean (3+3+3)/3 = 3.0
    };
}
