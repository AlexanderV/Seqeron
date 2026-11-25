// SPDX-License-Identifier: GPL-3.0-only
// Tests of the Port of ASCAT 3.2.0 and Battenberg 3.0.0; part of Seqeron.Genomics.AlleleSpecificCopyNumber.Tests, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the production project).
namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic tests for the Oncology area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-IMMUNE-001 — immune/stromal infiltration (ESTIMATE/ssGSEA) (Oncology).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 86.
///
/// API under test (ImmuneAnalyzer.EstimateInfiltration):
///   Computes ssGSEA enrichment of an immune (and stromal) signature gene set against a gene
///   expression profile. ssGSEA is RANK-based: it ranks genes by expression and integrates a
///   running enrichment over the ranked list, so only the ordering of expression matters.
///
/// Relations (derived from rank-based enrichment, NOT from output):
///   • INV  (scaling ⇒ same infiltration): multiplying every expression value by a positive
///          constant preserves the ranking, hence the scores are unchanged.
///   • INV/SYM (gene order independent): the score depends only on the gene→value map, not on
///          the order in which the profile was assembled.
///   • MON  (higher marker expression ⇒ higher score): raising the immune-signature genes'
///          expression lifts them up the ranking, raising the immune enrichment score.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class OncologyMetamorphicTests
{

    /// <summary>Replicates each (chrom, nA, nB) segment into <paramref name="lociPerSegment"/> adjacent loci 1 kb apart.</summary>
    private static List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> AscatSynthesiseLoci(
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
    [Description("INV (B24 F35): SegmentAlleleSpecific runs ascat.aspcf at ASCAT's default penalty 70, whose MAD-standardised penalised cost is translation-invariant, so a constant logR shift preserves the breakpoints and shifts each segment mean by c.")]
    public void Ascat_ConstantLogRShift_PreservesSegmentAlleleSpecificBreakpoints()
    {
        // Noisy two-level track (ASCAT's getMad is 0 on noise-free data, which places no breakpoint at all).
        // ascat.aspcf (R, penalty 70) gives the same breakpoints for shifts 0.3137, −0.8123, 1.9071 (Evidence).
        List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> loci =
            Seqeron.Genomics.Tests.Unit.Oncology.OncologyAnalyzer_AscatDerivation_Tests.AspcfStepTrack();

        IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary> baseline =
            AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.3);
        baseline.Count.Should().Be(2,
            because: "the noisy two-level track gives exactly one ASPCF breakpoint — the non-vacuity guard");

        foreach (double c in new[] { 0.3137, -0.8123, 1.9071 })
        {
            List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> shifted =
                Seqeron.Genomics.Tests.Unit.Oncology.OncologyAnalyzer_AscatDerivation_Tests.AspcfStepTrack(c);

            IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary> segs =
                AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(shifted, logRChangeThreshold: 0.3);

            segs.Count.Should().Be(baseline.Count,
                because: $"the standardised joint cost is translation-invariant, so a logR shift of {c} keeps the breakpoint count");
            for (int i = 0; i < baseline.Count; i++)
            {
                segs[i].LocusCount.Should().Be(baseline[i].LocusCount,
                    because: $"breakpoint {i} sits at the same locus after a constant logR shift of {c}");
                segs[i].Start.Should().Be(baseline[i].Start);
                segs[i].End.Should().Be(baseline[i].End);
                segs[i].MeanLogR.Should().BeApproximately(baseline[i].MeanLogR + c, 1e-9,
                    because: $"each segment's mean logR shifts by exactly the applied constant {c}");
                segs[i].MeanBAF.Should().BeApproximately(baseline[i].MeanBAF, 1e-12,
                    because: "a logR shift does not touch the BAF channel");
            }
        }
    }

    [Test]
    [Description("INV: ASCAT's ASPCF standardises each track by its MAD (translation-invariant) and minimises SSE/sd² + penalty (translation-invariant), so a constant logR shift preserves the breakpoints and shifts each segment mean by c.")]
    public void Ascat_ConstantLogRShift_PreservesAspcfBreakpoints()
    {
        // Noisy two-level track (ASCAT's getMad is 0 on noise-free data, which places no breakpoint at all).
        // ascat.aspcf (R, penalty 70) gives the same breakpoints for shifts 0.3137, −0.8123, 1.9071 (Evidence).
        List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> loci =
            Seqeron.Genomics.Tests.Unit.Oncology.OncologyAnalyzer_AscatDerivation_Tests.AspcfStepTrack();

        IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary> baseline =
            AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAspcf(loci, penalty: 70.0);
        baseline.Count.Should().Be(2,
            because: "the noisy two-level track gives exactly one breakpoint — the non-vacuity guard");

        // Shifts chosen so no shifted logR is exactly 0 (getMad drops exact zeros as imputed values).
        foreach (double c in new[] { 0.3137, -0.8123, 1.9071 })
        {
            List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> shifted =
                Seqeron.Genomics.Tests.Unit.Oncology.OncologyAnalyzer_AscatDerivation_Tests.AspcfStepTrack(c);

            IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary> segs =
                AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecificAspcf(shifted, penalty: 70.0);

            segs.Count.Should().Be(baseline.Count,
                because: $"the standardised joint cost is translation-invariant, so the optimum is unchanged by a logR shift of {c}");
            for (int i = 0; i < baseline.Count; i++)
            {
                segs[i].LocusCount.Should().Be(baseline[i].LocusCount,
                    because: $"the optimal breakpoint stays put under a constant logR shift of {c}");
                segs[i].MeanLogR.Should().BeApproximately(baseline[i].MeanLogR + c, 1e-9,
                    because: $"each segment's mean logR shifts by exactly the applied constant {c}");
                segs[i].MeanBAF.Should().BeApproximately(baseline[i].MeanBAF, 1e-12,
                    because: "the BAF channel is untouched by a logR shift");
            }
        }
    }

    [Test]
    [Description("INV: swapping the arbitrary A/B allele label maps BAF b→1−b, which sends nA↔nB (total CN invariant) and the segmenter's BAF-mirroring folds identically, so every segment's integer copy-number fit is unchanged.")]
    public void Ascat_AlleleSwap_PreservesTotalCopyNumber()
    {
        const double rho = 0.80, psi = 2.2;
        List<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus> loci = AscatSynthesiseLoci(AscatPlantedGenome, rho, psi);
        var swapped = loci
            .Select(l => new AlleleSpecificCopyNumberAnalyzer.AlleleSpecificLocus(l.Chromosome, l.Position, l.LogR, 1.0 - l.BAF))
            .ToList();

        IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary> baseSegs =
            AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(loci, logRChangeThreshold: 0.2);
        IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary> swapSegs =
            AlleleSpecificCopyNumberAnalyzer.SegmentAlleleSpecific(swapped, logRChangeThreshold: 0.2);

        // Non-vacuity: the planted LOH/gain segments are allele-imbalanced (mirrored BAF > 0.5),
        // so the swap b → 1−b is a genuine change of the raw input, not the identity.
        baseSegs.Any(s => s.MeanBAF > 0.5 + 1e-6).Should().BeTrue(
            because: "the planted LOH/gain segments are allele-imbalanced, so swapping A/B is non-trivial");

        AlleleSpecificCopyNumberAnalyzer.PurityPloidyFit baseFit = AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy(baseSegs, gamma: 1.0);
        AlleleSpecificCopyNumberAnalyzer.PurityPloidyFit swapFit = AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy(swapSegs, gamma: 1.0);

        swapFit.Segments.Count.Should().Be(baseFit.Segments.Count,
            because: "the allele swap does not change segmentation, so the fit emits the same number of segments");
        for (int i = 0; i < baseFit.Segments.Count; i++)
        {
            int baseTotal = baseFit.Segments[i].MajorCopyNumber + baseFit.Segments[i].MinorCopyNumber;
            int swapTotal = swapFit.Segments[i].MajorCopyNumber + swapFit.Segments[i].MinorCopyNumber;

            swapTotal.Should().Be(baseTotal,
                because: $"total copy number nA+nB is independent of the BAF, so swapping A/B preserves segment {i}'s total CN");
            // Stronger: BAF-mirroring folds b and 1−b identically, so even the major/minor split is preserved.
            swapFit.Segments[i].MajorCopyNumber.Should().Be(baseFit.Segments[i].MajorCopyNumber,
                because: $"the mirrored BAF is invariant under b→1−b, so segment {i}'s major CN is unchanged");
            swapFit.Segments[i].MinorCopyNumber.Should().Be(baseFit.Segments[i].MinorCopyNumber,
                because: $"the mirrored BAF is invariant under b→1−b, so segment {i}'s minor CN is unchanged");
        }
    }

    /// <summary>A planted diploid genome with balanced, copy-neutral-LOH and gain segments (some allele-imbalanced).</summary>
    // One chromosome per segment (5 loci < ASCAT kmin = 6), so ascat.aspcf — run by SegmentAlleleSpecific (B24 F35) —
    // emits one segment per planted segment.
    private static readonly (string Chrom, int NA, int NB)[] AscatPlantedGenome =
    {
        ("1", 1, 1), // balanced, b = 0.5
        ("2", 2, 0), // copy-neutral LOH, b ≠ 0.5
        ("3", 1, 1),
        ("4", 2, 1), // gain, b ≠ 0.5
        ("5", 1, 1),
    };


    // ═══════════════════════════════════════════════════════════════════
    //  ONCO-ASCAT-001 — allele-specific copy-number derivation (ASCAT)
    // ═══════════════════════════════════════════════════════════════════
    //
    // Theory (Van Loo et al. 2010 PNAS 107:16910; VanLoo-lab/ascat ascat.runAscat.R;
    //   docs/algorithms/Oncology/Allele_Specific_Copy_Number_Derivation.md):
    //   ASCAT maps a segment's total-coverage log-ratio (logR r) and B-allele frequency
    //   (BAF b) to allele-specific copy numbers (nA, nB). With platform parameter γ = 1:
    //       S      = 2^r · ((1−ρ)·2 + ρ·ψ)              (total-signal scale)
    //       nA     = (ρ−1 − (b−1)·S) / ρ
    //       nB     = (ρ−1 +  b   ·S) / ρ
    //       nA+nB  = (2(ρ−1) + S) / ρ                    ← does NOT depend on the BAF b
    //   Two metamorphic relations (checklist row 235) follow from these equations and from
    //   the way the upstream segmenter places breakpoints:
    //
    //   • INV (constant logR shift preserves breakpoints): ASCAT's ASPCF (run by both
    //     SegmentAlleleSpecificAspcf and, since B24 F35, SegmentAlleleSpecific at penalty 70) minimises the
    //     MAD-standardised within-segment logR SSE plus a fixed per-segment penalty. Adding a constant c to every
    //     locus's logR leaves every within-segment SSE and the MAD unchanged: the breakpoint set is invariant and
    //     each segment's mean logR simply shifts by c. The BAF channel is untouched.
    //
    //   • INV (A/B allele swap preserves total CN): "allele A" vs "allele B" is an arbitrary
    //     label. Swapping it maps the raw BAF b → 1−b, which the equations above send to
    //     nA ↔ nB exactly (so the total nA+nB is invariant), and which the segmenter's
    //     BAF-mirroring (ascat.aspcf max(b, 1−b)) folds to the *same* value. The whole integer
    //     copy-number fit — and in particular every segment's total copy number — is therefore
    //     invariant to the allele labelling.
    //
    // API under test: OncologyAnalyzer.SegmentAlleleSpecific, .SegmentAlleleSpecificAspcf,
    //   .FitPurityPloidy (AlleleSpecificLocus / AlleleSpecificSegmentSummary / PurityPloidyFit).


    /// <summary>ASCAT forward model (γ = 1): integer (nA, nB) at (ρ, ψ) → the (logR, BAF) a locus would show.</summary>
    private static (double LogR, double Baf) AscatForward(int nA, int nB, double rho, double psi)
    {
        int n = nA + nB;
        double denom = rho * n + 2.0 * (1.0 - rho);
        double d = rho * psi + 2.0 * (1.0 - rho);
        return (System.Math.Log2(denom / d), (rho * nB + (1.0 - rho)) / denom);
    }



    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ONCO-HETERO-001 — intratumour heterogeneity (MATH) (Oncology).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 116.
    //
    // API under test (OncologyAnalyzer.AnalyzeHeterogeneity):
    //   The MATH score is 100·(scaled MAD)/median of the VAFs (Mroz & Rocco 2013) — a
    //   coefficient-of-dispersion-like ratio.
    //
    // Relations (derived from the MAD/median ratio, NOT from output):
    //   • INV  (scaling all VAFs equally preserves MATH): MAD and median scale together, so the
    //          ratio (MATH) is unchanged.
    //   • MON  (wider VAF spread ⇒ ≥ heterogeneity): at a fixed median, a wider VAF spread raises
    //          the MAD and hence the MATH score.
    // ───────────────────────────────────────────────────────────────────────────

    private static double Math(double[] vafs) =>
        OncologyAnalyzer.AnalyzeHeterogeneity(vafs, vafs, 1).MathScore;
}
