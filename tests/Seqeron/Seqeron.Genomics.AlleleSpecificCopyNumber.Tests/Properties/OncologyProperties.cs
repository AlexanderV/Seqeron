// SPDX-License-Identifier: GPL-3.0-only
// Tests of the Port of ASCAT 3.2.0 and Battenberg 3.0.0; part of Seqeron.Genomics.AlleleSpecificCopyNumber.Tests, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the production project).
using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the Oncology algorithm group.
///
/// This file is the single home for the entire Oncology block of checklist 01
/// (rows #86–#120). Each test unit lives in its own <c>#region</c> so siblings can
/// be appended without disturbing the others. Oracles are derived INDEPENDENTLY
/// from the cited theory/doc (never routed through production constants), so a
/// self-consistent-but-wrong production constant is still caught.
///
/// Test Units: ONCO-IMMUNE-001, ONCO-SOMATIC-001, ONCO-VAF-001, ONCO-DRIVER-001,
/// ONCO-ARTIFACT-001, ONCO-ANNOT-001, ONCO-TMB-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Oncology")]
public class OncologyProperties
{



    // -------------------------------------------------------------------------
    // Theory (Van Loo et al. 2010, PNAS 107:16910 — ASCAT; ascat.runAscat.R):
    //   • FitPurityPloidy builds the ASCAT distance matrix, takes strict 7×7 local minima passing the runASCAT
    //     filter cascade, and maps each segment to integer allele-specific copy numbers (seg_raw rounding) ≥ 0.
    //   • The recovered purity ρ lies on the purity grid (values > 1 reported as 1); ψ lies on the ploidy grid.
    // The grid search is a deterministic pure function of its inputs.
    // -------------------------------------------------------------------------

    /// <summary>A non-empty set of allele-specific segment summaries with folded BAF ∈ [0.5,1] and finite logR.</summary>
    private static Arbitrary<IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary>> AscatSegmentsArbitrary() =>
        (from n in Gen.Choose(1, 5)
         from segs in (from logRcenti in Gen.Choose(-200, 200)
                       from bafCenti in Gen.Choose(0, 50)
                       from len in Gen.Choose(1, 100)
                       select new AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary(
                           "1", 0, len, logRcenti / 100.0, 0.5 + bafCenti / 100.0, len)).ArrayOf(n)
         select (IReadOnlyList<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary>)segs).ToArbitrary();

    /// <summary>
    /// R + P: whenever ASCAT finds an optimum (TryFitPurityPloidy = true) it returns a purity ρ ∈ (0,1] (grid points
    /// above 1 reported as 1), a ψ on the grid seq(min_ploidy − 0.5, max_ploidy + 0.5, 0.05), GoF ∈ (80, 100] and
    /// integer allele-specific segments with major ≥ minor ≥ 0 (ascat.runAscat.R); otherwise (rho = NA) the throwing
    /// FitPurityPloidy reports InvalidOperationException.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitPurityPloidy_PurityPloidyAndCopyNumbers_AreInValidRanges()
    {
        return Prop.ForAll(AscatSegmentsArbitrary(), segs =>
        {
            if (!AlleleSpecificCopyNumberAnalyzer.TryFitPurityPloidy(segs, out var fit))
            {
                bool throws;
                try { AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy(segs); throws = false; }
                catch (InvalidOperationException) { throws = true; }
                return throws.Label("no ASCAT optimum ⇒ FitPurityPloidy throws InvalidOperationException");
            }

            bool purityOk = fit.Purity > 0.0 && fit.Purity <= 1.0;
            bool psiOk = fit.Psi >= 1.0 - 1e-9 && fit.Psi <= 6.0 + 1e-9;
            bool gofOk = fit.GoodnessOfFit > 80.0 && fit.GoodnessOfFit <= 100.0 + 1e-9;
            bool cnOk = fit.Segments.All(s => s.MinorCopyNumber >= 0 && s.MajorCopyNumber >= s.MinorCopyNumber);
            return (purityOk && psiOk && gofOk && cnOk)
                .Label($"ρ={fit.Purity}, ψ={fit.Psi}, GoF={fit.GoodnessOfFit}, segments={fit.Segments.Count}");
        });
    }

    /// <summary>
    /// D (determinism): the ASCAT search is a pure function — identical segments yield an identical outcome
    /// (found / not found, purity, ψ, ploidy, goodness of fit, and every implied copy-number segment).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property FitPurityPloidy_IsDeterministic()
    {
        return Prop.ForAll(AscatSegmentsArbitrary(), segs =>
        {
            bool okA = AlleleSpecificCopyNumberAnalyzer.TryFitPurityPloidy(segs, out var a);
            bool okB = AlleleSpecificCopyNumberAnalyzer.TryFitPurityPloidy(segs, out var b);
            bool same = okA == okB && (!okA || (a.Purity == b.Purity && a.Ploidy == b.Ploidy && a.Psi == b.Psi
                        && a.GoodnessOfFit == b.GoodnessOfFit
                        && a.Segments.SequenceEqual(b.Segments)));
            return same.Label("FitPurityPloidy must be deterministic for identical input");
        });
    }

    /// <summary>
    /// Guards: a null or empty segment list is rejected (purity/ploidy is undefined with no data).
    /// </summary>
    [Test]
    [Category("Property")]
    public void FitPurityPloidy_NullOrEmpty_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy(null!));
        Assert.Throws<ArgumentException>(() =>
            AlleleSpecificCopyNumberAnalyzer.FitPurityPloidy(Array.Empty<AlleleSpecificCopyNumberAnalyzer.AlleleSpecificSegmentSummary>()));
    }


    /// <summary>
    /// Generates a variant together with random thresholds τ_t, τ_n ∈ [0, 1], so the status oracle is
    /// exercised against arbitrary cutoffs rather than only the defaults.
    /// </summary>
    private static Arbitrary<(OncologyAnalyzer.VariantObservation variant, double tauT, double tauN)>
        VariantWithThresholdsArbitrary() =>
        (from variant in VariantArbitrary().Generator
         from tt in Gen.Choose(0, 1000).Select(x => x / 1000.0)
         from tn in Gen.Choose(0, 1000).Select(x => x / 1000.0)
         select (variant, tt, tn))
        .ToArbitrary();

    /// <summary>An arbitrary candidate list paired with valid custom STAR-Fusion thresholds.</summary>
    private static Arbitrary<(OncologyAnalyzer.FusionCandidate[] candidates, OncologyAnalyzer.FusionDetectionThresholds thresholds)>
        FusionProblemArbitrary() =>
        (from candidates in FusionCandidateListArbitrary().Generator
         from minJunc in Gen.Choose(0, 3)
         from minSum in Gen.Choose(1, 6)
         from minSpan in Gen.Choose(1, 8)
         select (candidates, new OncologyAnalyzer.FusionDetectionThresholds(minJunc, minSum, minSpan))).ToArbitrary();

    private static Arbitrary<(OncologyAnalyzer.CopyNumberArmSegment[] segments, OncologyAnalyzer.FocalAmplificationThresholds thresholds)>
        FocalProblemArbitrary() =>
        (from segments in ArmSegmentGen().ArrayOf()
         from tampMilli in Gen.Choose(0, 500) // GISTIC2 -ta range [0, Inf] (gp_gistic2_from_seg.m)
         from cutoffMilli in Gen.Choose(500, 990)
         select (segments, new OncologyAnalyzer.FocalAmplificationThresholds(tampMilli / 1000.0, cutoffMilli / 1000.0)))
        .ToArbitrary();

    private static Arbitrary<(OncologyAnalyzer.CopyNumberArmSegment[] segments, double[] cutoffs, double ploidy)>
        DeletionProblemArbitrary() =>
        (from segments in DeletionSegmentGen().ArrayOf()
         from cutoffs in AscendingThresholdsArbitrary().Generator
         from ploidy in PloidyGen()
         select (segments, cutoffs, ploidy)).ToArbitrary(); // total ≥ alt ⇒ IMAF ≤ 1

    private static Arbitrary<(OncologyAnalyzer.TumorMarker[] panel, int posThresh, int minReads, int n)>
        MrdProblemArbitrary() =>
        (from count in Gen.Choose(1, 8)
         from panel in TumorMarkerGen().ArrayOf(count)
         from posThresh in Gen.Choose(1, 5)
         from minReads in Gen.Choose(1, 5)
         from n in Gen.Choose(0, 10_000)
         select (panel, posThresh, minReads, n)).ToArbitrary();


    /// <summary>
    /// Generates a contract-valid <see cref="OncologyAnalyzer.VariantObservation"/> (0 ≤ alt ≤ total in
    /// both tumor and normal), with tumor totals drawn ≥ 1 so f_t can land on either side of τ_t.
    /// </summary>
    private static Arbitrary<OncologyAnalyzer.VariantObservation> VariantArbitrary() =>
        (from tTotal in Gen.Choose(1, 500)
         from tAlt in Gen.Choose(0, tTotal)
         from nTotal in Gen.Choose(0, 500)
         from nAlt in Gen.Choose(0, Math.Max(0, nTotal))
         from pos in Gen.Choose(1, 1_000_000)
         select new OncologyAnalyzer.VariantObservation(
             "chr1", pos, "A", "T", tAlt, tTotal, nAlt, nTotal))
        .ToArbitrary();

    /// <summary>A list of 0..8 arbitrary fusion candidates.</summary>
    private static Arbitrary<OncologyAnalyzer.FusionCandidate[]> FusionCandidateListArbitrary() =>
        (from n in Gen.Choose(0, 8)
         from list in FusionCandidateGen().ArrayOf(n)
         select list).ToArbitrary();

    private static Gen<double> PloidyGen() => Gen.Choose(1000, 4000).Select(v => v / 1000.0); // 1.0 .. 4.0

    /// <summary>Four strictly ascending threshold cutoffs, built from positive increments.</summary>
    private static Arbitrary<double[]> AscendingThresholdsArbitrary() =>
        (from t0 in Gen.Choose(-4000, -1000)
         from d1 in Gen.Choose(1, 2000)
         from d2 in Gen.Choose(1, 2000)
         from d3 in Gen.Choose(1, 2000)
         select new[] { t0 / 1000.0, (t0 + d1) / 1000.0, (t0 + d1 + d2) / 1000.0, (t0 + d1 + d2 + d3) / 1000.0 })
        .ToArbitrary();

    private static Gen<OncologyAnalyzer.CopyNumberArmSegment> ArmSegmentGen() =>
        from arm in Gen.Elements(FocalArmPool)
        from armLen in Gen.Choose(100_000, 2_000_000)
        from segLen in Gen.Choose(1, 2_200_000)
        from log2Milli in Gen.Choose(-2000, 2000)
        select new OncologyAnalyzer.CopyNumberArmSegment(arm, 0, segLen, armLen, log2Milli / 1000.0);

    private static Gen<OncologyAnalyzer.CopyNumberArmSegment> DeletionSegmentGen() =>
        from arm in Gen.Elements(DeletionArmPool)
        from armLen in Gen.Choose(100_000, 2_000_000)
        from segLen in Gen.Choose(1, 2_000_000)
        from log2Milli in Gen.Choose(-2500, 1000) // biased negative so CN 0 occurs
        select new OncologyAnalyzer.CopyNumberArmSegment(arm, 0, segLen, armLen, log2Milli / 1000.0);


    /// <summary>An arbitrary fusion candidate: pooled gene symbols, non-negative read counts, optional frame info.</summary>
    private static Gen<OncologyAnalyzer.FusionCandidate> FusionCandidateGen() =>
        from g5 in Gen.Elements(FusionGenePool)
        from g3 in Gen.Elements(FusionGenePool)
        from split5 in Gen.Choose(0, 20)
        from split3 in Gen.Choose(0, 20)
        from disc in Gen.Choose(0, 20)
        from coding in Gen.Choose(-1, 400)
        from phase in Gen.Choose(-1, 2)
        select new OncologyAnalyzer.FusionCandidate(g5, g3, split5, split3, disc, coding, phase);




    // -------------------------------------------------------------------------
    // Theory (STAR-Fusion Haas 2017/2019; Arriba Uhrig 2021; exon-phase rule):
    //   • TotalSupport = split_reads1 + split_reads2 + discordant_mates.            (Arriba; INV-2)
    //   • JunctionReads = split_reads1 + split_reads2.                              (Arriba)
    //   • A candidate is CALLED iff gene5p ≠ gene3p (distinct genes; INV-1) AND     (STAR-Fusion;
    //       - junctionReads ≥ minJunctionReads AND total ≥ minSumFrags, OR           INV-3)
    //       - junctionReads == 0 AND discordant ≥ minSpanningFragsOnly.
    //   • Calls are ordered by descending TotalSupport (abundance of support).      (INV-4)
    //   • In-frame ⇔ (fivePrimeCodingBases − threePrimeStartPhase) mod 3 == 0.      (INV-5)
    //
    // The detection predicate, support sums and frame rule are reconstructed here
    // from the cited STAR-Fusion / Arriba defaults and the modulo-3 exon-phase rule
    // — NOT routed through DetectFusions / ComputeTotalSupport / IsInFrame — so a
    // self-consistent-but-wrong production threshold or sum is still caught.
    // -------------------------------------------------------------------------

    /// <summary>Gene-symbol pool for fusion candidates; small so same-gene (INV-1) pairs recur.</summary>
    private static readonly string[] FusionGenePool =
    {
        "EML4", "ALK", "ROS1", "CD74", "RET", "KIF5B", "TMPRSS2", "ERG", "BRAF", "NTRK1",
    };



    // -------------------------------------------------------------------------
    // Theory (Mermel 2011 GISTIC2.0; GISTIC2 t_amp / broad_len_cutoff; NCBI Gene):
    //   • Focal ⟺ segment length / arm length < broad_len_cutoff (default 0.98, strict).   (INV-1)
    //   • Amplified ⟺ log2 > t_amp (default 0.1, strict).                                    (INV-2)
    //   • DetectFocalAmplifications keeps the input subset (amplified ∧ focal), in order.    (INV-3)
    //   • Oncogene reported iff a focal amp lies on its arm: ERBB2 17q, MYC 8q, EGFR 7p,     (INV-4)
    //     CCND1 11q, MDM2 12q, CDK4 12q.
    //
    // The focal/amplitude predicate and the oncogene-arm map are reconstructed here from
    // GISTIC2 / NCBI Gene — NOT routed through production — so a wrong boundary sense or
    // gene-arm assignment is still caught.
    // -------------------------------------------------------------------------

    private static readonly string[] FocalArmPool = { "17q", "8q", "7p", "11q", "12q", "3p", "Xq", "5q" };



    // -------------------------------------------------------------------------
    // Theory (Cheng 2017; cBioPortal −2; CNVkit integer CN):
    //   • A segment is a homozygous (deep) deletion iff its hard-threshold integer CN is 0.   (INV-1)
    //   • A single-copy loss (integer CN 1) is heterozygous, never homozygous.                 (INV-2)
    //   • DetectHomozygousDeletions keeps the input subset (CN 0), in input order.             (INV-3)
    //   • A stricter (more negative) deletion cutoff yields ≤ deletions (subset).              (M)
    //   • Tumour suppressor reported iff a homozygous deletion lies on its arm: TP53 17p,      (INV-4)
    //     RB1 13q, CDKN2A 9p, PTEN 10q, BRCA1 17q, BRCA2 13q.
    //
    // The CN-0 predicate is reconstructed via the independent CNVkit threshold oracle
    // (OracleCallCopyNumber, shared with ONCO-CNA-001) and the arm→gene map is restated
    // from NCBI Gene — NOT routed through production.
    // -------------------------------------------------------------------------

    private static readonly string[] DeletionArmPool = { "17p", "13q", "9p", "10q", "17q", "8q", "3p", "Xp" };



    // -------------------------------------------------------------------------
    // Theory (Reinert 2019 / Signatera ≥2-variant rule; Wan 2020 INVAR IMAF; Avanzini 2020):
    //   • A marker is detected iff PlasmaAltReads ≥ minSupportingReads.
    //   • MRD-positive ⟺ detected variants ≥ positivity threshold (default 2).   (R, M)
    //   • IMAF = Σ alt / Σ total (depth-weighted plasma VAF) ∈ [0,1].             (R sensitivity)
    //   • Panel detection p = 1 − e^(−n·IMAF·m) ∈ [0,1], m = tracked markers.
    //   • Longitudinal: per-timepoint DetectMRD; first positive index = earliest positive.
    //
    // The detection count, status, IMAF and Poisson p are reconstructed independently
    // — NOT routed through production — so a wrong threshold or weighting is caught.
    // -------------------------------------------------------------------------

    private static Gen<OncologyAnalyzer.TumorMarker> TumorMarkerGen() =>
        from alt in Gen.Choose(0, 50)
        from extra in Gen.Choose(0, 100)
        select new OncologyAnalyzer.TumorMarker("1", 1, "A", "T", alt, alt + extra);
}
