// EPIGEN-DMR-001 — Differentially Methylated Regions (tiling window + Fisher's exact test)
// Evidence: docs/Evidence/EPIGEN-DMR-001-Evidence.md
// TestSpec: tests/TestSpecs/EPIGEN-DMR-001.md
// Source: Akalin A, et al. (2012). methylKit. Genome Biology 13:R87. doi:10.1186/gb-2012-13-10-r87
//         methylKit tileMethylCounts/getMethylDiff/calculateDiffMeth manuals (al2na/methylKit).
//         Fisher's exact test (hypergeometric 2x2; Fisher 1922, 1935).

namespace Seqeron.Genomics.Tests.Unit.Annotation;

using MethylationType = EpigeneticsAnalyzer.MethylationType;
using MethylationSite = EpigeneticsAnalyzer.MethylationSite;
using GeneAnnotation = EpigeneticsAnalyzer.GeneAnnotation;
using DMR = EpigeneticsAnalyzer.DifferentiallyMethylatedRegion;

[TestFixture]
public class EpigeneticsAnalyzer_DMR_Tests
{
    private static MethylationSite Site(int pos, double level, int coverage = 20) =>
        new(pos, MethylationType.CpG, "CG", level, coverage);

    private static IEnumerable<MethylationSite> Profile(double level, int count, int spacing = 1, int coverage = 20) =>
        Enumerable.Range(0, count).Select(i => Site(i * spacing, level, coverage));

    #region FindDMRs

    // M1 — Hyper: treatment fully methylated, control unmethylated → meanDiff +1, "Hypermethylated"
    // (Akalin 2012: hyper = higher methylation than control; getMethylDiff meth.diff > difference).
    [Test]
    public void FindDMRs_HyperMethylatedWindow_ReportedAsHypermethylated()
    {
        var control = Profile(level: 0.0, count: 3);
        var treatment = Profile(level: 1.0, count: 3);

        var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(dmrs, Has.Count.EqualTo(1), "one fully-differential window must yield one DMR");
            Assert.That(dmrs[0].MeanDifference, Is.EqualTo(1.0).Within(1e-10),
                "mean(level2 - level1) = 1.0 - 0.0 = 1.0 across all sites");
            Assert.That(dmrs[0].Annotation, Is.EqualTo("Hypermethylated"),
                "positive difference = treatment higher than control = hypermethylated (Akalin 2012)");
            Assert.That(dmrs[0].CpGCount, Is.EqualTo(3), "all three covered sites are in the region");
            // Pooled 2x2 table is [[0,60],[60,0]] (complete separation).
            // scipy.stats.fisher_exact([[0,60],[60,0]]) -> two-sided p = 2.070073888186964e-35.
            Assert.That(dmrs[0].PValue, Is.EqualTo(2.070073888186964e-35).Within(1e-44),
                "two-sided Fisher's exact p of the fully-separated pooled table (scipy reference)");
        });
    }

    // M2 — Hypo: control methylated, treatment unmethylated → meanDiff -1, "Hypomethylated"
    [Test]
    public void FindDMRs_HypoMethylatedWindow_ReportedAsHypomethylated()
    {
        var control = Profile(level: 1.0, count: 3);
        var treatment = Profile(level: 0.0, count: 3);

        var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(dmrs, Has.Count.EqualTo(1), "one fully-differential window must yield one DMR");
            Assert.That(dmrs[0].MeanDifference, Is.EqualTo(-1.0).Within(1e-10),
                "mean(level2 - level1) = 0.0 - 1.0 = -1.0 across all sites");
            Assert.That(dmrs[0].Annotation, Is.EqualTo("Hypomethylated"),
                "negative difference = treatment lower than control = hypomethylated (Akalin 2012)");
        });
    }

    // M3 — Strict cutoff: |meanDiff| == minDifference must NOT be reported
    // (getMethylDiff uses meth.diff > difference, strict).
    [Test]
    public void FindDMRs_DifferenceEqualToCutoff_NotReported()
    {
        // control 0.0, treatment 0.25 → meanDiff exactly 0.25 == minDifference.
        var control = Profile(level: 0.0, count: 3);
        var treatment = Profile(level: 0.25, count: 3);

        var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment, minDifference: 0.25).ToList();

        Assert.That(dmrs, Is.Empty,
            "meanDiff equal to the cutoff is excluded because the cutoff is strict (> difference)");
    }

    // M4 — Above cutoff is reported
    [Test]
    public void FindDMRs_DifferenceAboveCutoff_Reported()
    {
        var control = Profile(level: 0.0, count: 3);
        var treatment = Profile(level: 0.30, count: 3);

        var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment, minDifference: 0.25).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(dmrs, Has.Count.EqualTo(1), "meanDiff 0.30 > cutoff 0.25 → reported");
            Assert.That(dmrs[0].MeanDifference, Is.EqualTo(0.30).Within(1e-10),
                "mean(0.30 - 0.0) = 0.30");
        });
    }

    // M5 — Tiling: site clusters farther apart than windowSize are split into separate windows
    // (methylKit win.size tiling).
    [Test]
    public void FindDMRs_SitesBeyondWindowSize_SplitIntoSeparateWindows()
    {
        // Cluster A at 0,1,2 ; cluster B at 1000,1001,1002 ; windowSize 100.
        var control = new[]
        {
            Site(0, 0.0), Site(1, 0.0), Site(2, 0.0),
            Site(1000, 0.0), Site(1001, 0.0), Site(1002, 0.0),
        };
        var treatment = new[]
        {
            Site(0, 1.0), Site(1, 1.0), Site(2, 1.0),
            Site(1000, 1.0), Site(1001, 1.0), Site(1002, 1.0),
        };

        var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment, windowSize: 100).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(dmrs, Has.Count.EqualTo(2),
                "two clusters >windowSize apart form two separate windows");
            Assert.That(dmrs[0].Start, Is.EqualTo(0), "first window starts at the first cluster");
            Assert.That(dmrs[1].Start, Is.EqualTo(1000), "second window starts at the second cluster");
            Assert.That(dmrs[0].End, Is.EqualTo(2), "first window ends at its last site");
            Assert.That(dmrs[1].End, Is.EqualTo(1002), "second window ends at its last site");
        });
    }

    // M7 — Empty input → no DMRs
    [Test]
    public void FindDMRs_EmptyInput_ReturnsNoRegions()
    {
        var dmrs = EpigeneticsAnalyzer.FindDMRs(
            Enumerable.Empty<MethylationSite>(),
            Enumerable.Empty<MethylationSite>()).ToList();

        Assert.That(dmrs, Is.Empty, "no covered positions means no tiles and no DMRs");
    }

    // M8 — PValue within [0,1] AND equal to the externally-computed two-sided Fisher value.
    // control 0.1 cov 30 -> numC=3,numT=27 ; treatment 0.9 cov 30 -> numC=27,numT=3 ; 4 sites pooled
    // -> 2x2 table [[12,108],[108,12]]. scipy.stats.fisher_exact([[12,108],[108,12]]) -> p = 2.475428262210228e-39.
    [Test]
    public void FindDMRs_ReportedRegion_PValueWithinUnitInterval()
    {
        var control = Profile(level: 0.1, count: 4, coverage: 30);
        var treatment = Profile(level: 0.9, count: 4, coverage: 30);

        var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment).ToList();

        Assert.That(dmrs, Is.Not.Empty, "strongly differential window must be reported");
        Assert.Multiple(() =>
        {
            Assert.That(dmrs[0].PValue, Is.InRange(0.0, 1.0),
                "a Fisher's exact p-value is a probability in [0,1]");
            // Exact two-sided value (independent reference: scipy.stats.fisher_exact); lock it, do not
            // settle for a bounds check that any wrong implementation would also pass.
            Assert.That(dmrs[0].PValue, Is.EqualTo(2.475428262210228e-39).Within(1e-48),
                "two-sided Fisher's exact p of pooled table [[12,108],[108,12]] (scipy reference)");
        });
    }

    // S1 — Window with fewer covered sites than minCpGCount is not reported
    [Test]
    public void FindDMRs_FewerSitesThanMinCpGCount_NotReported()
    {
        var control = new[] { Site(0, 0.0), Site(1, 0.0) };   // only 2 sites
        var treatment = new[] { Site(0, 1.0), Site(1, 1.0) };

        var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment, minCpGCount: 3).ToList();

        Assert.That(dmrs, Is.Empty,
            "a region needs >= minCpGCount adjacent covered sites (Akalin 2012)");
    }

    // S3 — A position present only in sample1 is compared against an implicit level 0 in sample2
    [Test]
    public void FindDMRs_PositionMissingInSecondSample_TreatedAsZero()
    {
        // control has 3 methylated sites; treatment has none of those positions.
        var control = new[] { Site(0, 1.0), Site(1, 1.0), Site(2, 1.0) };
        var treatment = Enumerable.Empty<MethylationSite>();

        var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(dmrs, Has.Count.EqualTo(1),
                "missing treatment positions are compared against level 0");
            Assert.That(dmrs[0].MeanDifference, Is.EqualTo(-1.0).Within(1e-10),
                "mean(0 - 1.0) = -1.0 (treatment absent = level 0)");
        });
    }

    // S4 — Null sample throws ArgumentNullException
    [Test]
    public void FindDMRs_NullSample_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => EpigeneticsAnalyzer.FindDMRs(null!, Enumerable.Empty<MethylationSite>()).ToList(),
                NUnit.Framework.Throws.ArgumentNullException, "null sample1 must be rejected eagerly");
            Assert.That(() => EpigeneticsAnalyzer.FindDMRs(Enumerable.Empty<MethylationSite>(), null!).ToList(),
                NUnit.Framework.Throws.ArgumentNullException, "null sample2 must be rejected eagerly");
        });
    }

    // C1 — Determinism and ordering by Start ascending
    [Test]
    public void FindDMRs_RepeatedCalls_DeterministicAndOrderedByStart()
    {
        var control = new[]
        {
            Site(1000, 0.0), Site(1001, 0.0), Site(1002, 0.0),
            Site(0, 0.0), Site(1, 0.0), Site(2, 0.0),
        };
        var treatment = new[]
        {
            Site(1000, 1.0), Site(1001, 1.0), Site(1002, 1.0),
            Site(0, 1.0), Site(1, 1.0), Site(2, 1.0),
        };

        var first = EpigeneticsAnalyzer.FindDMRs(control, treatment, windowSize: 100).ToList();
        var second = EpigeneticsAnalyzer.FindDMRs(control, treatment, windowSize: 100).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(first.Select(d => d.Start), Is.EqualTo(new[] { 0, 1000 }),
                "DMRs are ordered by Start ascending regardless of input order");
            Assert.That(second.Select(d => (d.Start, d.End, d.MeanDifference)),
                Is.EqualTo(first.Select(d => (d.Start, d.End, d.MeanDifference))),
                "repeated calls produce identical output (deterministic)");
        });
    }

    #endregion

    #region FisherExactProbability

    // M6 — Single-table hypergeometric probability matches the published worked example.
    // Wikipedia Fisher's exact test: a=1,b=9,c=11,d=3 (n=24) → p ≈ 0.001346076.
    [Test]
    public void FisherExactProbability_WikipediaWorkedExample_MatchesPublishedValue()
    {
        double p = EpigeneticsAnalyzer.FisherExactProbability(1, 9, 11, 3);

        Assert.That(p, Is.EqualTo(0.001346076).Within(1e-7),
            "hypergeometric probability of the published 2x2 table (Fisher's exact test worked example)");
    }

    // FisherExactProbability: empty table (n == 0) is the degenerate everything-zero case → p = 1.0.
    [Test]
    public void FisherExactProbability_ZeroTotal_ReturnsOne()
    {
        Assert.That(EpigeneticsAnalyzer.FisherExactProbability(0, 0, 0, 0), Is.EqualTo(1.0),
            "an all-zero contingency table has no information → single-table probability defined as 1.0");
    }

    // FisherExactProbability: negative cell counts are invalid and rejected.
    [Test]
    public void FisherExactProbability_NegativeCell_Throws()
    {
        Assert.That(() => EpigeneticsAnalyzer.FisherExactProbability(-1, 9, 11, 3),
            NUnit.Framework.Throws.TypeOf<ArgumentOutOfRangeException>(),
            "contingency-table cells must be non-negative");
    }

    // FisherExactProbability: a balanced symmetric table — the single-table probability matches
    // the hypergeometric value computed independently (scipy hypergeom / C(.)·C(.)/C(.)).
    // [[5,5],[5,5]] -> C(10,5)*C(10,5)/C(20,10) = 252*252/184756 = 0.34371820130334063.
    [Test]
    public void FisherExactProbability_SymmetricTable_MatchesHypergeometric()
    {
        Assert.That(EpigeneticsAnalyzer.FisherExactProbability(5, 5, 5, 5),
            Is.EqualTo(0.34371820130334063).Within(1e-12),
            "single-table hypergeometric probability of [[5,5],[5,5]] (independent reference)");
    }

    // S2 — Degenerate margin (a zero row/column total) → only one feasible table → not differential.
    // Verified through FindDMRs: a window with zero coverage in one group is not reported.
    [Test]
    public void FisherExactProbability_DegenerateMargin_ReturnsOne()
    {
        // Sample2 has all-zero coverage within the window → numC2 = numT2 = 0 (zero row total).
        var control = new[] { Site(0, 0.0, 20), Site(1, 0.0, 20), Site(2, 0.0, 20) };
        var treatment = new[] { Site(0, 1.0, 0), Site(1, 1.0, 0), Site(2, 1.0, 0) };

        var dmrs = EpigeneticsAnalyzer.FindDMRs(control, treatment).ToList();

        // meanDiff = 1.0 (levels), so it passes the difference cutoff, but the 2x2 table is
        // degenerate (sample2 contributes no reads) → Fisher p = 1.0.
        Assert.Multiple(() =>
        {
            Assert.That(dmrs, Has.Count.EqualTo(1), "level difference still passes the cutoff");
            Assert.That(dmrs[0].PValue, Is.EqualTo(1.0).Within(1e-10),
                "a degenerate fixed margin admits only the observed table → p = 1.0");
        });
    }

    #endregion

    #region AnnotateDMRs

    // C2 — Delegation: a DMR overlapping a gene annotation is relabelled with the feature name.
    [Test]
    public void AnnotateDMRs_OverlappingAnnotation_LabelsRegion()
    {
        var dmrs = new[]
        {
            new DMR(Start: 100, End: 200, MeanDifference: 0.5, PValue: 0.001, CpGCount: 4, Annotation: "Hypermethylated"),
            new DMR(Start: 5000, End: 5100, MeanDifference: -0.5, PValue: 0.001, CpGCount: 4, Annotation: "Hypomethylated"),
        };
        var annotations = new[] { new GeneAnnotation("TP53", 150, 400) };

        var annotated = EpigeneticsAnalyzer.AnnotateDMRs(dmrs, annotations).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(annotated[0].Annotation, Is.EqualTo("TP53"),
                "DMR [100,200] overlaps feature [150,400) → labelled with the feature");
            Assert.That(annotated[1].Annotation, Is.EqualTo("Hypomethylated"),
                "non-overlapping DMR keeps its methylation annotation");
        });
    }

    // AnnotateDMRs — null inputs are rejected eagerly (input-validation contract).
    [Test]
    public void AnnotateDMRs_NullInput_Throws()
    {
        var dmrs = new[]
        {
            new DMR(Start: 100, End: 200, MeanDifference: 0.5, PValue: 0.001, CpGCount: 4, Annotation: "Hypermethylated"),
        };
        var annotations = new[] { new GeneAnnotation("TP53", 150, 400) };

        Assert.Multiple(() =>
        {
            Assert.That(() => EpigeneticsAnalyzer.AnnotateDMRs(null!, annotations),
                NUnit.Framework.Throws.ArgumentNullException, "null dmrs must be rejected");
            Assert.That(() => EpigeneticsAnalyzer.AnnotateDMRs(dmrs, null!),
                NUnit.Framework.Throws.ArgumentNullException, "null annotations must be rejected");
        });
    }

    #endregion
}
