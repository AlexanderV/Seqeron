using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Infrastructure;
using static Seqeron.Genomics.Annotation.TranscriptomeAnalyzer;
using TranscriptIsoform = Seqeron.Genomics.Annotation.TranscriptomeAnalyzer.TranscriptIsoform;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// TRANS-DIFF/EXPR/SPLICE-001 mutation killers: exact-value / boundary tests for the closed-form
/// transcriptome helpers the canonical fixtures only smoke-tested — log2 fold-change regulation and
/// the two-criterion DE significance rule, over-representation analysis (hypergeometric z), GSEA
/// running-sum enrichment, Pearson co-expression, alternative-splicing classification (Wang 2008
/// five-class taxonomy: SE/RI/A5SS/A3SS/MXE), differential PSI, isoform dominance/switching, RNA-seq
/// QC rates and the PCA top-variable-gene projection.
/// </summary>
[TestFixture]
public class TranscriptomeAnalyzerMutationTests
{
    private const double Tol = 1e-9;

    private static IReadOnlyList<double> L(params double[] v) => v;

    #region AnalyzeDifferentialExpression — regulation + two-criterion significance

    private static DifferentialExpression OneGene(
        IReadOnlyList<double> g1, IReadOnlyList<double> g2, double fcThr = 1.0, double pThr = 0.05)
        => AnalyzeDifferentialExpression(
            new[] { ("G", g1, g2) }, foldChangeThreshold: fcThr, pValueThreshold: pThr).Single();

    [Test]
    public void DiffExpr_RegulationUpDownUnchanged()
    {
        Assert.That(OneGene(L(1, 1, 1), L(5, 5, 5)).Regulation, Is.EqualTo("Upregulated"));   // mean2 > mean1
        Assert.That(OneGene(L(5, 5, 5), L(1, 1, 1)).Regulation, Is.EqualTo("Downregulated")); // mean2 < mean1
        // mean1 == 0 ⇒ log2FC forced to 0 ⇒ Unchanged (guard mean1 > 0).
        var unchanged = OneGene(L(0, 0, 0), L(5, 5, 5));
        Assert.That(unchanged.Log2FoldChange, Is.EqualTo(0).Within(Tol));
        Assert.That(unchanged.Regulation, Is.EqualTo("Unchanged"));
    }

    [Test]
    public void DiffExpr_Log2FoldChangeUsesHundredthPseudocount()
    {
        // Single replicate per group: log2((mean2 + 0.01)/(mean1 + 0.01)).
        var expected = Math.Log2((3 + 0.01) / (1 + 0.01));
        Assert.That(OneGene(L(1), L(3)).Log2FoldChange, Is.EqualTo(expected).Within(Tol));
    }

    [Test]
    public void DiffExpr_SignificantRequiresBothFoldChangeAndAdjustedP()
    {
        // Zero within-group variance + separated means ⇒ p = 0 ⇒ adjusted p = 0 (single gene).
        // |log2FC| large ⇒ both criteria met ⇒ significant.
        Assert.That(OneGene(L(1, 1, 1), L(5, 5, 5)).IsSignificant, Is.True);

        // Fold change below threshold ⇒ NOT significant even though adjusted p = 0 (kills && → ||).
        Assert.That(OneGene(L(1, 1, 1), L(1.05, 1.05, 1.05)).IsSignificant, Is.False);

        // Adjusted-p threshold 0 ⇒ 0 < 0 is false ⇒ NOT significant (kills < → <=).
        Assert.That(OneGene(L(1, 1, 1), L(5, 5, 5), pThr: 0.0).IsSignificant, Is.False);
    }

    #endregion

    #region FindDifferentiallyExpressed — exact threshold boundaries

    [Test]
    public void FindDiffExpr_FoldChangeThresholdIsInclusive()
    {
        // CalculateFoldChange uses pseudocount 1: log2((5+1)/(1+1)) = log2(3). At threshold == log2(3)
        // the inclusive |log2FC| >= threshold must hold (kills >= → >). Welch p = 0 (zero variance).
        var gene = FindDifferentiallyExpressed(
            new[] { ("G", L(1, 1, 1), L(5, 5, 5)) },
            alpha: 0.05, log2FoldChangeThreshold: Math.Log2(3.0)).Single();
        Assert.That(gene.Log2FoldChange, Is.EqualTo(Math.Log2(3.0)).Within(Tol));
        Assert.That(gene.IsSignificant, Is.True);
    }

    [Test]
    public void FindDiffExpr_AdjustedPIsStrictlyBelowAlpha()
    {
        // adjusted p = 0; alpha = 0 ⇒ 0 < 0 false ⇒ not significant (kills < → <=).
        var gene = FindDifferentiallyExpressed(
            new[] { ("G", L(1, 1, 1), L(5, 5, 5)) },
            alpha: 0.0, log2FoldChangeThreshold: 1.0).Single();
        Assert.That(gene.IsSignificant, Is.False);
    }

    #endregion

    #region PerformOverRepresentationAnalysis (hypergeometric z-approximation)

    [Test]
    public void Ora_EnrichmentScoreAndPValueMatchHypergeometricApproximation()
    {
        var de = new HashSet<string>(Enumerable.Range(1, 10).Select(i => $"g{i}")); // 10 DE genes
        var pathwayGenes = new HashSet<string>(
            Enumerable.Range(1, 5).Select(i => $"g{i}")            // 5 overlap with DE
            .Concat(Enumerable.Range(1, 15).Select(i => $"x{i}"))); // + 15 others ⇒ size 20
        int background = 100;

        var r = PerformOverRepresentationAnalysis(
            de, new[] { ("P1", "Pathway 1", (IReadOnlySet<string>)pathwayGenes) }, background).Single();

        // Re-derived ground truth (Wikipedia hypergeometric mean/variance + normal tail).
        double expected = 10.0 * 20 / background;                                  // 2.0
        double es = 5 / expected;                                                   // 2.5
        double variance = expected * (1 - 20.0 / background) * (background - 10.0) / (background - 1);
        double z = (5 - expected) / Math.Sqrt(variance);
        double p = 1 - StatisticsHelper.NormalCDF(z);

        Assert.That(r.GenesInPathway, Is.EqualTo(20));
        Assert.That(r.OverlappingGenes, Is.EqualTo(5));
        Assert.That(r.EnrichmentScore, Is.EqualTo(es).Within(Tol));
        Assert.That(r.PValue, Is.EqualTo(p).Within(Tol));
    }

    [Test]
    public void Ora_NonPositiveBackgroundYieldsNothing()
    {
        var de = new HashSet<string> { "g1" };
        var pw = new[] { ("P1", "Pathway 1", (IReadOnlySet<string>)new HashSet<string> { "g1" }) };
        Assert.That(PerformOverRepresentationAnalysis(de, pw, 0), Is.Empty);
    }

    #endregion

    #region CalculateEnrichmentScore (GSEA running sum)

    [Test]
    public void Gsea_RunningSumPeaksWhenHitsAreTopRanked()
    {
        // 2 hits / 2 misses; hits first ⇒ running sum climbs to +1.0 then falls.
        var ranked = new[] { "g1", "g2", "g3", "g4" };
        var set = new HashSet<string> { "g1", "g2" };
        Assert.That(CalculateEnrichmentScore(ranked, set), Is.EqualTo(1.0).Within(Tol));
    }

    [Test]
    public void Gsea_RunningSumTracksLargestNegativeDeviation()
    {
        // 1 hit / 2 misses; misses first ⇒ deepest deviation is −1.0 (kills miss-count and abs-deviation logic).
        var ranked = new[] { "m1", "m2", "h1" };
        var set = new HashSet<string> { "h1" };
        Assert.That(CalculateEnrichmentScore(ranked, set), Is.EqualTo(-1.0).Within(Tol));
    }

    [Test]
    public void Gsea_NoHitsScoresZero()
        => Assert.That(CalculateEnrichmentScore(new[] { "g1", "g2" }, new HashSet<string> { "x1" }),
            Is.EqualTo(0).Within(Tol));

    #endregion

    #region Pearson correlation + co-expression network

    [Test]
    public void Pearson_PerfectPositiveAndNegative()
    {
        Assert.That(CalculatePearsonCorrelation(L(1, 2, 3), L(2, 4, 6)), Is.EqualTo(1.0).Within(Tol));
        Assert.That(CalculatePearsonCorrelation(L(1, 2, 3), L(6, 4, 2)), Is.EqualTo(-1.0).Within(Tol));
    }

    [Test]
    public void Pearson_DegenerateInputsReturnZero()
    {
        Assert.That(CalculatePearsonCorrelation(L(5), L(5)), Is.EqualTo(0).Within(Tol));        // n < 2
        Assert.That(CalculatePearsonCorrelation(L(3, 3, 3), L(1, 2, 3)), Is.EqualTo(0).Within(Tol)); // zero variance
    }

    [Test]
    public void CoExpressionNetwork_ThresholdIsInclusive()
    {
        var genes = new[]
        {
            ("A", L(1, 2, 3)),
            ("B", L(2, 4, 6)), // r(A,B) = +1
            ("C", L(3, 2, 1)), // r(A,C) = r(B,C) = −1
        };
        // |r| == 1 for every pair; threshold 1.0 keeps all 3 edges (kills >= → >).
        var edges = BuildCoExpressionNetwork(genes, correlationThreshold: 1.0).ToList();
        Assert.That(edges, Has.Count.EqualTo(3));
    }

    #endregion

    #region Alternative-splicing classification (Wang 2008)

    private static TranscriptIsoform Iso(string id, string gene, params (int Start, int End)[] exons)
        => new(id, gene, exons.Sum(e => e.End - e.Start + 1), exons.Length, 0, true, exons);

    private static string SingleEventType(params TranscriptIsoform[] isoforms)
        => DetectAlternativeSplicing(isoforms).Single().EventType;

    [Test]
    public void AltSplicing_AlternativeFivePrimeSS_SameStartDifferentEnd()
        => Assert.That(SingleEventType(Iso("t1", "G", (100, 200)), Iso("t2", "G", (100, 250))),
            Is.EqualTo("AlternativeFivePrimeSS"));

    [Test]
    public void AltSplicing_AlternativeThreePrimeSS_SameEndDifferentStart()
        => Assert.That(SingleEventType(Iso("t1", "G", (150, 300)), Iso("t2", "G", (100, 300))),
            Is.EqualTo("AlternativeThreePrimeSS"));

    [Test]
    public void AltSplicing_MutuallyExclusiveExons_NonOverlappingAlternates()
        => Assert.That(SingleEventType(
                Iso("t1", "G", (100, 150), (200, 250), (400, 450)),
                Iso("t2", "G", (100, 150), (300, 350), (400, 450))),
            Is.EqualTo("MutuallyExclusiveExons"));

    [Test]
    public void AltSplicing_TouchingAlternatesAreNotMutuallyExclusive()
        // Alternates that abut (250..250) overlap ⇒ NOT MXE ⇒ skipped exon (kills the inclusive overlap boundary).
        => Assert.That(SingleEventType(
                Iso("t1", "G", (100, 150), (200, 250), (400, 450)),
                Iso("t2", "G", (100, 150), (250, 300), (400, 450))),
            Is.EqualTo("SkippedExon"));

    [Test]
    public void AltSplicing_RetainedIntron_SpanningExonCoversIntron()
        => Assert.That(SingleEventType(Iso("t1", "G", (100, 300)), Iso("t2", "G", (100, 200), (300, 400))),
            Is.EqualTo("RetainedIntron"));

    [Test]
    public void AltSplicing_RetainedIntron_SpanStartTouchesLeftExonEnd()
        // exon.Start == left.End is still a span (inclusive lower bound).
        => Assert.That(SingleEventType(Iso("t1", "G", (200, 400)), Iso("t2", "G", (100, 200), (300, 400))),
            Is.EqualTo("RetainedIntron"));

    #endregion

    #region Differential splicing (ΔPSI)

    [Test]
    public void DifferentialSplicing_DirectionAndInclusiveThreshold()
    {
        var events = DetectDifferentialSplicing(new (string, int, int, double, double)[]
        {
            ("Inc", 1, 100, 0.2, 0.5),   // ΔPSI +0.3 ⇒ IncreasedInclusion
            ("Skp", 1, 100, 0.5, 0.2),   // ΔPSI −0.3 ⇒ IncreasedSkipping
            ("Bnd", 1, 100, 0.1, 0.2),   // ΔPSI +0.1 == threshold ⇒ reported (kills >= → >)
            ("Sub", 1, 100, 0.40, 0.45), // ΔPSI +0.05 < threshold ⇒ dropped
        }, deltaPsiThreshold: 0.1).ToList();

        Assert.That(events.Select(e => e.GeneId).OrderBy(g => g), Is.EqualTo(new[] { "Bnd", "Inc", "Skp" }));
        Assert.That(events.Single(e => e.GeneId == "Inc").EventType, Is.EqualTo("IncreasedInclusion"));
        Assert.That(events.Single(e => e.GeneId == "Skp").EventType, Is.EqualTo("IncreasedSkipping"));
        Assert.That(events.Single(e => e.GeneId == "Inc").DeltaPSI, Is.EqualTo(0.3).Within(1e-9));
    }

    #endregion

    #region Isoform dominance + switching

    [Test]
    public void DominantIsoform_RatioIsTopOverTotalExpression()
    {
        var isoforms = new[]
        {
            new TranscriptIsoform("t1", "G", 100, 1, 6, true, Array.Empty<(int, int)>()),
            new TranscriptIsoform("t2", "G", 100, 1, 3, true, Array.Empty<(int, int)>()),
            new TranscriptIsoform("t3", "G", 100, 1, 1, true, Array.Empty<(int, int)>()),
        };
        var (geneId, dominant, ratio) = FindDominantIsoforms(isoforms).Single();
        Assert.That(geneId, Is.EqualTo("G"));
        Assert.That(dominant.TranscriptId, Is.EqualTo("t1"));
        Assert.That(ratio, Is.EqualTo(6.0 / 10.0).Within(Tol));
    }

    [Test]
    public void IsoformSwitching_ReportsReciprocalUsageSwitchWithScore()
    {
        var a = new TranscriptIsoform("tA", "G", 100, 1, 0, true, Array.Empty<(int, int)>());
        var b = new TranscriptIsoform("tB", "G", 100, 1, 0, true, Array.Empty<(int, int)>());
        var data = new[]
        {
            (a, 10.0, 1.0), // usage 10/11 → 1/11, Δ ≈ −0.818
            (b, 1.0, 10.0), // usage 1/11 → 10/11, Δ ≈ +0.818
        };
        var (geneId, t1, t2, score) = DetectIsoformSwitching(data, switchThreshold: 0.3).Single();

        double expectedScore = Math.Abs(10.0 / 11 - 1.0 / 11) + Math.Abs(1.0 / 11 - 10.0 / 11);
        Assert.That(geneId, Is.EqualTo("G"));
        Assert.That(t1, Is.EqualTo("tA")); // decreased isoform first
        Assert.That(t2, Is.EqualTo("tB")); // increased isoform second
        Assert.That(score, Is.EqualTo(expectedScore).Within(Tol));
    }

    #endregion

    #region QC metrics + PCA

    [Test]
    public void QualityMetrics_RatesAndDetectedGenes()
    {
        var m = CalculateQualityMetrics(1000, 800, 600, 40, L(5, 0, 3, 0, 1));
        Assert.That(m.MappingRate, Is.EqualTo(0.8).Within(Tol));   // mapped/total
        Assert.That(m.ExonicRate, Is.EqualTo(0.75).Within(Tol));   // exonic/mapped
        Assert.That(m.RRNARate, Is.EqualTo(0.05).Within(Tol));     // rRNA/mapped
        Assert.That(m.DetectedGenes, Is.EqualTo(3));               // counts > 0
    }

    [Test]
    public void QualityMetrics_ZeroMappedReadsGivesZeroRates()
    {
        var m = CalculateQualityMetrics(1000, 0, 600, 40, L(1, 1));
        Assert.That(m.ExonicRate, Is.EqualTo(0).Within(Tol));
        Assert.That(m.RRNARate, Is.EqualTo(0).Within(Tol));
    }

    [Test]
    public void Pca_SingleSampleProjectsToOrigin()
    {
        var single = PerformPCA(new[] { ("S1", L(1, 2, 3)) }).Single();
        Assert.That(single.PC1, Is.EqualTo(0).Within(Tol));
        Assert.That(single.PC2, Is.EqualTo(0).Within(Tol));
    }

    [Test]
    public void Pca_ProjectsTopGenesAsHalvedRunningSums()
    {
        // All 4 genes selected; PC1 = sum of first half of selected values, PC2 = sum of the second half.
        var pca = PerformPCA(new[] { ("S1", L(1, 2, 3, 4)), ("S2", L(5, 6, 7, 8)) }, topGenes: 500).ToList();
        var s1 = pca.Single(s => s.SampleId == "S1");
        var s2 = pca.Single(s => s.SampleId == "S2");
        Assert.That(s1.PC1, Is.EqualTo(1 + 2).Within(Tol));
        Assert.That(s1.PC2, Is.EqualTo(3 + 4).Within(Tol));
        Assert.That(s2.PC1, Is.EqualTo(5 + 6).Within(Tol));
        Assert.That(s2.PC2, Is.EqualTo(7 + 8).Within(Tol));
    }

    #endregion
}
