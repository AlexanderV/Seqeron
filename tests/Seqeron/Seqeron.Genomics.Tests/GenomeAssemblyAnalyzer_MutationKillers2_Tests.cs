using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Chromosome;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// ASSEMBLY-STATS-001 mutation killers (batch 2): exact-value tests for the remaining uncovered
/// analysis methods — scaffold structure, k-mer completeness/error estimation, tandem-repeat finding,
/// local-quality (GC / N / linguistic-complexity) windows, BUSCO-like completeness, repetitive- and
/// suspicious-region detection.
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_MutationKillers2_Tests
{
    private const double Tol = 1e-9;

    #region AnalyzeScaffolds

    [Test]
    public void AnalyzeScaffolds_SplitsContigsAndGapsWithExactSums()
    {
        // 3 A · 10 N · 3 A: two contigs (len 3 each) separated by a 10-bp gap.
        var s = AnalyzeScaffolds(new[] { ("sc", "AAANNNNNNNNNNAAA") }, minGapLength: 10).Single();

        Assert.That(s.ScaffoldId, Is.EqualTo("sc"));
        Assert.That(s.Contigs.Count, Is.EqualTo(2));
        Assert.That(s.Contigs[0], Is.EqualTo(("sc_contig1", 0, 2)));
        Assert.That(s.Contigs[1], Is.EqualTo(("sc_contig2", 13, 15)));
        Assert.That(s.Gaps.Count, Is.EqualTo(1));
        Assert.That(s.Gaps[0].Length, Is.EqualTo(10));
        Assert.That(s.TotalLength, Is.EqualTo(16));
        Assert.That(s.ContigLength, Is.EqualTo(6));
        Assert.That(s.GapLength, Is.EqualTo(10));
    }

    #endregion

    #region EstimateCompletenessFromKmers

    [Test]
    public void EstimateCompletenessFromKmers_PeakErrorAndSize()
    {
        // Spectrum: AAA×4, CCC×4, GGG×2, TTT×1.
        // Non-singleton mode → peak coverage 4; singleton ratio 1/4 = error 0.25;
        // total k-mers 11 / peak 4 = size 2; solid (≥peak/2=2) = 3 of (4−1) distinct ⇒ completeness 1.
        var (completeness, errorRate, size) = EstimateCompletenessFromKmers(new[]
        {
            ("AAA", 4), ("CCC", 4), ("GGG", 2), ("TTT", 1),
        });

        Assert.That(completeness, Is.EqualTo(1.0).Within(Tol));
        Assert.That(errorRate, Is.EqualTo(0.25).Within(Tol));
        Assert.That(size, Is.EqualTo(2L));
    }

    #endregion

    #region FindTandemRepeats

    [Test]
    public void FindTandemRepeats_PerfectDinucleotideRepeat()
    {
        // (AT)×4 = 8 bp, perfect ⇒ purity 1.0.
        var r = FindTandemRepeats(new[] { ("s", "ATATATAT") }).Single();
        Assert.That(r.SequenceId, Is.EqualTo("s"));
        Assert.That(r.Start, Is.EqualTo(0));
        Assert.That(r.End, Is.EqualTo(7));
        Assert.That(r.Unit, Is.EqualTo("AT"));
        Assert.That(r.Copies, Is.EqualTo(4));
        Assert.That(r.Purity, Is.EqualTo(1.0).Within(Tol));
    }

    #endregion

    #region CalculateLocalQuality (GC / N / linguistic complexity)

    [Test]
    public void CalculateLocalQuality_ExactWindowMetrics()
    {
        // "ACGTACGT": one window (len 8 < default 1000). GC = 4/8 = 0.5; N = 0;
        // linguistic complexity = distinct 4-mers (4) / min(validKmers 5, 4^4) = 4/5 = 0.8.
        var q = CalculateLocalQuality(new[] { ("s", "ACGTACGT") }).Single();
        Assert.That(q.SequenceId, Is.EqualTo("s"));
        Assert.That(q.Position, Is.EqualTo(0));
        Assert.That(q.WindowSize, Is.EqualTo(8));
        Assert.That(q.GcContent, Is.EqualTo(0.5).Within(Tol));
        Assert.That(q.NCount, Is.EqualTo(0));
        Assert.That(q.Complexity, Is.EqualTo(0.8).Within(Tol));
    }

    #endregion

    #region AssessCompleteness (BUSCO-like)

    [Test]
    public void AssessCompleteness_CompleteAndMissingMarkers()
    {
        // One marker is present (its sequence is a substring of the assembly) and one is absent.
        string present = new string('A', 40) + new string('C', 40); // 80-bp marker, 80-bp kmers
        string assembly = present;                                    // identical ⇒ full coverage
        string absent = new string('G', 80);                          // shares no k-mers

        var r = AssessCompleteness(
            new[] { ("scaf", assembly) },
            new[] { ("g1", present), ("g2", absent) });

        Assert.That(r.TotalGenes, Is.EqualTo(2));
        Assert.That(r.Complete, Is.EqualTo(1));
        Assert.That(r.Missing, Is.EqualTo(1));
        Assert.That(r.CompletenessPercent, Is.EqualTo(50.0).Within(Tol));

        // Empty marker set ⇒ all-zero result.
        var none = AssessCompleteness(new[] { ("scaf", assembly) }, Array.Empty<(string, string)>());
        Assert.That(none.TotalGenes, Is.EqualTo(0));
    }

    #endregion

    #region Region detectors (behavioural coverage)

    [Test]
    public void FindRepetitiveRegions_DetectsHighCopyKmerRun()
    {
        // A long tandem (CAG)n produces many high-copy k-mers ⇒ at least one repetitive region.
        string seq = string.Concat(Enumerable.Repeat("CAG", 200));
        var regions = FindRepetitiveRegions(new[] { ("s", seq) }, kmerSize: 6, minCopies: 3, windowSize: 60).ToList();
        Assert.That(regions, Is.Not.Empty);
        Assert.That(regions.All(r => r.SequenceId == "s"), Is.True);
        Assert.That(regions.All(r => r.Copies >= 3), Is.True);
    }

    [Test]
    public void FindSuspiciousRegions_FlagsLowComplexityHomopolymer()
    {
        // A 600-bp poly(A) run is low-complexity ⇒ a suspicious region citing "Low complexity".
        var regions = FindSuspiciousRegions(new[] { ("s", new string('A', 600)) }).ToList();
        Assert.That(regions, Is.Not.Empty);
        Assert.That(regions.Any(r => r.Reason.Contains("Low complexity")), Is.True);
    }

    #endregion
}
