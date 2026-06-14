// ONCO-FUSION-001 — Fusion Gene Detection
// Evidence: docs/Evidence/ONCO-FUSION-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-FUSION-001.md
// Source: Haas et al. (2017) STAR-Fusion (MIN_JUNCTION_READS=1, MIN_SUM_FRAGS=2, MIN_SPANNING_FRAGS_ONLY=5);
//         Uhrig et al. (2021) Arriba, Genome Research 31(3):448 (total support = split1+split2+discordant);
//         Genomics England + Wikipedia "Reading frame" (in-frame iff codon phase preserved, modulo 3).

using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Oncology;
using Candidate = Seqeron.Genomics.Oncology.OncologyAnalyzer.FusionCandidate;
using Frame = Seqeron.Genomics.Oncology.OncologyAnalyzer.FusionReadingFrame;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class OncologyAnalyzer_DetectFusions_Tests
{
    #region DetectFusions — minimum-support thresholds

    // M1 — STAR-Fusion: junction reads (5) ≥ 1 AND total support (9) ≥ 2 → detected; total = 3+2+4 = 9 (Arriba).
    [Test]
    public void DetectFusions_JunctionAndSumPass_Detected()
    {
        var candidates = new[] { new Candidate("EML4", "ALK", 3, 2, 4) };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.Multiple(() =>
        {
            Assert.That(calls, Has.Count.EqualTo(1),
                "EML4-ALK has 5 junction reads and total support 9, passing MIN_JUNCTION_READS=1 and MIN_SUM_FRAGS=2.");
            Assert.That(calls[0].Gene5Prime, Is.EqualTo("EML4"), "5' partner preserved.");
            Assert.That(calls[0].Gene3Prime, Is.EqualTo("ALK"), "3' partner preserved.");
            Assert.That(calls[0].JunctionReads, Is.EqualTo(5), "Junction reads = split1 + split2 = 3 + 2.");
            Assert.That(calls[0].TotalSupport, Is.EqualTo(9), "Total support = split1 + split2 + discordant = 3 + 2 + 4 (Arriba).");
        });
    }

    // M2 — STAR-Fusion: 0 junction reads, 5 discordant fragments → meets MIN_SPANNING_FRAGS_ONLY=5 → detected.
    [Test]
    public void DetectFusions_SpanningOnlyFiveFrags_Detected()
    {
        var candidates = new[] { new Candidate("CD74", "ROS1", 0, 0, 5) };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.Multiple(() =>
        {
            Assert.That(calls, Has.Count.EqualTo(1),
                "With no junction reads, 5 discordant fragments meets MIN_SPANNING_FRAGS_ONLY=5.");
            Assert.That(calls[0].TotalSupport, Is.EqualTo(5), "Total support = 0 + 0 + 5.");
        });
    }

    // M3 — STAR-Fusion: 0 junction reads, only 4 discordant fragments → below MIN_SPANNING_FRAGS_ONLY=5 → rejected.
    [Test]
    public void DetectFusions_SpanningOnlyFourFrags_Rejected()
    {
        var candidates = new[] { new Candidate("NCOA4", "RET", 0, 0, 4) };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.That(calls, Is.Empty,
            "4 discordant fragments with no junction reads is below MIN_SPANNING_FRAGS_ONLY=5 → not reported.");
    }

    // M4 — STAR-Fusion: 1 junction read but total support 1 < MIN_SUM_FRAGS=2 → rejected.
    [Test]
    public void DetectFusions_SumBelowTwo_Rejected()
    {
        var candidates = new[] { new Candidate("KIF5B", "RET", 1, 0, 0) };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.That(calls, Is.Empty,
            "Junction reads = 1 satisfies MIN_JUNCTION_READS but total support 1 < MIN_SUM_FRAGS=2 → not reported.");
    }

    // M5 — STAR-Fusion: junction=1 (≥1) AND total=2 (≥2) → detected (the minimal passing junction case).
    [Test]
    public void DetectFusions_JunctionOneSumTwo_Detected()
    {
        var candidates = new[] { new Candidate("TMPRSS2", "ERG", 1, 0, 1) };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.Multiple(() =>
        {
            Assert.That(calls, Has.Count.EqualTo(1),
                "1 junction read and total support 2 meet MIN_JUNCTION_READS=1 and MIN_SUM_FRAGS=2.");
            Assert.That(calls[0].TotalSupport, Is.EqualTo(2), "Total support = 1 + 0 + 1.");
        });
    }

    // M6 — Registry invariant INV-1: a gene fused with itself is not a fusion, even with very high support.
    [Test]
    public void DetectFusions_SameGene_Rejected()
    {
        var candidates = new[] { new Candidate("ALK", "ALK", 5, 5, 5) };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.That(calls, Is.Empty,
            "gene5p == gene3p is not a gene fusion (INV-1) regardless of supporting-read count.");
    }

    // M7 — Arriba: TotalSupport = split_reads1 + split_reads2 + discordant_mates.
    [Test]
    public void ComputeTotalSupport_SumOfThreeClasses()
    {
        int total = OncologyAnalyzer.ComputeTotalSupport(new Candidate("A", "B", 3, 2, 4));

        Assert.That(total, Is.EqualTo(9),
            "Total support is the sum of split_reads1 (3), split_reads2 (2) and discordant_mates (4) = 9 (Arriba).");
    }

    // M11 / INV-4 — results ordered by descending total support (abundance of supporting reads).
    [Test]
    public void DetectFusions_MixedCandidates_OrderedByDescendingSupport()
    {
        var candidates = new[]
        {
            new Candidate("TMPRSS2", "ERG", 1, 0, 1),  // total 2
            new Candidate("EML4", "ALK", 3, 2, 4),     // total 9
            new Candidate("CD74", "ROS1", 0, 0, 5),    // total 5
            new Candidate("NCOA4", "RET", 0, 0, 4),    // rejected (4 < 5)
        };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.Multiple(() =>
        {
            Assert.That(calls.Select(c => c.TotalSupport), Is.EqualTo(new[] { 9, 5, 2 }),
                "Detected fusions are ordered by descending total support (STAR-Fusion scores by abundance); NCOA4-RET is filtered.");
            Assert.That(calls[0].Gene5Prime, Is.EqualTo("EML4"), "Highest-support fusion (9) is first.");
            Assert.That(calls[2].Gene5Prime, Is.EqualTo("TMPRSS2"), "Lowest-support fusion (2) is last.");
        });
    }

    #endregion

    #region DetectFusions — input validation and boundaries

    // M12 — null candidates → ArgumentNullException (sibling-method convention).
    [Test]
    public void DetectFusions_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => OncologyAnalyzer.DetectFusions(null!),
            "Null candidate collection must throw ArgumentNullException.");
    }

    // M13 — a negative supporting-read count is invalid → ArgumentException.
    [Test]
    public void DetectFusions_NegativeCount_Throws()
    {
        var candidates = new[] { new Candidate("A", "B", -1, 0, 0) };

        Assert.Throws<ArgumentException>(() => OncologyAnalyzer.DetectFusions(candidates),
            "A negative supporting-read count is invalid input and must throw ArgumentException.");
    }

    // M14 — empty input → empty result (trivial).
    [Test]
    public void DetectFusions_EmptyInput_ReturnsEmpty()
    {
        IReadOnlyList<OncologyAnalyzer.FusionCall> calls =
            OncologyAnalyzer.DetectFusions(Array.Empty<Candidate>());

        Assert.That(calls, Is.Empty, "No candidates → no fusion calls.");
    }

    // S1 — STAR-Fusion thresholds are configurable: with min_spanning_frags_only=3, 4 discordant passes.
    [Test]
    public void DetectFusions_CustomSpanningThreshold_Detected()
    {
        var candidates = new[] { new Candidate("NCOA4", "RET", 0, 0, 4) };
        var thresholds = new OncologyAnalyzer.FusionDetectionThresholds(MinSpanningFragsOnly: 3);

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates, thresholds);

        Assert.That(calls, Has.Count.EqualTo(1),
            "With min_spanning_frags_only lowered to 3, 4 discordant fragments (no junction reads) now passes.");
    }

    // S2 — boundary: exactly 5 discordant fragments with no junction reads passes (≥, not >).
    [Test]
    public void DetectFusions_SpanningExactlyFive_Detected()
    {
        var candidates = new[] { new Candidate("CD74", "ROS1", 0, 0, 5) };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.That(calls, Has.Count.EqualTo(1),
            "Exactly MIN_SPANNING_FRAGS_ONLY=5 discordant fragments is sufficient (threshold is ≥).");
    }

    // S3 — boundary: junction=1 and total support exactly 2 passes (≥, not >).
    [Test]
    public void DetectFusions_SumExactlyTwo_Detected()
    {
        var candidates = new[] { new Candidate("TMPRSS2", "ERG", 2, 0, 0) };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.That(calls, Has.Count.EqualTo(1),
            "Total support exactly MIN_SUM_FRAGS=2 (2 junction reads) is sufficient (threshold is ≥).");
    }

    #endregion

    #region Reading frame (codon phase)

    // M8 — in-frame: 300 coding bases (multiple of 3), start phase 0 → preserved frame.
    [Test]
    public void IsInFrame_Phase0_True()
    {
        Assert.That(OncologyAnalyzer.IsInFrame(300, 0), Is.True,
            "(300 - 0) mod 3 == 0 → the 3' partner stays in codon phase → in-frame.");
    }

    // M9 — out-of-frame: 301 coding bases, start phase 0 → frame shifted by 1.
    [Test]
    public void IsInFrame_Phase1_False()
    {
        Assert.That(OncologyAnalyzer.IsInFrame(301, 0), Is.False,
            "(301 - 0) mod 3 == 1 → reading frame shifted by one base → out-of-frame.");
    }

    // M10 — in-frame with nonzero 3' start phase: 301 coding bases, start phase 1 → preserved.
    [Test]
    public void IsInFrame_NonzeroStartPhase_True()
    {
        Assert.That(OncologyAnalyzer.IsInFrame(301, 1), Is.True,
            "(301 - 1) mod 3 == 0 → matching the 3' partner's start phase 1 keeps codons in phase → in-frame.");
    }

    // S4 — out-of-frame: 302 coding bases, start phase 0 → frame shifted by 2.
    [Test]
    public void IsInFrame_Phase2_False()
    {
        Assert.That(OncologyAnalyzer.IsInFrame(302, 0), Is.False,
            "(302 - 0) mod 3 == 2 → reading frame shifted by two bases → out-of-frame.");
    }

    // Reading frame is surfaced on the FusionCall when coding-phase fields are supplied.
    [Test]
    public void DetectFusions_WithCodingPhase_ReportsReadingFrame()
    {
        var candidates = new[]
        {
            new Candidate("EML4", "ALK", 3, 2, 4, FivePrimeCodingBases: 300, ThreePrimeStartPhase: 0), // in-frame
            new Candidate("CD74", "ROS1", 3, 2, 4, FivePrimeCodingBases: 301, ThreePrimeStartPhase: 0), // out-of-frame
        };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.Multiple(() =>
        {
            Assert.That(calls.Single(c => c.Gene5Prime == "EML4").ReadingFrame, Is.EqualTo(Frame.InFrame),
                "EML4-ALK with (300,0) is in-frame.");
            Assert.That(calls.Single(c => c.Gene5Prime == "CD74").ReadingFrame, Is.EqualTo(Frame.OutOfFrame),
                "CD74-ROS1 with (301,0) is out-of-frame.");
        });
    }

    // When coding-phase fields are unset (-1), reading frame is Unknown (never guessed).
    [Test]
    public void DetectFusions_NoCodingPhase_ReadingFrameUnknown()
    {
        var candidates = new[] { new Candidate("EML4", "ALK", 3, 2, 4) };

        IReadOnlyList<OncologyAnalyzer.FusionCall> calls = OncologyAnalyzer.DetectFusions(candidates);

        Assert.That(calls[0].ReadingFrame, Is.EqualTo(Frame.Unknown),
            "Without coding-phase information the reading frame is reported Unknown, not guessed.");
    }

    // C1 — IsInFrame rejects invalid arguments (negative base count / phase outside {0,1,2}).
    [Test]
    public void IsInFrame_NegativeArgument_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.IsInFrame(-1, 0),
                "Negative coding-base count is invalid.");
            Assert.Throws<ArgumentOutOfRangeException>(() => OncologyAnalyzer.IsInFrame(10, 3),
                "Phase must be 0, 1, or 2; phase 3 is invalid.");
        });
    }

    #endregion
}
