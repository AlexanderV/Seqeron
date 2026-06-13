// SV-BREAKPOINT-001 — Breakpoint Detection from Split (soft-clipped) Reads
// Evidence: docs/Evidence/SV-BREAKPOINT-001-Evidence.md
// TestSpec: tests/TestSpecs/SV-BREAKPOINT-001.md
// Source: Suzuki S et al. (2011) ClipCrop, BMC Bioinformatics 12(S14):S7, doi:10.1186/1471-2105-12-S14-S7 (junction=breakpoint; cluster within 5 b);
//         Hart SN et al. (2013) SoftSearch, PLoS ONE 8(12):e83356, doi:10.1371/journal.pone.0083356 (support = clipped reads at a position; min support);
//         Tattini L et al. (2015) Front Bioeng Biotechnol 3:92, doi:10.3389/fbioe.2015.00092 (split-read signature; single-base resolution);
//         SAM/BAM Format Specification (samtools/hts-specs) — POS is contig-local; S consumes query not reference.

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.StructuralVariantAnalyzer;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class StructuralVariantAnalyzer_FindBreakpoints_Tests
{
    // Evidence defaults: junction cluster tolerance 5 b (ClipCrop), minimum support 2 (SoftSearch configurable minimum).
    private const int Tolerance = 5;
    private const int MinSupport = 2;

    // A split read whose breakpoint junction (SupplementaryPosition) is `junction` on `chromosome`.
    private static SplitRead Read(string chromosome, int junction, string readId = "r") =>
        new(
            ReadId: readId,
            Chromosome: chromosome,
            PrimaryPosition: junction - 500, // anchored position; irrelevant to junction clustering
            SupplementaryPosition: junction,
            ClipLength: 40,
            ClippedSequence: "ACGTACGTAC");

    #region FindBreakpoints

    // M1 — Junctions {5000,5002,5004} within tolerance 5, same chr, minSupport 2 => one breakpoint at mean 5002, support 3.
    [Test]
    public void FindBreakpoints_AgreeingJunctions_ReturnsOneBreakpointWithSupport()
    {
        var reads = new[]
        {
            Read("chr1", 5000, "r1"),
            Read("chr1", 5002, "r2"),
            Read("chr1", 5004, "r3"),
        };

        var breakpoints = FindBreakpoints(reads, Tolerance, MinSupport).ToList();

        Assert.That(breakpoints, Has.Count.EqualTo(1),
            "Three junctions agreeing within the 5-base cluster tolerance form one breakpoint (ClipCrop: clustered within 5-base differences).");
        Assert.Multiple(() =>
        {
            Assert.That(breakpoints[0].Position1, Is.EqualTo(5002),
                "The reported coordinate is the rounded mean of member junctions {5000,5002,5004} = 5002 (ASM-01).");
            Assert.That(breakpoints[0].SupportingReads, Is.EqualTo(3),
                "Support equals the number of clipped reads in the cluster (SoftSearch: reads beginning at the position).");
            Assert.That(breakpoints[0].Chromosome1, Is.EqualTo("chr1"),
                "The breakpoint is reported on the contig of its member reads (SAM POS is contig-local).");
        });
    }

    // M2 — A single isolated split read with minSupport 2 yields no breakpoint. SoftSearch "at least x soft-clipped reads".
    [Test]
    public void FindBreakpoints_BelowMinSupport_ReturnsEmpty()
    {
        var reads = new[] { Read("chr1", 5000) };

        var breakpoints = FindBreakpoints(reads, Tolerance, MinSupport).ToList();

        Assert.That(breakpoints, Is.Empty,
            "A clip position with fewer than the minimum supporting reads (default 2) is not reported (SoftSearch).");
    }

    // M3 — Junctions 5000 and 5100 are 100 b apart (> tolerance 5): two singleton groups, each below support => none. ClipCrop tolerance.
    [Test]
    public void FindBreakpoints_GapExceedsTolerance_ReturnsEmpty()
    {
        var reads = new[]
        {
            Read("chr1", 5000, "r1"),
            Read("chr1", 5100, "r2"),
        };

        var breakpoints = FindBreakpoints(reads, Tolerance, MinSupport).ToList();

        Assert.That(breakpoints, Is.Empty,
            "Junctions farther apart than the cluster tolerance are not merged; each singleton falls below min support (ClipCrop clusters only within tolerance).");
    }

    // M4 — Two clusters (~5000 x2 and ~9000 x2), minSupport 2 => two breakpoints, each support 2. ClipCrop sort/cluster; SoftSearch support.
    [Test]
    public void FindBreakpoints_TwoSeparateClusters_ReturnsTwoBreakpoints()
    {
        var reads = new[]
        {
            Read("chr1", 5000, "a1"),
            Read("chr1", 5002, "a2"),
            Read("chr1", 9000, "b1"),
            Read("chr1", 9002, "b2"),
        };

        var breakpoints = FindBreakpoints(reads, Tolerance, MinSupport)
            .OrderBy(b => b.Position1)
            .ToList();

        Assert.That(breakpoints, Has.Count.EqualTo(2),
            "Two clip stacks separated by more than the tolerance form two distinct breakpoints (ClipCrop sorts and clusters by position).");
        Assert.Multiple(() =>
        {
            Assert.That(breakpoints[0].Position1, Is.EqualTo(5001),
                "First breakpoint is the rounded mean of {5000,5002} = 5001 (ASM-01).");
            Assert.That(breakpoints[0].SupportingReads, Is.EqualTo(2), "First cluster has 2 supporting reads.");
            Assert.That(breakpoints[1].Position1, Is.EqualTo(9001),
                "Second breakpoint is the rounded mean of {9000,9002} = 9001 (ASM-01).");
            Assert.That(breakpoints[1].SupportingReads, Is.EqualTo(2), "Second cluster has 2 supporting reads.");
        });
    }

    // M5 — Same junction 5000 on chr1 (x2) and chr2 (x2) => two breakpoints, one per chromosome. SAM POS contig-local.
    [Test]
    public void FindBreakpoints_DifferentChromosomes_NotMerged()
    {
        var reads = new[]
        {
            Read("chr1", 5000, "a1"),
            Read("chr1", 5000, "a2"),
            Read("chr2", 5000, "b1"),
            Read("chr2", 5000, "b2"),
        };

        var breakpoints = FindBreakpoints(reads, Tolerance, MinSupport)
            .OrderBy(b => b.Chromosome1, StringComparer.Ordinal)
            .ToList();

        Assert.That(breakpoints, Has.Count.EqualTo(2),
            "Identical junction coordinates on different chromosomes are distinct breakpoints (SAM POS is contig-local).");
        Assert.Multiple(() =>
        {
            Assert.That(breakpoints[0].Chromosome1, Is.EqualTo("chr1"), "One breakpoint is reported on chr1.");
            Assert.That(breakpoints[1].Chromosome1, Is.EqualTo("chr2"), "The other breakpoint is reported on chr2.");
        });
    }

    // M6 — Empty input yields empty output (defined trivial behavior).
    [Test]
    public void FindBreakpoints_EmptyInput_ReturnsEmpty()
    {
        var reads = Array.Empty<SplitRead>();

        var breakpoints = FindBreakpoints(reads, Tolerance, MinSupport).ToList();

        Assert.That(breakpoints, Is.Empty, "No split reads means no breakpoints can be detected.");
    }

    // M7 — Null input throws ArgumentNullException (input-validation contract, eager before iteration).
    [Test]
    public void FindBreakpoints_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => FindBreakpoints(null!, Tolerance, MinSupport).ToList(),
            "Null split-read input violates the precondition and must throw ArgumentNullException.");
    }

    // S1 — Junctions exactly `tolerance` apart (5000, 5005) cluster (window is inclusive). ClipCrop within-5-base window.
    [Test]
    public void FindBreakpoints_JunctionsExactlyToleranceApart_Cluster()
    {
        var reads = new[]
        {
            Read("chr1", 5000, "r1"),
            Read("chr1", 5005, "r2"),
        };

        var breakpoints = FindBreakpoints(reads, Tolerance, MinSupport).ToList();

        Assert.That(breakpoints, Has.Count.EqualTo(1),
            "A gap exactly equal to the tolerance (5) is inside the cluster window; the two reads form one breakpoint.");
        Assert.That(breakpoints[0].SupportingReads, Is.EqualTo(2),
            "Both reads support the single clustered breakpoint.");
    }

    // S2 — Junctions tolerance+1 apart (5000, 5006) do not cluster => two singletons below support => none.
    [Test]
    public void FindBreakpoints_JunctionsBeyondTolerance_DoNotCluster()
    {
        var reads = new[]
        {
            Read("chr1", 5000, "r1"),
            Read("chr1", 5006, "r2"),
        };

        var breakpoints = FindBreakpoints(reads, Tolerance, MinSupport).ToList();

        Assert.That(breakpoints, Is.Empty,
            "A gap one base beyond the tolerance (6 > 5) splits the reads into separate singleton groups, each below min support.");
    }

    #endregion

    #region RefineBreakpoint

    // M8 — Region [4990,5010] over junctions {5000,5000,5004} => consensus = mode 5000. SoftSearch: reads accumulate at the breakpoint.
    [Test]
    public void RefineBreakpoint_ConsensusJunction_ReturnsMode()
    {
        var reads = new[]
        {
            Read("chr1", 5000, "r1"),
            Read("chr1", 5000, "r2"),
            Read("chr1", 5004, "r3"),
        };

        var refined = RefineBreakpoint("chr1", 4990, 5010, reads);

        Assert.That(refined, Is.EqualTo(5000),
            "The consensus is the junction shared by the most reads (mode 5000), the position where clipped reads accumulate (SoftSearch).");
    }

    // S3 — Region with no member junctions returns null (no support to form a consensus).
    [Test]
    public void RefineBreakpoint_NoReadsInRegion_ReturnsNull()
    {
        var reads = new[]
        {
            Read("chr1", 5000, "r1"),
            Read("chr1", 5004, "r2"),
        };

        var refined = RefineBreakpoint("chr1", 8000, 8100, reads);

        Assert.That(refined, Is.Null,
            "A region containing no split-read junction has no consensus breakpoint and must return null.");
    }

    // C1 — Null input throws ArgumentNullException (input-validation contract).
    [Test]
    public void RefineBreakpoint_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => RefineBreakpoint("chr1", 0, 100, null!),
            "Null split-read input violates the precondition and must throw ArgumentNullException.");
    }

    #endregion
}
