using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.StructuralVariantAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for split-read breakpoint detection — SV-BREAKPOINT-001. The unit
/// under test is the signature-then-cluster breakpoint localizer in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs:
///   • <see cref="StructuralVariantAnalyzer.FindBreakpoints"/> — sorts split-read
///     junctions per chromosome, groups adjacent junctions within a tolerance window
///     and emits one <see cref="StructuralVariantAnalyzer.Breakpoint"/> per cluster
///     that meets a minimum read-support gate (canonical);
///   • <see cref="StructuralVariantAnalyzer.RefineBreakpoint"/> — returns the modal
///     (tie → rounded-mean) junction inside a candidate region, or null.
///
/// SCOPE. This file is scoped strictly to SV-BREAKPOINT-001 — single-base breakpoint
/// localization from soft-clip-derived junctions. The per-read junction key is
/// <see cref="StructuralVariantAnalyzer.SplitRead.SupplementaryPosition"/> (the
/// aligned/clipped junction), NOT the anchored PrimaryPosition that the sibling
/// ClusterSplitReads groups on (docs §5.2). Paired-end SV typing
/// (DEL/DUP/INV/TRA — SV-DETECT-001) and junction-sequence assembly are SEPARATE
/// units and are not exercised here (docs §5.3 "Not implemented").
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary inputs to a unit and asserts that the code
/// NEVER fails in an undisciplined way: no hang, no nonsense output, no *unhandled*
/// runtime exception (e.g. an Average over an empty cluster, an overflow on a
/// MaxInt junction coordinate, a NaN/Infinity leaking into a reported position), and
/// no silent corruption of the support count. Every input must resolve to EITHER a
/// well-defined, theory-correct value OR a documented, intentional outcome
/// (ArgumentNullException for a null read sequence). The headline hazards for THIS
/// unit are:
///   • empty input → an empty, non-null breakpoint sequence; no Average-of-empty
///     DivideByZero / "Sequence contains no elements" (§3.3, §6.1, INV-05);
///   • a single split read → with the default minSupport 2 it is below the gate ⇒
///     NO breakpoint; no variance-of-one / one-element-cluster crash (§6.1, INV-01);
///   • identical junctions / no split reads → identical junctions collapse to ONE
///     breakpoint at that coordinate; zero split reads supporting a refine region ⇒
///     null, never a fabricated junction (§6.1, INV-02/INV-04);
///   • MaxInt / extreme coordinates → the rounded-mean position and the tolerance gap
///     |curr − prev| must not overflow into a negative gap that mis-merges clusters
///     (Boundary Exploitation MaxInt; INV-02/INV-03).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SV-BREAKPOINT-001 — Breakpoint detection from split reads (StructuralVar)
/// Checklist: docs/checklists/03_FUZZING.md, row 201.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — граничні значення: 0, -1, MaxInt, empty.
///     Targets (checklist row 201): "identical, single breakpoint, no split reads".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (Breakpoint_Detection.md)
/// (docs/algorithms/StructuralVar/Breakpoint_Detection.md)
/// ───────────────────────────────────────────────────────────────────────────
///   • per-read breakpoint = the junction coordinate SupplementaryPosition; single-
///     base resolution                                              (§2.2, §4.1[1])
///   • sort by chromosome (ordinal) then junction; new cluster on a chromosome change
///     OR a gap > clusterTolerance to the previous junction          (§4.1[3], INV-02)
///   • a breakpoint is emitted iff its cluster size ≥ minSupport     (§4.1[4], INV-01)
///   • reported Position = rounded mean of member junctions (ASM-01); it lies within
///     ≤ tolerance of every member                                  (§2.4 INV-03)
///   • SupportingReads = cluster size                               (§2.4 INV-04)
///   • #breakpoints ≤ #input split reads (a partition)             (§2.4 INV-05)
///   • defaults: clusterTolerance 5, minSupport 2                   (§3.1, §4.2)
///   • RefineBreakpoint: modal junction in [start,end] on chromosome (tie → rounded
///     mean); null if no member junction                            (§4.1[5], §6.1)
///   • null splitReads ⇒ ArgumentNullException (both methods)        (§3.3)
///   • worked example: junctions {5000,5002,5004} agree within tolerance 5 ⇒ one
///     breakpoint at 5002, support 3                                (§7.1)
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class StructuralVariantBreakpointFuzzTests
{
    // Documented defaults (§3.1/§4.2), mirrored locally so the test owns the gate
    // independently of the source constants.
    private const int DefaultTolerance = 5;   // ClipCrop "within 5-base differences"
    private const int DefaultMinSupport = 2;  // SoftSearch configurable minimum

    // Builds a split read whose junction (SupplementaryPosition) is the breakpoint
    // key. PrimaryPosition is deliberately set to an unrelated value so a test that
    // accidentally clustered on the anchored coordinate (the sibling key) would fail.
    private static SplitRead Read(string id, string chrom, int junction)
        => new(id, chrom, PrimaryPosition: junction + 500, SupplementaryPosition: junction,
               ClipLength: 20, ClippedSequence: "ACGT");

    // ── Well-formed-result assertion helper ──────────────────────────────────
    // Pins the documented structural invariants on EVERY accepted breakpoint set,
    // regardless of input. This is what stops a fuzz test from rubber-stamping a
    // green call: each reported position must be finite (no NaN/Infinity from a
    // rounded mean), support must equal a real cluster size ≥ the gate, and the
    // total count cannot exceed the input read count (partition).
    private static void AssertWellFormed(
        IReadOnlyList<Breakpoint> breakpoints, int inputReadCount, int minSupport)
    {
        breakpoints.Count.Should().BeLessThanOrEqualTo(
            inputReadCount, "each read joins exactly one cluster ⇒ a partition (INV-05)");

        foreach (var bp in breakpoints)
        {
            bp.SupportingReads.Should().BeGreaterThanOrEqualTo(
                minSupport, "every emitted breakpoint passed the minimum-support gate (INV-01)");

            bp.SupportingReads.Should().BeGreaterThan(0, "a cluster has at least one member");

            // Position is an int (no NaN/Infinity is representable), but assert the
            // record is intra-contig and the strand convention is the documented one.
            bp.Chromosome1.Should().Be(bp.Chromosome2,
                "intra-contig clustering ⇒ both ends share the contig (§3.2)");
            bp.Strand1.Should().Be('+');
            bp.Strand2.Should().Be('-');
            double.IsNaN(bp.Quality).Should().BeFalse("quality is a finite score");
            double.IsInfinity(bp.Quality).Should().BeFalse("quality is a finite score");
        }
    }

    #region SV-BREAKPOINT-001 — Positive sanity (documented worked example)

    [Test]
    public void FindBreakpoints_DocumentedWorkedExample_SingleBreakpointSupportThree()
    {
        // Docs §7.1: junctions {5000, 5002, 5004} all within tolerance 5 ⇒ ONE
        // cluster; rounded-mean position = round((5000+5002+5004)/3) = round(5002) =
        // 5002 (ASM-01); support = cluster size = 3 (INV-04). Hand-checked.
        var reads = new[]
        {
            Read("r1", "chr1", 5000),
            Read("r2", "chr1", 5002),
            Read("r3", "chr1", 5004),
        };

        var breakpoints = FindBreakpoints(reads, clusterTolerance: 5, minSupport: 2).ToList();

        breakpoints.Should().HaveCount(1, "all three junctions agree within tolerance 5 (§7.1)");
        breakpoints[0].Position1.Should().Be(5002, "rounded mean of {5000,5002,5004} (ASM-01)");
        breakpoints[0].Position2.Should().Be(5002);
        breakpoints[0].Chromosome1.Should().Be("chr1");
        breakpoints[0].SupportingReads.Should().Be(3, "cluster size = 3 (INV-04)");
        AssertWellFormed(breakpoints, inputReadCount: 3, minSupport: 2);
    }

    [Test]
    public void FindBreakpoints_TwoJunctionsBeyondTolerance_AreSeparateClusters()
    {
        // Junctions {1000, 1006} are 6 > tolerance 5 apart ⇒ two singleton clusters
        // (INV-02). With minSupport 1 both are reported, each at its own coordinate
        // with support 1. This pins the gap > tolerance split rule (§6.1).
        var reads = new[] { Read("a", "chr1", 1000), Read("b", "chr1", 1006) };

        var breakpoints = FindBreakpoints(reads, clusterTolerance: 5, minSupport: 1).ToList();

        breakpoints.Should().HaveCount(2, "gap 6 > tolerance 5 ⇒ separate clusters (INV-02)");
        breakpoints.Select(b => b.Position1).Should().BeEquivalentTo(new[] { 1000, 1006 });
        AssertWellFormed(breakpoints, inputReadCount: 2, minSupport: 1);
    }

    [Test]
    public void FindBreakpoints_SameJunctionDifferentChromosomes_AreSeparateBreakpoints()
    {
        // POS is contig-local (§2.2[1], §6.1): identical junction 2000 on chr1 and
        // chr2 must NOT merge despite a zero positional gap. Two singletons.
        var reads = new[]
        {
            Read("x", "chr1", 2000), Read("y", "chr1", 2000),
            Read("z", "chr2", 2000), Read("w", "chr2", 2000),
        };

        var breakpoints = FindBreakpoints(reads, clusterTolerance: 5, minSupport: 2).ToList();

        breakpoints.Should().HaveCount(2, "same coordinate on different contigs ⇒ separate (INV-02)");
        breakpoints.Select(b => b.Chromosome1).Should().BeEquivalentTo(new[] { "chr1", "chr2" });
        breakpoints.Should().OnlyContain(b => b.Position1 == 2000 && b.SupportingReads == 2);
        AssertWellFormed(breakpoints, inputReadCount: 4, minSupport: 2);
    }

    #endregion

    #region SV-BREAKPOINT-001 — BE boundary: empty input / no split reads

    [Test]
    public void FindBreakpoints_EmptyInput_YieldsEmptyNonNullSequence()
    {
        // §6.1: no junctions to cluster ⇒ empty result. The Average-of-empty hazard
        // must not surface as a DivideByZero / "no elements" — the early yield break
        // guards it (INV-05 vacuously).
        var breakpoints = FindBreakpoints(Array.Empty<SplitRead>(), clusterTolerance: 5, minSupport: 2);

        breakpoints.Should().NotBeNull().And.BeEmpty("no junctions ⇒ empty (§6.1)");
    }

    [Test]
    public void RefineBreakpoint_NoSplitReadsInRegion_ReturnsNull()
    {
        // "no split reads" boundary (checklist row 201): an entirely empty read set,
        // and a non-empty set whose junctions fall outside the region, both yield
        // null — never a fabricated consensus (§6.1).
        RefineBreakpoint("chr1", 100, 200, Array.Empty<SplitRead>())
            .Should().BeNull("no read supports the region ⇒ null (§6.1)");

        var outsideReads = new[] { Read("p", "chr1", 50), Read("q", "chr1", 300) };
        RefineBreakpoint("chr1", 100, 200, outsideReads)
            .Should().BeNull("no junction lies in [100,200] ⇒ null (§6.1)");

        // Right chromosome filter: a junction in range but on the wrong contig is
        // ignored (chromosome compare is ordinal/case-sensitive, §3.3).
        var wrongChrom = new[] { Read("r", "chr2", 150) };
        RefineBreakpoint("chr1", 100, 200, wrongChrom)
            .Should().BeNull("a region junction on chr2 does not support chr1 (§3.3)");
    }

    [Test]
    public void FindBreakpoints_NullInput_ThrowsArgumentNullException()
    {
        // §3.3: null splitReads ⇒ ArgumentNullException, thrown eagerly (before
        // iteration) per the sibling DetectSVs pattern (§5.2).
        Action act = () => FindBreakpoints(null!, clusterTolerance: 5, minSupport: 2);
        act.Should().Throw<ArgumentNullException>();

        Action refine = () => RefineBreakpoint("chr1", 0, 10, null!);
        refine.Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region SV-BREAKPOINT-001 — BE boundary: single breakpoint / below the support gate

    [Test]
    public void FindBreakpoints_SingleSplitRead_BelowDefaultGate_YieldsNothing()
    {
        // "single breakpoint" boundary: one split read is one cluster of size 1,
        // which is < the default minSupport 2 ⇒ NOT reported. No one-element-cluster
        // crash (§6.1, INV-01).
        var reads = new[] { Read("solo", "chr1", 7777) };

        FindBreakpoints(reads, clusterTolerance: DefaultTolerance, minSupport: DefaultMinSupport)
            .Should().BeEmpty("a single read is below minSupport 2 ⇒ no breakpoint (INV-01)");
    }

    [Test]
    public void FindBreakpoints_SingleSplitRead_MinSupportOne_YieldsOneBreakpointAtItsJunction()
    {
        // With minSupport relaxed to 1 the lone read IS a valid breakpoint at its own
        // junction with support 1 — the rounded mean of a singleton is the value
        // itself (ASM-01, INV-03/INV-04). Pins that the gate, not the clustering, is
        // what suppressed it above.
        var reads = new[] { Read("solo", "chr3", 12345) };

        var breakpoints = FindBreakpoints(reads, clusterTolerance: DefaultTolerance, minSupport: 1).ToList();

        breakpoints.Should().HaveCount(1);
        breakpoints[0].Position1.Should().Be(12345, "singleton rounded mean = the junction (ASM-01)");
        breakpoints[0].SupportingReads.Should().Be(1, "support = cluster size 1 (INV-04)");
        AssertWellFormed(breakpoints, inputReadCount: 1, minSupport: 1);
    }

    [Test]
    public void FindBreakpoints_ClusterExactlyBelowGate_IsDropped_AtGateIsKept()
    {
        // The gate is "≥ minSupport" (§4.1[4]): a 2-read cluster is dropped at
        // minSupport 3 but kept at minSupport 2 — boundary on the support threshold.
        var reads = new[] { Read("a", "chr1", 400), Read("b", "chr1", 401) };

        FindBreakpoints(reads, clusterTolerance: 5, minSupport: 3)
            .Should().BeEmpty("cluster size 2 < minSupport 3 ⇒ dropped (INV-01)");

        FindBreakpoints(reads, clusterTolerance: 5, minSupport: 2)
            .Should().HaveCount(1, "cluster size 2 ≥ minSupport 2 ⇒ kept (INV-01)");
    }

    #endregion

    #region SV-BREAKPOINT-001 — BE boundary: identical junctions

    [Test]
    public void FindBreakpoints_IdenticalJunctions_CollapseToOneBreakpoint()
    {
        // "identical" boundary (checklist row 201): many reads at the SAME junction
        // are one cluster (gap 0 ≤ any tolerance ≥ 0). Reported position = that exact
        // coordinate; support = the count (INV-03/INV-04).
        var reads = Enumerable.Range(0, 6).Select(i => Read($"r{i}", "chr1", 9000)).ToArray();

        var breakpoints = FindBreakpoints(reads, clusterTolerance: 5, minSupport: 2).ToList();

        breakpoints.Should().HaveCount(1, "identical junctions collapse to one breakpoint");
        breakpoints[0].Position1.Should().Be(9000, "mean of identical values is the value (INV-03)");
        breakpoints[0].SupportingReads.Should().Be(6, "support = 6 identical reads (INV-04)");
        AssertWellFormed(breakpoints, inputReadCount: 6, minSupport: 2);
    }

    [Test]
    public void FindBreakpoints_IdenticalJunctions_ToleranceZero_StillOneCluster()
    {
        // tolerance = 0 (the 0 boundary): identical junctions still cluster because
        // the gap is exactly 0 ≤ 0, but any +1 neighbour would split. Confirms the
        // window is inclusive of an exact match (INV-02).
        var same = Enumerable.Range(0, 4).Select(i => Read($"s{i}", "chr1", 555)).ToArray();
        FindBreakpoints(same, clusterTolerance: 0, minSupport: 2)
            .Should().ContainSingle().Which.SupportingReads.Should().Be(4);

        var adjacent = new[] { Read("a", "chr1", 555), Read("b", "chr1", 556) };
        FindBreakpoints(adjacent, clusterTolerance: 0, minSupport: 1)
            .Should().HaveCount(2, "gap 1 > tolerance 0 ⇒ separate clusters (INV-02)");
    }

    [Test]
    public void RefineBreakpoint_ModalJunctionWins_TieBrokenByRoundedMean()
    {
        // §4.1[5]: the consensus is the modal junction. Junctions {100,100,100,108}
        // in region ⇒ mode 100 (3 reads) beats 108 (1 read).
        var reads = new[]
        {
            Read("a", "chr1", 100), Read("b", "chr1", 100),
            Read("c", "chr1", 100), Read("d", "chr1", 108),
        };
        RefineBreakpoint("chr1", 90, 120, reads).Should().Be(100, "mode of the region junctions (§4.1)");

        // Tie between modes {100,110} (each twice) ⇒ rounded mean = 105 (§4.1).
        var tied = new[]
        {
            Read("a", "chr1", 100), Read("b", "chr1", 100),
            Read("c", "chr1", 110), Read("d", "chr1", 110),
        };
        RefineBreakpoint("chr1", 90, 120, tied).Should().Be(105, "tie ⇒ rounded mean of modal coords (§4.1)");
    }

    #endregion

    #region SV-BREAKPOINT-001 — BE boundary: extreme coordinates (0 / -1 / MaxInt) + sweep

    [Test]
    public void FindBreakpoints_ExtremeCoordinates_NoOverflowMisMerge()
    {
        // MaxInt / 0 / negative coordinates: junctions int.MaxValue and 0 are an
        // enormous gap apart. The gap |MaxInt − 0| must register as > tolerance (a
        // naive subtraction would NOT overflow here since both are non-negative, but
        // a clipped read at int.MinValue would). Assert no mis-merge and no throw.
        var reads = new[]
        {
            Read("lo", "chrX", 0),
            Read("loDup", "chrX", 0),
            Read("hi", "chrX", int.MaxValue),
            Read("hiDup", "chrX", int.MaxValue),
        };

        var breakpoints = FindBreakpoints(reads, clusterTolerance: 5, minSupport: 2).ToList();

        breakpoints.Should().HaveCount(2, "0 and int.MaxValue are far beyond tolerance ⇒ two clusters");
        breakpoints.Select(b => b.Position1).Should().BeEquivalentTo(new[] { 0, int.MaxValue });
        breakpoints.Should().OnlyContain(b => b.SupportingReads == 2);
        AssertWellFormed(breakpoints, inputReadCount: 4, minSupport: 2);
    }

    [Test]
    public void RefineBreakpoint_RegionWithoutOrdering_DoesNotThrow_StartAboveEnd()
    {
        // A degenerate region (start > end) cannot contain any junction, so the
        // inclusive filter [start,end] is empty ⇒ null, not a throw. Boundary on the
        // region precondition (§3.1 regionStart ≤ regionEnd treated as caller duty).
        var reads = new[] { Read("a", "chr1", 150) };
        RefineBreakpoint("chr1", 200, 100, reads).Should().BeNull("empty inclusive range ⇒ null");
    }

    [Test]
    [CancelAfter(30_000)]
    public void FindBreakpoints_RandomizedBoundarySweep_NeverThrows_WellFormed()
    {
        // Randomized BE sweep over the degenerate shapes: empty / singleton sets,
        // identical junctions, tolerance 0, extreme and negative coordinates, and
        // multi-contig mixes. Every accepted result must satisfy the documented
        // invariants and the call must never throw or hang (CancelAfter 30s).
        var rng = new Random(201_001);
        for (int trial = 0; trial < 3000; trial++)
        {
            int n = rng.Next(0, 12);                       // includes empty / singleton
            int tolerance = rng.Next(0, 8);                // includes the 0 boundary
            int minSupport = rng.Next(1, 5);
            int contigs = rng.Next(1, 4);

            var reads = new List<SplitRead>(n);
            for (int i = 0; i < n; i++)
            {
                string chrom = $"chr{rng.Next(0, contigs)}";
                // Mix of small, identical-prone, negative and MaxInt-adjacent coords.
                int junction = rng.Next(0, 5) switch
                {
                    0 => 0,
                    1 => -rng.Next(0, 50),                 // negative reference coords
                    2 => int.MaxValue - rng.Next(0, 6),    // MaxInt neighbourhood
                    3 => 1000,                             // identical-collision bait
                    _ => rng.Next(0, 2000),
                };
                reads.Add(Read($"r{i}", chrom, junction));
            }

            List<Breakpoint> breakpoints = null!;
            Action act = () => breakpoints = FindBreakpoints(reads, tolerance, minSupport).ToList();
            act.Should().NotThrow($"trial {trial}: degenerate split-read set must not throw");

            AssertWellFormed(breakpoints, reads.Count, minSupport);

            // Independent recomputation of the contract from the spec (§4.1): group by
            // contig, sort junctions, split on a gap > tolerance, keep clusters of size
            // ≥ minSupport. Compare ONLY counts/support against the implementation so a
            // wrong clustering would be caught (not echoing the code's own arrays).
            int expectedCount = 0;
            var expectedSupports = new List<int>();
            foreach (var grp in reads.GroupBy(r => r.Chromosome))
            {
                var js = grp.Select(r => r.SupplementaryPosition).OrderBy(p => p).ToList();
                int size = 1;
                for (int k = 1; k <= js.Count; k++)
                {
                    bool same = k < js.Count && (long)js[k] - js[k - 1] <= tolerance;
                    if (same)
                    {
                        size++;
                    }
                    else
                    {
                        if (size >= minSupport) { expectedCount++; expectedSupports.Add(size); }
                        size = 1;
                    }
                }
            }

            breakpoints.Count.Should().Be(expectedCount,
                $"trial {trial}: emitted breakpoint count must match the spec clustering (§4.1)");
            breakpoints.Select(b => b.SupportingReads).OrderBy(s => s)
                .Should().BeEquivalentTo(expectedSupports.OrderBy(s => s),
                    $"trial {trial}: per-cluster support must match cluster sizes (INV-04)");
        }
    }

    [Test]
    [CancelAfter(30_000)]
    public void RefineBreakpoint_RandomizedSweep_NeverThrows_NullOrModalInRange()
    {
        // Sweep RefineBreakpoint over empty/degenerate regions and junction sets. The
        // result is either null (no support) or a coordinate whose support is the
        // modal count among in-region junctions (§4.1[5]).
        var rng = new Random(201_002);
        for (int trial = 0; trial < 2000; trial++)
        {
            int n = rng.Next(0, 10);
            var reads = new List<SplitRead>(n);
            for (int i = 0; i < n; i++)
                reads.Add(Read($"r{i}", "chr1", rng.Next(-20, 120)));

            int a = rng.Next(-30, 130);
            int b = rng.Next(-30, 130);
            int start = Math.Min(a, b);
            int end = Math.Max(a, b);

            int? result = null;
            Action act = () => result = RefineBreakpoint("chr1", start, end, reads);
            act.Should().NotThrow($"trial {trial}: refine must not throw");

            var inRange = reads
                .Where(r => r.SupplementaryPosition >= start && r.SupplementaryPosition <= end)
                .Select(r => r.SupplementaryPosition)
                .ToList();

            if (inRange.Count == 0)
            {
                result.Should().BeNull($"trial {trial}: no in-range junction ⇒ null (§6.1)");
            }
            else
            {
                result.Should().NotBeNull($"trial {trial}: an in-range junction ⇒ a consensus");
                int modalCount = inRange.GroupBy(p => p).Max(g => g.Count());
                var modes = inRange.GroupBy(p => p).Where(g => g.Count() == modalCount)
                    .Select(g => g.Key).ToList();
                int expected = modes.Count == 1
                    ? modes[0]
                    : (int)Math.Round(modes.Average(p => (double)p));
                result.Should().Be(expected,
                    $"trial {trial}: consensus = modal junction (tie → rounded mean) (§4.1)");
            }
        }
    }

    #endregion
}
