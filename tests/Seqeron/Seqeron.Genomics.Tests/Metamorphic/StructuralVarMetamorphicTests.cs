using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the StructuralVar area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SV-BREAKPOINT-001 — split-read breakpoint detection (StructuralVar).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 201.
///
/// API under test (StructuralVariantAnalyzer.FindBreakpoints):
///   Clusters split-read junctions within a tolerance window and reports a breakpoint per cluster,
///   with support = cluster size and quality growing with support.
///
/// Relations (derived from the clustering/support model, NOT from output):
///   • SHIFT (prepend flank shifts breakpoints): shifting every junction coordinate shifts the
///          reported breakpoint position by the same offset.
///   • MON  (more split reads ⇒ ≥ confidence): a cluster with more supporting reads has a
///          greater-or-equal quality score.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class StructuralVarMetamorphicTests
{
    private static StructuralVariantAnalyzer.SplitRead Sr(string id, int junction) =>
        new(id, "chr1", 500, junction, 20, "ACGTACGTAC");

    #region SV-BREAKPOINT-001 MON — more split reads raise the confidence

    [Test]
    [Description("MON: a breakpoint cluster with more supporting split reads has a greater-or-equal quality score.")]
    public void Breakpoints_MoreSplitReads_HigherConfidence()
    {
        double previousQuality = -1;
        int previousSupport = -1;
        foreach (int k in new[] { 2, 3, 5 })
        {
            var reads = Enumerable.Range(0, k).Select(i => Sr($"r{i}", 1000)).ToList();
            var bp = StructuralVariantAnalyzer.FindBreakpoints(reads).Single();

            bp.SupportingReads.Should().Be(k, because: "all reads cluster at the same junction");
            bp.SupportingReads.Should().BeGreaterThan(previousSupport, because: "more reads were added");
            bp.Quality.Should().BeGreaterThanOrEqualTo(previousQuality, because: "more supporting reads cannot lower the breakpoint quality");
            previousQuality = bp.Quality;
            previousSupport = bp.SupportingReads;
        }
    }

    #endregion

    #region SV-BREAKPOINT-001 SHIFT — a coordinate shift shifts the breakpoint

    [Test]
    [Description("SHIFT: shifting every split-read junction coordinate by an offset shifts the reported breakpoint position by the same offset.")]
    public void Breakpoints_CoordinateShift_ShiftsPosition()
    {
        var reads = new[] { Sr("a", 1000), Sr("b", 1001), Sr("c", 1002) };
        int originalPosition = StructuralVariantAnalyzer.FindBreakpoints(reads).Single().Position1;

        foreach (int offset in new[] { 1000, 100000 })
        {
            var shifted = reads.Select(r => Sr(r.ReadId, r.SupplementaryPosition + offset)).ToList();
            StructuralVariantAnalyzer.FindBreakpoints(shifted).Single().Position1
                .Should().Be(originalPosition + offset,
                    because: $"shifting every junction by {offset} shifts the consensus breakpoint by {offset}");
        }
    }

    #endregion
}
