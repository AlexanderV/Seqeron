using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for structural-variant analysis (StructuralVariantAnalyzer): breakpoint
/// clustering, copy-number detection, and SV calling.
///
/// Test Units: SV-BREAKPOINT-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("StructuralVariant")]
public class StructuralVariantProperties
{
    #region SV-BREAKPOINT-001: R: positions valid; M: more split reads → ≥ confidence; D: deterministic

    // FindBreakpoints clusters split-read junctions within a tolerance; a cluster is reported when it
    // has ≥ minSupport reads, with Quality = min(support·15, 100) (SoftSearch/ClipCrop model).

    private static StructuralVariantAnalyzer.SplitRead Sr(string id, int junction) =>
        new(id, "chr1", PrimaryPosition: junction + 100, SupplementaryPosition: junction, ClipLength: 10, ClippedSequence: "ACGT");

    /// <summary>Generates 2..8 split reads on chr1 with junctions in [0,40] (so clusters can form).</summary>
    private static Arbitrary<StructuralVariantAnalyzer.SplitRead[]> SplitReadsArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int n = 2 + rng.Next(7);
            var reads = new StructuralVariantAnalyzer.SplitRead[n];
            for (int i = 0; i < n; i++) reads[i] = Sr($"r{i}", rng.Next(0, 40));
            return reads;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R): every reported breakpoint has valid (non-negative, single-junction) positions, meets
    /// the support threshold, and carries a quality in (0,100].
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Breakpoints_AreValid()
    {
        return Prop.ForAll(SplitReadsArbitrary(), reads =>
        {
            const int minSupport = 2;
            var bps = StructuralVariantAnalyzer.FindBreakpoints(reads, clusterTolerance: 5, minSupport: minSupport).ToList();
            return bps.All(b =>
                b.Position1 >= 0 && b.Position1 == b.Position2 &&
                b.SupportingReads >= minSupport &&
                b.Quality is > 0.0 and <= 100.0)
                .Label("a breakpoint had invalid position/support/quality");
        });
    }

    /// <summary>
    /// INV-2 (M): more supporting split reads give at least as high a breakpoint quality.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Breakpoints_MoreSupport_HigherQuality()
    {
        var few = Enumerable.Range(0, 2).Select(i => Sr($"r{i}", 100)).ToArray();
        var many = Enumerable.Range(0, 5).Select(i => Sr($"r{i}", 100)).ToArray();

        double qFew = StructuralVariantAnalyzer.FindBreakpoints(few, 5, 2).Single().Quality;
        double qMany = StructuralVariantAnalyzer.FindBreakpoints(many, 5, 2).Single().Quality;

        Assert.That(qMany, Is.GreaterThanOrEqualTo(qFew), "more split reads must not lower the quality");
        Assert.That(qMany, Is.GreaterThan(qFew), "more support should raise quality below the cap");
    }

    /// <summary>
    /// INV-3 (D): Breakpoint detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Breakpoints_AreDeterministic()
    {
        return Prop.ForAll(SplitReadsArbitrary(), reads =>
        {
            var a = StructuralVariantAnalyzer.FindBreakpoints(reads, 5, 2).Select(b => (b.Position1, b.SupportingReads)).ToList();
            var b = StructuralVariantAnalyzer.FindBreakpoints(reads, 5, 2).Select(x => (x.Position1, x.SupportingReads)).ToList();
            return a.SequenceEqual(b).Label("FindBreakpoints must be deterministic");
        });
    }

    #endregion
}
