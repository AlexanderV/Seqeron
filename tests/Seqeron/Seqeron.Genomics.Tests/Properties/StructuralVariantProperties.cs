using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;
using Seqeron.Genomics.Annotation;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for structural-variant analysis (StructuralVariantAnalyzer): breakpoint
/// clustering, copy-number detection, and SV calling.
///
/// Test Units: SV-BREAKPOINT-001, SV-CNV-001
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

    #region SV-CNV-001: R: copy number ≥ 0; M: higher coverage ratio → higher CN; D: deterministic

    // DetectCNV computes per-window log2(meanDepth/reference) → CN = round(2·2^log2) (read-depth CNV;
    // Yoon et al. 2009; CNVkit).

    /// <summary>Generates 10..40 per-base depths in [0,100].</summary>
    private static Arbitrary<int[]> DepthArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int n = 10 + rng.Next(31);
            var d = new int[n];
            for (int i = 0; i < n; i++) d[i] = rng.Next(0, 101);
            return d;
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R): every called segment has a non-negative copy number, ordered coordinates within the
    /// data, and a finite log-ratio.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Cnv_SegmentsAreValid()
    {
        return Prop.ForAll(DepthArbitrary(), depth =>
        {
            var segs = StructuralVariantAnalyzer.DetectCNV(depth, windowSize: 5).ToList();
            return segs.All(s =>
                s.CopyNumber >= 0 &&
                s.Start <= s.End && s.Start >= 0 && s.End < depth.Length &&
                double.IsFinite(s.LogRatio))
                .Label("a copy-number segment was invalid");
        });
    }

    /// <summary>
    /// INV-2 (M): with a fixed reference, a window of higher mean depth yields a higher copy number.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Cnv_HigherDepth_HigherCopyNumber()
    {
        // Window 0 at depth 10 (= reference → CN 2); window 1 at depth 20 (log2 1 → CN 4).
        var depth = Enumerable.Repeat(10, 5).Concat(Enumerable.Repeat(20, 5)).ToArray();
        var segs = StructuralVariantAnalyzer.DetectCNV(depth, windowSize: 5, referenceDepth: 10).ToList();

        Assert.That(segs, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            Assert.That(segs[0].CopyNumber, Is.EqualTo(2), "depth == reference → neutral CN 2");
            Assert.That(segs[1].CopyNumber, Is.GreaterThan(segs[0].CopyNumber), "double depth → higher CN");
        });
    }

    /// <summary>
    /// INV-3 (D): CNV detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Cnv_IsDeterministic()
    {
        return Prop.ForAll(DepthArbitrary(), depth =>
        {
            var a = StructuralVariantAnalyzer.DetectCNV(depth, 5).Select(s => (s.Start, s.CopyNumber)).ToList();
            var b = StructuralVariantAnalyzer.DetectCNV(depth, 5).Select(s => (s.Start, s.CopyNumber)).ToList();
            return a.SequenceEqual(b).Label("DetectCNV must be deterministic");
        });
    }

    #endregion
}
