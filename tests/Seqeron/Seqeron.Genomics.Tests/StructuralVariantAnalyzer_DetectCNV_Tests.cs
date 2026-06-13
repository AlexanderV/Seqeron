// SV-CNV-001 — Read-Depth Copy Number Variation Detection
// Evidence: docs/Evidence/SV-CNV-001-Evidence.md
// TestSpec: tests/TestSpecs/SV-CNV-001.md
// Source: Yoon S, Xuan Z, Makarov V, Ye K, Sebat J (2009). Genome Research 19(9):1586-1592,
//         doi:10.1101/gr.092981.109 (read depth proportional to copy number; windowed read counts);
//         CNVkit cnvlib/call.py (_log2_ratio_to_absolute_pure: n = r*2^v) and calling docs
//         (diploid: CN = 2 * 2^log2; log2(1/2)=-1.0, log2(3/2)=0.585), Talevich et al. (2016)
//         PLoS Comput Biol 12(4):e1004873, doi:10.1371/journal.pcbi.1004873.

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.StructuralVariantAnalyzer;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class StructuralVariantAnalyzer_DetectCNV_Tests
{
    // Reference (copy-number-neutral) read depth used as the log2 anchor throughout (Evidence dataset).
    private const double ReferenceDepth = 100.0;
    private const int Window = 4;

    private static int[] FlatWindow(int depth, int size = Window) =>
        Enumerable.Repeat(depth, size).ToArray();

    #region DetectCNV

    // M1 — Neutral: window mean RD == reference RD => log2 0 => CN 2 (CNVkit: round(2*2^0)=2; log2(2/2)=0). INV-01/INV-02.
    [Test]
    public void DetectCNV_NeutralWindow_CopyNumberTwo()
    {
        var depth = FlatWindow(100);

        var segments = DetectCNV(depth, Window, ReferenceDepth).ToList();

        Assert.That(segments, Has.Count.EqualTo(1), "One full window must yield exactly one segment.");
        Assert.Multiple(() =>
        {
            Assert.That(segments[0].LogRatio, Is.EqualTo(0.0).Within(1e-10),
                "Mean RD equal to the reference gives log2(100/100)=0 (CNVkit log2 ratio definition).");
            Assert.That(segments[0].CopyNumber, Is.EqualTo(2),
                "log2 0 maps to the diploid neutral copy number: round(2*2^0)=2 (CNVkit _log2_ratio_to_absolute_pure).");
        });
    }

    // M2 — Single-copy loss: mean RD 50, ref 100 => log2 -1.0 => CN 1 (CNVkit: "single-copy loss is log2(1/2) = -1.0"). INV-02.
    [Test]
    public void DetectCNV_HalfDepthWindow_CopyNumberOneDeletion()
    {
        var depth = FlatWindow(50);

        var segments = DetectCNV(depth, Window, ReferenceDepth).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(segments[0].LogRatio, Is.EqualTo(-1.0).Within(1e-10),
                "Half the reference depth gives log2(50/100)=log2(1/2)=-1.0 (CNVkit single-copy loss anchor).");
            Assert.That(segments[0].CopyNumber, Is.EqualTo(1),
                "log2 -1.0 maps to copy number 1: round(2*2^-1)=round(1)=1 (CNVkit), a single-copy deletion.");
        });
    }

    // M3 — Single-copy gain: mean RD 150, ref 100 => log2 0.585 => CN 3 (CNVkit: "single-copy gain ... log2(3/2)=0.585"). INV-02.
    [Test]
    public void DetectCNV_OnePointFiveDepthWindow_CopyNumberThreeDuplication()
    {
        var depth = FlatWindow(150);

        var segments = DetectCNV(depth, Window, ReferenceDepth).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(segments[0].LogRatio, Is.EqualTo(Math.Log2(1.5)).Within(1e-10),
                "1.5x the reference depth gives log2(150/100)=log2(3/2)=0.585 (CNVkit single-copy gain anchor).");
            Assert.That(segments[0].CopyNumber, Is.EqualTo(3),
                "log2 0.585 maps to copy number 3: round(2*2^0.585)=round(3.0)=3 (CNVkit), a single-copy duplication.");
        });
    }

    // M4 — Amplification: mean RD 200, ref 100 => log2 1.0 => CN 4 (CNVkit diploid: CN = 2*2^log2). INV-02.
    [Test]
    public void DetectCNV_DoubleDepthWindow_CopyNumberFour()
    {
        var depth = FlatWindow(200);

        var segments = DetectCNV(depth, Window, ReferenceDepth).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(segments[0].LogRatio, Is.EqualTo(1.0).Within(1e-10),
                "Double the reference depth gives log2(200/100)=log2(2)=1.0 (CNVkit log2 ratio definition).");
            Assert.That(segments[0].CopyNumber, Is.EqualTo(4),
                "log2 1.0 maps to copy number 4: round(2*2^1)=round(4)=4 (CNVkit diploid conversion).");
        });
    }

    // M5 — Windowing: 8 positions, window 4 => two non-overlapping windows; mean per window; ProbeCount = window size. INV-05.
    [Test]
    public void DetectCNV_EightPositionsWindowFour_TwoSegments()
    {
        // Window 0 mean = 100 (neutral), window 1 mean = 50 (loss).
        int[] depth = { 100, 100, 100, 100, 50, 50, 50, 50 };

        var segments = DetectCNV(depth, Window, ReferenceDepth).ToList();

        Assert.That(segments, Has.Count.EqualTo(2),
            "Eight positions with window size 4 form two non-overlapping full windows (Yoon et al. 2009 windowing).");
        Assert.Multiple(() =>
        {
            Assert.That(segments[0].Start, Is.EqualTo(0), "First window starts at position 0.");
            Assert.That(segments[0].End, Is.EqualTo(3), "First window is inclusive [0,3] for size 4.");
            Assert.That(segments[0].ProbeCount, Is.EqualTo(4), "ProbeCount equals the number of positions in the window.");
            Assert.That(segments[0].CopyNumber, Is.EqualTo(2), "First window (mean 100 == ref) is neutral CN 2.");
            Assert.That(segments[1].Start, Is.EqualTo(4), "Second window starts at position 4.");
            Assert.That(segments[1].CopyNumber, Is.EqualTo(1), "Second window (mean 50) is a single-copy loss CN 1.");
        });
    }

    // M6 — No baseline: reference RD defaults to the overall median of non-zero window means (Yoon overall median m). ASSUMPTION A1.
    [Test]
    public void DetectCNV_NoBaseline_UsesOverallMedianReference()
    {
        // Window means: 50, 100, 200. Overall median = 100 => median window is neutral CN 2.
        int[] depth =
        {
            50, 50, 50, 50,
            100, 100, 100, 100,
            200, 200, 200, 200,
        };

        var segments = DetectCNV(depth, Window /* referenceDepth: null */).ToList();

        Assert.That(segments, Has.Count.EqualTo(3), "Three full windows produce three segments.");
        Assert.Multiple(() =>
        {
            Assert.That(segments[1].LogRatio, Is.EqualTo(0.0).Within(1e-10),
                "With the overall median (100) as the reference, the median-depth window has log2 0 (Yoon overall-median baseline m).");
            Assert.That(segments[1].CopyNumber, Is.EqualTo(2),
                "The window at the median reference depth is copy-number neutral (CN 2).");
            Assert.That(segments[0].CopyNumber, Is.EqualTo(1),
                "The half-median window (mean 50) is a single-copy loss CN 1.");
            Assert.That(segments[2].CopyNumber, Is.EqualTo(4),
                "The double-median window (mean 200) is an amplification CN 4.");
        });
    }

    // M7 — Zero-depth window is undefined (log2(0)=-inf) => no-call, excluded from output. INV-06.
    [Test]
    public void DetectCNV_ZeroDepthWindow_ExcludedAsNoCall()
    {
        // Window 0 mean = 100 (neutral); window 1 mean = 0 (no-call).
        int[] depth = { 100, 100, 100, 100, 0, 0, 0, 0 };

        var segments = DetectCNV(depth, Window, ReferenceDepth).ToList();

        Assert.That(segments, Has.Count.EqualTo(1),
            "A zero-depth window has an undefined log2 ratio (log2(0)=-inf) and is reported as a no-call, not a finite call.");
        Assert.That(segments[0].Start, Is.EqualTo(0),
            "Only the non-zero window (positions 0-3) is emitted; the zero-depth window 4-7 is excluded.");
    }

    // S2 — Explicit baseline overrides the median: ref 50, window mean 50 => log2 0 => CN 2.
    [Test]
    public void DetectCNV_ExplicitBaseline_OverridesMedian()
    {
        var depth = FlatWindow(50);

        var segments = DetectCNV(depth, Window, referenceDepth: 50.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(segments[0].LogRatio, Is.EqualTo(0.0).Within(1e-10),
                "An explicit reference depth equal to the window mean gives log2 0, overriding the median default.");
            Assert.That(segments[0].CopyNumber, Is.EqualTo(2),
                "With the supplied reference the window is neutral CN 2 (CNVkit ratio-against-reference).");
        });
    }

    // S1 — Property: copy numbers are non-negative and non-decreasing as window mean RD increases. INV-03/INV-04.
    [Test]
    public void DetectCNV_NonNegativeDepths_CopyNumbersNonNegativeAndMonotone()
    {
        // Strictly increasing window means 25, 50, 100, 200, 400 against reference 100.
        int[] depth =
        {
            25, 25, 25, 25,
            50, 50, 50, 50,
            100, 100, 100, 100,
            200, 200, 200, 200,
            400, 400, 400, 400,
        };

        var copyNumbers = DetectCNV(depth, Window, ReferenceDepth)
            .Select(s => s.CopyNumber)
            .ToList();

        Assert.That(copyNumbers, Has.All.GreaterThanOrEqualTo(0),
            "Copy number is physically non-negative (CNVkit clamps max(0,n)). [INV-03]");
        for (int i = 1; i < copyNumbers.Count; i++)
        {
            Assert.That(copyNumbers[i], Is.GreaterThanOrEqualTo(copyNumbers[i - 1]),
                "CN = round(2*2^log2) is non-decreasing as window mean read depth increases (Yoon linear coverage<->CN). [INV-04]");
        }
    }

    // C1 — Empty input yields empty output (defined trivial behavior).
    [Test]
    public void DetectCNV_EmptyInput_ReturnsEmpty()
    {
        var segments = DetectCNV(Array.Empty<int>(), Window, ReferenceDepth).ToList();

        Assert.That(segments, Is.Empty, "No depth data means no window can be formed and no CNV can be called.");
    }

    // C2 — Null input throws ArgumentNullException (input-validation contract).
    [Test]
    public void DetectCNV_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => DetectCNV(null!, Window, ReferenceDepth).ToList(),
            "Null depth input violates the precondition and must throw ArgumentNullException.");
    }

    // C3 — Window larger than the data forms no full window => empty output. INV-05.
    [Test]
    public void DetectCNV_WindowLargerThanData_ReturnsEmpty()
    {
        int[] depth = { 100, 100, 100 }; // 3 positions, window 4 => no full window.

        var segments = DetectCNV(depth, windowSize: 4, ReferenceDepth).ToList();

        Assert.That(segments, Is.Empty,
            "A trailing partial window (3 positions < window size 4) is dropped, so no segment is emitted. [INV-05]");
    }

    // C4 — All-zero depth with no baseline: overall-median reference is 0 => no callable window => empty. INV-06.
    [Test]
    public void DetectCNV_AllZeroDepthNoBaseline_ReturnsEmpty()
    {
        int[] depth = { 0, 0, 0, 0, 0, 0, 0, 0 };

        var segments = DetectCNV(depth, Window /* referenceDepth: null */).ToList();

        Assert.That(segments, Is.Empty,
            "When every window is zero-depth the overall-median reference is 0 and no finite copy number can be called. [INV-06]");
    }

    // C2b — Non-positive window size throws ArgumentOutOfRangeException (input-validation contract).
    [Test]
    public void DetectCNV_NonPositiveWindow_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => DetectCNV(new[] { 100 }, windowSize: 0, ReferenceDepth).ToList(),
            "A window must have positive size; windowSize 0 must throw ArgumentOutOfRangeException.");
    }

    #endregion

    #region SegmentCopyNumber (delegate / segmentation variant)

    // M8 — log2 ratios [0,0,-1,-1] => two merged segments: CN 2 (len 2) then CN 1 (len 2). Same conversion as DetectCNV.
    [Test]
    public void SegmentCopyNumber_LogRatios_MergesAdjacentEqualCopyNumber()
    {
        double[] logRatios = { 0.0, 0.0, -1.0, -1.0 };

        var segments = SegmentCopyNumber(logRatios).ToList();

        Assert.That(segments, Has.Count.EqualTo(2),
            "Adjacent windows of equal copy number merge: a CN-2 run then a CN-1 run.");
        Assert.Multiple(() =>
        {
            Assert.That(segments[0].CopyNumber, Is.EqualTo(2), "log2 0 maps to CN 2 (CNVkit round(2*2^0)=2).");
            Assert.That(segments[0].ProbeCount, Is.EqualTo(2), "The neutral run merges its two windows.");
            Assert.That(segments[0].Start, Is.EqualTo(0), "First merged run starts at window index 0.");
            Assert.That(segments[0].End, Is.EqualTo(1), "First merged run ends at window index 1.");
            Assert.That(segments[1].CopyNumber, Is.EqualTo(1), "log2 -1.0 maps to CN 1 (CNVkit round(2*2^-1)=1).");
            Assert.That(segments[1].ProbeCount, Is.EqualTo(2), "The loss run merges its two windows.");
        });
    }

    // M8b — A NaN log2 ratio is a no-call: it breaks the current run and is not emitted as a segment.
    [Test]
    public void SegmentCopyNumber_NaNLogRatio_BreaksRunAndIsNotEmitted()
    {
        // CN-2 run, then a NaN no-call, then a CN-1 run.
        double[] logRatios = { 0.0, 0.0, double.NaN, -1.0, -1.0 };

        var segments = SegmentCopyNumber(logRatios).ToList();

        Assert.That(segments, Has.Count.EqualTo(2),
            "A NaN log2 ratio is a no-call window; it terminates the neutral run and starts a fresh loss run, yielding two segments.");
        Assert.Multiple(() =>
        {
            Assert.That(segments[0].CopyNumber, Is.EqualTo(2), "First run is the CN-2 neutral run before the no-call.");
            Assert.That(segments[0].End, Is.EqualTo(1), "The neutral run ends at index 1, before the NaN at index 2.");
            Assert.That(segments[1].CopyNumber, Is.EqualTo(1), "Second run is the CN-1 loss run after the no-call.");
            Assert.That(segments[1].Start, Is.EqualTo(3), "The loss run starts at index 3, after the NaN at index 2.");
        });
    }

    // C2c — Null log2 ratios throws ArgumentNullException (delegate smoke / input-validation).
    [Test]
    public void SegmentCopyNumber_NullInput_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => SegmentCopyNumber((System.Collections.Generic.IEnumerable<double>)null!).ToList(),
            "Null log2-ratio input violates the precondition and must throw ArgumentNullException.");
    }

    #endregion
}
