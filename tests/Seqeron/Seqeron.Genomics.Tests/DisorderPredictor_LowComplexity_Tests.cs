// DISORDER-LC-001 — Low-Complexity Region Detection in Protein Sequences (SEG)
// Evidence: docs/Evidence/DISORDER-LC-001-Evidence.md
// TestSpec: tests/TestSpecs/DISORDER-LC-001.md
// Source: Wootton J.C., Federhen S. (1993). Comput. Chem. 17(2):149-163;
//         NCBI BLAST blast_seg.c (kSegWindow=12, kSegLocut=2.2, kSegHicut=2.5);
//         GCG/NCBI seg help (complexity = Shannon entropy, bits/residue, max log2(20)=4.322).
using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// DISORDER-LC-001: Low-Complexity Region Detection (SEG).
/// Complexity H = -Σ pᵢ·log₂(pᵢ) over the window composition (bits/residue).
/// Defaults: triggerWindow W=12, K1=2.2, K2=2.5.
/// </summary>
[TestFixture]
public class DisorderPredictor_LowComplexity_Tests
{
    #region PredictLowComplexityRegions — MUST

    // M1 — Homopolymer ≥ W: each 12-window has H=0 ≤ K1=2.2 → one segment over whole run.
    [Test]
    public void PredictLowComplexityRegions_PolyQ26_SingleRegionWholeSequence()
    {
        string polyQ = new string('Q', 26);

        var regions = DisorderPredictor.PredictLowComplexityRegions(polyQ).ToList();

        Assert.That(regions, Has.Count.EqualTo(1),
            "26-residue homopolymer: every window has entropy 0 ≤ K1, so exactly one merged segment is expected");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "segment must start at index 0");
            Assert.That(regions[0].End, Is.EqualTo(25), "segment must span to the last index (0-based inclusive)");
            Assert.That(regions[0].Type, Is.EqualTo("Q-rich"), "single dominant residue Q (>50%) → 'Q-rich'");
        });
    }

    // M2 — Max-complexity sequence: each 12-window has ≥12 distinct residues,
    // H ≥ log₂(12) ≈ 3.585 bits > K2=2.5 → no segment triggers or extends.
    [Test]
    public void PredictLowComplexityRegions_AllTwentyAminoAcidsTwice_NoRegions()
    {
        string complex = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"; // 40 AA

        var regions = DisorderPredictor.PredictLowComplexityRegions(complex).ToList();

        Assert.That(regions, Is.Empty,
            "every 12-window has ≥12 distinct residues (H ≈ 3.585 bits > K2=2.5), so no low-complexity region exists");
    }

    // M3 — Four residue types × 3 each in a 12-window: H = 2.0 bits ≤ K1=2.2 → triggers.
    [Test]
    public void PredictLowComplexityRegions_FourTypesEntropy2_TriggersOneRegion()
    {
        string fourTypes = "AAABBBCCCDDD"; // L=12, H = -4·(0.25·log₂0.25) = 2.0 bits

        var regions = DisorderPredictor.PredictLowComplexityRegions(fourTypes).ToList();

        Assert.That(regions, Has.Count.EqualTo(1),
            "window entropy 2.0 bits ≤ K1=2.2 → exactly one triggered segment");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "segment starts at 0");
            Assert.That(regions[0].End, Is.EqualTo(11), "segment ends at last index 11");
        });
    }

    // M4 — Same input, strict K1=0.5: H=2.0 > 0.5 → no trigger.
    [Test]
    public void PredictLowComplexityRegions_FourTypesEntropy2_StrictTrigger_NoRegions()
    {
        string fourTypes = "AAABBBCCCDDD"; // H = 2.0 bits

        var regions = DisorderPredictor.PredictLowComplexityRegions(
            fourTypes, triggerThreshold: 0.5).ToList();

        Assert.That(regions, Is.Empty,
            "with K1=0.5, window entropy 2.0 bits > 0.5 → nothing triggers");
    }

    // M5 — Two-residue block 12×A + 12×L: each window has at most two types,
    // H ≤ log₂(2)=1.0 bits ≤ K1=2.2 → single merged segment over [0,23].
    [Test]
    public void PredictLowComplexityRegions_DipeptideBlock_SingleMergedRegion()
    {
        string sequence = new string('A', 12) + new string('L', 12); // 24 AA

        var regions = DisorderPredictor.PredictLowComplexityRegions(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(1),
            "every window has ≤2 residue types (H ≤ 1.0 bit ≤ K1), so spans merge into one segment");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "merged segment starts at 0");
            Assert.That(regions[0].End, Is.EqualTo(23), "merged segment ends at last index 23");
        });
    }

    // M6 — Two homopolymer runs separated by a high-complexity spacer (each pure-spacer
    // 12-window has ≥12 distinct residues, H ≈ 3.585 > K2) → exactly two separate segments
    // with a high-complexity gap between them. Boundaries (0,34) and (67,99) are the
    // SEG-defined result: junction windows that still contain the homopolymer trigger (≤K1)
    // and the segments then extend over neighbouring residues while segment H ≤ K2=2.5.
    [Test]
    public void PredictLowComplexityRegions_TwoRunsWithHighComplexitySpacer_TwoRegions()
    {
        string polyQ = new string('Q', 20);
        string spacer = string.Concat(Enumerable.Repeat("ACDEFGHIKLMNPQRSTVWY", 3)); // 60 AA, max complexity
        string polyA = new string('A', 20);
        string sequence = polyQ + spacer + polyA; // length 100

        var regions = DisorderPredictor.PredictLowComplexityRegions(sequence)
            .OrderBy(r => r.Start).ToList();

        Assert.That(regions, Has.Count.EqualTo(2),
            "two homopolymer runs separated by a max-complexity spacer (H > K2) → exactly two segments");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "first segment starts at 0 (poly-Q anchor)");
            Assert.That(regions[0].End, Is.EqualTo(34), "first segment extends to 34 while segment H ≤ K2=2.5");
            Assert.That(regions[0].Type, Is.EqualTo("Q-rich"), "Q is 20/35 ≈ 0.571 > 0.5 → Q-rich");
            Assert.That(regions[1].Start, Is.EqualTo(67), "second segment starts at 67 (extends left while H ≤ K2)");
            Assert.That(regions[1].End, Is.EqualTo(99), "second segment ends at last index 99 (poly-A anchor)");
            Assert.That(regions[1].Type, Is.EqualTo("A-rich"), "A is 20/33 ≈ 0.606 > 0.5 → A-rich");
            Assert.That(regions[0].End, Is.LessThan(regions[1].Start),
                "the two segments are separated by a high-complexity gap (not merged)");
        });
    }

    // M7 — Dominant-residue (>50%) label.
    [Test]
    public void PredictLowComplexityRegions_PolyA_TypeIsARich()
    {
        string polyA = new string('A', 20);

        var regions = DisorderPredictor.PredictLowComplexityRegions(polyA).ToList();

        Assert.That(regions, Has.Count.EqualTo(1), "homopolymer → one segment");
        Assert.That(regions[0].Type, Is.EqualTo("A-rich"),
            "single residue A at 100% (>50%) → 'A-rich' label");
    }

    #endregion

    #region PredictLowComplexityRegions — SHOULD / Corner

    // S1 — minLength filter removes a segment shorter than the threshold.
    [Test]
    public void PredictLowComplexityRegions_MinLengthAboveSegmentLength_FiltersRegion()
    {
        string seq = new string('Q', 12); // segment length 12

        var regions = DisorderPredictor.PredictLowComplexityRegions(seq, minLength: 15).ToList();

        Assert.That(regions, Is.Empty,
            "segment length 12 < minLength 15 → filtered out");
    }

    // S2 — Sequence shorter than W: no full trigger window.
    [Test]
    public void PredictLowComplexityRegions_ShorterThanWindow_Empty()
    {
        string seq = "QQQQQQ"; // 6 < W=12

        var regions = DisorderPredictor.PredictLowComplexityRegions(seq).ToList();

        Assert.That(regions, Is.Empty,
            "sequence shorter than the trigger window has no full window → empty");
    }

    // S3 — Exactly window length: one window at index 0.
    [Test]
    public void PredictLowComplexityRegions_ExactlyWindowLengthHomopolymer_OneRegion()
    {
        string exact = new string('Q', 12);

        var regions = DisorderPredictor.PredictLowComplexityRegions(exact).ToList();

        Assert.That(regions, Has.Count.EqualTo(1), "one full window at index 0 triggers");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "starts at 0");
            Assert.That(regions[0].End, Is.EqualTo(11), "ends at 11 (window length 12, 0-based)");
            Assert.That(regions[0].Type, Is.EqualTo("Q-rich"), "Q-rich");
        });
    }

    // S4 — Case-insensitivity: lowercase and uppercase give identical segments.
    [Test]
    public void PredictLowComplexityRegions_LowercaseInput_SameAsUppercase()
    {
        var lower = DisorderPredictor.PredictLowComplexityRegions(new string('q', 26)).ToList();
        var upper = DisorderPredictor.PredictLowComplexityRegions(new string('Q', 26)).ToList();

        Assert.That(lower, Has.Count.EqualTo(upper.Count), "same region count regardless of case");
        Assert.Multiple(() =>
        {
            Assert.That(lower[0].Start, Is.EqualTo(upper[0].Start), "same start");
            Assert.That(lower[0].End, Is.EqualTo(upper[0].End), "same end");
            Assert.That(lower[0].Type, Is.EqualTo(upper[0].Type), "same type (upper-cased)");
        });
    }

    // C1 — Empty string.
    [Test]
    public void PredictLowComplexityRegions_EmptyString_Empty()
    {
        var regions = DisorderPredictor.PredictLowComplexityRegions("").ToList();
        Assert.That(regions, Is.Empty, "empty input (length 0 < W) → empty");
    }

    // Null input → ArgumentNullException on enumeration.
    [Test]
    public void PredictLowComplexityRegions_Null_Throws()
    {
        Assert.That(() => DisorderPredictor.PredictLowComplexityRegions(null!).ToList(),
            NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
            "null sequence must raise ArgumentNullException");
    }

    #endregion

    #region Invariants

    // C2 / INV-01 — boundaries are valid for any input.
    [TestCase("QQQQQQQQQQQQQQQQQQQQQQQQQQ")]
    [TestCase("AAABBBCCCDDDAAABBBCCCDDD")]
    [TestCase("QQQQQQQQQQQQACDEFGHIKLMNQQQQQQQQQQQQQ")]
    public void PredictLowComplexityRegions_AnyInput_ValidBoundaries(string sequence)
    {
        var regions = DisorderPredictor.PredictLowComplexityRegions(sequence).ToList();

        foreach (var r in regions)
        {
            Assert.Multiple(() =>
            {
                Assert.That(r.Start, Is.GreaterThanOrEqualTo(0), "INV-01: Start ≥ 0");
                Assert.That(r.End, Is.LessThan(sequence.Length), "INV-01: End < sequence length");
                Assert.That(r.Start, Is.LessThanOrEqualTo(r.End), "INV-01: Start ≤ End");
                Assert.That(r.Type, Is.Not.Null.And.Not.Empty, "Type label must be non-empty");
            });
        }
    }

    // C3 / INV-02 — segments are non-overlapping and non-adjacent.
    [TestCase("QQQQQQQQQQQQQQQQQQQQQQQQQQ")]
    [TestCase("AAAAALLLLLAAAAALLLLLAAAAALLLLLL")]
    public void PredictLowComplexityRegions_AnyInput_NoOverlap(string sequence)
    {
        var regions = DisorderPredictor.PredictLowComplexityRegions(sequence)
            .OrderBy(r => r.Start).ToList();

        for (int i = 1; i < regions.Count; i++)
        {
            Assert.That(regions[i].Start, Is.GreaterThan(regions[i - 1].End),
                $"INV-02: region {i} must not overlap region {i - 1}");
        }
    }

    #endregion
}
