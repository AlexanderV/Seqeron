// PROTMOTIF-LC-001 — Low-Complexity Region Detection (SEG, Wootton & Federhen 1993)
// Evidence: docs/Evidence/PROTMOTIF-LC-001-Evidence.md
// TestSpec: tests/TestSpecs/PROTMOTIF-LC-001.md
// Source:   Wootton JC, Federhen S (1993) Comput. Chem. 17(2):149-163, eq.(3);
//           NCBI ncbi-seg man page (W=12, K1=2.2, K2=2.5, bits/residue, max log2(20)=4.322);
//           NCBI blast_seg.c (s_Entropy); SeqComplex `ce` (-Σ p·log2 p).

using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Canonical tests for PROTMOTIF-LC-001. SEG local complexity is the Shannon entropy of the window
/// composition in bits/residue, K = -Σ pᵢ·log₂(pᵢ). Expected values are computed independently from
/// the formula (not from the implementation output): e.g. an 11A/1B window = 0.413817 bits.
/// </summary>
[TestFixture]
[Category("PROTMOTIF-LC-001")]
public class ProteinMotifFinder_FindLowComplexityRegions_Tests
{
    private const double Log2Of20 = 4.321928094887363; // max amino-acid complexity (man page)

    // SEG complexity of a window, computed via the public API (single-window sequence).
    // A sequence of exactly W residues yields one window; with K2 large enough it is always emitted
    // (we pass triggerComplexity = extensionComplexity = Log2Of20 so any window triggers), letting us
    // read the per-window complexity back as the region's Complexity field.
    private static double WindowComplexity(string window)
    {
        var regions = ProteinMotifFinder
            .FindLowComplexityRegions(window, windowSize: window.Length,
                triggerComplexity: Log2Of20, extensionComplexity: Log2Of20)
            .ToList();
        Assert.That(regions, Has.Count.EqualTo(1), "A single full window with K1=K2=log2(20) must always be emitted.");
        return regions[0].Complexity;
    }

    #region CalculateSegComplexity (via single-window region complexity)

    // M1 — Homopolymer: 12×A -> p=1 -> K = -1·log2(1) = 0.0 bits exactly.
    [Test]
    public void FindLowComplexityRegions_HomopolymerWindow_ComplexityIsZero()
    {
        double k = WindowComplexity(new string('A', 12));

        Assert.That(k, Is.EqualTo(0.0).Within(1e-10),
            "A homopolymer window has a single symbol (p=1), so Shannon entropy K = -1·log2(1) = 0 bits.");
    }

    // M2 — 11×A + 1×B, L=12: K = -[(11/12)log2(11/12) + (1/12)log2(1/12)] = 0.413817 bits.
    [Test]
    public void FindLowComplexityRegions_ElevenToOneWindow_ComplexityMatchesEntropy()
    {
        double k = WindowComplexity(new string('A', 11) + "B");

        Assert.That(k, Is.EqualTo(0.413817).Within(1e-6),
            "K = -Σ pᵢ·log2(pᵢ) for {11/12, 1/12} = 0.413817 bits/residue (independently computed).");
    }

    // M3 — 6×A + 6×B, L=12: K = -2·(1/2)log2(1/2) = 1.0 bit exactly.
    [Test]
    public void FindLowComplexityRegions_EqualSplitWindow_ComplexityIsOneBit()
    {
        double k = WindowComplexity(new string('A', 6) + new string('B', 6));

        Assert.That(k, Is.EqualTo(1.0).Within(1e-10),
            "Two equiprobable symbols give K = -2·(1/2)log2(1/2) = 1 bit/residue.");
    }

    // M4 — 10×A + 2×B, L=12: K = 0.650022 bits.
    [Test]
    public void FindLowComplexityRegions_TenToTwoWindow_ComplexityMatchesEntropy()
    {
        double k = WindowComplexity(new string('A', 10) + "BB");

        Assert.That(k, Is.EqualTo(0.650022).Within(1e-6),
            "K = -Σ pᵢ·log2(pᵢ) for {10/12, 2/12} = 0.650022 bits/residue.");
    }

    // M5 — 12 distinct residues, L=12: K = log2(12) = 3.584963 bits.
    [Test]
    public void FindLowComplexityRegions_MaxDiversityWindow_ComplexityIsLog2Of12()
    {
        double k = WindowComplexity("ACDEFGHIKLMN"); // 12 distinct amino acids

        Assert.That(k, Is.EqualTo(3.584962500721156).Within(1e-9),
            "12 distinct residues, each p=1/12, give K = log2(12) = 3.584963 bits/residue.");
    }

    #endregion

    #region FindLowComplexityRegions (SEG two-pass detection)

    // M6 — Poly-Q tract embedded in diverse flanks: detected as exactly one region, complexity 0.
    [Test]
    public void FindLowComplexityRegions_PolyQTract_DetectedAsSingleRegion()
    {
        // 12-residue diverse flanks (high entropy) ensure flank windows exceed K2 and are not merged.
        string sequence = "ACDEFGHIKLMN" + new string('Q', 20) + "NMLKIHGFEDCA";

        var regions = ProteinMotifFinder.FindLowComplexityRegions(sequence).ToList();

        // Boundaries (independently computed from K = -Σ p·log2 p per window, W=12, K1=2.2, K2=2.5):
        // 44-residue sequence = flank[0..11] + Q[12..31] + flank[32..43]. Window K-values rise from the
        // left flank, dip to 0 inside the Q core, then rise again. The maximal run of windows with K ≤ K2
        // that contains a trigger window (K ≤ K1) is windows [6..26]; residue span = [6, 26+12-1] = [6,37].
        Assert.That(regions, Has.Count.EqualTo(1),
            "A single poly-Q tract between diverse flanks must yield exactly one low-complexity region.");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(6),
                "Run of windows ≤K2 containing a trigger begins at window 6 (residue 6).");
            Assert.That(regions[0].End, Is.EqualTo(37),
                "Run ends at window 26; span = 26+12-1 = 37 (0-based inclusive).");
            Assert.That(regions[0].Complexity, Is.EqualTo(0.0).Within(1e-10),
                "The poly-Q core contains homopolymer windows, so the region's minimum complexity is 0 bits.");
        });
    }

    // M6b — K1 trigger rule: a window with K1 < K ≤ K2 lies in the extension band but never triggers,
    // so on its own it must NOT be reported. "AAABBBCCDDEE" (counts 3,3,2,2,2) has
    // K = -Σ p·log2 p = 2.292481 bits/residue, which is > K1=2.2 and ≤ K2=2.5 (independently computed).
    // A correct SEG emits a region only if at least one window has K ≤ K1; this one has none -> empty.
    // This locks the trigger semantics: an implementation that emitted any run ≤ K2 would fail here.
    [Test]
    public void FindLowComplexityRegions_WindowAboveK1WithinK2_NotTriggered_ReturnsEmpty()
    {
        // Single 12-residue window, K=2.2925: extension-band but never a trigger.
        var regions = ProteinMotifFinder.FindLowComplexityRegions("AAABBBCCDDEE").ToList();

        Assert.That(regions, Is.Empty,
            "A window with K1(2.2) < K(2.2925) ≤ K2(2.5) is in the extension band but never triggers; "
            + "with no window ≤ K1 the run must not be emitted (SEG two-pass trigger rule).");
    }

    // M7 — Fully diverse protein: every window has K = log2(20+) > K2 = 2.5 -> no regions.
    [Test]
    public void FindLowComplexityRegions_DiverseProtein_ReturnsNoRegions()
    {
        string sequence = "ACDEFGHIKLMNPQRSTVWY"; // all 20 distinct amino acids

        var regions = ProteinMotifFinder.FindLowComplexityRegions(sequence).ToList();

        Assert.That(regions, Is.Empty,
            "Every window of 20 distinct residues has K ≈ log2(20)=4.32 > K2=2.5, so no region is reported.");
    }

    // M8 — Defaults: a 12-residue poly-A homopolymer (K=0 ≤ K1=2.2) inside diverse flanks is detected
    //      with default W=12, K1=2.2, K2=2.5, confirming the SEG default parameters are in force.
    [Test]
    public void FindLowComplexityRegions_DefaultParameters_DetectHomopolymerWithSegDefaults()
    {
        string sequence = "ACDEFGHIKLMN" + new string('A', 12) + "NMLKIHGFEDCA";

        var withDefaults = ProteinMotifFinder.FindLowComplexityRegions(sequence).ToList();
        var withExplicit = ProteinMotifFinder
            .FindLowComplexityRegions(sequence, windowSize: 12, triggerComplexity: 2.2, extensionComplexity: 2.5)
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(withDefaults, Has.Count.EqualTo(1),
                "Default call must detect the poly-A homopolymer (K=0 ≤ K1=2.2).");
            Assert.That(withDefaults.Select(r => (r.Start, r.End, r.Complexity)),
                Is.EqualTo(withExplicit.Select(r => (r.Start, r.End, r.Complexity))),
                "Defaults must equal the explicit SEG defaults W=12, K1=2.2, K2=2.5.");
        });
    }

    // M9 — Two homopolymer tracts (poly-G, poly-S) separated by a 12-residue all-distinct spacer.
    // The fully-diverse spacer window (start 14, K=log2(12)=3.585 > K2=2.5) breaks the run into two
    // regions. Boundaries derived from K = -Σ p·log2 p per window (W=12, K1=2.2, K2=2.5):
    //   poly-G run = windows [0..9] (window 9 K=2.2925 ≤ K2; window 10 K=2.6175 > K2) -> span [0,20];
    //   poly-S run = windows [20..28] (window 20 K=2.2925 ≤ K2) -> span [20,39].
    [Test]
    public void FindLowComplexityRegions_TwoTractsWithDiverseSpacer_ReturnsTwoRegions()
    {
        string spacer = "ACDEFGHIKLMN"; // 12 distinct residues; central window K=log2(12) > K2
        string sequence = new string('G', 14) + spacer + new string('S', 14); // length 40

        var regions = ProteinMotifFinder.FindLowComplexityRegions(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(2),
            "The fully-diverse central spacer window (K=3.585 > K2=2.5) splits the two tracts into two regions.");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0),
                "The poly-G region begins at residue 0 (homopolymer window K=0 triggers).");
            Assert.That(regions[0].End, Is.EqualTo(20),
                "Poly-G run = windows [0..9]; window 9 (11 G + 1 spacer) has K=2.29 ≤ K2, window 10 K=2.62 > K2. Span = 9+12-1 = 20.");
            Assert.That(regions[1].Start, Is.EqualTo(20),
                "Poly-S run starts at window 20 (K=2.29 ≤ K2), giving residue start 20.");
            Assert.That(regions[1].End, Is.EqualTo(39),
                "Poly-S run = windows [20..28]; span = 28+12-1 = 39 (last residue).");
            Assert.That(regions[0].Complexity, Is.EqualTo(0.0).Within(1e-10),
                "The poly-G run contains homopolymer windows, so its minimum complexity is 0.");
            Assert.That(regions[1].Complexity, Is.EqualTo(0.0).Within(1e-10),
                "The poly-S run contains homopolymer windows, so its minimum complexity is 0.");
        });
    }

    #endregion

    #region Edge cases and invariants

    // S1 — Sequence shorter than the window: no complete trigger window -> empty.
    [Test]
    public void FindLowComplexityRegions_SequenceShorterThanWindow_ReturnsEmpty()
    {
        var regions = ProteinMotifFinder.FindLowComplexityRegions("AAAAA", windowSize: 12).ToList();

        Assert.That(regions, Is.Empty,
            "A sequence shorter than the window contains no complete trigger window (NCBI SEG man page).");
    }

    // S1b — null and empty inputs -> empty.
    [Test]
    public void FindLowComplexityRegions_NullOrEmpty_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ProteinMotifFinder.FindLowComplexityRegions(null!).ToList(), Is.Empty,
                "Null sequence yields no regions.");
            Assert.That(ProteinMotifFinder.FindLowComplexityRegions(string.Empty).ToList(), Is.Empty,
                "Empty sequence yields no regions.");
        });
    }

    // S1c — Invalid window size throws.
    [Test]
    public void FindLowComplexityRegions_NonPositiveWindow_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => ProteinMotifFinder.FindLowComplexityRegions("ACDEFGHIKLMN", windowSize: 0).ToList(),
            "Window size must be positive.");
    }

    // S2 — Region boundaries stay within the sequence and are well-ordered.
    [Test]
    public void FindLowComplexityRegions_ReportedRegion_WithinBounds()
    {
        string sequence = "ACDEFGHIKLMN" + new string('P', 20) + "NMLKIHGFEDCA";

        var regions = ProteinMotifFinder.FindLowComplexityRegions(sequence).ToList();

        Assert.That(regions, Is.Not.Empty, "A poly-P tract must produce a region.");
        Assert.Multiple(() =>
        {
            foreach (var r in regions)
            {
                Assert.That(r.Start, Is.GreaterThanOrEqualTo(0).And.LessThanOrEqualTo(r.End),
                    "Start is non-negative and not after End.");
                Assert.That(r.End, Is.LessThan(sequence.Length),
                    "End is a valid 0-based inclusive index within the sequence.");
            }
        });
    }

    // S3 — Case-insensitive: lowercase poly tract is detected like uppercase.
    [Test]
    public void FindLowComplexityRegions_LowercaseInput_DetectedCaseInsensitively()
    {
        string upper = "ACDEFGHIKLMN" + new string('Q', 20) + "NMLKIHGFEDCA";
        string lower = upper.ToLowerInvariant();

        var upperRegions = ProteinMotifFinder.FindLowComplexityRegions(upper).ToList();
        var lowerRegions = ProteinMotifFinder.FindLowComplexityRegions(lower).ToList();

        Assert.That(lowerRegions.Select(r => (r.Start, r.End)),
            Is.EqualTo(upperRegions.Select(r => (r.Start, r.End))),
            "Lowercase input must be detected identically to uppercase (case-insensitive).");
    }

    // C1 — INV-01: every reported region complexity is within [0, log2(20)].
    [Test]
    public void FindLowComplexityRegions_RegionComplexity_WithinBounds()
    {
        string sequence = "ACDEFGHIKLMN" + new string('A', 10) + "BB" + new string('N', 14);

        var regions = ProteinMotifFinder.FindLowComplexityRegions(sequence).ToList();

        Assert.That(regions, Is.Not.Empty, "Biased tracts must produce regions.");
        Assert.Multiple(() =>
        {
            foreach (var r in regions)
                Assert.That(r.Complexity, Is.GreaterThanOrEqualTo(0.0).And.LessThanOrEqualTo(Log2Of20),
                    "SEG complexity is bounded by [0, log2(20)] for amino-acid windows.");
        });
    }

    // C2 — INV-05: deterministic, identical results across repeated runs.
    [Test]
    public void FindLowComplexityRegions_RepeatedRuns_AreDeterministic()
    {
        string sequence = new string('G', 14) + "ACDEFGHIKLMN" + new string('S', 14);

        var run1 = ProteinMotifFinder.FindLowComplexityRegions(sequence).ToList();
        var run2 = ProteinMotifFinder.FindLowComplexityRegions(sequence).ToList();

        Assert.That(run2.Select(r => (r.Start, r.End, r.Complexity)),
            Is.EqualTo(run1.Select(r => (r.Start, r.End, r.Complexity))),
            "Two runs on the same input must yield identical regions (deterministic).");
    }

    #endregion
}
