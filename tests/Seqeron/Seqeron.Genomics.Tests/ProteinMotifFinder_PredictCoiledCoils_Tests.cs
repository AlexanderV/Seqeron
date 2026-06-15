// PROTMOTIF-CC-001 — Coiled-Coil Prediction (heptad-repeat a/d hydrophobic-core detection)
// Evidence: docs/Evidence/PROTMOTIF-CC-001-Evidence.md
// TestSpec: tests/TestSpecs/PROTMOTIF-CC-001.md
// Source:   Mason JM, Arndt KM (2004) ChemBioChem 5(2):170-176; Lupas A, Van Dyke M, Stock J (1991)
//           Science 252:1162-1164; Wikipedia "Coiled coil" (a/d occupied by Ile, Leu, or Val).

using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for PROTMOTIF-CC-001. The window score is the fraction of heptad a/d core positions
/// occupied by a hydrophobic-core residue (I, L or V), maximised over the seven heptad registers.
/// Expected values are derived by closed-form occupancy counting, not from the implementation output.
/// </summary>
[TestFixture]
[Category("PROTMOTIF-CC-001")]
public class ProteinMotifFinder_PredictCoiledCoils_Tests
{
    // "LAALAAA": a=L (index 0), d=L (index 3) every heptad -> all a/d are hydrophobic core in register 0.
    private static string Repeat(string unit, int times) => string.Concat(Enumerable.Repeat(unit, times));

    #region PredictCoiledCoils

    // M1 — Perfect heptad: L at every a and d. W=28, thr=0.5. All 8 windows score 1.0 (register 0);
    // run [0..7] -> residues [0, 7+27] = [0,34], length 35 >= 21 -> exactly one region (0,34,1.0).
    [Test]
    public void PredictCoiledCoils_PerfectHeptadRepeat_ReturnsSingleFullRegionScoreOne()
    {
        string sequence = Repeat("LAALAAA", 5); // 35 aa

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(1),
            "A 5-heptad repeat with L at every a and d must yield exactly one coiled-coil region.");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0),
                "Region starts at residue 0 (first above-threshold window).");
            Assert.That(regions[0].End, Is.EqualTo(34),
                "Region end = last window index (7) + window (28) - 1 = 34 (inclusive).");
            Assert.That(regions[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "Every a/d position is L in {I,L,V}, so occupancy = 1.0 (Mason & Arndt; Wikipedia).");
        });
    }

    // M2 — No core residues: glycine has no I/L/V, so occupancy is 0 in every register/window.
    [Test]
    public void PredictCoiledCoils_NoCoreResidues_ReturnsNoRegions()
    {
        string sequence = new string('G', 40);

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence).ToList();

        Assert.That(regions, Is.Empty,
            "With no I/L/V at any a/d position the score is 0 < 0.5 threshold (core set {I,L,V}).");
    }

    // M3 — Below window length: 21 aa < window 28, no full window exists -> empty (Lupas window rule, INV-04).
    [Test]
    public void PredictCoiledCoils_SequenceShorterThanWindow_ReturnsEmpty()
    {
        string sequence = Repeat("LAALAAA", 3); // 21 aa < 28

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence).ToList();

        Assert.That(regions, Is.Empty,
            "A sequence shorter than the 28-residue window cannot be scored (INV-04).");
    }

    // M4 — Off-frame coiled coil: the perfect core is shifted by 2 (register 2). The 7-register max must
    // still find it. n=37, W=28 -> windows 0..9 all score 1.0; run [0..9] -> residues [0, 9+27] = [0,36].
    [Test]
    public void PredictCoiledCoils_OffFrameRepeat_FoundViaSevenRegisterMaximum()
    {
        string sequence = "AA" + Repeat("LAALAAA", 5); // 37 aa, core in register 2

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(1),
            "The off-frame coiled coil must be detected by maximising over the 7 heptad registers (INV-05).");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0),
                "All windows score 1.0 in register 2, so the region starts at residue 0.");
            Assert.That(regions[0].End, Is.EqualTo(36),
                "Region end = last window index (9) + window (28) - 1 = 36.");
            Assert.That(regions[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "Register 2 places L at every a/d position -> occupancy 1.0 (Lupas seven frames).");
        });
    }

    // M5 — Region exactly at the window length (28 >= MinRegion 21) is kept. n=28, W=28 -> one window,
    // score 1.0; run [0..0] -> residues [0,27], length 28 >= 21 -> region (0,27,1.0) (INV-02,03).
    [Test]
    public void PredictCoiledCoils_RegionEqualToWindow_IsKeptAtMinimumLength()
    {
        string sequence = Repeat("LAALAAA", 4); // 28 aa

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(1),
            "A single full window of perfect heptads (length 28 >= 21) yields exactly one region.");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "Region starts at residue 0.");
            Assert.That(regions[0].End, Is.EqualTo(27),
                "Region end = window index 0 + window 28 - 1 = 27, length 28 >= MinRegion 21.");
            Assert.That(regions[0].Score, Is.EqualTo(1.0).Within(1e-10), "All a/d are L -> occupancy 1.0.");
        });
    }

    // S1 — Half occupancy boundary: L only at a (index 0), A at d. Best register has 4/8 = 0.5 occupancy.
    // 0.5 >= 0.5 threshold -> region kept with score exactly 0.5.
    [Test]
    public void PredictCoiledCoils_HalfOccupancyAtThreshold_ReturnsRegionScoreOneHalf()
    {
        string sequence = Repeat("LAAAAAA", 5); // 35 aa, L only at position a

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(1),
            "Half-occupancy (0.5) meets the >= 0.5 threshold, so one region is reported.");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "Region starts at residue 0.");
            Assert.That(regions[0].End, Is.EqualTo(34), "Region end = window 7 + 28 - 1 = 34.");
            Assert.That(regions[0].Score, Is.EqualTo(0.5).Within(1e-10),
                "L occupies only a (4 of 8 a/d positions per window) -> occupancy exactly 0.5.");
        });
    }

    // S2 — Just below threshold: same 0.5 occupancy but threshold 0.5001 -> 0.5 < 0.5001 -> no region.
    [Test]
    public void PredictCoiledCoils_HalfOccupancyBelowThreshold_ReturnsEmpty()
    {
        string sequence = Repeat("LAAAAAA", 5);

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence, threshold: 0.5001).ToList();

        Assert.That(regions, Is.Empty,
            "Occupancy 0.5 is strictly below threshold 0.5001, so no region is reported.");
    }

    // S5 — Coiled core followed by a hydrophilic tail: the region must END inside the sequence (the
    // window score drops below threshold once it slides into the all-G tail). 5 perfect heptads (35 aa)
    // + 35 G. Closed-form occupancy scan -> single region (0, 48, 1.0); region ends at 48, not 69.
    [Test]
    public void PredictCoiledCoils_CoreFollowedByHydrophilicTail_RegionEndsInsideSequence()
    {
        string sequence = Repeat("LAALAAA", 5) + new string('G', 35); // 70 aa

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(1),
            "Only the leading heptad core scores above threshold; the G-tail does not.");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "Region starts at residue 0.");
            Assert.That(regions[0].End, Is.EqualTo(48),
                "Region ends at residue 48 (last window whose occupancy >= 0.5) + 28 - 1, before the tail end (69).");
            Assert.That(regions[0].Score, Is.EqualTo(1.0).Within(1e-10),
                "Peak occupancy in the run is 1.0 (perfect heptads).");
        });
    }

    // S3 — Null input -> empty (validation).
    [Test]
    public void PredictCoiledCoils_NullInput_ReturnsEmpty()
    {
        var regions = ProteinMotifFinder.PredictCoiledCoils(null!).ToList();

        Assert.That(regions, Is.Empty, "Null input must yield an empty result without throwing.");
    }

    // S4 — Empty input -> empty (validation).
    [Test]
    public void PredictCoiledCoils_EmptyInput_ReturnsEmpty()
    {
        var regions = ProteinMotifFinder.PredictCoiledCoils(string.Empty).ToList();

        Assert.That(regions, Is.Empty, "Empty input must yield an empty result.");
    }

    // S6 — Two separated cores -> two non-overlapping regions in increasing Start order (INV-03), and
    // exercises the mid-sequence "drop and emit, then continue" branch. core1 (LAALAAA x5) + G x40 +
    // core2 (LAALAAA x5) = 110 aa. The G gap drops the window score to 0 between the cores. Closed-form
    // occupancy scan (independently reproduced) -> region1 (0,48,1.0) [same drop as S5] and region2 starts
    // once windows re-enter core2 at window index 58 -> (58,109,1.0).
    [Test]
    public void PredictCoiledCoils_TwoSeparatedCores_ReturnsTwoNonOverlappingRegions()
    {
        string sequence = Repeat("LAALAAA", 5) + new string('G', 40) + Repeat("LAALAAA", 5); // 110 aa

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(2),
            "Two cores separated by a long hydrophilic gap must yield two distinct regions.");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "First region starts at residue 0.");
            Assert.That(regions[0].End, Is.EqualTo(48), "First region ends at 48 (window 21 + 28 - 1).");
            Assert.That(regions[0].Score, Is.EqualTo(1.0).Within(1e-10), "First core peak occupancy 1.0.");
            Assert.That(regions[1].Start, Is.EqualTo(58), "Second region starts at window index 58.");
            Assert.That(regions[1].End, Is.EqualTo(109), "Second region ends at the last residue (109).");
            Assert.That(regions[1].Score, Is.EqualTo(1.0).Within(1e-10), "Second core peak occupancy 1.0.");
            Assert.That(regions[0].End, Is.LessThan(regions[1].Start),
                "Regions are non-overlapping and in increasing Start order (INV-03).");
        });
    }

    // S7 — Custom windowSize that produces a scoring run shorter than MinRegion (21) -> BuildRegion rejects
    // it -> empty. With windowSize=7 the perfect single heptad scores 1.0 over windows [4..7]; that run maps
    // to residues [4, 7+7-1] = [4,13], length 10 < 21 -> rejected. Exercises the min-region filter branch and
    // the non-default windowSize parameter (independently reproduced by closed-form scan).
    [Test]
    public void PredictCoiledCoils_CustomWindowRunBelowMinRegion_ReturnsEmpty()
    {
        string sequence = new string('G', 7) + "LAALAAA" + new string('G', 7); // 21 aa

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence, windowSize: 7, threshold: 0.99).ToList();

        Assert.That(regions, Is.Empty,
            "A scoring run spanning fewer than 21 residues is rejected by the 3-heptad minimum (INV-02).");
    }

    // S8 — Custom windowSize (14 = 2 heptads) still detects a perfect repeat. LAALAAA x4 (28 aa): every
    // 14-residue window scores 1.0 over all 15 windows; run [0..14] -> residues [0, 14+14-1] = [0,27],
    // length 28 >= 21 -> one region (0,27,1.0). Confirms the windowSize parameter is honoured.
    [Test]
    public void PredictCoiledCoils_CustomWindowSize_DetectsPerfectRepeat()
    {
        string sequence = Repeat("LAALAAA", 4); // 28 aa

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence, windowSize: 14).ToList();

        Assert.That(regions, Has.Count.EqualTo(1), "A 2-heptad window must still detect the perfect repeat.");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "Region starts at residue 0.");
            Assert.That(regions[0].End, Is.EqualTo(27), "Region end = last window 14 + 14 - 1 = 27.");
            Assert.That(regions[0].Score, Is.EqualTo(1.0).Within(1e-10), "All a/d are L -> occupancy 1.0.");
        });
    }

    // C1 — Case-insensitive: lowercase perfect heptad must behave identically to M1.
    [Test]
    public void PredictCoiledCoils_LowercaseInput_TreatedSameAsUppercase()
    {
        string sequence = Repeat("laalaaa", 5);

        var regions = ProteinMotifFinder.PredictCoiledCoils(sequence).ToList();

        Assert.That(regions, Has.Count.EqualTo(1), "Lowercase residues must be recognised (input uppercased).");
        Assert.Multiple(() =>
        {
            Assert.That(regions[0].Start, Is.EqualTo(0), "Same as M1: region starts at 0.");
            Assert.That(regions[0].End, Is.EqualTo(34), "Same as M1: region ends at 34.");
            Assert.That(regions[0].Score, Is.EqualTo(1.0).Within(1e-10), "Same as M1: occupancy 1.0.");
        });
    }

    #endregion
}
