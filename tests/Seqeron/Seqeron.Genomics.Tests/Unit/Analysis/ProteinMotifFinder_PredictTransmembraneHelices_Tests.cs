// PROTMOTIF-TM-001 — Transmembrane Helix Prediction (Kyte-Doolittle hydropathy sliding window)
// Evidence: docs/Evidence/PROTMOTIF-TM-001-Evidence.md
// TestSpec: tests/TestSpecs/PROTMOTIF-TM-001.md
// Source:   Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132; Davidson KD background page;
//           QIAGEN CLC Hydrophobicity scales; Biopython Bio.SeqUtils.ProtParam.protein_scale.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// Canonical tests for PROTMOTIF-TM-001: Kyte-Doolittle (1982) transmembrane-helix prediction.
/// Expected values are derived from the published window/threshold rule (window 19, threshold 1.6)
/// and the Kyte-Doolittle hydropathy scale; the window profile is the arithmetic mean of the
/// window's per-residue values (Biopython protein_scale, edge weight 1.0).
/// </summary>
[TestFixture]
[Category("PROTMOTIF-TM-001")]
public class ProteinMotifFinder_PredictTransmembraneHelices_Tests
{
    // Kyte-Doolittle scale values used to derive expected scores (QIAGEN / Davidson tables).
    private const double LeuKd = 3.8;  // L
    private const double IleKd = 4.5;  // I
    private const double ValKd = 4.2;  // V

    #region PredictTransmembraneHelices

    // M1 — Single hydrophobic stretch (D×10 + L×20 + D×10), window 19, threshold 1.6.
    // Per the KD rule the run of above-threshold windows is profile indices 5..16 (hand-computed,
    // arithmetic mean over each 19-residue window). The first above-threshold window starts at
    // residue 5; the last starts at residue 16 and covers residues 16..34, so the segment spans
    // residues 5..34 (every residue lying within at least one above-threshold window). Peak = mean
    // of an all-Leu window = 3.8.
    [Test]
    public void PredictTransmembraneHelices_SingleHydrophobicStretch_ReturnsOneSegmentWithExactBounds()
    {
        string sequence = new string('D', 10) + new string('L', 20) + new string('D', 10);

        var segments = ProteinMotifFinder.PredictTransmembraneHelices(sequence).ToList();

        Assert.That(segments, Has.Count.EqualTo(1),
            "A single internal poly-Leu stretch must yield exactly one transmembrane segment.");
        Assert.Multiple(() =>
        {
            Assert.That(segments[0].Start, Is.EqualTo(5),
                "Segment start = first profile index whose window mean ≥ 1.6 (KD window 19 rule).");
            Assert.That(segments[0].End, Is.EqualTo(34),
                "Segment end = last covered residue = last above-threshold window start (16) + window (19) - 1 = 34.");
            Assert.That(segments[0].Score, Is.EqualTo(LeuKd).Within(1e-10),
                "Peak score = mean of an all-Leu window = 3.8 (KD value for L).");
        });
    }

    // M2 — All-hydrophilic sequence: D = -3.5 < 1.6, so no window mean crosses the threshold.
    [Test]
    public void PredictTransmembraneHelices_AllHydrophilic_ReturnsNoSegments()
    {
        string sequence = new string('D', 40);

        var segments = ProteinMotifFinder.PredictTransmembraneHelices(sequence).ToList();

        Assert.That(segments, Is.Empty,
            "Every window mean equals the Asp KD value (-3.5), which is below the 1.6 threshold.");
    }

    // M3 — Exactly one window of poly-Leu (length == window): single profile point of 3.8 ≥ 1.6.
    [Test]
    public void PredictTransmembraneHelices_ExactlyOneWindowOfLeucine_ReturnsSingleSpanningSegment()
    {
        string sequence = new string('L', 19);

        var segments = ProteinMotifFinder.PredictTransmembraneHelices(sequence).ToList();

        Assert.That(segments, Has.Count.EqualTo(1),
            "A 19-residue all-Leu sequence forms exactly one full-length transmembrane segment.");
        Assert.Multiple(() =>
        {
            Assert.That(segments[0].Start, Is.EqualTo(0), "Single window starts at residue 0.");
            Assert.That(segments[0].End, Is.EqualTo(18),
                "End clamped to the last residue index (length-1 = 18).");
            Assert.That(segments[0].Score, Is.EqualTo(LeuKd).Within(1e-10),
                "Peak score = KD value for Leu (3.8).");
        });
    }

    // M4 — Scale-value reproduction: a uniform 19-residue window reproduces the residue's KD value
    // as the peak score; Arg (-4.5 < 1.6) yields no segment.
    [Test]
    public void PredictTransmembraneHelices_UniformWindow_ReproducesKyteDoolittleValue()
    {
        var ileSegments = ProteinMotifFinder.PredictTransmembraneHelices(new string('I', 19)).ToList();
        var valSegments = ProteinMotifFinder.PredictTransmembraneHelices(new string('V', 19)).ToList();
        var argSegments = ProteinMotifFinder.PredictTransmembraneHelices(new string('R', 19)).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(ileSegments, Has.Count.EqualTo(1), "Poly-Ile (4.5 ≥ 1.6) is a TM segment.");
            Assert.That(ileSegments[0].Score, Is.EqualTo(IleKd).Within(1e-10),
                "Peak = KD value for Ile (4.5).");
            Assert.That(valSegments, Has.Count.EqualTo(1), "Poly-Val (4.2 ≥ 1.6) is a TM segment.");
            Assert.That(valSegments[0].Score, Is.EqualTo(ValKd).Within(1e-10),
                "Peak = KD value for Val (4.2).");
            Assert.That(argSegments, Is.Empty,
                "Poly-Arg window mean = -4.5 < 1.6, so no TM segment is reported.");
        });
    }

    // M5 — Null input returns an empty result (no window possible).
    [Test]
    public void PredictTransmembraneHelices_NullSequence_ReturnsEmpty()
    {
        var segments = ProteinMotifFinder.PredictTransmembraneHelices(null!).ToList();

        Assert.That(segments, Is.Empty, "Null input cannot form a window; result is empty.");
    }

    // M6 — Empty input returns an empty result.
    [Test]
    public void PredictTransmembraneHelices_EmptySequence_ReturnsEmpty()
    {
        var segments = ProteinMotifFinder.PredictTransmembraneHelices(string.Empty).ToList();

        Assert.That(segments, Is.Empty, "Empty input cannot form a window; result is empty.");
    }

    // M7 — Sequence shorter than the window returns an empty result.
    [Test]
    public void PredictTransmembraneHelices_ShorterThanWindow_ReturnsEmpty()
    {
        string sequence = new string('L', 18); // 18 < default window 19

        var segments = ProteinMotifFinder.PredictTransmembraneHelices(sequence).ToList();

        Assert.That(segments, Is.Empty,
            "A sequence shorter than the 19-residue window yields no profile and no segments.");
    }

    // S1 — A non-standard residue (X) carries no KD value and is excluded from the window mean;
    // the window of 18 Leu + 1 X averages to 3.8 (mean over the 18 scored residues).
    [Test]
    public void PredictTransmembraneHelices_NonStandardResidueInWindow_ExcludedFromMean()
    {
        string sequence = new string('L', 9) + "X" + new string('L', 9); // length 19

        var segments = ProteinMotifFinder.PredictTransmembraneHelices(sequence).ToList();

        Assert.That(segments, Has.Count.EqualTo(1),
            "One window with 18 Leu and one X still exceeds the threshold.");
        Assert.That(segments[0].Score, Is.EqualTo(LeuKd).Within(1e-10),
            "X has no scale value; mean is taken over the 18 Leu residues = 3.8.");
    }

    // S2 — Custom higher threshold suppresses a weak (Ala, 1.8) segment that passes at 1.6.
    [Test]
    public void PredictTransmembraneHelices_HigherThreshold_SuppressesWeakSegment()
    {
        string sequence = new string('D', 10) + new string('A', 20) + new string('D', 10);

        var atDefault = ProteinMotifFinder.PredictTransmembraneHelices(sequence, threshold: 1.6).ToList();
        var atRaised = ProteinMotifFinder.PredictTransmembraneHelices(sequence, threshold: 2.0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(atDefault, Has.Count.EqualTo(1),
                "Poly-Ala peak window mean (1.8) ≥ 1.6, so one segment is reported at the default threshold.");
            Assert.That(atRaised, Is.Empty,
                "Raising the threshold to 2.0 rejects the Ala stretch (peak 1.8 < 2.0).");
        });
    }

    // S3 — Lowercase input is accepted (case-insensitive) and gives the same segment as M1.
    [Test]
    public void PredictTransmembraneHelices_LowercaseInput_MatchesUppercase()
    {
        string sequence = (new string('d', 10) + new string('l', 20) + new string('d', 10));

        var segments = ProteinMotifFinder.PredictTransmembraneHelices(sequence).ToList();

        Assert.That(segments, Has.Count.EqualTo(1), "Lowercase sequence must be handled identically.");
        Assert.Multiple(() =>
        {
            Assert.That(segments[0].Start, Is.EqualTo(5), "Same start as the uppercase case (M1).");
            Assert.That(segments[0].End, Is.EqualTo(34), "Same end as the uppercase case (M1).");
            Assert.That(segments[0].Score, Is.EqualTo(LeuKd).Within(1e-10), "Same peak score as M1.");
        });
    }

    // S4 — A non-positive window size is guarded and returns an empty result.
    [Test]
    public void PredictTransmembraneHelices_NonPositiveWindow_ReturnsEmpty()
    {
        var segments = ProteinMotifFinder.PredictTransmembraneHelices(new string('L', 30), windowSize: 0).ToList();

        Assert.That(segments, Is.Empty, "A window size of 0 is invalid and yields no segments.");
    }

    // C1 — Property: every reported segment satisfies INV-01 (peak ≥ threshold) and INV-02 (Start ≤ End,
    // within bounds) on a real multi-segment sequence.
    [Test]
    public void PredictTransmembraneHelices_ReportedSegments_SatisfyDetectionInvariants()
    {
        const double threshold = 1.6;
        string sequence = new string('D', 10) + new string('L', 20) + new string('D', 10);

        var segments = ProteinMotifFinder.PredictTransmembraneHelices(sequence, threshold: threshold).ToList();

        Assert.That(segments, Is.Not.Empty, "Test input must produce at least one segment.");
        Assert.Multiple(() =>
        {
            foreach (var seg in segments)
            {
                Assert.That(seg.Score, Is.GreaterThanOrEqualTo(threshold),
                    "INV-01: a segment's peak window mean must be ≥ the threshold.");
                Assert.That(seg.Start, Is.LessThanOrEqualTo(seg.End),
                    "INV-02: Start ≤ End.");
                Assert.That(seg.Start, Is.GreaterThanOrEqualTo(0), "INV-02: Start ≥ 0.");
                Assert.That(seg.End, Is.LessThan(sequence.Length), "INV-02: End within bounds.");
            }
        });
    }

    #endregion
}
