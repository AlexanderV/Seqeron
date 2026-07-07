// DISORDER-MORF-001 — MoRF (Molecular Recognition Feature) Prediction
// Evidence: docs/Evidence/DISORDER-MORF-001-Evidence.md
// TestSpec: tests/TestSpecs/DISORDER-MORF-001.md
// Source: Mohan A et al. (2006) J Mol Biol 362(5):1043-1059 (PMID 16935303);
//         Cheng/Oldfield et al. PMC2570644 ("dips", threshold 0.5);
//         Campen A et al. (2008) Protein Pept Lett 15(9):956-963 (per-residue TOP-IDP scores).
//
// Definition under test: a MoRF is a short region of relative ORDER (a "dip",
// per-residue disorder score < 0.5) embedded WITHIN a longer disordered region
// (flanked by disorder on both sides), of length 10-70 residues.
//
// Expected coordinates and scores are derived independently of the implementation
// from the normalized TOP-IDP scale (Campen 2008 Table 2) and the disorder window
// (size 21) used by PredictDisorder. P (normalized 1.000) is disordered; L (0.298)
// and I (0.213) are ordered. Window smoothing makes the sub-0.5 "dip" narrower than
// the raw ordered run; the dip coordinates ARE the MoRF, per Cheng/Oldfield.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class DisorderPredictor_MoRF_Tests
{
    // 25 disordered P + 30 ordered L + 25 disordered P (length 80).
    // Smoothed disorder profile dips below 0.5 over residues 29-50 (length 22).
    private static readonly string DipInDisorderL =
        new string('P', 25) + new string('L', 30) + new string('P', 25);

    // Same flanks with a deeper-ordering I core (TOP-IDP -0.486 < L's -0.326).
    private static readonly string DipInDisorderI =
        new string('P', 25) + new string('I', 30) + new string('P', 25);

    #region PredictMoRFs — MUST

    [Test]
    public void M1_OrderedDipWithinDisorder_DetectedAtExactCoordinates()
    {
        // 25P + 30L + 25P: a single ordered dip embedded in disorder.
        // Derived: smoothed disorder < 0.5 over residues [29,50] (length 22, in 10-70);
        // mean disorder over the dip = 0.362033 → score = (0.5-0.362033)/0.5 = 0.275934.
        var morfs = DisorderPredictor.PredictMoRFs(DipInDisorderL).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(morfs, Has.Count.EqualTo(1),
                "exactly one dip-in-disorder MoRF expected (Cheng/Oldfield dip definition)");
            Assert.That(morfs[0].Start, Is.EqualTo(29),
                "dip start where smoothed disorder first drops below 0.5");
            Assert.That(morfs[0].End, Is.EqualTo(50),
                "dip end where smoothed disorder rises to >= 0.5");
            Assert.That(morfs[0].Score, Is.EqualTo(0.275934).Within(1e-6),
                "score = (0.5 - mean disorder 0.362033) / 0.5");
        });
    }

    [Test]
    public void M2_FullyOrderedSequence_NoMoRFs()
    {
        // All-L: ordered throughout, no surrounding disorder to embed a dip.
        var morfs = DisorderPredictor.PredictMoRFs(new string('L', 40)).ToList();

        Assert.That(morfs, Is.Empty,
            "no MoRF without disorder flanking the ordered region (Cheng/Oldfield)");
    }

    [Test]
    public void M3_FullyDisorderedSequence_NoMoRFs()
    {
        // All-P: disordered throughout, no ordered dip exists.
        var morfs = DisorderPredictor.PredictMoRFs(new string('P', 40)).ToList();

        Assert.That(morfs, Is.Empty,
            "no MoRF without an ordered dip (Cheng/Oldfield threshold 0.5)");
    }

    [Test]
    public void M4_DipShorterThanTenResidues_NotAMoRF()
    {
        // 25P + 16L + 25P: smoothed dip is only 8 residues (< 10) → not a MoRF.
        string seq = new string('P', 25) + new string('L', 16) + new string('P', 25);

        var morfs = DisorderPredictor.PredictMoRFs(seq).ToList();

        Assert.That(morfs, Is.Empty,
            "dip below the 10-residue minimum is not a MoRF (Mohan 2006 10-70 band)");
    }

    [Test]
    public void M5_DipLongerThanSeventyResidues_NotAMoRF()
    {
        // 25P + 95L + 25P: smoothed dip is 87 residues (> 70) → not a MoRF.
        string seq = new string('P', 25) + new string('L', 95) + new string('P', 25);

        var morfs = DisorderPredictor.PredictMoRFs(seq).ToList();

        Assert.That(morfs, Is.Empty,
            "dip above the 70-residue maximum is not a MoRF (Mohan 2006 10-70 band)");
    }

    [Test]
    public void M6_DipAtTerminus_NotAMoRF()
    {
        // 15L + 30P: the ordered dip touches the N-terminus, so it is not flanked
        // by disorder on both sides → not "within a longer region of disorder".
        string seq = new string('L', 15) + new string('P', 30);

        var morfs = DisorderPredictor.PredictMoRFs(seq).ToList();

        Assert.That(morfs, Is.Empty,
            "terminal dip is not embedded in disorder (Oldfield 2005 / Mohan 2006)");
    }

    [Test]
    public void M7_DeeperDipScoresHigher_BothBounded()
    {
        // I core (TOP-IDP -0.486) gives a deeper dip than L core (-0.326).
        // Derived: I-dip score = 0.399608 > L-dip score = 0.275934; both in [0,1].
        var morfL = DisorderPredictor.PredictMoRFs(DipInDisorderL).Single();
        var morfI = DisorderPredictor.PredictMoRFs(DipInDisorderI).Single();

        Assert.Multiple(() =>
        {
            Assert.That(morfI.Score, Is.EqualTo(0.399608).Within(1e-6),
                "deeper I-dip score = (0.5 - mean disorder 0.300196) / 0.5");
            Assert.That(morfL.Score, Is.EqualTo(0.275934).Within(1e-6),
                "shallower L-dip score = (0.5 - mean disorder 0.362033) / 0.5");
            Assert.That(morfI.Score, Is.GreaterThan(morfL.Score),
                "score is monotone in dip depth (deeper order = higher score)");
            Assert.That(morfI.Score, Is.InRange(0.0, 1.0), "score bounded in [0,1]");
            Assert.That(morfL.Score, Is.InRange(0.0, 1.0), "score bounded in [0,1]");
        });
    }

    #endregion

    #region PredictMoRFs — SHOULD

    [Test]
    public void S1_TwoSeparateDips_TwoNonOverlappingMoRFs()
    {
        // 25P + 30L + 30P + 30L + 25P: two ordered dips separated by disorder.
        // Derived: dips at [29,50] and [89,110].
        string seq = new string('P', 25) + new string('L', 30) + new string('P', 30)
                   + new string('L', 30) + new string('P', 25);

        var morfs = DisorderPredictor.PredictMoRFs(seq).OrderBy(m => m.Start).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(morfs, Has.Count.EqualTo(2), "two independent dips → two MoRFs");
            Assert.That(morfs[0].Start, Is.EqualTo(29), "first dip start");
            Assert.That(morfs[0].End, Is.EqualTo(50), "first dip end");
            Assert.That(morfs[1].Start, Is.EqualTo(89), "second dip start");
            Assert.That(morfs[1].End, Is.EqualTo(110), "second dip end");
            Assert.That(morfs[1].Start, Is.GreaterThan(morfs[0].End), "non-overlapping");
        });
    }

    [Test]
    public void S2_CaseInsensitive_SameResultAsUpperCase()
    {
        var upper = DisorderPredictor.PredictMoRFs(DipInDisorderL).ToList();
        var lower = DisorderPredictor.PredictMoRFs(DipInDisorderL.ToLowerInvariant()).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(lower, Has.Count.EqualTo(upper.Count), "same MoRF count");
            Assert.That(lower[0].Start, Is.EqualTo(upper[0].Start), "same start");
            Assert.That(lower[0].End, Is.EqualTo(upper[0].End), "same end");
            Assert.That(lower[0].Score, Is.EqualTo(upper[0].Score).Within(1e-12), "same score");
        });
    }

    [Test]
    public void S3_CustomMaxLength_ExcludesDipBeyondBound()
    {
        // The L-dip has length 22. With maxLength 20 it falls outside the band → excluded.
        var morfs = DisorderPredictor.PredictMoRFs(DipInDisorderL, minLength: 10, maxLength: 20).ToList();

        Assert.That(morfs, Is.Empty,
            "dip length 22 exceeds custom maxLength 20 → not reported (INV-2)");
    }

    [Test]
    public void S3_CustomMinLength_ExcludesDipBelowBound()
    {
        // The L-dip has length 22. With minLength 23 it falls below the band → excluded.
        var morfs = DisorderPredictor.PredictMoRFs(DipInDisorderL, minLength: 23, maxLength: 70).ToList();

        Assert.That(morfs, Is.Empty,
            "dip length 22 below custom minLength 23 → not reported (INV-2)");
    }

    #endregion

    #region PredictMoRFs — COULD (guards)

    [Test]
    public void C1_NullInput_ReturnsEmpty()
    {
        var morfs = DisorderPredictor.PredictMoRFs(null!).ToList();
        Assert.That(morfs, Is.Empty, "null input → empty");
    }

    [Test]
    public void C1_EmptyInput_ReturnsEmpty()
    {
        var morfs = DisorderPredictor.PredictMoRFs("").ToList();
        Assert.That(morfs, Is.Empty, "empty input → empty");
    }

    [Test]
    public void C2_SequenceShorterThanMinLength_ReturnsEmpty()
    {
        var morfs = DisorderPredictor.PredictMoRFs("EEKPP").ToList();
        Assert.That(morfs, Is.Empty, "5-residue sequence cannot contain a 10-residue MoRF");
    }

    #endregion

    #region Invariants

    [Test]
    public void INV_BoundariesScoreAndOrdering_HoldOnMultiDipSequence()
    {
        string seq = new string('P', 25) + new string('L', 30) + new string('P', 30)
                   + new string('L', 30) + new string('P', 25);

        var prediction = DisorderPredictor.PredictDisorder(seq);
        var morfs = DisorderPredictor.PredictMoRFs(seq).OrderBy(m => m.Start).ToList();

        Assert.Multiple(() =>
        {
            for (int i = 0; i < morfs.Count; i++)
            {
                var m = morfs[i];
                // INV-1: valid coordinates.
                Assert.That(m.Start, Is.GreaterThanOrEqualTo(0), "INV-1 start >= 0");
                Assert.That(m.End, Is.LessThan(seq.Length), "INV-1 end < length");
                Assert.That(m.Start, Is.LessThanOrEqualTo(m.End), "INV-1 start <= end");

                // INV-2: length within default band 10-70.
                int len = m.End - m.Start + 1;
                Assert.That(len, Is.InRange(10, 70), "INV-2 length in 10-70");

                // INV-5: score bounded.
                Assert.That(m.Score, Is.InRange(0.0, 1.0), "INV-5 score in [0,1]");

                // INV-3: flanked by disorder (score >= 0.5) on both sides.
                Assert.That(prediction.ResiduePredictions[m.Start - 1].DisorderScore,
                    Is.GreaterThanOrEqualTo(0.5), "INV-3 left flank disordered");
                Assert.That(prediction.ResiduePredictions[m.End + 1].DisorderScore,
                    Is.GreaterThanOrEqualTo(0.5), "INV-3 right flank disordered");

                // INV-4: non-overlapping, ordered.
                if (i > 0)
                    Assert.That(m.Start, Is.GreaterThan(morfs[i - 1].End),
                        "INV-4 non-overlapping and ordered by start");
            }
        });
    }

    #endregion
}
