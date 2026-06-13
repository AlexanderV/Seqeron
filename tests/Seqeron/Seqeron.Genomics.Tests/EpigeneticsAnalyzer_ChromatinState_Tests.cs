// EPIGEN-CHROM-001 — Chromatin State Prediction
// Evidence: docs/Evidence/EPIGEN-CHROM-001-Evidence.md
// TestSpec: tests/TestSpecs/EPIGEN-CHROM-001.md
// Source: Ernst & Kellis (2012). Nat Methods 9(3):215–216 (ChromHMM: marks modeled present/absent);
//         Roadmap Epigenomics chromatin state learning (15/18-state mark→state signatures);
//         Liang (2004) PNAS 101:7357 (H3K4me3→active promoter);
//         Rada-Iglesias (2018) Nat Genet 50:4 (H3K4me1→enhancer);
//         Creyghton (2010) PNAS 107:21931 (H3K27ac→active enhancer);
//         Ferrari (2014) Mol Cell 53:49 (H3K27me3→Polycomb); Nicetto (2019) Science 363:294 (H3K9me3→heterochromatin).

namespace Seqeron.Genomics.Tests;

using ChromatinState = EpigeneticsAnalyzer.ChromatinState;

[TestFixture]
public class EpigeneticsAnalyzer_ChromatinState_Tests
{
    // Helpers keep each test focused on the one mark combination under test.
    // Signals are deliberately well above (0.9) / below (0.0) the 0.5 presence call
    // so results do not depend on the exact default threshold (Evidence Assumption 1).
    private const double Present = 0.9;
    private const double Absent = 0.0;

    #region PredictChromatinState

    // M1 — H3K4me3 (+H3K27ac) present → active promoter (Roadmap TssA; Liang 2004)
    [Test]
    public void PredictChromatinState_K4me3AndK27ac_ActivePromoter()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Present, h3k4me1: Absent, h3k27ac: Present,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.ActivePromoter),
            "H3K4me3 is the canonical active-promoter (TssA) mark (Roadmap; Liang 2004)");
    }

    // M2 — H3K4me3 alone → active promoter (TssA only needs H3K4me3)
    [Test]
    public void PredictChromatinState_K4me3Only_ActivePromoter()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Present, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.ActivePromoter),
            "TssA signature is H3K4me3; H3K27ac is not required (Roadmap)");
    }

    // M3 — H3K4me1 + H3K27ac → active enhancer (active Enh; Creyghton 2010)
    [Test]
    public void PredictChromatinState_K4me1AndK27ac_ActiveEnhancer()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Absent, h3k4me1: Present, h3k27ac: Present,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.ActiveEnhancer),
            "H3K27ac on an H3K4me1 enhancer marks it active (Creyghton 2010)");
    }

    // M4 — H3K4me1 without H3K27ac → weak/poised enhancer (Roadmap Enh; Rada-Iglesias 2018)
    [Test]
    public void PredictChromatinState_K4me1NoK27ac_WeakEnhancer()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Absent, h3k4me1: Present, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.WeakEnhancer),
            "H3K4me1 without H3K27ac is a poised/weak enhancer (Roadmap Enh)");
    }

    // M5 — H3K36me3 → transcribed gene body (Roadmap Tx)
    [Test]
    public void PredictChromatinState_K36me3_Transcribed()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Absent, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Present, h3k27me3: Absent, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.Transcribed),
            "H3K36me3 marks transcribed gene bodies (Roadmap Tx)");
    }

    // M6 — H3K27me3 alone → Polycomb-repressed (Roadmap ReprPC; Ferrari 2014)
    [Test]
    public void PredictChromatinState_K27me3Only_Repressed()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Absent, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Present, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.Repressed),
            "H3K27me3 alone is the Polycomb-repressed state (Roadmap ReprPC)");
    }

    // M7 — H3K9me3 alone → heterochromatin (Roadmap Het; Nicetto 2019)
    [Test]
    public void PredictChromatinState_K9me3Only_Heterochromatin()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Absent, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Present);

        Assert.That(state, Is.EqualTo(ChromatinState.Heterochromatin),
            "H3K9me3 marks constitutive heterochromatin (Roadmap Het)");
    }

    // M8 — H3K4me3 + H3K27me3 → bivalent promoter (Roadmap TssBiv), not active or repressed
    [Test]
    public void PredictChromatinState_K4me3AndK27me3_BivalentPromoter()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Present, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Present, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.BivalentPromoter),
            "Co-occurring H3K4me3 and H3K27me3 is the bivalent/poised TSS (Roadmap TssBiv), not pure active/repressed");
    }

    // M9 — H3K4me1 + H3K27me3 → bivalent enhancer (Roadmap EnhBiv)
    [Test]
    public void PredictChromatinState_K4me1AndK27me3_BivalentEnhancer()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Absent, h3k4me1: Present, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Present, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.BivalentEnhancer),
            "Co-occurring H3K4me1 and H3K27me3 is the bivalent enhancer (Roadmap EnhBiv)");
    }

    // M10 — no mark present → quiescent/low (Roadmap Quies) — INV-02
    [Test]
    public void PredictChromatinState_NoMarks_LowSignal()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Absent, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.LowSignal),
            "No mark above the presence call → quiescent/low (Roadmap Quies)");
    }

    // M11 — binary invariance (INV-01): same present/absent pattern, different magnitudes → same state
    [Test]
    public void PredictChromatinState_BinaryInvariance_SameStateRegardlessOfMagnitude()
    {
        // Both have only H3K4me3 present (one just over, one near max) → same state.
        var low = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.6, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent);
        var high = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.99, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent);

        Assert.Multiple(() =>
        {
            Assert.That(low, Is.EqualTo(ChromatinState.ActivePromoter),
                "0.6 is above the 0.5 call → present");
            Assert.That(high, Is.EqualTo(low),
                "ChromHMM models marks present/absent; magnitude above the call does not change the state (INV-01)");
        });
    }

    // M12 — promoter precedence (INV via Roadmap TSS>Enh): H3K4me3 + H3K4me1 → active promoter
    [Test]
    public void PredictChromatinState_K4me3AndK4me1_PromoterWins()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: Present, h3k4me1: Present, h3k27ac: Present,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent);

        Assert.That(state, Is.EqualTo(ChromatinState.ActivePromoter),
            "When both promoter (H3K4me3) and enhancer (H3K4me1) marks are present, TSS ranks above Enh (Roadmap)");
    }

    // S1 — threshold boundary: a mark exactly at the threshold counts as present (inclusive >=)
    [Test]
    public void PredictChromatinState_MarkExactlyAtThreshold_CountsAsPresent()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.5, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent,
            presenceThreshold: 0.5);

        Assert.That(state, Is.EqualTo(ChromatinState.ActivePromoter),
            "Presence call is inclusive (signal >= threshold), so 0.5 at threshold 0.5 is present");
    }

    // S2 — negative/zero signals are absent → low signal
    [Test]
    public void PredictChromatinState_NegativeSignals_LowSignal()
    {
        var state = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: -1.0, h3k4me1: -0.5, h3k27ac: 0.0,
            h3k36me3: -2.0, h3k27me3: 0.0, h3k9me3: -0.1);

        Assert.That(state, Is.EqualTo(ChromatinState.LowSignal),
            "Signals below the presence call (incl. negatives) are absent → LowSignal");
    }

    // C1 — custom (high) threshold suppresses a marginal mark, changing the call
    [Test]
    public void PredictChromatinState_HighThreshold_SuppressesMarginalMark()
    {
        var withDefault = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.6, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent,
            presenceThreshold: 0.5);
        var withHigh = EpigeneticsAnalyzer.PredictChromatinState(
            h3k4me3: 0.6, h3k4me1: Absent, h3k27ac: Absent,
            h3k36me3: Absent, h3k27me3: Absent, h3k9me3: Absent,
            presenceThreshold: 0.7);

        Assert.Multiple(() =>
        {
            Assert.That(withDefault, Is.EqualTo(ChromatinState.ActivePromoter),
                "0.6 >= 0.5 → present → active promoter");
            Assert.That(withHigh, Is.EqualTo(ChromatinState.LowSignal),
                "0.6 < 0.7 → absent → low signal");
        });
    }

    #endregion

    #region AnnotateHistoneModifications

    // M13 — each region labelled by its single mark's canonical Roadmap state
    [Test]
    public void AnnotateHistoneModifications_PerMark_AssignsCanonicalState()
    {
        var mods = new List<(int, int, string, double)>
        {
            (0, 1000, "H3K4me3", 0.8),    // active promoter (TssA)
            (2000, 3000, "H3K27me3", 0.7), // repressed (ReprPC)
            (4000, 5000, "H3K36me3", 0.9), // transcribed (Tx)
            (6000, 7000, "H3K9me3", 0.9),  // heterochromatin (Het)
        };

        var ann = EpigeneticsAnalyzer.AnnotateHistoneModifications(mods).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(ann, Has.Count.EqualTo(4), "one annotation per input region");
            Assert.That(ann[0].PredictedState, Is.EqualTo(ChromatinState.ActivePromoter), "H3K4me3 → TssA");
            Assert.That(ann[1].PredictedState, Is.EqualTo(ChromatinState.Repressed), "H3K27me3 → ReprPC");
            Assert.That(ann[2].PredictedState, Is.EqualTo(ChromatinState.Transcribed), "H3K36me3 → Tx");
            Assert.That(ann[3].PredictedState, Is.EqualTo(ChromatinState.Heterochromatin), "H3K9me3 → Het");
            Assert.That(ann[0].Mark, Is.EqualTo("H3K4me3"), "mark identity is preserved");
        });
    }

    // S? — a mark below the presence call → LowSignal label
    [Test]
    public void AnnotateHistoneModifications_BelowThreshold_LowSignal()
    {
        var mods = new List<(int, int, string, double)> { (0, 1000, "H3K4me3", 0.1) };

        var ann = EpigeneticsAnalyzer.AnnotateHistoneModifications(mods).ToList();

        Assert.That(ann[0].PredictedState, Is.EqualTo(ChromatinState.LowSignal),
            "Signal 0.1 below the 0.5 presence call → LowSignal");
    }

    // M13b — enhancer/acetylation marks and unknown marks map per the Roadmap single-mark table
    [Test]
    public void AnnotateHistoneModifications_EnhancerAndUnknownMarks_MappedCorrectly()
    {
        var mods = new List<(int, int, string, double)>
        {
            (0, 100, "H3K4me1", 0.8),   // weak/poised enhancer (Enh; Rada-Iglesias 2018)
            (100, 200, "H3K27ac", 0.8), // active enhancer (Creyghton 2010)
            (200, 300, "H3K9ac", 0.8),  // H3K9 acetylation → activation (Wang 2008)
            (300, 400, "FOObar", 0.8),  // unrecognised mark → LowSignal
        };

        var ann = EpigeneticsAnalyzer.AnnotateHistoneModifications(mods).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(ann[0].PredictedState, Is.EqualTo(ChromatinState.WeakEnhancer), "H3K4me1 → poised/weak enhancer");
            Assert.That(ann[1].PredictedState, Is.EqualTo(ChromatinState.ActiveEnhancer), "H3K27ac → active enhancer");
            Assert.That(ann[2].PredictedState, Is.EqualTo(ChromatinState.ActivePromoter), "H3K9ac → active promoter");
            Assert.That(ann[3].PredictedState, Is.EqualTo(ChromatinState.LowSignal), "unrecognised mark → LowSignal");
        });
    }

    // S4 — empty input → empty result
    [Test]
    public void AnnotateHistoneModifications_Empty_ReturnsEmpty()
    {
        var ann = EpigeneticsAnalyzer
            .AnnotateHistoneModifications(new List<(int, int, string, double)>())
            .ToList();

        Assert.That(ann, Is.Empty, "no regions in → no annotations out");
    }

    #endregion

    #region FindAccessibleRegions

    // M14 — contiguous above-threshold positions merge into one region (INV-05)
    [Test]
    public void FindAccessibleRegions_ContiguousSignal_MergesOneRegion()
    {
        // 21 samples at 0.7, 10 bp apart → span 0..200 (>= minWidth 100), no gap > maxGap.
        var signal = Enumerable.Range(0, 21).Select(i => (i * 10, 0.7)).ToList();

        var regions = EpigeneticsAnalyzer.FindAccessibleRegions(
            signal, threshold: 0.5, minWidth: 100, maxGap: 50).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(regions, Has.Count.EqualTo(1), "contiguous above-threshold positions form one region");
            Assert.That(regions[0].Start, Is.EqualTo(0), "region starts at first above-threshold position");
            Assert.That(regions[0].End, Is.EqualTo(200), "region ends at last above-threshold position");
            Assert.That(regions[0].End, Is.GreaterThanOrEqualTo(regions[0].Start), "End >= Start (INV-05)");
            Assert.That(regions[0].AccessibilityScore, Is.EqualTo(0.7).Within(1e-10),
                "score is the max signal in the region, all 0.7");
        });
    }

    // M14b — peak strength label reflects the region's max score (Strong > 0.8)
    [Test]
    public void FindAccessibleRegions_HighScore_StrongPeakLabel()
    {
        var signal = Enumerable.Range(0, 21).Select(i => (i * 10, 0.9)).ToList();

        var regions = EpigeneticsAnalyzer.FindAccessibleRegions(
            signal, threshold: 0.5, minWidth: 100, maxGap: 50).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(regions, Has.Count.EqualTo(1), "one contiguous region");
            Assert.That(regions[0].PeakType, Is.EqualTo("Strong"), "max score 0.9 (> 0.8) → Strong label");
        });
    }

    // S5 — a region narrower than minWidth is excluded
    [Test]
    public void FindAccessibleRegions_NarrowerThanMinWidth_Excluded()
    {
        // Span 0..50 only (< minWidth 100).
        var signal = new List<(int, double)> { (0, 0.9), (10, 0.9), (20, 0.9), (30, 0.9), (40, 0.9), (50, 0.9) };

        var regions = EpigeneticsAnalyzer.FindAccessibleRegions(
            signal, threshold: 0.5, minWidth: 100).ToList();

        Assert.That(regions, Is.Empty, "a region narrower than minWidth (50 < 100) is filtered out");
    }

    // S3 — empty signal → empty result
    [Test]
    public void FindAccessibleRegions_EmptySignal_ReturnsEmpty()
    {
        var regions = EpigeneticsAnalyzer
            .FindAccessibleRegions(new List<(int, double)>())
            .ToList();

        Assert.That(regions, Is.Empty, "no signal in → no regions out");
    }

    // S? — a gap larger than maxGap splits into separate regions
    [Test]
    public void FindAccessibleRegions_GapLargerThanMaxGap_SplitsRegions()
    {
        var signal = new List<(int, double)>();
        for (int i = 0; i <= 150; i += 10) signal.Add((i, 0.8));   // region 1: 0..150
        for (int i = 300; i <= 450; i += 10) signal.Add((i, 0.8)); // region 2 after gap of 150

        var regions = EpigeneticsAnalyzer.FindAccessibleRegions(
            signal, threshold: 0.5, minWidth: 100, maxGap: 50).ToList();

        Assert.That(regions, Has.Count.EqualTo(2),
            "a gap (150) larger than maxGap (50) separates the track into two regions");
    }

    #endregion
}
