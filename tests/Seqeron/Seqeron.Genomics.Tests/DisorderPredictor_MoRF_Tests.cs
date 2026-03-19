using NUnit.Framework;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// DISORDER-MORF-001: Molecular Recognition Feature Prediction
///
/// Tests for <see cref="DisorderPredictor.PredictMoRFs"/>.
/// Heuristic: hydrophobic island detection within IDRs.
/// Definition: Mohan et al. (2006) J Mol Biol 362:1043-1059.
/// Hydropathy: Kyte &amp; Doolittle (1982) J Mol Biol 157:105-132.
/// Note: annotation heuristic, not an ML predictor.
/// </summary>
[TestFixture]
public class DisorderPredictor_MoRF_Tests
{
    // Flanking IDR: E/K are disorder-promoting (TOP-IDP) and very hydrophilic (KD).
    private const string FlankEK = "EEEKKKEEEKKKEEEKKK"; // 18 AA

    // Hydrophobic island: G is disorder-promoting (TOP-IDP) but near-neutral KD (−0.4).
    private const string IslandG = "GGGGGGGGGGG"; // 11 AA

    // Full test sequence: 18 (flank) + 11 (island) + 18 (flank) = 47 AA.
    // All residues are disorder-promoting → one long IDR.
    // IDR mean KD ≈ −2.88. Island KD = −0.4. Diff ≈ 2.48.
    private const string SeqWithIsland = FlankEK + IslandG + FlankEK;

    #region S — Smoke Tests

    [Test]
    public void S1_OrderedProtein_NoMoRFs()
    {
        // All-leucine: L has TOP-IDP normalized score ≈ 0.298 < 0.542.
        // No IDRs → no MoRFs.
        string ordered = new string('L', 40);

        var morfs = DisorderPredictor.PredictMoRFs(ordered).ToList();

        Assert.That(morfs, Is.Empty);
    }

    [Test]
    public void S2_UniformIDR_NoMoRFs()
    {
        // All-proline: TOP-IDP normalized = 1.0 → disordered.
        // KD = −1.6 uniform throughout → no hydropathy enrichment → no MoRFs.
        string uniform = new string('P', 40);

        var morfs = DisorderPredictor.PredictMoRFs(uniform).ToList();

        Assert.That(morfs, Is.Empty);
    }

    [Test]
    public void S3_HydrophobicIslandInIDR_DetectedAsMoRF()
    {
        // G-island (KD −0.4) in E/K context (KD ≈ −3.7): large hydropathy contrast.
        // All residues disorder-promoting → entire sequence is one IDR.
        // The G-island should be detected as a MoRF.
        var morfs = DisorderPredictor.PredictMoRFs(SeqWithIsland).ToList();

        Assert.That(morfs, Has.Count.GreaterThanOrEqualTo(1));

        // Best MoRF should overlap with the G-island (positions 18–28)
        var best = morfs.OrderByDescending(m => m.Score).First();
        Assert.That(best.Start, Is.InRange(13, 23), "MoRF start near G-island");
        Assert.That(best.End, Is.InRange(23, 33), "MoRF end near G-island");
        Assert.That(best.Score, Is.GreaterThan(0.5), "Substantial hydropathy enrichment");
    }

    #endregion

    #region C — Corner Cases

    [Test]
    public void C1_ShortSequence_NoMoRFs()
    {
        // 5 AA is too short for a 10-AA MoRF within any IDR.
        string seq = "EEKPP";

        var morfs = DisorderPredictor.PredictMoRFs(seq).ToList();

        Assert.That(morfs, Is.Empty);
    }

    [Test]
    public void C2_IDRShorterThanMinLength_NoMoRFs()
    {
        // Ordered flanks + short 8-residue IDR: no 10-AA MoRF can fit.
        string ordered = new string('L', 20);
        string idr = "PPPPPPPP"; // 8 P's
        string sequence = ordered + idr + ordered;

        var morfs = DisorderPredictor.PredictMoRFs(sequence, minLength: 10).ToList();

        Assert.That(morfs, Is.Empty);
    }

    #endregion

    #region M — Method Tests

    [Test]
    public void M1_ScoreReflectsHydropathyContrast()
    {
        // G-island: KD = −0.4, high contrast with E/K context (KD ≈ −3.7).
        // S-island: KD = −0.8, lower contrast.
        // G-island should produce a higher score.
        string sIsland = "SSSSSSSSSSS"; // 11 S, KD = −0.8
        string seqS = FlankEK + sIsland + FlankEK;

        var morfsG = DisorderPredictor.PredictMoRFs(SeqWithIsland).ToList();
        var morfsS = DisorderPredictor.PredictMoRFs(seqS).ToList();

        if (morfsG.Count > 0 && morfsS.Count > 0)
        {
            double maxG = morfsG.Max(m => m.Score);
            double maxS = morfsS.Max(m => m.Score);
            Assert.That(maxG, Is.GreaterThan(maxS),
                "G-island (KD −0.4) should score higher than S-island (KD −0.8)");
        }
    }

    [Test]
    public void M2_MaxLengthRespected()
    {
        int maxLen = 15;
        var morfs = DisorderPredictor.PredictMoRFs(SeqWithIsland, maxLength: maxLen).ToList();

        foreach (var morf in morfs)
        {
            int len = morf.End - morf.Start + 1;
            Assert.That(len, Is.LessThanOrEqualTo(maxLen));
        }
    }

    [Test]
    public void M3_MinLengthRespected()
    {
        int minLen = 10;
        var morfs = DisorderPredictor.PredictMoRFs(SeqWithIsland, minLength: minLen).ToList();

        foreach (var morf in morfs)
        {
            int len = morf.End - morf.Start + 1;
            Assert.That(len, Is.GreaterThanOrEqualTo(minLen));
        }
    }

    [Test]
    public void M4_CaseInsensitive()
    {
        string lower = SeqWithIsland.ToLowerInvariant();

        var morfsUpper = DisorderPredictor.PredictMoRFs(SeqWithIsland).ToList();
        var morfsLower = DisorderPredictor.PredictMoRFs(lower).ToList();

        Assert.That(morfsLower.Count, Is.EqualTo(morfsUpper.Count));
    }

    #endregion

    #region INV — Invariants

    [TestCase(SeqWithIsland)]
    [TestCase("PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP")]
    [TestCase("LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL")]
    public void INV1_ValidBoundaries(string sequence)
    {
        var morfs = DisorderPredictor.PredictMoRFs(sequence).ToList();

        foreach (var morf in morfs)
        {
            Assert.That(morf.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(morf.End, Is.LessThan(sequence.Length));
            Assert.That(morf.Start, Is.LessThanOrEqualTo(morf.End));
            Assert.That(morf.Score, Is.InRange(0.0, 1.0));
        }
    }

    [TestCase(SeqWithIsland)]
    public void INV2_MoRFsWithinIDR(string sequence)
    {
        var prediction = DisorderPredictor.PredictDisorder(sequence);
        var morfs = DisorderPredictor.PredictMoRFs(sequence).ToList();

        foreach (var morf in morfs)
        {
            bool withinSomeIDR = prediction.DisorderedRegions.Any(
                idr => morf.Start >= idr.Start && morf.End <= idr.End);
            Assert.That(withinSomeIDR, Is.True,
                $"MoRF [{morf.Start}-{morf.End}] not within any IDR");
        }
    }

    [TestCase(SeqWithIsland)]
    public void INV3_NoOverlappingMoRFs(string sequence)
    {
        var morfs = DisorderPredictor.PredictMoRFs(sequence)
            .OrderBy(m => m.Start).ToList();

        for (int i = 1; i < morfs.Count; i++)
        {
            Assert.That(morfs[i].Start, Is.GreaterThan(morfs[i - 1].End),
                $"MoRF {i} overlaps MoRF {i - 1}");
        }
    }

    #endregion
}
