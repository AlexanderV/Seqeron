// SEQ-REPLICATION-001 — finisher A1-2 / A1-3: all BA1F minimizers/maximizers, circular-genome option,
// Grigoriev (1998) windowed cumulative-skew prediction.
// Evidence: docs/Evidence/SEQ-REPLICATION-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-REPLICATION-001.md (A1..A5, C4..C8, W1..W3)
// Reference values: Rosalind BA1F sample + extra dataset (published answers); python brute force
// (per-base #G−#C walk, all argmin/argmax); numpy.cumsum(Bio.SeqUtils.GC_skew(seq, w)) (Biopython 1.88).

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GcSkewCalculator_ReplicationOriginExtensions_Tests
{
    private const string Ba1fSample =
        "CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG";

    private const string Ba1fExtraResource = "Seqeron.Genomics.Tests.TestData.Rosalind.ba1f_extra_dataset.txt";

    // Same 100 kb synthetic genome as GcSkewCalculator_PredictReplicationOrigin_Tests (total G−C = 0).
    private static string SyntheticGenome() =>
        string.Concat(Enumerable.Repeat("ACGTC", 6000))
        + string.Concat(Enumerable.Repeat("AGGTC", 10000))
        + string.Concat(Enumerable.Repeat("ACGTC", 4000));

    private static (string Genome, int[] Answer) LoadBa1fExtraDataset()
    {
        var asm = typeof(GcSkewCalculator_ReplicationOriginExtensions_Tests).Assembly;
        using Stream stream = asm.GetManifestResourceStream(Ba1fExtraResource)
            ?? throw new InvalidOperationException($"Embedded resource '{Ba1fExtraResource}' not found.");
        using var reader = new StreamReader(stream);
        string[] lines = reader.ReadToEnd().Split('\n').Select(l => l.Trim()).ToArray();
        // Layout: "Input" / genome / "Output" / answer.
        Assert.That(lines[0], Is.EqualTo("Input"));
        Assert.That(lines[2], Is.EqualTo("Output"));
        int[] answer = lines[3].Split(' ', StringSplitOptions.RemoveEmptyEntries).Select(int.Parse).ToArray();
        return (lines[1], answer);
    }

    private static string RotateLeft(string s, int r) => s[r..] + s[..r];

    #region A — all minimizers / maximizers (Rosalind BA1F)

    // A1 — BA1F sample: published output "53 97" (all minimizers, value −4).
    [Test]
    public void FindMinimumSkewPositions_Ba1fSample_ReturnsPublishedAnswer()
    {
        Assert.That(GcSkewCalculator.FindMinimumSkewPositions(new DnaSequence(Ba1fSample)), Is.EqualTo(new[] { 53, 97 }));
        Assert.That(GcSkewCalculator.FindMinimumSkewPositions(Ba1fSample), Is.EqualTo(new[] { 53, 97 }));
    }

    // A2 — BA1F extra dataset (93 523 bp): published answer "89969 89970 89971 90345 90346";
    // python brute force: min −184 at those indices, max +154 at [20377, 20378, 20379].
    [Test]
    public void FindMinimumSkewPositions_Ba1fExtraDataset_ReturnsPublishedAnswer()
    {
        var (genome, answer) = LoadBa1fExtraDataset();

        Assert.Multiple(() =>
        {
            Assert.That(genome, Has.Length.EqualTo(93523));
            Assert.That(answer, Is.EqualTo(new[] { 89969, 89970, 89971, 90345, 90346 }));
            Assert.That(GcSkewCalculator.FindMinimumSkewPositions(genome), Is.EqualTo(answer));
            Assert.That(GcSkewCalculator.FindMaximumSkewPositions(genome), Is.EqualTo(new[] { 20377, 20378, 20379 }));
            var p = GcSkewCalculator.PredictReplicationOrigin(genome);
            Assert.That(p.PredictedOrigin, Is.EqualTo(89969));
            Assert.That(p.OriginSkew, Is.EqualTo(-184.0));
            Assert.That(p.PredictedTerminus, Is.EqualTo(20377));
            Assert.That(p.TerminusSkew, Is.EqualTo(154.0));
        });
    }

    // A3 — All maximizers of the BA1F sample (python brute force: +2 at 16, 20, 21); the first equals
    // PredictedTerminus (16).
    [Test]
    public void FindMaximumSkewPositions_Ba1fSample_AllMaximizers()
    {
        Assert.That(GcSkewCalculator.FindMaximumSkewPositions(Ba1fSample), Is.EqualTo(new[] { 16, 20, 21 }));
        Assert.That(GcSkewCalculator.PredictReplicationOrigin(Ba1fSample).PredictedTerminus, Is.EqualTo(16));
    }

    // A4 — Small diagrams (python brute force): ties, endpoints, flat.
    [TestCase("CCGGCC", new[] { 2, 6 }, new[] { 0, 4 })]
    [TestCase("GGGCCC", new[] { 0, 6 }, new[] { 3 })]
    [TestCase("CGCG", new[] { 1, 3 }, new[] { 0, 2, 4 })]
    [TestCase("AATT", new[] { 0, 1, 2, 3, 4 }, new[] { 0, 1, 2, 3, 4 })]
    [TestCase("ccggcc", new[] { 2, 6 }, new[] { 0, 4 })]
    public void FindSkewPositions_SmallSequences_MatchBruteForce(string seq, int[] minimizers, int[] maximizers)
    {
        Assert.Multiple(() =>
        {
            Assert.That(GcSkewCalculator.FindMinimumSkewPositions(seq), Is.EqualTo(minimizers));
            Assert.That(GcSkewCalculator.FindMaximumSkewPositions(seq), Is.EqualTo(maximizers));
        });
    }

    // A5 — Input handling: null DnaSequence throws; null/empty string → [0] (diagram = Skew_0 only).
    [Test]
    public void FindSkewPositions_NullAndEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.FindMinimumSkewPositions((DnaSequence)null!));
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.FindMaximumSkewPositions((DnaSequence)null!));
            Assert.That(GcSkewCalculator.FindMinimumSkewPositions((string)null!), Is.EqualTo(new[] { 0 }));
            Assert.That(GcSkewCalculator.FindMaximumSkewPositions(string.Empty, circular: true), Is.EqualTo(new[] { 0 }));
        });
    }

    #endregion

    #region C — circular genomes

    // C4 — Circular mode identifies prefix n with 0: positions lie in [0, n−1] (python brute force over
    // Skew_0..Skew_{n−1}). GGGCCC: linear minimizers {0, 6} → circular {0}; G: linear max at 1 → circular 0.
    [TestCase("GGGCCC", new[] { 0 }, new[] { 3 })]
    [TestCase("CCGGCC", new[] { 2 }, new[] { 0, 4 })]
    [TestCase("CCGGGG", new[] { 2 }, new[] { 5 })]
    [TestCase("G", new[] { 0 }, new[] { 0 })]
    [TestCase("CGCG", new[] { 1, 3 }, new[] { 0, 2 })]
    public void FindSkewPositions_Circular_ReportsModuloN(string seq, int[] minimizers, int[] maximizers)
    {
        Assert.Multiple(() =>
        {
            Assert.That(GcSkewCalculator.FindMinimumSkewPositions(seq, circular: true), Is.EqualTo(minimizers));
            Assert.That(GcSkewCalculator.FindMaximumSkewPositions(seq, circular: true), Is.EqualTo(maximizers));
            Assert.That(GcSkewCalculator.PredictReplicationOrigin(seq, circular: true).PredictedOrigin, Is.EqualTo(minimizers[0]));
            Assert.That(GcSkewCalculator.PredictReplicationOrigin(seq, circular: true).PredictedTerminus, Is.EqualTo(maximizers[0]));
        });
    }

    // C5 — Total G−C = 0 ⇒ rotation-equivariance: rotating left by r shifts every minimizer/maximizer
    // to (p − r) mod n. Synthetic genome (python brute force): circular minimizers [29997, 30000, 30001]
    // (−6000), maximizers [79998, 79999] (+4001); rotation 12345 → [17652, 17655, 17656] / [67653, 67654].
    [TestCase(0)]
    [TestCase(12345)]
    [TestCase(50000)]
    [TestCase(99999)]
    public void FindSkewPositions_CircularBalancedGenome_RotationEquivariant(int r)
    {
        string g = SyntheticGenome();
        int n = g.Length;
        string rotated = RotateLeft(g, r);

        var mins = GcSkewCalculator.FindMinimumSkewPositions(rotated, circular: true);
        var maxs = GcSkewCalculator.FindMaximumSkewPositions(rotated, circular: true);

        Assert.Multiple(() =>
        {
            Assert.That(mins.Select(p => (p + r) % n).Order(), Is.EqualTo(new[] { 29997, 30000, 30001 }));
            Assert.That(maxs.Select(p => (p + r) % n).Order(), Is.EqualTo(new[] { 79998, 79999 }));
            Assert.That(mins.All(p => p >= 0 && p < n), Is.True);
        });
        if (r == 12345)
        {
            Assert.That(mins, Is.EqualTo(new[] { 17652, 17655, 17656 }));
            Assert.That(maxs, Is.EqualTo(new[] { 67653, 67654 }));
            var p = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(rotated), circular: true);
            Assert.That((p.PredictedOrigin, p.OriginSkew, p.PredictedTerminus, p.TerminusSkew),
                Is.EqualTo((17652, -3531.0, 67653, 6470.0)));
        }
    }

    // C6 — Total G−C ≠ 0 (documented): rotation is NOT equivariant. CCGGG (D = +1), python brute force:
    // circular minimizers {2}; rotated left by 3 ("GGCCG") gives {0, 4} → mapped back {3, 2} = {2, 3}.
    [Test]
    public void FindMinimumSkewPositions_CircularUnbalanced_RotationChangesSet()
    {
        Assert.That(GcSkewCalculator.FindMinimumSkewPositions("CCGGG", circular: true), Is.EqualTo(new[] { 2 }));
        var rotated = GcSkewCalculator.FindMinimumSkewPositions(RotateLeft("CCGGG", 3), circular: true);
        Assert.That(rotated, Is.EqualTo(new[] { 0, 4 }));
        Assert.That(rotated.Select(p => (p + 3) % 5).Order(), Is.EqualTo(new[] { 2, 3 }));
    }

    // C7 — circular = false is identical to the original overload; circular on BA1F (minimizers < n)
    // keeps 53/97; extra dataset circular = linear answer.
    [Test]
    public void PredictReplicationOrigin_CircularFalse_IdenticalToOriginal()
    {
        var (genome, answer) = LoadBa1fExtraDataset();
        Assert.Multiple(() =>
        {
            Assert.That(GcSkewCalculator.PredictReplicationOrigin(Ba1fSample, circular: false),
                Is.EqualTo(GcSkewCalculator.PredictReplicationOrigin(Ba1fSample)));
            Assert.That(GcSkewCalculator.PredictReplicationOrigin(new DnaSequence("CCGGGG"), circular: false),
                Is.EqualTo(GcSkewCalculator.PredictReplicationOrigin(new DnaSequence("CCGGGG"))));
            Assert.That(GcSkewCalculator.FindMinimumSkewPositions(Ba1fSample, circular: true), Is.EqualTo(new[] { 53, 97 }));
            Assert.That(GcSkewCalculator.FindMinimumSkewPositions(genome, circular: true), Is.EqualTo(answer));
        });
    }

    // C8 — Null/empty handling of the circular overload.
    [Test]
    public void PredictReplicationOrigin_Circular_NullAndEmpty()
    {
        Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.PredictReplicationOrigin((DnaSequence)null!, circular: true));
        Assert.That(GcSkewCalculator.PredictReplicationOrigin("", circular: true),
            Is.EqualTo(new ReplicationOriginPrediction(0, 0, 0, 0, false)));
    }

    #endregion

    #region W — windowed (Grigoriev 1998) prediction

    // W1 — numpy.cumsum(Bio.SeqUtils.GC_skew(seq, w)[:n//w]) argmin/argmax (first) → window centres.
    [TestCase("syn", 1000, 29500, -10.0, 79500, 6.666666666666664)]
    [TestCase("extra", 1000, 89500, -0.35203934153706795, 19500, 0.28817798587280524)]
    [TestCase("extra", 5000, 87500, -0.07055938354312184, 17500, 0.056616604707736475)]
    [TestCase("ba1f", 10, 45, -0.3095238095238095, 15, 0.5)]
    public void PredictReplicationOrigin_Windowed_MatchesBiopythonCumsum(
        string which, int w, int origin, double originSkew, int terminus, double terminusSkew)
    {
        string g = which switch
        {
            "syn" => SyntheticGenome(),
            "extra" => LoadBa1fExtraDataset().Genome,
            _ => Ba1fSample,
        };

        var p = GcSkewCalculator.PredictReplicationOrigin(new DnaSequence(g), w);

        Assert.Multiple(() =>
        {
            Assert.That(p.PredictedOrigin, Is.EqualTo(origin));
            Assert.That(p.OriginSkew, Is.EqualTo(originSkew).Within(1e-12));
            Assert.That(p.PredictedTerminus, Is.EqualTo(terminus));
            Assert.That(p.TerminusSkew, Is.EqualTo(terminusSkew).Within(1e-12));
            Assert.That(p.IsSignificant, Is.True);
        });
    }

    // W2 — Window 1: values are the per-nucleotide diagram without Skew_0 (GGGCCC: 1,2,3,2,1,0 →
    // min 0 at the last base, centre 5; max 3 at base 2).
    [Test]
    public void PredictReplicationOrigin_WindowOne_NoZeroBaseline()
    {
        var p = GcSkewCalculator.PredictReplicationOrigin("GGGCCC", 1);
        Assert.That((p.PredictedOrigin, p.OriginSkew, p.PredictedTerminus, p.TerminusSkew, p.IsSignificant),
            Is.EqualTo((5, 0.0, 2, 3.0, true)));
    }

    // W3 — Guards and degenerate input.
    [Test]
    public void PredictReplicationOrigin_Windowed_GuardsAndShortInput()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.PredictReplicationOrigin("ACGT", 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => GcSkewCalculator.PredictReplicationOrigin(new DnaSequence("ACGT"), -1));
            Assert.Throws<ArgumentNullException>(() => GcSkewCalculator.PredictReplicationOrigin((DnaSequence)null!, 10));
            Assert.That(GcSkewCalculator.PredictReplicationOrigin("ACG", 10), Is.EqualTo(new ReplicationOriginPrediction(0, 0, 0, 0, false)));
            Assert.That(GcSkewCalculator.PredictReplicationOrigin((string)null!, 10), Is.EqualTo(new ReplicationOriginPrediction(0, 0, 0, 0, false)));
        });
    }

    #endregion
}
