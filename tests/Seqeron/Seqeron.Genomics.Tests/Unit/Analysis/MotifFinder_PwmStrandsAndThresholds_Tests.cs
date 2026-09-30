namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// PAT-PWM-001 (review 2026-09, B05 follow-up): both-strand PWM search, per-window scores, count-matrix
/// PWMs with JASPAR pseudocounts, and the score distribution / thresholds.
/// All expected values are Biopython 1.88 outputs (Bio/motifs/matrix.py <c>calculate</c>, <c>search(both=True)</c>,
/// <c>mean</c>, <c>std</c>, <c>distribution</c>; Bio/motifs/thresholds.py <c>ScoreDistribution</c>;
/// Bio/motifs/jaspar <c>read</c> + <c>calculate_pseudocounts</c>). Biopython scores windows in float32 → 1e-5.
/// </summary>
[TestFixture]
[Category("PAT-PWM-001")]
public class MotifFinder_PwmStrandsAndThresholds_Tests
{
    private static readonly string[] WikipediaSequences =
    {
        "GAGGTAAAC", "TCCGTAAGT", "CAGGTTGGA", "ACAGTCAGT", "TAGGTCATT",
        "TAGGTACTG", "ATGGTAACT", "CAGGTATAC", "TGTGTGAGT", "AAGGTAAGT"
    };

    private const string WikipediaTarget = "CCTAGGTAAGTAACAGGTCAGTGG";

    private static PositionWeightMatrix WikipediaPwm() => MotifFinder.CreatePwm(WikipediaSequences, 0.25);

    private static void AssertHits(IReadOnlyList<PwmStrandMatch> actual, (int BioPos, double Score)[] expected, int n)
    {
        Assert.That(actual.Select(h => h.BiopythonPosition), Is.EqualTo(expected.Select(e => e.BioPos)), "Biopython positions");
        Assert.That(actual.Select(h => h.Score), Is.EqualTo(expected.Select(e => e.Score)).Within(1e-5), "scores");
        foreach (var h in actual)
        {
            Assert.That(h.Strand, Is.EqualTo(h.BiopythonPosition >= 0 ? '+' : '-'));
            Assert.That(h.Position, Is.EqualTo(h.Strand == '+' ? h.BiopythonPosition : h.BiopythonPosition + n));
        }
    }

    #region ScanWithPwmBothStrands (Biopython search both=True)

    [Test]
    [Description("Wikipedia example: search(target, 0.0, both=True) = [(2,11.7075),(6,4.7792),(-16,2.3877),(13,9.3161)]")]
    public void ScanWithPwmBothStrands_WikipediaExample_EqualsBiopython()
    {
        var hits = MotifFinder.ScanWithPwmBothStrands(new DnaSequence(WikipediaTarget), WikipediaPwm(), 0.0).ToList();

        AssertHits(hits, new[] { (2, 11.70753002166748), (6, 4.779160022735596), (-16, 2.387691020965576), (13, 9.316061019897461) },
            WikipediaTarget.Length);
        Assert.Multiple(() =>
        {
            Assert.That(hits.Select(h => h.Position), Is.EqualTo(new[] { 2, 6, 8, 13 }), "forward window starts, ascending");
            Assert.That(hits[2].MatchedSequence, Is.EqualTo("CCTGTTACT"), "minus site read 5'→3' = revcomp(AGTAACAGG)");
            Assert.That(hits.All(h => h.Pattern == "TAGGTAAGT"), Is.True);
        });
    }

    [Test]
    [Description("Threshold −3: Biopython adds minus hits at -21 (-0.79252) and -17 (-0.64793)")]
    public void ScanWithPwmBothStrands_LowerThreshold_EqualsBiopython()
    {
        var hits = MotifFinder.ScanWithPwmBothStrands(new DnaSequence(WikipediaTarget), WikipediaPwm(), -3.0).ToList();
        AssertHits(hits, new[]
        {
            (2, 11.70753002166748), (-21, -0.7925168871879578), (6, 4.779160022735596),
            (-17, -0.6479330062866211), (-16, 2.387691020965576), (13, 9.316061019897461)
        }, WikipediaTarget.Length);
    }

    [Test]
    [Description("Palindromic PWM (ACGT/ACGT/TCGA, rc == pwm): both strands of the same window are reported, '+' first")]
    public void ScanWithPwmBothStrands_PalindromicPwm_ReportsBothStrandsPlusFirst()
    {
        var pwm = MotifFinder.CreatePwm(new[] { "ACGT", "ACGT", "TCGA" }, 0.25);
        const string s = "GGACGTCCTCGAGG";
        var hits = MotifFinder.ScanWithPwmBothStrands(new DnaSequence(s), pwm, 0.0).ToList();

        // Biopython: [(2,5.74073),(-12,5.74073),(5,1.19229),(-9,1.19229),(8,4.04474),(-6,4.04474)]
        AssertHits(hits, new[]
        {
            (2, 5.740729331970215), (-12, 5.740729331970215), (5, 1.1922928094863892),
            (-9, 1.1922928094863892), (8, 4.044735431671143), (-6, 4.044735431671143)
        }, s.Length);
    }

    private static IEnumerable<TestCaseData> RandomBiopythonCases()
    {
        // random.seed(20260930); L = randint(5,8); 6 instances; 50-nt sequence; normalize(0.5); threshold 1.0
        yield return new TestCaseData(
            new[] { "TCGTAAGT", "TAATACGA", "TGATACAG", "TACCGGTG", "CCATATGG", "AGCGGGCA" },
            "CGGTTGGTAGGAAGCGTTAAGAGGAACAATTTTCGGAGGGAAGATCGCTA",
            new[] { 4, -46, -39, -36, 17, -21, 32 },
            new[] { 4.8678765296936035, 4.382450103759766, 1.6979516744613647, 1.8089829683303833, 4.505306720733643, 2.0605218410491943, 3.768341064453125 })
            .SetName("RandomSeed20260930_Case0");
        yield return new TestCaseData(
            new[] { "AGATT", "CGATA", "TACGA", "AGACA", "CAGCT", "TGCTT" },
            "CCATATTCCAATATGGCTTCATTTAGTCCCATATTCGACCGTGCCACGTG",
            new[] { 8, 13, -37, 14, -26, -24, -18, 41, -9, 44, -5 },
            new[] { 3.0659210681915283, 2.2060985565185547, 2.9430642127990723, 1.1065629720687866, 1.1065629720687866, 1.4691330194473267, 1.1065629720687866, 2.9430642127990723, 2.2060985565185547, 1.358101725578308, 1.358101725578308 })
            .SetName("RandomSeed20260930_Case1");
        yield return new TestCaseData(
            new[] { "GATTAG", "CTCGCT", "TCAGAA", "TTAGGT", "GACCCC", "AATACT" },
            "TCACAAGGTGTCGTCAGCTCCCTTTGATCAACCGCAGGACGGCGACTTGA",
            new[] { 3, -44, -41, 13, -35, 17, 22, -24, 37, -10 },
            new[] { 1.913917899131775, 1.4284909963607788, 2.165456533432007, 2.650883436203003, 1.0659209489822388, 1.4284909963607788, 2.9024221897125244, 2.9024221897125244, 1.4284909963607788, 2.165456533432007 })
            .SetName("RandomSeed20260930_Case2");
    }

    [TestCaseSource(nameof(RandomBiopythonCases))]
    [Description("Seeded random motifs/sequences: both-strand hits equal Biopython search(seq, 1.0, both=True)")]
    public void ScanWithPwmBothStrands_RandomCases_EqualBiopython(string[] instances, string sequence, int[] bioPositions, double[] scores)
    {
        var pwm = MotifFinder.CreatePwm(instances, 0.5);
        var hits = MotifFinder.ScanWithPwmBothStrands(new DnaSequence(sequence), pwm, 1.0).ToList();
        AssertHits(hits, bioPositions.Zip(scores).ToArray(), sequence.Length);
    }

    [Test]
    public void ScanWithPwmBothStrands_PlusStrandSubsetEqualsScanWithPwm()
    {
        var seq = new DnaSequence(WikipediaTarget);
        var pwm = WikipediaPwm();
        var plus = MotifFinder.ScanWithPwmBothStrands(seq, pwm, -5).Where(h => h.Strand == '+')
            .Select(h => (h.Position, h.MatchedSequence, h.Score));
        var forward = MotifFinder.ScanWithPwm(seq, pwm, -5).Select(m => (m.Position, m.MatchedSequence, m.Score));
        Assert.That(plus, Is.EqualTo(forward));
    }

    [Test]
    public void ScanWithPwmBothStrands_NullArgumentsAndShortSequence()
    {
        var pwm = WikipediaPwm();
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MotifFinder.ScanWithPwmBothStrands(null!, pwm));
            Assert.Throws<ArgumentNullException>(() => MotifFinder.ScanWithPwmBothStrands(new DnaSequence("ACGT"), null!));
            Assert.That(MotifFinder.ScanWithPwmBothStrands(new DnaSequence("ACGT"), pwm, double.NegativeInfinity), Is.Empty);
        });
    }

    #endregion

    #region CalculatePwmScores (Biopython calculate)

    [Test]
    [Description("calculate('ACGTNACgtacg') = [5.74073, NaN×4, 5.74073, -5.67807, -8, -5.67807] (N windows → NaN, mixed case)")]
    public void CalculatePwmScores_InvalidWindowsNaN_MixedCase_EqualsBiopython()
    {
        var pwm = MotifFinder.CreatePwm(new[] { "ACGT", "ACGA", "TCGT" }, 0.25);
        double[] scores = MotifFinder.CalculatePwmScores("ACGTNACgtacg", pwm);

        double[] expected = { 5.740729331970215, double.NaN, double.NaN, double.NaN, double.NaN, 5.740729331970215, -5.678071975708008, -8.0, -5.678071975708008 };
        Assert.That(scores, Has.Length.EqualTo(expected.Length));
        for (int i = 0; i < expected.Length; i++)
        {
            if (double.IsNaN(expected[i]))
                Assert.That(scores[i], Is.NaN, $"window {i}");
            else
                Assert.That(scores[i], Is.EqualTo(expected[i]).Within(1e-5), $"window {i}");
        }
    }

    [Test]
    public void CalculatePwmScores_DnaSequence_EqualsBiopythonCalculate_AndShortSequenceEmpty()
    {
        var pwm = WikipediaPwm();
        double[] expected =
        {
            -12.337397575378418, -5.2917890548706055, 11.70753002166748, -7.796898365020752,
            -14.033390998840332, -10.649341583251953, 4.779160022735596, -5.47497034072876,
            -8.048437118530273, -8.167141914367676, -9.097930908203125, -11.180948257446289,
            -6.4493303298950195, 9.316061019897461, -7.796898365020752, -12.810998916625977
        };
        Assert.Multiple(() =>
        {
            Assert.That(MotifFinder.CalculatePwmScores(new DnaSequence(WikipediaTarget), pwm), Is.EqualTo(expected).Within(1e-5));
            Assert.That(MotifFinder.CalculatePwmScores("ACGT", pwm), Is.Empty);
            Assert.Throws<ArgumentNullException>(() => MotifFinder.CalculatePwmScores((string)null!, pwm));
        });
    }

    #endregion

    #region Mean / Std / ScoreDistribution (Biopython thresholds.py)

    [Test]
    [Description("pssm.mean() = 5.522241422369563, pssm.std() = 3.2186083897634568, mean(bg .3/.2/.2/.3) = 5.975137019521005")]
    public void MeanStd_WikipediaExample_EqualBiopython()
    {
        var pwm = WikipediaPwm();
        Assert.Multiple(() =>
        {
            Assert.That(pwm.Mean(), Is.EqualTo(5.522241422369563).Within(1e-12));
            Assert.That(pwm.Std(), Is.EqualTo(3.2186083897634568).Within(1e-12));
            Assert.That(pwm.Mean(new[] { 0.3, 0.2, 0.2, 0.3 }), Is.EqualTo(5.975137019521005).Within(1e-12));
        });
    }

    [Test]
    [Description("distribution(): min −14.88138790352414, step 0.002954652535685415, 9000 points; thresholds equal Biopython")]
    public void ScoreDistribution_WikipediaExample_ThresholdsEqualBiopython()
    {
        var d = WikipediaPwm().ScoreDistribution();
        Assert.Multiple(() =>
        {
            Assert.That(d.MinScore, Is.EqualTo(-14.88138790352414).Within(1e-12));
            Assert.That(d.Step, Is.EqualTo(0.002954652535685415).Within(1e-15));
            Assert.That(d.PointCount, Is.EqualTo(9000));
            Assert.That(d.ThresholdFpr(0.01), Is.EqualTo(4.028388324862519).Within(1e-9));
            Assert.That(d.ThresholdFpr(0.001), Is.EqualTo(7.10122696197535).Within(1e-9));
            Assert.That(d.ThresholdFnr(0.1), Is.EqualTo(1.2303323735684302).Within(1e-9));
            Assert.That(d.ThresholdBalanced(), Is.EqualTo(0.1430202404361971).Within(1e-9));
            Assert.That(d.ThresholdBalanced(1.0, out double rate), Is.EqualTo(0.1430202404361971).Within(1e-9));
            Assert.That(rate, Is.EqualTo(0.06385040283203125).Within(1e-12));
            Assert.That(d.ThresholdPatser(), Is.EqualTo(2.5924271925194056).Within(1e-9));
            Assert.That(d.BackgroundDensity.Sum(), Is.EqualTo(1.0).Within(1e-9));
            Assert.That(d.MotifDensity.Sum(), Is.EqualTo(1.0).Within(1e-6));
        });
    }

    [Test]
    [Description("distribution(background={A:.3,C:.2,G:.2,T:.3}, precision=100) thresholds equal Biopython")]
    public void ScoreDistribution_BackgroundAndPrecision_EqualBiopython()
    {
        var d = WikipediaPwm().ScoreDistribution(new[] { 0.3, 0.2, 0.2, 0.3 }, precision: 100);
        Assert.Multiple(() =>
        {
            Assert.That(d.ThresholdFpr(0.01), Is.EqualTo(4.254351868562161).Within(1e-9));
            Assert.That(d.ThresholdFnr(0.05), Is.EqualTo(-0.4482487864018605).Within(1e-9));
            Assert.That(d.ThresholdBalanced(2.0), Is.EqualTo(3.3079165166197164).Within(1e-9));
            Assert.That(d.ThresholdPatser(), Is.EqualTo(3.4262209356125233).Within(1e-9));
            Assert.That(d.MeanScore, Is.EqualTo(5.975137019521005).Within(1e-12));
        });
    }

    [Test]
    public void ScoreDistribution_InvalidInputs_Throw()
    {
        var pwm = WikipediaPwm();
        var d = pwm.ScoreDistribution(precision: 10);
        Assert.Multiple(() =>
        {
            Assert.Throws<InvalidOperationException>(() => MotifFinder.CreatePwm(WikipediaSequences, 0.0).ScoreDistribution());
            Assert.Throws<ArgumentOutOfRangeException>(() => pwm.ScoreDistribution(precision: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => d.ThresholdFpr(-0.1));
            Assert.Throws<ArgumentOutOfRangeException>(() => d.ThresholdFnr(1.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => d.ThresholdFpr(double.NaN));
            Assert.Throws<ArgumentOutOfRangeException>(() => d.ThresholdBalanced(-1));
            Assert.That(d.ThresholdFpr(0.0), Is.EqualTo(d.MinScore + d.PointCount * d.Step).Within(1e-12), "fpr 0: one step above the grid (Biopython)");
        });
    }

    [TestCase(7.0, 2.0)]
    [TestCase(-7.0, 2.0)]
    [TestCase(7.0, -2.0)]
    [TestCase(-7.0, -2.0)]
    [TestCase(0.30000000000000004, 0.1)]
    [TestCase(-0.0, 3.0)]
    [TestCase(5.5, 0.5)]
    [Description("PythonFloorDiv reproduces CPython float //: 7//2=3, -7//2=-4, 7//-2=-4, -7//-2=3, 0.30000000000000004//0.1=3, -0.0//3=-0.0, 5.5//0.5=11")]
    public void PythonFloorDiv_MatchesCPython(double a, double b)
    {
        double expected = (a, b) switch
        {
            (7.0, 2.0) => 3, (-7.0, 2.0) => -4, (7.0, -2.0) => -4, (-7.0, -2.0) => 3,
            (0.30000000000000004, 0.1) => 3, (5.5, 0.5) => 11, _ => -0.0
        };
        Assert.That(PwmScoreDistribution.PythonFloorDiv(a, b), Is.EqualTo(expected));
    }

    #endregion

    #region FromCounts / JASPAR (Bucher 1990 matrices)

    [Test]
    [Description("JASPAR POL012.1: Biopython jaspar.calculate_pseudocounts = √389·0.25 = 4.926543751285817; full PSSM equals m.pssm")]
    public void FromCounts_JasparTataBox_EqualsBiopythonJasparPssm()
    {
        var tata = MotifFinder.BucherPromoterMatrices.TataBox;
        double[] pc = MotifFinder.JasparPseudocounts(tata.GetCounts());
        Assert.That(pc, Is.EqualTo(Enumerable.Repeat(4.926543751285817, 4)).Within(1e-12));

        double[,] expected =
        {
            { -0.6321326866402173, -2.2876580812831695, 1.8045631543459546, -3.688228292665394, 1.8126245782321924, 1.4174486715023265, 1.8722801429968936, 1.1511613311742428, 0.6463453687677879, -0.7459212343193731, -0.21669337572902317, -0.23319535222832682, -0.23319535222832702, -0.48654811435798495, -0.3186611827155647 },
            { 0.5531917871826982, -1.0045743269155591, -4.374344361085418, -2.775111988228752, -4.374344361085418, -4.374344361085418, -3.652490201348702, -3.8827845832880707, -1.0623747635055731, 0.45360562375401703, 0.5723099107046337, 0.36867082474523644, 0.26673243238218924, 0.13148817608621427, 0.052000113669752754 },
            { 0.6190253606574071, -2.155973255121007, -3.8827845832880707, -3.8827845832880707, -3.3636287464633328, -4.374344361085418, -2.73937389691206, -1.0623747635055731, 0.6642754575783802, 0.6005203032621591, 0.37956518063223604, 0.37956518063223604, 0.3795651806322358, 0.4942686459151711, 0.5042578124986132 },
            { -1.5079419854817468, 1.6193629790220367, -1.355643945397777, 1.8908541632509208, -1.548668251523926, 0.30151837658584674, -3.1894169197158244, 0.30151837658584674, -1.4297842308908504, -0.9490006905092385, -1.5079419854817468, -0.8438906252740562, -0.6321326866402175, -0.354317431218507, -0.428387798798603 },
        };
        var pwm = tata.Pwm;
        Assert.Multiple(() =>
        {
            for (int b = 0; b < 4; b++)
                for (int j = 0; j < 15; j++)
                    Assert.That(pwm.Matrix[b, j], Is.EqualTo(expected[b, j]).Within(1e-12), $"[{b},{j}]");
            Assert.That(pwm.MaxScore, Is.EqualTo(15.782082872405407).Within(1e-12));
            Assert.That(pwm.MinScore, Is.EqualTo(-34.61437864989872).Within(1e-12));
            Assert.That(pwm.Consensus, Is.EqualTo("GTATAAAAGGCGGGG"));
        });
    }

    private static IEnumerable<TestCaseData> BucherMatrixCases()
    {
        // name, JASPAR id, pseudocount, max, min, consensus, mean, fpr 1e-3, fpr 1e-4, patser, balanced (Biopython)
        yield return new TestCaseData("TATA Box", "POL012.1", 4.926543751285817, 15.782082872405407, -34.61437864989872, "GTATAAAAGGCGGGG", 8.967593531566934, 7.042753822501211, 10.506901554316364, 5.799558225244553, 0.269017811502799).SetName("POL012.1_TataBox");
        yield return new TestCaseData("Cap Signal", "POL002.1", 4.351723796382303, 7.359316438826737, -16.54820669256136, "TCAGTCTT", 3.97013071484069, 6.115969817666912, 6.863173315960076, 0.482055440536449, -0.06788633420731927).SetName("POL002.1_CapSignal");
        yield return new TestCaseData("CCAAT Box", "POL004.1", 3.307189138830738, 14.741533220189346, -32.33645571322996, "ACTAGCCAATCA", 8.857384140990742, 7.094645943646782, 10.53554904345583, 5.68218857086314, 0.27169213456134855).SetName("POL004.1_CcaatBox");
        yield return new TestCaseData("GC Box", "POL003.1", 4.138236339311712, 17.58949060499476, -41.790405797972014, "AGGGGGCGGGGCTG", 9.808703963688972, 6.836719879834341, 10.522778167867635, 6.633117120541364, 0.30446468585127917).SetName("POL003.1_GcBox");
    }

    [TestCaseSource(nameof(BucherMatrixCases))]
    [Description("Each Bucher/JASPAR matrix: pseudocounts, pssm max/min/consensus/mean and distribution thresholds equal Biopython")]
    public void BucherMatrix_PssmAndThresholds_EqualBiopython(string name, string id, double pseudo, double max, double min,
        string consensus, double mean, double fpr3, double fpr4, double patser, double balanced)
    {
        var e = MotifFinder.BucherPromoterMatrices.All.Single(x => x.MatrixId == id);
        var d = e.Pwm.ScoreDistribution();
        Assert.Multiple(() =>
        {
            Assert.That(e.Name, Is.EqualTo(name));
            Assert.That(MotifFinder.JasparPseudocounts(e.GetCounts())[0], Is.EqualTo(pseudo).Within(1e-12));
            Assert.That(e.Pwm.MaxScore, Is.EqualTo(max).Within(1e-11));
            Assert.That(e.Pwm.MinScore, Is.EqualTo(min).Within(1e-11));
            Assert.That(e.Pwm.Consensus, Is.EqualTo(consensus));
            Assert.That(e.Pwm.Mean(), Is.EqualTo(mean).Within(1e-11));
            Assert.That(d.ThresholdFpr(1e-3), Is.EqualTo(fpr3).Within(1e-9));
            Assert.That(d.ThresholdFpr(1e-4), Is.EqualTo(fpr4).Within(1e-9));
            Assert.That(d.ThresholdPatser(), Is.EqualTo(patser).Within(1e-9));
            Assert.That(d.ThresholdBalanced(), Is.EqualTo(balanced).Within(1e-9));
        });
    }

    // Adenovirus-major-late-like test promoter with an added CCAAT/GC-box block (reference sequence of the lock).
    private const string PromoterSequence =
        "GGGGCTATAAAAGGGGGTGGGGGCGCGTTCGTCCTCACTCTCTTCCGCATCGCTGTCTGCGAGGGCCAGCCAATCAGCGCCCCGCCCATTGGCTGGGCGGAGCC";

    [Test]
    [Description("FindPromoterElementsByMatrix(fpr 1e-3) = Biopython pssm.search at threshold_fpr(1e-3): TATA +4, cap +34, CCAAT +64/−(−19), GC +11/+17/−(−28)")]
    public void FindPromoterElementsByMatrix_EqualsBiopythonSearch()
    {
        var hits = MotifFinder.FindPromoterElementsByMatrix(new DnaSequence(PromoterSequence), 1e-3).ToList();
        int n = PromoterSequence.Length;
        var expected = new (string Id, int Pos, char Strand, double Score)[]
        {
            ("POL012.1", 4, '+', 14.674918174743652),
            ("POL002.1", 34, '+', 6.793403148651123),
            ("POL004.1", 64, '+', 12.077281951904297),
            ("POL004.1", n - 19, '-', 12.315014839172363),
            ("POL003.1", 11, '+', 13.369436264038086),
            ("POL003.1", 17, '+', 8.659268379211426),
            ("POL003.1", n - 28, '-', 13.07020378112793),
        };
        Assert.That(hits.Select(h => (h.MatrixId, h.Position, h.Strand)), Is.EqualTo(expected.Select(e => (e.Id, e.Pos, e.Strand))));
        Assert.That(hits.Select(h => h.Score), Is.EqualTo(expected.Select(e => e.Score)).Within(1e-5));
        Assert.That(hits.First(h => h.MatrixId == "POL012.1").Threshold, Is.EqualTo(7.042753822501211).Within(1e-9));

        // Strand-specific TATA / cap: given strand only even with bothStrands (Biopython both=False lists equal).
        var single = MotifFinder.FindPromoterElementsByMatrix(new DnaSequence(PromoterSequence), 1e-3, bothStrands: false).ToList();
        Assert.That(single.All(h => h.Strand == '+'), Is.True);
        Assert.That(single.Count, Is.EqualTo(5));
    }

    [Test]
    public void FromCounts_ValidationAndScalarPseudocount()
    {
        var counts = new double[,] { { 3, 0 }, { 0, 1 }, { 0, 1 }, { 0, 1 } };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => PositionWeightMatrix.FromCounts(null!, 0.5));
            Assert.Throws<ArgumentException>(() => PositionWeightMatrix.FromCounts(new double[3, 2], 0.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => PositionWeightMatrix.FromCounts(new double[,] { { -1 }, { 1 }, { 1 }, { 1 } }, 0.5));
            Assert.Throws<ArgumentOutOfRangeException>(() => PositionWeightMatrix.FromCounts(counts, -0.5));
            Assert.Throws<ArgumentException>(() => PositionWeightMatrix.FromCounts(new double[4, 1], 0.0), "empty column, no pseudocount");
            Assert.Throws<ArgumentException>(() => PositionWeightMatrix.FromCounts(counts, new[] { 0.1, 0.1 }));
            Assert.Throws<ArgumentException>(() => MotifFinder.JasparPseudocounts(new double[4, 0]));
        });

        // Counts of ACG×3/AAG... equal CreatePwm of the underlying sequences (same kernel).
        var seqs = new[] { "AC", "AG", "AT" };
        var fromSeqs = MotifFinder.CreatePwm(seqs, 0.25);
        var fromCounts = PositionWeightMatrix.FromCounts(counts, 0.25);
        Assert.That(fromCounts.Matrix, Is.EqualTo(fromSeqs.Matrix));
    }

    [Test]
    public void FindPromoterElementsByMatrix_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => MotifFinder.FindPromoterElementsByMatrix(null!, 1e-3));
            Assert.Throws<ArgumentOutOfRangeException>(() => MotifFinder.FindPromoterElementsByMatrix(new DnaSequence("ACGT"), 1.5));
        });
    }

    #endregion
}
