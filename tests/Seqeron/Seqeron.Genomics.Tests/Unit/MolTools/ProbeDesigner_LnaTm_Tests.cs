// PROBE-LNATM-001 — LNA-modified nearest-neighbour thermodynamics / Tm + citable MGB design rules
// Evidence: docs/Evidence/PROBE-DESIGN-001-LNA-Evidence.md
// TestSpec: tests/TestSpecs/PROBE-LNATM-001.md
// Sources:
//   McTigue PM, Peterson RJ, Kahn JD (2004). Biochemistry 43:5388-5405 — LNA-DNA NN increments.
//   Owczarzy R, You Y, Groth CL, Tataurov AV (2011). Biochemistry 50:9352-9367 — single-LNA, consecutive-LNA
//     and LNA-mismatch NN parameters; measured Tm (2 µM, 1 M Na⁺) of the LNA triplet duplexes below.
//   MELTING 5.2.0 (Dumousseau et al. 2012, BMC Bioinformatics 13:101; melting5.jar + data files shipped in
//     Bioconductor rmelting) — reference implementation; every expected ΔH°/ΔS°/Tm below was produced by
//     melting5.jar (Main.getMeltingResults, full double precision), e.g. rmelting test-method.locked.R:
//     CCATTLGCTACC mct04 63.61426, owc11 63.48299; GALCLC 12.94323.
//   Biopython 1.88 Bio.SeqUtils.MeltingTemp.Tm_NN(nn_table=DNA_NN3) — LNA-free reduction.
//   Kutyavin IV et al. (2000). Nucleic Acids Res 28(2):655-661 — 3'-MGB design rules.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Tests for the LNA-modified nearest-neighbour thermodynamics / Tm (McTigue 2004 and Owczarzy 2011 models on the
/// SantaLucia 1998 unified DNA set, as MELTING 5 implements them) and the citable 3'-MGB design-rule check.
/// </summary>
[TestFixture]
public class ProbeDesigner_LnaTm_Tests
{
    private const double Tol = 1e-9;
    private const double MeltingR = 1.99;   // MELTING NearestNeighborMode.computesMeltingTemperature
    private const PrimerDesigner.LnaNearestNeighborModel Owc = PrimerDesigner.LnaNearestNeighborModel.Owczarzy2011;
    private const PrimerDesigner.LnaNearestNeighborModel Mct = PrimerDesigner.LnaNearestNeighborModel.McTigue2004;

    // MELTING notation: "L" after a base marks it as LNA.
    private static (string Seq, int[] Lna) Parse(string melting)
    {
        var sb = new System.Text.StringBuilder();
        var lna = new List<int>();
        foreach (char c in melting)
        {
            if (c == 'L') lna.Add(sb.Length - 1);
            else sb.Append(c);
        }
        return (sb.ToString(), lna.ToArray());
    }

    private static void AssertMelting(string melting, string? target, PrimerDesigner.LnaNearestNeighborModel model,
        double c, double na, double mg, double hCal, double s, double tm)
    {
        var (seq, lna) = Parse(melting);
        var th = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, lna, model, target);
        var mode = mg > 0 ? PrimerDesigner.SaltCorrectionMode.Owczarzy2008Divalent : PrimerDesigner.SaltCorrectionMode.Owczarzy2004Monovalent;
        double t = PrimerDesigner.CalculateMeltingTemperatureNNLna(seq, lna, model, target, c, na, mg, 0, mode, MeltingR);
        Assert.That(th, Is.Not.Null, melting);
        Assert.Multiple(() =>
        {
            Assert.That(th!.Value.DeltaH * 1000, Is.EqualTo(hCal).Within(1e-6), $"{melting} ΔH° (MELTING)");
            Assert.That(th.Value.DeltaS, Is.EqualTo(s).Within(Tol), $"{melting} ΔS° (MELTING)");
            Assert.That(t, Is.EqualTo(tm).Within(1e-9), $"{melting} Tm (MELTING)");
        });
    }

    #region rmelting worked examples (test-method.locked.R) — C_T 1e-4 M, Na⁺ 1 M

    [Test]
    public void WorkedExample_McTigue2004_CCATTLGCTACC_MatchesMelting()
    {
        // MELTING -lck mct04: base all97 (−81.1/−222.5) + TTL/AA (+2.326/+8.1) + TLG/AC (−1.540/−3.0).
        AssertMelting("CCATTLGCTACC", null, Mct, 1e-4, 1, 0, -80314.0, -217.4, 63.614259353345176);
    }

    [Test]
    public void WorkedExample_Owczarzy2011_CCATTLGCTACC_MatchesMelting()
    {
        // MELTING default owc11: complete TTL/AA (−5.574/−14.149) and TLG/AC (−10.040/−25.744).
        AssertMelting("CCATTLGCTACC", null, Owc, 1e-4, 1, 0, -80314.0, -217.493, 63.48298667194416);
    }

    [Test]
    public void WorkedExample_Owczarzy2011_ConsecutiveLna_GALCLC_MatchesMelting()
    {
        AssertMelting("GALCLC", null, Owc, 1e-4, 1, 0, -24849.0, -65.76899999999999, 12.943226486909339);
    }

    [Test]
    public void DefaultOverloads_AreOwczarzy2011()
    {
        var (seq, lna) = Parse("CCATTLGCTACC");
        var a = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, lna);
        var b = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, lna, Owc);
        double ta = PrimerDesigner.CalculateMeltingTemperatureNNLna(seq, lna, 1e-4, 1.0,
            saltMode: PrimerDesigner.SaltCorrectionMode.None);
        // R = 1.9872 (class default): −80314 / (−217.493 + 1.9872·ln(2.5e-5)) − 273.15
        double expected = -80314.0 / (-217.493 + 1.9872 * Math.Log(1e-4 / 4)) - 273.15;
        Assert.Multiple(() =>
        {
            Assert.That(a, Is.EqualTo(b));
            Assert.That(ta, Is.EqualTo(expected).Within(1e-9));
        });
    }

    #endregion

    #region Owczarzy et al. (2011) LNA triplet duplexes — 2 µM, 1 M Na⁺ (MELTING test set)

    // (probe in MELTING notation, target 3'→5', MELTING ΔH cal/mol, ΔS, Tm, measured Tm (Owczarzy 2011))
    private static readonly object[] TripletCases =
    {
        new object[] { "TGACGGAGLCLGLATTCAGC", "ACTGCCTCGCTAAGTCG", -143605.0, -378.65899999999993, 79.2279033762257, 78.6 },
        new object[] { "CTATCCAGLGLCLATTCGCA", "GATAGGTCCGTAAGCGT", -142282.0, -376.9169999999999, 77.48030230203858, 77.4 },
        new object[] { "TTACTGTCLALALGGCAACT", "AATGACAGTTCCGTTGA", -133295.0, -356.491, 72.74444452341862, 73.0 },
        new object[] { "GCGTCAAGLCLGLACATCAT", "CGCAGTTCGCTGTAGTA", -143705.0, -381.0589999999999, 77.40880076697937, 75.0 },
        new object[] { "CGACTTGTLCLCLATACCTA", "GCTGAACAGGTATGGAT", -135851.0, -361.9319999999999, 74.46906334150503, 73.7 },
        new object[] { "CCATGCGTLALGLACAAGTG", "GGTACGCATCTGTTCAC", -139084.0, -368.41299999999995, 76.9360084459031, 74.4 },
        new object[] { "CTATCGCALTLCLTAATAAT", "GATAGCGTAGATTATTA", -131103.0, -360.6429999999999, 63.4299082370959, 63.8 },
        // central +X·mismatch (LNA-mismatch parameters)
        new object[] { "TGACGGAGLCLGLATTCAGC", "ACTGCCTCACTAAGTCG", -135441.0, -365.7269999999999, 70.08685927536197, 68.5 },
        new object[] { "CTATCCAGLGLCLATTCGCA", "GATAGGTCTGTAAGCGT", -142716.0, -384.1639999999999, 72.3790117774858, 70.7 },
        new object[] { "TTACTGTCLALALGGCAACT", "AATGACAGCTCCGTTGA", -116357.0, -316.32099999999997, 63.927874820717705, 61.4 },
        new object[] { "GCGTCAAGLCLGLACATCAT", "CGCAGTTCCCTGTAGTA", -128708.0, -352.0039999999999, 64.77605112645784, 57.4 },
        new object[] { "CGACTTGTLCLCLATACCTA", "GCTGAACACGTATGGAT", -106045.0, -292.738, 56.58142789274967, 53.6 },
        new object[] { "CCATGCGTLALGLACAAGTG", "GGTACGCAGCTGTTCAC", -113930.0, -306.07399999999996, 66.99414903034778, 65.6 },
        new object[] { "CTATCGCALTLCLTAATAAT", "GATAGCGTCGATTATTA", -106416.0, -302.4029999999999, 48.081383202006236, 50.6 },
    };

    [TestCaseSource(nameof(TripletCases))]
    public void Owczarzy2011_LnaTriplets_MatchMelting(string probe, string target, double h, double s, double tm, double measured)
    {
        AssertMelting(probe, target, Owc, 2e-6, 1, 0, h, s, tm);
    }

    [Test]
    public void Owczarzy2011_PerfectMatchTriplets_WithinMeasuredTm()
    {
        // The seven perfect-match triplet duplexes: model vs measured Tm (Owczarzy 2011), mean |error| 1.007 °C.
        double sum = 0;
        for (int i = 0; i < 7; i++)
        {
            var row = (object[])TripletCases[i];
            var (seq, lna) = Parse((string)row[0]);
            double t = PrimerDesigner.CalculateMeltingTemperatureNNLna(seq, lna, Owc, (string)row[1], 2e-6, 1.0,
                saltMode: PrimerDesigner.SaltCorrectionMode.Owczarzy2004Monovalent, gasConstant: MeltingR);
            sum += Math.Abs(t - (double)row[5]);
        }
        Assert.That(sum / 7, Is.EqualTo(1.0068).Within(0.001));
    }

    [Test]
    public void Owczarzy2011_CentralMismatch_LowersTmVsPerfectMatch()
    {
        for (int i = 0; i < 7; i++)
            Assert.That((double)((object[])TripletCases[i + 7])[4], Is.LessThan((double)((object[])TripletCases[i])[4]));
    }

    #endregion

    #region McTigue (2004) single-LNA duplexes — 5 µM, 1 M Na⁺, both models (MELTING)

    [TestCase("CCATTLGCTACC", -80314.0, -217.493, 55.276455530336534, -217.4, 55.40140463653762)]
    [TestCase("GGACLCTCGAC", -70614.0, -186.42799999999997, 57.63065896930431, -186.4, 57.674050443054966)]
    [TestCase("ACGTLCTTCG", -66442.0, -179.15699999999995, 49.06206516637917, -179.09999999999997, 49.15115657410496)]
    [TestCase("GTAGCGATLGTA", -80624.0, -220.181, 52.95954111103447, -220.1, 53.06641952486666)]
    [TestCase("CACLGGCTC", -57425.0, -151.11499999999998, 49.1657379120395, -151.09999999999997, 49.19287666227518)]
    public void McTigueSet_BothModels_MatchMelting(string probe, double h, double sOwc, double tmOwc, double sMct, double tmMct)
    {
        AssertMelting(probe, null, Owc, 5e-6, 1, 0, h, sOwc, tmOwc);
        AssertMelting(probe, null, Mct, 5e-6, 1, 0, h, sMct, tmMct);
    }

    [Test]
    public void McTigue2004_ConsecutiveRun_UsesOwczarzyTables_AsMelting()
    {
        // GCATT(L)G(L)CTA(L)CCAG: run 4–5 → Owczarzy single/tandem; isolated LNA 8 → McTigue increments.
        AssertMelting("GCATTLGLCTALCCAG", null, Mct, 5e-6, 1, 0, -96058.0, -253.39100000000002, 69.37626216810469);
        AssertMelting("GCATTLGLCTALCCAG", null, Owc, 5e-6, 1, 0, -96058.0, -253.44600000000003, 69.30909891329725);
    }

    #endregion

    #region Salt, Mg²⁺ and DNA mismatch (MELTING)

    [Test]
    public void Owczarzy2004Sodium_MatchesMelting()
    {
        AssertMelting("CCATTLGCTACC", null, Owc, 5e-6, 0.05, 0, -80314.0, -217.493, 41.57149460316782);
    }

    [Test]
    public void Owczarzy2008Magnesium_MatchesMelting()
    {
        AssertMelting("CCATTLGCTACC", null, Owc, 5e-6, 0.05, 0.003, -80314.0, -217.493, 49.48637755780163);
    }

    [Test]
    public void InternalDnaMismatchAwayFromLna_UsesAllawiSantaLuciaPeyret_MatchesMelting()
    {
        // T·T mismatch at position 7 (target 3'→5' GGTAACGTTGG).
        AssertMelting("CCATTLGCTACC", "GGTAACGTTGG", Owc, 5e-6, 1, 0, -70114.0, -192.493, 46.215135462730984);
    }

    #endregion

    #region Reduction, guards and not-computable cases

    [TestCase("CCATTGCTACC", 59.833634529845824)]
    [TestCase("GTGCATCGATGCAGC", 75.52391452117843)]
    [TestCase("GCATATGC", 46.231595611732246)] // self-complementary: symmetry term, k = C_T
    public void NoLna_ReducesToBiopythonTmNN_DnaNn3(string seq, double biopython)
    {
        // Biopython Tm_NN(seq, nn_table=DNA_NN3, dnac1=dnac2=C_T/2 (self-comp: dnac1=C_T, dnac2=0), Na=1000, saltcorr=0)
        foreach (var model in new[] { Owc, Mct })
        {
            double t = PrimerDesigner.CalculateMeltingTemperatureNNLna(seq, Array.Empty<int>(), model, null,
                1e-4, 1.0, saltMode: PrimerDesigner.SaltCorrectionMode.None, gasConstant: 1.987);
            Assert.That(t, Is.EqualTo(biopython).Within(1e-9), $"{seq} {model}");
        }
    }

    [Test]
    public void NoLna_SelfComplementary_At50mM_ReducesToBiopythonTmNN_Method6()
    {
        // Biopython Tm_NN("GCATATGC", nn_table=DNA_NN3, dnac1=500, dnac2=0, selfcomp=True, Na=50, saltcorr=6)
        double t = PrimerDesigner.CalculateMeltingTemperatureNNLna("GCATATGC", Array.Empty<int>(), Owc, null,
            0.5e-6, 0.05, gasConstant: 1.987);
        Assert.That(t, Is.EqualTo(16.621392992113726).Within(1e-9));
    }

    [Test]
    public void LnaModifiedSelfComplementarySequence_IsNotSymmetricDuplex()
    {
        // An LNA strand paired with an unmodified DNA complement is a heteroduplex (MELTING rejects -self with LNA).
        var th = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("GCATATGC", new[] { 3 });
        Assert.That(th!.Value.IsSelfComplementary, Is.False);
        var dna = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("GCATATGC", Array.Empty<int>());
        Assert.That(dna!.Value.IsSelfComplementary, Is.True);
    }

    [Test]
    public void NotComputable_Cases()
    {
        const string seq = "CCATTGCTACC";
        Assert.Multiple(() =>
        {
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 0 }), Is.Null, "terminal LNA");
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 10 }), Is.Null, "terminal LNA");
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 11 }), Is.Null, "out of range");
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { -1 }), Is.Null, "negative");
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(null!, Array.Empty<int>()), Is.Null);
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("", Array.Empty<int>()), Is.Null);
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("A", Array.Empty<int>()), Is.Null);
            Assert.That(PrimerDesigner.CalculateMeltingTemperatureNNLna("CCANTGCTACC", new[] { 4 }), Is.NaN, "non-ACGT");
            // Target of another length / terminal mismatch.
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 4 }, Owc, "GGTAACGATG"), Is.Null);
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 4 }, Owc, "AGTAACGATGG"), Is.Null);
            // Mismatch opposite an isolated LNA: no published parameter (MELTING: "GAL/CC … missing").
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("CCAGACAGG", new[] { 4 }, Owc, "GGTCCGTCC"), Is.Null);
            // Owczarzy mismatch table has no GLAL/CA (MELTING: "GLAL/C A … missing").
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("CCATGACAGG", new[] { 4, 5, 6 }, Owc, "GGTACAGTCC"), Is.Null);
            Assert.Throws<ArgumentNullException>(() => PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, null!));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 4 }, (PrimerDesigner.LnaNearestNeighborModel)7));
        });
    }

    [Test]
    public void DuplicateAndUnorderedPositions_SetSemantics()
    {
        const string seq = "CCATTGCTACC";
        var a = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 4, 4 });
        var b = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 4 });
        var c = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 6, 4 });
        var d = PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { 4, 6 });
        Assert.Multiple(() =>
        {
            Assert.That(a, Is.EqualTo(b));
            Assert.That(c, Is.EqualTo(d));
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna("ccattgctacc", new[] { 4 }), Is.EqualTo(b), "case-insensitive");
        });
    }

    [Test]
    public void EveryInternalLnaContext_IsParameterised_BothModels()
    {
        const string seq = "AACAGATCCGCTGGTTA"; // all 16 dinucleotides
        for (int i = 1; i < seq.Length - 1; i++)
        {
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { i }, Owc), Is.Not.Null, $"owc {i}");
            Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { i }, Mct), Is.Not.Null, $"mct {i}");
            if (i + 1 < seq.Length - 1)
                Assert.That(PrimerDesigner.CalculateNearestNeighborThermodynamicsLna(seq, new[] { i, i + 1 }, Owc), Is.Not.Null, $"tandem {i}");
        }
    }

    [Test]
    public void ConditionGuards_Throw()
    {
        var p = new[] { 4 };
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateMeltingTemperatureNNLna("CCATTGCTACC", p, 0.0));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateMeltingTemperatureNNLna("CCATTGCTACC", p, 1e-6, 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateMeltingTemperatureNNLna("CCATTGCTACC", p, 1e-6, 0.05, -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => PrimerDesigner.CalculateMeltingTemperatureNNLna("CCATTGCTACC", p, Owc, null, gasConstant: 0));
        });
    }

    #endregion

    #region MGB design rules (Kutyavin 2000)

    [Test]
    public void EvaluateMgbProbeDesign_LengthWindowAndThreePrimePlacement()
    {
        var fifteen = ProbeDesigner.EvaluateMgbProbeDesign("ACGTACGTACGTACG");
        var twentyFive = ProbeDesigner.EvaluateMgbProbeDesign("ACGTACGTACGTACGTACGTACGTA");
        var twelve = ProbeDesigner.EvaluateMgbProbeDesign("ACGTACGTACGT");
        var eleven = ProbeDesigner.EvaluateMgbProbeDesign("ACGTACGTACG");
        Assert.Multiple(() =>
        {
            Assert.That(fifteen.Length, Is.EqualTo(15));
            Assert.That(fifteen.LengthInMgbRange, Is.True);
            Assert.That(fifteen.MgbAttachmentEnd, Is.EqualTo("3'"));
            Assert.That(twelve.LengthInMgbRange, Is.True, "12mer is the lower bound");
            Assert.That(eleven.LengthInMgbRange, Is.False);
            Assert.That(twentyFive.LengthInMgbRange, Is.False);
            Assert.That(twentyFive.Guidance, Has.Some.Contains("12-20mer"));
        });
    }

    [Test]
    public void EvaluateMgbProbeDesign_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => ProbeDesigner.EvaluateMgbProbeDesign(null!));
    }

    #endregion
}
