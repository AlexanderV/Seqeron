using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// The canonical nearest-neighbour core <see cref="ThermoConstants.CalculateNearestNeighborDuplex"/> and the
/// basic-Tm helper <see cref="ThermoConstants.CalculateBasicTm"/> (B07 PRIMER-NNTM-001, cross-batch R9/R10).
/// Expected values are Biopython 1.88 <c>Bio.SeqUtils.MeltingTemp</c> outputs (<c>Tm_NN</c>,
/// <c>salt_correction</c>), printed with <c>repr</c>; the port is bit-exact on 23 000 random inputs (12 635 valid duplexes + 10 365 error cases)
/// (all four DNA tables, salt methods 0–7, K⁺/Tris/Mg²⁺/dNTPs, mismatches, dangling ends, shifts, inosine).
/// </summary>
[TestFixture]
public class ThermoConstants_NearestNeighborCore_Tests
{
    private const string Seq28 = "CGTTCCAAAGATGTGGGCATGAGCTTAC";
    private const string Seq16 = "CGTTCCAAAGATGTGG";
    private const double Exact = 1e-9;

    private static double Tm(string seq, string? cSeq = null, int shift = 0,
        NnParameterSet set = NnParameterSet.AllawiSantaLucia1997, double dnac1 = 25, double dnac2 = 25,
        bool self = false, double na = 50, double k = 0, double tris = 0, double mg = 0, double dntps = 0,
        NnSaltCorrection salt = NnSaltCorrection.SantaLucia1998Entropy, bool strict = true) =>
        ThermoConstants.CalculateNearestNeighborTm(seq, cSeq, shift, set, dnac1, dnac2, self,
            na, k, tris, mg, dntps, salt, strict: strict);

    // Biopython: Tm_NN(seq, nn_table=DNA_NN1..4) at the defaults (Na 50 mM, dnac1 = dnac2 = 25 nM, saltcorr 5).
    [TestCase(NnParameterSet.AllawiSantaLucia1997, Seq28, 60.32091949919129)]
    [TestCase(NnParameterSet.Breslauer1986, Seq28, 72.19499984118579)]
    [TestCase(NnParameterSet.Breslauer1986, "ATTATAATAAATTA", 20.52696567737297)] // init_allA/T branch
    [TestCase(NnParameterSet.Sugimoto1996, Seq28, 65.46933554780924)]
    [TestCase(NnParameterSet.SantaLuciaHicks2004, Seq28, 60.27151268472841)]
    public void Tm_EachParameterSet_MatchesBiopythonTmNN(NnParameterSet set, string seq, double expected)
    {
        Assert.That(Tm(seq, set: set), Is.EqualTo(expected).Within(Exact));
    }

    // Biopython: Tm_NN(Seq16, saltcorr=m, Na, K, Tris, Mg, dNTPs) — every salt_correction method,
    // incl. the von Ahsen sodium equivalent (methods 1–6) and the three Owczarzy 2008 regimes.
    [TestCase(NnSaltCorrection.SchildkrautLifson1965, 50, 10, 20, 1.5, 0.2, 46.45314133724665)]
    [TestCase(NnSaltCorrection.Wetmur1991, 100, 0, 20, 0, 0, 41.366606539326)]
    [TestCase(NnSaltCorrection.SantaLucia1996, 200, 0, 0, 0, 0, 49.077143586245725)]
    [TestCase(NnSaltCorrection.SantaLucia1998Tm, 20, 0, 0, 3, 0.8, 49.584963002528525)]
    [TestCase(NnSaltCorrection.SantaLucia1998Entropy, 50, 0, 0, 1.5, 0.6, 49.20942215244429)]
    [TestCase(NnSaltCorrection.Owczarzy2004, 50, 50, 0, 0, 0, 48.09656963636297)]
    [TestCase(NnSaltCorrection.Owczarzy2008, 500, 0, 0, 0.5, 0, 55.959566551580394)] // R < 0.22
    [TestCase(NnSaltCorrection.Owczarzy2008, 50, 0, 0, 1.5, 0.6, 48.27621274096333)] // 0.22 ≤ R < 6
    [TestCase(NnSaltCorrection.Owczarzy2008, 0, 0, 0, 10, 0, 53.293367200046305)]   // no monovalent
    public void Tm_EachSaltMethod_MatchesBiopythonTmNN(
        NnSaltCorrection salt, double na, double k, double tris, double mg, double dntps, double expected)
    {
        Assert.That(Tm(Seq16, na: na, k: k, tris: tris, mg: mg, dntps: dntps, salt: salt),
            Is.EqualTo(expected).Within(Exact));
    }

    [Test]
    public void SaltCorrection_MatchesBiopythonSaltCorrection()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ThermoConstants.CalculateNnSaltCorrection(NnSaltCorrection.SchildkrautLifson1965, sodium: 50),
                Is.EqualTo(-21.59709792802209).Within(1e-12));
            Assert.That(ThermoConstants.CalculateNnSaltCorrection(NnSaltCorrection.Wetmur1991, sodium: 100, tris: 20, magnesium: 1.5),
                Is.EqualTo(-10.988468412525648).Within(1e-12));
            Assert.That(ThermoConstants.CalculateNnSaltCorrection(NnSaltCorrection.SantaLucia1998Entropy, sodium: 50, sequence: "ACGTACGT"),
                Is.EqualTo(-7.71700633667508).Within(1e-12));
            Assert.That(ThermoConstants.CalculateNnSaltCorrection(NnSaltCorrection.Owczarzy2004, sodium: 50, sequence: "ACGTACGGG"),
                Is.EqualTo(0.00011701295321698038).Within(1e-18));
            Assert.That(ThermoConstants.CalculateNnSaltCorrection(NnSaltCorrection.Owczarzy2008, sodium: 50, magnesium: 1.5, dntps: 0.6, sequence: "ACGTACGGG"),
                Is.EqualTo(9.1205730064332e-05).Within(1e-18));
            Assert.That(ThermoConstants.CalculateNnSaltCorrection(NnSaltCorrection.Owczarzy2008, magnesium: 1.5, sequence: "ACGTACGGG"),
                Is.EqualTo(5.472338115589228e-05).Within(1e-18));
            Assert.That(ThermoConstants.CalculateNnSaltCorrection(NnSaltCorrection.None), Is.EqualTo(0));
        });
    }

    // Biopython Tm_NN with c_seq / shift: internal mismatch (DNA_IMM1), terminal mismatches (DNA_TMM1),
    // dangling ends (DNA_DE1) on either side, over-dangling trimming, negative shift, inosine.
    [TestCase("internal mismatch", Seq16, "GCAAGGTCTCTACACC", 0, NnParameterSet.AllawiSantaLucia1997, 32.081928634970666)]
    [TestCase("terminal mismatches", Seq16, "ACAAGGTTTCTACACA", 0, NnParameterSet.SantaLuciaHicks2004, 38.13300510881902)]
    [TestCase("left dangling (shift 1)", Seq16, "TGCAAGGTTTCTACACC", 1, NnParameterSet.AllawiSantaLucia1997, 44.41652011600104)]
    [TestCase("right dangling (padding)", "CGTTCCAAAGATGTGGA", "GCAAGGTTTCTACACC", 0, NnParameterSet.AllawiSantaLucia1997, 45.09243484077098)]
    [TestCase("over-dangling trimmed", Seq16, "AAAGCAAGGTTTCTACACCTT", 3, NnParameterSet.AllawiSantaLucia1997, 47.27826498960451)]
    [TestCase("negative shift", Seq16, "CAAGGTTTCTACACC", -1, NnParameterSet.AllawiSantaLucia1997, 39.13770984597102)]
    [TestCase("inosine", "CGTTCCAAIGATGTGG", "GCAAGGTTCCTACACC", 0, NnParameterSet.AllawiSantaLucia1997, 42.70983249202294)]
    public void Tm_MismatchAndDanglingDuplexes_MatchBiopythonTmNN(
        string label, string seq, string cSeq, int shift, NnParameterSet set, double expected)
    {
        Assert.That(Tm(seq, cSeq, shift, set), Is.EqualTo(expected).Within(Exact), label);
    }

    [Test]
    public void Tm_SelfComplementaryAndConcentrations_MatchBiopythonTmNN()
    {
        Assert.Multiple(() =>
        {
            // Tm_NN("GCGCATGCGC", selfcomp=True, dnac1=500): k = dnac1, symmetry term.
            Assert.That(Tm("GCGCATGCGC", self: true, dnac1: 500), Is.EqualTo(47.09214772582419).Within(Exact));
            // Tm_NN(Seq16, dnac1=250, dnac2=0): PCR template limit k = dnac1.
            Assert.That(Tm(Seq16, dnac1: 250, dnac2: 0), Is.EqualTo(48.707368868496474).Within(Exact));
        });
    }

    [Test]
    public void Tm_Check_NormalisesLikeBiopython()
    {
        // _check: whitespace removed, upper-cased, U back-transcribed → CGTTCCAAAGATGTGG.
        Assert.Multiple(() =>
        {
            Assert.That(Tm("cgu ucc aaa gau gug g"), Is.EqualTo(43.80238368804726).Within(Exact));
            Assert.That(Tm(Seq16), Is.EqualTo(43.80238368804726).Within(Exact));
            Assert.That(ThermoConstants.NormalizeForNearestNeighbor("a-cGu N i."), Is.EqualTo("ACGTI"));
        });
    }

    [Test]
    public void Tm_StrictAndNonStrict_MatchBiopython()
    {
        // 'AG/AA' (two adjacent mismatches) has no parameter: strict → error; non-strict → skipped.
        Assert.Multiple(() =>
        {
            Assert.That(() => Tm(Seq16, "GCAAGGTTAATACACC"), NUnit.Framework.Throws.ArgumentException);
            Assert.That(Tm(Seq16, "GCAAGGTTAATACACC", strict: false), Is.EqualTo(30.44830734549629).Within(Exact));
            Assert.That(Tm(Seq16, "GCAAGGTTTCTACAGG", strict: false), Is.EqualTo(33.76507074792329).Within(Exact));
        });
    }

    [Test]
    public void Tm_InvalidInputs_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => Tm("CGTT", na: 0), NUnit.Framework.Throws.ArgumentException, "zero total ion concentration (Biopython ValueError)");
            Assert.That(() => Tm("NNNN"), NUnit.Framework.Throws.ArgumentException, "nothing left after _check (Biopython IndexError)");
            Assert.That(() => Tm(Seq16, dnac1: 10, dnac2: 20), NUnit.Framework.Throws.ArgumentException, "k ≤ 0 (Biopython math domain error)");
            Assert.That(() => ThermoConstants.CalculateNnSaltCorrection(NnSaltCorrection.Owczarzy2004, sodium: 50),
                NUnit.Framework.Throws.ArgumentException, "sequence missing for method 6");
        });
    }

    [Test]
    public void Thermodynamics_AreTheSumsOfTheTable1Terms()
    {
        // NN4 (SantaLucia & Hicks 2004 Table 1), ACGT: init (0.2, −5.7) + 2 terminal A·T (2.2, 6.9)
        // + AC (−8.4, −22.4) + CG (−10.6, −27.2) + GT (−8.4, −22.4).
        var (dh, ds) = ThermoConstants.CalculateNearestNeighborThermodynamics(
            "ACGT", parameterSet: NnParameterSet.SantaLuciaHicks2004);
        var init = ThermoConstants.GetNearestNeighborInitiation(NnParameterSet.SantaLuciaHicks2004);
        Assert.Multiple(() =>
        {
            Assert.That(dh, Is.EqualTo(0.2 + 2 * 2.2 - 8.4 - 10.6 - 8.4).Within(1e-12));
            Assert.That(ds, Is.EqualTo(-5.7 + 2 * 6.9 - 22.4 - 27.2 - 22.4).Within(1e-12));
            Assert.That(init.Initiation, Is.EqualTo((0.2, -5.7)));
            Assert.That(init.TerminalAT, Is.EqualTo((2.2, 6.9)));
            Assert.That(init.Symmetry, Is.EqualTo((0.0, -1.4)));
            Assert.That(ThermoConstants.TryGetNearestNeighborStack(NnParameterSet.SantaLuciaHicks2004, "AA", out var aa), Is.True);
            Assert.That(aa, Is.EqualTo((-7.6, -21.3)), "2004 AA/TT (1998: −7.9/−22.2)");
            Assert.That(ThermoConstants.TryGetNearestNeighborStack(NnParameterSet.AllawiSantaLucia1997, "TT", out var tt), Is.True);
            Assert.That(tt, Is.EqualTo((-7.9, -22.2)));
            Assert.That(ThermoConstants.TryGetNearestNeighborStack(NnParameterSet.SantaLuciaHicks2004, "TA", out var ta), Is.True);
            Assert.That(ta, Is.EqualTo((-7.2, -21.3)), "TA/AT ΔS° −21.3 (ΔG°37 −0.58; MELTING's −20.4 is a typo)");
            Assert.That(ThermoConstants.TryGetNearestNeighborStack(NnParameterSet.SantaLuciaHicks2004, "AN", out _), Is.False);
        });
    }

    /// <summary>
    /// R9 routing contract for B03: the core with DNA_NN3 (per-terminal-pair initiation), salt method 5,
    /// dnac1 = dnac2 = C_T/2 (self-complementary: dnac1 = C_T) and ΔG°37 = ΔH° − 310.15·ΔS°/1000 reproduces
    /// <see cref="SequenceStatistics.CalculateThermodynamics(string, double, double, bool)"/> field by field.
    /// </summary>
    [Test]
    public void Core_NN3_SaltMethod5_ReproducesSequenceStatisticsCalculateThermodynamics()
    {
        var rng = new Random(9);
        const string alphabet = "ACGTACGTacgtUN";
        double[] sodium = { 0.01, 0.05, 0.1, 1.0 };
        double[] strands = { 2.5e-7, 1e-6, 5e-8 };
        int compared = 0;
        for (int i = 0; i < 3000; i++)
        {
            int n = rng.Next(2, 50);
            var chars = new char[n];
            for (int j = 0; j < n; j++) chars[j] = alphabet[rng.Next(alphabet.Length)];
            string s = new(chars);
            double na = sodium[rng.Next(sodium.Length)], ct = strands[rng.Next(strands.Length)];
            bool self = rng.Next(4) == 0;

            var expected = SequenceStatistics.CalculateThermodynamics(s, na, ct, self);
            if (ThermoConstants.NormalizeForNearestNeighbor(s).Replace("I", "").Length < 2)
                continue;
            var r = ThermoConstants.CalculateNearestNeighborDuplex(s,
                parameterSet: NnParameterSet.AllawiSantaLucia1997,
                dnac1: self ? ct * 1e9 : ct * 1e9 / 2, dnac2: self ? 0 : ct * 1e9 / 2,
                selfComplementary: self, sodium: na * 1000,
                saltCorrection: NnSaltCorrection.SantaLucia1998Entropy);
            double dg = r.DeltaH - 310.15 * r.SaltCorrectedDeltaS / 1000.0;
            Assert.Multiple(() =>
            {
                Assert.That(Math.Round(r.DeltaH, 2), Is.EqualTo(expected.DeltaH), s);
                Assert.That(Math.Round(r.SaltCorrectedDeltaS, 2), Is.EqualTo(expected.DeltaS), s);
                Assert.That(Math.Round(dg, 2), Is.EqualTo(expected.DeltaG), s);
                Assert.That(Math.Round(r.MeltingTemperature, 1), Is.EqualTo(expected.MeltingTemperature), s);
            });
            compared++;
        }
        Assert.That(compared, Is.GreaterThan(2500));
    }

    [TestCase(null, 0.0)]
    [TestCase("", 0.0)]
    [TestCase("NNNN", 0.0)]
    [TestCase("ACGT", 12.0)]            // Wallace 2·2 + 4·2
    [TestCase("ACGU", 12.0)]            // U read as T (Biopython _check)
    [TestCase("acgtacgtacgtac", 37.371428571428574)] // 14 counted bases: 64.9 + 41·(7 − 16.4)/14
    public void CalculateBasicTm_OligoCalcBasicFormulas(string? seq, double expected)
    {
        Assert.That(ThermoConstants.CalculateBasicTm(seq), Is.EqualTo(expected).Within(1e-12));
    }
}
