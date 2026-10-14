// SEQ-HYDRO-001 — profile options (FIN-B03 WP6/WP7, F29–F31): Biopython window convention,
// strict unknown-residue mode, ExPASy ProtScale 0..1 scale normalisation, caller-supplied scales.
// Oracle: Biopython 1.88 Bio.SeqUtils.ProtParam.ProteinAnalysis(seq).protein_scale(param_dict, window, edge)
// executed 2026-10-14 (script + output: docs/Evidence/SEQ-HYDRO-001-Evidence.md, "FIN-B03 WP6/WP7").
// Normalised cases: protein_scale({k: (v - min)/(max - min)}, window, edge) over the 20-entry ProtParamData table.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class SequenceStatistics_HydrophobicityProfileOptions_Tests
{
    private const double Tolerance = 1e-12;
    private const string Seq = "MKTAYIAKQRQISFVKSHFSRQ";

    private static readonly SequenceStatistics.HydrophobicityProfileOptions Bio =
        new() { Convention = SequenceStatistics.HydropathyWindowConvention.Biopython };

    // Shipped Biopython 1.88 ProtParamData scales (F52); the literal oracle copy lives in
    // ProteinScales_EqualBiopythonProtParamData below.
    private static readonly IReadOnlyDictionary<char, double> HoppWoods = ProteinScales.HoppWoods;
    private static readonly IReadOnlyDictionary<char, double> Emini = ProteinScales.Emini;
    private static readonly IReadOnlyDictionary<char, double> Janin = ProteinScales.Janin;
    private static readonly IReadOnlyDictionary<char, double> Flexibility = ProteinScales.Flexibility;

    private static void AssertProfile(IEnumerable<double> actual, double[] expected)
    {
        var profile = actual.ToList();
        Assert.That(profile, Has.Count.EqualTo(expected.Length));
        Assert.Multiple(() =>
        {
            for (int i = 0; i < expected.Length; i++)
                Assert.That(profile[i], Is.EqualTo(expected[i]).Within(Tolerance), $"window {i}");
        });
    }

    #region Biopython window convention (A11, F29)

    // protein_scale(kd, 4, 1.0): pairs (0,3),(1,2) + middle residue 2 again; / 5
    [Test]
    public void Biopython_EvenWindow4Edge1_CountsCentreTwiceOverWPlus1()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, 4, 1.0, Bio),
            new[] { -0.3199999999999999, -0.45999999999999996, 0.6, 2.2600000000000002, 0.5799999999999998, -1.0, -2.72, -3.9799999999999995, -2.1, 0.040000000000000036, 0.44000000000000006, 2.7, 1.3, -0.32, -0.9, -1.6600000000000001, 0.15999999999999986, -1.3, -2.1 });
    }

    // protein_scale(kd, 2, 1.0): (kd[i]+2·kd[i+1])/3
    [Test]
    public void Biopython_EvenWindow2Edge1_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, 2, 1.0, Bio),
            new[] { -1.9666666666666668, -1.7666666666666666, 0.9666666666666668, -0.26666666666666666, 2.566666666666667, 2.6999999999999997, -2.0, -3.6333333333333333, -4.166666666666667, -3.8333333333333335, 1.8333333333333333, 0.9666666666666668, 1.5999999999999999, 3.733333333333333, -1.2, -1.8333333333333333, -2.4, 0.7999999999999998, 0.3999999999999999, -3.266666666666667, -3.8333333333333335 });
    }

    // protein_scale(kd, 6, 0.5): even weighted window accepted by Biopython
    [Test]
    public void Biopython_EvenWindow6Edge05_AcceptedAndMatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, 6, 0.5, Bio),
            new[] { 0.451923076923077, 0.14615384615384625, 1.4615384615384615, 0.5730769230769232, -1.5365384615384616, -2.2211538461538463, -2.640384615384615, -2.375, 0.25769230769230755, 0.36346153846153834, 1.3442307692307691, 1.7423076923076928, -0.8134615384615383, -0.29423076923076913, -1.128846153846154, -0.553846153846154, -1.2576923076923074 });
    }

    // protein_scale(kd, 5, 1.0) on ACDXEFGHIK: X drops its pair partner; divisor unchanged
    [Test]
    public void Biopython_OffCentreUnknown_DropsSymmetricPartner()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile("ACDXEFGHIK", 5, 1.0, Bio),
            new[] { -1.04, -0.34, -1.48, -0.22000000000000003, 0.039999999999999925, -0.040000000000000126 });
    }

    // protein_scale(kd, 5, 0.4) on ACDXEFGHIK
    [Test]
    public void Biopython_OffCentreUnknownEdge04_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile("ACDXEFGHIK", 5, 0.4, Bio),
            new[] { -1.30625, -0.8687499999999998, -1.58125, 0.02187500000000009, -0.08750000000000004, -0.24062500000000014 });
    }

    // ProteinAnalysis upper-cases: same as ACDXEFGHIK
    [Test]
    public void Biopython_LowercaseInput_IsUpperCased()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile("acdxefghik", 5, 0.4, Bio),
            new[] { -1.30625, -0.8687499999999998, -1.58125, 0.02187500000000009, -0.08750000000000004, -0.24062500000000014 });
    }

    // AXCDE: pair (X,D) dropped → (A+E+C)/5 = 0.16
    [Test]
    public void Biopython_OuterPairUnknown_DropsAlanine()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile("AXCDE", 5, 1.0, Bio),
            new[] { 0.16 });
    }

    // AC-DEB, W=4: '-' and 'B' drop their pairs
    [Test]
    public void Biopython_GapAndAmbiguityCodes_MatchProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile("AC-DEB", 4, 1.0, Bio),
            new[] { -0.33999999999999997, -0.9, -2.1 });
    }

    // Standard convention keeps the off-centre unknown's partner: AXCDE → (1.8+2.5−3.5−3.5)/5 = −0.54
    // (≠ Biopython 0.16), and ACXDE (centre X) agrees with Biopython −0.54 in both conventions.
    [Test]
    public void Standard_VersusBiopython_DifferOnlyForOffCentreUnknown()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceStatistics.CalculateHydrophobicityProfile("AXCDE", 5, 1.0).Single(),
                Is.EqualTo(-0.54).Within(Tolerance));
            Assert.That(SequenceStatistics.CalculateHydrophobicityProfile("ACXDE", 5, 1.0, Bio).Single(),
                Is.EqualTo(-0.54).Within(Tolerance));
            Assert.That(SequenceStatistics.CalculateHydrophobicityProfile("ACXDE", 5, 1.0).Single(),
                Is.EqualTo(-0.54).Within(Tolerance));
        });
    }

    // Odd window, no unknowns: both conventions coincide (Biopython parity of the default path).
    [TestCase(9, 1.0)]
    [TestCase(5, 0.4)]
    public void Biopython_OddWindowStandardResidues_EqualsStandardConvention(int w, double edge)
    {
        var standard = SequenceStatistics.CalculateHydrophobicityProfile(Seq, w, edge).ToArray();
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, w, edge, Bio), standard);
    }

    // protein_scale(kd, 1, 1.0) raises ZeroDivisionError (unit = 2(1−edge)/(W−1)).
    [Test]
    public void Biopython_Window1_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(
            () => SequenceStatistics.CalculateHydrophobicityProfile("ACDE", 1, 1.0, Bio));
    }

    [Test]
    public void Standard_EvenWeightedWindow_StillThrows()
    {
        Assert.Throws<ArgumentException>(() => SequenceStatistics.CalculateHydrophobicityProfile(
            Seq, 6, 0.5, SequenceStatistics.HydrophobicityProfileOptions.Default));
    }

    [Test]
    public void DefaultOptions_EqualLegacyOverload()
    {
        var legacy = SequenceStatistics.CalculateHydrophobicityProfile(Seq, 5, 0.4).ToArray();
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, 5, 0.4, null), legacy);
    }

    #endregion

    #region Strict mode (A11, F29)

    // Biopython ProteinAnalysis("ACDX").gravy() raises KeyError 'X'; strict mode mirrors that.
    [TestCase("ACDX")]
    [TestCase("AC-DE")]
    [TestCase("ACDEB")]
    [TestCase("MKſ")]
    public void Strict_UnknownSymbol_Throws(string seq)
    {
        var strict = new SequenceStatistics.HydrophobicityProfileOptions { Strict = true };
        Assert.Throws<ArgumentException>(() => SequenceStatistics.CalculateHydrophobicityProfile(seq, 3, 1.0, strict));
        Assert.Throws<ArgumentException>(() => SequenceStatistics.CalculateHydrophobicityProfile(seq, 3, 1.0, strict with { Convention = SequenceStatistics.HydropathyWindowConvention.Biopython }));
    }

    [Test]
    public void Strict_StandardResiduesLowercase_EqualsLenient()
    {
        var strict = new SequenceStatistics.HydrophobicityProfileOptions { Strict = true };
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq.ToLowerInvariant(), 9, 1.0, strict),
            SequenceStatistics.CalculateHydrophobicityProfile(Seq, 9, 1.0).ToArray());
    }

    #endregion

    #region ProtScale 0..1 normalisation (A12, F30)

    // kd min −4.5 (R), max 4.5 (I): v' = (v + 4.5)/9
    [Test]
    public void Normalized_KyteDoolittleWindow9_MatchesProteinScaleOnMinMaxScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, 9, 1.0, new SequenceStatistics.HydrophobicityProfileOptions { NormalizeScale = true }),
            new[] { 0.4592592592592592, 0.3802469135802469, 0.3851851851851852, 0.4493827160493827, 0.41728395061728396, 0.4679012345679012, 0.46419753086419746, 0.39382716049382716, 0.43209876543209874, 0.4358024691358025, 0.5259259259259259, 0.5592592592592592, 0.44814814814814824, 0.4148148148148148 });
    }

    // normalised kd, W=5, edge 0.4
    [Test]
    public void Normalized_KyteDoolittleEdge04_MatchesProteinScaleOnMinMaxScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, 5, 0.4, new SequenceStatistics.HydrophobicityProfileOptions { NormalizeScale = true }),
            new[] { 0.43298611111111107, 0.5222222222222221, 0.6232638888888888, 0.6392361111111111, 0.5104166666666666, 0.3232638888888889, 0.15069444444444444, 0.18194444444444444, 0.31875, 0.528125, 0.6593749999999999, 0.6881944444444443, 0.5968749999999999, 0.44166666666666665, 0.396875, 0.3722222222222222, 0.42638888888888893, 0.3378472222222222 });
    }

    // normalised kd, Biopython convention, unknown X
    [Test]
    public void Normalized_UnknownCentre_BiopythonConvention()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile("ACDXEFGHIK", 5, 1.0, new SequenceStatistics.HydrophobicityProfileOptions { NormalizeScale = true, Convention = SequenceStatistics.HydropathyWindowConvention.Biopython }),
            new[] { 0.18444444444444447, 0.3622222222222223, 0.13555555555555557, 0.27555555555555555, 0.5044444444444445, 0.49555555555555547 });
    }

    // hw min −3.4 (W), max 3.0: v' = (v + 3.4)/6.4
    [Test]
    public void Normalized_HoppWoodsEdge04_MatchesProteinScaleOnMinMaxScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, 5, 0.4, new SequenceStatistics.HydrophobicityProfileOptions { NormalizeScale = true, Scale = HoppWoods }),
            new[] { 0.5268554687499999, 0.43798828125, 0.32275390624999994, 0.396484375, 0.5068359375, 0.69091796875, 0.7402343749999999, 0.71484375, 0.5917968749999999, 0.47021484374999994, 0.37353515625, 0.39160156249999994, 0.48681640624999994, 0.5781249999999999, 0.5532226562499999, 0.49609374999999994, 0.46679687499999994, 0.5571289062499999 });
    }

    [Test]
    public void Normalized_ConstantOrEmptyScale_Throws()
    {
        Assert.Throws<ArgumentException>(() => SequenceStatistics.CalculateHydrophobicityProfile(Seq, 3, 1.0,
            new SequenceStatistics.HydrophobicityProfileOptions { NormalizeScale = true, Scale = new Dictionary<char, double> { ['A'] = 1, ['C'] = 1 } }));
        Assert.Throws<ArgumentException>(() => SequenceStatistics.CalculateHydrophobicityProfile(Seq, 3, 1.0,
            new SequenceStatistics.HydrophobicityProfileOptions { NormalizeScale = true, Scale = new Dictionary<char, double>() }));
    }

    #endregion

    #region Caller-supplied scales (A13, F31)

    // protein_scale(ProtParamData.hw, 9, 1.0)
    [Test]
    public void Scale_HoppWoodsWindow9_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, HoppWoods, 9),
            new[] { -0.06666666666666662, 0.4111111111111112, 0.10000000000000002, -0.055555555555555455, 0.033333333333333395, 0.01111111111111112, 0.04444444444444444, 0.43333333333333335, 0.1333333333333334, 0.05555555555555555, -0.5555555555555556, -0.5444444444444445, -0.011111111111111151, -0.02222222222222224 });
    }

    // protein_scale(ProtParamData.hw, 5, 0.4)
    [Test]
    public void Scale_HoppWoodsWindow5Edge04_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, HoppWoods, 5, 0.4),
            new[] { -0.02812500000000013, -0.5968749999999998, -1.3343749999999999, -0.8624999999999999, -0.15625000000000006, 1.0218749999999999, 1.3374999999999997, 1.1749999999999998, 0.38749999999999984, -0.390625, -1.009375, -0.8937499999999999, -0.2843750000000001, 0.3, 0.14062499999999983, -0.22499999999999992, -0.41249999999999987, 0.16562499999999986 });
    }

    // protein_scale(ProtParamData.em (Emini surface probability), 9, 1.0)
    [Test]
    public void Scale_EminiWindow9_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, Emini, 9),
            new[] { 1.0731111111111111, 1.1576666666666666, 1.1357777777777778, 1.0712222222222223, 1.1045555555555557, 1.0607777777777776, 1.0611111111111111, 1.1422222222222222, 1.0944444444444446, 1.0757777777777777, 0.989111111111111, 0.9632222222222223, 1.0601111111111112, 1.086 });
    }

    // protein_scale(ProtParamData.em (Emini surface probability), 5, 0.4)
    [Test]
    public void Scale_EminiWindow5Edge04_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, Emini, 5, 0.4),
            new[] { 1.1116249999999999, 1.02040625, 0.9003749999999999, 0.8999374999999998, 1.0291875, 1.21571875, 1.3522499999999997, 1.3191875, 1.1836874999999998, 0.99846875, 0.876625, 0.8621562499999998, 0.9581249999999999, 1.0936562499999998, 1.1071562499999998, 1.0971875, 1.04296875, 1.139125 });
    }

    // protein_scale(ProtParamData.ja, 9, 1.0)
    [Test]
    public void Scale_JaninWindow9_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, Janin, 9),
            new[] { -0.31222222222222223, -0.48666666666666664, -0.38333333333333336, -0.2811111111111111, -0.33333333333333326, -0.26555555555555554, -0.26555555555555554, -0.4766666666666667, -0.3177777777777777, -0.27555555555555555, -0.09777777777777778, -0.04222222222222224, -0.23555555555555557, -0.2911111111111111 });
    }

    // protein_scale(ProtParamData.ja, 5, 0.4)
    [Test]
    public void Scale_JaninWindow5Edge04_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, Janin, 5, 0.4),
            new[] { -0.35812499999999997, -0.1428125, 0.14062499999999997, 0.04843749999999995, -0.24062499999999992, -0.6634375, -0.870625, -0.785625, -0.44375, -0.08999999999999994, 0.16125, 0.10593749999999996, -0.11374999999999996, -0.3978125, -0.3490625, -0.26406250000000003, -0.131875, -0.33312499999999995 });
    }

    // protein_scale(ProtParamData.Flex, 9, 1.0)
    [Test]
    public void Scale_FlexibilityWindow9_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, Flexibility, 9),
            new[] { 1.0015555555555555, 1.0077777777777779, 1.0005555555555556, 0.9927777777777776, 0.9996666666666667, 0.9981111111111112, 0.9985555555555558, 1.0116666666666667, 1.0054444444444444, 0.9957777777777778, 0.9854444444444445, 0.9864444444444445, 0.9954444444444444, 0.9944444444444444 });
    }

    // protein_scale(ProtParamData.Flex, 5, 0.4)
    [Test]
    public void Scale_FlexibilityWindow5Edge04_MatchesProteinScale()
    {
        AssertProfile(SequenceStatistics.CalculateHydrophobicityProfile(Seq, Flexibility, 5, 0.4),
            new[] { 1.003, 0.9824374999999999, 0.9559687499999999, 0.9689062499999999, 0.99709375, 1.02834375, 1.03825, 1.0223125, 1.00771875, 0.98571875, 0.9758125000000001, 0.97203125, 0.9936562499999999, 1.00996875, 1.0065000000000002, 0.9943437500000002, 0.9793124999999999, 0.9959062499999999 });
    }

    // F52: the shipped tables equal Biopython 1.88 Bio.SeqUtils.ProtParamData hw / em / ja / Flex
    // (literal copy printed from the installed module 2026-10-14) and are read-only.
    [Test]
    public void ProteinScales_EqualBiopythonProtParamData()
    {
        var hw = new Dictionary<char, double>
        {
            ['A'] = -0.5, ['C'] = -1.0, ['D'] = 3.0, ['E'] = 3.0, ['F'] = -2.5, ['G'] = 0.0, ['H'] = -0.5, ['I'] = -1.8, ['K'] = 3.0, ['L'] = -1.8, ['M'] = -1.3, ['N'] = 0.2, ['P'] = 0.0, ['Q'] = 0.2, ['R'] = 3.0, ['S'] = 0.3, ['T'] = -0.4, ['V'] = -1.5, ['W'] = -3.4, ['Y'] = -2.3,
        };

        var em = new Dictionary<char, double>
        {
            ['A'] = 0.815, ['C'] = 0.394, ['D'] = 1.283, ['E'] = 1.445, ['F'] = 0.695, ['G'] = 0.714, ['H'] = 1.18, ['I'] = 0.603, ['K'] = 1.545, ['L'] = 0.603, ['M'] = 0.714, ['N'] = 1.296, ['P'] = 1.236, ['Q'] = 1.348, ['R'] = 1.475, ['S'] = 1.115, ['T'] = 1.184, ['V'] = 0.606, ['W'] = 0.808, ['Y'] = 1.089,
        };

        var ja = new Dictionary<char, double>
        {
            ['A'] = 0.28, ['C'] = 0.97, ['D'] = -0.52, ['E'] = -1.01, ['F'] = 0.46, ['G'] = 0.43, ['H'] = -0.31, ['I'] = 0.6, ['K'] = -1.62, ['L'] = 0.6, ['M'] = 0.43, ['N'] = -0.55, ['P'] = -0.42, ['Q'] = -0.69, ['R'] = -1.14, ['S'] = -0.19, ['T'] = -0.32, ['V'] = 0.6, ['W'] = 0.29, ['Y'] = -0.15,
        };

        var flex = new Dictionary<char, double>
        {
            ['A'] = 0.984, ['C'] = 0.906, ['D'] = 1.068, ['E'] = 1.094, ['F'] = 0.915, ['G'] = 1.031, ['H'] = 0.95, ['I'] = 0.927, ['K'] = 1.102, ['L'] = 0.935, ['M'] = 0.952, ['N'] = 1.048, ['P'] = 1.049, ['Q'] = 1.037, ['R'] = 1.008, ['S'] = 1.046, ['T'] = 0.997, ['V'] = 0.931, ['W'] = 0.904, ['Y'] = 0.929,
        };

        Assert.Multiple(() =>
        {
            Assert.That(ProteinScales.HoppWoods, Is.EquivalentTo(hw), "hw");
            Assert.That(ProteinScales.Emini, Is.EquivalentTo(em), "em");
            Assert.That(ProteinScales.Janin, Is.EquivalentTo(ja), "ja");
            Assert.That(ProteinScales.Flexibility, Is.EquivalentTo(flex), "Flex");
            Assert.That(ProteinScales.KyteDoolittle, Is.SameAs(ProteinPhysicochemistry.KyteDoolittle), "kd is canonical in PPC");
        });

        foreach (var scale in new[] { ProteinScales.HoppWoods, ProteinScales.Emini, ProteinScales.Janin, ProteinScales.Flexibility })
            Assert.Throws<NotSupportedException>(() => ((IDictionary<char, double>)scale)['A'] = 0, "read-only");
    }

    [Test]
    public void Scale_Null_Throws()
    {
        Assert.Throws<ArgumentNullException>(
            () => SequenceStatistics.CalculateHydrophobicityProfile(Seq, (IReadOnlyDictionary<char, double>)null!, 9));
    }

    #endregion
}
