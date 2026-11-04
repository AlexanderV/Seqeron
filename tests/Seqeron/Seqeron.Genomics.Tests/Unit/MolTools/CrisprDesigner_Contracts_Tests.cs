namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Argument and range contracts of <see cref="CrisprDesigner"/> (audit round 4, B08 H14-H20): what the public
/// methods accept, which exception they raise, and in which order. Numeric behaviour is pinned elsewhere; these tests
/// pin that malformed input fails loudly and uniformly instead of returning a plausible-looking number.
/// TestSpec: TestSpecs/CRISPR-OFF-001.md and CRISPR-GUIDE-001.md ("Round 4 contracts").
/// </summary>
[TestFixture]
public class CrisprDesigner_Contracts_Tests
{
    private const string Guide = "GACGCATAAAGATGAGACGC";
    private static readonly DnaSequence Genome = new("CTTGCGTCTCATCTTTATGCGTAAGGTTTTTT");

    // ------------------------------------------------------------------------------------------------
    // H15 / H14: per-hit scores are finite and non-negative; the Parkinson error blames the right thing.
    // Reference: crispor.py calcMitGuideScore / calcCfdGuideScore are module-private and only ever fed sums of
    // non-negative hit scores; fed -2, calcCfdGuideScore returns -5000.0 and calcMitGuideScore(-100) divides by 0.
    // ------------------------------------------------------------------------------------------------

    [TestCase(-150.0)]
    [TestCase(-100.0)]
    [TestCase(-0.5)]
    [TestCase(double.NaN)]
    [TestCase(double.PositiveInfinity)]
    [TestCase(double.NegativeInfinity)]
    public void MitSpecificity_List_RefusesAScoreThatIsNotFiniteAndNonNegative(double score)
    {
        // Before: [-150] -> -200, [-100] -> Infinity, [NaN] -> NaN, all returned from a method documented as [0, 100].
        var e = Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.CalculateMitSpecificityScore(new[] { 10.0, score }));
        Assert.That(e!.ParamName, Is.EqualTo("singleHitScores"));
    }

    [Test]
    public void MitSpecificity_List_KeepsTheDocumentedRangeForValidScores()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateMitSpecificityScore(Array.Empty<double>()), Is.EqualTo(100.0));
            Assert.That(CrisprDesigner.CalculateMitSpecificityScore(new[] { 0.0 }), Is.EqualTo(100.0));
            Assert.That(CrisprDesigner.CalculateMitSpecificityScore(new[] { 100.0 }), Is.EqualTo(50.0));
            Assert.That(CrisprDesigner.CalculateMitSpecificityScore(new[] { 1e300, 1e300 }), Is.InRange(0.0, 100.0));
        });
    }

    [TestCase(CfdAggregation.MitStyle)]
    [TestCase(CfdAggregation.ParkinsonNormalized)]
    public void CfdSpecificity_List_RefusesNegativeAndNonFiniteScores_WithAMessageThatIsTrue(CfdAggregation aggregation)
    {
        foreach (double bad in new[] { -2.0, double.NaN, double.PositiveInfinity })
        {
            var e = Assert.Throws<ArgumentOutOfRangeException>(() =>
                CrisprDesigner.CalculateCfdSpecificityScore(new[] { bad }, aggregation));
            Assert.That(e!.ParamName, Is.EqualTo("cfdScores"));
            // The old message asserted "undefined for a CFD sum of 0" for a sum of -2 (and for NaN).
            Assert.That(e.Message, Does.Contain("finite and non-negative"));
            Assert.That(e.Message, Does.Not.Contain("sum of 0"));
        }
    }

    [Test]
    public void CfdSpecificity_ParkinsonNormalized_ZeroSum_BlamesTheArgumentTheCallerPassed()
    {
        var emptyList = Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.CalculateCfdSpecificityScore(Array.Empty<double>(), CfdAggregation.ParkinsonNormalized));
        Assert.That(emptyList!.ParamName, Is.EqualTo("cfdScores"));
        Assert.That(emptyList.Message, Does.Contain("sum of 0"));

        // A genome with no off-target has CFD sum 0 as well; "cfdScores" is not a parameter of the genome methods.
        var noPam = new DnaSequence(new string('T', 40));
        var genomeOverload = Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.CalculateCfdSpecificityScore(Guide, noPam, 4, CrisprSystemType.SpCas9,
                CfdAggregation.ParkinsonNormalized));
        Assert.That(genomeOverload!.ParamName, Is.EqualTo("aggregation"));
        var alt = Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(Guide, noPam, 4, CrisprSystemType.SpCas9,
                CfdAggregation.ParkinsonNormalized));
        Assert.That(alt!.ParamName, Is.EqualTo("aggregation"));
    }

    // ------------------------------------------------------------------------------------------------
    // H16: the four genome aggregates share one argument prologue, so one malformed call fails the same way in all.
    // Order: null/empty guide, null genome, maxMismatches 0-5, the model refusal (PAM-before-target), then the guide
    // length and alphabet; the MIT pair used to run the refusal and the interval check before everything else.
    // ------------------------------------------------------------------------------------------------

    private static IEnumerable<(string Name, Action Call)> Aggregates(
        string? guide, DnaSequence? genome, int maxMismatches, CrisprSystemType system) =>
    [
        ("MIT", () => CrisprDesigner.CalculateMitSpecificityScore(guide!, genome!, maxMismatches, system)),
        ("MIT+alt", () => CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(guide!, genome!, maxMismatches, system)),
        ("CFD", () => CrisprDesigner.CalculateCfdSpecificityScore(guide!, genome!, maxMismatches, system)),
        ("CFD+alt", () => CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(guide!, genome!, maxMismatches, system)),
    ];

    private static (Type Type, string? Param) Outcome(Action call)
    {
        var ex = Assert.Catch(() => call());
        return (ex!.GetType(), (ex as ArgumentException)?.ParamName);
    }

    [Test]
    public void GenomeAggregates_FailIdenticallyOnTheSameMalformedCall()
    {
        var cases = new (string Label, string? Guide, DnaSequence? Genome, int MaxMm, CrisprSystemType System, Type Expected, string Param)[]
        {
            ("null guide, Cas12a", null, Genome, 4, CrisprSystemType.Cas12a, typeof(ArgumentNullException), "guideSequence"),
            ("empty guide, Cas12a", "", Genome, 4, CrisprSystemType.Cas12a, typeof(ArgumentException), "guideSequence"),
            ("null genome, Cas12a", Guide, null, 4, CrisprSystemType.Cas12a, typeof(ArgumentNullException), "genome"),
            ("maxMismatches 9, Cas12a", Guide, Genome, 9, CrisprSystemType.Cas12a, typeof(ArgumentOutOfRangeException), "maxMismatches"),
            ("valid arguments, Cas12a", Guide + "ACG", Genome, 4, CrisprSystemType.Cas12a, typeof(ArgumentException), "systemType"),
            ("wrong length, SpCas9", "ACGT", Genome, 4, CrisprSystemType.SpCas9, typeof(ArgumentException), "guideSequence"),
            ("junk alphabet, SpCas9", "GACGCATAAAGATGAGACGN", Genome, 4, CrisprSystemType.SpCas9, typeof(ArgumentException), "guideSequence"),
        };

        Assert.Multiple(() =>
        {
            foreach (var c in cases)
                foreach (var (name, call) in Aggregates(c.Guide, c.Genome, c.MaxMm, c.System))
                {
                    var outcome = Outcome(call);
                    Assert.That(outcome, Is.EqualTo((c.Expected, (string?)c.Param)), $"{c.Label} / {name}");
                }
        });
    }

    // ------------------------------------------------------------------------------------------------
    // H17: EvaluateGuideRna enforces the system's guide length, as the off-target scans do.
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void EvaluateGuideRna_WrongLength_IsRefused_NotScoredAsPerfect()
    {
        var short4 = Assert.Throws<ArgumentException>(() => CrisprDesigner.EvaluateGuideRna("ACGT"));
        Assert.That(short4!.ParamName, Is.EqualTo("guideSequence"));
        Assert.That(short4.Message, Does.Contain("Guide length (4)").And.Contain("expected length (20)"));
        Assert.Throws<ArgumentException>(() => CrisprDesigner.EvaluateGuideRna("ACGTACGT", CrisprSystemType.LbCas12a));
        Assert.Throws<ArgumentException>(() => CrisprDesigner.EvaluateGuideRna(new string('A', 21)));
    }

    [Test]
    public void EvaluateGuideRna_AcceptsExactlyTheGuideLengthOfEverySystem()
    {
        foreach (var type in Enum.GetValues<CrisprSystemType>())
        {
            int length = CrisprDesigner.GetSystem(type).GuideLength;
            string guide = string.Concat(Enumerable.Repeat("ACGGACTAGC", 3))[..length];

            Assert.DoesNotThrow(() => CrisprDesigner.EvaluateGuideRna(guide, type), type.ToString());
            Assert.Throws<ArgumentException>(() => CrisprDesigner.EvaluateGuideRna(guide + "A", type), type.ToString());
            Assert.Throws<ArgumentException>(() => CrisprDesigner.EvaluateGuideRna(guide[..^1], type), type.ToString());
        }
    }

    // ------------------------------------------------------------------------------------------------
    // H18: GuideRnaParameters are validated; the documented -1 transcription and default(...) are refused loudly.
    // ------------------------------------------------------------------------------------------------

    public static IEnumerable<TestCaseData> InvalidParameters()
    {
        var d = GuideRnaParameters.Default;
        yield return new TestCaseData(d with { MinGcContent = 80, MaxGcContent = 20 }).SetName("MinGc 80 > MaxGc 20");
        yield return new TestCaseData(d with { MinGcContent = -1 }).SetName("MinGc below 0");
        yield return new TestCaseData(d with { MaxGcContent = 101 }).SetName("MaxGc above 100");
        yield return new TestCaseData(d with { MinGcContent = double.NaN }).SetName("MinGc NaN");
        yield return new TestCaseData(d with { MinScore = double.NaN }).SetName("MinScore NaN");
        yield return new TestCaseData(d with { MaxSelfComplementaryStems = -1 }).SetName("stem limit -1 (CHOPCHOP's sentinel)");
        yield return new TestCaseData(default(GuideRnaParameters)).SetName("default(GuideRnaParameters)");
    }

    [TestCaseSource(nameof(InvalidParameters))]
    public void GuideRnaParameters_Invalid_IsRefusedByDesignAndEvaluation(GuideRnaParameters parameters)
    {
        const string text = "ACGTGACGCATAAAGATGAGACGCAGGTGA";
        var design = Assert.Catch<ArgumentException>(() =>
            CrisprDesigner.DesignGuideRnas(new DnaSequence(text), 0, text.Length - 1, CrisprSystemType.SpCas9, parameters));
        var evaluate = Assert.Catch<ArgumentException>(() =>
            CrisprDesigner.EvaluateGuideRna(Guide, CrisprSystemType.SpCas9, parameters));

        Assert.That(design!.ParamName, Is.EqualTo("parameters"));
        Assert.That(evaluate!.ParamName, Is.EqualTo("parameters"));
    }

    [Test]
    public void GuideRnaParameters_ParameterlessConstructor_IsTheDocumentedDefault()
    {
        var p = new GuideRnaParameters();

        Assert.That(p, Is.EqualTo(GuideRnaParameters.Default));
        Assert.That(p.ScoreGcContent, Is.True);
        // 0 and null are legitimate stem limits; a window equal at both ends is legitimate too.
        Assert.DoesNotThrow(() => CrisprDesigner.EvaluateGuideRna(Guide, CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { MaxSelfComplementaryStems = 0 }));
        Assert.DoesNotThrow(() => CrisprDesigner.EvaluateGuideRna(Guide, CrisprSystemType.SpCas9,
            GuideRnaParameters.Default with { MinGcContent = 50, MaxGcContent = 50 }));
    }

    // ------------------------------------------------------------------------------------------------
    // H19: a null reference is a missing argument, an empty string is a wrong value.
    // ------------------------------------------------------------------------------------------------

    public static IEnumerable<TestCaseData> StringEntryPoints()
    {
        var genome = Genome;
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.EvaluateGuideRna(g!)), "guideSequence").SetName("EvaluateGuideRna");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.FindOffTargets(g!, genome).ToList()), "guideSequence").SetName("FindOffTargets");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.FindOffTargetsWithPam(g!, genome).ToList()), "guideSequence").SetName("FindOffTargetsWithPam");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.FindOffTargetsWithBulges(g!, genome).ToList()), "guideSequence").SetName("FindOffTargetsWithBulges");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.CalculateHeuristicSpecificityScore(g!, genome)), "guideSequence").SetName("CalculateHeuristicSpecificityScore");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.CalculateMitSpecificityScore(g!, genome)), "guideSequence").SetName("CalculateMitSpecificityScore");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.CalculateCfdSpecificityScore(g!, genome)), "guideSequence").SetName("CalculateCfdSpecificityScore");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.CalculateCfdScore(g!, Guide, "GG")), "sgRna20").SetName("CalculateCfdScore");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.CalculateMitHitScore(g!, Guide)), "guide20").SetName("CalculateMitHitScore");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.CalculateSaCas9HitScore(g!, new string('A', 21))), "saGuide21").SetName("CalculateSaCas9HitScore");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.GetGrafMotif(g!)), "guideSequence").SetName("GetGrafMotif");
        yield return new TestCaseData((Func<string?, object?>)(g => CrisprDesigner.CountSelfComplementaryStems(g!)), "guideSequence").SetName("CountSelfComplementaryStems");
    }

    [TestCaseSource(nameof(StringEntryPoints))]
    public void StringGuideEntryPoints_NullIsAMissingArgument_EmptyIsAWrongValue(Func<string?, object?> call, string param)
    {
        var nullEx = Assert.Throws<ArgumentNullException>(() => call(null));
        Assert.That(nullEx!.ParamName, Is.EqualTo(param));
        var emptyEx = Assert.Throws<ArgumentException>(() => call(""));   // exact type: not ArgumentNullException
        Assert.That(emptyEx!.ParamName, Is.EqualTo(param));
    }

    // ------------------------------------------------------------------------------------------------
    // H20 / H53: FindPamSites(string) validates its system at the call, like the DnaSequence overload; a null, empty or
    // shorter-than-the-PAM string is "no sites" (PAM_Site_Detection.md section 3.3).
    // ------------------------------------------------------------------------------------------------

    [Test]
    public void FindPamSites_UnknownSystem_ThrowsAtTheCall_OnBothOverloads()
    {
        var unknown = (CrisprSystemType)99;

        // No ToList(): an iterator would defer the throw to the first MoveNext.
        Assert.Throws<ArgumentException>(() => CrisprDesigner.FindPamSites("ACGTACGTACGTACGTACGTAGG", unknown));
        Assert.Throws<ArgumentException>(() => CrisprDesigner.FindPamSites(new DnaSequence("ACGTACGTACGTACGTACGTAGG"), unknown));
        Assert.Throws<ArgumentException>(() => CrisprDesigner.FindPamSites((string)null!, unknown), "validated before the null shortcut");
        Assert.Throws<ArgumentException>(() => CrisprDesigner.FindPamSites("", unknown));
    }

    [Test]
    public void FindPamSites_StringOverload_NullEmptyAndTooShort_YieldNoSites()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.FindPamSites((string)null!, CrisprSystemType.SpCas9), Is.Empty);
            Assert.That(CrisprDesigner.FindPamSites("", CrisprSystemType.SpCas9), Is.Empty);
            foreach (var type in Enum.GetValues<CrisprSystemType>())
            {
                var system = CrisprDesigner.GetSystem(type);
                string tooShort = new string('G', system.PamSequence.Length - 1);
                Assert.That(CrisprDesigner.FindPamSites(tooShort, type), Is.Empty, $"{type}: shorter than the PAM motif");
                Assert.That(CrisprDesigner.FindPamSites(new string('G', system.GuideLength + system.PamSequence.Length - 1), type),
                    Is.Empty, $"{type}: one base short of PAM + guide");
            }
        });
    }
}
