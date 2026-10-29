namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Findings of audit round 6 (the confirmation audit of round 5) on the CRISPR slice: B08 A6-1.
/// The alternative-PAM aggregates validated <c>minAlternativePamHitScore</c> inside the hit collection, which round 5 moved
/// behind the repeat-guide test, so whether an invalid cutoff was refused came to depend on the genome's content.
/// Each test failed on the code the auditor read (bec761ac).
/// TestSpec: TestSpecs/CRISPR-OFF-001.md ("Round 6").
/// </summary>
[TestFixture]
public class CrisprDesigner_Round6_Tests
{
    private const string Guide = "GACGCATAAAGATGAGACGC";

    // The guide three times: with maxOccurrences = 2 the repeat rule (CRISPOR MAXOCC) returns 0 before any hit is collected.
    private static readonly DnaSequence RepeatGenome = new(string.Join("TTTTT", Enumerable.Repeat(Guide + "TGG", 3)));

    private static readonly DnaSequence UniqueGenome = new("TTTTT" + Guide + "TGG" + "CCCCCACACACACACACAGTGTG");

    public static IEnumerable<TestCaseData> InvalidCutoffs()
    {
        foreach (var (cutoff, label) in new[] { (double.NaN, "NaN"), (-1.0, "minus_one"), (double.NegativeInfinity, "minus_infinity"), (-double.Epsilon, "minus_epsilon") })
            yield return new TestCaseData(cutoff).SetName($"cutoff_{label}");
    }

    [TestCaseSource(nameof(InvalidCutoffs))]
    public void MitAlternativePams_InvalidCutoff_IsRefused_ForARepeatGuideAsForAUniqueOne(double cutoff)
    {
        var onRepeat = Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(Guide, RepeatGenome, 4, CrisprSystemType.SpCas9, cutoff, maxOccurrences: 2));
        var onUnique = Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(Guide, UniqueGenome, 4, CrisprSystemType.SpCas9, cutoff, maxOccurrences: 2));

        Assert.Multiple(() =>
        {
            Assert.That(onRepeat!.ParamName, Is.EqualTo("minAlternativePamHitScore"));
            Assert.That(onUnique!.ParamName, Is.EqualTo("minAlternativePamHitScore"));
        });
    }

    [TestCaseSource(nameof(InvalidCutoffs))]
    public void CfdAlternativePams_InvalidCutoff_IsRefused_ForARepeatGuideAsForAUniqueOne(double cutoff)
    {
        var onRepeat = Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(Guide, RepeatGenome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, cutoff, maxOccurrences: 2));
        var onUnique = Assert.Throws<ArgumentOutOfRangeException>(() =>
            CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(Guide, UniqueGenome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, cutoff, maxOccurrences: 2));

        Assert.Multiple(() =>
        {
            Assert.That(onRepeat!.ParamName, Is.EqualTo("minAlternativePamHitScore"));
            Assert.That(onUnique!.ParamName, Is.EqualTo("minAlternativePamHitScore"));
        });
    }

    [TestCase(0.0)]
    [TestCase(1.0)]
    [TestCase(double.PositiveInfinity)]
    public void AlternativePams_ValidCutoff_StillScores_AndARepeatGuideStillReturnsZero(double cutoff)
    {
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(Guide, RepeatGenome, 4, CrisprSystemType.SpCas9, cutoff, maxOccurrences: 2), Is.EqualTo(0.0));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(Guide, RepeatGenome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, cutoff, maxOccurrences: 2), Is.EqualTo(0.0));
            Assert.That(CrisprDesigner.CalculateMitSpecificityScoreWithAlternativePams(Guide, UniqueGenome, 4, CrisprSystemType.SpCas9, cutoff), Is.InRange(0.0, 100.0));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(Guide, UniqueGenome, 4, CrisprSystemType.SpCas9, CfdAggregation.MitStyle, cutoff), Is.InRange(0.0, 100.0));
        });
    }
}
