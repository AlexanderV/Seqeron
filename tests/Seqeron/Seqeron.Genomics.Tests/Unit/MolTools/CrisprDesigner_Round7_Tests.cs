namespace Seqeron.Genomics.Tests.Unit.MolTools;

/// <summary>
/// Findings of audit round 7 (the confirmation audit of round 6) on the CRISPR slice: B08 A7-1.
/// Round 6 made the alternative-PAM cutoff independent of the genome's content (A6-1); the CFD aggregates still validated their
/// <c>aggregation</c> argument only after the repeat-guide shortcut, in the last line of <c>AggregateCfd</c>.
/// Each test failed on the code the auditor read (cbcecfd3).
/// TestSpec: TestSpecs/CRISPR-OFF-001.md ("Round 7").
/// </summary>
[TestFixture]
public class CrisprDesigner_Round7_Tests
{
    private const string Guide = "GACGCATAAAGATGAGACGC";

    // The guide three times: with maxOccurrences = 2 the repeat rule (CRISPOR MAXOCC) returns 0 before any hit is collected.
    private static readonly DnaSequence RepeatGenome = new(string.Join("TTTTT", Enumerable.Repeat(Guide + "TGG", 3)));

    private static readonly DnaSequence UniqueGenome = new("TTTTT" + Guide + "TGG" + "CCCCCACACACACACACAGTGTG");

    [TestCase(99)]
    [TestCase(-1)]
    [TestCase(2)]
    [TestCase(int.MaxValue)]
    [TestCase(int.MinValue)]
    public void CfdAggregates_UndefinedAggregation_IsRefused_ForARepeatGuideAsForAUniqueOne(int aggregationValue)
    {
        var aggregation = (CfdAggregation)aggregationValue;

        Assert.Multiple(() =>
        {
            foreach (var (label, genome) in new[] { ("repeat", RepeatGenome), ("unique", UniqueGenome) })
            {
                var plain = Assert.Throws<ArgumentOutOfRangeException>(() =>
                    CrisprDesigner.CalculateCfdSpecificityScore(Guide, genome, 4, CrisprSystemType.SpCas9, aggregation, maxOccurrences: 2), $"plain, {label}");
                var alt = Assert.Throws<ArgumentOutOfRangeException>(() =>
                    CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(Guide, genome, 4, CrisprSystemType.SpCas9, aggregation, maxOccurrences: 2), $"alt-PAM, {label}");

                Assert.That(plain!.ParamName, Is.EqualTo("aggregation"), $"plain, {label}");
                Assert.That(alt!.ParamName, Is.EqualTo("aggregation"), $"alt-PAM, {label}");
            }
        });
    }

    [TestCase(CfdAggregation.MitStyle)]
    [TestCase(CfdAggregation.ParkinsonNormalized)]
    public void CfdAggregates_DefinedAggregation_ARepeatGuideStillReturnsZero(CfdAggregation aggregation)
    {
        Assert.Multiple(() =>
        {
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScore(Guide, RepeatGenome, 4, CrisprSystemType.SpCas9, aggregation, maxOccurrences: 2), Is.EqualTo(0.0));
            Assert.That(CrisprDesigner.CalculateCfdSpecificityScoreWithAlternativePams(Guide, RepeatGenome, 4, CrisprSystemType.SpCas9, aggregation, maxOccurrences: 2), Is.EqualTo(0.0));
        });
    }
}
