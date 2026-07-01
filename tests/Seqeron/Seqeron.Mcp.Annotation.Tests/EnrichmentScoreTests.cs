using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class EnrichmentScoreTests
{
    private static List<string> Ranked() => new() { "A", "B", "C", "D" };

    [Test]
    public void EnrichmentScore_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.EnrichmentScore(Ranked(), new List<string> { "A" }));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.EnrichmentScore(null!, new List<string> { "A" }));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.EnrichmentScore(Ranked(), null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.EnrichmentScore(new List<string>(), new List<string> { "A" }));
        Assert.Throws<ArgumentException>(() => AnnotationTools.EnrichmentScore(Ranked(), new List<string>()));
    }

    [Test]
    public void EnrichmentScore_Binding_InvokesSuccessfully()
    {
        // TranscriptomeAnalyzer.CalculateEnrichmentScore: running sum +1/hits on hit, -1/misses on miss;
        // returns max deviation. Set {A,B} at the top of [A,B,C,D]: +0.5,+1.0,+0.5,0 -> ES = 1.0.
        var result = AnnotationTools.EnrichmentScore(Ranked(), new List<string> { "A", "B" });
        Assert.That(result.Score, Is.EqualTo(1.0).Within(1e-9));
    }

    [Test]
    public void EnrichmentScore_BottomRankedSet_NegativeScore()
    {
        // Set {C,D} at the bottom: -0.5,-1.0,-0.5,0 -> max deviation is -1.0.
        var result = AnnotationTools.EnrichmentScore(Ranked(), new List<string> { "C", "D" });
        Assert.That(result.Score, Is.EqualTo(-1.0).Within(1e-9));
    }

    [Test]
    public void EnrichmentScore_AllHitsOrNoHits_ReturnsZero()
    {
        // No gene of the set is in the ranked list (missCount == n, hitCount 0) -> 0.
        var noHit = AnnotationTools.EnrichmentScore(Ranked(), new List<string> { "Z" });
        Assert.That(noHit.Score, Is.EqualTo(0.0));
    }
}
