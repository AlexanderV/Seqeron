using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindRibosomeBindingSitesTests
{
    // Mirrors GenomeAnnotator_Gene_Tests: minimal ORF = ATG + (aa-1)*3 A's + TAA.
    private static string CreateMinimalOrf(int aminoAcidCount = 100)
        => "ATG" + new string('A', (aminoAcidCount - 1) * 3) + "TAA";

    // Mirrors GenomeAnnotator_Gene_Tests.CreateSequenceWithSd: padding + SD + spacer + ORF.
    private static string CreateSequenceWithSd(string sdMotif, int distanceToStart, int orfLength = 100)
        => new string('C', 10) + sdMotif + new string('C', distanceToStart) + CreateMinimalOrf(orfLength);

    [Test]
    public void FindRibosomeBindingSites_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindRibosomeBindingSites(CreateSequenceWithSd("AGGAGG", 8)));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindRibosomeBindingSites(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindRibosomeBindingSites(null!));
    }

    [Test]
    public void FindRibosomeBindingSites_DetectsConsensusAggagg()
    {
        // Mirrors GenomeAnnotator_Gene_Tests.FindRibosomeBindingSites_DetectsConsensusAggagg.
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 8, orfLength: 100);

        var result = AnnotationTools.FindRibosomeBindingSites(
            sequence, upstreamWindow: 20, minDistance: 4, maxDistance: 15);

        var aggagg = result.Sites.Single(s => s.Sequence == "AGGAGG");
        Assert.Multiple(() =>
        {
            // SD begins right after the 10-C padding.
            Assert.That(aggagg.Position, Is.EqualTo(10));
            // score = motif.Length / 6.0 (GenomeAnnotator.cs L389); AGGAGG -> 1.0.
            Assert.That(aggagg.Score, Is.EqualTo(1.0).Within(1e-9));
        });
    }

    [Test]
    public void FindRibosomeBindingSites_TooClose_FullConsensusFiltered()
    {
        // Mirrors FindRibosomeBindingSites_TooClose_NotDetected: AGGAGG at aligned spacing 2 (< minDistance=4).
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 2, orfLength: 100);

        var result = AnnotationTools.FindRibosomeBindingSites(
            sequence, upstreamWindow: 20, minDistance: 4, maxDistance: 15);

        var fullConsensusAtInvalid = result.Sites.Where(s => s.Sequence == "AGGAGG" && s.Position == 10);
        Assert.That(fullConsensusAtInvalid.Count(), Is.EqualTo(0));
    }

    [Test]
    public void FindRibosomeBindingSites_AtMinDistance_Detected()
    {
        // Mirrors FindRibosomeBindingSites_AtMinDistance_Detected.
        string sequence = CreateSequenceWithSd("AGGAGG", distanceToStart: 4, orfLength: 100);

        var result = AnnotationTools.FindRibosomeBindingSites(
            sequence, upstreamWindow: 25, minDistance: 4, maxDistance: 15);

        Assert.That(result.Sites.Any(s => s.Sequence == "AGGAGG"), Is.True);
    }
}
