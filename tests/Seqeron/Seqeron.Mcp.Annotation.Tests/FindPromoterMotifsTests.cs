using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

/// <summary>find_promoter_motifs — delegation to GenomeAnnotator.FindPromoterMotifs (σ70 pairs, Harley &amp; Reynolds 1987).</summary>
[TestFixture]
public class FindPromoterMotifsTests
{
    private const string Consensus = "GG" + "TTGACA" + "CCCCCCCCCCCCCCCCC" + "TATAAT" + "GG";

    [Test]
    public void FindPromoterMotifs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindPromoterMotifs(Consensus));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindPromoterMotifs(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindPromoterMotifs(null!));
    }

    [Test]
    public void FindPromoterMotifs_ConsensusPromoter_ReturnsBothBoxes()
    {
        var result = AnnotationTools.FindPromoterMotifs(Consensus);

        Assert.That(result.Motifs.Select(m => (m.Position, m.Type, m.Sequence, m.Score)), Is.EqualTo(new[]
        {
            (2, "-35 box", "TTGACA", 1.0),
            (25, "-10 box", "TATAAT", 1.0),
        }));
    }

    [Test]
    public void FindPromoterMotifs_UnpairedBox_ReturnsEmpty()
    {
        Assert.That(AnnotationTools.FindPromoterMotifs("GGGGGTTGACAGGGGG").Motifs, Is.Empty);
        Assert.That(AnnotationTools.FindPromoterMotifs("CCCCCTATAATCCCCC").Motifs, Is.Empty);
    }

    [Test]
    public void FindPromoterMotifs_DelegatesToGenomeAnnotator()
    {
        const string seq = "CCGATATTTTAAAAGTTTTTTTAAAGACTAACTATAAATTTGACATTATG";
        var expected = GenomeAnnotator.FindPromoterMotifs(seq).ToList();
        var actual = AnnotationTools.FindPromoterMotifs(seq).Motifs.Select(m => (m.Position, m.Type, m.Sequence, m.Score)).ToList();
        Assert.That(actual, Is.EqualTo(expected));
        Assert.That(actual, Has.Count.EqualTo(6));
    }

    [Test]
    public void FindPromoterMotifs_NoMotif_ReturnsEmpty()
    {
        var result = AnnotationTools.FindPromoterMotifs("GGGGGGGGGGGGGGGG");
        Assert.That(result.Motifs, Is.Empty);
    }
}
