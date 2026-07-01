using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindPromoterMotifsTests
{
    [Test]
    public void FindPromoterMotifs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindPromoterMotifs("GGGGGTTGACAGGGGG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindPromoterMotifs(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindPromoterMotifs(null!));
    }

    [Test]
    public void FindPromoterMotifs_FullMinus35Consensus_ReturnsCorrectHit()
    {
        // Mirrors GenomeAnnotator_PromoterMotif_Tests.FindPromoterMotifs_FullMinus35Consensus_ReturnsCorrectHit.
        var result = AnnotationTools.FindPromoterMotifs("GGGGGTTGACAGGGGG");

        var minus35 = result.Motifs.Where(m => m.Type == "-35 box").ToList();
        Assert.That(minus35, Has.Count.EqualTo(4)); // TTGACA, TTGAC, TGACA, TTGA

        var full = minus35.Single(h => h.Sequence == "TTGACA");
        Assert.Multiple(() =>
        {
            Assert.That(full.Position, Is.EqualTo(5));
            Assert.That(full.Type, Is.EqualTo("-35 box"));
            Assert.That(full.Score, Is.EqualTo(1.0).Within(0.001));
        });
    }

    [Test]
    public void FindPromoterMotifs_FullMinus10Consensus_ReturnsCorrectHit()
    {
        // Mirrors FindPromoterMotifs_FullMinus10Consensus_ReturnsCorrectHit (TATAAT after 5 C's).
        var result = AnnotationTools.FindPromoterMotifs("CCCCCTATAATCCCCC");

        var minus10 = result.Motifs.Where(m => m.Type == "-10 box").ToList();
        Assert.That(minus10, Has.Count.EqualTo(4)); // TATAAT, TATAA, ATAAT, TATA

        var full = minus10.Single(h => h.Sequence == "TATAAT");
        Assert.Multiple(() =>
        {
            Assert.That(full.Position, Is.EqualTo(5));
            Assert.That(full.Type, Is.EqualTo("-10 box"));
            Assert.That(full.Score, Is.EqualTo(1.0).Within(0.001));
        });
    }

    [Test]
    public void FindPromoterMotifs_NoMotif_ReturnsEmpty()
    {
        var result = AnnotationTools.FindPromoterMotifs("GGGGGGGGGGGGGGGG");
        Assert.That(result.Motifs, Is.Empty);
    }
}
