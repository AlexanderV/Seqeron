using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PredictVariantEffectTests
{
    private static VariantDto Snp(int pos, string refA, string altA) => new(pos, refA, altA, "SNP", pos);

    [Test]
    public void PredictVariantEffect_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PredictVariantEffect(Snp(2, "A", "G"), "TTA", 2));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.PredictVariantEffect(null!, "TTA", 2));
        Assert.Throws<ArgumentException>(() => AnnotationTools.PredictVariantEffect(Snp(2, "A", "G"), "", 2));
    }

    [Test]
    public void PredictVariantEffect_SynonymousMissenseNonsense()
    {
        // Mirrors VariantCallerTests.PredictEffect_* cases.
        Assert.Multiple(() =>
        {
            // TTA (Leu) -> TTG (Leu): synonymous.
            Assert.That(AnnotationTools.PredictVariantEffect(Snp(2, "A", "G"), "TTA", 2).Effect, Is.EqualTo("Synonymous"));
            // ATG (Met) -> ACG (Thr): missense.
            Assert.That(AnnotationTools.PredictVariantEffect(Snp(1, "T", "C"), "ATG", 1).Effect, Is.EqualTo("Missense"));
            // TAT (Tyr) -> TAA (Stop): nonsense.
            Assert.That(AnnotationTools.PredictVariantEffect(Snp(2, "T", "A"), "TAT", 2).Effect, Is.EqualTo("Nonsense"));
        });
    }

    [Test]
    public void PredictVariantEffect_Indel_ReturnsFrameshift()
    {
        // Mirrors PredictEffect_Indel_ReturnsFrameshift.
        var result = AnnotationTools.PredictVariantEffect(new VariantDto(1, "-", "T", "Insertion", 1), "ATGC", 1);
        Assert.That(result.Effect, Is.EqualTo("Frameshift"));
    }

    [Test]
    public void PredictVariantEffect_OutOfRangePosition_ReturnsUnknown()
    {
        var result = AnnotationTools.PredictVariantEffect(Snp(99, "A", "G"), "ATG", 99);
        Assert.That(result.Effect, Is.EqualTo("Unknown"));
    }
}
