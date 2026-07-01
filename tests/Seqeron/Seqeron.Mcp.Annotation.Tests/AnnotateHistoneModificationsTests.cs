using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AnnotateHistoneModificationsTests
{
    private static HistoneModificationInputDto Mod(string mark, double signal) =>
        new(0, 100, mark, signal);

    [Test]
    public void AnnotateHistoneModifications_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.AnnotateHistoneModifications(null!));
        // Empty list is valid (yields no annotations).
        Assert.DoesNotThrow(() => AnnotationTools.AnnotateHistoneModifications(new List<HistoneModificationInputDto>()));
    }

    [Test]
    public void AnnotateHistoneModifications_Binding_InvokesSuccessfully()
    {
        // Present-mark Roadmap mapping (presence threshold = 0.5), per EpigeneticsAnalyzer.InferStateFromMark.
        var input = new List<HistoneModificationInputDto>
        {
            Mod("H3K4me3", 0.9),   // TssA -> ActivePromoter
            Mod("H3K4me1", 0.9),   // Enh (weak, no H3K27ac) -> WeakEnhancer
            Mod("H3K27ac", 0.9),   // active enhancer -> ActiveEnhancer
            Mod("H3K36me3", 0.9),  // Tx -> Transcribed
            Mod("H3K27me3", 0.9),  // ReprPC -> Repressed
            Mod("H3K9me3", 0.9),   // Het -> Heterochromatin
            Mod("H3K4me3", 0.3),   // below 0.5 threshold -> LowSignal
        };

        var result = AnnotationTools.AnnotateHistoneModifications(input);

        Assert.That(result.Annotations, Has.Count.EqualTo(7));
        Assert.Multiple(() =>
        {
            Assert.That(result.Annotations[0].PredictedState, Is.EqualTo("ActivePromoter"));
            Assert.That(result.Annotations[1].PredictedState, Is.EqualTo("WeakEnhancer"));
            Assert.That(result.Annotations[2].PredictedState, Is.EqualTo("ActiveEnhancer"));
            Assert.That(result.Annotations[3].PredictedState, Is.EqualTo("Transcribed"));
            Assert.That(result.Annotations[4].PredictedState, Is.EqualTo("Repressed"));
            Assert.That(result.Annotations[5].PredictedState, Is.EqualTo("Heterochromatin"));
            Assert.That(result.Annotations[6].PredictedState, Is.EqualTo("LowSignal"));
            // Coordinates and signal are preserved.
            Assert.That(result.Annotations[0].Start, Is.EqualTo(0));
            Assert.That(result.Annotations[0].End, Is.EqualTo(100));
            Assert.That(result.Annotations[0].Signal, Is.EqualTo(0.9));
        });
    }
}
