using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class PredictChromatinStateTests
{
    // Signature order: h3k4me3, h3k4me1, h3k27ac, h3k36me3, h3k27me3, h3k9me3.

    [Test]
    public void PredictChromatinState_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.PredictChromatinState(1, 0, 0, 0, 0, 0));
        // Signals must be in [0,1]; ValidateSignal throws ArgumentOutOfRangeException (an ArgumentException).
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.PredictChromatinState(-0.1, 0, 0, 0, 0, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.PredictChromatinState(0, 1.5, 0, 0, 0, 0));
    }

    [Test]
    public void PredictChromatinState_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.PredictChromatinState (ChromHMM present/absent, threshold 0.5).
        Assert.Multiple(() =>
        {
            // H3K4me3 alone -> active promoter.
            Assert.That(AnnotationTools.PredictChromatinState(1, 0, 0, 0, 0, 0).State,
                Is.EqualTo("ActivePromoter"));
            // H3K4me1 + H3K27ac -> active enhancer.
            Assert.That(AnnotationTools.PredictChromatinState(0, 1, 1, 0, 0, 0).State,
                Is.EqualTo("ActiveEnhancer"));
            // H3K4me1 without H3K27ac -> weak enhancer.
            Assert.That(AnnotationTools.PredictChromatinState(0, 1, 0, 0, 0, 0).State,
                Is.EqualTo("WeakEnhancer"));
            // H3K36me3 -> transcribed.
            Assert.That(AnnotationTools.PredictChromatinState(0, 0, 0, 1, 0, 0).State,
                Is.EqualTo("Transcribed"));
            // H3K27me3 -> Polycomb-repressed.
            Assert.That(AnnotationTools.PredictChromatinState(0, 0, 0, 0, 1, 0).State,
                Is.EqualTo("Repressed"));
            // H3K9me3 -> heterochromatin.
            Assert.That(AnnotationTools.PredictChromatinState(0, 0, 0, 0, 0, 1).State,
                Is.EqualTo("Heterochromatin"));
        });
    }

    [Test]
    public void PredictChromatinState_BivalentAndLowSignal()
    {
        // H3K4me3 + H3K27me3 -> bivalent promoter (checked before ActivePromoter).
        Assert.That(AnnotationTools.PredictChromatinState(1, 0, 0, 0, 1, 0).State,
            Is.EqualTo("BivalentPromoter"));
        // H3K4me1 + H3K27me3 -> bivalent enhancer.
        Assert.That(AnnotationTools.PredictChromatinState(0, 1, 0, 0, 1, 0).State,
            Is.EqualTo("BivalentEnhancer"));
        // No mark present -> low signal.
        Assert.That(AnnotationTools.PredictChromatinState(0, 0, 0, 0, 0, 0).State,
            Is.EqualTo("LowSignal"));
    }
}
