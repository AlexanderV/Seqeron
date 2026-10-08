using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class ComplexityKmerEntropyTests
{
    [Test]
    public void ComplexityKmerEntropy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.ComplexityKmerEntropy("ATGCGATCGATCG", 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityKmerEntropy("", 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityKmerEntropy(null!, 2));
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityKmerEntropy("ATGC", 0)); // k < 1
        Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityKmerEntropy("ZZZCGA", 2)); // invalid DNA
    }

    [Test]
    public void ComplexityKmerEntropy_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.ComplexityKmerEntropy("ATGCGATCGATCG", 2);
        Assert.That(result.Entropy, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.K, Is.EqualTo(2));

        // Low complexity (repetitive) sequence has lower entropy
        var lowComplexity = SequenceTools.ComplexityKmerEntropy("AAAAAAAAAAAAAAAA", 2);
        var highComplexity = SequenceTools.ComplexityKmerEntropy("ATGCATGCATGCATGC", 2);
        Assert.That(highComplexity.Entropy, Is.GreaterThan(lowComplexity.Entropy));
    }

    // R entropy 1.3.2 entropy.MillerMadow(c(3,2), unit="log2") = 1.1152200985435652; Grassberger (2003) eq. via mpmath
    // = 1.2692841903863027; BBTools EntropyTracker.calcEntropy("ATATAT", k=2) = 0.41816565 (= 0.9709505944546686 / log2 5).
    [Test]
    public void ComplexityKmerEntropy_CorrectionAndNormalize_MatchReferences()
    {
        var mm = SequenceTools.ComplexityKmerEntropy("ATATAT", 2, "millerMadow");
        var gr = SequenceTools.ComplexityKmerEntropy("ATATAT", 2, "Grassberger");
        var norm = SequenceTools.ComplexityKmerEntropy("ATATAT", 2, normalize: true);
        var plain = SequenceTools.ComplexityKmerEntropy("ATATAT", 2);
        Assert.Multiple(() =>
        {
            Assert.That(mm.Entropy, Is.EqualTo(1.1152200985435652).Within(1e-12));
            Assert.That(mm.Correction, Is.EqualTo("millerMadow"));
            Assert.That(gr.Entropy, Is.EqualTo(1.2692841903863027).Within(1e-12));
            Assert.That(gr.Correction, Is.EqualTo("grassberger"));
            Assert.That((float)norm.Entropy, Is.EqualTo(0.41816565f));
            Assert.That(norm.Normalized, Is.True);
            Assert.That(plain.Entropy, Is.EqualTo(0.9709505944546686).Within(1e-15));
            Assert.That(plain.Correction, Is.EqualTo("none"));
            Assert.Throws<ArgumentException>(() => SequenceTools.ComplexityKmerEntropy("ATATAT", 2, "nsb"));
        });
    }
}
