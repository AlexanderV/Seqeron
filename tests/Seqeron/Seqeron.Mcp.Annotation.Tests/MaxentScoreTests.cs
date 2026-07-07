using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class MaxentScoreTests
{
    [Test]
    public void MaxentScore_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.MaxentScore("CAGGUAAGU", "Donor"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.MaxentScore("", "Donor"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.MaxentScore("CAGGUAAGU", ""));
        // Unknown splice-site type is rejected.
        Assert.Throws<ArgumentException>(() => AnnotationTools.MaxentScore("CAGGUAAGU", "Enhancer"));
    }

    [Test]
    public void MaxentScore_Binding_InvokesSuccessfully()
    {
        // SpliceSitePredictor.CalculateMaxEntScore for a Donor sums log2(weight + 0.01) over the 9-mer PWM.
        // For a perfect CAGGUAAGU (every position weight 1.0): score = 9 * log2(1.01).
        var result = AnnotationTools.MaxentScore("CAGGUAAGU", "Donor");
        Assert.That(result.Score, Is.EqualTo(9.0 * Math.Log2(1.01)).Within(1e-9));
    }

    [Test]
    public void MaxentScore_NonConsensusScoresLower()
    {
        // A donor breaking the invariant positions scores far lower than the perfect consensus.
        var perfect = AnnotationTools.MaxentScore("CAGGUAAGU", "Donor").Score;
        var broken = AnnotationTools.MaxentScore("GGGGUAAGU", "Donor").Score;
        Assert.That(broken, Is.LessThan(perfect));
    }

    [Test]
    public void MaxentScore_TypeIsCaseInsensitive()
    {
        Assert.That(AnnotationTools.MaxentScore("CAGGUAAGU", "donor").Score,
            Is.EqualTo(AnnotationTools.MaxentScore("CAGGUAAGU", "Donor").Score));
    }
}
