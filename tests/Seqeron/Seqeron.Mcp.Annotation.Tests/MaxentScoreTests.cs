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
        // Donor = MaxEntScan score5 (Yeo & Burge 2004): maxentpy score5('CAGGTAAGT') = 10.858313101356437
        // (score5.pl prints 10.86).
        var result = AnnotationTools.MaxentScore("CAGGUAAGU", "Donor");
        Assert.That(result.Score, Is.EqualTo(10.858313101356437).Within(1e-12));
    }

    [Test]
    public void MaxentScore_DonorNot9Nt_Throws()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.MaxentScore("CAGGUAAGUA", "Donor"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.MaxentScore("CAGGNAAGU", "Donor"));
    }

    [Test]
    public void MaxentScore_NonConsensusScoresLower()
    {
        // A donor breaking the invariant positions scores far lower than the perfect consensus.
        var perfect = AnnotationTools.MaxentScore("CAGGUAAGU", "Donor").Score;
        var broken = AnnotationTools.MaxentScore("CAGGCAAGU", "Donor").Score; // GC donor: 3.10 (maxentpy)
        Assert.That(broken, Is.EqualTo(3.104161570996358).Within(1e-12));
        Assert.That(broken, Is.LessThan(perfect));
    }

    [Test]
    public void MaxentScore_TypeIsCaseInsensitive()
    {
        Assert.That(AnnotationTools.MaxentScore("CAGGUAAGU", "donor").Score,
            Is.EqualTo(AnnotationTools.MaxentScore("CAGGUAAGU", "Donor").Score));
    }
}
