using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ImpactLevelTests
{
    [Test]
    public void ImpactLevel_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ImpactLevel("MissenseVariant"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ImpactLevel(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ImpactLevel(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ImpactLevel("NotAConsequence"));
    }

    [Test]
    public void ImpactLevel_MapsConsequenceToImpact()
    {
        // Mirrors VariantAnnotator.GetImpactLevel switch mapping.
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.ImpactLevel("StopGained").Impact, Is.EqualTo("High"));
            Assert.That(AnnotationTools.ImpactLevel("FrameshiftVariant").Impact, Is.EqualTo("High"));
            Assert.That(AnnotationTools.ImpactLevel("MissenseVariant").Impact, Is.EqualTo("Moderate"));
            Assert.That(AnnotationTools.ImpactLevel("SynonymousVariant").Impact, Is.EqualTo("Low"));
            Assert.That(AnnotationTools.ImpactLevel("IntronVariant").Impact, Is.EqualTo("Modifier"));
        });
    }

    [Test]
    public void ImpactLevel_IsCaseInsensitive()
    {
        Assert.That(AnnotationTools.ImpactLevel("missensevariant").Impact, Is.EqualTo("Moderate"));
    }
}
