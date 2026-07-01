using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ClassifyVariantTests
{
    [Test]
    public void ClassifyVariant_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ClassifyVariant("A", "G"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyVariant("", "G"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyVariant("A", ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyVariant(null!, "G"));
    }

    [Test]
    public void ClassifyVariant_MatchesVariantAnnotatorClassification()
    {
        // Mirrors VariantAnnotatorTests ClassifyVariant cases.
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.ClassifyVariant("A", "G").VariantType, Is.EqualTo("SNV"));
            Assert.That(AnnotationTools.ClassifyVariant("A", "ACGT").VariantType, Is.EqualTo("Insertion"));
            Assert.That(AnnotationTools.ClassifyVariant("ACGT", "A").VariantType, Is.EqualTo("Deletion"));
            Assert.That(AnnotationTools.ClassifyVariant("AC", "GT").VariantType, Is.EqualTo("MNV"));
            Assert.That(AnnotationTools.ClassifyVariant("ACG", "TT").VariantType, Is.EqualTo("Indel"));
        });
    }
}
