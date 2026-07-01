using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ClassifyMutationTests
{
    private static VariantDto Snp(string refA, string altA) => new(0, refA, altA, "SNP", 0);

    [Test]
    public void ClassifyMutation_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ClassifyMutation(Snp("A", "G")));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.ClassifyMutation(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyMutation(new VariantDto(0, "", "G", "SNP", 0)));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ClassifyMutation(new VariantDto(0, "A", "", "SNP", 0)));
    }

    [Test]
    public void ClassifyMutation_PurineOrPyrimidineSwap_IsTransition()
    {
        // Mirrors VariantCaller_CallVariants_Tests: A->G and C->T are transitions.
        Assert.Multiple(() =>
        {
            Assert.That(AnnotationTools.ClassifyMutation(Snp("A", "G")).MutationType, Is.EqualTo("Transition"));
            Assert.That(AnnotationTools.ClassifyMutation(Snp("C", "T")).MutationType, Is.EqualTo("Transition"));
            Assert.That(AnnotationTools.ClassifyMutation(Snp("G", "A")).MutationType, Is.EqualTo("Transition"));
        });
    }

    [Test]
    public void ClassifyMutation_PurinePyrimidineSwap_IsTransversion()
    {
        // A->C crosses purine<->pyrimidine, a transversion.
        Assert.That(AnnotationTools.ClassifyMutation(Snp("A", "C")).MutationType, Is.EqualTo("Transversion"));
    }

    [Test]
    public void ClassifyMutation_NonSnp_IsOther()
    {
        // Non-SNP variant types classify as Other.
        var result = AnnotationTools.ClassifyMutation(new VariantDto(3, "-", "T", "Insertion", 3));
        Assert.That(result.MutationType, Is.EqualTo("Other"));
    }
}
