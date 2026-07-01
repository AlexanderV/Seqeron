using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class GenerateSeedVariantsTests
{
    [Test]
    public void GenerateSeedVariants_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.GenerateSeedVariants("GAGGUAG"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.GenerateSeedVariants(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.GenerateSeedVariants(null!));
    }

    [Test]
    public void GenerateSeedVariants_Binding_InvokesSuccessfully()
    {
        // MiRnaAnalyzer.GenerateSeedVariants yields the original, then every single-position
        // substitution over A/C/G/U. For "AC": original + 2 positions * 3 substitutions = 7.
        var result = AnnotationTools.GenerateSeedVariants("AC");

        Assert.That(result.Variants, Is.EqualTo(new[]
        {
            "AC",                 // original
            "CC", "GC", "UC",     // position 0 substitutions (A -> C/G/U)
            "AA", "AG", "AU",     // position 1 substitutions (C -> A/G/U)
        }));
    }

    [Test]
    public void GenerateSeedVariants_CountFormula()
    {
        // A length-7 seed yields 1 + 7*3 = 22 variants (all unique for a distinct seed).
        var result = AnnotationTools.GenerateSeedVariants("GAGGUAG");
        Assert.Multiple(() =>
        {
            Assert.That(result.Variants, Has.Count.EqualTo(22));
            Assert.That(result.Variants[0], Is.EqualTo("GAGGUAG"));
        });
    }
}
