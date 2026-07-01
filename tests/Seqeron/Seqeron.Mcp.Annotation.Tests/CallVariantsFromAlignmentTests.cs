using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class CallVariantsFromAlignmentTests
{
    [Test]
    public void CallVariantsFromAlignment_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.CallVariantsFromAlignment("AT-GC", "ATXGC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CallVariantsFromAlignment("", "ATXGC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CallVariantsFromAlignment(null!, "ATXGC"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.CallVariantsFromAlignment("AT-GC", ""));
        // Different aligned lengths -> ArgumentException (surfaced when the tool materializes the list).
        Assert.Throws<ArgumentException>(() => AnnotationTools.CallVariantsFromAlignment("ATGC", "ATG"));
    }

    [Test]
    public void CallVariantsFromAlignment_Insertion_DetectsIt()
    {
        // Reference gap at column 2 -> Insertion of 'X' at reference position 2.
        var result = AnnotationTools.CallVariantsFromAlignment("AT-GC", "ATXGC");

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        var v = result.Variants[0];
        Assert.Multiple(() =>
        {
            Assert.That(v.Type, Is.EqualTo("Insertion"));
            Assert.That(v.ReferenceAllele, Is.EqualTo("-"));
            Assert.That(v.AlternateAllele, Is.EqualTo("X"));
            Assert.That(v.Position, Is.EqualTo(2));
        });
    }

    [Test]
    public void CallVariantsFromAlignment_Deletion_DetectsIt()
    {
        // Query gap at column 2 -> Deletion of 'G' at reference position 2.
        var result = AnnotationTools.CallVariantsFromAlignment("ATGC", "AT-C");

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        var v = result.Variants[0];
        Assert.Multiple(() =>
        {
            Assert.That(v.Type, Is.EqualTo("Deletion"));
            Assert.That(v.ReferenceAllele, Is.EqualTo("G"));
            Assert.That(v.AlternateAllele, Is.EqualTo("-"));
            Assert.That(v.Position, Is.EqualTo(2));
        });
    }
}
