using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindRepetitiveElementsTests
{
    [Test]
    public void FindRepetitiveElements_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindRepetitiveElements("ATTCGATTCGATTCG", minRepeatLength: 5));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindRepetitiveElements(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindRepetitiveElements(null!));
    }

    [Test]
    public void FindRepetitiveElements_ThreeAdjacentCopies_ReturnsSingleTandemRepeat()
    {
        // Mirrors GenomeAnnotator_FindRepetitiveElements_Tests: ATTCG x3 -> one tandem_repeat [0,15).
        var result = AnnotationTools.FindRepetitiveElements("ATTCGATTCGATTCG", minRepeatLength: 5, minCopies: 2);

        var tr = result.Repeats.Single(e => e.Type == "tandem_repeat" && e.Start == 0 && e.End == 15);
        Assert.That(tr.Sequence, Is.EqualTo("ATTCGATTCGATTCG"));
    }

    [Test]
    public void FindRepetitiveElements_Palindrome_ReturnsInvertedRepeat()
    {
        // Mirrors FindRepetitiveElements_ReverseComplementPalindrome_ReturnsInvertedRepeat: GAATTC (EcoRI).
        var result = AnnotationTools.FindRepetitiveElements("GAATTC", minRepeatLength: 3, minCopies: 2);

        var ir = result.Repeats.Single(e => e.Type == "inverted_repeat" && e.Start == 0 && e.End == 6);
        Assert.That(ir.Sequence, Is.EqualTo("GAATTC"));
    }

    [Test]
    public void FindRepetitiveElements_NonRepetitive_NoTandemRepeat()
    {
        // Mirrors FindRepetitiveElements_MotifAppearsOnce_NoTandemRepeatForThatMotif.
        var result = AnnotationTools.FindRepetitiveElements("ATTCGAAAAA", minRepeatLength: 5, minCopies: 2);

        Assert.That(result.Repeats.Any(e => e.Type == "tandem_repeat" && e.Sequence.Contains("ATTCGATTCG")), Is.False);
    }

    [Test]
    public void FindRepetitiveElements_SubStemsAndRunsOfN_OnlyMaximalAcgtElementsDelegated()
    {
        // EMBOSS palindrome -minpallen 3 -gaplimit 100 -nummismatches 0: one maximal stem AGGATCCT [0,8);
        // the NNNNNNNN tail forms neither a tandem unit nor stem pairs (ACGT only, B11 F14/F15).
        var result = AnnotationTools.FindRepetitiveElements("AGGATCCTNNNNNNNN", minRepeatLength: 3, minCopies: 2);

        Assert.That(result.Repeats.Select(e => (e.Start, e.End, e.Type, e.Sequence)),
            Is.EqualTo(new[] { (0, 8, "inverted_repeat", "AGGATCCT") }));
    }
}
