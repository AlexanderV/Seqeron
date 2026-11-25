using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AssessHairpinByMfeTests
{
    private const string Mir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";
    private const string Unstructured = "GAUCAGUCAUGACUAGCUAGAUUGACAAGCUUACGAUCAGGACUUAGCAUCGUAGCAUGACUAG";

    [Test]
    public void AssessHairpinByMfe_Schema_ValidatesCorrectly()
    {
        Assert.Throws<ArgumentException>(() => AnnotationTools.AssessHairpinByMfe(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AssessHairpinByMfe(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.AssessHairpinByMfe(Mir21, double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.AssessHairpinByMfe(Mir21, -0.1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.AssessHairpinByMfe(Mir21, 0.85, -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.AssessHairpinByMfe(new string('A', 1001)));
    }

    [Test]
    public void AssessHairpinByMfe_Binding_HsaMir21_MatchesLibraryLockedValues()
    {
        // Library MiRnaAnalyzer_PreMiRna_Tests CL1: hsa-mir-21 (miRBase MI0000077) MFE -35.13 (Turner 2004), MFEI 1.0037142857142858 (Zhang 2006).
        var r = AnnotationTools.AssessHairpinByMfe(Mir21);
        Assert.That(r.Found, Is.True);
        var h = r.Hairpin!;
        var lib = MiRnaAnalyzer.AssessHairpinByMfe(Mir21)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(h.FreeEnergy, Is.EqualTo(-35.13).Within(1e-9));
            Assert.That(h.Mfei, Is.EqualTo(1.0037142857142858).Within(1e-12));
            Assert.That(h.Start, Is.EqualTo(0));
            Assert.That(h.End, Is.EqualTo(Mir21.Length - 1));
            Assert.That(h.StemBasePairs, Is.GreaterThanOrEqualTo(16)); // Ambros 2003
            Assert.That(h.TerminalLoopSize, Is.InRange(3, 25));        // Bartel 2004
            // delegation equality
            Assert.That(h.DotBracket, Is.EqualTo(lib.DotBracket));
            Assert.That(h.StemBasePairs, Is.EqualTo(lib.StemBasePairs));
            Assert.That(h.TerminalLoopStart, Is.EqualTo(lib.TerminalLoopStart));
            Assert.That(h.TerminalLoopSize, Is.EqualTo(lib.TerminalLoopSize));
            Assert.That(h.Amfe, Is.EqualTo(lib.Amfe));
            Assert.That(h.Sequence, Is.EqualTo(lib.Sequence));
        });
    }

    [Test]
    public void AssessHairpinByMfe_Rejected_ReturnsNotFound()
    {
        Assert.That(AnnotationTools.AssessHairpinByMfe(Unstructured).Found, Is.False);
        Assert.That(AnnotationTools.AssessHairpinByMfe(Unstructured).Hairpin, Is.Null);
        // MFEI 1.0037 < 1.5 -> rejected by the Zhang threshold only
        Assert.That(AnnotationTools.AssessHairpinByMfe(Mir21, 1.5).Found, Is.False);
    }

    [Test]
    public void AssessHairpinByMfe_DnaInput_EqualsRna()
    {
        var a = AnnotationTools.AssessHairpinByMfe(Mir21.Replace('U', 'T').ToLowerInvariant());
        Assert.That(a.Hairpin!.FreeEnergy, Is.EqualTo(-35.13).Within(1e-9));
    }
}
