using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class SimulateBisulfiteConversionTests
{
    [Test]
    public void SimulateBisulfiteConversion_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.SimulateBisulfiteConversion("ACGT"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.SimulateBisulfiteConversion(""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.SimulateBisulfiteConversion(null!));
    }

    [Test]
    public void SimulateBisulfiteConversion_FullyUnmethylated_ConvertsAllC()
    {
        // EpigeneticsAnalyzer.SimulateBisulfiteConversion (Frommer 1992): unmethylated C -> T.
        // ACGCGT: C at index 1 and 3 both convert -> ATGTGT.
        var result = AnnotationTools.SimulateBisulfiteConversion("ACGCGT");
        Assert.That(result.Converted, Is.EqualTo("ATGTGT"));
    }

    [Test]
    public void SimulateBisulfiteConversion_MethylatedPositionsProtected()
    {
        // Protecting the cytosine at index 1 keeps it a C; the C at index 3 still converts.
        var result = AnnotationTools.SimulateBisulfiteConversion("ACGCGT", new List<int> { 1 });
        Assert.That(result.Converted, Is.EqualTo("ACGTGT"));
    }

    [Test]
    public void SimulateBisulfiteConversion_PreservesCaseAndNonCytosine()
    {
        // Lower-case c -> t; non-cytosine bases pass through unchanged.
        var result = AnnotationTools.SimulateBisulfiteConversion("acgtGA");
        Assert.That(result.Converted, Is.EqualTo("atgtGA"));
    }
}
