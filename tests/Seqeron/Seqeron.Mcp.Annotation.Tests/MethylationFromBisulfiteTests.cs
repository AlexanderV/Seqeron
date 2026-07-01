using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class MethylationFromBisulfiteTests
{
    private const string Reference = "ACGTACGT"; // CpG sites at 1 and 5

    private static List<BisulfiteReadDto> Reads() =>
        new()
        {
            new BisulfiteReadDto("ACGTACGT", 0), // C at both CpG sites (methylated)
            new BisulfiteReadDto("ACGTATGT", 0), // C at site 1, T at site 5 (converted)
        };

    [Test]
    public void MethylationFromBisulfite_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.MethylationFromBisulfite(Reference, Reads()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.MethylationFromBisulfite("", Reads()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.MethylationFromBisulfite(null!, Reads()));
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.MethylationFromBisulfite(Reference, new List<BisulfiteReadDto>()));
        Assert.Throws<ArgumentException>(
            () => AnnotationTools.MethylationFromBisulfite(Reference, null!));
    }

    [Test]
    public void MethylationFromBisulfite_Binding_InvokesSuccessfully()
    {
        // EpigeneticsAnalyzer.CalculateMethylationFromBisulfite (Bismark): at each reference CpG,
        // read C = methylated, read T = unmethylated; level = meth/(meth+unmeth), coverage = C/T calls.
        // Site 1: two C calls -> level 1.0, coverage 2. Site 5: one C + one T -> level 0.5, coverage 2.
        var result = AnnotationTools.MethylationFromBisulfite(Reference, Reads());

        Assert.That(result.Sites, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            Assert.That(result.Sites[0].Position, Is.EqualTo(1));
            Assert.That(result.Sites[0].Type, Is.EqualTo("CpG"));
            Assert.That(result.Sites[0].Context, Is.EqualTo("CGT"));
            Assert.That(result.Sites[0].MethylationLevel, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(result.Sites[0].Coverage, Is.EqualTo(2));

            Assert.That(result.Sites[1].Position, Is.EqualTo(5));
            Assert.That(result.Sites[1].MethylationLevel, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(result.Sites[1].Coverage, Is.EqualTo(2));
        });
    }

    [Test]
    public void MethylationFromBisulfite_UncoveredSites_Excluded()
    {
        // A read covering only site 1 leaves site 5 at zero coverage -> excluded.
        var reads = new List<BisulfiteReadDto> { new BisulfiteReadDto("ACG", 0) };
        var result = AnnotationTools.MethylationFromBisulfite(Reference, reads);
        Assert.That(result.Sites, Has.Count.EqualTo(1));
        Assert.That(result.Sites[0].Position, Is.EqualTo(1));
    }
}
