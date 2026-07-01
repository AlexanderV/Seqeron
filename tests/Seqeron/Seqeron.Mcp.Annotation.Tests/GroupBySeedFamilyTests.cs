using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class GroupBySeedFamilyTests
{
    private static MiRnaDto Make(string name, string seq) =>
        AnnotationTools.CreateMiRna(name, seq).MiRna;

    // Seeds: let-7a/let-7c share GAGGUAG; other has GAGGUCG.
    private static List<MiRnaDto> MiRnas() =>
        new()
        {
            Make("let-7a", "UGAGGUAGUAGGUUGUAUAGUU"),  // seed GAGGUAG
            Make("let-7c", "UGAGGUAGUAGGUUGUGUGGUU"),  // seed GAGGUAG (same family)
            Make("other",  "UGAGGUCGUAGGUUGUAUAGUU"),  // seed GAGGUCG (different)
        };

    [Test]
    public void GroupBySeedFamily_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.GroupBySeedFamily(MiRnas()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.GroupBySeedFamily(null!));
    }

    [Test]
    public void GroupBySeedFamily_Binding_InvokesSuccessfully()
    {
        // MiRnaAnalyzer.GroupBySeedFamily groups by identical seed sequence.
        var result = AnnotationTools.GroupBySeedFamily(MiRnas());

        Assert.That(result.Families, Has.Count.EqualTo(2));
        Assert.Multiple(() =>
        {
            var gaggUag = result.Families.Single(f => f.SeedFamily == "GAGGUAG");
            Assert.That(gaggUag.Members.Select(m => m.Name),
                Is.EquivalentTo(new[] { "let-7a", "let-7c" }));

            var gaggUcg = result.Families.Single(f => f.SeedFamily == "GAGGUCG");
            Assert.That(gaggUcg.Members.Select(m => m.Name), Is.EquivalentTo(new[] { "other" }));
        });
    }

    [Test]
    public void GroupBySeedFamily_EmptyInput_ReturnsEmpty()
    {
        var result = AnnotationTools.GroupBySeedFamily(new List<MiRnaDto>());
        Assert.That(result.Families, Is.Empty);
    }
}
