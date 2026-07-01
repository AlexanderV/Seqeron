using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindSimilarMiRnasTests
{
    private static MiRnaDto Make(string name, string seq) =>
        AnnotationTools.CreateMiRna(name, seq).MiRna;

    // Query seed GAGGUAG. Database seeds: GAGGUAG (0), GAGGUCG (1), GAGGCCG (2).
    private static MiRnaDto Query() => Make("query", "UGAGGUAGAAAAAAAAAAAAAA");

    private static List<MiRnaDto> Database() =>
        new()
        {
            Make("same-seed",  "UGAGGUAGCCCCCCCCCCCCCC"), // seed GAGGUAG, 0 mismatches
            Make("one-off",    "UGAGGUCGAAAAAAAAAAAAAA"), // seed GAGGUCG, 1 mismatch
            Make("two-off",    "UGAGGCCGAAAAAAAAAAAAAA"), // seed GAGGCCG, 2 mismatches
            Make("query",      "UGAGGUAGAAAAAAAAAAAAAA"), // same name -> excluded
        };

    [Test]
    public void FindSimilarMiRnas_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.FindSimilarMiRnas(Query(), Database()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FindSimilarMiRnas(null!, Database()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.FindSimilarMiRnas(Query(), null!));
    }

    [Test]
    public void FindSimilarMiRnas_Binding_InvokesSuccessfully()
    {
        // MiRnaAnalyzer.FindSimilarMiRnas: exclude same Name; include seed Hamming distance <= maxMismatches.
        // Default maxMismatches=1: same-seed (0) and one-off (1) included; two-off (2) and self excluded.
        var result = AnnotationTools.FindSimilarMiRnas(Query(), Database());

        Assert.That(result.Matches.Select(m => m.Name),
            Is.EquivalentTo(new[] { "same-seed", "one-off" }));
    }

    [Test]
    public void FindSimilarMiRnas_MaxMismatchesTwo_IncludesTwoOff()
    {
        var result = AnnotationTools.FindSimilarMiRnas(Query(), Database(), maxMismatches: 2);
        Assert.That(result.Matches.Select(m => m.Name),
            Is.EquivalentTo(new[] { "same-seed", "one-off", "two-off" }));
    }

    [Test]
    public void FindSimilarMiRnas_ExactMatchOnly()
    {
        var result = AnnotationTools.FindSimilarMiRnas(Query(), Database(), maxMismatches: 0);
        Assert.That(result.Matches.Select(m => m.Name), Is.EquivalentTo(new[] { "same-seed" }));
    }
}
