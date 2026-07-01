using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

[TestFixture]
public class CalculateCoverageTests
{
    [Test]
    public void CalculateCoverage_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.CalculateCoverage("ACGTACGTAC", new[] { "ACGTA" }, 5));
        // Coverage array length always equals reference length.
        var r = AlignmentTools.CalculateCoverage("ACGTACGTAC", Array.Empty<string>(), 5);
        Assert.That(r.Coverage, Has.Length.EqualTo(10));
        Assert.That(r.Coverage, Is.All.EqualTo(0));
    }

    [Test]
    public void CalculateCoverage_Binding_InvokesSuccessfully()
    {
        // Reference "ACGTACGTAC" (len 10). Read "ACGTA" best-places at 0 (covers [0,5));
        // read "GTACG" best-places at 2 (covers [2,7)). Per-base depth follows.
        var r = AlignmentTools.CalculateCoverage("ACGTACGTAC", new[] { "ACGTA", "GTACG" }, 5);
        Assert.That(r.Coverage, Is.EqualTo(new[] { 1, 1, 2, 2, 2, 1, 1, 0, 0, 0 }));
    }
}
