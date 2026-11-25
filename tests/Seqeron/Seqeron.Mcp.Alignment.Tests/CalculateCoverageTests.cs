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

    [Test]
    public void CalculateCoverage_StrandAndMultiMapOptions_DelegateToCore()
    {
        // samtools depth -a oracle (Evidence): TTGACC = revcomp(GGTCAA at [2,8)).
        Assert.That(AlignmentTools.CalculateCoverage("ACGGTCAATGCC", new[] { "TTGACC" }, 6).Coverage,
            Is.EqualTo(new[] { 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 }), "A6-1: both strands by default");
        Assert.That(AlignmentTools.CalculateCoverage("ACGGTCAATGCC", new[] { "TTGACC" }, 6, includeReverseComplement: false).Coverage,
            Is.EqualTo(new int[12]), "explicit forward-only (Bowtie 2 --norc)");
        // ACGTA hits loci 0 and 4 → skipped with excludeMultiMapped (samtools depth -Q 2).
        Assert.That(AlignmentTools.CalculateCoverage("ACGTACGTAC", new[] { "ACGTA", "GTACG" }, 5,
                includeReverseComplement: false, excludeMultiMapped: true).Coverage,
            Is.EqualTo(new[] { 0, 0, 1, 1, 1, 1, 1, 0, 0, 0 }));
        // Both strands (default): GTACG also matches its reverse complement CGTAC at loci 1 and 5 -> multi-mapped too
        // (samtools depth -a -Q 2 oracle: all zeros).
        Assert.That(AlignmentTools.CalculateCoverage("ACGTACGTAC", new[] { "ACGTA", "GTACG" }, 5, excludeMultiMapped: true).Coverage,
            Is.EqualTo(new int[10]));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.CalculateCoverage("ACGT", new[] { "ACGT" }, 0));
    }
}
