using System.Collections.Generic;
using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for the <c>analyze_karyotype</c> MCP tool.
/// Expected values come from ChromosomeAnalyzer.AnalyzeKaryotype documented behavior
/// (Seqeron.Genomics.Tests/Unit/Chromosome/ChromosomeAnalyzer_Karyotype_Tests.cs, CHROM-KARYO-001) —
/// NOT from the wrapper's own output.
/// </summary>
[TestFixture]
public class AnalyzeKaryotypeTests
{
    private static ChromosomeInput C(string name, long len, bool sex) => new(name, len, sex);

    [Test]
    public void AnalyzeKaryotype_Schema_ValidatesCorrectly()
    {
        var valid = new List<ChromosomeInput> { C("chr1_1", 100, false), C("chr1_2", 100, false) };

        Assert.DoesNotThrow(() => ChromosomeTools.AnalyzeKaryotype(valid));
        // Empty list is documented-valid (returns empty karyotype), not an error.
        Assert.DoesNotThrow(() => ChromosomeTools.AnalyzeKaryotype(new List<ChromosomeInput>()));

        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeKaryotype(null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeKaryotype(valid, expectedPloidyLevel: 0));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeKaryotype(valid, expectedPloidyLevel: -1));
    }

    [Test]
    public void AnalyzeKaryotype_Binding_InvokesSuccessfully()
    {
        // Normal diploid set: 2 autosome pairs + XY sex chromosomes (CHROM-KARYO-001).
        var chromosomes = new List<ChromosomeInput>
        {
            C("chr1_1", 248956422, false),
            C("chr1_2", 248956422, false),
            C("chr2_1", 242193529, false),
            C("chr2_2", 242193529, false),
            C("chrX", 156040895, true),
            C("chrY", 57227415, true),
        };

        var result = ChromosomeTools.AnalyzeKaryotype(chromosomes);

        Assert.Multiple(() =>
        {
            Assert.That(result.TotalChromosomes, Is.EqualTo(6));
            Assert.That(result.AutosomeCount, Is.EqualTo(4));
            Assert.That(result.SexChromosomes, Is.EquivalentTo(new[] { "chrX", "chrY" }));
            Assert.That(result.PloidyLevel, Is.EqualTo(2));
            Assert.That(result.HasAneuploidy, Is.False);
            Assert.That(result.Abnormalities, Is.Empty);
            Assert.That(result.TotalGenomeSize, Is.EqualTo(248956422L * 2 + 242193529L * 2 + 156040895L + 57227415L));
        });
    }

    [Test]
    public void AnalyzeKaryotype_Trisomy_IsDetected()
    {
        // Trisomy 21: three copies where diploid (2) is expected.
        var chromosomes = new List<ChromosomeInput>
        {
            C("chr21_1", 46709983, false),
            C("chr21_2", 46709983, false),
            C("chr21_3", 46709983, false),
        };

        var result = ChromosomeTools.AnalyzeKaryotype(chromosomes);

        Assert.Multiple(() =>
        {
            Assert.That(result.HasAneuploidy, Is.True);
            Assert.That(result.AutosomeCount, Is.EqualTo(3));
            Assert.That(result.Abnormalities, Is.EqualTo(new[] { "Trisomy chr21" }));
        });
    }
}
