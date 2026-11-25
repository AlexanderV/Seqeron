using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Models;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for <c>estimate_completeness_from_kmers</c> (GenomeAssemblyAnalyzer.EstimateCompletenessFromKmers):
/// Merqury solid threshold, GenomeScope error rate / genome size, and with an assembly Merqury k-mer
/// completeness and QV. Values from tools/review/oracles/B10/completeness_oracle.py (ASSEMBLY-COMPLETE-001).
/// </summary>
[TestFixture]
public class EstimateCompletenessFromKmersTests
{
    // 20 canonical 3-mers, histogram 1:5, 2:2, 3:1, 4:3, 5:6, 6:3.
    private static List<KmerCount> ValleySpectrum() => new()
    {
        new("AAA", 1), new("AAC", 1), new("AAG", 1), new("AAT", 1), new("ACA", 1), new("ACC", 2), new("ACG", 2),
        new("ACT", 3), new("AGA", 4), new("AGC", 4), new("AGG", 4), new("ATA", 5), new("ATC", 5), new("ATG", 5),
        new("CAA", 5), new("CAC", 5), new("CAG", 5), new("CCA", 6), new("CCC", 6), new("CCG", 6),
    };

    [Test]
    public void EstimateCompletenessFromKmers_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.EstimateCompletenessFromKmers(ValleySpectrum(), 30));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateCompletenessFromKmers(null!, 30));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateCompletenessFromKmers(ValleySpectrum(), expectedCoverage: -1));
        // Mixed k-mer lengths and non-canonical keys (with an assembly) are rejected by the library.
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateCompletenessFromKmers(
            new List<KmerCount> { new("AAA", 2), new("AAAA", 2) }));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.EstimateCompletenessFromKmers(
            new List<KmerCount> { new("TTT", 2) }, assembly: new List<NamedSequence> { new("a", "AAAA") }));
    }

    [Test]
    public void EstimateCompletenessFromKmers_SpectrumOnly_ErrorRateAndSize_NoCompleteness()
    {
        // filt 3, E = 12 of T = 72 ⇒ error 1 − (60/72)^(1/3); peak 5 ⇒ size 12; no assembly ⇒ null completeness.
        var r = ChromosomeTools.EstimateCompletenessFromKmers(ValleySpectrum());
        Assert.Multiple(() =>
        {
            Assert.That(r.Completeness, Is.Null);
            Assert.That(r.ErrorRate, Is.EqualTo(0.05896397111897156).Within(1e-12));
            Assert.That(r.EstimatedGenomeSize, Is.EqualTo(12));
            Assert.That((r.K, r.SolidThreshold, r.PeakCoverage, r.TotalKmers, r.ErrorKmers, r.SolidKmers),
                Is.EqualTo((3, 3, 5, 72L, 12L, 12L)));
            Assert.That(r.QualityValue, Is.Null);
            Assert.That(r.AssemblyKmers, Is.Null);
        });
    }

    [Test]
    public void EstimateCompletenessFromKmers_WithAssembly_MerquryCompletenessAndQv()
    {
        var asm = new List<NamedSequence> { new("a", "ACGTTGCAAC"), new("b", "GGNNCATTAC") };
        var r = ChromosomeTools.EstimateCompletenessFromKmers(ValleySpectrum(), assembly: asm);
        Assert.Multiple(() =>
        {
            Assert.That(r.Completeness, Is.EqualTo(2.0 / 12).Within(1e-12));
            Assert.That((r.AssemblySolidKmers, r.AssemblyKmers, r.AssemblyOnlyKmers), Is.EqualTo(((long?)2, (long?)12, (long?)4)));
            Assert.That(r.AssemblyErrorRate, Is.EqualTo(0.1264195352637011).Within(1e-12));
            Assert.That(r.QualityValue, Is.EqualTo(8.981858105339732).Within(1e-9));
        });

        // Expected coverage overrides the peak: (72 − 12)/4 = 15.
        Assert.That(ChromosomeTools.EstimateCompletenessFromKmers(ValleySpectrum(), 4).EstimatedGenomeSize, Is.EqualTo(15));
    }
}
