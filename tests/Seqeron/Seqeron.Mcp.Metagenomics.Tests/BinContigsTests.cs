using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.BinContigs (deterministic GC-sorted k-means init).
// Reference values from Seqeron.Genomics.Tests MetagenomicsAnalyzer_GenomeBinning_Tests
// (CheckM completeness/contamination; Teeling 2004 TETRA; Parks 2014).
[TestFixture]
public class BinContigsTests
{
    // Repeating "GCTA" -> ~50% GC content.
    private static string MidGc(int length)
    {
        const string pattern = "GCTA";
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(pattern[i % 4]);
        return sb.ToString();
    }

    private static ContigInput[] MidGcContigs(int count, int len, double coverage)
        => Enumerable.Range(0, count)
            .Select(i => new ContigInput($"contig_{i}", MidGc(len), coverage))
            .ToArray();

    [Test]
    public void BinContigs_Schema_ValidatesCorrectly()
    {
        // Empty input is defined (no bins), not an error.
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.BinContigs(System.Array.Empty<ContigInput>()));

        // A well-formed call does not throw.
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.BinContigs(MidGcContigs(10, 200000, 15.0), numBins: 1, minBinSize: 100000));
    }

    [Test]
    public void BinContigs_Binding_InvokesSuccessfully()
    {
        // Below minBinSize -> no bin reported.
        var tooSmall = MetagenomicsTools.BinContigs(
            new[] { new ContigInput("c1", MidGc(1000), 10.0) });
        Assert.That(tooSmall.Items, Is.Empty,
            "1000 bp < 500000 default minBinSize -> no bins.");

        // 10 x 200kb mid-GC contigs, numBins=1 -> a single bin.
        // TotalLength 2,000,000; Completeness = 2M/4M * 100 = 50.0; uniform GC -> Contamination 0;
        // GcContent = 0.5.
        var result = MetagenomicsTools.BinContigs(
            MidGcContigs(10, 200000, 15.0), numBins: 1, minBinSize: 100000);

        Assert.Multiple(() =>
        {
            Assert.That(result.Items, Has.Count.EqualTo(1));
            var bin = result.Items[0];
            Assert.That(bin.TotalLength, Is.EqualTo(2_000_000.0));
            Assert.That(bin.Completeness, Is.EqualTo(50.0),
                "Completeness = min(2M/4M * 100, 100) = 50.0 (CheckM/Parks 2014).");
            Assert.That(bin.Contamination, Is.EqualTo(0.0),
                "Uniform within-bin GC -> zero GC stddev -> contamination 0.");
            Assert.That(bin.GcContent, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(bin.Coverage, Is.EqualTo(15.0).Within(1e-9));
            Assert.That(bin.ContigIds, Has.Count.EqualTo(10));
        });
    }
}
