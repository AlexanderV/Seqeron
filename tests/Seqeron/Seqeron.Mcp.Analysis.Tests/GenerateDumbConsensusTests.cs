using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// generate_dumb_consensus — delegation to <c>MotifFinder.GenerateDumbConsensus</c>.
/// Expected values: Biopython 1.85 dumb_consensus outputs locked in MotifFinder_AlignmentConsensus_Tests.
/// </summary>
[TestFixture]
public class GenerateDumbConsensusTests
{
    [Test]
    public void GenerateDumbConsensus_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.GenerateDumbConsensus(new[] { "ACGT" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDumbConsensus(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDumbConsensus(Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDumbConsensus(new[] { "ACGT" }, ambiguous: "XY"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.GenerateDumbConsensus(new[] { "ACGT" }, double.NaN));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDumbConsensus(new[] { "ACG", "AC" }));
    }

    [Test]
    public void GenerateDumbConsensus_EqualsBiopython()
    {
        // Dumb_000: threshold 0.7, ambiguous 'N'
        Assert.That(AnalysisTools.GenerateDumbConsensus(new[] { "ACGT", "ATGT", "ATGT" }, 0.7, "N").Consensus, Is.EqualTo("ANGT"));
        // Dumb_002: require_multiple
        Assert.That(AnalysisTools.GenerateDumbConsensus(new[] { "A-", "-A", "--" }, 0.7, "X", true).Consensus, Is.EqualTo("XX"));
        // Dumb_003: protein, threshold 0.5
        Assert.That(AnalysisTools.GenerateDumbConsensus(new[] { "MKV", "MKV", "mRV" }, 0.5).Consensus, Is.EqualTo("MKV"));
    }

    [Test]
    public void GenerateDumbConsensus_Defaults_EqualLibraryDefaults()
    {
        var rows = new[] { "ACGT", "AT-T", "CT-T", "GT-T" };
        Assert.That(AnalysisTools.GenerateDumbConsensus(rows).Consensus,
            Is.EqualTo(global::Seqeron.Genomics.Analysis.MotifFinder.GenerateDumbConsensus(rows)));
    }
}
