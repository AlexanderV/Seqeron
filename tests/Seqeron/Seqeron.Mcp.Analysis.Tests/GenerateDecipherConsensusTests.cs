using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// generate_decipher_consensus — delegation to <c>MotifFinder.GenerateDecipherConsensus</c>.
/// Expected values: DECIPHER man/ConsensusSequence.Rd examples, reproduced by DECIPHER's own
/// R/C source (locked in MotifFinder_DecipherConsensus_Tests).
/// </summary>
[TestFixture]
public class GenerateDecipherConsensusTests
{
    [Test]
    public void GenerateDecipherConsensus_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.GenerateDecipherConsensus(new[] { "ACGT" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDecipherConsensus(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDecipherConsensus(Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDecipherConsensus(new[] { "A" }, "xyz"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDecipherConsensus(new[] { "A" }, noConsensusChar: "NN"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateDecipherConsensus(new[] { "A" }, noConsensusChar: "Z"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.GenerateDecipherConsensus(new[] { "A" }, threshold: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.GenerateDecipherConsensus(new[] { "A" }, minInformation: 0));
    }

    [Test]
    public void GenerateDecipherConsensus_ManualExamples_EqualDecipher()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.GenerateDecipherConsensus(new[] { "A", "A", "A", "T" }).Consensus, Is.EqualTo("W"));
            Assert.That(AnalysisTools.GenerateDecipherConsensus(new[] { "A", "A", "A", "T" }, threshold: 0.3, minInformation: 0.8, noConsensusChar: "N").Consensus, Is.EqualTo("N"));
            Assert.That(AnalysisTools.GenerateDecipherConsensus(new[] { "ANGCT-", "-ACCT-" }, includeTerminalGaps: true).Consensus, Is.EqualTo("+NSCT-"));
            Assert.That(AnalysisTools.GenerateDecipherConsensus(new[] { "ANQIH-", "ADELW." }, "protein").Consensus, Is.EqualTo("ABZJX-"));
            Assert.That(AnalysisTools.GenerateDecipherConsensus(new[] { "A-+.A", "AAAAA" }, includeNonLetters: true).Consensus, Is.EqualTo("AAAAA"));
            Assert.That(AnalysisTools.GenerateDecipherConsensus(new[] { "AWNDA", "AAAAA" }, ambiguity: false).Consensus, Is.EqualTo("AAAAA"));
            Assert.That(AnalysisTools.GenerateDecipherConsensus(new[] { "ACGU", "ACGU", "UCGA" }, "RNA").Consensus, Is.EqualTo("WCGW"));
        });
    }

    [Test]
    public void GenerateDecipherConsensus_EqualsLibrary()
    {
        string[] rows = { "GTT", "GAA", "CTG" };
        Assert.That(AnalysisTools.GenerateDecipherConsensus(rows, threshold: 0.5).Consensus,
            Is.EqualTo(global::Seqeron.Genomics.Analysis.MotifFinder.GenerateDecipherConsensus(rows, threshold: 0.5)));
    }
}
