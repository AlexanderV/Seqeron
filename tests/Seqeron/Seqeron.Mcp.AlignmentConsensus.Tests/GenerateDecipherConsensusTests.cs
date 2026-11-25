// SPDX-License-Identifier: GPL-3.0-only
// Wrapper of the DECIPHER ConsensusSequence port (GPL-3.0-only).
// Part of Seqeron.Mcp.AlignmentConsensus, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Mcp.AlignmentConsensus.Tools;
using Seqeron.Genomics.AlignmentConsensus;
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.AlignmentConsensus.Tests;

/// <summary>
/// generate_decipher_consensus — delegation to <c>DecipherConsensus.GenerateDecipherConsensus</c>.
/// Expected values: DECIPHER man/ConsensusSequence.Rd examples, reproduced by DECIPHER's own
/// R/C source (locked in MotifFinder_DecipherConsensus_Tests).
/// </summary>
[TestFixture]
public class GenerateDecipherConsensusTests
{
    [Test]
    public void GenerateDecipherConsensus_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "ACGT" }));
        Assert.Throws<ArgumentException>(() => AlignmentConsensusTools.GenerateDecipherConsensus(null!));
        Assert.Throws<ArgumentException>(() => AlignmentConsensusTools.GenerateDecipherConsensus(Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "A" }, "xyz"));
        Assert.Throws<ArgumentException>(() => AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "A" }, noConsensusChar: "NN"));
        Assert.Throws<ArgumentException>(() => AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "A" }, noConsensusChar: "Z"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "A" }, threshold: 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "A" }, minInformation: 0));
    }

    [Test]
    public void GenerateDecipherConsensus_ManualExamples_EqualDecipher()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "A", "A", "A", "T" }).Consensus, Is.EqualTo("W"));
            Assert.That(AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "A", "A", "A", "T" }, threshold: 0.3, minInformation: 0.8, noConsensusChar: "N").Consensus, Is.EqualTo("N"));
            Assert.That(AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "ANGCT-", "-ACCT-" }, includeTerminalGaps: true).Consensus, Is.EqualTo("+NSCT-"));
            Assert.That(AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "ANQIH-", "ADELW." }, "protein").Consensus, Is.EqualTo("ABZJX-"));
            Assert.That(AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "A-+.A", "AAAAA" }, includeNonLetters: true).Consensus, Is.EqualTo("AAAAA"));
            Assert.That(AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "AWNDA", "AAAAA" }, ambiguity: false).Consensus, Is.EqualTo("AAAAA"));
            Assert.That(AlignmentConsensusTools.GenerateDecipherConsensus(new[] { "ACGU", "ACGU", "UCGA" }, "RNA").Consensus, Is.EqualTo("WCGW"));
        });
    }

    [Test]
    public void GenerateDecipherConsensus_EqualsLibrary()
    {
        string[] rows = { "GTT", "GAA", "CTG" };
        Assert.That(AlignmentConsensusTools.GenerateDecipherConsensus(rows, threshold: 0.5).Consensus,
            Is.EqualTo(global::Seqeron.Genomics.AlignmentConsensus.DecipherConsensus.GenerateDecipherConsensus(rows, threshold: 0.5)));
    }
}
