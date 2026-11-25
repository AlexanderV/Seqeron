// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossMotifs, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.EmbossMotifs.Tests;

/// <summary>
/// generate_emboss_consensus — delegation to <c>EmbossConsensus.GenerateEmbossConsensus</c>.
/// Expected values: EMBOSS 6.6.0 cons outputs locked in MotifFinder_AlignmentConsensus_Tests (Emboss_000/003/010/013).
/// </summary>
[TestFixture]
public class GenerateEmbossConsensusTests
{
    private static readonly string[] Dna = { "ACGTAC-T", "ACGTTCAT", "AGGTAC-T", "tCGAAG-T" };
    private static readonly string[] Protein = { "MKVLAAGIVG", "MKVLSAGIVA", "MRVLAAG-VG", "MKILTAGLVG" };

    [Test]
    public void GenerateEmbossConsensus_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => EmbossMotifsTools.GenerateEmbossConsensus(Dna));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.GenerateEmbossConsensus(null!));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.GenerateEmbossConsensus(new[] { "ACGT" }));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.GenerateEmbossConsensus(Dna, "rna"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.GenerateEmbossConsensus(new[] { "ACG", "AC" }));
        Assert.Throws<ArgumentOutOfRangeException>(() => EmbossMotifsTools.GenerateEmbossConsensus(Dna, identity: -1));
    }

    [Test]
    public void GenerateEmbossConsensus_Defaults_EqualsCons()
    {
        // cons -sequence a.fa -snucleotide
        Assert.That(EmbossMotifsTools.GenerateEmbossConsensus(Dna).Consensus, Is.EqualTo("ACGTACnT"));
        // cons -sequence a.fa -sprotein
        Assert.That(EmbossMotifsTools.GenerateEmbossConsensus(Protein, "protein").Consensus, Is.EqualTo("MKVLAAGIVG"));
    }

    [Test]
    public void GenerateEmbossConsensus_PluralitySetcase_EqualsCons_AndEqualsLibrary()
    {
        // cons -plurality 1.0 -setcase 3.0
        Assert.That(EmbossMotifsTools.GenerateEmbossConsensus(Dna, plurality: 1.0, setcase: 3.0).Consensus, Is.EqualTo("acGtacaT"));
        Assert.That(EmbossMotifsTools.GenerateEmbossConsensus(Protein, "Protein", 1.0, 0, 3.0).Consensus, Is.EqualTo("MKVLaAGiVg"));
        Assert.That(EmbossMotifsTools.GenerateEmbossConsensus(Dna, identity: 4).Consensus,
            Is.EqualTo(EmbossConsensus.GenerateEmbossConsensus(Dna, identity: 4)));
    }

    [Test]
    public void GenerateEmbossConsensus_AutoAndPadRaggedRows_EqualCons()
    {
        // cons -sequence a.fa (no -snucleotide/-sprotein): ragged rows padded by ajSeqsetFill.
        Assert.That(EmbossMotifsTools.GenerateEmbossConsensus(new[] { "ACGTAC", "ACG", "AC" }, "auto", padRaggedRows: true).Consensus, Is.EqualTo("ACGnnn"));
        // cons -sequence a.fa -setcase 2.15: N/X decided by the first sequence.
        Assert.That(EmbossMotifsTools.GenerateEmbossConsensus(new[] { "V", "x", "~" }, "auto", setcase: 2.15).Consensus, Is.EqualTo("n"));
        Assert.That(EmbossMotifsTools.GenerateEmbossConsensus(Protein, "auto").Consensus, Is.EqualTo("MKVLAAGIVG"));
        Assert.Throws<ArgumentException>(() => EmbossMotifsTools.GenerateEmbossConsensus(new[] { "ACG", "AC" }, "auto"));
    }
}
