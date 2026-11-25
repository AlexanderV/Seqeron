// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.EmbossSequenceFormats.Tools;

namespace Seqeron.Mcp.EmbossSequenceFormats.Tests;

/// <summary>supermatcher_align / needleall_align / ksw2_extend and the `format` parameter of the pairwise and multiple-alignment tools (audit A2-27b / A2-28c / A2-29b, B09 F124 - F126). The tools delegate to
/// EmbossPairwisePrograms.Supermatcher / NeedleAll, Ksw2Aligner and AlignmentFormats / EmbossAlignmentWriter; the oracle comparisons live in the library tests.</summary>
[TestFixture]
public class EmbossProgramToolTests
{

    [Test]
    public void FormatParameter_PairwiseTools_WriteTheCanonicalFormats()
    {
        var plain = AlignmentTools.GlobalAlign("ACGT", "ACGT");
        Assert.That(plain.Formatted, Is.Null);
        var fasta = AlignmentTools.GlobalAlign("ACGT", "AC-T".Replace("-", "G"), format: "fasta", names: ["x", "y"]);
        Assert.That(fasta.Formatted, Is.EqualTo(">x\nACGT\n>y\nACGT\n"));
        var clustal = AlignmentTools.GlobalAlign("ACGT", "ACGT", format: "clustal");
        Assert.That(clustal.Formatted, Does.StartWith("Biopython 1.88 multiple sequence alignment\n\n\nseq1"));
        var pair = AlignmentTools.LocalAlignAffine("ACGTACGT", "ACGTACGT", format: "pair", names: ["a", "b"]);
        Assert.That(pair.Formatted, Does.StartWith("########################################\n# Program: water\n"), "a local result is a water report (R14-1)");
        Assert.That(AlignmentTools.GlobalAlignAffine("ACGT", "ACGT", format: "pair").Formatted, Does.StartWith("########################################\n# Program: needle\n"));
        var emboss = AlignmentTools.GlobalAlign("ACGT", "ACGT", format: "emboss:fasta", names: ["a", "b"]);
        Assert.That(emboss.Formatted, Is.EqualTo(EmbossSequenceFormatsTools.WriteAlignment(["ACGT", "ACGT"], "fasta", ["a", "b"]).Text));
        Assert.Throws<ArgumentException>(() => AlignmentTools.GlobalAlign("ACGT", "ACGT", format: "bogus"));
    }

    [Test]
    public void WriteAlignment_BiopythonDialect_WritesTheBioAlignFiles()
    {
        var r = EmbossSequenceFormatsTools.WriteAlignment(["AC-T", "ACGT"], "fasta", ["s1", "s2"], descriptions: ["first", ""], dialect: "biopython");
        Assert.That(r.Dialect, Is.EqualTo("biopython"));
        Assert.That(r.Text, Is.EqualTo(">s1 first\nAC-T\n>s2\nACGT\n"));
        var nexus = EmbossSequenceFormatsTools.WriteAlignment(["AC-T", "ACGT"], "nexus", ["s1", "s2"], dialect: "biopython");
        Assert.That(nexus.Text, Is.EqualTo("#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\nformat datatype=dna missing=? gap=-;\nmatrix\ns1 AC-T\ns2 ACGT\n;\nend;\n"));
        Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.WriteAlignment(["AC-T"], "pair", dialect: "biopython"));
        Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.WriteAlignment(["AC-T", "ACGT"], "fasta", dialect: "other"));
    }}
