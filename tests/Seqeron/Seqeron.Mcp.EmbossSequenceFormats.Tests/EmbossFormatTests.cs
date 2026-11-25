// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossSequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.EmbossSequenceFormats.Tools;

namespace Seqeron.Mcp.EmbossSequenceFormats.Tests;

/// <summary>stretcher_align / write_alignment: EMBOSS 6.6 stretcher and the -aformat output formats (audit A2-27 / A2-28, B09 F120 / F121). The tools delegate to EmbossPairwisePrograms.StretcherAlign and EmbossAlignmentWriter; the oracle comparisons live in the library tests.</summary>
[TestFixture]
public class EmbossFormatTests
{

    [Test]
    public void WriteAlignment_PairFormat_IsTheNeedleLayout()
    {
        var r = EmbossSequenceFormatsTools.WriteAlignment(["ACGTACGTAGCTAGCTAGCATCG-ATCGATCGAT", B], "pair", ["a", "b"], score: 146, matrix: "EDNAFULL", gapOpen: 10, gapExtend: 0.5,
            runDate: "2026-11-19T21:43:53", commandLine: ["-asequence a.fa", "-bsequence b.fa", "-aformat3 pair", "-auto", "-stdout"]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Format, Is.EqualTo("pair"));
            Assert.That(r.Dialect, Is.EqualTo("emboss"));
            Assert.That(r.Text, Does.StartWith("########################################\n# Program: needle\n# Rundate: Thu 19 Nov 2026 21:43:53\n# Commandline: needle\n#    -asequence a.fa\n"));
            Assert.That(r.Text, Does.Contain("# Identity:      32/34 (94.1%)\n"));
            Assert.That(r.Text, Does.Contain("# Score: 146.0\n"));
            Assert.That(r.Text, Does.Contain("a                  1 ACGTACGTAGCTAGCTAGCATCG-ATCGATCGAT     33\n"));
        });
    }

    [Test]
    public void WriteAlignment_SequenceSetFormats_TakeAnyNumberOfRows()
    {
        var rows = new[] { "ACGT-ACGTA", "AC-TTACG-A", "ACGTTAC--A" };
        var fasta = EmbossSequenceFormatsTools.WriteAlignment(rows, "fasta", ["s1", "s2", "s3"], descriptions: ["desc one", "desc one", "desc one"]);
        Assert.That(fasta.Text, Is.EqualTo(">s1 desc one\nACGT-ACGTA\n>s2 desc one\nAC-TTACG-A\n>s3 desc one\nACGTTAC--A\n"));
        var clustal = EmbossSequenceFormatsTools.WriteAlignment(rows, "clustal", ["s1", "s2", "s3"]);
        Assert.That(clustal.Text, Does.StartWith("CLUSTAL W (1.83) multiple sequence alignment\n\n\ns1              ACGT-ACGTA\n"));
        var phylip = EmbossSequenceFormatsTools.WriteAlignment(rows, "phylip", ["s1", "s2", "s3"]);
        Assert.That(phylip.Text, Is.EqualTo(" 3 10\ns1        ACGT-ACGTA\ns2        AC-TTACG-A\ns3        ACGTTAC--A\n"));
    }

    [Test]
    public void WriteAlignment_Validation()
    {
        Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.WriteAlignment(["ACGT", "ACGT"], "bam"));
        Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.WriteAlignment([], "fasta"));
        Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.WriteAlignment(["ACGT", "ACGT", "ACGT"], "pair"), "pairwise formats take exactly two rows");
        Assert.Throws<ArgumentNullException>(() => EmbossSequenceFormatsTools.WriteAlignment(null!, "fasta"));
        Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.WriteAlignment(["ACGT"], "score"), "one row: score needs two");
    }

    [Test]
    public void WriteAlignment_ThreeRows_SimpleHasTheConsensusMarkup_PairFormatsRefuseLikeTheProgram()
    {
        string[] rows = ["MKVLAAGIVGLLLAQP-", "MKVLSAGIVA-LLAQPP", "MRVLAAGLVGLLLVQP-"];
        var simple = EmbossSequenceFormatsTools.WriteAlignment(rows, "simple", names: ["a", "b", "c"], program: "aligncopy", runDate: "2026-11-21T19:14:37", commandLine: ["-sequences a3.fa", "-aformat simple", "-outfile stdout", "-auto"]);
        var ex = Assert.Throws<ArgumentException>(() => EmbossSequenceFormatsTools.WriteAlignment(rows, "pair", names: ["a", "b", "c"]));
        Assert.Multiple(() =>
        {
            Assert.That(simple.Text, Does.Contain("# Aligned_sequences: 3\n# 1: a\n# 2: b\n# 3: c\n# Matrix: EBLOSUM62\n#\n# Length: 17\n# Identity:      10/17 (58.8%)\n# Similarity:    16/17 (94.1%)\n# Gaps:           2/17 (11.8%)\n"));
            Assert.That(simple.Text, Does.Contain("a                  1 MKVLAAGIVGLLLAQP-     16\nb                  1 MKVLSAGIVA-LLAQPP     16\nc                  1 MRVLAAGLVGLLLVQP-     16\n                     |:||:||:|::||:|| \n"));
            Assert.That(ex!.Message, Does.StartWith("Alignment format pair specifies at most 2 sequences, alignment has 3"));
        });
    }

    [Test]
    public void WriteAlignment_ThreeRows_TCoffeeListsEveryPairOfRows()
    {
        var r = EmbossSequenceFormatsTools.WriteAlignment(["ACGT-ACGTTAGC", "ACGTTACG-TAGC", "AC-TTACGTTAGA"], "tcoffee", names: ["s1", "s2", "s3"]);
        Assert.Multiple(() =>
        {
            Assert.That(r.Text, Does.StartWith("3\ns1 12 ACGTACGTTAGC\ns2 12 ACGTTACGTAGC\ns3 12 ACTTACGTTAGA\n"));
            Assert.That(r.Text, Does.Contain("#1 2\n1 1 0 1 0\n2 2 0 1 0\n"));
            Assert.That(r.Text, Does.Contain("#1 3\n1 1 0 1 0\n2 2 0 1 0\n4 3 0 1 0\n"));
            Assert.That(r.Text, Does.Contain("#2 3\n"));
            Assert.That(r.Text, Does.EndWith("! SEQ_1_TO_N\n"));
        });
    }

    private const string B = "ACGTACGTAGCTAGATAGCATCGAATCGATCGAT";
}
