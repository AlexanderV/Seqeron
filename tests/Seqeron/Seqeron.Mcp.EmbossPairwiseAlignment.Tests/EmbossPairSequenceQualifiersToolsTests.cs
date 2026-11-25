// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>B09 A5-7: the <c>sequenceQualifiers</c> parameter of the EMBOSS pairwise tools. Values are what the installed EMBOSS 6.6.0 prints (also locked by the 180-row oracles of the library).</summary>
[TestFixture]
public class EmbossPairSequenceQualifiersToolsTests
{
    private const string A = "ACGTACGTTAGCCGATTAGGCATCGATCGGATCG";
    private const string B = "CGGATCGATCAGGACTATCGGCTAACGTACGTTAGC";

    [Test]
    public void NeedleAlign_RangeAndReverse_ReportsTheOriginalNumbering()
    {
        // needle a.fa b.fa -sbegin1 3 -send1 20 -sreverse2 -slower1 -aformat3 pair
        var r = EmbossPairwiseAlignmentTools.NeedleAlign(A, B, sequenceQualifiers: "-sbegin1 3 -send1 20 -sreverse2 -slower1");
        Assert.That(r.Formatted, Does.Contain("a                  3 ------gtacgttagccgattagg-------------     20"));
        Assert.That(r.Formatted, Does.Contain("b                 36 GCTAACGTACGTTAGCCGA-TAGTCCTGATCGATCCG      1"));
        Assert.That(r.Score, Is.EqualTo(66));
    }

    [Test]
    public void NeedleAlign_UnnumberedQualifiers_BelongToTheSecondSequence_LikeTheProgramAfterItsFiles()
    {
        // A6-2: needle a.fa b.fa -sbegin 3 -send 8 -aformat3 pair (installed 6.6.0, Commandline echo "-sbegin2 3 -send2 8"): only b is cut
        var r = EmbossPairwiseAlignmentTools.NeedleAlign(A, B, sequenceQualifiers: "-sbegin 3 -send 8");
        Assert.That(r.Formatted, Does.Contain("a                  1 ACGTACGTTAGCCGATTAGGCATCGATCGGATCG-     34\n"));
        Assert.That(r.Formatted, Does.Contain("b                  3 -----------------------------GATCGA      8\n"));
        Assert.That(r.Score, Is.EqualTo(25));
        // the same with the numbers, the = form and unique abbreviations; -sbe (no number) is -sbegin of EVERY sequence (the program resolves it without the tie to the second file)
        foreach (string text in new[] { "-sbegin2 3 -send2 8", "-sbegin2=3 --SEND2=8", "-sbeg2 3 -sen2 8" })
            Assert.That(EmbossPairwiseAlignmentTools.NeedleAlign(A, B, sequenceQualifiers: text).Formatted, Is.EqualTo(r.Formatted), text);
        var both = EmbossPairwiseAlignmentTools.NeedleAlign(A, B, sequenceQualifiers: "-sbe 3");
        Assert.That(both.Formatted, Does.Contain("a                  3 ").And.Contain("b                  3 "));
        Assert.That(Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.NeedleAlign(A, B, sequenceQualifiers: "-s 3"))!.Message, Does.Contain("ambiguous qualifier"));
    }

    [Test]
    public void NeedleAlign_WithoutQualifiers_IsUnchanged()
    {
        var plain = EmbossPairwiseAlignmentTools.NeedleAlign(A, B);
        var empty = EmbossPairwiseAlignmentTools.NeedleAlign(A, B, sequenceQualifiers: "");
        Assert.That(empty.Score, Is.EqualTo(plain.Score));
        Assert.That(empty.AlignedSequence1, Is.EqualTo(plain.AlignedSequence1));
    }

    [Test]
    public void ForcedTypeContradictingTheInput_IsRefused()
    {
        // needle: the second sequence of a nucleotide first sequence must be nucleotide ("Sequence is not nucleic")
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.NeedleAlign("ACGTACGT", "ACDEFGAIK", sequenceQualifiers: "-sbegin1 1"));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.WaterAlign(A, B, sequenceQualifiers: "-gapopen 3"));
    }

    [Test]
    public void Dotmatcher_Qualifiers_ScansTheWholeLengthOverTheCutText()
    {
        // dotmatcher da.fa db.fa -windowsize 5 -threshold 15 -matrixfile EDNAFULL -sbegin1 10: the first segment of the plain run (32 32 40 40) stays: -sbegin shifts the range to the front
        var r = EmbossPairwiseAlignmentTools.Dotmatcher("ACGTTGCATGCATTAGCGATCGATCGGATTACAGGCATCG", "TTGCATGCATTAGCGATCGATAAGGATTACAGGGGCATCG", 5, 15, "EDNAFULL", sequenceQualifiers: "-sbegin1 10");
        // the whole plot of the program (dotmatcher ... -graph data, installed 6.6.0; the program's x / y are the second / first sequence): 22 segments, the main-diagonal one twice
        var expected = new[] { (32, 32, 40, 40), (16, 6, 23, 13), (26, 14, 32, 20), (20, 1, 24, 5), (22, 3, 30, 11), (35, 14, 40, 19), (24, 1, 30, 7), (32, 4, 39, 11), (33, 1, 39, 7), (32, 32, 40, 40), (17, 18, 22, 23), (1, 3, 5, 7),
            (9, 11, 18, 20), (1, 7, 16, 22), (15, 22, 27, 34), (23, 32, 31, 40), (5, 15, 12, 22), (2, 19, 7, 24), (1, 24, 7, 30), (12, 35, 17, 40), (8, 35, 13, 40), (1, 35, 5, 39) };
        Assert.That(r.Segments.Select(s => (s.Start1, s.Start2, s.End1, s.End2)).ToArray(), Is.EqualTo(expected));
    }

    [Test]
    public void Merger_RangesAreIgnored_AsTheProgramDoes()
    {
        // merger a.fa b.fa -send1 -39 -send2 34 ... : merger does not trim, the whole sequences are merged (installed 6.6.0, oracle row of emboss_seqopts_oracle.tsv)
        const string a = "TCCGCCAGGTACACCCCTTGTGGGAAGAGTGCCTAGACC", b = "TTGTGGGAAGAGTGCCTAGACCGCCCATCTTGCG";
        var plain = EmbossPairwiseAlignmentTools.Merger(a, b);
        var ranged = EmbossPairwiseAlignmentTools.Merger(a, b, sequenceQualifiers: "-send1 -39 -send2 34");
        Assert.That(ranged.MergedSequence, Is.EqualTo("tccgccaggtacaccccttgtgggaagagtgcctagaccgcccatcttgcg"));
        Assert.That(ranged.MergedSequence, Is.EqualTo(plain.MergedSequence));
        Assert.That(ranged.Conflicts, Is.Empty);
    }

    [Test]
    public void Merger_ReverseAndCase_AreApplied()
    {
        // merger -sreverse1 -slower1 -sbegin2 52 -send2 24 (installed 6.6.0): the reverse complement of a is merged; the first conflict is a 24 'G' / b 2 'c' using 'G'
        var r = EmbossPairwiseAlignmentTools.Merger("TCCTTTCCAAACCCCCGCGCTCTTGTTAAGGTCAAATG", "GCGCTCTTGTTAAGGTCAAATGCTACTGATTGACAG", sequenceQualifiers: "-sreverse1 -slower1 -sbegin2 52 -send2 24");
        Assert.That(r.MergedSequence, Is.EqualTo("catttgaccttaacaagagcgcgGgGGTttgTTaaggTcaaatgctactgattgacag"));
        Assert.That(r.Conflicts[0], Is.EqualTo(new MergerConflictDto(24, "G", 2, "c", "G")));
        Assert.Throws<ArgumentException>(() => EmbossPairwiseAlignmentTools.Merger("ACGT", "ACGT", sequenceQualifiers: "-sbegin1"));
    }
}
