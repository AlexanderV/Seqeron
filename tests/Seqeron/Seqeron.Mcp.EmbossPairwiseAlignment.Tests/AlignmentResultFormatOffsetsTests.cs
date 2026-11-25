// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>
/// R14-1 (B09 F288, WP80): the EMBOSS alignment formats take the residues before and after each row from the ROWS, not from an end coordinate whose convention differed per producer. The core aligners (global_align /
/// local_align / semi_global_align, MIT) return the rows; <c>write_alignment</c> of the EMBOSS server writes them. Truth: the installed EMBOSS 6.6.0 `water` (same alignment ACGT / ACGT of a = TTACGTTT, b = ACGT)
/// and the fitting rows of the query inside the reference. water -aformat3 markx10 prints ">>>a, 8 nt vs b, 4 nt", al_start 3 al_stop 6 / 1 4; pair "a  3 ACGT  6" / "b  1 ACGT  4"; sam "b 0 a 3 255 4M".
/// </summary>
[TestFixture]
public class AlignmentResultFormatOffsetsTests
{
    private const string A = "TTACGTTT";
    private const string B = "ACGT";
    private static readonly string[] Names = ["a", "b"];

    // the alignment of a core aligner written as an EMBOSS format: the offsets come from the rows (EmbossAlignedSequence.Flanks); water for a local alignment, needle otherwise
    private static string Emboss(AlignmentResultDto d, string format, int length1, int length2, string[] names)
    {
        var (offset1, offend1) = EmbossAlignedSequence.Flanks(d.AlignedSequence1, d.StartPosition1, length1);
        var (offset2, offend2) = EmbossAlignedSequence.Flanks(d.AlignedSequence2, d.StartPosition2, length2);
        string program = string.Equals(d.AlignmentType, "Local", StringComparison.OrdinalIgnoreCase) ? "water" : "needle";
        return EmbossSequenceFormatsTools.WriteAlignment([d.AlignedSequence1, d.AlignedSequence2], format, names, null, [offset1, offset2], [offend1, offend2], d.ExactScore != 0 ? d.ExactScore : d.Score, program: program).Text;
    }

    [Test]
    public void LocalAlign_Markx10_HasTheLengthsAndCoordinatesOfWater()
    {
        string text = Emboss(AlignmentTools.LocalAlign(A, B), "markx10", A.Length, B.Length, Names);
        Assert.That(text, Does.Contain(">>>a, 8 nt vs b, 4 nt"));
        Assert.That(text, Does.Contain(">a ..\n; sq_len: 8\n; sq_type: D\n; al_start: 3\n; al_stop: 6\n; al_display_start: 3\nACGT\n"));
        Assert.That(text, Does.Contain(">b ..\n; sq_len: 4\n; sq_type: D\n; al_start: 1\n; al_stop: 4\n; al_display_start: 1\nACGT\n"));
    }

    [Test]
    public void LocalAlign_Pair_HasTheCoordinatesOfWater()
    {
        string text = Emboss(AlignmentTools.LocalAlign(A, B), "pair", A.Length, B.Length, Names);
        Assert.That(text, Does.Contain("a                  3 ACGT      6\n"));
        Assert.That(text, Does.Contain("b                  1 ACGT      4\n"));
    }

    [Test]
    public void LocalAlign_Sam_HasTheStartOfWater()
    {
        string text = Emboss(AlignmentTools.LocalAlign(A, B), "sam", A.Length, B.Length, Names);
        Assert.That(text, Does.Contain("b\t0\ta\t3\t255\t4M\t"));
    }

    [Test]
    public void LocalAlign_NoPositiveRegion_HasNoNegativeOffset()
    {
        string text = Emboss(AlignmentTools.LocalAlign("AAAA", "CCCC"), "markx10", 4, 4, Names);
        Assert.That(text, Does.Contain(">>>a, 4 nt vs b, 4 nt"));
        Assert.That(text, Does.Not.Contain("-1"));
    }

    [Test]
    public void SemiGlobalAlign_Pair_ListsTheWholeReferenceFromOne()
    {
        // fitting: the rows carry the whole reference (free flanks as gaps in the query), so the reference runs 1 .. 8, not 3 .. 10
        string text = Emboss(AlignmentTools.SemiGlobalAlign(B, A), "pair", B.Length, A.Length, ["q", "ref"]);
        Assert.That(text, Does.Contain("ref                1 TTACGTTT      8\n"));
        Assert.That(text, Does.Contain("q                  1 --ACGT--      4\n"));
    }

    [Test]
    public void SemiGlobalAlign_Markx10_HasTheWholeLengths()
    {
        string text = Emboss(AlignmentTools.SemiGlobalAlign(B, A), "markx10", B.Length, A.Length, ["q", "ref"]);
        Assert.That(text, Does.Contain(">>>q, 4 nt vs ref, 8 nt"));
        Assert.That(text, Does.Contain("; al_start: 1\n; al_stop: 8\n"));
    }

    [Test]
    public void DualAffineAlign_Local_Markx10_HasTheCoordinatesOfWater()
    {
        string text = Emboss(AlignmentTools.DualAffineAlign(A, B, mode: "local"), "markx10", A.Length, B.Length, Names);
        Assert.That(text, Does.Contain(">>>a, 8 nt vs b, 4 nt"));
        Assert.That(text, Does.Contain("; al_start: 3\n; al_stop: 6\n"));
    }

    [Test]
    public void LocalAlign_EmbossSamIsTheEmbossRecord_NotTheBiopythonOne()
    {
        string text = Emboss(AlignmentTools.LocalAlign(A, B), "sam", A.Length, B.Length, Names);
        Assert.That(text, Does.Contain("@HD").Or.Contain("b\t0\ta\t3\t255\t4M\t"));
        Assert.That(text, Is.Not.EqualTo(AlignmentTools.LocalAlign(A, B, format: "sam", names: Names).Formatted));
    }

    [Test]
    public void GlobalAlign_Pair_Needle_Stretcher_AgreeOnTheFullRows()
    {
        foreach (string text in new[]
        {
            Emboss(AlignmentTools.GlobalAlign(A, B), "pair", A.Length, B.Length, Names),
            EmbossPairwiseAlignmentTools.NeedleAlign(A, B, format: "pair", names: Names).Formatted!,
            EmbossPairwiseAlignmentTools.StretcherAlign(A, B, format: "pair", names: Names).Formatted!,
        })
        {
            Assert.That(text, Does.Contain("a                  1 "), text);
            Assert.That(text, Does.Contain("      8\n"), text);
            Assert.That(text, Does.Contain("b                  1 ").And.Contain("      4\n"), text);
        }
    }
}
