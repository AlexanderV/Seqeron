// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>
/// Audit A3-5 (B09 F133): the MCP wrappers of the ported programs hand the request token to the library (<c>SequenceAligner.BeginCancellationScope</c>), so the row / iteration polls that A3-5
/// added to needle / water / stretcher / supermatcher / needleall / ksw2 / dual affine / T-Coffee / MUSCLE refinement / anchors take effect through MCP. A token that is already cancelled
/// must stop every one of these tools with <see cref="OperationCanceledException"/>.
/// </summary>
[TestFixture]
public class CancellationTokenTests
{
    private static readonly CancellationToken Cancelled = new(true);
    private const string A = "ACGTACGTTAGCACGTAGCTAGCTAGGATC";
    private const string B = "ACGTCGTTAGCACGTAGCTTAGCTAGGATC";

    [Test]
    public void AnUncancelledToken_GivesTheSameResultAsNoToken()
    {
        using var cts = new CancellationTokenSource();
        Assert.Multiple(() =>
        {
            Assert.That(EmbossPairwiseAlignmentTools.StretcherAlign(A, B, cancellationToken: cts.Token).AlignedSequence1, Is.EqualTo(EmbossPairwiseAlignmentTools.StretcherAlign(A, B).AlignedSequence1));
            Assert.That(AlignmentTools.Ksw2Extend(A, B, cancellationToken: cts.Token).Score, Is.EqualTo(AlignmentTools.Ksw2Extend(A, B).Score));
        });
    }

    /// <summary>A4-1 (B09 F157): dotmatcher, minimap2_chains and merger took no token; they now hand the request token to the ports' ambient polls.</summary>
    [Test]
    public void A4_1_DotmatcherMinimap2ChainsMerger_StopWithACancelledToken()
    {
        const string a = "AAAAACGTACGTTAGCACGTAGCTAGCTAGGATCAGGTCAAGTCAGTCGATTC";
        string b = a.Replace('G', 'T');
        Assert.Multiple(() =>
        {
            Assert.Catch<OperationCanceledException>(() => EmbossPairwiseAlignmentTools.Dotmatcher(a, b, cancellationToken: Cancelled), "dotmatcher");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.Minimap2Chains([a], b, cancellationToken: Cancelled), "minimap2_chains");
            Assert.Catch<OperationCanceledException>(() => EmbossPairwiseAlignmentTools.Merger(a, b, cancellationToken: Cancelled), "merger");
        });
    }
}
