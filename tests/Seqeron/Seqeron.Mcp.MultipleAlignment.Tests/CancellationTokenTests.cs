// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.MultipleAlignment.Tools;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>
/// Audit A3-5 (B09 F133): the MCP wrappers of the ported programs hand the request token to the library (<c>SequenceAligner.BeginCancellationScope</c>), so the row / iteration polls that A3-5
/// added to needle / water / stretcher / supermatcher / needleall / ksw2 / dual affine / T-Coffee / MUSCLE refinement / anchors take effect through MCP. A token that is already cancelled
/// must stop every one of these tools with <see cref="OperationCanceledException"/>.
/// </summary>
[TestFixture]
public class CancellationTokenTests
{
    private const string A = "ACGTACGTTAGCACGTAGCTAGCTAGGATC";
    private static readonly CancellationToken Cancelled = new(true);
    private const string B = "ACGTCGTTAGCACGTAGCTTAGCTAGGATC";

    [Test]
    public void EveryPortedProgramTool_StopsWithACancelledToken()
    {
        Assert.Multiple(() =>
        {
            Assert.Catch<OperationCanceledException>(() => EmbossPairwiseAlignmentTools.StretcherAlign(A, B, cancellationToken: Cancelled), "stretcher_align");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.DualAffineAlign(A, B, cancellationToken: Cancelled), "dual_affine_align");
            Assert.Catch<OperationCanceledException>(() => EmbossPairwiseAlignmentTools.SupermatcherAlign([A], [A], cancellationToken: Cancelled), "supermatcher_align");
            Assert.Catch<OperationCanceledException>(() => EmbossPairwiseAlignmentTools.NeedleallAlign([A], [B], cancellationToken: Cancelled), "needleall_align");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.Ksw2Extend(A, B, cancellationToken: Cancelled), "ksw2_extend");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.Ksw2Extend(A, B, variant: "extz", cancellationToken: Cancelled), "ksw2_extend extz");
            Assert.Catch<OperationCanceledException>(() => MultipleAlignmentTools.TCoffeeSlowPairAlign(A, B, cancellationToken: Cancelled), "tcoffee_slow_pair_align");
            Assert.Catch<OperationCanceledException>(() => MultipleAlignmentTools.TCoffeePrimaryLibrary([A, B, A[..20]], cancellationToken: Cancelled), "tcoffee_primary_library");
            Assert.Catch<OperationCanceledException>(() => MultipleAlignmentTools.MultipleAlignTCoffee([A, B, A[..20]], cancellationToken: Cancelled), "multiple_align_tcoffee");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.MultipleAlignMuscle([A, B, A[..20]], maxIters: 4, cancellationToken: Cancelled), "multiple_align_muscle (refinement)");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.MultipleAlignAnchored([A, B, A[..20]], cancellationToken: Cancelled), "multiple_align_anchored");
        });
    }}
