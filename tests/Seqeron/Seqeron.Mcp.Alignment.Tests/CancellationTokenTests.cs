using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

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
    public void EveryPortedProgramTool_StopsWithACancelledToken()
    {
        Assert.Multiple(() =>
        {
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.StretcherAlign(A, B, cancellationToken: Cancelled), "stretcher_align");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.DualAffineAlign(A, B, cancellationToken: Cancelled), "dual_affine_align");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.SupermatcherAlign([A], [A], cancellationToken: Cancelled), "supermatcher_align");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.NeedleallAlign([A], [B], cancellationToken: Cancelled), "needleall_align");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.Ksw2Extend(A, B, cancellationToken: Cancelled), "ksw2_extend");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.Ksw2Extend(A, B, variant: "extz", cancellationToken: Cancelled), "ksw2_extend extz");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.TCoffeeSlowPairAlign(A, B, cancellationToken: Cancelled), "tcoffee_slow_pair_align");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.TCoffeePrimaryLibrary([A, B, A[..20]], cancellationToken: Cancelled), "tcoffee_primary_library");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.MultipleAlignTCoffee([A, B, A[..20]], cancellationToken: Cancelled), "multiple_align_tcoffee");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.MultipleAlignMuscle([A, B, A[..20]], maxIters: 4, cancellationToken: Cancelled), "multiple_align_muscle (refinement)");
            Assert.Catch<OperationCanceledException>(() => AlignmentTools.MultipleAlignAnchored([A, B, A[..20]], cancellationToken: Cancelled), "multiple_align_anchored");
        });
    }

    [Test]
    public void AnUncancelledToken_GivesTheSameResultAsNoToken()
    {
        using var cts = new CancellationTokenSource();
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.StretcherAlign(A, B, cancellationToken: cts.Token).AlignedSequence1, Is.EqualTo(AlignmentTools.StretcherAlign(A, B).AlignedSequence1));
            Assert.That(AlignmentTools.Ksw2Extend(A, B, cancellationToken: cts.Token).Score, Is.EqualTo(AlignmentTools.Ksw2Extend(A, B).Score));
        });
    }
}
