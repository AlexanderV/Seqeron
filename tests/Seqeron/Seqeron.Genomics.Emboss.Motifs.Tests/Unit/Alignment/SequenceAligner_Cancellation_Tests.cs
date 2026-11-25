// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 / ALIGN-MULTI-001 / GENOMIC-* (audit A2-40, B09 F86): cooperative cancellation of the long-running APIs.
///
/// Contract: a call whose token is cancelled (before it starts, or while it runs) throws <see cref="OperationCanceledException"/> promptly - the dynamic
/// programs poll the token at their row / iteration boundaries - and returns no partial result; the token never leaks into a later call. The token is
/// supplied through <c>ConsistencyMsaOptions / ClustalWOptions / ClustalWProteinOptions / MuscleOptions / TCoffeeProgramOptions.CancellationToken</c> (the convention of
/// <see cref="System.Threading.Tasks.ParallelOptions.CancellationToken"/>), through <see cref="SequenceAligner.WithCancellation{T}(CancellationToken, Func{T})"/> for
/// every other aligner, and through the trailing <c>cancellationToken</c> parameter of the <see cref="GenomicAnalyzer"/> scans.
/// </summary>
[TestFixture]
// Audit A3-15: runs alone. The fixture builds 3 M-base suffix trees (~0.5 GB each) and measures how long cancelled calls take to stop; next to seven other workers (LevelOfParallelism 8 on 4 cores) the
// uninterruptible tree builds took longer than the limit once (flake "the cancelled call did not stop within the limit") and the sum of the fixtures' working sets once ended in an OOM kill.
[NonParallelizable]
public class SequenceAligner_Cancellation_Tests
{

    // ------------------------------------------------------------------ A4-1: the ports added in WP36 - WP38 poll the ambient token (one test per class: cancelled while running, the poll named)

    [Test]
    public void PreCancelled_A4_1_EveryPortAddedInWp36ToWp38_Throws()
    {
        string a = Random(300, 71), b = Random(300, 72);
        var anchors = new[] { new MinimapAnchor(0, false, 100, 100, 15), new MinimapAnchor(0, false, 120, 120, 15) };
        var matches = new (long, long, long)[] { (1, 1, 20), (31, 31, 20) };
        Assert.Multiple(() =>
        {
            Assert.Catch<OperationCanceledException>(() => DotMatcher.Compute(a, b, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => Minimap2Chaining.Chain(anchors, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => Minimap2Chaining.CollectAnchors([a], b, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => Minimap2Chaining.MapChains([a], b, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => NucmerClusters.Cluster(matches, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => EmbossPairwisePrograms.Merger(a, b, cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => EmbossPairwisePrograms.LocalAlignAffineTopN(a, b, new ScoringMatrix(5, -4, -10, -4), 3, 0, TopNTieOrder.Sim)));
            Assert.Catch<OperationCanceledException>(() => EmbossMotifPrograms.FindTranslatedPatternHits(a, "L-x-S", "F", 0, Cancelled));
            Assert.Catch<OperationCanceledException>(() => MultipleAlignmentDisplay.ShowAlign(["a", "b"], [a, b], new ShowAlignOptions { CancellationToken = Cancelled }));
            Assert.Catch<OperationCanceledException>(() => EmbossMotifPrograms.GenerateEmbossAmbiguousConsensus([a, b], cancellationToken: Cancelled));
            Assert.Catch<OperationCanceledException>(() => GenomicAnalyzer.FindNucmerClusters(a, b, 8, cancellationToken: Cancelled));
        });
    }

    [Test]
    public void CancelledWhileRunning_FuzztranTranslatedPatternHits_StopsPromptly()
    {
        string nt = Random(600_000, 81);
        AssertCancelsPromptly(ct => EmbossMotifPrograms.FindTranslatedPatternHits(nt, "L-x-S", "F", 0, ct), pollInside: "FindTranslatedPatternHits");
    }

    [Test]
    public void CancelledWhileRunning_FuzzproProteinPatternHits_StopsPromptly_InsideThePositionLoop()
    {
        string protein = Random(300_000, 86, "ACDEFGHIKLMNPQRSTVWY");
        AssertCancelsPromptly(ct => EmbossMotifPrograms.FindProteinPatternHits(protein, "L-x-S", ct), pollInside: "ProteinPatternHitsCore");
    }

    [Test]
    public void CancelledWhileRunning_Consambig_StopsPromptly()
    {
        var rows = Enumerable.Range(0, 20).Select(i => Random(40_000, 90 + i)).ToList();
        AssertCancelsPromptly(ct => EmbossMotifPrograms.GenerateEmbossAmbiguousConsensus(rows, cancellationToken: ct), pollInside: "GenerateEmbossAmbiguousConsensus");
    }

    /// <summary>
    /// A5-4: the call runs ONCE with a token that is cancelled at the third poll (of <paramref name="pollInside"/> when named) - <see cref="PollCancellation"/>: deterministic, no sleep, no thread, no time bound (the
    /// earlier version cancelled after a wall-clock delay, which under CPU load landed in another stage). The call must stop with <see cref="OperationCanceledException"/> raised by that poll and run no later poll.
    /// </summary>
    private static void AssertCancelsPromptly(Func<CancellationToken, object?> work, string pollInside, int nth = 3)
    {
        string stack = PollCancellation.CancelInside(work, pollInside, nth);
        Assert.That(stack, Does.Contain(pollInside), "the cancelling poll is inside " + pollInside);
    }

    private static readonly CancellationToken Cancelled = new(true);

    private static string Random(int n, int seed, string alphabet = "ACGT")
    {
        var rng = new System.Random(seed);
        return new string(Enumerable.Range(0, n).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
    }
}
