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

    // ------------------------------------------------------------------ A3-5: the pairwise programs, ksw2, T-Coffee DPs, MUSCLE refinement, anchors, distmat (B09 F133)

    [Test]
    public void PreCancelled_A3_5_EveryPortedProgram_Throws()
    {
        const string a = "ACGTACGTTAGCACGTAGCTAGCTAGGATC", b = "ACGTCGTTAGCACGTAGCTTAGCTAGGATC";
        var one = new[] { new EmbossSequence("a", a) };
        var q = Ksw2Aligner.EncodeNucleotides(a);
        var t = Ksw2Aligner.EncodeNucleotides(b);
        var mat = Ksw2Aligner.SimpleMatrix(5, 2, 4, 1);
        Assert.Multiple(() =>
        {
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => EmbossPairwisePrograms.NeedleAlign(a, b)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => EmbossPairwisePrograms.WaterAlign(a, b)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => EmbossPairwisePrograms.StretcherAlign(a, b)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => EmbossPairwisePrograms.Supermatcher(one, one)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => EmbossPairwisePrograms.NeedleAll(one, one)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => Ksw2Aligner.ExtendZ(q, t, 5, mat, 4, 2, -1, -1)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => Ksw2Aligner.ExtendD(q, t, 5, mat, 4, 2, 24, 1, -1, -1)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => Ksw2Aligner.ExtendZ2(q, t, 5, mat, 4, 2, -1, -1, 0)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => Ksw2Aligner.ExtendD2(q, t, 5, mat, 4, 2, 24, 1, -1, -1, 0)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => SequenceAligner.GlobalAlignDualAffine(a, b, DualAffineScoring.Minimap2)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeSlowPairAlign(a, b)));
            Assert.Catch<OperationCanceledException>(() => SequenceAligner.WithCancellation(Cancelled, () => EmbossAlignmentStatistics.CalculateDistanceMatrix(new[] { "ACGT", "ACGA", "ACCA" }, new MsaDistanceOptions())));
        });
    }

    [Test]
    public void CancelledWhileRunning_Needle_Water_StopPromptly()
    {
        string a = Random(600, 51), b = Random(600, 52);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => EmbossPairwisePrograms.NeedleAlign(a, b)), pollInside: "NeedleCore");
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => EmbossPairwisePrograms.WaterAlign(a, b)), pollInside: "WaterCore");
    }

    [Test]
    public void CancelledWhileRunning_Stretcher_StopsPromptly()
    {
        string a = Random(3000, 53), b = Random(3000, 54);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => EmbossPairwisePrograms.StretcherAlign(a, b)), pollInside: "StretcherEngine");
    }

    [Test]
    public void CancelledWhileRunning_Supermatcher_StopsPromptly()
    {
        var seq = new[] { new EmbossSequence("a", Random(3000, 55)) };
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => EmbossPairwisePrograms.Supermatcher(seq, seq, width: 600)), pollInside: "SwPathCalcFast");
    }

    [Test]
    public void CancelledWhileRunning_NeedleAll_StopsPromptly()
    {
        var set = Enumerable.Range(0, 60).Select(i => new EmbossSequence("s" + i, Random(600, 200 + i))).ToArray();
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => EmbossPairwisePrograms.NeedleAll(set, set)), pollInside: "NeedleAll");
    }

    [Test]
    public void CancelledWhileRunning_DistanceMatrix_StopsPromptly()
    {
        var rows = Enumerable.Range(0, 40).Select(i => Random(800, 400 + i)).ToArray();
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => EmbossAlignmentStatistics.CalculateDistanceMatrix(rows, new MsaDistanceOptions())), pollInside: "CalcMatch");
    }

    [Test]
    public void CancelledWhileRunning_DotMatcher_StopsPromptly()
    {
        string a = Random(3000, 73), b = Random(3000, 74);
        AssertCancelsPromptly(ct => DotMatcher.Compute(a, b, cancellationToken: ct), pollInside: "Diagonal");
    }

    [Test]
    public void CancelledWhileRunning_Merger_StopsPromptly()
    {
        string a = Random(600, 77), b = Random(600, 78);
        AssertCancelsPromptly(ct => EmbossPairwisePrograms.Merger(a, b, cancellationToken: ct), pollInside: "MergerPathCalc");
    }

    [Test]
    public void CancelledWhileRunning_TopNSim_StopsPromptly()
    {
        string a = Random(600, 79), b = Random(600, 80);
        AssertCancelsPromptly(ct => SequenceAligner.WithCancellation(ct, () => EmbossPairwisePrograms.LocalAlignAffineTopN(a, b, new ScoringMatrix(5, -4, -10, -4), 40, 0, TopNTieOrder.Sim)), pollInside: "SimMatcher");
    }

    [Test]
    public void CancelledWhileRunning_MegaMergeMinimalMatchSet_StopsPromptly()
    {
        // 3 000 matches that all survive (the same diagonal, 100 apart): every kept match is compared with every later one
        var list = Enumerable.Range(0, 3_000).Select(i => new EmbossPairwisePrograms.MegaMatch(100 * i, 100 * i, 30)).ToList();
        AssertCancelsPromptly(ct => { EmbossPairwisePrograms.MinimalMatchSet(list, 20, ct); return list; }, pollInside: "MinimalMatchSet");
    }

    [Test]
    public void CancelledWhileRunning_ShowAlign_StopsPromptly()
    {
        // the consensus of 150 000 columns is computed first (its own poll), then 150 000 one-column blocks are written (the poll of Output)
        string row1 = Random(150_000, 82), row2 = Random(150_000, 83);
        AssertCancelsPromptly(ct => MultipleAlignmentDisplay.ShowAlign(["a", "b"], [row1, row2], new ShowAlignOptions { Width = 1, CancellationToken = ct }), pollInside: "Output");
        AssertCancelsPromptly(ct => MultipleAlignmentDisplay.ShowAlign(["a", "b"], [row1, row2], new ShowAlignOptions { Width = 1, CancellationToken = ct }), pollInside: "EmbossConsensusCore");
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
