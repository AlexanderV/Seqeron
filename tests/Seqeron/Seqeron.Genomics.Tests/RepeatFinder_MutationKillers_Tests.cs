// Mutation killers for RepeatFinder.cs (rows 13 REP-STR-001, 17 REP-PALIN-001, 256 REP-APPROX-001).
//
// Stryker baseline 76.55% (< 80% target). The surviving mutants clustered in code paths whose existing
// fixtures only assert "a repeat was found" rather than the exact tract geometry, plus the cancellable
// FindMicrosatellites overload and the TRF Bernoulli / tandem-repeat-summary helpers. These tests pin the
// documented business properties — exact repeat unit / copy count / span, exact loop geometry of a hairpin,
// exact Bernoulli match/indel probabilities, and the per-type repeat counts of the summary — so that the
// arithmetic, relational, boundary and null-coalescing mutants on those lines diverge from the truth.
//
// Evidence / theory:
//   * Microsatellites (STRs) are 1-6 bp motifs repeated >= minRepeats times consecutively (Ellegren 2004,
//     Nat Rev Genet 5:435-445). A perfect run of n bases of a k-mer is exactly floor(n/k) copies.
//   * Inverted repeat / hairpin: the right arm is the reverse complement of the left arm; the intervening
//     loop forms a stem-loop only when it is long enough (>= 3 nt here).
//   * Tandem Repeats Finder Bernoulli model (Benson 1999, NAR 27(2):573-580): two adjacent copies are
//     modelled as independent Bernoulli trials; PM (matching probability) is the fraction of matching
//     columns between adjacent copies and PI the fraction of indel columns; a perfect tract has PM = 1.

using System;
using System.Linq;
using System.Threading;
using NUnit.Framework;
using Seqeron.Genomics;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class RepeatFinder_MutationKillers_Tests
{
    private const double Tol = 1e-9;

    #region FindMicrosatellites — cancellable overload, exact tract (kills extension/threshold/length mutants)

    // A perfect (AT) x3 run that ENDS exactly at the sequence end. The last copy is only counted when the
    // extension bound is inclusive (j + unitLen <= seq.Length); a strict "<" stops one copy short and reports
    // only 2 copies (below minRepeats=3), so the repeat disappears. Pins copy count, span and total length.
    [Test]
    public void FindMicrosatellites_Cancellable_PerfectDinucleotide_ReportsExactTract()
    {
        var results = RepeatFinder
            .FindMicrosatellites("ATATAT", 1, 6, 3, CancellationToken.None)
            .ToList();

        Assert.That(results, Has.Count.EqualTo(1), "exactly one (AT)x3 tract, no shorter sub-repeats");
        var r = results[0];
        Assert.Multiple(() =>
        {
            Assert.That(r.Position, Is.EqualTo(0), "tract starts at index 0");
            Assert.That(r.RepeatUnit, Is.EqualTo("AT"), "period-2 motif AT");
            Assert.That(r.RepeatCount, Is.EqualTo(3), "AT AT AT = 3 consecutive copies");
            Assert.That(r.TotalLength, Is.EqualTo(6), "TotalLength = RepeatCount * unitLength = 3 * 2 = 6");
            Assert.That(r.RepeatType, Is.EqualTo(RepeatType.Dinucleotide), "2 bp unit -> Dinucleotide");
        });
    }

    // A homopolymer A x6 is the maximal mononucleotide tract; every shorter run starting later (A x5, x4, x3)
    // is CONTAINED in it and must be suppressed by the containment filter. A single result must remain. This
    // pins both the end-coordinate arithmetic (end = i + repeats*unitLen - 1) and the containment comparison
    // (r.Start <= i && r.End >= end): corrupting either re-admits the contained sub-runs.
    [Test]
    public void FindMicrosatellites_Cancellable_Homopolymer_SuppressesContainedSubRuns()
    {
        var results = RepeatFinder
            .FindMicrosatellites("AAAAAA", 1, 6, 3, CancellationToken.None)
            .ToList();

        Assert.That(results, Has.Count.EqualTo(1), "the maximal A x6 run subsumes all shorter contained runs");
        var r = results[0];
        Assert.Multiple(() =>
        {
            Assert.That(r.Position, Is.EqualTo(0));
            Assert.That(r.RepeatUnit, Is.EqualTo("A"));
            Assert.That(r.RepeatCount, Is.EqualTo(6), "6 consecutive A copies");
            Assert.That(r.TotalLength, Is.EqualTo(6));
            Assert.That(r.RepeatType, Is.EqualTo(RepeatType.Mononucleotide));
        });
    }

    // A non-repetitive sequence has NO microsatellite at minRepeats=3. If the report threshold were inverted
    // to "repeats < minRepeats", every single position (1 copy) would be emitted, so the result would be
    // non-empty. Pins the "repeats >= minRepeats" reporting rule.
    [Test]
    public void FindMicrosatellites_Cancellable_NonRepetitive_ReturnsEmpty()
    {
        var results = RepeatFinder
            .FindMicrosatellites("ATGC", 1, 6, 3, CancellationToken.None)
            .ToList();

        Assert.That(results, Is.Empty, "no motif repeats >= 3 times in ATGC");
    }

    // The cancellable string overload validates its numeric parameters eagerly (these guards were uncovered).
    [Test]
    public void FindMicrosatellites_Cancellable_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RepeatFinder.FindMicrosatellites("ACGT", 0, 6, 3, CancellationToken.None).ToList(),
                "minUnitLength < 1 is rejected");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RepeatFinder.FindMicrosatellites("ACGT", 4, 2, 3, CancellationToken.None).ToList(),
                "maxUnitLength < minUnitLength is rejected");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RepeatFinder.FindMicrosatellites("ACGT", 1, 6, 1, CancellationToken.None).ToList(),
                "minRepeats < 2 is rejected");
        });
    }

    #endregion

    #region FindMicrosatellites — non-cancellable core, contained-run suppression

    [Test]
    public void FindMicrosatellites_Core_Homopolymer_SuppressesContainedSubRuns()
    {
        var results = RepeatFinder.FindMicrosatellites("AAAAAA", 1, 6, 3).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].RepeatCount, Is.EqualTo(6));
    }

    #endregion

    #region FindInvertedRepeats — exact hairpin geometry (kills loop-length / CanFormHairpin mutants)

    // GGGGG ... CCCCC where CCCCC = revcomp(GGGGG). With a 3 nt loop the stem-loop is hairpin-capable.
    [Test]
    public void FindInvertedRepeats_StemLoop_ReportsExactArmsAndLoop()
    {
        var results = RepeatFinder
            .FindInvertedRepeats("GGGGGTTTCCCCC", minArmLength: 5, maxLoopLength: 50, minLoopLength: 3)
            .ToList();

        var hairpin = results.SingleOrDefault(r => r.LeftArmStart == 0 && r.ArmLength == 5);
        Assert.That(hairpin, Is.Not.EqualTo(default(InvertedRepeatResult)), "the 5 bp stem must be detected");
        Assert.Multiple(() =>
        {
            Assert.That(hairpin.RightArmStart, Is.EqualTo(8), "right arm begins after the 5 bp stem + 3 nt loop");
            Assert.That(hairpin.LeftArm, Is.EqualTo("GGGGG"));
            Assert.That(hairpin.RightArm, Is.EqualTo("CCCCC"), "right arm is the reverse complement of the left");
            Assert.That(hairpin.LoopLength, Is.EqualTo(3), "loop = RightArmStart - (LeftArmStart + ArmLength) = 8 - 5");
            Assert.That(hairpin.Loop, Is.EqualTo("TTT"));
            Assert.That(hairpin.CanFormHairpin, Is.True, "loop length 3 >= 3 -> hairpin-capable");
            Assert.That(hairpin.TotalLength, Is.EqualTo(13), "2 * arm(5) + loop(3) = 13");
        });
    }

    // A 1 nt loop is too short to fold: CanFormHairpin must be FALSE. This kills the "loopLength >= 3"
    // boundary being relaxed to "loopLength >= 0".
    [Test]
    public void FindInvertedRepeats_LoopTooShort_CannotFormHairpin()
    {
        var results = RepeatFinder
            .FindInvertedRepeats("GGGGGTCCCCC", minArmLength: 5, maxLoopLength: 50, minLoopLength: 1)
            .ToList();

        var stem = results.Single(r => r.LeftArmStart == 0 && r.ArmLength == 5);
        Assert.Multiple(() =>
        {
            Assert.That(stem.LoopLength, Is.EqualTo(1), "single-nt loop");
            Assert.That(stem.CanFormHairpin, Is.False, "loop length 1 < 3 -> cannot form a hairpin");
        });
    }

    #endregion

    #region ComputeBernoulliStatistics — TRF Bernoulli model (Benson 1999)

    // A perfect tandem tract: every column between adjacent copies is a match, so PM = 1, PI = 0 and the
    // Bernoulli-mean expected matches equal the trial count.
    [Test]
    public void ComputeBernoulliStatistics_PerfectTract_PmIsOne()
    {
        var stats = RepeatFinder.ComputeBernoulliStatistics("ACGTACGTACGT", period: 4);

        Assert.Multiple(() =>
        {
            Assert.That(stats.AdjacentCopyPairs, Is.EqualTo(2), "3 copies -> 2 adjacent pairs");
            Assert.That(stats.Matches, Is.EqualTo(8), "2 pairs x 4 matching columns");
            Assert.That(stats.Mismatches, Is.EqualTo(0));
            Assert.That(stats.Indels, Is.EqualTo(0));
            Assert.That(stats.BernoulliTrials, Is.EqualTo(8), "total alignment columns between adjacent copies");
            Assert.That(stats.MatchProbability, Is.EqualTo(1.0).Within(Tol), "PM = matches / trials = 8/8");
            Assert.That(stats.IndelProbability, Is.EqualTo(0.0).Within(Tol));
            Assert.That(stats.PercentMatches, Is.EqualTo(100.0).Within(Tol));
            Assert.That(stats.ExpectedMatches, Is.EqualTo(8.0).Within(Tol), "E[heads] = PM * trials = 1 * 8");
            Assert.That(stats.MeetsExpectedMatchProbability, Is.True, "PM 1.0 >= default 0.80");
        });
    }

    // The inclusive [0, 1] bounds of expectedMatchProbability must be ACCEPTED (the boundary mutants relax
    // "< 0 || > 1" to "<= 0 || >= 1", which would wrongly reject the valid endpoints 0.0 and 1.0). period = 1
    // is likewise valid (>= 1) and must not throw.
    [Test]
    public void ComputeBernoulliStatistics_BoundaryArguments_AreAccepted()
    {
        Assert.Multiple(() =>
        {
            Assert.DoesNotThrow(() => RepeatFinder.ComputeBernoulliStatistics("ACGTACGT", 4, 0.0),
                "expectedMatchProbability = 0.0 is a valid lower bound");
            Assert.DoesNotThrow(() => RepeatFinder.ComputeBernoulliStatistics("ACGTACGT", 4, 1.0),
                "expectedMatchProbability = 1.0 is a valid upper bound");
            Assert.DoesNotThrow(() => RepeatFinder.ComputeBernoulliStatistics("AA", 1),
                "period = 1 is valid");
        });
    }

    [Test]
    public void ComputeBernoulliStatistics_InvalidArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RepeatFinder.ComputeBernoulliStatistics("ACGTACGT", 0), "period < 1 is rejected");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RepeatFinder.ComputeBernoulliStatistics("ACGTACGT", 4, -0.1), "probability below 0 is rejected");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RepeatFinder.ComputeBernoulliStatistics("ACGTACGT", 4, 1.1), "probability above 1 is rejected");
        });
    }

    #endregion

    #region GetTandemRepeatSummary — per-type counts (kills null-coalescing remove-left mutants)

    // A sequence containing mononucleotide, dinucleotide and trinucleotide tracts (but no tetranucleotide).
    // Each per-type count must reflect the PRESENT group's real, non-zero count, while the absent
    // tetranucleotide group reports 0. The "byType.GetValueOrDefault(type)?.Count ?? 0" mutants that drop
    // the left operand would collapse every present count to 0, diverging from the asserted values.
    // The di/tri counts are 2, not 1, because overlapping-but-not-contained frames are each reported
    // (e.g. (CG)x3 at index 6 and the offset (GC)x3 at index 7 are distinct, non-contained tracts).
    [Test]
    public void GetTandemRepeatSummary_MixedTypes_CountsPresentTypesNonZero()
    {
        // A x6 (mono) | CGCGCG (di frames) | CAGCAGCAG (tri frames)
        var summary = RepeatFinder.GetTandemRepeatSummary(new DnaSequence("AAAAAACGCGCGCAGCAGCAG"), 3);

        Assert.Multiple(() =>
        {
            Assert.That(summary.MononucleotideRepeats, Is.EqualTo(1), "one A x6 mononucleotide tract");
            Assert.That(summary.DinucleotideRepeats, Is.EqualTo(2), "two non-contained dinucleotide frames");
            Assert.That(summary.TrinucleotideRepeats, Is.EqualTo(2), "two non-contained trinucleotide frames");
            Assert.That(summary.TetranucleotideRepeats, Is.EqualTo(0), "no tetranucleotide tract present");
        });
    }

    #endregion
}
