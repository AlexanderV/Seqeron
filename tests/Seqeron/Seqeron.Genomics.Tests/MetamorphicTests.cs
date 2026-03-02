using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for algorithms changed in March 2026 commits:
///   PAT-IUPAC-001, PAT-PWM-001, PAT-APPROX-002,
///   REP-STR-001, REP-TANDEM-001, REP-INV-001, REP-DIRECT-001.
///
/// Each test encodes a metamorphic relation (MR) — a property that relates
/// the outputs of multiple executions under input transformations,
/// without requiring a test oracle.
///
/// Relations verified:
///   - Subset/hierarchy (degeneracy, tolerance, range widening)
///   - Monotonicity (threshold tightening/relaxation)
///   - Invariance under input transformation (shift, permutation, case)
///   - Symmetry (reverse complement)
///   - Metric axioms (non-negativity)
///   - Composition (exact ⊆ approximate)
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class MetamorphicTests
{
    #region Helpers

    private static readonly Random Rng = new(42);

    /// <summary>Generates a random DNA string of given length.</summary>
    private static string RandomDna(int length)
    {
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[Rng.Next(4)];
        return new string(chars);
    }

    /// <summary>Generates a DNA string that contains a guaranteed microsatellite.</summary>
    private static string DnaWithMicrosatellite(string unit, int repeats, int flankLen = 10)
    {
        return RandomDna(flankLen) + string.Concat(Enumerable.Repeat(unit, repeats)) + RandomDna(flankLen);
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-IUPAC-001 — IUPAC Degenerate Matching
    // ═══════════════════════════════════════════════════════════════════

    #region MR1: IUPAC degeneracy hierarchy — more degenerate code ⊇ less degenerate

    /// <summary>
    /// MR1-a: For any base b, MatchesIupac(b, N) = true whenever MatchesIupac(b, X) = true
    /// for any less-degenerate code X.
    /// N matches all 4 bases, so its match set is a superset of every other code's match set.
    /// </summary>
    [Test]
    public void MR1a_Iupac_N_SupersetOfEveryCode()
    {
        char[] allBases = { 'A', 'C', 'G', 'T' };
        char[] allCodes = { 'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N' };

        foreach (char code in allCodes)
        {
            foreach (char b in allBases)
            {
                if (IupacHelper.MatchesIupac(b, code))
                {
                    IupacHelper.MatchesIupac(b, 'N').Should().BeTrue(
                        because: $"N is the most degenerate code and must match '{b}' if '{code}' matches it");
                }
            }
        }
    }

    /// <summary>
    /// MR1-b: Degeneracy hierarchy for purines: matches(R) ⊇ matches(A) and matches(R) ⊇ matches(G).
    /// R = {A, G}, so any base matching A alone or G alone must also match R.
    /// </summary>
    [Test]
    public void MR1b_Iupac_R_SupersetOfA_And_G()
    {
        char[] allBases = { 'A', 'C', 'G', 'T' };

        foreach (char b in allBases)
        {
            if (IupacHelper.MatchesIupac(b, 'A'))
                IupacHelper.MatchesIupac(b, 'R').Should().BeTrue(
                    because: $"R ⊇ A, so '{b}' matching A must also match R");

            if (IupacHelper.MatchesIupac(b, 'G'))
                IupacHelper.MatchesIupac(b, 'R').Should().BeTrue(
                    because: $"R ⊇ G, so '{b}' matching G must also match R");
        }
    }

    /// <summary>
    /// MR1-c: Three-base codes are supersets of their constituent two-base codes.
    /// B = not-A = {C,G,T}, so B ⊇ Y(={C,T}), B ⊇ S(={G,C}), B ⊇ K(={G,T}).
    /// </summary>
    [Test]
    public void MR1c_Iupac_ThreeBase_SupersetOfTwoBase()
    {
        // B = {C,G,T} ⊇ Y={C,T}, S={G,C}, K={G,T}
        // D = {A,G,T} ⊇ R={A,G}, W={A,T}, K={G,T}
        // H = {A,C,T} ⊇ Y={C,T}, W={A,T}, M={A,C}
        // V = {A,C,G} ⊇ R={A,G}, S={G,C}, M={A,C}
        var supersets = new (char ThreeBase, char[] TwoBases)[]
        {
            ('B', new[] { 'Y', 'S', 'K' }),
            ('D', new[] { 'R', 'W', 'K' }),
            ('H', new[] { 'Y', 'W', 'M' }),
            ('V', new[] { 'R', 'S', 'M' }),
        };

        char[] allBases = { 'A', 'C', 'G', 'T' };

        foreach (var (threeBase, twoBases) in supersets)
        {
            foreach (char twoBase in twoBases)
            {
                foreach (char b in allBases)
                {
                    if (IupacHelper.MatchesIupac(b, twoBase))
                    {
                        IupacHelper.MatchesIupac(b, threeBase).Should().BeTrue(
                            because: $"'{threeBase}' ⊇ '{twoBase}', so '{b}' matching '{twoBase}' must match '{threeBase}'");
                    }
                }
            }
        }
    }

    #endregion

    #region MR2: IUPAC degenerate motif — more degenerate pattern yields ≥ matches

    /// <summary>
    /// MR2: Replacing a concrete base in a pattern with a more degenerate IUPAC code
    /// can only increase (or maintain) the number of matches.
    /// Example: pattern "ATG" → "RTG" (A→R). R matches {A,G}, so matches("RTG") ⊇ matches("ATG").
    /// </summary>
    [Test]
    public void MR2_DegenerateMotif_MoreDegenerate_MoreOrEqualMatches()
    {
        var sequence = new DnaSequence("ATGATGGTGCTGATG");

        // Concrete pattern
        var concreteMatches = MotifFinder.FindDegenerateMotif(sequence, "ATG").ToList();

        // More degenerate: A→R (R matches A and G)
        var degenerateMatches = MotifFinder.FindDegenerateMotif(sequence, "RTG").ToList();

        degenerateMatches.Count.Should().BeGreaterThanOrEqualTo(concreteMatches.Count,
            because: "replacing A with R (⊇ A) can only add matches, not remove them");

        // Every position from concrete must appear in degenerate
        var degeneratePositions = degenerateMatches.Select(m => m.Position).ToHashSet();
        foreach (var m in concreteMatches)
        {
            degeneratePositions.Should().Contain(m.Position,
                because: "a wider IUPAC code must still match everywhere the narrower one did");
        }
    }

    /// <summary>
    /// MR2-b: Replacing any position with N always gives ≥ matches than the original pattern.
    /// </summary>
    [Test]
    public void MR2b_DegenerateMotif_NAtAnyPosition_SupersetOfOriginal()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        string pattern = "ACG";

        var originalMatches = MotifFinder.FindDegenerateMotif(sequence, pattern).ToList();
        var originalPositions = originalMatches.Select(m => m.Position).ToHashSet();

        // Replace each position with N and verify superset
        for (int i = 0; i < pattern.Length; i++)
        {
            var chars = pattern.ToCharArray();
            chars[i] = 'N';
            var nPattern = new string(chars);

            var nMatches = MotifFinder.FindDegenerateMotif(sequence, nPattern).ToList();
            var nPositions = nMatches.Select(m => m.Position).ToHashSet();

            nPositions.Should().BeSubsetOf(nPositions, // self but check superset:
                because: $"pattern with N at pos {i} must match everywhere original did");

            foreach (int pos in originalPositions)
            {
                nPositions.Should().Contain(pos,
                    because: $"N at position {i}: original match at {pos} must be preserved");
            }
        }
    }

    #endregion

    #region MR3: IUPAC — prepending flanking DNA shifts positions, preserves count

    /// <summary>
    /// MR3: If we prepend a flanking sequence F to sequence S, then every match position p
    /// in S appears at position p + |F| in (F + S), and the total count is ≥.
    /// </summary>
    [Test]
    public void MR3_DegenerateMotif_PrependFlank_ShiftsPositionsPreservesCount()
    {
        string seq = "ACGTACGTACGT";
        string flank = "TTTT";
        string pattern = "RCG"; // R = A or G

        var original = MotifFinder.FindDegenerateMotif(new DnaSequence(seq), pattern).ToList();
        var shifted = MotifFinder.FindDegenerateMotif(new DnaSequence(flank + seq), pattern).ToList();

        shifted.Count.Should().BeGreaterThanOrEqualTo(original.Count,
            because: "prepending cannot destroy matches in the original part");

        var shiftedPositions = shifted.Select(m => m.Position).ToHashSet();
        foreach (var m in original)
        {
            shiftedPositions.Should().Contain(m.Position + flank.Length,
                because: $"match at {m.Position} must shift to {m.Position + flank.Length}");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-PWM-001 — Position Weight Matrix
    // ═══════════════════════════════════════════════════════════════════

    #region MR4: PWM — permutation invariance of training sequences

    /// <summary>
    /// MR4: The PWM is built from frequency counts, so reordering training sequences
    /// must produce an identical matrix (same consensus, same scores).
    /// </summary>
    [Test]
    public void MR4_CreatePwm_PermutedTrainingSet_SameMatrix()
    {
        var sequences = new[] { "ACGT", "ACGT", "AGGT", "ATGT" };
        var shuffled = new[] { "ATGT", "ACGT", "AGGT", "ACGT" };

        var pwm1 = MotifFinder.CreatePwm(sequences);
        var pwm2 = MotifFinder.CreatePwm(shuffled);

        pwm1.Consensus.Should().Be(pwm2.Consensus,
            because: "PWM is order-invariant: same multiset of training sequences → same consensus");
        pwm1.Length.Should().Be(pwm2.Length);
        pwm1.MaxScore.Should().BeApproximately(pwm2.MaxScore, 1e-10);
        pwm1.MinScore.Should().BeApproximately(pwm2.MinScore, 1e-10);

        // Compare all matrix entries
        for (int b = 0; b < 4; b++)
        {
            for (int i = 0; i < pwm1.Length; i++)
            {
                pwm1.Matrix[b, i].Should().BeApproximately(pwm2.Matrix[b, i], 1e-10,
                    because: $"matrix[{b},{i}] must be identical for permuted training sets");
            }
        }
    }

    #endregion

    #region MR5: PWM — duplicating training set doesn't change consensus

    /// <summary>
    /// MR5: Doubling every training sequence preserves relative frequencies,
    /// so log-odds scores (with pseudocount) should converge but consensus stays the same.
    /// </summary>
    [Test]
    public void MR5_CreatePwm_DuplicatedTrainingSet_SameConsensus()
    {
        var sequences = new[] { "ACGT", "AGGT", "ATGT" };
        var doubled = sequences.Concat(sequences).ToArray();

        var pwm1 = MotifFinder.CreatePwm(sequences);
        var pwm2 = MotifFinder.CreatePwm(doubled);

        pwm1.Consensus.Should().Be(pwm2.Consensus,
            because: "doubling all training sequences preserves relative base frequencies");
    }

    #endregion

    #region MR6: PWM scan — lower threshold ⊇ higher threshold results

    /// <summary>
    /// MR6: For any two thresholds t₁ < t₂, the set of matches at t₁ must be
    /// a superset of matches at t₂ (more permissive → more results).
    /// </summary>
    [Test]
    public void MR6_ScanWithPwm_LowerThreshold_SupersetOfHigherThreshold()
    {
        var training = new[] { "ACGT", "ACGT", "AGGT" };
        var pwm = MotifFinder.CreatePwm(training);
        var seq = new DnaSequence("ACGTACGTAGGTTTTT");

        double lowThreshold = -5.0;
        double highThreshold = 0.0;

        var lowResults = MotifFinder.ScanWithPwm(seq, pwm, lowThreshold).ToList();
        var highResults = MotifFinder.ScanWithPwm(seq, pwm, highThreshold).ToList();

        lowResults.Count.Should().BeGreaterThanOrEqualTo(highResults.Count,
            because: "lower threshold is more permissive, so it yields ≥ matches");

        var lowPositions = lowResults.Select(m => m.Position).ToHashSet();
        foreach (var m in highResults)
        {
            lowPositions.Should().Contain(m.Position,
                because: "every match above high threshold must also be above low threshold");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  PAT-APPROX-002 — Edit Distance & Approximate Matching
    // ═══════════════════════════════════════════════════════════════════

    #region MR7: Edit distance — non-negativity

    /// <summary>
    /// MR7: Edit distance is a metric: d(a, b) ≥ 0 for all strings a, b.
    /// </summary>
    [Test]
    public void MR7_EditDistance_NonNegative()
    {
        var pairs = new[]
        {
            ("ACGT", "TGCA"),
            ("", "ACGT"),
            ("ACGT", ""),
            ("A", "A"),
            ("GATTACA", "GCATGCU"),
            (RandomDna(50), RandomDna(50)),
            (RandomDna(100), RandomDna(30)),
        };

        foreach (var (a, b) in pairs)
        {
            ApproximateMatcher.EditDistance(a, b).Should().BeGreaterThanOrEqualTo(0,
                because: $"edit distance is a metric and must be non-negative for ({a.Length},{b.Length})");
        }
    }

    #endregion

    #region MR8: FindWithEdits — monotonicity wrt maxEdits

    /// <summary>
    /// MR8: Increasing maxEdits can only add matches, never remove them.
    /// FindWithEdits(seq, pat, k) ⊆ FindWithEdits(seq, pat, k+1) by position.
    /// </summary>
    [Test]
    public void MR8_FindWithEdits_HigherMaxEdits_SupersetOfLower()
    {
        string sequence = "ACGTACGTACGTAAAA";
        string pattern = "ACGT";

        var results0 = ApproximateMatcher.FindWithEdits(sequence, pattern, 0).ToList();
        var results1 = ApproximateMatcher.FindWithEdits(sequence, pattern, 1).ToList();
        var results2 = ApproximateMatcher.FindWithEdits(sequence, pattern, 2).ToList();

        // Positions should form a chain of supersets
        var pos0 = results0.Select(r => r.Position).ToHashSet();
        var pos1 = results1.Select(r => r.Position).ToHashSet();
        var pos2 = results2.Select(r => r.Position).ToHashSet();

        pos0.Should().BeSubsetOf(pos1,
            because: "every exact match (d≤0) is also a match with d≤1");
        pos1.Should().BeSubsetOf(pos2,
            because: "every match at d≤1 is also a match at d≤2");

        results1.Count.Should().BeGreaterThanOrEqualTo(results0.Count);
        results2.Count.Should().BeGreaterThanOrEqualTo(results1.Count);
    }

    #endregion

    #region MR9: FindWithMismatches — monotonicity wrt maxMismatches

    /// <summary>
    /// MR9: Increasing maxMismatches can only add matches.
    /// FindWithMismatches(seq, pat, k) ⊆ FindWithMismatches(seq, pat, k+1).
    /// </summary>
    [Test]
    public void MR9_FindWithMismatches_HigherMax_SupersetOfLower()
    {
        string sequence = "ACGTACGTACGTAAAA";
        string pattern = "ACGT";

        var pos0 = ApproximateMatcher.FindWithMismatches(sequence, pattern, 0)
            .Select(r => r.Position).ToHashSet();
        var pos1 = ApproximateMatcher.FindWithMismatches(sequence, pattern, 1)
            .Select(r => r.Position).ToHashSet();
        var pos2 = ApproximateMatcher.FindWithMismatches(sequence, pattern, 2)
            .Select(r => r.Position).ToHashSet();

        pos0.Should().BeSubsetOf(pos1,
            because: "every match with 0 mismatches is also a match with ≤1 mismatches");
        pos1.Should().BeSubsetOf(pos2,
            because: "every match with ≤1 mismatches is also a match with ≤2 mismatches");
    }

    #endregion

    #region MR10: Exact match ⊆ approximate match

    /// <summary>
    /// MR10: Every exact match of a pattern must also appear as an approximate match
    /// with distance 0 (composition relation).
    /// </summary>
    [Test]
    public void MR10_ExactMatch_IsSubsetOf_ApproximateMatch()
    {
        string sequence = "ACGTACGTACGT";
        string pattern = "ACGT";

        var exactPositions = ApproximateMatcher.FindWithMismatches(sequence, pattern, 0)
            .Select(r => r.Position).ToHashSet();
        var approxResults = ApproximateMatcher.FindWithEdits(sequence, pattern, 1).ToList();
        var approxPositions = approxResults.Select(r => r.Position).ToHashSet();

        exactPositions.Should().BeSubsetOf(approxPositions,
            because: "every exact match is by definition within edit distance 1");

        // All exact matches should have distance 0 in the approximate results
        foreach (int pos in exactPositions)
        {
            var r = approxResults.Where(r => r.Position == pos && r.Distance == 0);
            r.Should().NotBeEmpty(
                because: $"position {pos} is an exact match, so it must appear with distance=0");
        }
    }

    #endregion

    #region MR11: Edit distance — prepending same prefix to both strings doesn't change distance

    /// <summary>
    /// MR11: d(P+A, P+B) = d(A, B) when P is a common prefix prepended to both.
    /// This is NOT always true for edit distance in general (it IS true for suffix-free
    /// comparison), but adding the same prefix should yield d(PA, PB) ≤ d(A,B) + 0
    /// since optimal alignment can match P exactly. In fact d(PA, PB) = d(A, B).
    /// </summary>
    [Test]
    public void MR11_EditDistance_CommonPrefix_DoesNotChangeDistance()
    {
        var pairs = new[]
        {
            ("ACG", "ATG"),
            ("GATTACA", "GCATGCA"),
            ("AA", "TT"),
        };

        string prefix = "CCCC";

        foreach (var (a, b) in pairs)
        {
            int original = ApproximateMatcher.EditDistance(a, b);
            int withPrefix = ApproximateMatcher.EditDistance(prefix + a, prefix + b);

            withPrefix.Should().Be(original,
                because: $"prepending identical prefix '{prefix}' to both strings doesn't change distance");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-STR-001 — Microsatellite (STR) Detection
    // ═══════════════════════════════════════════════════════════════════

    #region MR12: Microsatellites — lower minRepeats ⊇ higher minRepeats

    /// <summary>
    /// MR12: Relaxing minRepeats threshold can only increase the result set.
    /// Positions found with minRepeats=5 must also be found with minRepeats=3.
    /// </summary>
    [Test]
    public void MR12_Microsatellites_LowerMinRepeats_SupersetOfHigher()
    {
        var seq = new DnaSequence("ACACACACACACACACGCATATATATAT");

        var results3 = RepeatFinder.FindMicrosatellites(seq, 1, 6, 3).ToList();
        var results5 = RepeatFinder.FindMicrosatellites(seq, 1, 6, 5).ToList();

        results3.Count.Should().BeGreaterThanOrEqualTo(results5.Count,
            because: "lower minRepeats is more permissive");

        // Every result from the stricter set must appear in the relaxed set
        var positions3 = results3.Select(r => (r.Position, r.RepeatUnit)).ToHashSet();
        foreach (var r in results5)
        {
            positions3.Should().Contain((r.Position, r.RepeatUnit),
                because: $"STR '{r.RepeatUnit}' at {r.Position} found with minRepeats=5 must also be found with minRepeats=3");
        }
    }

    #endregion

    #region MR13: Microsatellites — wider unit range ⊇ narrower

    /// <summary>
    /// MR13: Widening the unit length range [minUnit, maxUnit] can only add results.
    /// Results from [2,4] must be a subset of results from [1,6].
    /// </summary>
    [Test]
    public void MR13_Microsatellites_WiderUnitRange_SupersetOfNarrower()
    {
        var seq = new DnaSequence("AAAAAACACACACACATGATGATGATG");

        var narrow = RepeatFinder.FindMicrosatellites(seq, 2, 4, 3).ToList();
        var wide = RepeatFinder.FindMicrosatellites(seq, 1, 6, 3).ToList();

        wide.Count.Should().BeGreaterThanOrEqualTo(narrow.Count,
            because: "wider unit range [1,6] ⊇ [2,4]");

        var wideSet = wide.Select(r => (r.Position, r.RepeatUnit)).ToHashSet();
        foreach (var r in narrow)
        {
            wideSet.Should().Contain((r.Position, r.RepeatUnit),
                because: $"STR '{r.RepeatUnit}' at {r.Position} from [2,4] must appear in [1,6]");
        }
    }

    #endregion

    #region MR14: Microsatellites — doubling the repeat region increases repeat count

    /// <summary>
    /// MR14: If a sequence S contains a microsatellite with unit U repeated k times,
    /// then a sequence with unit U repeated 2k times must yield repeat count ≥ 2k.
    /// </summary>
    [Test]
    public void MR14_Microsatellites_DoubledRepeat_IncreasedCount()
    {
        string unit = "AC";
        int k = 5;
        string single = string.Concat(Enumerable.Repeat(unit, k));
        string doubled = string.Concat(Enumerable.Repeat(unit, 2 * k));
        string flank = "TTTTTTTTTT";

        var seq1 = new DnaSequence(flank + single + flank);
        var seq2 = new DnaSequence(flank + doubled + flank);

        var results1 = RepeatFinder.FindMicrosatellites(seq1, 2, 2, 3).ToList();
        var results2 = RepeatFinder.FindMicrosatellites(seq2, 2, 2, 3).ToList();

        var acResult1 = results1.FirstOrDefault(r => r.RepeatUnit == "AC");
        var acResult2 = results2.FirstOrDefault(r => r.RepeatUnit == "AC");

        acResult1.RepeatCount.Should().Be(k);
        acResult2.RepeatCount.Should().BeGreaterThanOrEqualTo(2 * k,
            because: "doubling the repeat region should double the repeat count");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-TANDEM-001 — Tandem Repeat Detection (GenomicAnalyzer)
    // ═══════════════════════════════════════════════════════════════════

    #region MR15: Tandem repeats — prepending shifts positions, preserves count

    /// <summary>
    /// MR15: Prepending a non-repeating flank F to sequence S shifts all tandem repeat
    /// positions by |F|, but the number of detected repeats in the original region is ≥.
    /// </summary>
    [Test]
    public void MR15_TandemRepeats_PrependFlank_ShiftsPositions()
    {
        string seq = "ATGATGATGATGATG";
        string flank = "CCCCCC";

        var original = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(seq), 2, 2).ToList();
        var shifted = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(flank + seq), 2, 2).ToList();

        // Every tandem repeat in original should appear shifted in the new sequence
        foreach (var r in original)
        {
            shifted.Should().Contain(s =>
                s.Unit == r.Unit &&
                s.Position == r.Position + flank.Length &&
                s.Repetitions == r.Repetitions,
                because: $"tandem '{r.Unit}' at {r.Position} should appear at {r.Position + flank.Length} after prepending flank");
        }
    }

    #endregion

    #region MR16: Tandem repeats — higher minRepetitions ⊆ lower

    /// <summary>
    /// MR16: Raising minRepetitions threshold removes results.
    /// Results(minRep=4) ⊆ Results(minRep=2) by (Unit, Position).
    /// </summary>
    [Test]
    public void MR16_TandemRepeats_HigherMinReps_SubsetOfLower()
    {
        var seq = new DnaSequence("ATGATGATGATGATGCGCGCG");

        var low = GenomicAnalyzer.FindTandemRepeats(seq, 2, 2).ToList();
        var high = GenomicAnalyzer.FindTandemRepeats(seq, 2, 4).ToList();

        high.Count.Should().BeLessThanOrEqualTo(low.Count,
            because: "stricter minRepetitions yields fewer or equal results");

        var lowSet = low.Select(r => (r.Unit, r.Position)).ToHashSet();
        foreach (var r in high)
        {
            lowSet.Should().Contain((r.Unit, r.Position),
                because: $"tandem '{r.Unit}' at {r.Position} found with minRep=4 must exist with minRep=2");
        }
    }

    #endregion

    #region MR17: Tandem repeats — higher minUnitLength ⊆ lower

    /// <summary>
    /// MR17: Raising minUnitLength removes shorter repeats.
    /// Results(minUnit=3) ⊆ Results(minUnit=2) by (Unit, Position).
    /// </summary>
    [Test]
    public void MR17_TandemRepeats_HigherMinUnit_SubsetOfLower()
    {
        var seq = new DnaSequence("ACACACACATGATGATG");

        var low = GenomicAnalyzer.FindTandemRepeats(seq, 2, 2).ToList();
        var high = GenomicAnalyzer.FindTandemRepeats(seq, 3, 2).ToList();

        high.Count.Should().BeLessThanOrEqualTo(low.Count,
            because: "higher minUnitLength is more restrictive");

        var lowSet = low.Select(r => (r.Unit, r.Position)).ToHashSet();
        foreach (var r in high)
        {
            lowSet.Should().Contain((r.Unit, r.Position),
                because: $"tandem '{r.Unit}' at {r.Position} with minUnit=3 must also appear with minUnit=2");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-INV-001 — Inverted Repeat Detection
    // ═══════════════════════════════════════════════════════════════════

    #region MR18: Inverted repeats — right arm is reverse complement of left arm

    /// <summary>
    /// MR18: By definition, for every inverted repeat result the right arm must be
    /// the reverse complement of the left arm. This is a consistency invariant.
    /// </summary>
    [Test]
    public void MR18_InvertedRepeats_RightArm_IsReverseComplementOfLeftArm()
    {
        var seq = new DnaSequence("ACGTAAATTTTTACGT");
        var results = RepeatFinder.FindInvertedRepeats(seq, 4, 50, 0).ToList();

        results.Should().NotBeEmpty("test sequence has a known inverted repeat");

        foreach (var r in results)
        {
            string expectedRevComp = DnaSequence.GetReverseComplementString(r.LeftArm);
            r.RightArm.Should().Be(expectedRevComp,
                because: $"inverted repeat right arm must be reverse complement of left arm '{r.LeftArm}'");
        }
    }

    #endregion

    #region MR19: Inverted repeats — reverse complement of sequence preserves count

    /// <summary>
    /// MR19: If sequence S has inverted repeats, its reverse complement revcomp(S)
    /// should have the same number of inverted repeats (with mirrored positions).
    /// This follows from the symmetry of Watson-Crick base pairing.
    /// </summary>
    [Test]
    public void MR19_InvertedRepeats_ReverseComplement_PreservesCount()
    {
        // Sequence with a clear inverted repeat: ACGT...loop...ACGT(revcomp)
        var seq = new DnaSequence("TTACGTAAACCCACGTTT");
        var revComp = seq.ReverseComplement();

        var original = RepeatFinder.FindInvertedRepeats(seq, 4, 50, 3).ToList();
        var reversed = RepeatFinder.FindInvertedRepeats(revComp, 4, 50, 3).ToList();

        reversed.Count.Should().Be(original.Count,
            because: "reverse complement preserves inverted repeat structures");
    }

    #endregion

    #region MR20: Inverted repeats — relaxing arm length adds results

    /// <summary>
    /// MR20: Lowering minArmLength can only add new results (shorter arms detected).
    /// Results(minArm=6) ⊆ Results(minArm=4).
    /// </summary>
    [Test]
    public void MR20_InvertedRepeats_LowerMinArm_SupersetOfHigher()
    {
        var seq = new DnaSequence("ACGTAAAAATTTTTACGTCCCCGGGGCCCC");

        var strict = RepeatFinder.FindInvertedRepeats(seq, 6, 50, 3).ToList();
        var relaxed = RepeatFinder.FindInvertedRepeats(seq, 4, 50, 3).ToList();

        relaxed.Count.Should().BeGreaterThanOrEqualTo(strict.Count,
            because: "lower minArmLength is more permissive");

        var relaxedSet = relaxed
            .Select(r => (r.LeftArmStart, r.RightArmStart, r.ArmLength)).ToHashSet();
        foreach (var r in strict)
        {
            relaxedSet.Should().Contain((r.LeftArmStart, r.RightArmStart, r.ArmLength),
                because: $"IR at ({r.LeftArmStart},{r.RightArmStart},arm={r.ArmLength}) "
                       + $"found with minArm=6 must exist with minArm=4");
        }
    }

    #endregion

    #region MR21: Inverted repeats — wider loop range ⊇ narrower

    /// <summary>
    /// MR21: Widening the loop range [minLoop, maxLoop] can only add results.
    /// Results(minLoop=3, maxLoop=20) ⊆ Results(minLoop=0, maxLoop=50).
    /// </summary>
    [Test]
    public void MR21_InvertedRepeats_WiderLoopRange_SupersetOfNarrower()
    {
        var seq = new DnaSequence("ACGTAAATTTTTACGTCCCCGGGGCCCC");

        var narrow = RepeatFinder.FindInvertedRepeats(seq, 4, 20, 3).ToList();
        var wide = RepeatFinder.FindInvertedRepeats(seq, 4, 50, 0).ToList();

        wide.Count.Should().BeGreaterThanOrEqualTo(narrow.Count,
            because: "wider loop range is more permissive");

        var wideSet = wide
            .Select(r => (r.LeftArmStart, r.RightArmStart, r.ArmLength)).ToHashSet();
        foreach (var r in narrow)
        {
            wideSet.Should().Contain((r.LeftArmStart, r.RightArmStart, r.ArmLength),
                because: "narrower loop range results must appear in wider range results");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-DIRECT-001 — Direct Repeat Detection
    // ═══════════════════════════════════════════════════════════════════

    #region MR22: Direct repeats — wider length range ⊇ narrower

    /// <summary>
    /// MR22: Widening [minLength, maxLength] can only add results.
    /// Results(min=5, max=10) ⊆ Results(min=3, max=20).
    /// </summary>
    [Test]
    public void MR22_DirectRepeats_WiderLengthRange_SupersetOfNarrower()
    {
        var seq = new DnaSequence("ACGTACGTTTTTACGTACGT");

        var narrow = RepeatFinder.FindDirectRepeats(seq, 5, 10, 1).ToList();
        var wide = RepeatFinder.FindDirectRepeats(seq, 3, 20, 1).ToList();

        wide.Count.Should().BeGreaterThanOrEqualTo(narrow.Count,
            because: "wider length range is more permissive");

        var wideSet = wide
            .Select(r => (r.FirstPosition, r.SecondPosition, r.RepeatSequence)).ToHashSet();
        foreach (var r in narrow)
        {
            wideSet.Should().Contain((r.FirstPosition, r.SecondPosition, r.RepeatSequence),
                because: $"direct repeat '{r.RepeatSequence}' at ({r.FirstPosition},{r.SecondPosition}) "
                       + $"from [5,10] must appear in [3,20]");
        }
    }

    #endregion

    #region MR23: Direct repeats — prepend flank shifts positions

    /// <summary>
    /// MR23: Prepending a non-matching flank F shifts all direct repeat positions by |F|.
    /// </summary>
    [Test]
    public void MR23_DirectRepeats_PrependFlank_ShiftsPositions()
    {
        string seq = "ACGTACGTTTTTTACGTACGT";
        string flank = "NNNNNN".Replace('N', 'T'); // TTTTTT as flank (simple, won't create new repeats with seq)

        // Use a flank that definitely won't create new repeats matching the interior
        flank = "GGGCCC";

        var original = RepeatFinder.FindDirectRepeats(new DnaSequence(seq), 5, 15, 1).ToList();
        var shifted = RepeatFinder.FindDirectRepeats(new DnaSequence(flank + seq), 5, 15, 1).ToList();

        // Every original repeat should appear shifted in the new sequence
        foreach (var r in original)
        {
            shifted.Should().Contain(s =>
                s.RepeatSequence == r.RepeatSequence &&
                s.FirstPosition == r.FirstPosition + flank.Length &&
                s.SecondPosition == r.SecondPosition + flank.Length,
                because: $"direct repeat '{r.RepeatSequence}' at ({r.FirstPosition},{r.SecondPosition}) "
                       + $"should shift to ({r.FirstPosition + flank.Length},{r.SecondPosition + flank.Length})");
        }
    }

    #endregion

    #region MR24: Direct repeats — lower minSpacing ⊇ higher

    /// <summary>
    /// MR24: Lowering minSpacing can only add results (allows closer repeats).
    /// Results(minSpacing=5) ⊆ Results(minSpacing=1).
    /// </summary>
    [Test]
    public void MR24_DirectRepeats_LowerMinSpacing_SupersetOfHigher()
    {
        var seq = new DnaSequence("ACGTACGTACGTACGTACGT");

        var strict = RepeatFinder.FindDirectRepeats(seq, 4, 10, 5).ToList();
        var relaxed = RepeatFinder.FindDirectRepeats(seq, 4, 10, 1).ToList();

        relaxed.Count.Should().BeGreaterThanOrEqualTo(strict.Count,
            because: "lower minSpacing is more permissive");

        var relaxedSet = relaxed
            .Select(r => (r.FirstPosition, r.SecondPosition, r.RepeatSequence)).ToHashSet();
        foreach (var r in strict)
        {
            relaxedSet.Should().Contain((r.FirstPosition, r.SecondPosition, r.RepeatSequence),
                because: "results with minSpacing=5 must be included in minSpacing=1");
        }
    }

    #endregion

    #region MR25: Direct repeats — symmetry of repeat pair detection

    /// <summary>
    /// MR25: For every direct repeat (pos₁, pos₂, seq), the repeat sequence at pos₁
    /// must equal the repeat sequence at pos₂. Cross-referencing both positions against
    /// the original sequence must yield identical substrings.
    /// This verifies internal consistency across all detected pairs.
    /// </summary>
    [Test]
    public void MR25_DirectRepeats_BothPositions_MatchRepeatSequence()
    {
        string seqStr = "ACGTACGTTTTACGTACGTGGGACGTACGT";
        var seq = new DnaSequence(seqStr);

        var results = RepeatFinder.FindDirectRepeats(seq, 4, 15, 1).ToList();
        results.Should().NotBeEmpty();

        foreach (var r in results)
        {
            string atFirst = seqStr.Substring(r.FirstPosition, r.Length);
            string atSecond = seqStr.Substring(r.SecondPosition, r.Length);

            atFirst.Should().Be(r.RepeatSequence,
                because: $"sequence at FirstPosition={r.FirstPosition} must match RepeatSequence");
            atSecond.Should().Be(r.RepeatSequence,
                because: $"sequence at SecondPosition={r.SecondPosition} must match RepeatSequence");
            atFirst.Should().Be(atSecond,
                because: "direct repeat means identical sequence at both positions");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  Cross-module metamorphic relations
    // ═══════════════════════════════════════════════════════════════════

    #region MR26: Exact motif count ≤ degenerate motif count for same pattern

    /// <summary>
    /// MR26: FindExactMotif with pattern P finds positions that are a subset
    /// of FindDegenerateMotif with the same pattern P (since all concrete bases
    /// are valid IUPAC codes, the degenerate matcher must find at least the same hits).
    /// </summary>
    [Test]
    public void MR26_ExactMotif_SubsetOf_DegenerateMotif_ForConcretePattern()
    {
        var seq = new DnaSequence("ACGTACGTACGT");
        string pattern = "ACG";

        var exactPositions = MotifFinder.FindExactMotif(seq, pattern).ToHashSet();
        var degeneratePositions = MotifFinder.FindDegenerateMotif(seq, pattern)
            .Select(m => m.Position).ToHashSet();

        exactPositions.Should().BeSubsetOf(degeneratePositions,
            because: "for a concrete pattern, exact and degenerate matching must agree; degenerate ⊇ exact");

        // In fact they should be equal for a pattern with no ambiguity codes
        degeneratePositions.Should().BeEquivalentTo(exactPositions,
            because: "concrete IUPAC pattern matches exactly like exact matcher");
    }

    #endregion

    #region MR27: Hamming distance ≤ Edit distance for equal-length strings

    /// <summary>
    /// MR27: Hamming distance only allows substitutions, while edit distance allows
    /// substitutions + insertions + deletions. For equal-length strings:
    /// EditDistance(a, b) ≤ HammingDistance(a, b).
    /// </summary>
    [Test]
    public void MR27_EditDistance_LessThanOrEqual_HammingDistance_ForEqualLength()
    {
        var pairs = new[]
        {
            ("ACGT", "TGCA"),
            ("AAAA", "AAAA"),
            ("ACGT", "ACGT"),
            ("GATTACA", "GCATGCA"),
            ("AACCGGTT", "TTGGCCAA"),
        };

        foreach (var (a, b) in pairs)
        {
            int hamming = ApproximateMatcher.HammingDistance(a, b);
            int edit = ApproximateMatcher.EditDistance(a, b);

            edit.Should().BeLessThanOrEqualTo(hamming,
                because: $"edit distance (allowing indels) ≤ Hamming distance for '{a}' vs '{b}'");
        }
    }

    #endregion
}
