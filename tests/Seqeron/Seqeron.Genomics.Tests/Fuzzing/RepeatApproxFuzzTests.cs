namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Repeats area — approximate (imperfect / interrupted) tandem-repeat
/// detection (REP-APPROX-001), the Tandem Repeats Finder (Benson 1999) alignment-model
/// detector <see cref="RepeatFinder.FindApproximateTandemRepeats(DnaSequence,int,int,int)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain values to a unit and asserts the
/// code NEVER fails in an undisciplined way: it must not HANG / infinite-loop, must not
/// throw an *unhandled* runtime exception (IndexOutOfRange, DivideByZero from a period of
/// 0), and must not emit OUT-OF-CONTRACT output (a reported repeat whose bounds exceed the
/// sequence, a copy number below the documented minimum, or a percentage outside [0,100]).
/// Every input must resolve to EITHER a well-defined, theory-correct result (incl. the empty
/// result) OR a documented validation exception (ArgumentNullException / ArgumentException /
/// ArgumentOutOfRangeException). The headline hazard here is the classic approximate-scan
/// infinite-loop / DivByZero trap (a zero / negative period), so EVERY test is
/// <c>[CancelAfter]</c>-guarded — a hang manifests as the timeout firing rather than a
/// non-terminating materialization.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: REP-APPROX-001 — approximate tandem-repeat detection
/// Checklist: docs/checklists/03_FUZZING.md, row 256.
/// Fuzz strategy exercised for THIS unit (docs/checklists/03_FUZZING.md §Description):
///   • BE = Boundary Exploitation — the degenerate boundaries the checklist row calls out:
///          minReps 0 (here the analogous PERIOD floor — a 0/negative period is the
///          DivByZero / infinite-loop trap), a unit/period LONGER than the sequence (no
///          repeat can fit → empty), the empty sequence, and a single-character sequence.
///   • MC = Malformed Content — non-ACGT content (all-N) fed to BOTH documented surfaces:
///          the typed DnaSequence surface (validates and REJECTS non-ACGT at construction)
///          and the raw-string surface (does NOT validate; uppercases and scans 'N' as an
///          ordinary symbol — an all-N homopolymer is therefore legal input that must
///          produce only in-contract output, never a crash).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The approximate-tandem contract under test (REP-APPROX-001, TRF 4.10.0 model)
/// ───────────────────────────────────────────────────────────────────────────
/// FindApproximateTandemRepeats examines candidate distances d ≤ maxPeriod where a k-tuple match run
/// meets Benson's sum-of-heads criterion, aligns the sequence by wraparound DP against the candidate
/// pattern, takes the majority consensus, realigns and reports TRF's statistics (score ≥ minScore,
/// ≥ 1.9 copies for consensus ≤ 50, 1.8 above 100), then applies TRF's redundancy elimination.
/// Validation is eager on both overloads: minPeriod ≥ 1, maxPeriod ≥ minPeriod, maxPeriod ≤ 2000,
/// minScore ≥ 1 (ArgumentOutOfRangeException); a null DnaSequence throws ArgumentNullException; a
/// null / empty string yields no repeats. TRF scoring matches only identical A/C/G/T, so N (raw
/// surface) never matches and an all-N input yields no repeat.
///
/// Documented invariants pinned on every positive result:
///   INV-bounds  : 0 ≤ Start and Start + SpanLength ≤ sequence.Length.
///   INV-period  : minPeriod ≤ Period ≤ maxPeriod; |Consensus| = ConsensusSize ≥ 1 (TRF: the consensus
///                 size "may differ slightly from the period size").
///   INV-copies  : CopyNumber ≥ 1.8 (TRF minimum copy number; 1.9 for consensus ≤ 50).
///   INV-percent : PercentMatches, PercentIndels, PercentA..T ∈ [0,100]; PercentMatches + PercentIndels ≤ 100;
///                 Entropy ∈ [0,2].
///   INV-score   : AlignmentScore ≥ minScore.
///
/// Every test forces enumeration (`.ToList()`) so the in-Core validation surfaces and any hang
/// would manifest as the [CancelAfter] timeout firing.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RepeatApproxFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomDna(int length, int seed)
    {
        const string bases = "ACGT";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>Asserts the full documented output contract on a single approximate-repeat result.</summary>
    private static void AssertInContract(ApproximateTandemRepeatResult r, int seqLen, int minPeriod, int maxPeriod, int minScore)
    {
        r.Start.Should().BeGreaterThanOrEqualTo(0, "INV-bounds: Start is a real index into the sequence");
        (r.Start + r.SpanLength).Should().BeLessThanOrEqualTo(seqLen,
            "INV-bounds: the reported window never extends past the end of the sequence");
        r.Period.Should().BeInRange(minPeriod, maxPeriod, "INV-period: the period stays within the searched range");
        r.ConsensusSize.Should().BeGreaterThanOrEqualTo(1, "INV-period: a consensus has at least one base");
        r.Consensus.Length.Should().Be(r.ConsensusSize, "INV-period: the consensus string has ConsensusSize bases");
        r.CopyNumber.Should().BeGreaterThanOrEqualTo(r.ConsensusSize <= 50 ? 1.9 : 1.8,
            "INV-copies: TRF minimum copy number (1.9 copies, 1.8 for large patterns)");
        r.PercentMatches.Should().BeInRange(0.0, 100.0, "INV-percent: percent matches is a percentage");
        r.PercentIndels.Should().BeInRange(0.0, 100.0, "INV-percent: percent indels is a percentage");
        (r.PercentMatches + r.PercentIndels).Should().BeLessThanOrEqualTo(100.0 + 1e-9,
            "INV-percent: matches and indels are disjoint outcomes of the adjacent-copy trials");
        (r.PercentA + r.PercentC + r.PercentG + r.PercentT).Should().BeLessThanOrEqualTo(100.0 + 1e-9,
            "INV-percent: composition over the region (N excluded from the four columns)");
        r.Entropy.Should().BeInRange(0.0, 2.0, "INV-percent: entropy of a 4-letter composition is 0..2 bits");
        r.AlignmentScore.Should().BeGreaterThanOrEqualTo(minScore, "INV-score: only repeats reaching minScore are emitted");
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-APPROX-001 — approximate tandem repeat detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region REP-APPROX-001 — approximate tandem repeat detection

    #region BE — Boundary: degenerate period floor (minReps 0 / period 0)

    /// <summary>
    /// BE: minPeriod = 0 is the row's "minReps 0" analogue and the KEY DivByZero / hang trap.
    /// A period of 0 would make every distance / copy-number computation divide by zero. The contract
    /// REJECTS minPeriod &lt; 1 with ArgumentOutOfRangeException(nameof(minPeriod)) eagerly, before any
    /// scan — never a DivideByZeroException and never a hang. Pinned on
    /// BOTH the typed and the raw-string surface; enumeration is forced so a regression to late
    /// validation would still be caught.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindApproximate_PeriodZero_ThrowsArgumentOutOfRange_NeverDividesByZero()
    {
        var typed = () => RepeatFinder.FindApproximateTandemRepeats(new DnaSequence("ATGATGATGATG"), minPeriod: 0, maxPeriod: 6).ToList();
        var raw = () => RepeatFinder.FindApproximateTandemRepeats("ATGATGATGATG", minPeriod: 0, maxPeriod: 6).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>(
                "a period of 0 would divide by zero in the copy count; the contract rejects minPeriod < 1")
            .Which.ParamName.Should().Be("minPeriod");
        raw.Should().Throw<ArgumentOutOfRangeException>(
                "the raw-string surface enforces the same minPeriod >= 1 floor")
            .Which.ParamName.Should().Be("minPeriod");
    }

    /// <summary>
    /// BE: a NEGATIVE minPeriod is nonsensical and must be rejected just like 0 — pinning the
    /// rejection boundary is at minPeriod &lt; 1, not merely at == 0. Both surfaces throw
    /// ArgumentOutOfRangeException(nameof(minPeriod)).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindApproximate_NegativePeriod_ThrowsArgumentOutOfRange()
    {
        var typed = () => RepeatFinder.FindApproximateTandemRepeats(new DnaSequence("ATGATGATGATG"), minPeriod: -3, maxPeriod: 6).ToList();
        var raw = () => RepeatFinder.FindApproximateTandemRepeats("ATGATGATGATG", minPeriod: -3, maxPeriod: 6).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>("a negative period is below the documented floor of 1")
            .Which.ParamName.Should().Be("minPeriod");
        raw.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("minPeriod");
    }

    /// <summary>
    /// BE: maxPeriod &lt; minPeriod is an inverted range — there is no valid period to search. The
    /// contract REJECTS it with ArgumentOutOfRangeException(nameof(maxPeriod)) rather than silently scanning an empty range, on both surfaces.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindApproximate_InvertedPeriodRange_ThrowsArgumentOutOfRange()
    {
        var typed = () => RepeatFinder.FindApproximateTandemRepeats(new DnaSequence("ATGATGATGATG"), minPeriod: 6, maxPeriod: 2).ToList();
        var raw = () => RepeatFinder.FindApproximateTandemRepeats("ATGATGATGATG", minPeriod: 6, maxPeriod: 2).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>("maxPeriod < minPeriod is an empty/inverted range")
            .Which.ParamName.Should().Be("maxPeriod");
        raw.Should().Throw<ArgumentOutOfRangeException>().Which.ParamName.Should().Be("maxPeriod");
    }

    #endregion

    #region BE — Boundary: period / unit longer than the sequence

    /// <summary>
    /// BE: a period LONGER than the sequence cannot hold even the two contiguous copies a tandem
    /// repeat requires (no distance d can have a k-tuple match), so the result is cleanly EMPTY, never an
    /// out-of-range index. Here a 5-base sequence is searched with periods 6..10, all of
    /// which exceed the length. Pinned on both surfaces.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindApproximate_PeriodLongerThanSequence_IsEmptyAndDoesNotThrow()
    {
        var typed = () => RepeatFinder.FindApproximateTandemRepeats(new DnaSequence("ACGTA"), minPeriod: 6, maxPeriod: 10).ToList();
        var raw = () => RepeatFinder.FindApproximateTandemRepeats("ACGTA", minPeriod: 6, maxPeriod: 10).ToList();

        typed.Should().NotThrow("an oversized period makes the scan bound false; the loop simply never runs");
        raw.Should().NotThrow();

        RepeatFinder.FindApproximateTandemRepeats(new DnaSequence("ACGTA"), 6, 10).Should().BeEmpty(
            "no period larger than the sequence can fit two contiguous copies; the result is empty, not a crash");
        RepeatFinder.FindApproximateTandemRepeats("ACGTA", 6, 10).Should().BeEmpty();
    }

    /// <summary>
    /// BE: the exact fitting boundary. A period equal to HALF the sequence length is the largest
    /// period for which exactly two copies fit; one base longer fits zero copies. With "ATATATAT"
    /// (8 bases) a period of 4 fits exactly two copies and a period of 5 fits none — searching only
    /// period 5 must yield empty without crashing, pinning the off-by-one at the fitting edge.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindApproximate_PeriodJustOverHalfLength_IsEmpty()
    {
        RepeatFinder.FindApproximateTandemRepeats("ATATATAT", minPeriod: 5, maxPeriod: 5).Should().BeEmpty(
            "a period of 5 needs 10 bases for two copies but only 8 are present; nothing fits");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The typed surface materialises an empty
    /// DnaSequence and the scan has no position to visit; the raw surface short-circuits null/empty to
    /// the empty result (after eager parameter validation). Neither path divides, indexes, or hangs. Pinned for the default and a minimal
    /// period range.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindApproximate_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var typed = () => RepeatFinder.FindApproximateTandemRepeats(new DnaSequence(string.Empty), 1, 6).ToList();
        var rawEmpty = () => RepeatFinder.FindApproximateTandemRepeats(string.Empty, 1, 6).ToList();
        var rawNull = () => RepeatFinder.FindApproximateTandemRepeats((string)null!, 1, 6).ToList();

        typed.Should().NotThrow("an empty sequence has no region long enough to hold a tandem repeat");
        rawEmpty.Should().NotThrow("the raw-string surface short-circuits empty input to an empty result");
        rawNull.Should().NotThrow("the raw-string surface treats null input as empty, not as an error");

        RepeatFinder.FindApproximateTandemRepeats(new DnaSequence(string.Empty), 1, 6).Should().BeEmpty();
        RepeatFinder.FindApproximateTandemRepeats(string.Empty, 1, 6).Should().BeEmpty();
        RepeatFinder.FindApproximateTandemRepeats((string)null!, 1, 6).Should().BeEmpty();
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". The typed overload guards it
    /// with an explicit ArgumentNullException (ThrowIfNull) raised eagerly
    /// at the call — never a NullReferenceException.
    /// </summary>
    [Test]
    public void FindApproximate_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => RepeatFinder.FindApproximateTandemRepeats((DnaSequence)null!, 1, 6);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload null-guards its sequence; null is rejected, never dereferenced");
    }

    #endregion

    #region BE — Boundary: single-character sequence

    /// <summary>
    /// BE: a single-character sequence cannot hold a tandem repeat — a tandem needs ≥ 2 copies, and
    /// one base is shorter than even the minimal period-1 ×2 repeat (which needs 2 bases). The
    /// scan needs a match at some distance d ≥ 1 between two positions, which one base cannot provide.
    /// The detector returns empty with no crash and no hang, on both surfaces.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindApproximate_SingleCharSequence_IsEmptyAndDoesNotThrow()
    {
        var typed = () => RepeatFinder.FindApproximateTandemRepeats(new DnaSequence("A"), 1, 6).ToList();
        var raw = () => RepeatFinder.FindApproximateTandemRepeats("A", 1, 6).ToList();

        typed.Should().NotThrow("a single base cannot hold two consecutive copies of any period");
        raw.Should().NotThrow();

        RepeatFinder.FindApproximateTandemRepeats(new DnaSequence("A"), 1, 6).Should().BeEmpty(
            "one base is too short for even a period-1 unit repeated twice");
        RepeatFinder.FindApproximateTandemRepeats("A", 1, 6).Should().BeEmpty();
    }

    #endregion

    #region MC — Malformed Content: all-N (non-ACGT)

    /// <summary>
    /// MC: all-N input on the TYPED surface. The DnaSequence constructor validates its content and
    /// REJECTS any non-ACGT base with ArgumentException ("Invalid nucleotide 'N'…",
    /// DnaSequence.cs lines 112–124), so an all-N sequence never even reaches the approximate scan —
    /// it is rejected at construction, a documented validation exception, never an undisciplined
    /// crash inside the scan.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindApproximate_AllN_TypedSurface_RejectedAtConstruction()
    {
        var construct = () => new DnaSequence("NNNNNNNNNNNN");

        construct.Should().Throw<ArgumentException>(
            "the typed surface validates nucleotides; an all-N sequence is rejected before any approximate scan");
    }

    /// <summary>
    /// MC: all-N input on the RAW-string surface. This surface does NOT validate nucleotides, but TRF's
    /// scoring matrix gives +2 only to identical A/C/G/T pairs ("changed to use Similarity Matrix to avoid
    /// N matching itself", TRF source) and k-tuples containing N are never formed — so an all-N run is
    /// never a repeat (compiled TRF 4.10.0 reports nothing). The scan must complete promptly and return
    /// no result.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void FindApproximate_AllN_RawSurface_CompletesWithInContractOutput()
    {
        const string allN = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"; // 30 'N'
        const int minScore = RepeatFinder.DefaultApproximateMinScore;

        var act = () => RepeatFinder.FindApproximateTandemRepeats(allN, 1, 6, minScore).ToList();
        act.Should().NotThrow("the raw surface accepts N; an all-N run never crashes the scan");

        RepeatFinder.FindApproximateTandemRepeats(allN, 1, 6, minScore).Should().BeEmpty(
            "N never matches under TRF scoring, so an all-N run holds no tandem repeat");
        RepeatFinder.FindApproximateTandemRepeats(allN, 1, 6, minScore: 1).Should().BeEmpty(
            "even the lowest score threshold reports nothing on all-N input");
    }

    #endregion

    #region Positive sanity — a planted perfect array and a planted single mismatch

    /// <summary>
    /// Positive sanity: a PERFECT tandem array is found with the correct period and consensus, and
    /// fully in-contract. "(ATG)×10" = 30 bases of a perfect period-3 repeat. A perfect alignment
    /// scores match-weight (+2) per base = 2·30 = 60 ≥ the default minScore 50, so it is reported.
    /// We pin: a result with Period 3 and Consensus "ATG" exists, its CopyNumber ≈ 10, its score is
    /// at least 50, and every result is in-contract. This is the core-function anchor that the
    /// degenerate-boundary probes must not silently break.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void FindApproximate_PerfectTrinucleotideArray_FoundWithCorrectPeriodAndConsensus()
    {
        string array = string.Concat(Enumerable.Repeat("ATG", 10)); // 30 bp, perfect period-3
        const int minScore = RepeatFinder.DefaultApproximateMinScore;

        var results = RepeatFinder.FindApproximateTandemRepeats(new DnaSequence(array), 1, 6, minScore).ToList();

        results.Should().NotBeEmpty("a perfect 10-copy period-3 array scores 2·30 = 60 ≥ the default minScore 50");
        foreach (var r in results)
            AssertInContract(r, array.Length, 1, 6, minScore);

        var atg = results.FirstOrDefault(r => r.Period == 3);
        atg.Period.Should().Be(3, "the period-3 interpretation of a perfect (ATG)n array must be reported");
        atg.Consensus.Should().Be("ATG", "the majority-rule consensus of a perfect (ATG)n array is exactly 'ATG'");
        atg.CopyNumber.Should().BeApproximately(10.0, 0.5, "ten contiguous copies span the 30-base array");
        atg.PercentMatches.Should().BeApproximately(100.0, 0.0001, "a perfect array aligns with 100% matches");
        atg.PercentIndels.Should().BeApproximately(0.0, 0.0001, "a perfect array has no indels");
    }

    /// <summary>
    /// Positive sanity / identity threshold: a planted SINGLE mismatch in an otherwise-perfect array
    /// is STILL found when the alignment score reaches the threshold, and is NOT found when the
    /// threshold is raised above the achievable score — re-derived from the TRF scoring (match +2,
    /// mismatch −7), not hardcoded.
    /// "(ATG)×12" with one base flipped (one column mismatched) over 36 bases: a perfect tiling would
    /// score 2·36 = 72; flipping one base turns one +2 into a −7, costing 9, so the best achievable
    /// alignment score is ~63 (still well above the default 50). We pin: with minScore 50 the repeat
    /// IS reported (one mismatch is within tolerance and its score ≥ 50); with an unreachable
    /// minScore of 1000 NOTHING is reported (below threshold → suppressed). Both runs are in-contract.
    /// </summary>
    [Test]
    [CancelAfter(20000)]
    public void FindApproximate_SingleMismatch_FoundWithinThreshold_RejectedAboveThreshold()
    {
        var chars = string.Concat(Enumerable.Repeat("ATG", 12)).ToCharArray(); // 36 bp
        chars[16] = chars[16] == 'A' ? 'C' : 'A'; // flip one interior base → exactly one mismatch column
        string oneMismatch = new string(chars);

        // Within tolerance: one mismatch costs 9 off a perfect 72 → best ≈ 63 ≥ 50, so it is reported.
        var found = RepeatFinder.FindApproximateTandemRepeats(oneMismatch, 1, 6, minScore: 50).ToList();
        found.Should().NotBeEmpty("a single mismatch keeps the score well above the default minScore of 50");
        foreach (var r in found)
            AssertInContract(r, oneMismatch.Length, 1, 6, 50);
        found.Should().Contain(r => r.Period == 3 && r.PercentMatches < 100.0,
            "the imperfect period-3 array is detected with below-100% identity (an approximate repeat)");

        // Above any achievable score: an unreachable minScore suppresses every candidate.
        var none = RepeatFinder.FindApproximateTandemRepeats(oneMismatch, 1, 6, minScore: 1000).ToList();
        none.Should().BeEmpty(
            "no alignment over 36 bases can reach a score of 1000; below threshold → nothing reported, INV-score holds");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and produce ONLY
    /// in-contract results — no out-of-range bounds, no sub-minimum copy number, no out-of-range
    /// percentage, no hang — so the degenerate-boundary guards never corrupt the scan on ordinary input.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindApproximate_RandomSequence_ProducesOnlyInContractResults()
    {
        const int minScore = RepeatFinder.DefaultApproximateMinScore;
        string seq = RandomDna(200, seed: 256_001);

        var results = RepeatFinder.FindApproximateTandemRepeats(new DnaSequence(seq), 1, 6, minScore).ToList();

        foreach (var r in results)
            AssertInContract(r, seq.Length, 1, 6, minScore);
    }

    /// <summary>
    /// RB/MC: seeded random sequences with planted imperfect repeats, N and lowercase, across period
    /// ranges up to 100 — every result is in-contract, the output is ordered by (Start, end, Period),
    /// and the string and case-folded inputs agree.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void FindApproximate_RandomPlantedRepeatsWithNoise_AlwaysInContract()
    {
        var rng = new Random(20260930);
        for (int t = 0; t < 60; t++)
        {
            var sb = new System.Text.StringBuilder(RandomDna(rng.Next(0, 60), seed: t));
            int period = rng.Next(1, 40);
            string unit = RandomDna(period, seed: 1000 + t);
            int copies = rng.Next(2, 10);
            for (int c = 0; c < copies; c++)
                foreach (char ch in unit)
                {
                    double x = rng.NextDouble();
                    if (x < 0.05) continue;                                 // deletion
                    sb.Append(x < 0.12 ? "ACGTN"[rng.Next(5)] : ch);        // substitution / N
                    if (x > 0.97) sb.Append("ACGT"[rng.Next(4)]);           // insertion
                }
            sb.Append(RandomDna(rng.Next(0, 60), seed: 5000 + t));
            string seq = sb.ToString();
            int maxPeriod = rng.Next(1, 101);
            int minScore = rng.Next(1, 80);

            var results = RepeatFinder.FindApproximateTandemRepeats(seq, 1, maxPeriod, minScore).ToList();
            foreach (var r in results)
                AssertInContract(r, seq.Length, 1, maxPeriod, minScore);
            results.Should().BeInAscendingOrder(r => r.Start, "results are ordered by start position");
            RepeatFinder.FindApproximateTandemRepeats(seq.ToLowerInvariant(), 1, maxPeriod, minScore)
                .Should().Equal(results, "case-insensitive");
        }
    }

    /// <summary>
    /// Complexity guard: 100 kb of random DNA with TRF's recommended MaxPeriod 500 completes quickly
    /// (the k-tuple trigger is O(n·maxPeriod); WDP runs only on candidates), with in-contract output.
    /// </summary>
    [Test]
    [CancelAfter(60000)]
    public void FindApproximate_100kbRandom_MaxPeriod500_CompletesInContract()
    {
        string seq = RandomDna(100_000, seed: 424_242);

        var results = RepeatFinder.FindApproximateTandemRepeats(seq, 1, 500).ToList();

        foreach (var r in results)
            AssertInContract(r, seq.Length, 1, 500, RepeatFinder.DefaultApproximateMinScore);
    }

    #endregion

    #endregion
}
