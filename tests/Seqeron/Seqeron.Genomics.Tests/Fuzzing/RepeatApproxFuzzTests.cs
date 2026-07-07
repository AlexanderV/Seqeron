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
/// The approximate-tandem contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// FindApproximateTandemRepeats enumerates each period p in [minPeriod, maxPeriod], and at
/// each start grows a tandem window, derives the majority-rule consensus, globally aligns the
/// window against a whole number of tandem copies of the consensus with TRF scoring
/// (match +2, mismatch −7, indel −7; RepeatFinder.cs lines 280–303), and reports the repeat
/// when its alignment score reaches <c>minScore</c> (default
/// <see cref="RepeatFinder.DefaultApproximateMinScore"/> = 50, per Benson 1999). Benson's
/// definition — "two or more contiguous, approximate copies of a pattern" — fixes the
/// MINIMUM at two copies: the start loop requires <c>start + period·2 ≤ seq.Length</c> and
/// the window starts at <c>spanLen = period·2</c> (RepeatFinder.cs lines 367–368, 412), so
/// nothing with fewer than two copies (CopyNumber ≥ 2) is ever reported. The result record
/// (ApproximateTandemRepeatResult, RepeatFinder.cs line 1059) carries Start, SpanLength,
/// Period, ConsensusSize, Consensus, CopyNumber, PercentMatches, PercentIndels and
/// AlignmentScore.
///
/// Documented parameter / boundary contract (RepeatFinder.cs lines 320–391):
///   • sequence == null            → the typed overload throws ArgumentNullException
///     (ThrowIfNull, line 326) BEFORE touching the scan — never a NullReferenceException.
///   • minPeriod &lt; 1            → ArgumentOutOfRangeException(nameof(minPeriod)) (line 353).
///     THIS is the row's "minReps 0" analogue: a period of 0 would make the inner window
///     loop step `for (spanLen = period·2 = 0; …; spanLen++)` start at 0 and the copy-count
///     `copies = (spanLen + period − 1) / period` divide by zero — the classic DivByZero /
///     non-terminating trap. The contract REJECTS minPeriod &lt; 1 so the scan can never
///     reach a zero period. A negative period is likewise rejected.
///   • maxPeriod &lt; minPeriod    → ArgumentOutOfRangeException(nameof(maxPeriod)) (line 354).
///   • period &gt; sequence length → the start-loop bound `start + period·2 ≤ seq.Length`
///     is false at every start, so no window is ever grown for that period → empty result,
///     no out-of-range Substring (RepeatFinder.cs line 368).
///   • empty sequence             → the typed surface materialises an empty DnaSequence and
///     the Core short-circuits `IsNullOrEmpty(seq)` to the empty list (line 357); the raw
///     surface short-circuits null/empty to the empty enumerable (line 341). No division,
///     no indexing, no hang.
///   • single character           → `start + period·2 ≤ 1` is false for every period ≥ 1, so
///     no window is grown → empty result.
///   • all-N (non-ACGT)           → the TYPED surface rejects it at DnaSequence construction
///     (ArgumentException "Invalid nucleotide", DnaSequence.cs lines 112–124); the RAW
///     surface does NOT validate — it uppercases and treats 'N' as an ordinary symbol, so
///     "NNNNNN…" is a legal homopolymer-like input that may legitimately score a repeat, and
///     every such result must still be in-contract (bounds, copy number, percentages).
///
/// Documented invariants pinned on every positive result (RepeatFinder.cs lines 320–535,
/// 1059–1068; Benson 1999):
///   INV-bounds  : 0 ≤ Start and Start + SpanLength ≤ sequence.Length (the window is a real
///                 Substring(start, spanLen) inside the sequence).
///   INV-period  : minPeriod ≤ Period ≤ maxPeriod, ConsensusSize = Period, |Consensus| = Period.
///   INV-copies  : CopyNumber ≥ 2 (Benson's "two or more contiguous copies" minimum).
///   INV-percent : PercentMatches ∈ [0,100] and PercentIndels ∈ [0,100].
///   INV-score   : AlignmentScore ≥ minScore (only repeats reaching the threshold are emitted).
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
        r.ConsensusSize.Should().Be(r.Period, "INV-period: consensus size equals the period");
        r.Consensus.Length.Should().Be(r.Period, "INV-period: the consensus string has exactly Period bases");
        r.CopyNumber.Should().BeGreaterThanOrEqualTo(2.0,
            "INV-copies: Benson's minimum is two or more contiguous copies");
        r.PercentMatches.Should().BeInRange(0.0, 100.0, "INV-percent: percent matches is a percentage");
        r.PercentIndels.Should().BeInRange(0.0, 100.0, "INV-percent: percent indels is a percentage");
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
    /// A period of 0 would make the window-grow loop start at `spanLen = period·2 = 0` and the
    /// copy count `copies = (spanLen + period − 1) / period` divide by zero. The contract REJECTS
    /// minPeriod &lt; 1 with ArgumentOutOfRangeException(nameof(minPeriod)) (RepeatFinder.cs line
    /// 353) BEFORE any window is grown — never a DivideByZeroException and never a hang. Pinned on
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
    /// contract REJECTS it with ArgumentOutOfRangeException(nameof(maxPeriod)) (RepeatFinder.cs line
    /// 354) rather than silently scanning an empty range, on both surfaces.
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
    /// repeat requires, so the start-loop bound `start + period·2 ≤ seq.Length` is false at every
    /// start and no window is ever grown — a clean EMPTY result, never an out-of-range Substring
    /// (RepeatFinder.cs line 368). Here a 5-base sequence is searched with periods 6..10, all of
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
    /// DnaSequence and the Core short-circuits `IsNullOrEmpty(seq)` to the empty list
    /// (RepeatFinder.cs line 357); the raw surface short-circuits null/empty to the empty enumerable
    /// (line 341). Neither path divides, indexes, or hangs. Pinned for the default and a minimal
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
    /// with an explicit ArgumentNullException (ThrowIfNull, RepeatFinder.cs line 326) raised eagerly
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
    /// start-loop bound `start + period·2 ≤ 1` is false for every period ≥ 1, so no window is grown.
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
    /// MC: all-N input on the RAW-string surface. This surface does NOT validate nucleotides — it
    /// uppercases and scans 'N' as an ordinary symbol (RepeatFinder.cs line 344), so an all-N run is
    /// a legal homopolymer-like input. It is the maximal-repeat watch point: every window of any
    /// period is a perfect tandem of an all-N consensus, so a repeat MAY legitimately be reported —
    /// but the scan must complete promptly (no hang) and every result must be fully in-contract
    /// (bounds within [0,len], CopyNumber ≥ 2, percentages in [0,100], score ≥ minScore, consensus
    /// of all 'N'). We pin in-contract output, not absence of output.
    /// </summary>
    [Test]
    [CancelAfter(15000)]
    public void FindApproximate_AllN_RawSurface_CompletesWithInContractOutput()
    {
        const string allN = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"; // 30 'N'
        const int minScore = RepeatFinder.DefaultApproximateMinScore;

        var act = () => RepeatFinder.FindApproximateTandemRepeats(allN, 1, 6, minScore).ToList();
        act.Should().NotThrow("the raw surface treats 'N' as an ordinary symbol; an all-N run never crashes the scan");

        var results = RepeatFinder.FindApproximateTandemRepeats(allN, 1, 6, minScore).ToList();
        foreach (var r in results)
        {
            AssertInContract(r, allN.Length, 1, 6, minScore);
            r.Consensus.Should().MatchRegex("^N+$",
                "the majority-rule consensus of an all-N window is itself all 'N'");
        }
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
    /// mismatch −7; RepeatFinder.cs lines 280–283), not hardcoded.
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
    /// percentage, no hang — so the degenerate-boundary guards never corrupt the scan on ordinary
    /// input. Length kept modest because the detector is super-linear; a hang would trip [CancelAfter].
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

    #endregion

    #endregion
}
