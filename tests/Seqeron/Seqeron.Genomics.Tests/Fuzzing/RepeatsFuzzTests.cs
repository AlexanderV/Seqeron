using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Repeats area — short tandem repeat (microsatellite / STR)
/// detection (REP-STR-001) and general tandem repeat detection (REP-TANDEM-001).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain parameter values to a
/// unit and asserts that the code NEVER fails in an undisciplined way: no hang or
/// infinite loop, no state corruption, no nonsense output, and no *unhandled*
/// runtime exception (IndexOutOfRangeException, ArgumentOutOfRangeException
/// leaking from internal indexing, DivideByZero, OutOfMemory). Every input must
/// resolve to EITHER a well-defined, theory-correct result, OR a *documented,
/// intentional* validation exception (ArgumentException / ArgumentNullException /
/// ArgumentOutOfRangeException). A raw runtime exception, a hang, or a blow-up of
/// spurious "repeats" on garbage parameters is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: REP-STR-001 — short tandem repeat (microsatellite / STR) detection
/// Checklist: docs/checklists/03_FUZZING.md, row 13.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate parameter boundaries called out
///          in the checklist row: minRepeats = 0, minRepeats = 1, maxUnitLen >
///          seqLen, and the empty sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The STR-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A microsatellite is a motif U of length m (1 ≤ m ≤ 6 by default) repeated
/// consecutively k times: S[p .. p + k·m) = Uᵏ. Detection enumerates unit lengths
/// minUnitLength..maxUnitLength, skips redundant motifs (motifs that are themselves
/// a repetition of a smaller motif, e.g. "ATAT"), counts consecutive copies, and
/// emits a MicrosatelliteResult only when the count reaches minRepeats.
///   — docs/algorithms/Repeat_Analysis/Microsatellite_Detection.md §2.2, §2.4.
///
/// API entry: RepeatFinder.FindMicrosatellites(...)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs).
///
/// Documented parameter contract (Microsatellite_Detection.md §3.1, §3.3, §6.1;
/// RepeatFinder.cs lines 25–99) — the SAME validation gate is wired on ALL four
/// overloads (typed DnaSequence, raw string, and both cancellable variants):
///   • minUnitLength &lt; 1            → ArgumentOutOfRangeException;
///   • maxUnitLength &lt; minUnitLength → ArgumentOutOfRangeException;
///   • minRepeats &lt; 2              → ArgumentOutOfRangeException. This is THE
///     boundary the checklist row probes: a *literal* minRepeats = 0 would make
///     every single position a "repeat" (k = 1 ≥ 0) and explode the output, and
///     minRepeats = 1 would report every substring as a trivial 1-copy repeat —
///     so the contract REJECTS both with the documented exception rather than
///     emitting nonsense. The accepted floor is minRepeats = 2
///     (Microsatellite_Detection.md §6.1). We pin that 0 and 1 throw, and that 2
///     is accepted, so the floor cannot silently drift.
///
/// Two documented surfaces with DIFFERENT null/empty handling — fuzzing pins both
/// and the boundary between them (Microsatellite_Detection.md §3.1, §3.3):
///   (1) The TYPED overload — FindMicrosatellites(DnaSequence, …)
///       (RepeatFinder.cs lines 25–37): a null DnaSequence → ArgumentNullException
///       (explicit ThrowIfNull guard, never a NullReferenceException). An empty
///       DnaSequence (the ctor's IsNullOrEmpty short-circuit materializes one) →
///       the empty enumerable: the scan loop bound `i ≤ len − unitLen·minRepeats`
///       is negative, so the loop never runs (Microsatellite_Detection.md §6.1).
///   (2) The RAW-STRING overload — FindMicrosatellites(string, …)
///       (RepeatFinder.cs lines 65–80): null OR empty input → the empty enumerable
///       via `yield break` — NOT an exception. Non-empty input is upper-cased
///       before scanning. Parameter validation (the three OutOfRange checks) runs
///       BEFORE the null/empty short-circuit, so a bad parameter throws even on a
///       null/empty sequence.
///
/// The maxUnitLength &gt; seqLen boundary: when the largest searched unit cannot
/// even fit minRepeats copies inside the sequence, the inner loop bound
/// `i ≤ len − unitLen·minRepeats` is negative for that unit length, so the loop
/// body never executes — no Substring out-of-range, no crash, an empty (or shorter)
/// result. We pin that an oversized maxUnitLength can never index past the end.
///
/// Note these overloads return a LAZY IEnumerable (iterator with `yield`); the
/// OutOfRange parameter validation on the string overloads is itself inside the
/// iterator body, so it only fires on enumeration. Every test therefore forces
/// enumeration (`.ToList()`) so the documented exception actually surfaces and any
/// hang would manifest as a non-terminating materialization.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: REP-TANDEM-001 — general tandem repeat detection
/// Checklist: docs/checklists/03_FUZZING.md, row 14.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: minRepetitions = 0, minUnitLength = 0, a unit length
///          ceiling of 1 (only homopolymer runs count), the empty sequence, and a
///          single-character sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The tandem-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A tandem repeat is a region S[p .. p + k·|U|) = Uᵏ where a unit U of length m
/// occurs consecutively k ≥ 2 times with no gap (Tandem_Repeat_Detection.md §2.1,
/// §2.2). The CANONICAL exact detector is
///   GenomicAnalyzer.FindTandemRepeats(DnaSequence sequence,
///                                     int minUnitLength = 2,
///                                     int minRepetitions = 2)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs).
/// This is the surface the checklist row probes (minReps / minUnitLen are its
/// parameters); it is also the primary surface named in the differential checklist
/// (08_DIFFERENTIAL_TESTING.md row 14: GenomicAnalyzer vs RepeatFinder). The sibling
/// aggregator RepeatFinder.GetTandemRepeatSummary delegates to FindMicrosatellites
/// and is not the parameterised detector under fuzz here.
///
/// Documented parameter contract — BEFORE this unit's work the doc recorded that
/// FindTandemRepeats performed NO validation and that "callers can trigger
/// undefined or exception-driven behavior with nonsensical values" was an *accepted*
/// assumption (Tandem_Repeat_Detection.md §3.3, §5.4). Fuzzing proved that the two
/// degenerate boundaries the checklist row targets were UNDISCIPLINED failures, not
/// merely "exception-driven":
///   • minRepetitions = 0 → the unit-length bound `unitLen ≤ seq.Length / minReps`
///     divides by zero → a raw DivideByZeroException (a crash, not a documented
///     ArgumentException). Per the fuzzing doctrine (ADVANCED §8) a raw runtime
///     exception is a bug.
///   • minUnitLength = 0 → a zero-length unit: `unit = Substring(start, 0) = ""`,
///     and the extension `while (… && Substring(pos, 0) == "")` is ALWAYS true while
///     `pos += 0` never advances → a NON-TERMINATING loop (a hang). A hang is the
///     single worst fuzzing outcome.
/// Both were fixed at the source by adding the same validation gate the sibling
/// repeat finders already use (FindMicrosatellites / FindPalindromes):
///   • sequence == null            → ArgumentNullException (ThrowIfNull);
///   • minUnitLength &lt; 1         → ArgumentOutOfRangeException (a 0-length unit
///     would hang; negative is nonsense);
///   • minRepetitions &lt; 2        → ArgumentOutOfRangeException (k ≥ 2 is the
///     definition of a tandem repeat; 0 would divide-by-zero, 1 would mark every
///     substring a trivial 1-copy "repeat").
/// The validation is hoisted into an eager wrapper (FindTandemRepeatsCore holds the
/// `yield` body) so the exception surfaces at the call, not only on enumeration.
/// The tests below PIN these contracts so the floors cannot silently drift, and
/// pin that the homopolymer-only (unitLen ceiling 1), empty, and single-char
/// boundaries return clean empty/correct results with no crash and no hang.
///
/// Documented invariants pinned on every positive result (Tandem_Repeat_Detection.md
/// §2.4): INV-01 Repetitions ≥ minRepetitions and |Unit| ≥ minUnitLength;
/// INV-02 TotalLength = |Unit| × Repetitions; INV-03 Position + |Unit|×Repetitions
/// ≤ sequence.Length.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: REP-INV-001 — inverted repeat (stem-loop / hairpin) detection
/// Checklist: docs/checklists/03_FUZZING.md, row 15.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: minLen (= minArmLength) = 0, minLen &gt; seqLen/2, the empty
///          sequence, and a sequence with "no complement possibilities" (a homopolymer
///          whose reverse complement never appears downstream).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The inverted-repeat-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// An inverted repeat is a left arm L followed downstream — optionally across a loop
/// X — by a right arm R = ReverseComplement(L); with loop length 0 it is exactly a
/// palindrome, and with a loop it is a stem-loop / hairpin (Inverted_Repeat_Detection.md
/// §2.1, §2.2). The canonical exact detector is
///   RepeatFinder.FindInvertedRepeats(sequence,
///                                    int minArmLength = 4,
///                                    int maxLoopLength = 50,
///                                    int minLoopLength = 3)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs lines 275–352)
/// in two overloads — typed DnaSequence and raw string. The checklist's "minLen" is
/// THIS detector's minArmLength (the minimum length of each stem arm).
///
/// Documented parameter contract (Inverted_Repeat_Detection.md §3.1, §3.3; RepeatFinder.cs
/// lines 281–286, 297–302):
///   • sequence == null      → the typed overload throws ArgumentNullException
///     (ThrowIfNull); the raw-string overload treats null/empty as the empty result.
///   • minArmLength &lt; 2    → ArgumentOutOfRangeException. This is THE boundary the
///     checklist row probes: a literal minArmLength = 0 means a zero-length arm
///     (`Substring(i, 0) = ""`), whose reverse complement is also "", so EVERY
///     downstream gap position would "match" the empty arm and the detector would emit
///     a blow-up of nonsense zero-length-arm "repeats". The contract REJECTS &lt; 2 so a
///     repeat must have at least a 2-bp stem to be reported.
///   • minLoopLength &lt; 0   → ArgumentOutOfRangeException.
/// BEFORE this unit's work the raw-string overload did NOT replicate the numeric
/// validation of the typed overload (recorded as accepted Deviation #1,
/// Inverted_Repeat_Detection.md §5.4): a degenerate minArmLength = 0 fed through the
/// raw-string surface (the surface the MCP `find_inverted_repeats` tool forwards raw
/// user input to, AnalysisTools.cs line 351) emitted spurious empty-arm results instead
/// of being rejected. Fuzzing this row PROVED that asymmetry was an undisciplined failure
/// (nonsense output on a degenerate parameter, ADVANCED §8). It was fixed at the source
/// by mirroring the typed overload's two guards onto the raw-string overload
/// (`minArmLength < 2` and `minLoopLength < 0` now both throw on BOTH surfaces). The
/// tests below PIN these floors on both surfaces so they cannot silently drift.
///
/// The minLen &gt; seqLen/2 boundary: two arms of length minArmLength plus the minimum
/// loop cannot fit when `seq.Length &lt; 2·minArmLength + minLoopLength`, so the outer
/// scan bound `i ≤ seq.Length − 2·minArmLength − minLoopLength` is negative and the loop
/// body never runs — an empty result, never an out-of-range Substring
/// (Inverted_Repeat_Detection.md §6.1). "No complement possibilities" (e.g. an all-A
/// homopolymer, whose reverse complement TTTT… never appears in the same run) likewise
/// yields nothing without crashing. Every test forces enumeration (`.ToList()`) so the
/// in-iterator validation surfaces and any hang would manifest as a non-terminating
/// materialization.
///
/// Documented invariants pinned on positive results (Inverted_Repeat_Detection.md §2.4):
/// INV-01 ReverseComplement(LeftArm) = RightArm; INV-02 TotalLength = 2·ArmLength +
/// LoopLength; INV-03 LoopLength = RightArmStart − (LeftArmStart + ArmLength);
/// INV-04 CanFormHairpin ⇔ LoopLength ≥ 3.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: REP-DIRECT-001 — direct repeat (same-orientation recurring subsequence) detection
/// Checklist: docs/checklists/03_FUZZING.md, row 16.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: minLen (= minLength) = 0, minLen = 1, the empty sequence,
///          an all-unique-characters sequence (no recurrence), and an all-same-character
///          homopolymer (maximal recurrence — the combinatorial blow-up watch point).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The direct-repeat-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A direct repeat is a subsequence that recurs VERBATIM downstream in the SAME 5'→3'
/// orientation (unlike an inverted repeat, the downstream copy is the literal sequence,
/// not its reverse complement), optionally across a spacer of zero or more bases
/// (Direct_Repeat_Detection.md §2.1, §2.2). For length L, first position i and second
/// position j a pair is reported when S[i..i+L) = S[j..j+L) and j &gt; i + L − 1 + minSpacing.
/// The canonical exact detector is
///   RepeatFinder.FindDirectRepeats(sequence,
///                                  int minLength = 5,
///                                  int maxLength = 50,
///                                  int minSpacing = 1)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs lines 369–433)
/// in two overloads — typed DnaSequence and raw string. The checklist's "minLen" is THIS
/// detector's minLength (the minimum repeat length tested).
///
/// Documented parameter contract (Direct_Repeat_Detection.md §3.1, §3.3; RepeatFinder.cs
/// lines 375–377, 391–392):
///   • sequence == null     → the typed overload throws ArgumentNullException (ThrowIfNull);
///     the raw-string overload treats null/empty as the empty result (yield break).
///   • minLength &lt; 2       → ArgumentOutOfRangeException. This is THE boundary the
///     checklist row probes: a literal minLength = 0 means a zero-length candidate
///     (`Substring(i, 0) = ""`), whose suffix-tree occurrence lookup
///     (SuffixTree.FindAllOccurrences("") returns EVERY start position) matches at every
///     downstream gap → an O(n²) blow-up of spurious empty-/single-base "repeats" — the
///     nonsense-output failure mode ADVANCED §8 forbids. The contract REJECTS minLength
///     &lt; 2 so a direct repeat must be at least a 2-base subsequence. minLength = 1 is
///     likewise below the floor (every single base that recurs would be a "repeat") and is
///     rejected too.
///   • maxLength &lt; minLength → ArgumentOutOfRangeException.
/// BEFORE this unit's work the raw-string overload did NOT replicate the numeric validation
/// of the typed overload (recorded as accepted Deviation #1, Direct_Repeat_Detection.md §5.4):
/// a degenerate minLength = 0 fed through the raw-string surface (the surface the MCP
/// `find_direct_repeats` tool forwards raw user input to, AnalysisTools.cs line 368) emitted
/// the spurious empty-arm blow-up instead of being rejected. Fuzzing this row PROVED that
/// asymmetry was an undisciplined failure (nonsense O(n²) output on a degenerate parameter,
/// ADVANCED §8). It was fixed at the source by mirroring the typed overload's two numeric
/// guards onto the raw-string overload (`minLength < 2` and `maxLength < minLength` now both
/// throw on BOTH surfaces, hoisted into an eager wrapper so the exception surfaces at the call,
/// not only on enumeration). The tests below PIN these floors on both surfaces so they cannot
/// silently drift.
///
/// The all-unique-characters boundary (e.g. "ACGT"): no subsequence of length ≥ 2 recurs, so
/// the result is cleanly empty — never a crash. The all-same-character homopolymer (e.g.
/// "AAAA…") is the maximal-recurrence / combinatorial-blow-up watch point: every length-L
/// window equals every other length-L window, so with minSpacing = 0 many overlapping pairs
/// ARE legitimately reported, but the count stays polynomially bounded (deduped by the
/// (i, j, len) hash set, INV-04) and the scan completes promptly — no hang. The
/// minLength &gt; available-room boundary makes the position loop bound
/// `i ≤ len·2 + minSpacing` exceed seq.Length so the loop body never runs (empty result,
/// no out-of-range Substring). Every test forces enumeration (`.ToList()`) so the in-iterator
/// scan runs and any hang would manifest as a non-terminating materialization.
///
/// Documented invariants pinned on positive results (Direct_Repeat_Detection.md §2.4):
/// INV-01 RepeatSequence is identical at FirstPosition and SecondPosition;
/// INV-02 Spacing = SecondPosition − FirstPosition − Length;
/// INV-03 with minSpacing &gt; 0 the copies do not overlap;
/// INV-04 each (FirstPosition, SecondPosition, Length) tuple is unique.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: REP-PALIN-001 — DNA palindrome (restriction-site / self-complementary) detection
/// Checklist: docs/checklists/03_FUZZING.md, row 17.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row: minLen = 0, an ODD minLen, maxLen &gt; seqLen, and the empty
///          sequence.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The palindrome-detection contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// A DNA palindrome is a sequence equal to its own reverse complement —
/// S = ReverseComplement(S) — the biological (not textual) notion that underlies
/// Type-II restriction recognition sites such as EcoRI 'GAATTC'
/// (Palindrome_Detection.md §2.1, §2.2). Because each base must pair with a
/// complementary partner across the symmetry axis (A↔T, G↔C), a biological DNA
/// palindrome is ALWAYS EVEN length (Palindrome_Detection.md §2.2, INV-02). The
/// canonical validating detector is
///   RepeatFinder.FindPalindromes(sequence,
///                                int minLength = 4,
///                                int maxLength = 12)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/RepeatFinder.cs lines 553–608)
/// in two overloads — typed DnaSequence and raw string. The checklist's "minLen" is
/// THIS detector's minLength (the minimum palindrome window scanned).
///
/// Documented parameter contract (Palindrome_Detection.md §3.1, §3.3, §6.1;
/// RepeatFinder.cs lines 558–562, 575–578) — the SAME validation gate is wired on
/// BOTH overloads (unlike the sibling inverted/direct finders, here the raw-string
/// overload ALREADY mirrors the typed overload's numeric checks, so no source fix is
/// needed for this unit):
///   • sequence == null        → the typed overload throws ArgumentNullException
///     (ThrowIfNull); the raw-string overload treats null/empty as the empty result
///     (yield break).
///   • minLength &lt; 4 OR minLength ODD → ArgumentOutOfRangeException
///     (nameof(minLength), "Must be even and >= 4"). These are THE two boundaries the
///     checklist row probes:
///       – minLength = 0 is below the floor of 4 (a 0-length window would compare ""
///         to its empty reverse complement and "match" at every position — a nonsense
///         blow-up) → rejected by the `< 4` clause.
///       – an ODD minLength (e.g. 5) is THE distinctive target for this unit: a DNA
///         reverse-complement palindrome must have EVEN length, so an odd request is
///         biologically impossible. The validating API does NOT round or silently
///         drop it — it REJECTS it outright with ArgumentOutOfRangeException
///         (Palindrome_Detection.md §6.1, "Odd minLength → throws"). We pin that an
///         odd minLength throws on BOTH surfaces, the contract-correct behavior.
///   • maxLength &lt; minLength  → ArgumentOutOfRangeException (nameof(maxLength)).
///
/// The maxLength &gt; seqLen boundary: the inner scan bound `i ≤ seq.Length − len`
/// goes negative for every candidate length larger than the sequence, so the loop
/// body never runs for those lengths and no Substring is taken past the end — an
/// empty (or shorter) result, never an out-of-range crash, no palindrome longer than
/// the sequence (Palindrome_Detection.md §2.4 INV-03, §6.1). We pin that an oversized
/// maxLength can never index past the end and never invents an over-long palindrome.
///
/// The empty sequence: the raw-string overload short-circuits null/empty to the
/// empty enumerable via `yield break` (RepeatFinder.cs line 580–581) AFTER the numeric
/// gate, so valid lengths over empty input yield nothing; the typed overload over an
/// empty DnaSequence has a negative scan bound and likewise yields nothing — no
/// division, no indexing past the end, no hang (Palindrome_Detection.md §6.1).
///
/// Note both overloads are LAZY iterators (`yield`); the typed overload's numeric
/// validation runs eagerly inside FindPalindromes BEFORE handing off to the iterator
/// core, but the raw-string overload's checks sit INSIDE the iterator body, so they
/// only fire on enumeration. Every test therefore forces enumeration (`.ToList()`) so
/// the documented exception surfaces and any hang would manifest as a non-terminating
/// materialization.
///
/// Documented invariants pinned on positive results (Palindrome_Detection.md §2.4):
/// INV-01 Sequence = ReverseComplement(Sequence); INV-02 every reported Length is
/// EVEN; INV-03 Position is within sequence bounds; INV-04 Length equals the reported
/// palindrome's actual length.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class RepeatsFuzzTests
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

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-STR-001 — short tandem repeat detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region REP-STR-001 — short tandem repeat detection

    #region BE — Boundary: minRepeats = 0

    /// <summary>
    /// BE: minRepeats = 0 is the degenerate floor. A literal 0 would make k = 1 ≥ 0
    /// true at every position, so EVERY base would be a "microsatellite" — a
    /// nonsense blow-up of the result set. The documented contract REJECTS it with
    /// ArgumentOutOfRangeException (minRepeats &lt; 2 → throw, RepeatFinder.cs lines
    /// 34/73/95; Microsatellite_Detection.md §3.3, §6.1) on BOTH the typed and the
    /// raw-string surface, BEFORE producing any result. The string overload's check
    /// lives inside the iterator, so we force enumeration to surface it.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MinRepeatsZero_ThrowsArgumentOutOfRange()
    {
        var typed = () => RepeatFinder.FindMicrosatellites(new DnaSequence("ACACACAC"), 1, 6, 0).ToList();
        var raw = () => RepeatFinder.FindMicrosatellites("ACACACAC", 1, 6, 0).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>(
            "minRepeats = 0 is below the documented floor of 2; a literal 0 would mark every position a repeat");
        raw.Should().Throw<ArgumentOutOfRangeException>(
            "the raw-string overload enforces the same minRepeats >= 2 floor before scanning");
    }

    #endregion

    #region BE — Boundary: minRepeats = 1

    /// <summary>
    /// BE: minRepeats = 1 is the trivial-repeat boundary. With minRepeats = 1 every
    /// substring is a 1-copy "repeat", which would produce a defined-but-useless
    /// quadratic blow-up. The contract REJECTS it with ArgumentOutOfRangeException
    /// (minRepeats &lt; 2 → throw, RepeatFinder.cs lines 34/73/95;
    /// Microsatellite_Detection.md §6.1) on both surfaces — there is no hang and no
    /// blow-up because the value never reaches the scan. The accepted floor is
    /// exactly 2.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MinRepeatsOne_ThrowsArgumentOutOfRange()
    {
        var typed = () => RepeatFinder.FindMicrosatellites(new DnaSequence("ACACACAC"), 1, 6, 1).ToList();
        var raw = () => RepeatFinder.FindMicrosatellites("ACACACAC", 1, 6, 1).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>(
            "minRepeats = 1 is below the documented floor of 2; every substring would be a trivial 1-copy repeat");
        raw.Should().Throw<ArgumentOutOfRangeException>(
            "the raw-string overload enforces the same minRepeats >= 2 floor before scanning");
    }

    /// <summary>
    /// BE: the accepted floor minRepeats = 2 must NOT throw and must behave —
    /// pinning that the rejection boundary is exactly at &lt; 2, not at ≤ 2. A clean
    /// di-nucleotide with two copies ("ATAT" is redundant as a unit, so we use a
    /// non-redundant 2-mer repeated twice: "ATAT" → use "AC" twice = "ACAC") is
    /// reported with count 2.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MinRepeatsTwo_IsAcceptedAndDetects()
    {
        var act = () => RepeatFinder.FindMicrosatellites("ACAC", 2, 2, 2).ToList();
        act.Should().NotThrow("minRepeats = 2 is the accepted floor, not below it");

        var results = RepeatFinder.FindMicrosatellites("ACAC", 2, 2, 2).ToList();

        results.Should().ContainSingle(
            "exactly one non-redundant 2-mer 'AC' repeats twice across 'ACAC'");
        results[0].RepeatUnit.Should().Be("AC");
        results[0].RepeatCount.Should().Be(2);
    }

    #endregion

    #region BE — Boundary: maxUnitLength > sequence length

    /// <summary>
    /// BE: maxUnitLength far larger than the sequence length must NOT crash. For a
    /// unit length that cannot fit minRepeats copies in the sequence, the scan
    /// bound `i ≤ len − unitLen·minRepeats` is negative, so the inner loop never
    /// runs and no Substring is taken past the end (RepeatFinder.cs lines 184–195).
    /// The result is simply the repeats found by the unit lengths that DO fit —
    /// here, the short sequence "AC" has no room for any 2+ copy repeat, so the
    /// result is empty. Pinned on both the typed and raw surfaces.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MaxUnitLengthExceedsSequence_IsEmptyAndDoesNotThrow()
    {
        // 2-char sequence, maxUnitLength = 1000 ≫ length. No unit length 1..1000 can
        // fit ≥ 2 copies in 2 chars except a 1-mer homopolymer, which "AC" is not.
        var typedAct = () => RepeatFinder.FindMicrosatellites(new DnaSequence("AC"), 1, 1000, 2).ToList();
        var rawAct = () => RepeatFinder.FindMicrosatellites("AC", 1, 1000, 2).ToList();

        typedAct.Should().NotThrow(
            "an oversized maxUnitLength makes the scan bound negative for the big units; the loop simply never runs");
        rawAct.Should().NotThrow(
            "the raw-string overload is equally guarded against indexing past the sequence end");

        RepeatFinder.FindMicrosatellites(new DnaSequence("AC"), 1, 1000, 2).Should().BeEmpty(
            "no unit can repeat 2× inside a 2-char heteropolymer; the oversized max yields nothing, not a crash");
        RepeatFinder.FindMicrosatellites("AC", 1, 1000, 2).Should().BeEmpty();
    }

    /// <summary>
    /// BE: maxUnitLength &gt; seqLen on a sequence that DOES contain a short repeat
    /// must still find that short repeat and must not be derailed by the oversized
    /// upper bound. "ACACAC" has a 'AC' di-repeat (3 copies); a maxUnitLength of 100
    /// (far past the 6-char length) leaves that detection intact while the unit
    /// lengths > 2 simply contribute nothing.
    /// </summary>
    [Test]
    public void FindMicrosatellites_MaxUnitLengthExceedsSequence_StillFindsShortRepeat()
    {
        var results = RepeatFinder.FindMicrosatellites("ACACAC", 1, 100, 2).ToList();

        results.Should().Contain(r => r.RepeatUnit == "AC" && r.RepeatCount == 3,
            "the 'AC' di-repeat (3 copies) is found despite maxUnitLength ≫ sequence length");
        results.Should().OnlyContain(r => r.RepeatCount >= 2,
            "INV-01: every reported repeat satisfies RepeatCount >= minRepeats, even at the boundary");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The raw-string overload
    /// short-circuits null/empty to the empty enumerable via `yield break`
    /// (RepeatFinder.cs lines 75–76) — no exception. The typed overload over an
    /// empty DnaSequence (the ctor materializes an empty sequence from "" via its
    /// IsNullOrEmpty short-circuit) yields nothing because the scan bound is
    /// negative. Neither path divides, indexes, or hangs on empty input
    /// (Microsatellite_Detection.md §6.1).
    /// </summary>
    [Test]
    public void FindMicrosatellites_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var typedAct = () => RepeatFinder.FindMicrosatellites(new DnaSequence(string.Empty), 1, 6, 3).ToList();
        var rawEmptyAct = () => RepeatFinder.FindMicrosatellites(string.Empty, 1, 6, 3).ToList();
        var rawNullAct = () => RepeatFinder.FindMicrosatellites((string)null!, 1, 6, 3).ToList();

        typedAct.Should().NotThrow("an empty sequence has no region long enough to hold a repeat; it yields nothing");
        rawEmptyAct.Should().NotThrow("the raw-string overload short-circuits empty input to an empty result");
        rawNullAct.Should().NotThrow("the raw-string overload treats null input as empty, not as an error");

        RepeatFinder.FindMicrosatellites(new DnaSequence(string.Empty), 1, 6, 3).Should().BeEmpty();
        RepeatFinder.FindMicrosatellites(string.Empty, 1, 6, 3).Should().BeEmpty();
        RepeatFinder.FindMicrosatellites((string)null!, 1, 6, 3).Should().BeEmpty();
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". The typed
    /// overload guards it with an explicit ArgumentNullException
    /// (ArgumentNullException.ThrowIfNull, RepeatFinder.cs line 31) — never a
    /// NullReferenceException. The exception is raised eagerly at the call (it is
    /// outside the iterator body), so no enumeration is required to surface it.
    /// </summary>
    [Test]
    public void FindMicrosatellites_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => RepeatFinder.FindMicrosatellites((DnaSequence)null!, 1, 6, 3);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload null-guards its sequence; null is rejected, never dereferenced");
    }

    #endregion

    #region Positive sanity — a clear microsatellite is detected correctly

    /// <summary>
    /// Positive sanity: alongside the boundary probes, a textbook microsatellite
    /// must be detected with the CORRECT unit and copy count — fuzzing the degenerate
    /// edges must not come at the cost of the core function silently breaking.
    /// "CACACACA" is a di-nucleotide STR: unit 'CA' repeated 4 times. Pinned per
    /// INV-01..INV-04 (Microsatellite_Detection.md §2.4): RepeatCount = 4,
    /// unit length 2 → Dinucleotide, TotalLength = 2 × 4 = 8.
    /// </summary>
    [Test]
    public void FindMicrosatellites_DinucleotideStr_DetectedWithCorrectUnitAndCount()
    {
        var results = RepeatFinder.FindMicrosatellites("CACACACA", 1, 6, 3).ToList();

        var ca = results.Should().ContainSingle(r => r.RepeatUnit == "CA").Subject;
        ca.RepeatCount.Should().Be(4, "'CA' tiles 'CACACACA' exactly four times");
        ca.RepeatType.Should().Be(RepeatType.Dinucleotide, "a 2-bp unit classifies as a dinucleotide repeat");
        ca.TotalLength.Should().Be(8, "INV-03: TotalLength = unitLength (2) × RepeatCount (4)");
        ca.Position.Should().Be(0, "the repeat tract starts at the first base");
    }

    /// <summary>
    /// Positive sanity: a trinucleotide STR — the clinically important class
    /// (Huntington CAG, Fragile-X CGG, …). "(GAT)5" = "GATGATGATGATGAT" is the
    /// motif 'GAT' repeated 5 times. Pinned: unit 'GAT', count 5, Trinucleotide,
    /// TotalLength = 3 × 5 = 15 (Microsatellite_Detection.md §2.1, §2.4).
    /// </summary>
    [Test]
    public void FindMicrosatellites_TrinucleotideStr_DetectedWithCorrectUnitAndCount()
    {
        var results = RepeatFinder.FindMicrosatellites("GATGATGATGATGAT", 1, 6, 3).ToList();

        var gat = results.Should().ContainSingle(r => r.RepeatUnit == "GAT").Subject;
        gat.RepeatCount.Should().Be(5, "'GAT' tiles '(GAT)5' exactly five times");
        gat.RepeatType.Should().Be(RepeatType.Trinucleotide, "a 3-bp unit classifies as a trinucleotide repeat");
        gat.TotalLength.Should().Be(15, "INV-03: TotalLength = unitLength (3) × RepeatCount (5)");
    }

    /// <summary>
    /// Positive sanity at scale: a long but modest STR tract must complete promptly
    /// (no hang / no quadratic blow-up at the accepted floor) and report a count
    /// consistent with the tract length. 'CA' repeated 500 times → a single
    /// dinucleotide result with RepeatCount 500. Bounds kept modest to guard against
    /// hang while still exercising the extension loop.
    /// </summary>
    [Test]
    public void FindMicrosatellites_LongHomogeneousTract_CompletesAndCountsCorrectly()
    {
        const int copies = 500;
        string tract = string.Concat(Enumerable.Repeat("CA", copies)); // 1000 bases

        var results = RepeatFinder.FindMicrosatellites(tract, 2, 2, 3).ToList();

        var ca = results.Should().ContainSingle(r => r.RepeatUnit == "CA").Subject;
        ca.RepeatCount.Should().Be(copies, "the extension loop counts every consecutive 'CA' copy without overflow or hang");
        ca.TotalLength.Should().Be(copies * 2);
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must produce only
    /// well-formed results and must complete — no spurious negative counts, no
    /// out-of-range positions, no hang. Every result must satisfy the documented
    /// invariants INV-01..INV-03 (Microsatellite_Detection.md §2.4) regardless of
    /// what random content appears.
    /// </summary>
    [Test]
    public void FindMicrosatellites_RandomSequence_ProducesOnlyWellFormedResults()
    {
        string seq = RandomDna(2000, seed: 13_001);

        var results = RepeatFinder.FindMicrosatellites(seq, 1, 6, 3).ToList();

        results.Should().OnlyContain(r =>
            r.RepeatCount >= 3 &&                                   // INV-01
            r.RepeatUnit.Length >= 1 && r.RepeatUnit.Length <= 6 && // INV-02
            r.TotalLength == r.RepeatUnit.Length * r.RepeatCount &&  // INV-03
            r.Position >= 0 && r.Position + r.TotalLength <= seq.Length,
            "every result on random input is well-formed: count ≥ minRepeats, unit in range, total = unit×count, in bounds");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-TANDEM-001 — general tandem repeat detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region REP-TANDEM-001 — tandem repeat detection

    #region BE — Boundary: minRepetitions = 0

    /// <summary>
    /// BE: minRepetitions = 0 is the degenerate floor and was the KEY crash target.
    /// In the pre-fix detector the unit-length bound `unitLen ≤ seq.Length / minReps`
    /// divides by zero → a raw DivideByZeroException (a crash, not a documented
    /// validation error). A tandem repeat is defined as k ≥ 2 copies
    /// (Tandem_Repeat_Detection.md §2.2), so the contract now REJECTS 0 with
    /// ArgumentOutOfRangeException at the call site (eager validation in the wrapper),
    /// never dividing by zero. We force enumeration so a regression to lazy/late
    /// validation would still be caught.
    /// </summary>
    [Test]
    public void FindTandemRepeats_MinRepetitionsZero_ThrowsArgumentOutOfRange()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats(new DnaSequence("ATGATGATG"), 2, 0).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
                "minRepetitions = 0 would divide by zero in the unit-length bound; the contract rejects it")
            .Which.ParamName.Should().Be("minRepetitions");
    }

    #endregion

    #region BE — Boundary: minUnitLength = 0

    /// <summary>
    /// BE: minUnitLength = 0 is the degenerate zero-length-unit boundary and was the
    /// KEY HANG target. A 0-length unit makes `unit = Substring(start, 0) = ""`, and
    /// the extension loop `while (… && Substring(pos, 0) == "")` is always true while
    /// `pos += 0` never advances — a non-terminating loop (the worst fuzzing outcome).
    /// The contract now REJECTS minUnitLength &lt; 1 with ArgumentOutOfRangeException
    /// so the scan can never enter that infinite loop. Pinned eagerly; the guarded
    /// call must return promptly, not hang.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindTandemRepeats_MinUnitLengthZero_ThrowsArgumentOutOfRange()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats(new DnaSequence("ATGATGATG"), 0, 2).ToList();

        act.Should().Throw<ArgumentOutOfRangeException>(
                "a zero-length unit tiles every position infinitely and never advances the scan; the contract rejects it")
            .Which.ParamName.Should().Be("minUnitLength");
    }

    #endregion

    #region BE — Boundary: unit-length ceiling of 1 (homopolymer-only)

    /// <summary>
    /// BE: a unit-length ceiling of 1 — i.e. searching ONLY length-1 units
    /// (minUnitLength = 1) — means only homopolymer runs (AAAA…, the mononucleotide
    /// class) can ever be a tandem repeat; no di-/tri-/longer unit is even considered.
    /// "AAATGC" has a single 'A' run of 3 and nothing else repeats, so exactly one
    /// result — unit "A", 3 copies at position 0 — must be reported, and every result
    /// must be a 1-bp unit. This pins that the smallest legal unit length behaves and
    /// does not over- or under-report (no crash, no spurious multi-bp units).
    /// </summary>
    [Test]
    public void FindTandemRepeats_UnitLengthOne_FindsOnlyHomopolymerRuns()
    {
        var results = GenomicAnalyzer.FindTandemRepeats(new DnaSequence("AAATGC"), minUnitLength: 1, minRepetitions: 2).ToList();

        results.Should().ContainSingle("only the 'A' homopolymer run repeats ≥ 2× in 'AAATGC'");
        results[0].Unit.Should().Be("A");
        results[0].Repetitions.Should().Be(3, "the run 'AAA' is three consecutive 'A' copies");
        results[0].Position.Should().Be(0);
        results.Should().OnlyContain(r => r.Unit.Length == 1,
            "with a unit-length ceiling of 1 only mononucleotide (homopolymer) units are searched");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. With valid thresholds the
    /// scan bounds `unitLen ≤ 0/minReps = 0` and `start ≤ 0 − unitLen·minReps` are
    /// degenerate, so the detector yields nothing — no division, no indexing past the
    /// end, no hang (Tandem_Repeat_Detection.md §6.1). Pinned for the default and a
    /// minimal unit length.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindTandemRepeats_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats(new DnaSequence(string.Empty), 2, 2).ToList();
        act.Should().NotThrow("an empty sequence has no region long enough to hold a tandem repeat");

        GenomicAnalyzer.FindTandemRepeats(new DnaSequence(string.Empty), 2, 2).Should().BeEmpty();
        GenomicAnalyzer.FindTandemRepeats(new DnaSequence(string.Empty), 1, 2).Should().BeEmpty(
            "even searching 1-bp units, the empty sequence yields no homopolymer run");
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no input". The pre-fix detector
    /// dereferenced `sequence.Sequence` immediately → a raw NullReferenceException.
    /// The contract now guards it with ArgumentNullException (ThrowIfNull), raised
    /// eagerly at the call — never a NullReferenceException
    /// (Tandem_Repeat_Detection.md §3.3, now enforced).
    /// </summary>
    [Test]
    public void FindTandemRepeats_NullSequence_ThrowsArgumentNullException()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats((DnaSequence)null!, 2, 2);

        act.Should().Throw<ArgumentNullException>(
            "a null sequence is null-guarded, never dereferenced into a NullReferenceException");
    }

    #endregion

    #region BE — Boundary: single-character sequence

    /// <summary>
    /// BE: a single-character sequence cannot hold a tandem repeat — a tandem needs
    /// ≥ 2 copies (Tandem_Repeat_Detection.md §2.2), and one base is half of even the
    /// shortest possible 1-bp ×2 repeat. The detector must return empty with no crash
    /// and no hang, whether the unit-length floor is the default 2 or the minimal 1.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindTandemRepeats_SingleCharSequence_IsEmptyAndDoesNotThrow()
    {
        var act = () => GenomicAnalyzer.FindTandemRepeats(new DnaSequence("A"), 1, 2).ToList();
        act.Should().NotThrow("a single base cannot hold two consecutive copies of any unit");

        GenomicAnalyzer.FindTandemRepeats(new DnaSequence("A"), 1, 2).Should().BeEmpty(
            "one base is too short for even a 1-bp unit repeated twice");
        GenomicAnalyzer.FindTandemRepeats(new DnaSequence("A"), 2, 2).Should().BeEmpty();
    }

    #endregion

    #region Positive sanity — a clear tandem repeat is detected correctly

    /// <summary>
    /// Positive sanity: alongside the degenerate probes, the textbook tandem from the
    /// algorithm doc — "(ATG)3" = "ATGATGATG" — must be detected with the CORRECT unit
    /// and copy count, so the boundary hardening never silently breaks the core
    /// function. Pinned per INV-01..INV-03 (Tandem_Repeat_Detection.md §2.4): unit
    /// "ATG", 3 copies at position 0, TotalLength = 3 × 3 = 9.
    /// </summary>
    [Test]
    public void FindTandemRepeats_ClearTrinucleotideTandem_DetectedWithCorrectUnitAndCount()
    {
        var results = GenomicAnalyzer.FindTandemRepeats(new DnaSequence("ATGATGATG"), minUnitLength: 3, minRepetitions: 3).ToList();

        var atg = results.Should().ContainSingle(r => r.Unit == "ATG").Subject;
        atg.Repetitions.Should().Be(3, "'ATG' tiles '(ATG)3' exactly three times");
        atg.Position.Should().Be(0, "the tandem block starts at the first base");
        atg.TotalLength.Should().Be(9, "INV-02: TotalLength = |Unit| (3) × Repetitions (3)");
        (atg.Position + atg.TotalLength).Should().BeLessThanOrEqualTo("ATGATGATG".Length,
            "INV-03: the reported tandem stays within the sequence bounds");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// produce only well-formed results — no spurious counts, no out-of-range
    /// positions, no hang — so the degenerate-boundary guards do not corrupt the scan
    /// on ordinary input. Every result must satisfy INV-01..INV-03
    /// (Tandem_Repeat_Detection.md §2.4). Length kept modest because the detector is
    /// O(n²·m); a hang here would trip the timeout.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindTandemRepeats_RandomSequence_ProducesOnlyWellFormedResults()
    {
        const int minUnit = 1;
        const int minReps = 2;
        string seqStr = RandomDna(500, seed: 14_001);
        var seq = new DnaSequence(seqStr);

        var results = GenomicAnalyzer.FindTandemRepeats(seq, minUnit, minReps).ToList();

        results.Should().OnlyContain(r =>
            r.Repetitions >= minReps &&                          // INV-01 (count)
            r.Unit.Length >= minUnit &&                          // INV-01 (unit length)
            r.TotalLength == r.Unit.Length * r.Repetitions &&    // INV-02
            r.Position >= 0 && r.Position + r.TotalLength <= seqStr.Length, // INV-03
            "every result on random input is well-formed: count ≥ minReps, unit ≥ minUnit, total = unit×count, in bounds");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-INV-001 — inverted repeat detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region REP-INV-001 — inverted repeat detection

    #region BE — Boundary: minLen (minArmLength) = 0

    /// <summary>
    /// BE: minArmLength = 0 is the degenerate floor and the KEY nonsense-output target.
    /// A literal 0 means a zero-length arm — `leftArm = Substring(i, 0) = ""`, whose
    /// reverse complement is also "" — so EVERY downstream gap position would "match"
    /// the empty arm and the detector would blow up the result set with spurious
    /// zero-length-arm "repeats". A tandem stem needs ≥ 2 paired bases, so the contract
    /// REJECTS minArmLength &lt; 2 with ArgumentOutOfRangeException
    /// (RepeatFinder.cs line 282/300; Inverted_Repeat_Detection.md §3.3) — and, after
    /// this unit's fix, on BOTH the typed and the raw-string surface (the raw-string
    /// overload previously skipped this guard, Deviation #1). The string overload's
    /// check precedes its `yield break`, so we force enumeration to surface it.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindInvertedRepeats_MinArmLengthZero_ThrowsArgumentOutOfRange()
    {
        var typed = () => RepeatFinder.FindInvertedRepeats(new DnaSequence("GCGCAAAAGCGC"), minArmLength: 0).ToList();
        var raw = () => RepeatFinder.FindInvertedRepeats("GCGCAAAAGCGC", minArmLength: 0).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>(
                "minArmLength = 0 is below the documented floor of 2; a 0-length arm matches every gap as a spurious empty repeat")
            .Which.ParamName.Should().Be("minArmLength");
        raw.Should().Throw<ArgumentOutOfRangeException>(
                "the raw-string overload now enforces the same minArmLength >= 2 floor, never emitting nonsense empty-arm results")
            .Which.ParamName.Should().Be("minArmLength");
    }

    /// <summary>
    /// BE: minArmLength = 1 is the trivial single-base-arm boundary, still below the
    /// stem floor of 2 (a 1-bp "arm" is a single complementary base, not a stem). Both
    /// surfaces reject it with ArgumentOutOfRangeException, pinning that the rejection
    /// boundary is exactly &lt; 2 — and that minLoopLength &lt; 0 is rejected too.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_DegenerateParameters_RejectedOnBothSurfaces()
    {
        var typedArm1 = () => RepeatFinder.FindInvertedRepeats(new DnaSequence("GCGCAAAAGCGC"), minArmLength: 1).ToList();
        var rawArm1 = () => RepeatFinder.FindInvertedRepeats("GCGCAAAAGCGC", minArmLength: 1).ToList();
        var typedLoopNeg = () => RepeatFinder.FindInvertedRepeats(new DnaSequence("GCGCAAAAGCGC"), minArmLength: 4, minLoopLength: -1).ToList();
        var rawLoopNeg = () => RepeatFinder.FindInvertedRepeats("GCGCAAAAGCGC", minArmLength: 4, minLoopLength: -1).ToList();

        typedArm1.Should().Throw<ArgumentOutOfRangeException>("a 1-bp arm is below the documented stem floor of 2");
        rawArm1.Should().Throw<ArgumentOutOfRangeException>("the raw-string overload enforces the same minArmLength >= 2 floor");
        typedLoopNeg.Should().Throw<ArgumentOutOfRangeException>("a negative minLoopLength is nonsensical and rejected");
        rawLoopNeg.Should().Throw<ArgumentOutOfRangeException>("the raw-string overload now rejects a negative minLoopLength too");

        // The accepted floor minArmLength = 2 must NOT throw — pinning the boundary is at < 2, not ≤ 2.
        var arm2 = () => RepeatFinder.FindInvertedRepeats("GCGCAAAAGCGC", minArmLength: 2).ToList();
        arm2.Should().NotThrow("minArmLength = 2 is the accepted stem floor, not below it");
    }

    #endregion

    #region BE — Boundary: minLen > seqLen / 2

    /// <summary>
    /// BE: minArmLength &gt; seqLen/2 means two arms of that length cannot fit in the
    /// sequence even before counting the loop. The outer scan bound
    /// `i ≤ seq.Length − 2·minArmLength − minLoopLength` is negative, so the loop body
    /// never runs and no Substring is taken past the end — an empty result, never a
    /// crash (Inverted_Repeat_Detection.md §6.1). Probed with minArmLength = 8 on an
    /// 12-base sequence (2·8 = 16 ≫ 12) on both surfaces.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindInvertedRepeats_MinArmLengthExceedsHalfLength_IsEmptyAndDoesNotThrow()
    {
        const string seq = "GCGCAAAAGCGC"; // length 12; minArmLength 8 → 2·8 = 16 > 12
        var typedAct = () => RepeatFinder.FindInvertedRepeats(new DnaSequence(seq), minArmLength: 8).ToList();
        var rawAct = () => RepeatFinder.FindInvertedRepeats(seq, minArmLength: 8).ToList();

        typedAct.Should().NotThrow("an arm longer than half the sequence cannot fit twice; the scan bound is negative and the loop never runs");
        rawAct.Should().NotThrow("the raw-string overload is equally guarded against indexing past the sequence end");

        RepeatFinder.FindInvertedRepeats(new DnaSequence(seq), minArmLength: 8).Should().BeEmpty(
            "no two 8-bp arms plus a loop fit inside 12 bases; the oversized minArmLength yields nothing, not a crash");
        RepeatFinder.FindInvertedRepeats(seq, minArmLength: 8).Should().BeEmpty();
    }

    /// <summary>
    /// BE: the exact-fit edge — when 2·minArmLength + minLoopLength equals the sequence
    /// length the structure JUST fits and a genuine inverted repeat must still be found,
    /// pinning that the boundary guard rejects "too big" without also dropping the
    /// largest legal arm. "GCGCTTTGCGC" (length 11) holds arms 'GCGC' (revcomp 'GCGC')
    /// with a 3-base 'TTT' loop: 2·4 + 3 = 11. minArmLength = 4 must find it.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_ArmLengthAtExactFit_StillFindsRepeat()
    {
        var results = RepeatFinder.FindInvertedRepeats("GCGCTTTGCGC", minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        results.Should().Contain(
            r => r.LeftArm == "GCGC" && r.ArmLength == 4 && r.LoopLength == 3,
            "the 'GCGC…GCGC' hairpin fits exactly (2·4 + 3 = 11) and is found at the fit boundary");
    }

    #endregion

    #region BE — Boundary: no complement possibilities (homopolymer)

    /// <summary>
    /// BE: a sequence with "no complement possibilities" — an all-A homopolymer — can
    /// never contain an inverted repeat, because the reverse complement of any all-A
    /// arm is an all-T arm that is absent from the same run (Inverted_Repeat_Detection.md
    /// §6.1, homopolymer edge case). The detector must return empty with no crash and no
    /// hang. Pinned on both surfaces and for a long run to exercise the full scan.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void FindInvertedRepeats_HomopolymerNoComplement_IsEmptyAndDoesNotThrow()
    {
        string allA = new string('A', 60);

        var typedAct = () => RepeatFinder.FindInvertedRepeats(new DnaSequence(allA), minArmLength: 4).ToList();
        var rawAct = () => RepeatFinder.FindInvertedRepeats(allA, minArmLength: 4).ToList();

        typedAct.Should().NotThrow("revcomp of an all-A arm is all-T, which never appears in an all-A run; the scan finds nothing");
        rawAct.Should().NotThrow();

        RepeatFinder.FindInvertedRepeats(new DnaSequence(allA), minArmLength: 4).Should().BeEmpty(
            "a homopolymer has no reverse-complement arm downstream — no inverted repeat, no crash");
        RepeatFinder.FindInvertedRepeats(allA, minArmLength: 4).Should().BeEmpty();
    }

    /// <summary>
    /// BE: a heteropolymer that simply contains no downstream reverse-complement arm —
    /// a strictly increasing-then-non-complementary run — must also return empty without
    /// crashing (Inverted_Repeat_Detection.md §6.1, "no complementary regions"). Pins
    /// that "found nothing" is a clean empty result, distinct from the rejection cases.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_NoComplementaryArm_IsEmpty()
    {
        // 'AAAACCCC' — no arm's revcomp appears downstream within the loop window.
        RepeatFinder.FindInvertedRepeats("AAAACCCC", minArmLength: 4, maxLoopLength: 10, minLoopLength: 3)
            .Should().BeEmpty("no downstream substring equals a reverse-complement candidate; the result is cleanly empty");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The raw-string overload
    /// short-circuits null/empty to the empty enumerable via `yield break`
    /// (RepeatFinder.cs lines 300–301) — no exception, even though the parameter
    /// validation now runs first (valid params here, so it passes through). The typed
    /// overload over an empty DnaSequence has a negative outer scan bound, so it yields
    /// nothing. Neither path divides, indexes, or hangs on empty input
    /// (Inverted_Repeat_Detection.md §6.1).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindInvertedRepeats_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var typedAct = () => RepeatFinder.FindInvertedRepeats(new DnaSequence(string.Empty), minArmLength: 4).ToList();
        var rawEmptyAct = () => RepeatFinder.FindInvertedRepeats(string.Empty, minArmLength: 4).ToList();
        var rawNullAct = () => RepeatFinder.FindInvertedRepeats((string)null!, minArmLength: 4).ToList();

        typedAct.Should().NotThrow("an empty sequence has no room for two arms and a loop; it yields nothing");
        rawEmptyAct.Should().NotThrow("the raw-string overload short-circuits empty input to an empty result");
        rawNullAct.Should().NotThrow("the raw-string overload treats null input as empty, not as an error");

        RepeatFinder.FindInvertedRepeats(new DnaSequence(string.Empty), minArmLength: 4).Should().BeEmpty();
        RepeatFinder.FindInvertedRepeats(string.Empty, minArmLength: 4).Should().BeEmpty();
        RepeatFinder.FindInvertedRepeats((string)null!, minArmLength: 4).Should().BeEmpty();
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". The typed overload
    /// guards it with an explicit ArgumentNullException (ThrowIfNull, RepeatFinder.cs
    /// line 281), raised eagerly at the call — never a NullReferenceException.
    /// </summary>
    [Test]
    public void FindInvertedRepeats_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => RepeatFinder.FindInvertedRepeats((DnaSequence)null!, minArmLength: 4);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload null-guards its sequence; null is rejected, never dereferenced");
    }

    #endregion

    #region Positive sanity — a clear inverted repeat / hairpin is detected correctly

    /// <summary>
    /// Positive sanity: a textbook self-complementary hairpin — "GAATTC" (the EcoRI
    /// palindrome, revcomp(GAATTC) = GAATTC) repeated across a 4-base loop in
    /// "GAATTCAAAAGAATTC" — must be detected with the CORRECT arm coordinates, so the
    /// boundary hardening never silently breaks the core function. Pinned per
    /// INV-01..INV-04 (Inverted_Repeat_Detection.md §2.4): left arm 'GAATTC' at 0,
    /// right arm 'GAATTC' at 10, arm length 6, loop length 4, TotalLength = 2·6 + 4 = 16,
    /// CanFormHairpin true (loop ≥ 3).
    /// </summary>
    [Test]
    public void FindInvertedRepeats_EcoRiHairpin_DetectedWithCorrectArmCoords()
    {
        var results = RepeatFinder.FindInvertedRepeats("GAATTCAAAAGAATTC", minArmLength: 4, maxLoopLength: 10, minLoopLength: 3).ToList();

        var hairpin = results.Should().ContainSingle(r => r.ArmLength == 6).Subject;
        hairpin.LeftArmStart.Should().Be(0, "the left arm starts at the first base");
        hairpin.RightArmStart.Should().Be(10, "the right arm starts after the 6-bp arm and 4-bp loop");
        hairpin.LeftArm.Should().Be("GAATTC");
        hairpin.RightArm.Should().Be("GAATTC", "INV-01: revcomp(GAATTC) = GAATTC (self-complementary EcoRI site)");
        DnaSequence.GetReverseComplementString(hairpin.LeftArm).Should().Be(hairpin.RightArm,
            "INV-01: the right arm equals the reverse complement of the left arm");
        hairpin.LoopLength.Should().Be(4, "INV-03: loop = rightStart(10) − (leftStart(0) + armLen(6))");
        hairpin.TotalLength.Should().Be(16, "INV-02: TotalLength = 2·ArmLength(6) + LoopLength(4)");
        hairpin.CanFormHairpin.Should().BeTrue("INV-04: a loop of 4 ≥ 3 is hairpin-viable");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// produce only well-formed results — every reported arm pair is a true reverse
    /// complement, every coordinate is in bounds, and the loop/total invariants hold —
    /// so the degenerate-boundary guards do not corrupt the scan on ordinary input.
    /// Pinned per INV-01..INV-04 (Inverted_Repeat_Detection.md §2.4). Length kept modest
    /// because the detector is O(n·A²·L); a hang would trip the timeout.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindInvertedRepeats_RandomSequence_ProducesOnlyWellFormedResults()
    {
        string seq = RandomDna(400, seed: 15_001);

        var results = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4, maxLoopLength: 20, minLoopLength: 3).ToList();

        results.Should().OnlyContain(r =>
            DnaSequence.GetReverseComplementString(r.LeftArm) == r.RightArm &&        // INV-01
            r.TotalLength == 2 * r.ArmLength + r.LoopLength &&                          // INV-02
            r.LoopLength == r.RightArmStart - (r.LeftArmStart + r.ArmLength) &&         // INV-03
            r.CanFormHairpin == (r.LoopLength >= 3) &&                                  // INV-04
            r.ArmLength >= 4 && r.LeftArmStart >= 0 &&
            r.RightArmStart + r.ArmLength <= seq.Length,
            "every result on random input is well-formed: arms are reverse complements, total/loop invariants hold, coords in bounds");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-DIRECT-001 — direct repeat detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region REP-DIRECT-001 — direct repeat detection

    #region BE — Boundary: minLen (minLength) = 0

    /// <summary>
    /// BE: minLength = 0 is the degenerate floor and the KEY nonsense-output target. A
    /// literal 0 makes a zero-length candidate — `repeat = Substring(i, 0) = ""` — and the
    /// suffix-tree lookup `FindAllOccurrences("")` returns EVERY start position, so the empty
    /// candidate "matches" at every downstream gap and the detector blows the result set up
    /// with O(n²) spurious empty-length "repeats". A direct repeat needs ≥ 2 bases, so the
    /// contract REJECTS minLength &lt; 2 with ArgumentOutOfRangeException
    /// (RepeatFinder.cs line 376/391; Direct_Repeat_Detection.md §3.3) — and, after this
    /// unit's fix, on BOTH the typed and the raw-string surface (the raw-string overload
    /// previously skipped this guard, Deviation #1). Both surfaces validate eagerly, so the
    /// throw surfaces at the call even before enumeration; we still force enumeration so a
    /// regression to lazy/late validation would also be caught.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindDirectRepeats_MinLengthZero_ThrowsArgumentOutOfRange()
    {
        var typed = () => RepeatFinder.FindDirectRepeats(new DnaSequence("ACGTACGTACGT"), minLength: 0).ToList();
        var raw = () => RepeatFinder.FindDirectRepeats("ACGTACGTACGT", minLength: 0).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>(
                "minLength = 0 is below the documented floor of 2; a 0-length candidate matches every gap as a spurious empty repeat")
            .Which.ParamName.Should().Be("minLength");
        raw.Should().Throw<ArgumentOutOfRangeException>(
                "the raw-string overload now enforces the same minLength >= 2 floor, never emitting the O(n^2) empty-repeat blow-up")
            .Which.ParamName.Should().Be("minLength");
    }

    #endregion

    #region BE — Boundary: minLen (minLength) = 1

    /// <summary>
    /// BE: minLength = 1 is the trivial single-base-repeat boundary, still below the floor of
    /// 2 (a length-1 "repeat" is a single base that happens to recur — every common base would
    /// qualify, a defined-but-useless blow-up). Both surfaces REJECT it with
    /// ArgumentOutOfRangeException, pinning that the rejection boundary is exactly &lt; 2 and
    /// cannot drift to ≤ 1. The accepted floor minLength = 2 must NOT throw — pinned alongside
    /// so the boundary is fixed at exactly 2.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MinLengthOne_RejectedOnBothSurfaces_FloorIsExactlyTwo()
    {
        var typedOne = () => RepeatFinder.FindDirectRepeats(new DnaSequence("ACGTACGTACGT"), minLength: 1).ToList();
        var rawOne = () => RepeatFinder.FindDirectRepeats("ACGTACGTACGT", minLength: 1).ToList();

        typedOne.Should().Throw<ArgumentOutOfRangeException>(
            "minLength = 1 is below the documented floor of 2; every recurring single base would be a trivial repeat");
        rawOne.Should().Throw<ArgumentOutOfRangeException>(
            "the raw-string overload enforces the same minLength >= 2 floor");

        // The accepted floor minLength = 2 must NOT throw — pinning the boundary is at < 2, not ≤ 2.
        var typedTwo = () => RepeatFinder.FindDirectRepeats(new DnaSequence("ACGTACGTACGT"), minLength: 2, maxLength: 4).ToList();
        var rawTwo = () => RepeatFinder.FindDirectRepeats("ACGTACGTACGT", minLength: 2, maxLength: 4).ToList();
        typedTwo.Should().NotThrow("minLength = 2 is the accepted floor, not below it");
        rawTwo.Should().NotThrow("the raw-string overload accepts the floor of 2 just as the typed overload does");
    }

    /// <summary>
    /// BE: a maxLength below minLength is the inverted-range boundary. Both surfaces reject it
    /// with ArgumentOutOfRangeException (RepeatFinder.cs line 377/392), pinning that the range
    /// guard mirrored onto the raw-string overload fires there too and not only on the typed
    /// surface.
    /// </summary>
    [Test]
    public void FindDirectRepeats_MaxLengthBelowMinLength_RejectedOnBothSurfaces()
    {
        var typed = () => RepeatFinder.FindDirectRepeats(new DnaSequence("ACGTACGTACGT"), minLength: 8, maxLength: 4).ToList();
        var raw = () => RepeatFinder.FindDirectRepeats("ACGTACGTACGT", minLength: 8, maxLength: 4).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>("maxLength < minLength is an inverted range and rejected")
            .Which.ParamName.Should().Be("maxLength");
        raw.Should().Throw<ArgumentOutOfRangeException>("the raw-string overload now rejects the inverted range too")
            .Which.ParamName.Should().Be("maxLength");
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The raw-string overload short-circuits
    /// null/empty to the empty enumerable via `yield break` (RepeatFinder.cs lines 391–392) — no
    /// exception (with valid params, which pass the now-eager numeric gate). The typed overload
    /// over an empty DnaSequence has a position-loop bound that is negative, so it yields nothing.
    /// Neither path indexes past the end or hangs on empty input (Direct_Repeat_Detection.md §6.1).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindDirectRepeats_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var typedAct = () => RepeatFinder.FindDirectRepeats(new DnaSequence(string.Empty), minLength: 5).ToList();
        var rawEmptyAct = () => RepeatFinder.FindDirectRepeats(string.Empty, minLength: 5).ToList();
        var rawNullAct = () => RepeatFinder.FindDirectRepeats((string)null!, minLength: 5).ToList();

        typedAct.Should().NotThrow("an empty sequence has no room for two repeat copies; it yields nothing");
        rawEmptyAct.Should().NotThrow("the raw-string overload short-circuits empty input to an empty result");
        rawNullAct.Should().NotThrow("the raw-string overload treats null input as empty, not as an error");

        RepeatFinder.FindDirectRepeats(new DnaSequence(string.Empty), minLength: 5).Should().BeEmpty();
        RepeatFinder.FindDirectRepeats(string.Empty, minLength: 5).Should().BeEmpty();
        RepeatFinder.FindDirectRepeats((string)null!, minLength: 5).Should().BeEmpty();
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". The typed overload guards
    /// it with an explicit ArgumentNullException (ThrowIfNull, RepeatFinder.cs line 375), raised
    /// eagerly at the call — never a NullReferenceException.
    /// </summary>
    [Test]
    public void FindDirectRepeats_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => RepeatFinder.FindDirectRepeats((DnaSequence)null!, minLength: 5);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload null-guards its sequence; null is rejected, never dereferenced");
    }

    #endregion

    #region BE — Boundary: all-unique characters (no recurrence)

    /// <summary>
    /// BE: an all-unique-characters sequence — every base distinct, e.g. "ACGT" — can contain no
    /// direct repeat: no subsequence of length ≥ 2 recurs anywhere downstream. The detector must
    /// return a CLEAN empty result with no crash and no hang, distinct from the rejection cases
    /// above. Pinned on both surfaces and at the smallest legal repeat length (2).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindDirectRepeats_AllUniqueChars_IsEmptyAndDoesNotThrow()
    {
        const string unique = "ACGT"; // four distinct bases — nothing recurs

        var typedAct = () => RepeatFinder.FindDirectRepeats(new DnaSequence(unique), minLength: 2, maxLength: 4).ToList();
        var rawAct = () => RepeatFinder.FindDirectRepeats(unique, minLength: 2, maxLength: 4).ToList();

        typedAct.Should().NotThrow("no subsequence recurs in an all-unique sequence; the scan finds nothing");
        rawAct.Should().NotThrow();

        RepeatFinder.FindDirectRepeats(new DnaSequence(unique), minLength: 2, maxLength: 4).Should().BeEmpty(
            "no length ≥ 2 subsequence recurs in 'ACGT' — no direct repeat, no crash");
        RepeatFinder.FindDirectRepeats(unique, minLength: 2, maxLength: 4).Should().BeEmpty();
    }

    #endregion

    #region BE — Boundary: all-same character (homopolymer — maximal recurrence)

    /// <summary>
    /// BE: an all-same-character homopolymer (e.g. "AAAA…") is the MAXIMAL-recurrence /
    /// combinatorial-blow-up watch point — every length-L window equals every other, so direct
    /// repeats genuinely abound. The detector must NOT hang or blow up unboundedly: the (i, j, len)
    /// dedup hash set (INV-04) keeps the count polynomial, and the scan must complete promptly.
    /// We pin (a) it completes within the timeout, (b) every reported pair is a real same-character
    /// repeat satisfying the documented invariants, and (c) for a homopolymer with minSpacing = 0 a
    /// known adjacent pair (e.g. "AA" at 0 abutting "AA" at 2) IS found — the boundary produces
    /// correct, bounded output, not nonsense.
    /// </summary>
    [Test]
    [CancelAfter(10000)]
    public void FindDirectRepeats_HomopolymerMaximalRecurrence_CompletesWithBoundedWellFormedResults()
    {
        string allA = new string('A', 40);

        var act = () => RepeatFinder.FindDirectRepeats(allA, minLength: 2, maxLength: 6, minSpacing: 0).ToList();
        act.Should().NotThrow("a homopolymer is the maximal-recurrence case but the (i,j,len) dedup keeps it bounded — no hang, no crash");

        var results = RepeatFinder.FindDirectRepeats(allA, minLength: 2, maxLength: 6, minSpacing: 0).ToList();

        results.Should().NotBeEmpty("a homopolymer trivially contains many same-orientation direct repeats");
        results.Should().OnlyContain(r =>
            r.RepeatSequence == allA.Substring(r.FirstPosition, r.Length) &&         // INV-01 (copy 1)
            r.RepeatSequence == allA.Substring(r.SecondPosition, r.Length) &&        // INV-01 (copy 2)
            r.RepeatSequence.All(c => c == 'A') &&                                   // every base is the homopolymer char
            r.Length >= 2 && r.Length <= 6 &&                                        // within the tested length band
            r.Spacing == r.SecondPosition - r.FirstPosition - r.Length &&            // INV-02
            r.SecondPosition > r.FirstPosition &&                                    // downstream copy
            r.SecondPosition + r.Length <= allA.Length,                             // in bounds
            "every homopolymer result is a well-formed same-character direct-repeat pair, not a spurious or out-of-range entry");

        // INV-04: each (FirstPosition, SecondPosition, Length) tuple is reported at most once.
        results.Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).Should()
            .OnlyHaveUniqueItems("INV-04: the (i, j, len) dedup hash set suppresses duplicate result keys");

        // A concrete known pair: 'AA' at 0 and 'AA' at 2 (adjacent, spacing 0) must be present.
        results.Should().Contain(
            r => r.Length == 2 && r.FirstPosition == 0 && r.SecondPosition == 2 && r.Spacing == 0,
            "with minSpacing = 0 the adjacent 'AA'…'AA' pair is a legitimate direct repeat and is reported");
    }

    #endregion

    #region Positive sanity — a clear direct repeat is detected correctly

    /// <summary>
    /// Positive sanity: a textbook direct repeat — the motif "ATCG" recurring verbatim across a
    /// 4-base spacer in "ATCGTTTTATCG" — must be detected with the CORRECT coordinates, so the
    /// boundary hardening never silently breaks the core function. Pinned per INV-01..INV-04
    /// (Direct_Repeat_Detection.md §2.4): first copy at 0, second copy at 8, length 4, the two
    /// copies identical, Spacing = 8 − 0 − 4 = 4.
    /// </summary>
    [Test]
    public void FindDirectRepeats_ClearDirectRepeat_DetectedWithCorrectCoords()
    {
        const string seq = "ATCGTTTTATCG"; // 'ATCG' at 0 and at 8, spacer 'TTTT'

        var results = RepeatFinder.FindDirectRepeats(seq, minLength: 4, maxLength: 4, minSpacing: 1).ToList();

        var pair = results.Should().ContainSingle(r => r.RepeatSequence == "ATCG").Subject;
        pair.FirstPosition.Should().Be(0, "the first copy of 'ATCG' starts at the first base");
        pair.SecondPosition.Should().Be(8, "the second copy of 'ATCG' starts after the 4-bp motif and 4-bp spacer");
        pair.Length.Should().Be(4);
        seq.Substring(pair.FirstPosition, pair.Length).Should().Be(seq.Substring(pair.SecondPosition, pair.Length),
            "INV-01: the two copies are identical at both positions");
        pair.Spacing.Should().Be(4, "INV-02: Spacing = SecondPosition(8) − FirstPosition(0) − Length(4)");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and produce only
    /// well-formed results — every reported pair is a verbatim same-orientation recurrence, the
    /// spacing invariant holds, copies are non-overlapping under minSpacing &gt; 0, and all keys are
    /// unique — so the degenerate-boundary guards do not corrupt the scan on ordinary input. Pinned
    /// per INV-01..INV-04 (Direct_Repeat_Detection.md §2.4). Length kept modest because detection is
    /// O(r·n·(m+k)); a hang would trip the timeout.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindDirectRepeats_RandomSequence_ProducesOnlyWellFormedResults()
    {
        string seq = RandomDna(400, seed: 16_001);

        var results = RepeatFinder.FindDirectRepeats(seq, minLength: 5, maxLength: 20, minSpacing: 1).ToList();

        results.Should().OnlyContain(r =>
            r.RepeatSequence == seq.Substring(r.FirstPosition, r.Length) &&          // INV-01 (copy 1)
            r.RepeatSequence == seq.Substring(r.SecondPosition, r.Length) &&         // INV-01 (copy 2)
            r.Spacing == r.SecondPosition - r.FirstPosition - r.Length &&            // INV-02
            r.SecondPosition >= r.FirstPosition + r.Length + 1 &&                    // INV-03 (minSpacing = 1 → no overlap)
            r.Length >= 5 && r.Length <= 20 &&
            r.FirstPosition >= 0 && r.SecondPosition + r.Length <= seq.Length,
            "every result on random input is well-formed: copies are verbatim recurrences, spacing/bounds invariants hold, no overlap");

        results.Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).Should()
            .OnlyHaveUniqueItems("INV-04: every (FirstPosition, SecondPosition, Length) key is unique");
    }

    #endregion

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  REP-PALIN-001 — DNA palindrome detection : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region REP-PALIN-001 — palindrome detection

    #region BE — Boundary: minLen (minLength) = 0

    /// <summary>
    /// BE: minLength = 0 is the degenerate floor and the KEY nonsense-output target. A
    /// literal 0 means a zero-length window — `candidate = Substring(i, 0) = ""`, whose
    /// reverse complement is also "" — so EVERY position would "match" the empty
    /// palindrome and the detector would blow up the result set with spurious
    /// zero-length "palindromes". A biological palindrome needs ≥ 4 paired bases
    /// (the validating API excludes trivial sub-4 windows), so the contract REJECTS
    /// minLength &lt; 4 with ArgumentOutOfRangeException (RepeatFinder.cs lines 559/575;
    /// Palindrome_Detection.md §3.3, §6.1) on BOTH the typed and the raw-string surface.
    /// The raw-string overload's check lives inside the iterator (before its
    /// `yield break`), so we force enumeration to surface it.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPalindromes_MinLengthZero_ThrowsArgumentOutOfRange()
    {
        var typed = () => RepeatFinder.FindPalindromes(new DnaSequence("GAATTCGAATTC"), minLength: 0).ToList();
        var raw = () => RepeatFinder.FindPalindromes("GAATTCGAATTC", minLength: 0).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>(
                "minLength = 0 is below the documented even floor of 4; a 0-length window matches every position as a spurious empty palindrome")
            .Which.ParamName.Should().Be("minLength");
        raw.Should().Throw<ArgumentOutOfRangeException>(
                "the raw-string overload enforces the same minLength >= 4 even floor before scanning")
            .Which.ParamName.Should().Be("minLength");
    }

    #endregion

    #region BE — Boundary: minLen (minLength) is ODD

    /// <summary>
    /// BE — THE distinctive target for this unit. A DNA reverse-complement palindrome
    /// is self-complementary across a central symmetry axis, so every base pairs with
    /// a complementary partner and the length is ALWAYS EVEN (Palindrome_Detection.md
    /// §2.2, INV-02). An ODD minLength is therefore biologically impossible. The
    /// validating API does NOT round the request up/down and does NOT silently emit
    /// only even results — it REJECTS the odd value OUTRIGHT with
    /// ArgumentOutOfRangeException (the `minLength % 2 != 0` clause, RepeatFinder.cs
    /// lines 559/575; Palindrome_Detection.md §6.1 "Odd minLength → throws"). We pin
    /// that behavior on BOTH surfaces, across several odd values (one below and several
    /// above the floor of 4), so the contract cannot silently drift to rounding.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    [TestCase(3)]
    [TestCase(5)]
    [TestCase(7)]
    [TestCase(11)]
    public void FindPalindromes_MinLengthOdd_ThrowsArgumentOutOfRange(int oddMinLength)
    {
        var typed = () => RepeatFinder.FindPalindromes(new DnaSequence("GAATTCGAATTC"), minLength: oddMinLength, maxLength: 12).ToList();
        var raw = () => RepeatFinder.FindPalindromes("GAATTCGAATTC", minLength: oddMinLength, maxLength: 12).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>(
                "a DNA reverse-complement palindrome must have EVEN length; an odd minLength is rejected, never rounded")
            .Which.ParamName.Should().Be("minLength");
        raw.Should().Throw<ArgumentOutOfRangeException>(
                "the raw-string overload enforces the same even-length requirement on minLength")
            .Which.ParamName.Should().Be("minLength");
    }

    /// <summary>
    /// BE: the accepted floor minLength = 4 — even and at the boundary — must NOT throw,
    /// pinning that the rejection boundary is exactly "&lt; 4 OR odd", not "≤ 4". A clean
    /// even request of 4 is accepted and behaves. Paired with the odd-rejection cases
    /// this fixes the boundary at exactly even-and-≥-4.
    /// </summary>
    [Test]
    public void FindPalindromes_MinLengthFourEven_IsAcceptedOnBothSurfaces()
    {
        var typed = () => RepeatFinder.FindPalindromes(new DnaSequence("GAATTC"), minLength: 4, maxLength: 6).ToList();
        var raw = () => RepeatFinder.FindPalindromes("GAATTC", minLength: 4, maxLength: 6).ToList();

        typed.Should().NotThrow("minLength = 4 is the accepted even floor, not below it");
        raw.Should().NotThrow("the raw-string overload accepts the even floor of 4 just as the typed overload does");
    }

    #endregion

    #region BE — Boundary: maxLength > sequence length

    /// <summary>
    /// BE: maxLength far larger than the sequence length must NOT crash and must never
    /// invent a palindrome longer than the sequence. For every candidate length larger
    /// than the sequence the inner scan bound `i ≤ seq.Length − len` is negative, so the
    /// loop body never runs and no Substring is taken past the end (RepeatFinder.cs lines
    /// 592–606; Palindrome_Detection.md §2.4 INV-03, §6.1). "GAATTC" (length 6) with
    /// maxLength = 1000 ≫ 6 must still find ONLY its own 6-bp palindrome and nothing
    /// longer. Pinned on both surfaces.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPalindromes_MaxLengthExceedsSequence_DoesNotCrashAndNeverExceedsSeqLength()
    {
        const string seq = "GAATTC"; // length 6; maxLength 1000 ≫ 6

        var typedAct = () => RepeatFinder.FindPalindromes(new DnaSequence(seq), minLength: 4, maxLength: 1000).ToList();
        var rawAct = () => RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 1000).ToList();

        typedAct.Should().NotThrow("an oversized maxLength makes the scan bound negative for the big lengths; the loop simply never runs");
        rawAct.Should().NotThrow("the raw-string overload is equally guarded against indexing past the sequence end");

        var results = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 1000).ToList();

        results.Should().Contain(r => r.Sequence == "GAATTC" && r.Length == 6,
            "the genuine 6-bp EcoRI palindrome is still found despite maxLength ≫ sequence length");
        results.Should().OnlyContain(r => r.Length <= seq.Length,
            "INV-03: no reported palindrome can be longer than the sequence, even when maxLength far exceeds it");
        results.Should().OnlyContain(r => r.Length % 2 == 0,
            "INV-02: every reported palindrome length is even, even at the oversized-maxLength boundary");
    }

    /// <summary>
    /// BE: maxLength &gt; seqLen on a sequence too SHORT to hold even the minimum window
    /// must yield a clean empty result, not a crash. "GC" (length 2) with the default
    /// even floor minLength = 4 and an oversized maxLength has no window of length ≥ 4
    /// that fits, so the loop never runs — empty, no out-of-range Substring.
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPalindromes_MaxLengthExceedsShortSequence_IsEmptyAndDoesNotThrow()
    {
        var typedAct = () => RepeatFinder.FindPalindromes(new DnaSequence("GC"), minLength: 4, maxLength: 500).ToList();
        var rawAct = () => RepeatFinder.FindPalindromes("GC", minLength: 4, maxLength: 500).ToList();

        typedAct.Should().NotThrow("no length-≥-4 window fits a 2-base sequence; the scan bound is negative and the loop never runs");
        rawAct.Should().NotThrow();

        RepeatFinder.FindPalindromes(new DnaSequence("GC"), minLength: 4, maxLength: 500).Should().BeEmpty(
            "a 2-base sequence cannot hold a ≥ 4-bp palindrome; the oversized max yields nothing, not a crash");
        RepeatFinder.FindPalindromes("GC", minLength: 4, maxLength: 500).Should().BeEmpty();
    }

    #endregion

    #region BE — Boundary: empty sequence

    /// <summary>
    /// BE: the empty sequence is the lower size boundary. The raw-string overload
    /// short-circuits null/empty to the empty enumerable via `yield break`
    /// (RepeatFinder.cs lines 580–581) — no exception, even though the numeric gate
    /// runs first (valid lengths here, so it passes through). The typed overload over
    /// an empty DnaSequence (the ctor materializes an empty sequence from "") has a
    /// negative scan bound, so it yields nothing. Neither path divides, indexes, or
    /// hangs on empty input (Palindrome_Detection.md §6.1).
    /// </summary>
    [Test]
    [CancelAfter(5000)]
    public void FindPalindromes_EmptySequence_IsEmptyAndDoesNotThrow()
    {
        var typedAct = () => RepeatFinder.FindPalindromes(new DnaSequence(string.Empty), minLength: 4, maxLength: 12).ToList();
        var rawEmptyAct = () => RepeatFinder.FindPalindromes(string.Empty, minLength: 4, maxLength: 12).ToList();
        var rawNullAct = () => RepeatFinder.FindPalindromes((string)null!, minLength: 4, maxLength: 12).ToList();

        typedAct.Should().NotThrow("an empty sequence has no window long enough to hold a palindrome; it yields nothing");
        rawEmptyAct.Should().NotThrow("the raw-string overload short-circuits empty input to an empty result");
        rawNullAct.Should().NotThrow("the raw-string overload treats null input as empty, not as an error");

        RepeatFinder.FindPalindromes(new DnaSequence(string.Empty), minLength: 4, maxLength: 12).Should().BeEmpty();
        RepeatFinder.FindPalindromes(string.Empty, minLength: 4, maxLength: 12).Should().BeEmpty();
        RepeatFinder.FindPalindromes((string)null!, minLength: 4, maxLength: 12).Should().BeEmpty();
    }

    /// <summary>
    /// BE/INJ: a null DnaSequence is the boundary of "no typed input". The typed
    /// overload guards it with an explicit ArgumentNullException (ThrowIfNull,
    /// RepeatFinder.cs line 558), raised eagerly at the call — never a
    /// NullReferenceException.
    /// </summary>
    [Test]
    public void FindPalindromes_NullDnaSequence_ThrowsArgumentNullException()
    {
        var act = () => RepeatFinder.FindPalindromes((DnaSequence)null!, minLength: 4, maxLength: 12);

        act.Should().Throw<ArgumentNullException>(
            "the typed overload null-guards its sequence; null is rejected, never dereferenced");
    }

    /// <summary>
    /// BE: maxLength below minLength is the inverted-range boundary. Both surfaces reject
    /// it with ArgumentOutOfRangeException (RepeatFinder.cs lines 561/577), pinning that
    /// the ordered-bounds guard fires on both surfaces.
    /// </summary>
    [Test]
    public void FindPalindromes_MaxLengthBelowMinLength_RejectedOnBothSurfaces()
    {
        var typed = () => RepeatFinder.FindPalindromes(new DnaSequence("GAATTCGAATTC"), minLength: 8, maxLength: 4).ToList();
        var raw = () => RepeatFinder.FindPalindromes("GAATTCGAATTC", minLength: 8, maxLength: 4).ToList();

        typed.Should().Throw<ArgumentOutOfRangeException>("maxLength < minLength is an inverted range and rejected")
            .Which.ParamName.Should().Be("maxLength");
        raw.Should().Throw<ArgumentOutOfRangeException>("the raw-string overload rejects the inverted range too")
            .Which.ParamName.Should().Be("maxLength");
    }

    #endregion

    #region Positive sanity — a known restriction-site palindrome is detected correctly

    /// <summary>
    /// Positive sanity: the textbook restriction-site palindrome — "GAATTC" (the EcoRI
    /// recognition sequence, length 6) — must be detected with the CORRECT sequence,
    /// position and length, so the boundary hardening never silently breaks the core
    /// function. We first PROVE the biological-palindrome property the detector relies
    /// on — revcomp(GAATTC) == GAATTC — then pin INV-01..INV-04
    /// (Palindrome_Detection.md §2.4): Sequence 'GAATTC', Position 0, Length 6 (even),
    /// and the reported sequence equals its own reverse complement.
    /// </summary>
    [Test]
    public void FindPalindromes_EcoRiSite_DetectedAsSelfComplementaryPalindrome()
    {
        // Prove the biological-palindrome property the detector hinges on.
        DnaSequence.GetReverseComplementString("GAATTC").Should().Be("GAATTC",
            "EcoRI 'GAATTC' is self-complementary: revcomp(GAATTC) = GAATTC");

        var results = RepeatFinder.FindPalindromes("GAATTC", minLength: 4, maxLength: 12).ToList();

        var ecori = results.Should().ContainSingle(r => r.Length == 6).Subject;
        ecori.Sequence.Should().Be("GAATTC", "the 6-bp EcoRI site is the palindrome present");
        ecori.Position.Should().Be(0, "the palindrome starts at the first base");
        ecori.Length.Should().Be(6, "INV-04: Length equals the actual palindrome length");
        (ecori.Length % 2).Should().Be(0, "INV-02: a biological DNA palindrome has even length");
        DnaSequence.GetReverseComplementString(ecori.Sequence).Should().Be(ecori.Sequence,
            "INV-01: every reported palindrome equals its own reverse complement");
    }

    /// <summary>
    /// Positive sanity: a 4-bp restriction-site palindrome — "AGCT" (the AluI site,
    /// revcomp(AGCT) = AGCT) — at the smallest legal even floor must be detected, pinning
    /// that the floor-4 boundary still finds genuine palindromes. Embedded with flanks
    /// ("TTAGCTAA") so the detector locates it at an interior position.
    /// </summary>
    [Test]
    public void FindPalindromes_AluiSiteAtFloor_DetectedAtCorrectPosition()
    {
        DnaSequence.GetReverseComplementString("AGCT").Should().Be("AGCT",
            "AluI 'AGCT' is self-complementary: revcomp(AGCT) = AGCT");

        var results = RepeatFinder.FindPalindromes("TTAGCTAA", minLength: 4, maxLength: 4).ToList();

        var alui = results.Should().ContainSingle(r => r.Sequence == "AGCT").Subject;
        alui.Position.Should().Be(2, "the 'AGCT' palindrome sits after the 2-base 'TT' flank");
        alui.Length.Should().Be(4, "INV-04: the 4-bp AluI site is reported at the even floor");
    }

    /// <summary>
    /// Positive sanity / RB: a fixed-seed random sequence must complete promptly and
    /// produce only well-formed results — every reported window genuinely equals its own
    /// reverse complement, every length is even, and every coordinate is in bounds — so
    /// the degenerate-boundary guards do not corrupt the scan on ordinary input. Pinned
    /// per INV-01..INV-04 (Palindrome_Detection.md §2.4). Length kept modest because the
    /// scan is O(n·r·m); a hang would trip the timeout.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void FindPalindromes_RandomSequence_ProducesOnlyWellFormedResults()
    {
        string seq = RandomDna(2000, seed: 17_001);

        var results = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12).ToList();

        results.Should().OnlyContain(r =>
            DnaSequence.GetReverseComplementString(r.Sequence) == r.Sequence &&     // INV-01
            r.Length % 2 == 0 &&                                                    // INV-02
            r.Length >= 4 && r.Length <= 12 &&                                      // within tested band
            r.Position >= 0 && r.Position + r.Length <= seq.Length &&               // INV-03
            r.Sequence.Length == r.Length &&                                        // INV-04
            r.Sequence == seq.Substring(r.Position, r.Length),
            "every result on random input is well-formed: it equals its own reverse complement, has even length, and sits in bounds");
    }

    #endregion

    #endregion
}
